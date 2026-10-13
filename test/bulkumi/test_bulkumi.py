#!/usr/bin/env python3
"""End-to-end test for ESAT -umiMode {gene,position} and -umiDropN.

Builds a synthetic genome (chrT), a 3-gene mapping file and a coordinate-sorted BAM whose
reads are designed so that the expected number of collapsed molecules per gene is known for
every mode. Runs the ESAT jar under four configurations and checks
  (1) the number of collapsed UMI keys per gene in <exp>.umi.distributions.txt
  (2) the gene counts in <out>.gene.txt
against the expectations.

Usage: python test_bulkumi.py <esat.jar> [parent_dir_for_a_new_workdir]
"""
import os, re, subprocess, sys, collections, shutil

import tempfile
JAR = os.path.abspath(sys.argv[1])
# Always work in a NEW directory: inside the optional parent dir if given, else the system temp dir.
# The test never deletes a path it did not create.
parent = sys.argv[2] if len(sys.argv) > 2 and sys.argv[2] else None
if parent:
    os.makedirs(parent, exist_ok=True)
WD = tempfile.mkdtemp(prefix="esat_bulkumi_test_", dir=parent)
os.chdir(WD)
BC = "AAAAAA"
L = 30
CHR_LEN = 100000

# ---------------- annotation (0-based starts, right-open ends, UCSC style) --------------
genes = [  # name, strand, start, end
    ("GENEA", "+", 10000, 12000),
    ("GENEB", "-", 30000, 32000),
    ("GENEC", "+", 50000, 52000),
]
with open("genes.meta", "w") as f:
    f.write("\t".join(["name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd",
                       "exonCount", "exonStarts", "exonEnds", "name2"]) + "\n")
    for g, s, a, b in genes:
        f.write("\t".join([g + ".1", "chrT", s, str(a), str(b), str(a), str(b), "1",
                           f"{a},", f"{b},", g]) + "\n")

# ---------------- reads ----------------------------------------------------------------
reads = []   # (umi, flag, 1-based leftmost pos, cigar, NH)
def add(umi, flag, pos, cigar=f"{L}M", nh=1, n=1):
    for _ in range(n):
        reads.append((umi, flag, pos, cigar, nh))

# GENEA (+ strand; 5' end of a + read = alignment start)
add("AAAAAA", 0, 11900, n=5)               # fragment 1: 5 PCR duplicates
add("AAAAAA", 0, 11920, n=3)               # same UMI, different start -> new key only in position mode
add("CCCCCC", 0, 11910, n=2)               # 1 key in both modes
add("ACGTNA", 0, 11930, n=1)               # N-containing UMI: removed by -umiDropN
add("GGGGGG", 0, 11940, f"{L}M")           # soft-clipped duplicate shares the alignment start ->
add("GGGGGG", 0, 11940, f"5S{L-5}M")       #   one key in both modes
# GENEB (- strand; 5' end of a - read = alignment END)
add("TTTTTT", 16, 30050, f"{L}M", n=2)       # end = 30079
add("TTTTTT", 16, 30055, f"{L-5}M5S")        # different start, same end 30079 -> same key in position mode
add("TTTTTT", 16, 30070, f"{L}M")            # end = 30099 -> new key in position mode
add("ACACAC", 16, 30060, f"{L}M")
# GENEC: multimapper ignored (-multimap ignore), one unique read
add("CATCAT", 0, 51900, nh=2, n=4)
add("GTGTGT", 0, 51910)

EXPECTED = {   # (umiMode, dropN) -> {gene: collapsed keys}
    ("gene", False):     {"GENEA": 4, "GENEB": 2, "GENEC": 1},
    ("gene", True):      {"GENEA": 3, "GENEB": 2, "GENEC": 1},
    ("position", False): {"GENEA": 5, "GENEB": 3, "GENEC": 1},
    ("position", True):  {"GENEA": 4, "GENEB": 3, "GENEC": 1},
}
EXPECTED_B_POSITION_KEYS = ["ACACAC@30089", "TTTTTT@30079", "TTTTTT@30099"]

reads.sort(key=lambda r: r[2])
with open("reads.sam", "w") as f:
    f.write("@HD\tVN:1.4\tSO:coordinate\n")
    f.write(f"@SQ\tSN:chrT\tLN:{CHR_LEN}\n")
    for i, (umi, flag, pos, cig, nh) in enumerate(reads):
        f.write("\t".join([f"read{i}:{BC}:{umi}", str(flag), "chrT", str(pos), "255", cig, "*", "0", "0",
                           "A" * L, "I" * L, f"NH:i:{nh}"]) + "\n")
subprocess.run(["samtools", "view", "-b", "-o", "reads.bam", "reads.sam"], check=True)
subprocess.run(["samtools", "index", "reads.bam"], check=True)
with open("filelist.txt", "w") as f:
    f.write(f"test\t{WD}/reads.bam\n")

# ---------------- run ESAT -------------------------------------------------------------
failures = 0
for (mode, dropN), exp in EXPECTED.items():
    tag = f"{mode}{'_dropN' if dropN else ''}"
    os.makedirs(tag, exist_ok=True)
    cmd = ["java", "-Xmx2g", "-jar", JAR, "-alignments", "../filelist.txt", "-out", f"{tag}_esat.txt",
           "-geneMapping", "../genes.meta", "-task", "score3p", "-wLen", "100", "-wOlap", "50",
           "-wExt", "1000", "-sigTest", "1", "-multimap", "ignore", "-scPrep", "-umiMode", mode]
    if dropN:
        cmd.append("-umiDropN")
    with open(f"{tag}/esat.log", "w") as log:
        rc = subprocess.run(cmd, cwd=tag, stdout=log, stderr=subprocess.STDOUT).returncode
    if rc != 0:
        print(f"[FAIL] {tag}: ESAT exit code {rc}"); failures += 1; continue
    keys, reads_per_gene, gkeys = collections.Counter(), collections.Counter(), collections.defaultdict(list)
    for line in open(f"{tag}/test.umi.distributions.txt"):
        bc, g, umi, n = line.rstrip("\n").split("\t")
        keys[g] += 1; reads_per_gene[g] += int(n); gkeys[g].append(umi)
    gene_counts = {}
    with open(f"{tag}/{tag}_esat.txt.gene.txt") as fh:
        fh.readline()
        for line in fh:
            p = line.rstrip("\n").split("\t")
            gene_counts[p[0]] = float(p[-1])
    for g, e in exp.items():
        ok = keys.get(g, 0) == e and gene_counts.get(g, 0) == e
        failures += not ok
        print(f"[{'PASS' if ok else 'FAIL'}] {tag:15s} {g}: expected {e}  umi.distributions keys={keys.get(g,0)}"
              f"  gene.txt={gene_counts.get(g,0):g}  reads={reads_per_gene.get(g,0)}")
    if mode == "position":
        ok = sorted(gkeys["GENEB"]) == EXPECTED_B_POSITION_KEYS
        failures += not ok
        print(f"[{'PASS' if ok else 'FAIL'}] {tag:15s} minus-strand keys use alignment end: {sorted(gkeys['GENEB'])}")
print("ALL TESTS PASSED" if failures == 0 else f"{failures} FAILURE(S)")
sys.exit(1 if failures else 0)
