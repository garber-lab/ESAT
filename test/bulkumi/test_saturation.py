#!/usr/bin/env python3
"""End-to-end test for ESAT -umiSaturation (UMI-collision correction) with -umiMode position.

Builds a synthetic genome and BAM with 4-nt UMIs drawn from a skewed nucleotide composition (256 UMIs, effective
number ~140), so collisions are frequent at small read numbers (SKEW: 1200 molecules on 120 positions, one of
them carrying 60% of the molecules):
  BG01..BG30  150 molecules each at 150 distinct positions (no collisions; used to estimate UMI frequencies)
  MID         200 molecules on 3 positions (moderate saturation)
  SAT         1500 molecules on 3 positions (strong saturation: ~2x more molecules than keys)
Checks
  (1) <exp>.umi.saturation.txt lists every gene, background factors ~1;
  (2) the exact estimator recovers the simulated molecule numbers of MID and SAT within 10%;
  (3) <out>.gene.corrected.txt = gene.txt x the chosen estimator's factor, for every mode; factors >= 1;
  (4) -umiSaturation none (default) writes no corrected file; an illegal mode fails.

Usage: python test_saturation.py <esat.jar> [parent_dir_for_a_new_workdir]
"""
import csv, os, random, subprocess, sys, tempfile

JAR = os.path.abspath(sys.argv[1])
parent = sys.argv[2] if len(sys.argv) > 2 and sys.argv[2] else None
if parent:
    os.makedirs(parent, exist_ok=True)
WD = tempfile.mkdtemp(prefix="esat_saturation_test_", dir=parent)   # never deletes anything
os.chdir(WD)
rng = random.Random(7)
BC, L, GLEN, GAP = "AAAAAA", 30, 2000, 3000
NT, NT_P = "ACGT", [0.35, 0.15, 0.15, 0.35]

genes = [f"BG{i:02d}" for i in range(1, 31)] + ["MID", "SAT", "SKEW"]
TRUE_N = {"MID": 200, "SAT": 1500, "SKEW": 1200}
start = {g: 10000 + i * (GLEN + GAP) for i, g in enumerate(genes)}
CHR_LEN = 10000 + len(genes) * (GLEN + GAP) + 10000

with open("genes.meta", "w") as f:
    f.write("\t".join(["name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount",
                       "exonStarts", "exonEnds", "name2"]) + "\n")
    for g in genes:
        a, b = start[g], start[g] + GLEN
        f.write("\t".join([g + ".1", "chrT", "+", str(a), str(b), str(a), str(b), "1", f"{a},", f"{b},", g]) + "\n")

def umi():
    return "".join(rng.choices(NT, NT_P, k=4))

reads = []   # (umi, 1-based pos)
for g in genes[:30]:
    for p in rng.sample(range(start[g] + 500, start[g] + 1900), 150):
        reads += [(umi(), p)] * rng.randint(1, 3)          # 1-3 PCR duplicates per molecule
for g in ["MID", "SAT"]:
    pos = [start[g] + 1700, start[g] + 1750, start[g] + 1800]
    for _ in range(TRUE_N[g]):
        reads += [(umi(), rng.choices(pos, [0.6, 0.3, 0.1])[0])] * rng.randint(1, 3)
# SKEW: 120 positions (more than the 64 bins used for the UMI probabilities), two of them dominant; the busiest
# positions must not be averaged together (a quantile-binning bug under-corrected such genes)
skew_pos = [start["SKEW"] + 500 + 10 * i for i in range(120)]
skew_w = [0.6, 0.05] + [0.35 / 118] * 118
for _ in range(TRUE_N["SKEW"]):
    reads += [(umi(), rng.choices(skew_pos, skew_w)[0])] * rng.randint(1, 3)
reads.sort(key=lambda r: r[1])
with open("reads.sam", "w") as f:
    f.write("@HD\tVN:1.4\tSO:coordinate\n")
    f.write(f"@SQ\tSN:chrT\tLN:{CHR_LEN}\n")
    for i, (u, pos) in enumerate(reads):
        f.write("\t".join([f"read{i}:{BC}:{u}", "0", "chrT", str(pos), "255", f"{L}M", "*", "0", "0",
                           "A" * L, "I" * L, "NH:i:1"]) + "\n")
subprocess.run(["samtools", "view", "-b", "-o", "reads.bam", "reads.sam"], check=True)
subprocess.run(["samtools", "index", "reads.bam"], check=True)
with open("filelist.txt", "w") as f:
    f.write(f"test\t{WD}/reads.bam\n")


def run(tag, mode):
    os.makedirs(tag, exist_ok=True)
    cmd = ["java", "-Xmx2g", "-jar", JAR, "-alignments", "../filelist.txt", "-out", f"{tag}_esat.txt",
           "-geneMapping", "../genes.meta", "-task", "score3p", "-wLen", "100", "-wOlap", "50", "-wExt", "1000",
           "-sigTest", "1", "-multimap", "ignore", "-scPrep", "-umiMode", "position"]
    if mode is not None:
        cmd += ["-umiSaturation", mode]
    with open(f"{tag}/esat.log", "w") as log:
        return subprocess.run(cmd, cwd=tag, stdout=log, stderr=subprocess.STDOUT).returncode


def gene_file(path):
    with open(path) as fh:
        fh.readline()
        return {p[0]: float(p[-1]) for p in (l.rstrip("\n").split("\t") for l in fh)}


failures = 0
def check(ok, msg):
    global failures
    failures += not ok
    print(f"[{'PASS' if ok else 'FAIL'}] {msg}")

for mode in ["exact", "uniform", "reads", "original"]:
    rc = run(mode, mode)
    check(rc == 0, f"{mode}: ESAT exit code {rc}")
    if rc:
        continue
    sat = {r["gene"]: r for r in csv.DictReader(open(f"{mode}/test.umi.saturation.txt"), delimiter="\t")}
    check(set(sat) == set(genes), f"{mode}: saturation file lists all {len(genes)} genes ({len(sat)})")
    g, c = gene_file(f"{mode}/{mode}_esat.txt.gene.txt"), gene_file(f"{mode}/{mode}_esat.txt.gene.corrected.txt")
    fac = {k: float(v[f"factor_{mode}"]) for k, v in sat.items()}
    bad = [k for k in genes if abs(c[k] - g[k] * fac[k]) > 1e-3 * max(1.0, c[k])]
    check(not bad, f"{mode}: gene.corrected.txt = gene.txt x factor_{mode} for all genes {bad[:3]}")
    low = min(min(float(v["factor_uniform"]), float(v["factor_exact"])) for v in sat.values())
    check(low >= 0.999, f"{mode}: uniform and exact factors >= 1 for all genes (min {low:.3f})")
    if mode == "exact":
        bg = max(float(sat[k]["factor_exact"]) for k in genes[:30])
        check(bg < 1.05, f"exact: background genes factor <= {bg:.3f} (< 1.05)")
        for k, n in TRUE_N.items():
            est, keys = float(sat[k]["n_exact"]), int(sat[k]["keys"])
            check(abs(est / n - 1) < 0.10, f"exact: {k} keys={keys} true n={n} estimate={est:.0f} "
                                            f"(uniform {float(sat[k]['n_uniform']):.0f})")

rc = run("none", None)
check(rc == 0 and not os.path.exists("none/none_esat.txt.gene.corrected.txt")
      and not os.path.exists("none/test.umi.saturation.txt"), "default (none): no corrected or saturation file")
check(run("illegal", "bogus") != 0, "illegal -umiSaturation value fails")
print("ALL TESTS PASSED" if failures == 0 else f"{failures} FAILURE(S)")
sys.exit(1 if failures else 0)
