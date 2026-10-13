#!/usr/bin/env python3
"""Deterministic synthetic CEL-seq2-like fixture for ESAT regression testing.

Writes into <outdir>:
  genes.meta      gene-to-transcript mapping (UCSC genePred columns + name2)
  reads.sam       coordinate-sorted alignments, read names <id>:<barcode>:<umi>
The fixture exercises: both strands, multi-exon genes, multi-isoform genes (meta-transcript
collapsing), spliced reads, 3' extension (-wExt) incl. truncation at a neighbouring gene,
antisense overlap, a second chromosome, multimappers (NH>1), soft clips, PCR duplicates,
UMI collisions in a highly expressed gene, N-containing UMIs and intergenic reads.
The output is fully determined by SEED; do not change it without regenerating the goldens.

Usage: python make_fixture.py <outdir>
"""
import os, random, sys

SEED = 20261013
L = 50
BC = "TCACAG"
CHR = {"chr1": 200000, "chr2": 100000}
rng = random.Random(SEED)

# name, chrom, strand, isoforms (each a list of (start0, end) exons), molecules
GENES = [
    ("PLUS1",  "chr1", "+", [[(10000, 10400), (11000, 11600)]], 40),
    ("MINUS1", "chr1", "-", [[(20000, 20600), (21500, 21800)]], 60),
    ("ISO2",   "chr1", "+", [[(30000, 30300), (31000, 31500)],
                              [(30000, 30300), (30800, 30900), (31200, 31700)]], 80),
    ("HIGH",   "chr1", "+", [[(40000, 40500), (41000, 42000)]], 6000),   # saturates the 4,096 UMIs
    ("NEIGHA", "chr1", "+", [[(50000, 51000)]], 30),                      # extension truncated by NEIGHB
    ("NEIGHB", "chr1", "+", [[(51300, 52300)]], 30),
    ("SENSE",  "chr1", "+", [[(60000, 61500)]], 50),                      # antisense pair
    ("ANTI",   "chr1", "-", [[(61000, 62500)]], 50),
    ("LOW1",   "chr1", "-", [[(70000, 70800)]], 3),
    ("CHR2A",  "chr2", "+", [[(5000, 5400), (6000, 6800)]], 120),
    ("CHR2B",  "chr2", "-", [[(20000, 21200)]], 200),
]

def meta_exons(isoforms):
    iv = sorted(e for iso in isoforms for e in iso)
    out = []
    for s, e in iv:
        if out and s <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return out

def rand_umi():
    u = "".join(rng.choice("ACGT") for _ in range(6))
    if rng.random() < 0.02:            # ~2% N-containing UMIs
        i = rng.randrange(6); u = u[:i] + "N" + u[i + 1:]
    return u

def rand_seq(n):
    return "".join(rng.choice("ACGT") for _ in range(n))

def read_from_5p(exons, strand, five_p):
    """(pos1, cigar, aligned length) for a read of length L whose 5' end is at genomic 0-based
    coordinate five_p, following the exon structure (spliced when it crosses a junction)."""
    bases = []
    if strand == "+":
        for s, e in exons:
            for x in range(max(s, five_p), e):
                bases.append(x)
                if len(bases) >= L: break
            if len(bases) >= L: break
    else:
        for s, e in reversed(exons):
            for x in range(min(e - 1, five_p), s - 1, -1):
                bases.append(x)
                if len(bases) >= L: break
            if len(bases) >= L: break
        bases = bases[::-1]
    if len(bases) < 20:
        return None
    cig, run, prev = [], 1, bases[0]
    for b in bases[1:]:
        if b == prev + 1:
            run += 1
        else:
            cig.append(f"{run}M"); cig.append(f"{b - prev - 1}N"); run = 1
        prev = b
    cig.append(f"{run}M")
    return bases[0] + 1, "".join(cig), len(bases)

records = []
def emit(chrom, flag, pos1, cigar, seqlen, umi, nh=1, clip=0):
    if clip:
        cigar = (f"{clip}S" + cigar) if flag == 0 else (cigar + f"{clip}S")
        seqlen += clip
    records.append((chrom, pos1, flag, cigar, seqlen, umi, nh))

for name, chrom, strand, isos, nmol in GENES:
    ex = meta_exons(isos)
    tx = [x for s, e in ex for x in range(s, e)]          # transcript bases 5'->3'
    if strand == "-":
        tx = tx[::-1]
    for m in range(nmol):
        umi = rand_umi()
        nfrag = 1 + (rng.random() < 0.35) + (rng.random() < 0.10)   # IVT copies fragmented differently
        for _ in range(nfrag):
            if rng.random() < 0.06:
                # fragment whose 5' end lies in the 3' extension (past txEnd for +, before txStart for -)
                off = rng.randrange(10, 900)
                if strand == "+":
                    five_p = ex[-1][1] + off; r = read_from_5p([(five_p, five_p + L)], strand, five_p)
                else:
                    five_p = ex[0][0] - off; r = read_from_5p([(five_p - L + 1, five_p + 1)], strand, five_p)
            else:
                # 5' end of the read ~Exp(mean 180) upstream of the 3' end, within the last 600 nt
                d = int(rng.expovariate(1 / 180.0)) % min(len(tx), 600)
                five_p = tx[max(0, len(tx) - 1 - d - L)]
                r = read_from_5p(ex, strand, five_p)
            if r is None:
                continue
            pos1, cig, sl = r
            ndup = 1 + int(rng.expovariate(1 / 3.0))                      # PCR duplicates
            flag = 0 if strand == "+" else 16
            nh = 2 if rng.random() < 0.04 else 1                          # multimappers
            for _d in range(ndup):
                clip = 3 if rng.random() < 0.1 else 0                     # occasional soft clip
                emit(chrom, flag, pos1, cig, sl, umi, nh, clip)
for i in range(40):                                                       # intergenic reads
    chrom = rng.choice(list(CHR)); pos1 = rng.randrange(1000, CHR[chrom] - 1000)
    emit(chrom, rng.choice([0, 16]), pos1, f"{L}M", L, rand_umi())

order = {c: i for i, c in enumerate(CHR)}
records.sort(key=lambda r: (order[r[0]], r[1]))

out = sys.argv[1]
os.makedirs(out, exist_ok=True)
with open(os.path.join(out, "genes.meta"), "w") as f:
    f.write("\t".join(["name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd",
                       "exonCount", "exonStarts", "exonEnds", "name2"]) + "\n")
    for name, chrom, strand, isos, _ in GENES:
        for j, iso in enumerate(isos):
            a, b = iso[0][0], iso[-1][1]
            f.write("\t".join([f"{name}.{j+1}", chrom, strand, str(a), str(b), str(a), str(b), str(len(iso)),
                               ",".join(str(s) for s, _ in iso) + ",", ",".join(str(e) for _, e in iso) + ",",
                               name]) + "\n")
with open(os.path.join(out, "reads.sam"), "w") as f:
    f.write("@HD\tVN:1.4\tSO:coordinate\n")
    for c, n in CHR.items():
        f.write(f"@SQ\tSN:{c}\tLN:{n}\n")
    for i, (chrom, pos1, flag, cig, sl, umi, nh) in enumerate(records):
        f.write("\t".join([f"r{i}:{BC}:{umi}", str(flag), chrom, str(pos1), "255", cig, "*", "0", "0",
                           rand_seq(sl), "I" * sl, f"NH:i:{nh}"]) + "\n")
print(f"{len(records)} alignments, {len(GENES)} genes")
