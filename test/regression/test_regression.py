#!/usr/bin/env python3
"""Golden regression test: the default ESAT code path must reproduce the production jar.

Golden outputs in test/regression/golden/<config>/ were produced ONCE by the production
binary esat.v0.1_09.09.16_24.18.umihack.jar (UMass HPC, /share/data/umw_biocore/dnext_data/
singleCell/singleCellScripts/) on the fixture from make_fixture.py; see golden/PROVENANCE.txt.

  python test_regression.py <esat.jar>                 # compare against goldens
  python test_regression.py <esat.jar> --make-golden   # (re)create goldens; use ONLY with the production jar
  python test_regression.py <esat.jar> --expect-fail -umiMode position
                                                       # sensitivity check: extra args must change outputs

Each config is run with the default UMI mode and explicitly with -umiMode gene; both must match
the goldens exactly (gene.txt, window.txt, sorted umi.distributions).
Work happens in a fresh temporary directory created by this script (and removed by it on success).
"""
import os, sys, shutil, subprocess, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
GOLD = os.path.join(HERE, "golden")
CONFIGS = {
    # the parameters used by the CEL-seq2 DolphinNext pipeline
    "production": "-task score3p -wLen 100 -wOlap 50 -wExt 1000 -sigTest .01 -multimap ignore -scPrep",
    # alternative windowing / thresholds to cover more of the counting code
    "alt": "-task score3p -wLen 50 -wOlap 25 -wExt 500 -sigTest .05 -multimap ignore -scPrep -umiMin 2",
}

def run(jar, fix, workdir, params, extra):
    os.makedirs(workdir, exist_ok=True)
    with open(os.path.join(workdir, "filelist.txt"), "w") as f:
        f.write(f"lib1\t{fix}/reads.bam\n")
    cmd = ["java", "-Xmx2g", f"-Djava.io.tmpdir={workdir}", "-jar", jar, "-alignments", "filelist.txt",
           "-out", "lib1_esat.txt", "-geneMapping", f"{fix}/genes.meta"] + params.split() + extra
    with open(os.path.join(workdir, "esat.log"), "w") as log:
        rc = subprocess.run(cmd, cwd=workdir, stdout=log, stderr=subprocess.STDOUT).returncode
    outs = {}
    for name, fn, sort in [("gene.txt", "lib1_esat.txt.gene.txt", False),
                           ("window.txt", "lib1_esat.txt.window.txt", False),
                           ("umi.distributions.sorted.txt", "lib1.umi.distributions.txt", True)]:
        p = os.path.join(workdir, fn)
        if not os.path.exists(p):
            outs[name] = None; continue
        lines = open(p).read().splitlines()
        outs[name] = "\n".join(sorted(lines) if sort else lines) + "\n"
    with open(os.path.join(workdir, "esat.log"), errors="replace") as fh:
        outs["log_lines"] = sum(1 for _ in fh)
    return rc, outs

def main():
    args = sys.argv[1:]
    jar = os.path.abspath(args.pop(0))
    make = "--make-golden" in args
    expect_fail = "--expect-fail" in args
    user_extra = [a for a in args if a not in ("--make-golden", "--expect-fail")]
    tmp = tempfile.mkdtemp(prefix="esat_regress_")
    fix = os.path.join(tmp, "fixture")
    subprocess.run([sys.executable, os.path.join(HERE, "make_fixture.py"), fix], check=True)
    subprocess.run(["samtools", "view", "-b", "-o", f"{fix}/reads.bam", f"{fix}/reads.sam"], check=True)
    subprocess.run(["samtools", "index", f"{fix}/reads.bam"], check=True)
    if make or expect_fail:
        modes = {"custom": user_extra}
    else:
        modes = {"default": [], "explicit_gene": ["-umiMode", "gene"]}
    failures, mismatches = 0, 0
    for cfg, params in CONFIGS.items():
        for mode, extra in modes.items():
            rc, outs = run(jar, fix, os.path.join(tmp, cfg, mode), params, extra)
            if rc != 0:
                print(f"[FAIL] {cfg}/{mode}: ESAT exit {rc} (log {tmp}/{cfg}/{mode}/esat.log)"); failures += 1; continue
            gdir = os.path.join(GOLD, cfg)
            if make:
                os.makedirs(gdir, exist_ok=True)
                for name, txt in outs.items():
                    if name == "log_lines":
                        open(os.path.join(gdir, "log_lines.txt"), "w").write(f"{txt}\n"); continue
                    open(os.path.join(gdir, name), "w").write(txt or "")
                print(f"[GOLDEN] {cfg}: " + ", ".join(f"{k} {v if k == 'log_lines' else (v or '').count(chr(10))} lines" for k, v in outs.items()))
                continue
            # log volume guard: per-read/per-base logging slows ESAT by >10x on real data
            gl = int(open(os.path.join(gdir, "log_lines.txt")).read())
            okl = outs["log_lines"] <= 2 * gl + 50
            if not expect_fail:
                failures += not okl
            print(f"[{'MATCH' if okl else 'FAIL '}] {cfg:10s} {mode:13s} log volume {outs['log_lines']} lines (production {gl})")
            for name, txt in outs.items():
                if name == "log_lines":
                    continue
                gold = open(os.path.join(gdir, name)).read()
                ok = txt == gold
                mismatches += not ok
                msg = ""
                if not ok:
                    g, t = gold.splitlines(), (txt or "").splitlines()
                    diff = [(a, b) for a, b in zip(g, t) if a != b][:2]
                    msg = f"  ({len(g)} golden vs {len(t)} lines; first differences: {diff})"
                if not expect_fail:
                    failures += not ok
                print(f"[{'MATCH' if ok else 'DIFF '}] {cfg:10s} {mode:13s} {name}{msg if not expect_fail else ''}")
    if expect_fail:
        ok = mismatches > 0 and failures == 0
        print("SENSITIVITY CHECK PASSED (outputs changed)" if ok else "SENSITIVITY CHECK FAILED (outputs unchanged)")
        failures = 0 if ok else 1
    elif not make:
        print("ALL REGRESSION TESTS PASSED" if failures == 0 else f"{failures} REGRESSION FAILURE(S)")
    if failures == 0:
        shutil.rmtree(tmp, ignore_errors=True)   # only the directory this script created
    sys.exit(1 if failures else 0)

if __name__ == "__main__":
    main()
