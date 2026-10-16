"""UMI-collision (saturation) correction of ESAT position-mode counts: Python version of ESAT -umiSaturation.

ESAT position mode collapses reads into keys (UMI, 5' position) per gene; different molecules that draw the same UMI
at the same position are merged, so a gene's key count k under-counts its molecules n. Two methods (ESAT
-umiSaturation original | exact; these functions reproduce ESAT's numbers to ~1e-6):

exact (recommended; Agyemang, Irizarry & Baharav, J Comput Biol 2026, with the 5' position as additional barcode
    letters): a molecule draws position p with probability q_p (the gene's reads, which collisions do not merge) and
    UMI u with probability p_u (the library's UMI frequencies), so E[k | n] = sum_p sum_u [1 - (1 - q_p p_u)^n];
    n solves E[k | n] = k (method of moments, one equation per gene).
original (the current count table's rule, scripts/qc3/build_position_tables.py, kept for reference):
    n = -K ln(1 - k/K) with all K = U_eff x P_eff keys equally likely (U_eff = 1 / sum_u p_u^2, P_eff = R^2 / sum_p
    r_p^2; k/K capped at 0.95), but reads of the gene / reads per key of genes with 5-30 keys when that correction
    exceeds 1.5x. (The current table's step 3 used STAR whole-gene reads instead, see notebook 07.)

Also: a simulator with known truth, and the real-data validation with shortened UMIs (truncation_test,
self_consistency) used in notebook 07 and tools/umi_saturation/validate_saturation.py of the ESAT repository.
The full exploration (uniform, reads, per-position variants) is archived in
analysis/07_saturation/archive_2026-10-16/.
"""
import gzip
import itertools
import numpy as np
import pandas as pd

ORIGINAL_FALLBACK = 1.5
METHOD_COLORS = {"none": "#8a8984", "current table": "#e34948", "exact": "#2a78d6"}


# ----------------------------------------------------------------------------- input
def read_keys(path: str) -> pd.DataFrame:
    """ESAT <exp>.umi.distributions.txt(.gz) in position mode: barcode, gene, UMI@pos, reads ->
    DataFrame gene, umi, pos, reads (one row per key)."""
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as fh:
        d = pd.read_csv(fh, sep="\t", header=None, names=["bc", "gene", "key", "reads"],
                        dtype={"bc": "category", "gene": "category", "key": str, "reads": np.int32})
    up = d["key"].str.split("@", n=1, expand=True)
    return pd.DataFrame({"gene": d["gene"], "umi": up[0].astype("category"), "pos": up[1].astype(np.int64),
                         "reads": d["reads"]})


def umi_frequencies(keys: pd.DataFrame, max_keys: int = 200) -> pd.Series:
    """Library UMI frequencies p_u from the keys of genes with <= max_keys keys (collisions negligible there, so common
    UMIs are not under-counted); all 4^L UMIs of the observed length, pseudo-count 0.5."""
    kpg = keys.groupby("gene", observed=True).size()
    c = keys[keys["gene"].isin(kpg.index[kpg <= max_keys])]["umi"].value_counts()
    L = int(pd.Series(c.index.astype(str)).str.len().mode()[0])
    allu = ["".join(t) for t in itertools.product("ACGT", repeat=L)]
    c = c.reindex(allu, fill_value=0).astype(float) + 0.5
    return c / c.sum()


# ----------------------------------------------------------------------------- exact estimator
def bin_probs(p: np.ndarray, nbins: int = 64):
    """Quantile bins of a probability vector: (bin mean probability, number of elements)."""
    s = np.sort(np.asarray(p, float))
    idx = np.array_split(np.arange(len(s)), min(nbins, len(s)))
    return np.array([s[i].mean() for i in idx]), np.array([len(i) for i in idx], float)


def bin_positions(q: np.ndarray, max_exact: int = 512, keep_top: int = 256, tail_bins: int = 64):
    """Position probabilities as (probabilities, multiplicities). A gene's busiest positions dominate its collisions
    and are never averaged together: up to max_exact positions as they are; beyond that the keep_top largest stay
    individual and only the tail of rare positions is binned."""
    q = np.sort(np.asarray(q, float))[::-1]
    if len(q) <= max_exact:
        return q, np.ones(len(q))
    tb, tw = bin_probs(q[keep_top:], tail_bins)
    return np.concatenate([q[:keep_top], tb]), np.concatenate([np.ones(keep_top), tw])


def solve_exact(k: float, qb: np.ndarray, qw: np.ndarray, pb: np.ndarray, wb: np.ndarray, n0: float | None = None,
                max_factor: float = 20.0, tol: float = 1e-7, max_iter: int = 100):
    """n with E[k | n] = sum_ij qw_i wb_j [1 - (1 - qb_i pb_j)^n] = k: safeguarded Newton from n0 in [k, max_factor k].
    Returns (n, hit_cap)."""
    if k <= 0:
        return 0.0, False
    L = np.log1p(-np.outer(qb, pb))
    W = np.outer(qw, wb)
    f = lambda n: float((W * (1.0 - np.exp(n * L))).sum())
    df = lambda n: float((-W * L * np.exp(n * L)).sum())
    lo, hi = float(k), max_factor * float(k)
    if f(hi) < k:
        return hi, True
    n = min(max(n0 if n0 else k, lo), hi)
    for _ in range(max_iter):
        fn = f(n) - k
        if fn < 0:
            lo = n
        else:
            hi = n
        d = df(n)
        nn = n - fn / d if d > 0 else np.inf
        if not (lo < nn < hi):
            nn = 0.5 * (lo + hi)
        if abs(nn - n) < tol * n:
            return nn, False
        n = nn
    return n, False


def estimate(keys: pd.DataFrame, p_u: pd.Series | None = None, cap: float = 0.95, exact_min_frac: float = 0.002,
             low_keys=(5, 30)) -> pd.DataFrame:
    """Per gene: k (keys), R (reads), positions, P_eff, U_eff, k/K, n_exact, n_original (+ whether original used
    reads) and the factors n/k. Below k/K = exact_min_frac exact is set to the formula (identical to < 0.1%)."""
    p_u = umi_frequencies(keys) if p_u is None else p_u
    U_eff = 1.0 / float((p_u ** 2).sum())
    pb, wb = bin_probs(p_u.to_numpy(), 64)
    gp = keys.groupby(["gene", "pos"], observed=True).agg(reads=("reads", "sum"), keys=("reads", "size")).reset_index()
    g = gp.groupby("gene", observed=True).agg(k=("keys", "sum"), R=("reads", "sum"), positions=("pos", "size"))
    ssr = gp.assign(r2=gp["reads"].astype(float) ** 2).groupby("gene", observed=True)["r2"].sum()
    g["P_eff"] = g["R"].astype(float) ** 2 / ssr
    g["U_eff"] = U_eff
    g["k_over_K"] = g["k"] / (g["P_eff"] * U_eff)
    frac = g["k_over_K"].clip(upper=cap)
    formula = pd.Series(np.where(frac < 1e-9, g["k"], g["k"] * (-np.log1p(-frac) / frac.where(frac > 0, 1.0))),
                        index=g.index)
    lo, hi = low_keys
    low = g[(g["k"] >= lo) & (g["k"] <= hi)]
    reads_based = g["R"] / (low["R"].sum() / low["k"].sum())
    g["original_used_reads"] = formula / g["k"] > ORIGINAL_FALLBACK
    g["n_original"] = np.where(g["original_used_reads"], reads_based, formula)
    g["n_exact"] = formula
    g["exact_capped"] = False
    todo = g.index[g["k_over_K"] >= exact_min_frac]
    for gene, d in gp[gp["gene"].isin(todo)].groupby("gene", observed=True):
        q = d["reads"].to_numpy(float)
        qb, qw = bin_positions(q / q.sum())
        n, capped = solve_exact(float(g.at[gene, "k"]), qb, qw, pb, wb, n0=float(formula[gene]))
        g.at[gene, "n_exact"], g.at[gene, "exact_capped"] = n, capped
    g["factor_exact"], g["factor_original"] = g["n_exact"] / g["k"], g["n_original"] / g["k"]
    g.attrs = {"U_eff": U_eff}
    return g


# ----------------------------------------------------------------------------- validation
def nucleotide_umi_probs(L: int = 6, comp: dict | None = None) -> np.ndarray:
    """UMI probabilities from independent bases (the paper's model; default its composition A .23 C .24 G .21 T .32)."""
    comp = comp or {"A": 0.23, "C": 0.24, "G": 0.21, "T": 0.32}
    p = np.array([np.prod([comp[b] for b in t]) for t in itertools.product("ACGT", repeat=L)])
    return p / p.sum()


def simulate_grid(p_u: np.ndarray, peff_targets=(2, 5, 20, 80), n_grid=None, reps: int = 3, seed: int = 0,
                  reads_per_mol_extra: float = 0.6, cap: float = 0.95) -> pd.DataFrame:
    """Simulation with known truth: a Dirichlet position profile with about the target effective position number;
    n molecules draw positions and UMIs; each molecule gets 1 + Poisson(extra) reads. The methods see only the
    observed keys and reads (position profile and reads per molecule estimated from them). One row per simulation:
    n (truth), k, k/K, n_original, n_exact."""
    rng = np.random.default_rng(seed)
    n_grid = n_grid if n_grid is not None else np.unique(np.logspace(2, 5.3, 16).astype(int))
    U_eff = 1.0 / (p_u ** 2).sum()
    pb, wb = bin_probs(p_u, 64)
    rows = []
    for target in peff_targets:
        npos = max(int(target * 3), 3)
        for rep in range(reps):
            for _ in range(50):
                q = rng.dirichlet(np.ones(npos) * 0.7)
                if abs(1 / (q ** 2).sum() / target - 1) < 0.25:
                    break
            for n in n_grid:
                pos = rng.choice(npos, size=n, p=q)
                k = np.unique(pos.astype(np.int64) * len(p_u) + rng.choice(len(p_u), size=n, p=p_u)).size
                r = 1 + rng.poisson(reads_per_mol_extra, size=n)
                rp = np.bincount(pos, weights=r, minlength=npos)
                rp = rp[rp > 0]
                P_eff = rp.sum() ** 2 / (rp ** 2).sum()
                K = P_eff * U_eff
                fr = min(k / K, cap)
                formula = k * (-np.log1p(-fr) / fr)
                n_orig = rp.sum() / (1 + reads_per_mol_extra) if formula / k > ORIGINAL_FALLBACK else formula
                qb, qw = bin_positions(rp / rp.sum())
                n_e, _ = solve_exact(k, qb, qw, pb, wb, n0=formula)
                rows.append(dict(P_eff_target=target, rep=rep, n=n, k=k, k_over_K=k / K, n_original=n_orig,
                                 n_exact=n_e))
    return pd.DataFrame(rows)


def truncate_umis(keys: pd.DataFrame, L: int) -> pd.DataFrame:
    """Synthetically shortened UMIs (the paper's validation): keep the first L bases, re-collapse (gene, UMI, pos)."""
    t = keys.assign(umi=keys["umi"].astype(str).str[:L])
    return t.groupby(["gene", "umi", "pos"], observed=True, as_index=False)["reads"].sum().astype({"umi": "category"})


def truncation_test(keys: pd.DataFrame, lengths=(5, 4, 3), truth_max_frac: float = 0.2) -> pd.DataFrame:
    """Genes nearly unsaturated with the full UMIs (k/K < truth_max_frac): truth = their full-length exact estimate.
    With UMIs shortened to each length: keys, k/K, original and exact estimates."""
    full = estimate(keys)
    ok = full.index[(full["k_over_K"] < truth_max_frac) & (full["k"] >= 20)]
    out = []
    for L in lengths:
        e = estimate(truncate_umis(keys, L))
        g = ok.intersection(e.index)
        out.append(pd.DataFrame({"L": L, "truth": full.loc[g, "n_exact"], "k": e.loc[g, "k"],
                                 "k_over_K": e.loc[g, "k_over_K"], "original": e.loc[g, "n_original"],
                                 "exact": e.loc[g, "n_exact"]}, index=g))
    return pd.concat(out)


SELF_METHODS = {"none": "k", "current table": "n_original", "exact": "n_exact"}


def self_consistency(keys: pd.DataFrame, lengths=(5, 4, 3), min_keys: int = 20) -> pd.DataFrame:
    """Robustness to saturation: each method with the full UMIs and with UMIs shortened to each length (collisions
    4x, 16x, 64x more likely; same molecules). Per gene x length: <method>_full, <method>, k/K after shortening."""
    full = estimate(keys)
    full = full[full["k"] >= min_keys]
    out = []
    for L in lengths:
        e = estimate(truncate_umis(keys, L)).reindex(full.index)
        d = pd.DataFrame({"L": L, "k_over_K": e["k_over_K"], "k_over_K_full": full["k_over_K"]}, index=full.index)
        for lab, c in SELF_METHODS.items():
            d[f"{lab}_full"], d[lab] = full[c].astype(float), e[c].astype(float)
        out.append(d)
    return pd.concat(out)


def self_consistency_summary(sc: pd.DataFrame, sat_min: float = 0.2) -> pd.DataFrame:
    """Per UMI length and method: R^2 of log2(shortened) vs log2(full) around the identity line, for all genes and for
    genes saturated after shortening (k/K > sat_min); bias (median log2 shortened / full) and spread (robust SD)."""
    rows = []
    for L, d in sc.groupby("L"):
        sat = d["k_over_K"] > sat_min
        for lab in SELF_METHODS:
            x, y = np.log2(d[f"{lab}_full"]), np.log2(d[lab])
            ok = np.isfinite(x) & np.isfinite(y)
            r2 = lambda m: 1 - ((y[m] - x[m]) ** 2).sum() / ((y[m] - y[m].mean()) ** 2).sum()
            lr = (y - x)[ok & sat]
            rows.append({"UMI length": L, "UMIs": 4 ** L, "method": lab, "genes": int(ok.sum()),
                         "R² all genes": round(float(r2(ok)), 4), f"saturated genes (k/K > {sat_min})":
                         int((ok & sat).sum()), "R² saturated": round(float(r2(ok & sat)), 4),
                         "bias saturated (median log2)": round(float(lr.median()), 3),
                         "spread saturated (robust SD log2)": round(float(1.4826 * (lr - lr.median()).abs().median()), 3)})
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------- all libraries, plotting
def load_exact_vs_current(min_count: float = 10.0):
    """All 509 libraries: the current table (counts_pos_final), the exact table (ESAT -umiSaturation exact,
    data/position_mode/exact_509/) and the uncorrected position counts, aligned; plus a long table of gene x library
    entries with >= min_count in either table (gene, library, current, exact, position, log2 (exact+1)/(current+1))."""
    from rqtl_utils import ROOT
    d = ROOT / "data" / "position_mode"
    cur = pd.read_csv(d / "counts_pos_final.tsv.gz", sep="\t", index_col=0)
    ex = pd.read_csv(d / "exact_509" / "counts_exact.tsv.gz", sep="\t", index_col=0).loc[cur.index, cur.columns]
    pos = pd.read_csv(d / "position_mode_gene_counts.tsv.gz", sep="\t", index_col=0).loc[cur.index, cur.columns]
    keep = ((cur >= min_count) | (ex >= min_count)).to_numpy()
    gi, li = np.nonzero(keep)
    c, e, p = cur.to_numpy()[gi, li], ex.to_numpy()[gi, li], pos.to_numpy()[gi, li]
    long = pd.DataFrame({"gene": cur.index[gi], "library": cur.columns[li], "current": c, "exact": e, "position": p,
                         "log2_ratio": np.log2((e + 1) / (c + 1))})
    return cur, ex, pos, long


def summary_box(x_labels, values, order, name: str, color: str, showlegend: bool = True, legendgroup=None):
    """A plotly Box from precomputed statistics (quartiles, whiskers at the 5th / 95th percentiles), so the page does
    not carry every point."""
    import plotly.graph_objects as go
    d = pd.DataFrame({"x": np.asarray(x_labels), "v": np.asarray(values, float)}).dropna()
    q = d.groupby("x")["v"].quantile([0.05, 0.25, 0.5, 0.75, 0.95]).unstack()
    cats = [c for c in order if c in q.index]
    q = q.loc[cats]
    return go.Box(x=cats, q1=q[0.25].round(4), median=q[0.5].round(4), q3=q[0.75].round(4),
                  lowerfence=q[0.05].round(4), upperfence=q[0.95].round(4), name=name, marker_color=color,
                  showlegend=showlegend, legendgroup=legendgroup or name, boxpoints=False)
