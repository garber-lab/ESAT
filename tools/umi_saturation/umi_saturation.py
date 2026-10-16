"""UMI-collision (saturation) corrections for ESAT position-mode counts: prototype of the ESAT -umiSaturation flag.

ESAT position mode collapses reads into keys (UMI, 5' position) per gene. Different molecules that draw the same UMI
and start at the same position are merged, so a gene's key count k under-counts its fragments n. Two estimators of n:

uniform (our original handling, scripts/qc3/build_position_tables.py): all K = U_eff x P_eff keys equally likely,
    U_eff = 1 / sum_u p_u^2 (effective number of UMIs of the library), P_eff = R^2 / sum_p r_p^2 (effective number of
    5' positions of the gene, from reads r_p, which collisions do not saturate). E[k | n] = K (1 - (1 - 1/K)^n)
    ~ K (1 - exp(-n / K)), so n = -K ln(1 - k/K); k/K capped at 0.95 (factor 3.15).
exact (Agyemang, Irizarry & Baharav, J Comput Biol 2026, extended to positions): molecules draw a position p with
    probability q_p (the gene's read distribution) and a UMI u with probability p_u (the library's UMI frequencies),
    independently, so E[k | n] = sum_p sum_u [1 - (1 - q_p p_u)^n]; n = the solution of E[k | n] = k (method of
    moments; f is increasing in n). The paper uses UMIs only (single cell, gene x cell) and a fixed nucleotide
    composition; here p_u is measured from the library and positions are added.
reads (fallback of the original handling): reads of the gene / reads per key of low-count genes (5-30 keys), on the
    key scale; unaffected by collisions but noisier (PCR / IVT amplification varies by molecule).
"""
import gzip
import numpy as np
import pandas as pd


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
    """Library UMI frequencies p_u from the keys of genes with <= max_keys keys (where collisions are negligible, so
    common UMIs are not under-counted); all 4^L UMIs of the observed length, unseen ones get a pseudo-count 0.5."""
    kpg = keys.groupby("gene", observed=True).size()
    low = keys[keys["gene"].isin(kpg.index[kpg <= max_keys])]
    c = low["umi"].value_counts()
    L = int(pd.Series(c.index.astype(str)).str.len().mode()[0])
    import itertools
    allu = ["".join(t) for t in itertools.product("ACGT", repeat=L)]
    c = c.reindex(allu, fill_value=0).astype(float) + 0.5
    return c / c.sum()


def bin_probs(p: np.ndarray, nbins: int = 128):
    """Compress the UMI distribution into nbins quantile bins: (bin mean probability, number of UMIs in the bin)."""
    s = np.sort(np.asarray(p, float))
    idx = np.array_split(np.arange(len(s)), nbins)
    return np.array([s[i].mean() for i in idx]), np.array([len(i) for i in idx], float)


def bin_positions(q: np.ndarray, max_exact: int = 512, keep_top: int = 256, tail_bins: int = 64):
    """Position probabilities for the exact model, as (probabilities, multiplicities). A gene's few busiest positions
    dominate its collisions and must never be averaged together: up to max_exact positions are used as they are;
    beyond that the keep_top largest stay individual and only the tail of rare positions is put into tail_bins
    quantile bins. (Plain quantile bins over all positions averaged busy positions in pairs, e.g. RPL26 with 103
    positions, and under-corrected by up to 8%.)"""
    q = np.sort(np.asarray(q, float))[::-1]
    if len(q) <= max_exact:
        return q, np.ones(len(q))
    tb, tw = bin_probs(q[keep_top:], tail_bins)
    return np.concatenate([q[:keep_top], tb]), np.concatenate([np.ones(keep_top), tw])


def gene_positions(keys: pd.DataFrame) -> pd.DataFrame:
    """Per gene x 5' position: reads and keys."""
    return keys.groupby(["gene", "pos"], observed=True).agg(reads=("reads", "sum"), keys=("reads", "size")).reset_index()


def forward_exact(n: float, q: np.ndarray, pb: np.ndarray, wb: np.ndarray) -> float:
    """E[k | n] = sum_p sum_b w_b [1 - (1 - q_p p_b)^n]."""
    L = np.log1p(-np.outer(q, pb))
    return float((wb * (1.0 - np.exp(n * L))).sum())


def solve_exact(k: float, qb: np.ndarray, qw: np.ndarray, pb: np.ndarray, wb: np.ndarray, n0: float | None = None,
                max_factor: float = 20.0, tol: float = 1e-7, max_iter: int = 100):
    """Fast version used by ESAT: positions and UMIs both binned (qb, qw = bin probability and number of positions;
    pb, wb = same for UMIs); safeguarded Newton on f(n) - k from n0 (default: k), bracket [k, max_factor k].
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
        step = fn / d if d > 0 else np.inf
        nn = n - step
        if not (lo < nn < hi):
            nn = 0.5 * (lo + hi)
        if abs(nn - n) < tol * n:
            return nn, False
        n = nn
    return n, False


def invert_exact(k: float, q: np.ndarray, pb: np.ndarray, wb: np.ndarray, max_factor: float = 20.0, tol=1e-6):
    """n with E[k | n] = k (bisection on n in [k, max_factor k]); returns (n, hit_cap)."""
    if k <= 0:
        return 0.0, False
    L = np.log1p(-np.outer(q, pb))
    f = lambda n: float((wb * (1.0 - np.exp(n * L))).sum())
    lo, hi = float(k), float(k) * 2
    while f(hi) < k:
        lo, hi = hi, hi * 2
        if hi > max_factor * k:
            return max_factor * k, True
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if f(mid) < k:
            lo = mid
        else:
            hi = mid
        if hi - lo < tol * hi:
            break
    return 0.5 * (lo + hi), False


def estimate(keys: pd.DataFrame, p_u: pd.Series | None = None, nbins: int = 64, cap: float = 0.95,
             exact_min_frac: float = 0.002, low_keys=(5, 30)) -> pd.DataFrame:
    """Per gene: k (keys), R (reads), positions, P_eff, U_eff, K, k/K, n_uniform, n_exact, n_reads and the factors
    n/k. exact is computed for genes with k/K >= exact_min_frac (below it, exact = uniform to < 0.1%), UMIs in nbins
    quantile bins, positions via bin_positions (solve_exact, the algorithm ESAT uses)."""
    p_u = umi_frequencies(keys) if p_u is None else p_u
    U_eff = 1.0 / float((p_u ** 2).sum())
    pb, wb = bin_probs(p_u.to_numpy(), nbins)
    gp = gene_positions(keys)
    g = gp.groupby("gene", observed=True).agg(k=("keys", "sum"), R=("reads", "sum"), positions=("pos", "size"))
    g["ssr"] = gp.assign(r2=gp["reads"].astype(float) ** 2).groupby("gene", observed=True)["r2"].sum()
    g["P_eff"] = g["R"].astype(float) ** 2 / g["ssr"]
    g["U_eff"] = U_eff
    g["K"] = g["P_eff"] * U_eff
    g["k_over_K"] = g["k"] / g["K"]
    frac = g["k_over_K"].clip(upper=cap)
    # factor -ln(1 - f)/f applied to k (as in build_position_tables.py): = -K ln(1 - k/K) below the cap, 3.15 k above
    g["n_uniform"] = np.where(frac < 1e-9, g["k"], g["k"] * (-np.log1p(-frac) / frac.where(frac > 0, 1.0)))
    g["uniform_capped"] = g["k_over_K"] > cap
    g["n_exact"] = g["n_uniform"]
    g["exact_capped"] = False
    todo = g.index[g["k_over_K"] >= exact_min_frac]
    grp = gp[gp["gene"].isin(todo)].groupby("gene", observed=True)
    pb, wb = bin_probs(p_u.to_numpy(), nbins)
    for gene, d in grp:
        q = d["reads"].to_numpy(float)
        q /= q.sum()
        qb, qw = bin_positions(q)
        n, capped = solve_exact(float(g.at[gene, "k"]), qb, qw, pb, wb, n0=float(g.at[gene, "n_uniform"]))
        g.at[gene, "n_exact"], g.at[gene, "exact_capped"] = n, capped
    lo, hi = low_keys
    low = g[(g["k"] >= lo) & (g["k"] <= hi)]
    rpk = low["R"].sum() / low["k"].sum()
    g["n_reads"] = g["R"] / rpk
    for m in ["uniform", "exact", "reads"]:
        g[f"factor_{m}"] = g[f"n_{m}"] / g["k"]
    g.attrs = {"U_eff": U_eff, "reads_per_key_low": rpk}
    return g.drop(columns="ssr")


def simulate_keys(n: int, q: np.ndarray, p_u: np.ndarray, rng) -> int:
    """Observed distinct (position, UMI) keys when n molecules draw positions from q and UMIs from p_u."""
    pos = rng.choice(len(q), size=n, p=q)
    umi = rng.choice(len(p_u), size=n, p=p_u)
    return int(np.unique(pos.astype(np.int64) * len(p_u) + umi).size)


def nucleotide_umi_probs(L: int = 6, comp: dict | None = None) -> np.ndarray:
    """UMI probabilities from independent bases (the paper's model; default its composition A .23 C .24 G .21 T .32)."""
    import itertools
    comp = comp or {"A": 0.23, "C": 0.24, "G": 0.21, "T": 0.32}
    p = np.array([np.prod([comp[b] for b in t]) for t in itertools.product("ACGT", repeat=L)])
    return p / p.sum()


def simulate_grid(p_u: np.ndarray, peff_targets=(2, 5, 20, 80), n_grid=None, reps: int = 3, seed: int = 0,
                  reads_per_mol_extra: float = 0.6, cap: float = 0.95) -> pd.DataFrame:
    """Simulation with known truth. For each target effective position number, a Dirichlet position profile q
    (concentration tuned to the target); n molecules draw positions and UMIs; each molecule gets 1 + Poisson(extra)
    reads (PCR / IVT duplicates). Observed: keys k and reads per position (q estimated from reads, as in ESAT).
    Returns one row per simulation: P_eff, n, k, k/K, n_uniform, n_exact (both from the observed data)."""
    rng = np.random.default_rng(seed)
    n_grid = n_grid if n_grid is not None else np.unique(np.logspace(2, 5.3, 16).astype(int))
    U_eff = 1.0 / (p_u ** 2).sum()
    pb, wb = bin_probs(p_u, 64)
    rows = []
    for target in peff_targets:
        npos = max(int(target * 3), 3)
        for rep in range(reps):
            for _ in range(50):  # profile with an effective position number near the target
                q = rng.dirichlet(np.ones(npos) * 0.7)
                if abs(1 / (q ** 2).sum() / target - 1) < 0.25:
                    break
            for n in n_grid:
                pos = rng.choice(npos, size=n, p=q)
                umi = rng.choice(len(p_u), size=n, p=p_u)
                key = pos.astype(np.int64) * len(p_u) + umi
                k = np.unique(key).size
                r = 1 + rng.poisson(reads_per_mol_extra, size=n)
                rp = np.bincount(pos, weights=r, minlength=npos)
                rp = rp[rp > 0]
                qh = rp / rp.sum()
                P_eff = rp.sum() ** 2 / (rp ** 2).sum()
                K = P_eff * U_eff
                fr = min(k / K, cap)
                n_u = k * (-np.log1p(-fr) / fr)
                qb, qw = bin_positions(qh)
                n_e, capped = solve_exact(k, qb, qw, pb, wb, n0=n_u)
                rows.append(dict(P_eff_target=target, rep=rep, n=n, k=k, P_eff=P_eff, k_over_K=k / K,
                                 n_uniform=n_u, n_exact=n_e, exact_capped=capped, uniform_capped=k / K > cap))
    return pd.DataFrame(rows)


PROTO_LIBS = {  # the five prototype libraries (most saturated genes per cell type, plus the PBS of F108)
    "trimmed56bp_TGCAAC_F117_FRB_TNF_3ctk_S1mod": "F117 FRB TNF",
    "trimmed56bp_AGTGTC_F124_KRT_PBS_3ctk_S1mod": "F124 KRT PBS",
    "trimmed56bp_AGCTTC_F108_KRT_PBS_3ctk_S1mod": "F108 KRT PBS",
    "trimmed56bp_ACAGTG_F108_KRT_IFNG_3ctk_S1mod": "F108 KRT IFNG",
    "trimmed56bp_ACTCGA_F113_MEL_TNF_3ctk_S1mod": "F113 MEL TNF",
}
METHODS = ["uniform", "exact", "reads", "original"]
METHOD_COLORS = {"none": "#8a8984", "uniform": "#eda100", "exact": "#2a78d6", "reads": "#1baf7a",
                 "original": "#e34948", "current table": "#e34948", "per-position": "#4a3aa7"}


def load_esat_saturation(path) -> pd.DataFrame:
    """ESAT <exp>.umi.saturation.txt -> DataFrame indexed by gene (one barcode per bulk library)."""
    d = pd.read_csv(path, sep="\t")
    return d.drop_duplicates("gene").set_index("gene")


def load_gene_txt(path) -> pd.Series:
    d = pd.read_csv(path, sep="\t")
    return d.groupby("Symbol")[d.columns[-1]].sum()


def truncate_umis(keys: pd.DataFrame, L: int) -> pd.DataFrame:
    """Synthetically shortened UMIs (the paper's validation): keep the first L bases and re-collapse keys
    (gene, UMI[:L], pos), summing reads."""
    t = keys.assign(umi=keys["umi"].astype(str).str[:L])
    return t.groupby(["gene", "umi", "pos"], observed=True, as_index=False)["reads"].sum() \
            .astype({"umi": "category"})


def truncation_test(keys: pd.DataFrame, lengths=(5, 4), truth_max_frac: float = 0.05) -> pd.DataFrame:
    """For genes nearly unsaturated with the full UMIs (k/K < truth_max_frac; truth = full-length exact estimate),
    estimate the molecule number again from keys with UMIs shortened to each length; returns per gene and length
    the truth, the shortened keys and the uniform / exact / reads estimates."""
    full = estimate(keys)
    ok = full.index[(full["k_over_K"] < truth_max_frac) & (full["k"] >= 20)]
    out = []
    for L in lengths:
        tk = truncate_umis(keys, L)
        e = estimate(tk)
        npos = estimate_per_position(tk)
        g = ok.intersection(e.index)
        out.append(pd.DataFrame({"L": L, "truth": full.loc[g, "n_exact"], "k": e.loc[g, "k"],
                                 "k_over_K": e.loc[g, "k_over_K"], "uniform": e.loc[g, "n_uniform"],
                                 "exact": e.loc[g, "n_exact"], "per_position": npos.reindex(g),
                                 "reads": e.loc[g, "n_reads"],
                                 "U_eff": e.attrs["U_eff"]}, index=g))
    return pd.concat(out)


def umi_only_table(p_u: np.ndarray, y_max: int | None = None) -> np.ndarray:
    """The paper's UMI-only estimator as a lookup table: n_hat[y] for y = 0..U (U = number of UMIs) with
    E[Y | n] = sum_u [1 - (1 - p_u)^n], inverted by bisection (monotone); y = U (all UMIs seen, not estimable)
    is extrapolated as n_hat[U-1] + (n_hat[U-1] - n_hat[U-2]) * 2, a conservative stand-in for the paper's
    quadratic extrapolation."""
    U = len(p_u)
    y_max = U if y_max is None else min(y_max, U)
    L = np.log1p(-np.asarray(p_u, float))
    f = lambda n: float((1.0 - np.exp(n * L)).sum())
    out = np.zeros(y_max + 1)
    for y in range(1, min(y_max, U - 1) + 1):
        lo, hi = float(y), float(y) * 2
        while f(hi) < y:
            lo, hi = hi, hi * 2
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if f(mid) < y else (lo, mid)
        out[y] = 0.5 * (lo + hi)
    if y_max == U:
        out[U] = out[U - 1] + 2 * (out[U - 1] - out[U - 2])
    return out


def estimate_per_position(keys: pd.DataFrame, p_u: pd.Series | None = None) -> pd.Series:
    """Per-position collision correction: collisions only happen within a 5' position, where the key is the UMI
    alone, so the paper's UMI-only estimator applies to each position's key count k_p; n = sum_p n_hat(k_p).
    Needs no position profile from reads."""
    p_u = umi_frequencies(keys) if p_u is None else p_u
    kp = keys.groupby(["gene", "pos"], observed=True).size()
    tab = umi_only_table(p_u.to_numpy(), int(kp.max()))
    return pd.Series(tab[kp.to_numpy()], index=kp.index).groupby(level="gene", observed=True).sum()


def summary_box(x_labels, values, order, name: str, color: str, showlegend: bool = True, legendgroup=None):
    """A plotly Box trace from precomputed statistics (quartiles, whiskers at the 5th / 95th percentiles), so the
    page does not carry every point: x_labels and values are aligned arrays; order = categories to draw."""
    import plotly.graph_objects as go
    d = pd.DataFrame({"x": np.asarray(x_labels), "v": np.asarray(values, float)}).dropna()
    q = d.groupby("x")["v"].quantile([0.05, 0.25, 0.5, 0.75, 0.95]).unstack()
    cats = [c for c in order if c in q.index]
    q = q.loc[cats]
    return go.Box(x=cats, q1=q[0.25].round(4), median=q[0.5].round(4), q3=q[0.75].round(4),
                  lowerfence=q[0.05].round(4), upperfence=q[0.95].round(4), name=name, marker_color=color,
                  showlegend=showlegend, legendgroup=legendgroup or name, boxpoints=False)


def current_rule(e: pd.DataFrame, fallback: float = 1.5) -> pd.Series:
    """The current table's rule applied to an estimate() table: uniform when its factor <= fallback, else reads."""
    return pd.Series(np.where(e["n_uniform"] / e["k"] <= fallback, e["n_uniform"], e["n_reads"]), index=e.index)


SELF_METHODS = {"none (keys)": "k", "uniform": "n_uniform", "current table": "current", "exact": "n_exact",
                "per-position": "per_position"}


def self_consistency(keys: pd.DataFrame, lengths=(5, 4, 3), min_keys: int = 20) -> pd.DataFrame:
    """Robustness to saturation: every method is applied to the full UMIs and to UMIs shortened to each length
    (collisions become 4x, 16x, 64x more likely; the molecules are the same). Returns one row per gene x length with,
    for each method, its full-UMI estimate (<m>_full) and its shortened-UMI estimate (<m>), plus k/K after shortening.
    Genes: >= min_keys keys with the full UMIs."""
    def all_methods(k):
        e = estimate(k)
        e["current"] = current_rule(e)
        e["per_position"] = estimate_per_position(k).reindex(e.index)
        return e
    full = all_methods(keys)
    full = full[full["k"] >= min_keys]
    out = []
    for L in lengths:
        e = all_methods(truncate_umis(keys, L)).reindex(full.index)
        d = pd.DataFrame({"L": L, "k_over_K": e["k_over_K"], "k_over_K_full": full["k_over_K"]}, index=full.index)
        for lab, c in SELF_METHODS.items():
            d[f"{lab}_full"], d[lab] = full[c].astype(float), e[c].astype(float)
        out.append(d)
    return pd.concat(out)


def self_consistency_summary(sc: pd.DataFrame, sat_min: float = 0.2) -> pd.DataFrame:
    """Per UMI length and method: R^2 of log2(shortened) against log2(full) around the identity line
    (1 - sum (y - x)^2 / sum (y - mean y)^2), for all genes and for genes saturated after shortening (k/K > sat_min);
    median log2(shortened / full) (bias) and robust SD (1.4826 MAD) for the saturated genes."""
    rows = []
    for L, d in sc.groupby("L"):
        sat = d["k_over_K"] > sat_min
        for lab in SELF_METHODS:
            x, y = np.log2(d[f"{lab}_full"]), np.log2(d[lab])
            ok = np.isfinite(x) & np.isfinite(y)
            r2 = lambda m: 1 - ((y[m] - x[m]) ** 2).sum() / ((y[m] - y[m].mean()) ** 2).sum()
            lr = (y - x)[ok & sat]
            rows.append({"UMI length": L, "UMIs": 4 ** L, "method": lab, "genes": int(ok.sum()),
                         "R² all genes": round(float(r2(ok)), 4), "saturated genes (k/K > %.1f)" % sat_min:
                         int((ok & sat).sum()), "R² saturated": round(float(r2(ok & sat)), 4),
                         "bias saturated (median log2)": round(float(lr.median()), 3),
                         "spread saturated (robust SD log2)": round(float(1.4826 * (lr - lr.median()).abs().median()), 3)})
    return pd.DataFrame(rows)
