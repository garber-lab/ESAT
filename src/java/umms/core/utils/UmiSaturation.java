package umms.core.utils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.log4j.Logger;

/**
 * UMI-collision (saturation) correction of collapsed UMI counts (option -umiSaturation).
 *
 * After collapsing, a gene's number of keys k (UMI, or UMI + 5' position with -umiMode position) under-counts its
 * molecules n, because different molecules that draw the same key are merged. Estimators of n:
 * <ul>
 * <li><b>exact</b>: a molecule draws position p with probability q_p (the gene's reads) and UMI u with probability p_u
 *     (the experiment's UMI frequencies) independently, so E[k | n] = sum_p sum_u [1 - (1 - q_p p_u)^n]; n solves
 *     E[k | n] = k (method of moments; Agyemang, Irizarry and Baharav, J Comput Biol 2026, with the 5' position as
 *     additional barcode letters). UMI probabilities in 64 quantile bins; positions exact up to 512 per gene, beyond
 *     that the 256 busiest kept individually and the rest binned (binPositions); safeguarded Newton iterations.</li>
 * <li><b>original</b>: the handling used for the skin response-QTL counts before this option, kept for reference:
 *     n = -K ln(1 - k/K) with all K = U_eff x P_eff keys equally likely (U_eff = 1 / sum_u p_u^2, P_eff = R^2 /
 *     sum_p r_p^2 from reads; k/K capped at 0.95), but, when that correction exceeds 1.5x, reads of the gene / reads
 *     per key of genes with 5-30 keys. Over-corrects below the switch, inherits read noise above it.</li>
 * </ul>
 * p_u is estimated from the keys of genes with at most 200 keys (collisions negligible), pooled over the barcodes
 * of the experiment; all 4^L UMIs get a pseudo-count of 0.5.
 */
public class UmiSaturation {

	static Logger logger = Logger.getLogger(UmiSaturation.class.getName());

	public static final List<String> MODES = Arrays.asList("none", "original", "exact");

	public static int maxKeysUmiFreq = 200;		// genes used to estimate the UMI frequencies
	public static int umiBins = 64;				// quantile bins for UMI probabilities
	public static int posMaxExact = 512;		// genes with up to this many positions: no position binning
	public static int posKeepTop = 256;			// otherwise the busiest positions stay individual ...
	public static int posTailBins = 64;			// ... and the remaining ones go into quantile bins
	public static double uniformCap = 0.95;		// cap on k/K for the uniform estimator
	public static double exactMinFrac = 0.002;	// below this k/K the exact estimate is set to the uniform one
	public static double maxFactor = 20.0;		// largest correction factor of the exact estimator
	public static int lowKeysMin = 5, lowKeysMax = 30;	// genes used for reads per key
	public static double originalFallback = 1.5;	// "original": reads when the formula's factor exceeds this

	/** One gene x barcode. */
	public static class Estimate {
		public String strand, chr, gene, bc;
		public int keys, positions;
		public long reads;
		public double pEff, uEff, kOverK, nUniform, nExact, nReads;
		public boolean uniformCapped, exactCapped;

		/* "original" switches to reads above the fallback factor */
		public boolean originalUsesReads() {
			return keys > 0 && nUniform / keys > originalFallback;
		}

		public double n(String mode) {
			if (mode.equals("exact")) return nExact;
			if (mode.equals("original")) return originalUsesReads() ? nReads : nUniform;
			return keys;
		}

		public double factor(String mode) {
			return keys > 0 ? n(mode) / keys : 1.0;
		}
	}

	/* key "UMI@pos" (position mode) or "UMI" (gene mode) */
	static String umiOf(String key) {
		int i = key.indexOf('@');
		return i < 0 ? key : key.substring(0, i);
	}

	static String posOf(String key) {
		int i = key.indexOf('@');
		return i < 0 ? "" : key.substring(i + 1);
	}

	/** Quantile bins of a probability vector: {bin mean probability[], number of elements[]}. */
	static double[][] bin(double[] p, int nbins) {
		double[] s = p.clone();
		Arrays.sort(s);
		int nb = Math.min(nbins, s.length);
		double[] pb = new double[nb], wb = new double[nb];
		for (int b = 0; b < nb; b++) {
			int from = (int) ((long) b * s.length / nb), to = (int) ((long) (b + 1) * s.length / nb);
			double sum = 0;
			for (int i = from; i < to; i++) sum += s[i];
			wb[b] = to - from;
			pb[b] = sum / wb[b];
		}
		return new double[][] {pb, wb};
	}

	/** Position probabilities for the exact model, {probability[], multiplicity[]}: a gene's busiest positions dominate
	 *  its collisions and are never averaged together. Up to posMaxExact positions are used as they are; beyond that
	 *  the posKeepTop largest stay individual and the rest go into posTailBins quantile bins. (Quantile bins over all
	 *  positions averaged busy positions in pairs and under-corrected genes such as RPL26, 103 positions, by up to 8%.) */
	static double[][] binPositions(double[] q) {
		double[] s = q.clone();
		Arrays.sort(s);		// ascending
		int n = s.length;
		if (n <= posMaxExact) {
			double[] w = new double[n];
			Arrays.fill(w, 1.0);
			return new double[][] {s, w};
		}
		double[][] tail = bin(Arrays.copyOfRange(s, 0, n - posKeepTop), posTailBins);
		int nt = tail[0].length;
		double[] p = new double[nt + posKeepTop], w = new double[nt + posKeepTop];
		System.arraycopy(tail[0], 0, p, 0, nt);
		System.arraycopy(tail[1], 0, w, 0, nt);
		for (int i = 0; i < posKeepTop; i++) { p[nt + i] = s[n - posKeepTop + i]; w[nt + i] = 1.0; }
		return new double[][] {p, w};
	}

	/** Experiment UMI frequencies from genes with <= maxKeysUmiFreq keys; all 4^L UMIs, pseudo-count 0.5. */
	static double[] umiFrequencies(HashMap<String, HashMap<String, HashMap<String, HashMap<String, HashMap<String, Integer>>>>> umiCount) {
		HashMap<String, Double> c = new HashMap<String, Double>();
		int L = -1;
		for (HashMap<String, HashMap<String, HashMap<String, HashMap<String, Integer>>>> s : umiCount.values())
			for (HashMap<String, HashMap<String, HashMap<String, Integer>>> chr : s.values())
				for (HashMap<String, HashMap<String, Integer>> g : chr.values())
					for (HashMap<String, Integer> bc : g.values()) {
						if (bc.size() > maxKeysUmiFreq) continue;
						for (String key : bc.keySet()) {
							String u = umiOf(key);
							if (L < 0) L = u.length();
							if (u.length() != L || u.indexOf('N') >= 0) continue;
							Double v = c.get(u);
							c.put(u, v == null ? 1.0 : v + 1.0);
						}
					}
		if (L <= 0) return new double[] {1.0};
		int K = 1 << (2 * L);
		double[] p = new double[K];
		Arrays.fill(p, 0.5);
		String nt = "ACGT";
		for (Map.Entry<String, Double> e : c.entrySet()) {
			String u = e.getKey();
			int idx = 0;
			boolean ok = true;
			for (int i = 0; i < u.length(); i++) {
				int b = nt.indexOf(u.charAt(i));
				if (b < 0) { ok = false; break; }
				idx = idx * 4 + b;
			}
			if (ok) p[idx] += e.getValue();
		}
		double tot = 0;
		for (double v : p) tot += v;
		for (int i = 0; i < K; i++) p[i] /= tot;
		return p;
	}

	/** n with E[k | n] = k; returns {n, capped (1/0)}. */
	static double[] solveExact(double k, double[] qb, double[] qw, double[] pb, double[] wb, double n0) {
		if (k <= 0) return new double[] {0, 0};
		int P = qb.length, B = pb.length;
		double[] L = new double[P * B], W = new double[P * B];
		for (int i = 0; i < P; i++)
			for (int j = 0; j < B; j++) {
				L[i * B + j] = Math.log1p(-qb[i] * pb[j]);
				W[i * B + j] = qw[i] * wb[j];
			}
		double lo = k, hi = maxFactor * k;
		if (f(hi, L, W) < k) return new double[] {hi, 1};
		double n = Math.min(Math.max(n0 > 0 ? n0 : k, lo), hi);
		for (int it = 0; it < 100; it++) {
			double fn = 0, d = 0;
			for (int t = 0; t < L.length; t++) {
				double e = Math.exp(n * L[t]);
				fn += W[t] * (1.0 - e);
				d -= W[t] * L[t] * e;
			}
			fn -= k;
			if (fn < 0) lo = n; else hi = n;
			double nn = d > 0 ? n - fn / d : Double.NaN;
			if (!(nn > lo && nn < hi)) nn = 0.5 * (lo + hi);
			if (Math.abs(nn - n) < 1e-7 * n) return new double[] {nn, 0};
			n = nn;
		}
		return new double[] {n, 0};
	}

	static double f(double n, double[] L, double[] W) {
		double s = 0;
		for (int t = 0; t < L.length; t++) s += W[t] * (1.0 - Math.exp(n * L[t]));
		return s;
	}

	/** All estimators for every gene x barcode of one experiment (umiCount: strand -> chr -> gene -> bc -> key -> reads). */
	public static List<Estimate> estimate(HashMap<String, HashMap<String, HashMap<String, HashMap<String, HashMap<String, Integer>>>>> umiCount) {
		double[] pu = umiFrequencies(umiCount);
		double sq = 0;
		for (double v : pu) sq += v * v;
		double uEff = 1.0 / sq;
		double[][] ub = bin(pu, umiBins);
		List<Estimate> out = new ArrayList<Estimate>();
		HashMap<String, long[]> low = new HashMap<String, long[]>();	// bc -> {reads, keys} of low-count genes
		for (String strand : umiCount.keySet())
			for (String chr : umiCount.get(strand).keySet())
				for (String g : umiCount.get(strand).get(chr).keySet())
					for (Map.Entry<String, HashMap<String, Integer>> bce : umiCount.get(strand).get(chr).get(g).entrySet()) {
						Estimate e = new Estimate();
						e.strand = strand; e.chr = chr; e.gene = g; e.bc = bce.getKey(); e.uEff = uEff;
						HashMap<String, Long> rp = new HashMap<String, Long>();
						for (Map.Entry<String, Integer> ke : bce.getValue().entrySet()) {
							String pos = posOf(ke.getKey());
							Long v = rp.get(pos);
							rp.put(pos, v == null ? (long) ke.getValue() : v + ke.getValue());
							e.reads += ke.getValue();
						}
						e.keys = bce.getValue().size();
						e.positions = rp.size();
						double ss = 0;
						double[] q = new double[rp.size()];
						int i = 0;
						for (long r : rp.values()) { ss += (double) r * r; q[i++] = (double) r / e.reads; }
						e.pEff = (double) e.reads * e.reads / ss;
						double K = e.pEff * uEff;
						e.kOverK = e.keys / K;
						double frac = Math.min(e.kOverK, uniformCap);
						e.uniformCapped = e.kOverK > uniformCap;
						// factor -ln(1 - f)/f applied to k: = -K ln(1 - k/K) below the cap, 3.15 k above it
						e.nUniform = frac < 1e-9 ? e.keys : e.keys * (-Math.log1p(-frac) / frac);
						e.nExact = e.nUniform;
						if (e.kOverK >= exactMinFrac) {
							double[][] qbin = binPositions(q);
							double[] r = solveExact(e.keys, qbin[0], qbin[1], ub[0], ub[1], e.nUniform);
							e.nExact = r[0];
							e.exactCapped = r[1] > 0;
						}
						if (e.keys >= lowKeysMin && e.keys <= lowKeysMax) {
							long[] a = low.get(e.bc);
							if (a == null) { a = new long[2]; low.put(e.bc, a); }
							a[0] += e.reads; a[1] += e.keys;
						}
						out.add(e);
					}
		for (Estimate e : out) {
			long[] a = low.get(e.bc);
			double rpk = (a == null || a[1] == 0) ? 1.0 : (double) a[0] / a[1];
			e.nReads = e.reads / rpk;
		}
		logger.info("UMI saturation: U_eff = " + String.format("%.1f", uEff) + " of " + pu.length + " UMIs; " + out.size()
				+ " gene x barcode estimates");
		return out;
	}

	/** Writes <exp>.umi.saturation.txt with all estimators. */
	public static void write(List<Estimate> est, File f) throws IOException {
		BufferedWriter w = new BufferedWriter(new FileWriter(f));
		w.write("barcode\tgene\tchr\tstrand\tkeys\treads\tpositions\tP_eff\tU_eff\tk_over_K\tn_exact\tn_original"
				+ "\tfactor_exact\tfactor_original\toriginal_used_reads\texact_capped");
		w.newLine();
		for (Estimate e : est) {
			w.write(e.bc + "\t" + e.gene + "\t" + e.chr + "\t" + e.strand + "\t" + e.keys + "\t" + e.reads + "\t" + e.positions
					+ "\t" + fmt(e.pEff) + "\t" + fmt(e.uEff) + "\t" + fmt(e.kOverK) + "\t" + fmt(e.nExact) + "\t"
					+ fmt(e.n("original")) + "\t" + fmt(e.factor("exact")) + "\t" + fmt(e.factor("original")) + "\t"
					+ (e.originalUsesReads() ? 1 : 0) + "\t" + (e.exactCapped ? 1 : 0));
			w.newLine();
		}
		w.close();
	}

	static String fmt(double v) {
		return String.format("%.6g", v);
	}
}
