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
 * molecules n, because different molecules that draw the same key are merged.
 *
 * <b>exact</b> (method of moments; Agyemang, Irizarry and Baharav, J Comput Biol 2026, with the 5' position as
 * additional barcode letters): a molecule draws position p with probability q_p (the gene's reads, which collisions
 * do not merge) and UMI u with probability p_u (the experiment's UMI frequencies) independently, so
 * E[k | n] = sum_p sum_u [1 - (1 - q_p p_u)^n]; n solves E[k | n] = k. UMI probabilities in 64 quantile bins;
 * positions exact up to 512 per gene, beyond that the 256 busiest kept individually and the rest binned
 * (binPositions); safeguarded Newton iterations started from the uniform-key-space value -K ln(1 - k/K),
 * K = U_eff x P_eff (U_eff = 1 / sum_u p_u^2, P_eff = R^2 / sum_p r_p^2), which is also used for genes with
 * k/K below exactMinFrac, where both agree to better than 0.1%.
 * p_u is estimated from the keys of genes with at most 200 keys (collisions negligible), pooled over the barcodes
 * of the experiment; all 4^L UMIs get a pseudo-count of 0.5.
 */
public class UmiSaturation {

	static Logger logger = Logger.getLogger(UmiSaturation.class.getName());

	public static final List<String> MODES = Arrays.asList("none", "exact");

	public static int maxKeysUmiFreq = 200;		// genes used to estimate the UMI frequencies
	public static int umiBins = 64;				// quantile bins for UMI probabilities
	public static int posMaxExact = 512;		// genes with up to this many positions: no position binning
	public static int posKeepTop = 256;			// otherwise the busiest positions stay individual ...
	public static int posTailBins = 64;			// ... and the remaining ones go into quantile bins
	public static double startCap = 0.95;		// cap on k/K for the Newton starting value
	public static double exactMinFrac = 0.002;	// below this k/K the exact estimate is set to the starting value (< 0.1% apart)
	public static double maxFactor = 20.0;		// largest correction factor of the exact estimator

	/** One gene x barcode. */
	public static class Estimate {
		public String strand, chr, gene, bc;
		public int keys, positions;
		public long reads;
		public double pEff, uEff, kOverK, nExact;
		public boolean exactCapped;

		public double n(String mode) {
			return mode.equals("exact") ? nExact : keys;
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
						// starting value: uniform key space, k [-ln(1 - f) / f] with f = k/K (capped)
						double frac = Math.min(e.kOverK, startCap);
						double nStart = frac < 1e-9 ? e.keys : e.keys * (-Math.log1p(-frac) / frac);
						e.nExact = nStart;
						if (e.kOverK >= exactMinFrac) {
							double[][] qbin = binPositions(q);
							double[] r = solveExact(e.keys, qbin[0], qbin[1], ub[0], ub[1], nStart);
							e.nExact = r[0];
							e.exactCapped = r[1] > 0;
						}
						out.add(e);
					}
		logger.info("UMI saturation: U_eff = " + String.format("%.1f", uEff) + " of " + pu.length + " UMIs; " + out.size()
				+ " gene x barcode estimates");
		return out;
	}

	/** Writes <exp>.umi.saturation.txt. */
	public static void write(List<Estimate> est, File f) throws IOException {
		BufferedWriter w = new BufferedWriter(new FileWriter(f));
		w.write("barcode\tgene\tchr\tstrand\tkeys\treads\tpositions\tP_eff\tU_eff\tk_over_K\tn_exact\tfactor_exact\texact_capped");
		w.newLine();
		for (Estimate e : est) {
			w.write(e.bc + "\t" + e.gene + "\t" + e.chr + "\t" + e.strand + "\t" + e.keys + "\t" + e.reads + "\t" + e.positions
					+ "\t" + fmt(e.pEff) + "\t" + fmt(e.uEff) + "\t" + fmt(e.kOverK) + "\t" + fmt(e.nExact) + "\t"
					+ fmt(e.factor("exact")) + "\t" + (e.exactCapped ? 1 : 0));
			w.newLine();
		}
		w.close();
	}

	static String fmt(double v) {
		return String.format("%.6g", v);
	}
}
