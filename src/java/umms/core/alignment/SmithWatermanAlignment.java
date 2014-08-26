package umms.core.alignment;

import jaligner.Sequence;
import jaligner.SmithWatermanGotoh;
import jaligner.matrix.Matrix;
import jaligner.matrix.MatrixGenerator;

public class SmithWatermanAlignment {
	
	private Sequence sequence;
	private Matrix scoringMatrix;
	private float gapOpen;
	private float gapExtend;
	
	/**
	 * Default match score for Smith Waterman
	 */
	public static float DEFAULT_MATCH_SCORE = 5;
	
	/**
	 * Default mismatch score for Smith Waterman
	 */
	public static float DEFAULT_MISMATCH_SCORE = -4;
	
	/**
	 * Default gap open penalty for Smith Waterman
	 */
	public static float DEFAULT_GAP_OPEN_PENALTY = 8;
	
	/**
	 * Default gap extend penalty for Smith Waterman
	 */
	public static float DEFAULT_GAP_EXTEND_PENALTY = 2;

	
	public SmithWatermanAlignment(String seq) {
		this(seq, DEFAULT_MATCH_SCORE, DEFAULT_MISMATCH_SCORE, DEFAULT_GAP_OPEN_PENALTY, DEFAULT_GAP_EXTEND_PENALTY);
	}
	
	public SmithWatermanAlignment(String seq, float matchScore, float mismatchScore, float gapOpenPenalty, float gapExtendPenalty) {
		sequence = new Sequence(seq);
		scoringMatrix = MatrixGenerator.generate(matchScore, mismatchScore);
		gapOpen = gapOpenPenalty;
		gapExtend = gapExtendPenalty;
	}
	
	public jaligner.Alignment align(String seq) {
		return SmithWatermanGotoh.align(sequence, new Sequence(seq), scoringMatrix, gapOpen, gapExtend);
	}
	
	public static jaligner.Alignment align(String seq1, String seq2, float matchScore, float mismatchScore, float gapOpenPenalty, float gapExtendPenalty) {
		Matrix matrix = MatrixGenerator.generate(matchScore, mismatchScore);
		return SmithWatermanGotoh.align(new Sequence(seq1), new Sequence(seq2), matrix, gapOpenPenalty, gapExtendPenalty);
	}
	
	
	public static String getFullPrintableAlignment(jaligner.Alignment alignment) {
		int start1 = alignment.getStart1();
		int start2 = alignment.getStart2();
		int maxStart = Math.max(start1, start2);
		int minStart = Math.min(start1, start2);
		String m = new String(alignment.getMarkupLine());
		String markup = "";
		for(int i = 0; i < maxStart; i++) {
			markup += " ";
		}
		markup += m;
		int padding = maxStart - minStart;
		String extra = "";
		for(int i = 0; i < padding; i++) {
			extra += " ";
		}
		String orig1 = alignment.getOriginalSequence1().getSequence();
		String orig2 = alignment.getOriginalSequence2().getSequence();
		String rtrn = (maxStart == start1 ? orig1 : (extra + orig1)) + "\n";
		rtrn += markup + "\n";
		rtrn += maxStart == start2 ? orig2 : (extra + orig2);
		return rtrn;
	}
	
	
	/**
	 * Get start position on first sequence of an ungapped Smith-Waterman alignment
	 * First sequence must be at least as long as second sequence
	 * @param seq1 Longer sequence
	 * @param seq2 Shorter sequence
	 * @param matchScore Match score
	 * @param mismatchScore Mismatch score
	 * @param gapOpenPenalty Gap open penalty
	 * @param gapExtendPenalty Gap extend penalty
	 * @param minPctIdentity Min percent identity for alignment
	 * @return Start position of alignment on first sequence, or -1 if no suitable alignment exists
	 */
	public static int ungappedMatchStartOnFirstSequence(String seq1, String seq2, float matchScore, float mismatchScore, float minPctIdentity) {
		if(seq1.length() < seq2.length()) {
			throw new IllegalArgumentException("Seq1 must be at least as long as seq2");
		}
		return matchStartOnFirstSequence(seq1, seq2, matchScore, mismatchScore, Float.MAX_VALUE, Float.MAX_VALUE, minPctIdentity);
	}
	
	/**
	 * Get start position on first sequence of the Smith-Waterman alignment
	 * @param seq1 First sequence
	 * @param seq2 Second sequence
	 * @param matchScore Match score
	 * @param mismatchScore Mismatch score
	 * @param gapOpenPenalty Gap open penalty
	 * @param gapExtendPenalty Gap extend penalty
	 * @param minPctIdentity Min percent identity for alignment
	 * @return Start position of alignment on first sequence, or -1 if no suitable alignment exists
	 */
	public static int matchStartOnFirstSequence(String seq1, String seq2, float matchScore, float mismatchScore, float gapOpenPenalty, float gapExtendPenalty, float minPctIdentity) {
		jaligner.Alignment align = align(seq1, seq2, matchScore, mismatchScore, gapOpenPenalty, gapExtendPenalty);
		float pctId = align.getIdentity() / (float) seq2.length();
		if(pctId < minPctIdentity) {
			return -1;
		}
		return align.getStart1();
	}
	
	/**
	 * Whether there is a full length ungapped match of sequence 2 to sequence 1, starting at specified position of sequence 1, allowing mismatches
	 * Uses default match and mismatch scores
	 * @param seq1 Sequence 1
	 * @param seq2 Sequence 2
	 * @param seq1start Required start position of match on sequence 1
	 * @return Whether the best ungapped alignment is full length and has at most the max number of mismatches
	 */
	public static boolean containsFullLengthUngappedMatch(String seq1, String seq2, int seq1start, int maxMismatches) {
		return containsFullLengthUngappedMatch(seq1, seq2, seq1start, DEFAULT_MATCH_SCORE, DEFAULT_MISMATCH_SCORE, maxMismatches);
	}

	/**
	 * Whether there is a full length ungapped match of sequence 2 to sequence 1, starting at specified position of sequence 1, allowing mismatches
	 * @param seq1 Sequence 1
	 * @param seq2 Sequence 2
	 * @param seq1start Required start position of match on sequence 1
	 * @param matchScore Match score
	 * @param mismatchScore Mismatch score
	 * @param maxMismatches Max allowable number of mismatches
	 * @return Whether the best ungapped alignment is full length and has at most the max number of mismatches
	 */
	public static boolean containsFullLengthUngappedMatch(String seq1, String seq2, int seq1start, float matchScore, float mismatchScore, int maxMismatches) {
		int seq1end = seq1start + seq2.length();
		if(seq1end > seq1.length()) {
			throw new IllegalArgumentException("Sequence is too short (" + seq1.length() + ") for a match of length " + seq2.length() + " starting at position " + seq1start + ".");
		}
		String seq1substring = seq1.substring(seq1start, seq1end);
		return containsFullLengthUngappedMatch(seq1substring, seq2, matchScore, mismatchScore, maxMismatches);
	}

	
	/**
	 * Whether there is a full length ungapped match of the shorter sequence to the longer sequence, allowing mismatches
	 * Uses default match and mismatch scores
	 * @param seq1 Sequence 1
	 * @param seq2 Sequence 2
	 * @param maxMismatches Max allowable number of mismatches
	 * @return Whether the best ungapped alignment is full length and has at most the max number of mismatches
	 */
	public static boolean containsFullLengthUngappedMatch(String seq1, String seq2, int maxMismatches) {
		return containsFullLengthUngappedMatch(seq1, seq2, DEFAULT_MATCH_SCORE, DEFAULT_MISMATCH_SCORE, maxMismatches);
	}
	
	/**
	 * Whether there is a full length ungapped match of the shorter sequence to the longer sequence, allowing mismatches
	 * @param seq1 Sequence 1
	 * @param seq2 Sequence 2
	 * @param matchScore Match score
	 * @param mismatchScore Mismatch score
	 * @param maxMismatches Max allowable number of mismatches
	 * @return Whether the best ungapped alignment is full length and has at most the max number of mismatches
	 */
	public static boolean containsFullLengthUngappedMatch(String seq1, String seq2, float matchScore, float mismatchScore, int maxMismatches) {
		jaligner.Alignment align = align(seq1, seq2, matchScore, mismatchScore, Float.MAX_VALUE, Float.MAX_VALUE);
		int fullLen = Math.min(seq1.length(), seq2.length());
		int mismatches = fullLen - getNumberOfMatches(align);
		return mismatches <= maxMismatches;
	}
	
	private static int getNumberOfMatches(jaligner.Alignment align) {
		//System.err.println("In here");
		int matches = 0;
		char c1, c2; // the next character
		for (int i = 0; i <align.getSequence1().length; i++) {
			c1 = align.getSequence1()[i];
			c2 = align.getSequence2()[i];
			if (c1 != jaligner.Alignment.GAP && c2!=jaligner.Alignment.GAP) {
				if(c1==c2){matches++;}
			}
		}
		return matches;
	}

	
}
