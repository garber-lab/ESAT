package umms.core.sequence;

public class SequenceUtils {
	/**
	 * Computes the hamming distance between two sequences
	 * @param seq1
	 * @param seq2
	 * @return The distance between the two sequences
	 * @throws IllegalArgumentException If the sequences have different lengths or one of the Strings is null.
	 */
	public static int hamming(String seq1, String seq2) throws IllegalArgumentException{
		if (seq1 == null || seq2 == null || seq1.length() != seq2.length()) {
			throw new IllegalArgumentException("seq1 " + seq1 + " and seq2 " + seq2 + " have different lengths");
		}
		
		// compute hamming distance
		int dist = 0;
		for (int i = 0; i < seq1.length(); i++) {
			if (seq1.charAt(i) != seq2.charAt(i)) {
				dist++;
			}
		}
		
		return dist;
	}
}
