package broad.pda.seq.segmentation;


/**
 * Implements a strategy for read count normalization, so that counts estimated from a given alignment can be compared to counts
 * estimated in a different alignment.
 * @author mgarber
 *
 */
public interface ReadCountNormalizer {
	
	/**
	 * Normalizes a given read count
	 * @param count -- raw count estimated from this data model or alignment
	 * @return normalized count
	 */
	double normalize(double count);
}
