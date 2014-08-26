package broad.core.annotation;

public interface GenomicAnnotationFilter<T extends LightweightGenomicAnnotation> {
	
	/**
	 * Evaluates whether current annotation meets the filtering criteria
	 * @param annotation
	 * @return true if the record meets criteria
	 */
	boolean accept(T annotation);
	
	/**
	 * Evaluates whether annotation should stop sequencial read in of further annotations
	 * @param annotaiton
	 * @return true if annotation determines that no further annotations will meet criteria.
	 * 
	 */
	boolean isEnough(T annotation);

}
