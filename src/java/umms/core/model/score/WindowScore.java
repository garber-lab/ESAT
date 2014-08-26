package umms.core.model.score;

import umms.core.annotation.Annotation;

/**
 * Represents a window scoring object.
 * All of the code for calculating scores should be included in implementers of this class
 */
public interface WindowScore {
	
	/**
	 * Returns the underlying window
	 */
	public Annotation getAnnotation();
	
	/**
	 * Returns the number of reads overlapping this region
	 * @return number of reds overlapping this window
	 */
	//public double getCount();
	
	/**
	 * Returns the "score" as defined by each subclass (e.g. count, ratio, sum of annotation scores, etc.)
	 */
	public double getScore();
	
	public abstract class AbstractWindowScore implements WindowScore {
		protected Annotation annotation;
		public AbstractWindowScore(Annotation t) {
			annotation = t;
		}
		
		public Annotation getAnnotation() { return annotation; }
	}
}
