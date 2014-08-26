package umms.core.feature;

import java.util.Collection;

import umms.core.annotation.Annotation;


public interface Window extends Annotation{

	/**
	 * Add a pointer to the annotation from which this Window originated
	 * @param annotation
	 */
	void addSourceAnnotation(Annotation annotation);

	/**
	 * @return a collection of annotations from which this originated
	 */
	Collection<? extends Annotation> getSourceAnnotations();
	
	/**
	 * Get collection of windows spanning the annotation
	 * @param windowSize Window size
	 * @param stepSize Step size
	 * @return The collection of windows
	 */
	public Collection<Window> getWindows(int windowSize, int stepSize);
	
}
