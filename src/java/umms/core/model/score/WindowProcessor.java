package umms.core.model.score;

import umms.core.annotation.Annotation;

/**
 * Every WindowScore superclass can also implement a WindowProcessor class
 * that will wrap the WindowScore class so that it can be passed to the
 * WindowScoringIterator.
 */
public interface WindowProcessor<S extends WindowScore> {
	public S processWindow(Annotation annotation);
	public S processWindow(Annotation annotation, S previousScore);

	public void initRegion(Annotation region);
	public void finishedRegion();

	public abstract class AbstractProcessor<S extends WindowScore> implements WindowProcessor<S> {
		public void initRegion(Annotation region) {}
		public void finishedRegion() {}		
	}
}
