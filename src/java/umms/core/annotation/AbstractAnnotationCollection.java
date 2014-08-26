package umms.core.annotation;

import java.util.Iterator;

import umms.core.feature.Window;
import umms.core.model.score.CountScore;
import umms.core.model.score.WindowProcessor;
import umms.core.model.score.WindowScore;
import umms.core.model.score.WindowScoreIterator;
import net.sf.samtools.util.CloseableIterator;

/**
 * @author engreitz
 * Abstract class for functions common to both AlignmentModel and AnnotationSet.
 * Mostly here for consistency in function between the two classes
 * @param <T>
 */
public abstract class AbstractAnnotationCollection<T extends Annotation> implements AnnotationCollection<T> {
	

	@Override
	public int getBasesCovered(Annotation region) {
		return getBasesCovered(region, false);
	}
	
	
	@Override
	public double getTotalCoverage(Annotation region) {
		WindowScoreIterator<CountScore> itr = scan(region, 1, 0);
		double total = 0;
		while (itr.hasNext()) {
			total += itr.next().getCount();
		}
		return total;
	}
	
	@Override
	public double getAverageCoverage(Annotation region) {
		return getTotalCoverage(region) / region.length();
	}
	
	@Override
	public double getCount(Annotation region) { 
		return getCount(region, false); 
	}
	
	
	@Override
	public double getCount(AnnotationList<? extends Annotation> set) {
		return getCount(set, false);
	}
	
	
	@Override
	public double getCount(AnnotationList<? extends Annotation> set, boolean fullyContained) {
		double count = 0;
		CloseableIterator<? extends Annotation> itr = set.iterator();
		while (itr.hasNext()) {
			Annotation next = itr.next();
			count += getCount(next, fullyContained);
		}
		itr.close();
		return count;
	}
	
	
	@Override
	public CloseableIterator<T> getOverlappingAnnotations(Annotation region) {
		return getOverlappingAnnotations(region, false);
	}
	
	
	@Override
	public Annotation collapse(Annotation region, boolean fullyContained) {
		CloseableIterator<T> overlappers = getOverlappingAnnotations(region, fullyContained);
		if (!overlappers.hasNext()) return null;
		Annotation result = overlappers.next();
		while (overlappers.hasNext()) {
			result = result.union(overlappers.next());
		}
		result = result.intersect(region);
		return result;
	}

	
	public WindowScoreIterator<CountScore> scan(int windowSize, int overlap) {
		return scan(windowSize, overlap, getCountProcessor());
	}
	
	/**
	 * Iterate through the whole coordinate space in windows and score them using the provided window processor
	 * @param windowSize size of the windows to score
	 */
	public <W extends WindowScore> WindowScoreIterator<W> scan(int windowSize, int overlap, WindowProcessor<W> processor) {
		//call the coordinate space window iterator
		Iterator<? extends Window> windowIterator = getCoordinateSpace().getWindowIterator(windowSize, overlap); 
		
		//score each window
		return new WindowScoreIterator<W>(windowIterator, processor, null);
	}

	public <W extends WindowScore> WindowScoreIterator<W> scan(Annotation region, int windowSize, int overlap, WindowProcessor<W> processor) {
		Iterator<? extends Window> windowIterator = getCoordinateSpace().getWindowIterator(region, windowSize, overlap);
		return new WindowScoreIterator<W>(windowIterator, processor, region);
	}
	
	
	public WindowScoreIterator<CountScore> scan(Annotation region, int windowSize, int overlap) {
		return scan(region, windowSize, overlap, getCountProcessor());
	}
	
	
	protected WindowProcessor<CountScore> getCountProcessor() {
		return new CountScore.Processor(this);
	}

}
