package umms.core.model.score;

import java.util.*;
import umms.core.annotation.*;
import net.sf.samtools.util.CloseableIterator;

/**
 * @author engreitz
 * This class calculates statistics (sum, average, stddev, median, etc.) on the scores 
 * of annotations in the given region, in addition to counting the number of annotations.
 */
public class SumScore extends AbstractListScore {
	
	public SumScore(Annotation a) {
		super(a);
	}
	
	public SumScore(AnnotationCollection<? extends Annotation> model, Annotation region) {
		this(region);
		
		List<Double> scores = new ArrayList<Double>();
		CloseableIterator<? extends Annotation> itr = model.getOverlappingAnnotations(region);
		while (itr.hasNext()) {
			Annotation curr = itr.next();
			scores.add(curr.getScore());
		}
		itr.close();
		setScores(scores);
	}
	
	@Override
	public double getScore() {
		return getSum();
	}
	
	public static class Processor extends WindowProcessor.AbstractProcessor<SumScore> {
		protected AnnotationCollection<? extends Annotation> model;
		
		public Processor(AnnotationCollection<? extends Annotation> model) {
			this.model = model;
		}
		
		public SumScore processWindow(Annotation annotation) {
			return new SumScore(model, annotation);
		}
		
		public SumScore processWindow(Annotation annotation, SumScore previousScore) {
			return processWindow(annotation);
		}
	}
}
