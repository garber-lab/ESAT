package umms.core.model.score;

import java.util.*;

import umms.core.annotation.*;

public class CoverageScore extends AbstractListScore {

	//private Map<Integer,CountScore> countMap = new TreeMap<Integer,CountScore>();
	
	public CoverageScore(Annotation a) {
		super(a);
	}
	
	public CoverageScore(AnnotationCollection<? extends Annotation> model, Annotation a) {
		this(a);
		
		List<Double> scores = new ArrayList<Double>();
		CountScore.Processor processor = new CountScore.Processor(model, true);
		WindowScoreIterator<CountScore> itr = model.scan(a, 1, 0, processor);
		while (itr.hasNext()) {
			CountScore curr = itr.next();
			scores.add(curr.getCount());
			//countMap.put(Integer.valueOf(curr.getAnnotation().getStart()), curr);
		}		
		setScores(scores);
	}

	
	public static class Processor extends WindowProcessor.AbstractProcessor<CoverageScore> {
		protected AnnotationCollection<? extends Annotation> model;
		
		public Processor(AnnotationCollection<? extends Annotation> model) {
			this.model = model;
		}
		
		public CoverageScore processWindow(Annotation annotation) {
			return new CoverageScore(model, annotation);
		}

		public CoverageScore processWindow(Annotation annotation, CoverageScore previousScore) {
			return processWindow(annotation);
		}
	}
}
