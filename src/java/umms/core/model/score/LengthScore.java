package umms.core.model.score;

import java.util.*;
import umms.core.annotation.*;
import umms.core.alignment.Alignment;
import net.sf.samtools.util.CloseableIterator;

public class LengthScore extends AbstractListScore {
	
	public LengthScore(Annotation a) {
		super(a);
	}
	
	public LengthScore(AnnotationCollection<? extends Annotation> model, Annotation a) {
		this(a);
		
		List<Double> lengths = new ArrayList<Double>();
		CloseableIterator<? extends Annotation> itr = model.getOverlappingAnnotations(a);
		while (itr.hasNext()) {
			Annotation curr = itr.next();
			if (curr instanceof Alignment) {
				lengths.add(((Alignment) curr).getFragmentSize(model.getCoordinateSpace()).iterator().next().doubleValue());
			} else {
				lengths.add(Double.valueOf(itr.next().size()));
			}
		}
		itr.close();
		setScores(lengths);
	}
	
	
	public static class Processor extends WindowProcessor.AbstractProcessor<LengthScore> {
		protected AnnotationCollection<? extends Annotation> model;
		
		public Processor(AnnotationCollection<? extends Annotation> model) {
			this.model = model;
		}
		
		public LengthScore processWindow(Annotation annotation) {
			return new LengthScore(model, annotation);
		}

		public LengthScore processWindow(Annotation annotation, LengthScore previousScore) {
			return processWindow(annotation);
		}
	}
}
