package umms.core.annotation.filter;

import java.util.ArrayList;
import java.util.List;
import umms.core.annotation.*;

import org.apache.commons.collections15.Predicate;

/**
 * @author engreitz
 * Return annotations that overlap the input region
 */
public class OverlapFilter implements Predicate<Annotation> {
	
	 // Query annotations pass if they overlap any of the annotations in this list
	private List<Annotation> regions = new ArrayList<Annotation>(); 
	
	public OverlapFilter(Annotation region) {
		regions.add(region);
	}
	
	public OverlapFilter(List<? extends Annotation> regions) {
		this.regions.addAll(regions);
	}
	
	public OverlapFilter(AnnotationList<? extends Annotation> regions) {
		this.regions.addAll(regions.toList());
	}
	
	@Override
	public boolean evaluate(Annotation query) {
		for (Annotation region : regions) {
			if (region.overlaps(query)) return true;
		}
		return false;
	}
}