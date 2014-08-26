package umms.core.annotation.filter;

import java.util.ArrayList;
import java.util.List;

import umms.core.annotation.Annotation;
import umms.core.annotation.AnnotationList;

import org.apache.commons.collections15.Predicate;

/**
 * @author engreitz
 * Return only annotations that are fully contained by provided regions
 */
public class FullyContainedFilter implements Predicate<Annotation> {
	
	 // Query annotations pass if they are fully contained in any of the annotations in this list
	private List<Annotation> regions = new ArrayList<Annotation>(); 
	
	public FullyContainedFilter(Annotation region) {
		regions.add(region);
	}
	
	public FullyContainedFilter(List<? extends Annotation> regions) {
		this.regions.addAll(regions);
	}
	
	public FullyContainedFilter(AnnotationList<? extends Annotation> regions) {
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
