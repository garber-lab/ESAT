package umms.core.readFilters;

import umms.core.alignment.Alignment;
import umms.core.annotation.Annotation;

import org.apache.commons.collections15.Predicate;

public class ReadStartsInFilter implements Predicate<Alignment> {
	
	private Annotation window;
	
	public ReadStartsInFilter(Annotation w){
		window = w;
	}
	
	@Override
	public boolean evaluate(Alignment read) {
		
		if(window.getStart()<read.getOrientedStart() && window.getEnd()>read.getOrientedStart())
			return true;
		else
			return false;
	}
}