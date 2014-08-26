package umms.core.readFilters;

import umms.core.alignment.Alignment;
import umms.core.annotation.Annotation;

import org.apache.commons.collections15.Predicate;

public class SplicedReadFilter implements Predicate<Alignment>{

	@Override
	public boolean evaluate(Alignment align) {
		//if has splice junctions --> return true
		if(!align.getSpliceConnections().isEmpty()){
			return true;
		}
		return false;
	}

	public boolean equals(Object other){
		if(other instanceof SplicedReadFilter){
			return true;
		}
		return false;
	}

	@Override
	public int hashCode() {
		throw new UnsupportedOperationException("Class should override hashCode() because it overrides equals()");
	}
	
}
