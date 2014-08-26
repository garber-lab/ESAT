package umms.core.readFilters;

import umms.core.alignment.Alignment;

import org.apache.commons.collections15.Predicate;

public class PCRDuplicateFilter implements Predicate<Alignment> {
		
	@Override
	public boolean evaluate(Alignment align) {
		return !align.isDuplicate();
	}
}
