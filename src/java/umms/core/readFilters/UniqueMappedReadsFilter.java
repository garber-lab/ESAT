package umms.core.readFilters;

import umms.core.alignment.Alignment;
import org.apache.commons.collections15.Predicate;


public class UniqueMappedReadsFilter implements Predicate<Alignment>{

	@Override
	public boolean evaluate(Alignment align) {
		if(align.getWeight()<1.0)
			return false;
		return true;
	}
}
