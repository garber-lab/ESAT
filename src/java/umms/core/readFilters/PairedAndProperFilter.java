package nextgen.core.readFilters;

import org.apache.commons.collections15.Predicate;

import nextgen.core.alignment.Alignment;

public class PairedAndProperFilter implements Predicate<Alignment> {
	@Override
	public boolean evaluate(Alignment read) {
		//check if the read is paired
		//if so, then return isProperPair
		if(read.isPaired())
			if(read.isProperPair())
				return true;
		return false;
	}
}
