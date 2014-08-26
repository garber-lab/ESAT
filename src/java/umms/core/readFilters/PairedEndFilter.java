package umms.core.readFilters;

import java.util.Iterator;

import org.apache.commons.collections15.Predicate;

import umms.core.alignment.Alignment;
import umms.core.annotation.Annotation;

public class PairedEndFilter implements Predicate<Alignment> {
	@Override
	public boolean evaluate(Alignment read) {
		//check if the read is paired and both mates on same chromosome
		if(read.isPaired()){
			Iterator<Annotation> mates = read.getReadAlignments(null).iterator();
			
			if(mates.next().getReferenceName().equals(mates.next().getReferenceName()))
				return true;
		}
/*			String refName = null;
			while(mates.hasNext()){
				if(refName==null)
					refName = mates.next().getReferenceName();
				else
					if(refName.equals(mates.next().getReferenceName())){
						return true;
					}
			}
		}*/
		return false;
	}
}
