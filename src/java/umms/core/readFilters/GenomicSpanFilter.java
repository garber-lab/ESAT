/**
 * 
 */
package nextgen.core.readFilters;

import nextgen.core.alignment.Alignment;

import org.apache.commons.collections15.Predicate;

/**
 * @author prussell
 *
 */
public class GenomicSpanFilter implements Predicate<Alignment> {
	
	private int maxSpan;
	
	/**
	 * Instantiate the filter
	 * @param maxGenomicSpan Maximum genomic span for a fragment (inclusive)
	 */
	public GenomicSpanFilter(int maxGenomicSpan) {
		maxSpan = maxGenomicSpan;
	}
	
	@Override
	public boolean evaluate(Alignment align) {

		if(align.getFragmentEnd() - align.getFragmentStart() + 1 > maxSpan) return false;
		return true;
		
	}

}
