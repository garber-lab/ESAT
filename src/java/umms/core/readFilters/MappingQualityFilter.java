package umms.core.readFilters;

import umms.core.alignment.Alignment;
import org.apache.commons.collections15.Predicate;

public class MappingQualityFilter implements Predicate<Alignment> {
	
	private double minSingleEndQuality, minPairedEndQuality;
	
	public MappingQualityFilter(double singleEndQuality, double pairedEndQuality) {
		minSingleEndQuality = singleEndQuality;
		minPairedEndQuality = pairedEndQuality;
	}
	
	public MappingQualityFilter(double minQuality) {
		this(minQuality, minQuality);
	}
	
	@Override
	public boolean evaluate(Alignment align) {
		return ((align.isPaired() && align.getMappingQuality() >= minPairedEndQuality) ||
				(!align.isPaired() && align.getMappingQuality() >= minSingleEndQuality));
	}
}
