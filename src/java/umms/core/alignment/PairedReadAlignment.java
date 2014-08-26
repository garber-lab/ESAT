package umms.core.alignment;

import java.util.HashMap;
import org.apache.log4j.Logger;

/**
 * @author prussell
 * Alignment of two paired reads without the fragment in between
 */
public class PairedReadAlignment extends AbstractPairedEndAlignment {

	private static Logger logger = Logger.getLogger(PairedReadAlignment.class.getName());

	
    /**
     * Constructs a paired end alignment object from two alignments.
     * @param read1
     * @param read2
     */
    public PairedReadAlignment(SingleEndAlignment read1, SingleEndAlignment read2) {
    	this(read1, read2, TranscriptionRead.UNSTRANDED);
    }
    
    /**
     * Constructs a paired end alignment object from two alignments and provides the read that is in the direction of transcription
     * @param read1
     * @param read2
     * @param strand
     */
    public PairedReadAlignment(SingleEndAlignment read1, SingleEndAlignment read2,TranscriptionRead strand) {
    	super(asAnnotation(read1, read2, false));   //TODO How to deal with pairs that have different chromosomes?

    	this.firstMate = read1;
        this.secondMate = read2;
        setFragmentStrand(strand);
        firstMate.setOrientation(this.getFragmentStrand());
        secondMate.setOrientation(this.getFragmentStrand());
        refreshAttributeMap();
    }
    
    @Override
	public boolean equals(Object o) {
		if(!o.getClass().equals(PairedReadAlignment.class)) {
			return false;
		}
		PairedReadAlignment f = (PairedReadAlignment)o;
		return f.getFirstMate().equals(getFirstMate()) && f.getSecondMate().equals(getSecondMate());
	}

	@Override
	public int hashCode() {
		return (getFirstMate().toString() + getSecondMate().toString()).hashCode();
	}
	
	
	

}
