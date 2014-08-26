package umms.core.alignment;

import java.util.HashMap;

import umms.core.annotation.BasicAnnotation;

import org.apache.log4j.Logger;

/**
 * @author prussell
 * A full fragment alignment constructed from two paired reads
 */
public class FragmentAlignment extends AbstractPairedEndAlignment {

	private static Logger logger = Logger.getLogger(PairedReadAlignment.class.getName());

	
    /**
     * Constructs a paired end alignment object from two alignments.
     * @param read1
     * @param read2
     */
    public FragmentAlignment(SingleEndAlignment read1, SingleEndAlignment read2) {
    	this(read1, read2, TranscriptionRead.UNSTRANDED);
    }
    
    /**
     * Constructs a paired end alignment object from two alignments and provides the read that is in the direction of transcription
     * @param read1
     * @param read2
     * @param transcriptionRead
     */
    public FragmentAlignment(SingleEndAlignment read1, SingleEndAlignment read2,TranscriptionRead transcriptionRead) {
    	super(asAnnotation(read1, read2, true));   //TODO How to deal with pairs that have different chromosomes?

    	this.firstMate = read1;
        this.secondMate = read2;
        setFragmentStrand(transcriptionRead);
        refreshAttributeMap();
    }
    
    @Override
	public boolean equals(Object o) {
		if(!o.getClass().equals(FragmentAlignment.class)) {
			return false;
		}
		FragmentAlignment f = (FragmentAlignment)o;
		return f.getFirstMate().equals(getFirstMate()) && f.getSecondMate().equals(getSecondMate());
	}

	@Override
	public int hashCode() {
		return (getFirstMate().toString() + getSecondMate().toString()).hashCode();
	}

}
