package umms.core.alignment;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.log4j.Logger;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.StringUtil;
import umms.core.coordinatesystem.CoordinateSpace;
import umms.core.feature.GenomeWindow;
import umms.core.feature.Window;
import umms.core.annotation.*;
import umms.core.utils.AnnotationUtils;
import umms.core.writers.PairedEndWriter;

/**
 * @author prussell
 *
 */
public abstract class AbstractPairedEndAlignment extends BasicAnnotation implements Alignment,java.io.Serializable {


    /**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	//First in pair alignment
    SingleEndAlignment firstMate;
    //Second in pair alignment
    SingleEndAlignment secondMate;
    //Whether first of pair or second of pair are in direction of transcription. Unstranded option in case of an unstranded library
    Map<String,String> attributeMap;
    private Strand orientation = Strand.UNKNOWN;
    boolean isProperPair;
	private static Logger logger = Logger.getLogger(FragmentAlignment.class.getName());

	/**
	 * @param b
	 */
	public AbstractPairedEndAlignment(BasicAnnotation b) {
		super(b);
	}
     
	/**
	 * Populate the attribute map with combined attribute values based on the two mates
	 */
	public void refreshAttributeMap() {
		attributeMap = new HashMap<String,String>();
		List<SAMRecord.SAMTagAndValue> firstMateAttributes = firstMate.toSAMRecord().getAttributes();
		List<SAMRecord.SAMTagAndValue> secondMateAttributes = secondMate.toSAMRecord().getAttributes();
		Map<String, String> firstMateAttributeMap = new TreeMap<String, String>();
		Map<String, String> secondMateAttributeMap = new TreeMap<String, String>();
		TreeSet<String> allTags = new TreeSet<String>();
		for(SAMRecord.SAMTagAndValue tv : firstMateAttributes) {
			String tag = tv.tag;
			String value = tv.value.toString();
			firstMateAttributeMap.put(tag, value);
			allTags.add(tag);
			//logger.info("First mate tag " + tag + " value " + value);
		}
		for(SAMRecord.SAMTagAndValue tv : secondMateAttributes) {
			String tag = tv.tag;
			String value = tv.value.toString();
			secondMateAttributeMap.put(tag, value);
			allTags.add(tag);
			//logger.info("Second mate tag " + tag + " value " + value);
		}
		for(String tag : allTags) {
			String combinedAttribute = combineAttributes(tag, firstMateAttributeMap.get(tag), secondMateAttributeMap.get(tag));
			//logger.info("Combined " + combinedAttribute);
			if(combinedAttribute != null) attributeMap.put(tag, combinedAttribute);
		}
	}
	
	/**
	 * Combine sam tag values for the two mates
	 * Currently sums integer values and ignores string values
	 * @param tag The tag (not currently used but maybe support different tags differently in the future)
	 * @param firstMateValue First mate value for tag
	 * @param secondMateValue Second mate value for tag
	 * @return The combined value or null if can't combine
	 */
	private static String combineAttributes(String tag, String firstMateValue, String secondMateValue) {
		
		//TODO support string values
		
		if(firstMateValue == null && secondMateValue == null) {
			return null;
		}
		
		if(firstMateValue == null) {
			try {
				int value2 = Integer.parseInt(secondMateValue);
				return Integer.valueOf(value2).toString();
			} catch (NumberFormatException e) {
				return null;
			}			
		}
		
		if(secondMateValue == null) {
			try {
				int value1 = Integer.parseInt(firstMateValue);
				return Integer.valueOf(value1).toString();
			} catch (NumberFormatException e) {
				return null;
			}			
		}
		
		try {
			int value1 = Integer.parseInt(firstMateValue);
			int value2 = Integer.parseInt(secondMateValue);
			return Integer.valueOf(value1 + value2).toString();
		} catch (NumberFormatException e) {}
		
		return null;
		
	}
	
   	/**
     * Returns the strand for the first of pair
     * @return
     */
    Strand getFirstOfPairStrand(){
    	//Mate will always exist because object is only created if both mates exist.
    	return this.firstMate.getFragmentStrand();
    }

    /**
     * Returns the strand for the second of pair
     * @return
     */
    Strand getSecondOfPairStrand(){
    	//Mate will always exist because object is only created if both mates exist.
    	return this.secondMate.getFragmentStrand();
    }

	/**
	 * Returns the name of this read pair.
	 */
	@Override
	public String getReadName() {
		return firstMate.getReadName();
	}

	/**
	 * Returns the chromosome for this alignment
	 * @return Reference name
	 */
	public String getChromosome() {
		return this.getReferenceName();
	}

	/**
	 * Returns the sum of the mapping qualities of both mates
	 */
	@Override
	public int getMappingQuality() {
		
		return this.firstMate.getMappingQuality()+this.secondMate.getMappingQuality();
	}

	@Override
	public int[] getIntervalBetweenReads() {
		if(firstMate.overlaps(secondMate)) return null;
		int[] rtrn = new int[2];
		if(firstMate.getEnd() < secondMate.getStart()) {
			rtrn[0] = firstMate.getEnd();
			rtrn[1] = secondMate.getStart();
			return rtrn;
		}
		rtrn[0] = secondMate.getEnd();
		rtrn[1] = firstMate.getStart();
		return rtrn;
	}
	
	/**
	 * Always returns true for objects of this class.
	 */
	@Override
	public boolean isPaired() {
		return true;
	}

	/**
     * For paired end, returns true if the read in direction of transcription is negative stranded
     * @return true if the read in direction of transcription is negative stranded
     */
	@Override
	public boolean isNegativeStrand() {
		if ("-".equals(this.getFragmentStrand())) return true;
		//Even unstranded would be false
		return false;
	}

	/**
	 * Returns true is the alignment pair is a duplicate
	 */
	@Override
	public boolean isDuplicate() {
		return firstMate.isDuplicate();
	}	

	/**
	 * Sets the duplicate flag
	 */
	@Override
	public void setDuplicateFlag(boolean duplicateFlag) {
		firstMate.setDuplicateFlag(duplicateFlag);
		secondMate.setDuplicateFlag(duplicateFlag);
	}
	
	/**
	 * Get the first read
	 * @return The first read in pair
	 */
	public Alignment getFirstMate() {
		return this.firstMate;
	}
	
	/**
	 * Get the second read
	 * @return The second read in pair
	 */
	public Alignment getSecondMate() {
		return this.secondMate;
	}
	
	@Override
	public int getFragmentMidpoint(Annotation annot) {
		return AnnotationUtils.getSubAnnotationMidpointWithinAnnotation(annot, this);
	}
	
	/**
	 * Returns an object for the fragment between the read pair in the specified coordinate space
	 */
	@Override
	public Collection<? extends Window> getFragment(CoordinateSpace C) {
		if (C==null) {
			Collection<Window> rtrn=new TreeSet<Window>();
			rtrn.add(new GenomeWindow(this.getChr(), this.getFragmentStart(), this.getFragmentEnd()));
			return rtrn;
		}
		return (C.getFragment(this.getChr(), this.getFragmentStart(), this.getFragmentEnd()));
	}

	/**
	 * Returns the strand for the fragment, depending on the transcription read
	 */
	@Override
	public Strand getFragmentStrand() {
		return orientation;
	}
	
	/**
	 * Sets the strand for the fragment to the strand of the transcription read
	 * @param transcriptionRead
	 */
	@Override
	public void setFragmentStrand(TranscriptionRead transcriptionRead) {
		if(transcriptionRead.equals(TranscriptionRead.FIRST_OF_PAIR)){
			orientation = firstMate.getFragmentStrand();
		}
		else if(transcriptionRead.equals(TranscriptionRead.SECOND_OF_PAIR)){
			orientation = secondMate.getFragmentStrand();
		}
		else
			orientation = Strand.UNKNOWN;
			
	}
	
	/**
	 * Returns the size of the fragment
	 */
	@Override
	public Collection<Integer> getFragmentSize(CoordinateSpace C) {
		Collection<Integer> rtrn=new ArrayList<Integer>();
		Collection<? extends Window> fragments = this.getFragment(C);
		
		if(fragments == null) {
			rtrn.add(Integer.valueOf(this.getEnd() - this.getStart()));
			return rtrn;
		}
		
		for(Window w: fragments){
			rtrn.add(Integer.valueOf(w.getSize()));
		}
		
		return rtrn;
	}

	/**
	 * Returns the start of the fragment
	 */
	@Override
	public int getFragmentStart() {
		return Math.min(this.firstMate.getFragmentStart(), this.secondMate.getFragmentStart());
	}

	/**
	 * Returns the end of the fragment
	 */
	@Override
	public int getFragmentEnd() {
		return Math.max(this.firstMate.getFragmentEnd(), this.secondMate.getFragmentEnd());
	}

	/**
	 * Sets the specified attribute value
	 * @param attribute 
	 * @param value 
	 */
	public void setAttribute(String attribute, String value) {
		attributeMap.put(attribute, value);
	}

	/**
	 * Returns the value of the specified attribute
	 */
	@Override
	public String getAttribute(String attribute) {
		return attributeMap.get(attribute);
	}
	
	/**
	 * This enumerator represents with of the mates in the read pair (or neither for unstranded library) is in the direction of transcription
	 * @author skadri
	 *
	 */
	public static enum TranscriptionRead implements java.io.Serializable{
		
		/**
		 * 
		 */
		FIRST_OF_PAIR,
		/**
		 * 
		 */
		SECOND_OF_PAIR,
		/**
		 * 
		 */
		UNSTRANDED;
		
		private static String FIRST_NAME = "first";
		private static String SECOND_NAME = "second";
		private static String UNSTRANDED_NAME = "unstranded";
		
		public static String commaSeparatedList() {
			String rtrn = values()[0].toString();
			for(int i=1; i<values().length; i++) {
				rtrn += ", " + values()[i].toString();
			}
			return rtrn;
		}
		
		public String toString() {
			switch(this) {
			case FIRST_OF_PAIR:
				return FIRST_NAME;
			case SECOND_OF_PAIR:
				return SECOND_NAME;
			case UNSTRANDED:
				return UNSTRANDED_NAME;
			default:
				throw new IllegalArgumentException("Transcription read " + this + " not recognized.");
			}
		}

		public static TranscriptionRead fromString(String name) {
			if(name.equals(FIRST_NAME)) {
				return FIRST_OF_PAIR;
			} else if(name.equals(SECOND_NAME)) {
				return SECOND_OF_PAIR;
			} else if(name.equals(UNSTRANDED_NAME)) {
				return UNSTRANDED;
			} else {
				throw new IllegalArgumentException("Transcription read " + name + " not recognized. Options: " + commaSeparatedList());
			}

		}
		
	}

	@Override
	public Annotation getReadAlignmentBlocks(CoordinateSpace C) {
		/*Collection<Annotation> blocks=new TreeSet<Annotation>();
		blocks.addAll(this.firstMate.getReadAlignmentBlocks(C).getBlocks());
		blocks.addAll(this.secondMate.getReadAlignmentBlocks(C).getBlocks());
		Annotation alignmentBlock=new Gene(blocks, this.firstMate.getName());
		return alignmentBlock;*/
		return firstMate.union(secondMate);
	}

	public Collection<Alignment> getReadMates(){
		Collection<Alignment> mates=new TreeSet<Alignment>();
		mates.add(getFirstMate());
		mates.add(getSecondMate());
		return mates;
	}
	@Override
	public String toString(){
		/*
		 * add name into the bed string by @zhuxp
		 */
		String bed = this.getReadAlignmentBlocks(null).toBED();
		String[] a=bed.split("\t");
		a[3]=this.getName();
		return StringUtil.join("\t",a);
		
	}

	@Override
	public int getStart() {
		return this.getAlignmentStart();
	}

	@Override
	public int getEnd() {
		return this.getAlignmentEnd();
	}

	@Override
	public String getChr() {
		return this.getReferenceName();
	}

	@Override
	public int getAlignmentStart() {
		return Math.min(this.firstMate.getAlignmentStart(), this.secondMate.getAlignmentStart());
	}

	@Override
	public int getAlignmentEnd() {
		return Math.max(this.firstMate.getAlignmentEnd(), this.secondMate.getAlignmentEnd());
	}

	@Override
	public boolean isMapped() {
		return true;
	}
	
	@Override
	public boolean isChimera() {
		return !(firstMate.getReferenceName().equals(secondMate.getReferenceName()));
	}

	@Override
	public String getReadSequence() {
		return this.firstMate.getReadSequence()+this.secondMate.getReadSequence();
	}

	@Override
	public double getWeight() {
		if(this.attributeMap.containsKey("NH")){
			return 1.0/new Double(this.attributeMap.get("NH").toString()).doubleValue();
		}
		return 1.0;
	}

	@Override
	public boolean isProperPair() {
		return this.isProperPair;
	}

	@Override
	public void setProperPairFlag(boolean properPairFlag) {
		this.isProperPair=properPairFlag;
	}
	
	@Override
	public void shift(int delta) {
		super.shift(delta);
		firstMate.shift(delta);
		secondMate.shift(delta);
	}
	
	@Override
	public void moveToCoordinate(int coord) {
		shift(coord - getStart());
	}
	
	//TODO should this be implemented separately in the two subclasses?
	@Override
	public final SAMRecord toSAMRecord() {
		SAMRecord record = firstMate.toSAMRecord();
		
		record.setAttribute(PairedEndWriter.readStartFlag, Integer.valueOf(firstMate.getSAMStart()));
		record.setAttribute(PairedEndWriter.readCigarFlag, firstMate.getCigarString());
		record.setAttribute(PairedEndWriter.mateSequenceFlag, secondMate.getReadSequence());
		record.setAttribute(PairedEndWriter.mateCigarFlag, secondMate.getCigarString());
		
		record.setMateAlignmentStart(secondMate.getSAMStart());
		// add by @zhuxp
        record.setAlignmentStart(this.getAlignmentStart()+1);
     
		// end of add (test version)
        
        Annotation fragment = getReadAlignmentBlocks(null);
		record.setCigarString(fragment.getLengthOnReference() + "M");  // NOTE: losing information about indels in the SingleEndAlignments	
		record.setInferredInsertSize(fragment.getLengthOnReference());

		return record;
	}
	
	@Override
	public Strand getOrientation(){
		return getFragmentStrand();
	}

	@Override
	public Collection<Annotation> getReadAlignments(CoordinateSpace space) {
		Collection<Annotation> rtrn=new ArrayList<Annotation>();
		rtrn.add(firstMate.getReadAlignmentBlocks(space));
		rtrn.add(secondMate.getReadAlignmentBlocks(space));
		//set orientation to strand orientation
		for(Annotation m:rtrn){
			m.setOrientation(getFragmentStrand());
		}
		return rtrn;
	}

	@Override
	public boolean hasIndel() {
		return this.firstMate.hasIndel() || this.secondMate.hasIndel();
	}
	
	/**
	 * Get the beginning position of the fragment considering strand
	 * @return Highest position if negative strand, lowest position otherwise
	 */
	@Override
	public int getFirstFragmentPositionStranded() {
		Strand strand = getFragmentStrand();
		if(strand.equals(Strand.UNKNOWN)) {
			throw new IllegalStateException("Fragment strand is unknown: " + toSAMRecord().getSAMString());
		}
		if(strand.equals(Strand.NEGATIVE)) return getFragmentEnd()-1;
		return getFragmentStart();
	}

	/**
	 * Get the ending position of the fragment considering strand
	 * @return Lowest position if negative strand, highest position otherwise
	 */
	@Override
	public int getLastFragmentPositionStranded() {
		Strand strand = getFragmentStrand();
		if(strand.equals(Strand.UNKNOWN)) {
			throw new IllegalStateException("Fragment strand is unknown: " + toSAMRecord().getSAMString());
		}
		if(strand.equals(Strand.NEGATIVE)) return getFragmentStart();
		return getFragmentEnd()-1;
	}
	
	/**
	 * Get the alignment as an Annotation object
	 * @param firstMate First mate
	 * @param secondMate Second mate
	 * @param fullFragment Whether to get full fragment as one block or separate reads as blocks
	 * @return An annotation consisting of the aligned coordinates
	 */
	protected static BasicAnnotation asAnnotation(SingleEndAlignment firstMate, SingleEndAlignment secondMate, boolean fullFragment) {
		
		String chr = firstMate.getChr();
		if(!chr.equals(secondMate.getChr())) {
			throw new IllegalArgumentException("Can't make annotation from alignments on different chromosomes:\n" + firstMate.toBED() + "\n" + secondMate.toBED());
		}
		
		if(fullFragment) {
			int start = Math.min(firstMate.getStart(), secondMate.getStart());
			int end = Math.max(firstMate.getEnd(), secondMate.getEnd());
			BasicAnnotation ann = new BasicAnnotation(chr, start, end);
			ann.setName(firstMate.getName());
			return ann;
		}
		
		Collection<Annotation> blocks = new ArrayList<Annotation>();
		blocks.add(firstMate);
		blocks.add(secondMate);
		return new BasicAnnotation(blocks);
		
	}
	
	@Override
	public void setHeader(SAMFileHeader header){
		getFirstMate().setHeader(header);
		getSecondMate().setHeader(header);
	}
	
	public SAMFileHeader getHeader(){
		return getFirstMate().getHeader();
	}

	
}
