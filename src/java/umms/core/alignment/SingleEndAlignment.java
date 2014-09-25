package umms.core.alignment;

import java.io.EOFException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;

import org.apache.log4j.Logger;


import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.TextCigarCodec;
import net.sf.samtools.SAMRecord.SAMTagAndValue;
import umms.core.coordinatesystem.CoordinateSpace;
import umms.core.feature.GenomeWindow;
import umms.core.feature.Window;
import umms.core.utils.AnnotationUtils;
import umms.core.alignment.AbstractPairedEndAlignment.TranscriptionRead;
import umms.core.annotation.*;

public class SingleEndAlignment extends BasicAnnotation implements Alignment,java.io.Serializable {
	
    /**
	 * 
	 */
	private static final long serialVersionUID = 2L;

	static Logger logger = Logger.getLogger(SingleEndAlignment.class.getName());
	
    TranscriptionRead txnRead;
	// SAMRecord will keep track of everything EXCEPT chromosome, start, end, orientation, name, blocks
	// BasicAnnotation will handle these functions.
	private SAMRecord record;
	//WrapSamRecord record;
    Collection<Annotation> splicedEdges;
	boolean hasIndel;
    /**
     * 	true: this unpaired mate was the first mate
	 *	false: this unpaired mate was second mate
	 */
	boolean isFirstMate;
	
    public SingleEndAlignment(SAMRecord read) {
    	super(read.getReferenceName(), read.getAlignmentStart()-1, read.getAlignmentStart()); //This is a dummy setup
    	parseCigar(read.getCigarString(), read.getReferenceName(), read.getAlignmentStart()-1);
    	long startT = System.nanoTime();
    	setName(read.getReadName());
		long cTime = System.nanoTime() - startT;
		cTime = Math.round(cTime/(double)1000000);
		if(cTime > 50) {
			logger.debug("parse set read name  took " + cTime);
		}
		
		startT = System.nanoTime();
    	if (read.getReadNegativeStrandFlag()) 
    		setOrientation(Strand.NEGATIVE);
    	else
    		setOrientation(Strand.POSITIVE);
		cTime = System.nanoTime() - startT;
		cTime = Math.round(cTime/(double)1000000);
		if(cTime > 50) {
			logger.debug("Set orientation took " + cTime);
		}

    	startT = System.nanoTime();
    	try {
    		record = (SAMRecord) read.clone();
    	} catch (CloneNotSupportedException e) {
    		logger.warn("Caught exception on record " + read.getReadName());
    		e.printStackTrace();
    	}
    	//record=new WrapSamRecord(read);
    	assert(read.getAlignmentEnd() == getEnd()); // sanity check
		cTime = System.nanoTime() - startT;
		cTime = Math.round(cTime/(double)1000000);
		if(cTime > 50) {
			logger.debug("Cloning object took " + cTime);
		}
    }
    
     /**
     * This is to populate our Alignment from the old legacy IGV alignments
     * @param read
     */
    public SingleEndAlignment(Alignment read) {
    	super(read.getChr(), read.getAlignmentStart(), read.getAlignmentStart()+1); //This is a dummy setup
    	parseCigar(read.toCigar(), read.getChr(), read.getAlignmentStart());
    	setName(read.getReadName());
    	
    	if (read.isNegativeStrand()) 
    		setOrientation(Strand.NEGATIVE);
    	else
    		setOrientation(Strand.POSITIVE);
    }
    
    public SingleEndAlignment(SAMRecord read,boolean isFirstMateFlag){
    	this(read);
    	setIsFirstMate(isFirstMateFlag);
    }

    public SingleEndAlignment(){
    	
    }
    
    public void setSplicedEdges(Collection<Annotation> se){
    	splicedEdges=se;
    }
    
    public void setHasIndel(boolean h){
    	this.hasIndel=h;
    }
    
    private void writeObject(java.io.ObjectOutputStream out) throws IOException{

   	//logger.info("Writing to disk");
    	out.writeObject(txnRead);
    	out.writeObject(splicedEdges);
    	out.writeBoolean(hasIndel);
    	out.writeBoolean(isFirstMate);
    	
    	out.writeObject(record);
    }
    
	private void readObject(java.io.ObjectInputStream in) throws IOException, ClassNotFoundException {

//		logger.info("Read "+this.getName()+" is being read");
//		logger.info("Attempt to read");
    	try{
    		this.setFragmentStrand((TranscriptionRead)in.readObject());
    	    this.setSplicedEdges((Collection<Annotation>)in.readObject());
    	    this.setHasIndel((boolean)in.readBoolean());
    	    this.setIsFirstMate((boolean)in.readBoolean());
//    	    System.out.println("In read object "+this.getName());
    	    //Object r=in.readO;
    	    this.record=(SAMRecord)in.readObject();
    	    
    	}
    	catch(EOFException e){
    		logger.debug(e.getMessage());
    	}
    }
    //Cigar string is used to populate the alignment blocks and read length fields
    private void parseCigar(String cigarString, String chr, int start) {
    	Cigar cigar = TextCigarCodec.getSingleton().decode(cigarString);
    	long startT = System.nanoTime();	
    	this.splicedEdges=new TreeSet<Annotation>();
    	this.hasIndel=false;
		List<CigarElement> elements=cigar.getCigarElements();
		
		int currentOffset = start;
		
		for(CigarElement element: elements){
			CigarOperator op=element.getOperator();
			int length=element.getLength();
			
			//then lets create a block from this
			if(op.equals(CigarOperator.MATCH_OR_MISMATCH)){
				long startT2 = System.nanoTime();
				int blockStart=currentOffset;
				int blockEnd=blockStart+length;
				addBlocks(new BasicAnnotation(chr, blockStart, blockEnd));
				currentOffset=blockEnd;
				long cTime = System.nanoTime() - startT2;
				cTime = Math.round(cTime/(double)1000000);
				if(cTime > 50) {
					logger.debug("parseCigar match or missmatch took " + cTime + " for cigar " + cigarString);
				}
			}
			else if(op.equals(CigarOperator.N)){
				long startT2 = System.nanoTime();
				int blockStart=currentOffset;
				int blockEnd=blockStart+length;
				splicedEdges.add(new BasicAnnotation(chr, blockStart, blockEnd));
				currentOffset=blockEnd;
				long cTime = System.nanoTime() - startT2;
				cTime = Math.round(cTime/(double)1000000);
				if(cTime > 50) {
					logger.debug("parseCigar splice edges took " + cTime + " for cigar " + cigarString);
				}
			}
			else if(op.equals(CigarOperator.INSERTION) ||  op.equals(CigarOperator.H) || op.equals(CigarOperator.DELETION)|| op.equals(CigarOperator.SKIPPED_REGION)){
				currentOffset+=length;
				this.hasIndel=true;
			}
			else{
				//TODO This needs to handle all cigar strings
				//logger.warn("We arent accounting for Cigar operator "+op.toString()+" so we are skipping read "+this.getChr()+" "+this.getAlignmentStart()+" "+this.getAlignmentEnd()+" "+cigarString+" "+this.getReadName());
			}
			
		}
		long cTime = System.nanoTime() - startT;
		cTime = Math.round(cTime/(double)1000000);
		if(cTime > 50) {
			logger.debug("parseCigar took " + cTime + " for cigar " + cigarString);
		}
		
	}

    
    public String toString(){
		return this.getReadAlignmentBlocks(null).toBED();
	}
  
    /**
     * Will set the value of the isFirstMate flag 
     * @param value true: this unpaired mate was the first mate
     * 				false: this unpaired mate was second mate
     */
    public void setIsFirstMate(boolean value){
    	isFirstMate = value;
    }
    
    public boolean getIsFirstMate(){
    	return isFirstMate;
    }

	/**
     * Returns the length of the read
     * @return
     */
    public int getReadLength(){
    	return length();
    }
    
	/**
	 * Returns the name of this read pair.
	 */
	public String getReadName() {
		return getName();
	}
	
	/**
	 * Returns the chromosome for this alignment
	 */
	public String getChromosome() {
		return getReferenceName();
	}

	/**
	 * Returns the start of the fragment
	 */
	public int getFragmentStart() {
		return getStart();
	}

	/**
	 * Returns the mapping quality of the alignment
	 */
	public int getMappingQuality() {
		return record.getMappingQuality();
	}

	/**
	 * Always returns false for objects of this class.
	 */
	public boolean isPaired() {
		return false;
	}
	
	public boolean isChimera() {
		return false;
	}

	/**
	 * Returns the end of the fragment
	 */
	public int getFragmentEnd() {
		return getEnd();
	}

	/**
	 * Returns the size of the fragment
	 */
	public Collection<Integer> getFragmentSize(CoordinateSpace C) {
		Collection<Integer> rtrn=new ArrayList<Integer>();
		
		Collection<? extends Window> fragments = this.getFragment(C);
		
		if (fragments == null) {
			rtrn.add(Integer.valueOf(this.getEnd() - this.getStart()));
			return rtrn;
		}
		
		for(Window w: fragments){
			rtrn.add(w.getSize());
		}
		
		return rtrn;
	}
	
	/**
	 * Returns an object for the fragment from start to end of alignment in the specified coordinate space
	 */
	//TODO Add optional extension factor into fragment
	public Collection<? extends Window> getFragment(CoordinateSpace C) {
		if(C==null){
			Collection<Window> rtrn=new TreeSet<Window>();
			rtrn.add(new GenomeWindow(getReferenceName(), getStart(), getEnd()));
			return rtrn;
		}
		else {return C.getFragment(getReferenceName(), getStart(), getEnd());}
	}
	

    //TODO: Add a flag for "is first of read or second of read through paired end writer" then set fragment strand based on that
	/**
	 * Returns the strand for the fragment, which is same as read strand
	 */
	public Strand getFragmentStrand() {
		return this.getOrientation();
	}
	
	/**
	 * Get the ending position of the read considering strand
	 * @return Lowest position if negative strand, highest position otherwise
	 */
	@Override
	public int getLastFragmentPositionStranded() {
		Strand strand = getFragmentStrand();
		if(strand.equals(Strand.NEGATIVE)) return getFragmentStart();
		return getFragmentEnd()-1;
	}

	/**
	 * Get the beginning position of the read considering strand
	 * @return Highest position if negative strand, lowest position otherwise
	 */
	@Override
	public int getFirstFragmentPositionStranded() {
		Strand strand = getFragmentStrand();
		if(strand.equals(Strand.NEGATIVE)) return getFragmentEnd()-1;
		return getFragmentStart();
	}

	
	/**
	 * Returns true is the alignment pair is a duplicate
	 */
	public boolean isDuplicate() {
		return record.getDuplicateReadFlag();
	}	

	/**
	 * Sets the duplicate flag
	 */
	public void setDuplicateFlag(boolean duplicateFlag) {
		record.setDuplicateReadFlag(duplicateFlag);
	}

	/**
	 * Returns the alignment Blocks 
	 */
	public Annotation getReadAlignmentBlocks(CoordinateSpace C) {
		//return this.alignmentBlocks;
		return this.copy();
	}

	@Override
	public String getChr() {
		return getReferenceName();
	}

	@Override
	public int getAlignmentStart() {
		return getStart();
	}

	@Override
	public int getAlignmentEnd() {
		return getEnd();
	}

	@Override
	public boolean isMapped() {
		return true;
	}


	@Override
	public Object getAttribute(String string) {
		return record.getAttribute(string);
		/*
		for(SAMTagAndValue tag: this.attributeMap){
			if(tag.tag.equals(string)){return tag.value;}
		}
		return null;*/
	}

	@Override
	public String getReadSequence() {
		return record.getReadString();
	}

	@Override
	public double getWeight() {
		Integer value = record.getIntegerAttribute("NH");
		if (value != null) {
			return 1.0/new Double(value);
		}
		return 1.0;
	}

	@Override
	public boolean isProperPair() {
		return false;
	}

	@Override
	public void setProperPairFlag(boolean properPairFlag) {
		throw new UnsupportedOperationException("SingleEndAlignment must be unpaired in the current implementation");
	}

	
	@Override
	public final SAMRecord toSAMRecord() {
		// Update the stored record with the properties of the Annotation, in case they have been altered since the creation
		// e.g. this happens when shuffling reads
		
		record.setAlignmentStart(getSAMStart());
		record.setReadName(getReadName());
		record.setReadNegativeStrandFlag(isNegativeStrand());
		record.setReferenceName(getReferenceName());
		
		// Return the updated record
		return record;
	}
	
	
	public String getCigarString() {
		return record.getCigarString();
	}

	@Override
	public Collection<Annotation> getReadAlignments(CoordinateSpace space) {
		Collection<Annotation> rtrn=new ArrayList<Annotation>();
		Annotation annotation=this.getReadAlignmentBlocks(space);
		rtrn.add(annotation);
		return rtrn;
	}

	@Override
	public Collection<? extends Annotation> getSpliceConnections() {
		return this.splicedEdges;
	}
	
	@Override
	public boolean equals (Annotation other){
		if(other instanceof SingleEndAlignment){
			//if name equals
			if(this.getName().equals(other.getName())){
				return super.equals(other);
			}
		}
		return false;
	}
	
	/**
	 * Returns the alignment object
	 * @return
	 */
	public Collection<Alignment> getReadMates(){
		Collection<Alignment> mates=new TreeSet<Alignment>();
		mates.add(this);
		return mates;
	}
	
	/**
	 * Sets the strand for the fragment to the strand passed as argument
	 * @param strand
	 * TODO: CHECK IF THIS IS CORRECT
	 */
	public void setFragmentStrand(TranscriptionRead strand) {
		txnRead = strand;
	}

	@Override
	public boolean hasIndel() {
		return this.hasIndel;
	}

	@Override
	public int[] getIntervalBetweenReads() {
		return null;
	}

	@Override
	public int getFragmentMidpoint(Annotation annot) {
		return AnnotationUtils.getSubAnnotationMidpointWithinAnnotation(annot, this);
	}
	
	@Override
	public void setHeader(SAMFileHeader header){
		this.record.setHeader(header);
	}
	
	public SAMFileHeader getHeader(){
		return record.getHeader();
	}

	/**
	 * This class is a wrapper for the SAMRecord class so that this class can be correctly serialized
	 * @author skadri
	 *
	 */
	public class WrapSamRecord implements java.io.Serializable{
		
		int mappingQuality;
		boolean isDuplicate;
		String readSequence;
		
		public WrapSamRecord(SAMRecord record){
			
			this.mappingQuality=record.getMappingQuality();
			this.isDuplicate=record.getDuplicateReadFlag();
			this.readSequence=record.getReadString();
		}
		
		public int getMappingQuality(){
			return this.mappingQuality;
		}
		
		public boolean getDuplicateReadFlag(){
			return this.isDuplicate;
		}
		
		public void setDuplicateReadFlag(boolean duplicateFlag) {
			this.isDuplicate = duplicateFlag;
		}
		
		public String getReadString() {
			return this.readSequence;
		}
	}
}
