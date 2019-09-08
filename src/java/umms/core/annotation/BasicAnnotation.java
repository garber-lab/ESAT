package umms.core.annotation;

import broad.core.parser.StringParser;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Stack;
import java.util.TreeSet;

import org.apache.log4j.Logger;

import com.sleepycat.persist.model.Persistent;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.TextCigarCodec;
import umms.core.general.TabbedReader;
import umms.core.annotation.Annotation;
import broad.core.error.ParseException;


/**
 * @author engreitz
 * Coordinate system:  0-based, including first base but not the last
 */
@Persistent
public class BasicAnnotation extends AbstractAnnotation implements java.io.Serializable {
	protected CompoundInterval blocks = new CompoundInterval();
	private String referenceName;
	private Strand orientation = Strand.UNKNOWN;
	private String name = null;
	private double score;
	private static Logger logger = Logger.getLogger(BasicAnnotation.class.getName());
	
	/********************************************************************************
	 * CONSTRUCTORS
	 */
	
	/**
	 * For Berkeley DB only
	 * Do not use this constructor
	 */
	public BasicAnnotation() {}
	
	/**
	 * Create an annotation from a sam record
	 * The blocks are the mapped blocks
	 * This read only; does not consider the mate
	 * @param samRecord Annotation with name equal to the read name, orientation and blocks as specified in the alignment
	 */
	public BasicAnnotation(SAMRecord samRecord) {
		
		if(samRecord.getReadUnmappedFlag()) {
			throw new IllegalArgumentException("Can't make annotation from unmapped read");
		}
		
		String cigarString = samRecord.getCigarString();
		int start = samRecord.getAlignmentStart();
		String chr = samRecord.getReferenceName();
		boolean neg = samRecord.getReadNegativeStrandFlag();
		
		setName(samRecord.getReadName());
		setReferenceName(chr);
		setOrientation(neg ? Strand.NEGATIVE : Strand.POSITIVE);
		
    	Cigar cigar = TextCigarCodec.getSingleton().decode(cigarString);
 		List<CigarElement> elements=cigar.getCigarElements();
		
		int currentOffset = start;
		
		for(CigarElement element: elements){
			CigarOperator op=element.getOperator();
			int length=element.getLength();
			
			if(op.equals(CigarOperator.MATCH_OR_MISMATCH)){
				int blockStart=currentOffset;
				int blockEnd=blockStart+length;
				addBlocks(new BasicAnnotation(chr, blockStart, blockEnd));
				currentOffset=blockEnd;
			}
			else if(op.equals(CigarOperator.N)){
				int blockStart=currentOffset;
				int blockEnd=blockStart+length;
				currentOffset=blockEnd;
			}
			else if(op.equals(CigarOperator.INSERTION) ||  op.equals(CigarOperator.H) || op.equals(CigarOperator.DELETION)|| op.equals(CigarOperator.SKIPPED_REGION)){
				currentOffset+=length;
			}
			else{
				//TODO This needs to handle all cigar strings
				//logger.warn("We arent accounting for Cigar operator "+op.toString()+" so we are skipping read "+this.getChr()+" "+this.getAlignmentStart()+" "+this.getAlignmentEnd()+" "+cigarString+" "+this.getReadName());
			}
			
		}

	}
	
	public BasicAnnotation(String ucsc) {
		String ref=ucsc.split(":")[0];
		String start=ucsc.split(":")[1].split("-")[0];
		String end=ucsc.split(":")[1].split("-")[1];
		setReferenceName(ref);
		blocks.addInterval(new Integer(start), new Integer(end));
		score = 0;
	}
	
	public BasicAnnotation(String referenceName, int start, int end, Strand orientation, String name) {
		long startT = System.nanoTime();
		score = 0;
		setName(name);
		setReferenceName(referenceName);
		setOrientation(orientation);
		try {
			blocks.addInterval(start, end);
		} catch(Exception e) {
			e.printStackTrace();
			throw new IllegalArgumentException(toBED());
		}
		long cTime = System.nanoTime() - startT;
		cTime = Math.round(cTime/(double)1000000);
		if(cTime > 50) {
			logger.debug("Building BasicAnnotation took " + cTime);
		}
		
	}
	
	public BasicAnnotation(String referenceName, Strand orientation, String name, Collection<? extends Annotation> blocks) {
		score = 0;
		setName(name);
		setReferenceName(referenceName);
		setOrientation(orientation);
		addBlocks(blocks);
	}
	
	public BasicAnnotation(String referenceName, int start, int end) {
		this(referenceName, start, end, Strand.UNKNOWN);
	}
	
	public BasicAnnotation(String referenceName, int start, int end, Strand orientation) {
		this(referenceName, start, end, orientation, "");
	}
	
	public BasicAnnotation(String referenceName, int start, int end, String orientation) {
		this(referenceName, start, end, AbstractAnnotation.getStrand(orientation));
	}
	
	public BasicAnnotation(BasicAnnotation other) {
		setReferenceName(other.getReferenceName());
		setOrientation(other.getOrientation());
		setName(other.name);
		this.blocks = new CompoundInterval(other.blocks);
		setScore(other.getScore());
	}
	
	public BasicAnnotation(Annotation other) {
		setReferenceName(other.getReferenceName());
		setOrientation(other.getOrientation());
		setName(other.getName());
		addBlocks(other);
		setScore(other.getScore());
	}
	
	public BasicAnnotation(Collection<? extends Annotation> blocks) {
		if (blocks.size() == 0) throw new IllegalArgumentException("cannot create empty BasicAnnotation");
		Iterator<? extends Annotation> itr = blocks.iterator();
		Annotation first = itr.next();
		setReferenceName(first.getReferenceName());
		setOrientation(first.getOrientation());
		score = 0;
		setName(first.getName());
		addBlocks(first);
		while (itr.hasNext()) addBlocks(itr.next());
	}
	
	public BasicAnnotation(Collection<? extends Annotation> blocks, Strand orientation, String name) {
		this(blocks);
		setName(name);
		setOrientation(orientation);
	}
	
	public BasicAnnotation(Collection<? extends Annotation> blocks, String name) {
		this(blocks);
		setName(name);
	}
	
	public BasicAnnotation(String referenceName, CompoundInterval blocks, Strand orientation, String name) {
		setReferenceName(referenceName);
		this.blocks = blocks;
		setOrientation(orientation);
		score = 0;
		setName(name);
	}
	
	public BasicAnnotation(String referenceName, CompoundInterval blocks, Strand orientation) {
		this(referenceName, blocks, orientation, "");
	}
	
	public BasicAnnotation(String referenceName, CompoundInterval blocks) {
		this(referenceName, blocks, Strand.UNKNOWN);
	}
	

	public Annotation copy() {
		return new BasicAnnotation(this);
	}
	
	public static Annotation createFromUCSC(String ucsc) {
		String [] firstSplit = ucsc.split(":");
		String [] secondSplit = firstSplit[1].split("-");
		return new BasicAnnotation(firstSplit[0], Integer.valueOf(secondSplit[0]), Integer.valueOf(secondSplit[1]));
	}
	
	/**
	 * Create from string produced by AbstractAnnotation.getFullInfoString();
	 * @param fullInfoString Full info string
	 * @return Corresponding annotation object
	 */
	public static BasicAnnotation fromFullInfoString(String fullInfoString) {
		StringParser s1 = new StringParser();
		StringParser s2 = new StringParser();
		s1.parse(fullInfoString, "_");
		String name = s1.asString(0);
		String chr = s1.asString(1);
		Strand orientation = Strand.fromString(s1.asString(2));
		String block1 = s1.asString(3);
		s2.parse(block1, "-");
		int start1 = s2.asInt(0);
		int end1 = s2.asInt(1);
		BasicAnnotation rtrn = new BasicAnnotation(chr, start1, end1, orientation);
		rtrn.setName(name);
		for(int i = 5; i < s1.getFieldCount(); i++) {
			String block = s1.asString(i);
			s2.parse(block, "-");
			int start = s2.asInt(0);
			int end = s2.asInt(1);
			BasicAnnotation bl = new BasicAnnotation(chr, start, end, orientation);
			rtrn.addBlocks(bl);
		}
		return rtrn;
	}
	
	public static class Factory implements TabbedReader.Factory<BasicAnnotation> {
		@Override
		public BasicAnnotation create(String[] rawFields) throws ParseException {
			// Default is to read in BED foramt
			
			if (rawFields.length < 3) {
				throw new IllegalArgumentException("Cannot create BasicAnnotation from less than 3 fields");
			}
			
			String chr = rawFields[0];
			int start = Integer.parseInt(rawFields[1]);
			int end = Integer.parseInt(rawFields[2]);
			BasicAnnotation a = null;
			
			if (rawFields.length > 9) {
				// Assemble the annotation from the block information
				int nBlocks = Integer.parseInt(rawFields[9]);
				if (nBlocks > 0) {
					rawFields[10] = rawFields[10].replaceAll("\"", "");
					rawFields[11] = rawFields[11].replaceAll("\"", "");
					rawFields[10] = rawFields[10].replaceAll(" ", "");
					rawFields[11] = rawFields[11].replaceAll(" ", "");
				
					String [] sizes = rawFields[10].split(",");
					String [] starts   = rawFields[11].split(",");
					if (sizes.length < nBlocks || starts.length < nBlocks) {
						throw new ParseException("BAD BED FORMAT apparently the number of start ("+rawFields[10] +") and end ("+rawFields[11] +
							") items does not agree with te blockCount " + rawFields[9]);
					}
					
					for (int i = 0; i < nBlocks; i++) {
						int blockSize = Integer.parseInt(sizes[i]);
						int blockStart = Integer.parseInt(starts[i]);
						if (i == 0) {
							a = new BasicAnnotation(chr, start, start + blockSize);
						} else {
							a.addBlocks(new BasicAnnotation(chr, start + blockStart, start + blockStart + blockSize));
						}
					}
					if (a.getEnd() != end) throw new IllegalArgumentException("End specified by blocks does not match BED end");
				}
			}
			
			if (a == null) {
				// If we have no block information or block parsing failed, create the annotation from the start and end coordinates
				a = new BasicAnnotation(chr, start, end);
			}
			
			if (rawFields.length > 3) a.setName(rawFields[3]);
			if (rawFields.length > 4) a.setScore(Double.parseDouble(rawFields[4]));
			if (rawFields.length > 5) a.setOrientation(Annotation.Strand.fromString(rawFields[5]));
			
			// don't care about thickStart, thickEnd, rgb (rawFields[6], [7], and [8])
			
			return a;
		}
	}
	
	/********************************************************************************
	 * QUERY METHODS
	 */
	public int getStart() {
		return blocks.getStart();
	}
	
	public int getEnd() {
		return blocks.getEnd();
	}

	public int getMidpoint(){
		//get midpoint
		int mid=length()/2;
		//convert to reference space
		return getReferenceCoordinateAtPosition(mid);
	}
	
	public String getReferenceName() {
		return referenceName;
	}
	
	public String getName() {
		return name;  // even if it's null! if we override this when null it makes things complicated
	}
	
	public Strand getOrientation(){
		return orientation;
	}
	
	public boolean isUnoriented() {
		return !Strand.NEGATIVE.equals(getOrientation()) && !Strand.POSITIVE.equals(getOrientation()) ;
	}
	
	/**
	 * Returns unoriented blocks for this gene
	 * 
	 */
	public List<? extends Annotation> getBlocks() {
		return getBlocks(false);
	}
		
	public List<? extends Annotation> getBlocks(boolean oriented) {
		List<Annotation> list = new LinkedList<Annotation>();
		for (SingleInterval block : blocks.getBlocks()) {
			list.add(new BasicAnnotation(referenceName, block.getStart(), block.getEnd(), orientation));
		}
		if (oriented && orientation == Strand.NEGATIVE) {
			Collections.reverse(list);
		}
		return list;
	}

	/**
	 * This function returns the blocks on either ends of a given splicejunction
	 * @param spliceJunction
	 * @return
	 */
	public Annotation[] getFlankingBlocks(Annotation spliceJunction) {
		Annotation[] list = new Annotation[2];
		//List<? extends Annotation> blks = getBlocks();
		for (SingleInterval block : blocks.getBlocks()) {
			if(block.getEnd()==spliceJunction.getStart())
				list[0] = new BasicAnnotation(referenceName, block.getStart(), block.getEnd(), orientation);
			else if(block.getStart()==spliceJunction.getEnd()){
				list[1] = new BasicAnnotation(referenceName, block.getStart(), block.getEnd(), orientation);
			}
		}
		return list;
	}
	
	public int numBlocks() {
		return blocks.numBlocks();
	}
	
	public int length() {
		return blocks.length();
	}
	
	public int getReferenceCoordinateAtPosition(int positionInAnnotation, boolean ignoreOrientation) {
		if (getStrand() == Strand.NEGATIVE && !ignoreOrientation) {
			positionInAnnotation = length() - positionInAnnotation; //TODO who took out the -1? 
		}
		return blocks.getCoordinateAtPosition(positionInAnnotation);
	}
	
	public int getPositionAtReferenceCoordinate(int referenceCoordinate, boolean ignoreOrientation) {
		int positionInAnnotation = blocks.getPositionAtCoordinate(referenceCoordinate);
		if (getStrand() == Strand.NEGATIVE && !ignoreOrientation) {
			positionInAnnotation = length() - positionInAnnotation; //TODO who took out the -1?  
		}
		return positionInAnnotation;
	}
	
	@Override
	public double getScore() { 
		return score;
	}
	
	public String toString() {
		return toBED();
	}


	/********************************************************************************
	 * SETTING METHODS
	 */
	
	public void setStart(int start) {
		blocks.setStart(start);
	}
	
	public void setEnd(int end) {
		blocks.setEnd(end);
	}

	
	public void setOrientation(Strand orientation) {
		if(orientation==null){orientation=Strand.UNKNOWN;}
		this.orientation = orientation;
	}
	
	public void setOrientedStart(int orientedStart) {
		// TODO is this the correct 0-based behavior?
		if (orientation == Strand.NEGATIVE) {
			setEnd(orientedStart);
		} else {
			setStart(orientedStart);
		}
	}
	
	public void setOrientedEnd(int orientedEnd) {
		if (orientation == Strand.NEGATIVE) {
			setStart(orientedEnd);
		} else {
			setEnd(orientedEnd);
		}
	}
	
	public void setReferenceName(String refName) {
		this.referenceName = refName != null ? refName.intern() : null;
	}
	
	public void setName(String n) {
		this.name =  n;
	}
	
	public void setScore(double s) {
		score = s;
	}
	
	public void addBlocks(Annotation block) {
		if (numBlocks() > 0 && !block.getReferenceName().equals(getReferenceName())) {
			throw new IllegalArgumentException("Reference names must match: "+block.getReferenceName()+" "+getReferenceName()+" "+getName());
		}
		// strand of the new block is ignored
		
		// recursively add blocks if the Annotation has multiple blocks
		if (block.numBlocks() > 1) {
			for (Annotation a : block.getBlocks()) {
				addBlocks(a);
			}
		} else {
			if (numBlocks() == 0) {
				setReferenceName(block.getReferenceName());
			}
			blocks.addInterval(block.getStart(), block.getEnd());
		}
	}

	public void addBlocks(Collection<? extends Annotation> blocks) {
		for(Annotation block: blocks){
			addBlocks(block);			
		}
	}
	
	public void shift(int delta) {
		blocks.shift(delta);
	}
	
	public void moveToCoordinate(int coordinateInReference) {
		blocks.moveToCoordinate(coordinateInReference);
	}
	
	/********************************************************************************
	 * COMPARISON METHODS
	 */
	
	@Override
	public boolean overlaps(Annotation other, int buffer, boolean considerOrientation) {
		if (considerOrientation && (getOrientation() != other.getOrientation())) return false;
		if (!getReferenceName().equalsIgnoreCase(other.getReferenceName())) return false;
		return overlaps(other.getBlocks(),buffer);
		/*if (BasicAnnotation.class.isInstance(other)) {
			return overlaps((BasicAnnotation) other, buffer, considerOrientation);
		} else {
			BasicAnnotation basic = new BasicAnnotation(other);
			return overlaps(basic, buffer, considerOrientation);
		}*/
	}
	
	public boolean overlaps(List<? extends Annotation> otherBlocks,int buffer){
		//For each of this block
		for(Annotation block:getBlocks()){
			//For each of the other's blocks
			SingleInterval X = new SingleInterval(block.getStart(),block.getEnd());
			
			for(Annotation other: otherBlocks){
				//getBlocks() is already returning single intervals
				SingleInterval Y = new SingleInterval(other.getStart(),other.getEnd());
				//If this block overlaps the other block
				if(X.overlaps(Y))
					return true;
			}
		}
		return false;
	}
	
	@Override
	public int getOverlap(Annotation other) {
		int overlap = 0;
		if (getReferenceName().equals(other.getReferenceName())) {
			BasicAnnotation basic = new BasicAnnotation(other);
			overlap = blocks.intersect(basic.blocks).length();
		}
		return overlap; 
	}
	
	public boolean contains(Annotation other) {
		boolean result = false;
		if (getReferenceName().equals(other.getReferenceName())) {
			BasicAnnotation basic = new BasicAnnotation(other);
			result = blocks.contains(basic.blocks);
		}
		return result;
	}
	
	public Annotation union(Annotation other) {
		if(!getReferenceName().equals(other.getReferenceName())) {
			throw new IllegalArgumentException("Cannot merge annotations on different reference sequences.");
		}
		BasicAnnotation basic = new BasicAnnotation(other);
		CompoundInterval newBlocks = blocks.union(basic.blocks);
		Annotation result = new BasicAnnotation(getReferenceName(), newBlocks, getStrand());
		if (getOrientation() == other.getOrientation()) result.setOrientation(getOrientation());
		return result;
	}
	
	public Annotation intersect(Annotation other) {
		Annotation result = null;
		if (getReferenceName().equals(other.getReferenceName())) {
			BasicAnnotation basic = new BasicAnnotation(other);
			CompoundInterval newBlocks = blocks.intersect(basic.blocks);
			
			// JE 1/17/13 changed to return null if empty
			if (newBlocks.length() > 0) result = new BasicAnnotation(getReferenceName(), newBlocks, getStrand());
		}
		return result;
	}


	@Override
	public List<Annotation> disect(Annotation a) {
		/*List<Annotation> rtrn=new ArrayList<Annotation>();
		
		//We will split the current annotation by the block of the other
		List<? extends Annotation> blocks=a.getBlocks();
		for(Annotation block: blocks){
			if(getStart()<block.getStart()){
				rtrn.addAll(split(getStart(), block.getStart()).getBlocks());
			}
			if(block.getStart()<block.getEnd()){
				rtrn.addAll(split(block.getStart(), block.getEnd()).getBlocks());
			}
			if(block.getEnd()<getEnd()){
				rtrn.addAll(split(block.getEnd(), getEnd()).getBlocks());
			}
		}
		
		return rtrn;*/
		
		throw new UnsupportedOperationException("TODO");
	}


	@Override
	public List<Annotation> disect(List<? extends Annotation> disectors) {
		List<Annotation> rtrn=new ArrayList<Annotation>();
		
		for(Annotation annotation: disectors){
			rtrn.addAll(disect(annotation));
		}
		return rtrn;
	}


	@Override
	public Annotation minus(Annotation other) {
		BasicAnnotation basic=new BasicAnnotation(other);
		CompoundInterval newBlocks =this.blocks.minus(basic.blocks);
		BasicAnnotation rtrn= new BasicAnnotation(getReferenceName(), newBlocks, getStrand());
		//List<Annotation> list=new ArrayList<Annotation>();
		//list.add(rtrn);
		return rtrn;
		
		//throw new UnsupportedOperationException("TODO");
	}


	@Override
	public Annotation minus(Collection<? extends Annotation> others) {
		BasicAnnotation basic=new BasicAnnotation(others);
		return minus(basic);
	}


	@Override
	public void stitchTo(Annotation next) {
		setStart(Math.min(getStart(), next.getStart()));
		setEnd(Math.max(getEnd(), next.getEnd()));
		//throw new UnsupportedOperationException("TODO");
	}

	@Override
	public Annotation complement() {
		CompoundInterval interval=this.blocks.complement();
		BasicAnnotation rtrn=new BasicAnnotation(this.referenceName, interval, this.orientation, this.name);
		return rtrn;
	}

	@Override
	/**
	 * The BasicAnnotation will return the gaps themselves
	 */
	public Collection<? extends Annotation> getSpliceConnections() {
		Collection<Annotation> rtrn=new TreeSet<Annotation>();
		//if has exons
		if(this.getBlocks().size()>1){
			rtrn.addAll(this.complement().getBlocks());
		}
		return rtrn;
	}
	
	
	/* EXTRACTED FROM BROAD BasicLightweightAnnotation */
	public static List<Annotation> stitchList(Collection<? extends Annotation> sortedList, int maxDistanceToStitch) {

		Stack<Annotation> result = new Stack<Annotation>();

		if(sortedList.size() == 0 ) {
			return result;
		}
		
		Iterator<? extends Annotation> it = sortedList.iterator();
		result.push(it.next());
		while(it.hasNext()) {
			Annotation next = it.next();
			Annotation curr = result.pop();
			if(curr.overlaps(next,maxDistanceToStitch)) {
				curr.stitchTo(next);
				curr.setName(curr.getName()+"-"+next.getName());
				result.push(curr);
			} else {
				result.push(curr);
				result.push(next);
			}
		}	
		return result;
	}

	public static void main(String[] args)throws IOException{
		
		Collection<BasicAnnotation> blocks = new ArrayList<BasicAnnotation>();
		blocks.add(new BasicAnnotation("chr1",1000,2000,Strand.POSITIVE));
		blocks.add(new BasicAnnotation("chr1",3000,4000,Strand.POSITIVE));
		blocks.add(new BasicAnnotation("chr1",5000,6000,Strand.POSITIVE));
		BasicAnnotation a = new BasicAnnotation(blocks);
		
		Collection<BasicAnnotation> blocks1 = new ArrayList<BasicAnnotation>();
		blocks1.add(new BasicAnnotation("chr1",1000,2000,Strand.POSITIVE));
		blocks1.add(new BasicAnnotation("chr1",3000,4500,Strand.POSITIVE));
		blocks1.add(new BasicAnnotation("chr1",5500,7000,Strand.POSITIVE));
		BasicAnnotation b = new BasicAnnotation(blocks1);
		
		System.out.println(a.toBED());
		System.out.println(b.toBED());
		
		for(Annotation exon1: b.getBlocks()){
			for(Annotation intron2:a.getSpliceConnections()){
				if(exon1.overlaps(intron2)){
					Annotation p = a.minus(exon1);
					Annotation q = b.minus(exon1);
					System.out.println(p.toBED());
					System.out.println(q.toBED());
					//System.out.println(BuildScriptureCoordinateSpace.compatible(p,q));
				}
			}
		}
		Annotation c = b.intersect(a);
		System.out.println(c.toBED());
		//System.out.println(BuildScriptureCoordinateSpace.compatible(b,c));
	}
	
	
	/********************************************************************************
	 * PRIVATE UTILITY METHODS
	 */
	
}
