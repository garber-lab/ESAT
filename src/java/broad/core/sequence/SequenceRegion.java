package broad.core.sequence;

import java.util.Collection;
import java.util.List;

import nextgen.core.annotation.Annotation;
import broad.core.annotation.BasicGenomicAnnotation;
import broad.core.annotation.GenomicAnnotation;
import broad.core.annotation.LightweightGenomicAnnotation;
import broad.core.annotation.TwoSubjectAnnotation;

public class SequenceRegion extends Sequence implements Annotation{
	BasicGenomicAnnotation annotation;
	private String containingSequenceId;
	public static final int INF = Integer.MAX_VALUE;
	
	public SequenceRegion(String containingSequenceId) {
		super(null);
		this.containingSequenceId = containingSequenceId;
		annotation = new BasicGenomicAnnotation(containingSequenceId);
		annotation.setSequence(this);
		annotation.setEnd(INF);
	}
	
	public SequenceRegion(String containingSequenceId, String chr, int start, int end) {
		super(null);
		this.containingSequenceId = containingSequenceId;
		annotation = new BasicGenomicAnnotation(containingSequenceId, chr, start, end);
		annotation.setSequence(this);
	}
	
	public SequenceRegion(String containingSequenceId, Annotation annotation) {
		super(annotation.getName());
		this.containingSequenceId = containingSequenceId;
		this.annotation = new BasicGenomicAnnotation(annotation);
		
	}
	public int getRegionEnd() {
		return annotation.getEnd() == INF && super.getSequenceBases() != null && getSequenceBases().length() > 0 
			? getSequenceBases().length() + annotation.getStart() - 1
			: annotation.getEnd();
	}
	public void setRegionEnd(int regionEnd) {
		annotation.setEnd(regionEnd);
	}
	public int getRegionStart() {
		return annotation.getStart();
	}
	public void setRegionStart(int regionStart) {
		annotation.setStart(regionStart);
	}
	
	public WindowSlider getSlider(int windowSize, int overlap) {
		return WindowSlider.getSlider(this, windowSize, overlap);
	}
	
	public int length() {
		if(getEnd() == INF && getSequenceBases() == null && getSequenceBases().length() == 0) {
			return INF;
		} else if(getEnd() < INF) {
			return getEnd() - getStart();
		} else {
			return super.getLength();
		}
	}
	
	public Strand getOrientation() { return annotation.getOrientation(); }
	public void setReversedOrientation(boolean reversed) { annotation.setOrientation(!reversed); }
	
	public String getContainingSequenceId() {
		return containingSequenceId;
	}
	
	public String getId() {
		return super.getId() == null ? containingSequenceId + ":" + getStart() + "-" + getEnd() : super.getId();
	}
	
	public int getStart() {
		return annotation.getStart();
	}
	
	public int getEnd() {
		return annotation.getEnd();
	}
	
	public double getScore() {
		return annotation.getScore();
	}
	
	public void setScore(double score) {
		annotation.setScore(score);
	}
	
	public double getExtraScore(int i ) {
		return annotation.getExtraScore(i);
	}
	
	
	public String getName() {
		return getId();
	}
	
	public void setId(String id) {
		super.setId(id);
		annotation.setName(id);
	}
	
	public void setName(String name) {
		setId(name);
	}
	
	public boolean inReversedOrientation() {
		return annotation.inReversedOrientation();
	}
	
	public long getMiddle() {
		return annotation.getMiddle();
	}
	
	public void setStart(int start) {
		setRegionStart(start);
	}
	
	public void setEnd(int end) {
		setRegionEnd(end);
	}
	public Sequence getSequence() {
		return this;
	}
	public void setSequence(Sequence seq) {
		super.setSequenceBases(seq.getSequenceBases());
	}
	
	public String getChromosome() {
		return annotation.getChromosome();
	}
	
	public void setChromosome(String chr) {
		annotation.setChromosome(chr);
	}
	public boolean overlaps(LightweightGenomicAnnotation other, int buffer) {
		return annotation.overlaps(other, buffer);
	}
	public boolean overlaps(LightweightGenomicAnnotation other) {
		return annotation.overlaps(other);
	}
	
	public int getOverlap(LightweightGenomicAnnotation other) { return annotation.getOverlap(other);}
	
	public Annotation minus(Annotation other) {
		return annotation.minus(other);
	}
	
	public Annotation minus(List<? extends Annotation> others) {
		return annotation.minus(others);
	}
	
	public void takeIntersection(LightweightGenomicAnnotation other) {
		annotation.takeIntersection(other);
	}
	
	public void takeUnion(LightweightGenomicAnnotation other) {
		annotation.takeUnion(other);
	}
	public void stitchTo(LightweightGenomicAnnotation other) {
		annotation.stitchTo(other);
	}
	
	public boolean isFlankedBy(TwoSubjectAnnotation twoSubjectAnnotation, int buffer) {
		return annotation.isFlankedBy(twoSubjectAnnotation, buffer);
	}
	public int getFivePrimeBases() {
		return annotation.getFivePrimeBases();
	}
	public int getThreePrimeBases() {
		return annotation.getThreePrimeBases();
	}
	
	public String toString() {
		return getContainingSequenceId() + ":" + getRegionStart() + "-" + getRegionEnd();
	}
	
	public int getAbsoluteStart(int absoluteStart){
		return getRegionStart()+absoluteStart;
	}
	
	public int getAbsoluteEnd(int absoluteStart){
		return getRegionEnd()+absoluteStart;
	}
	
	public String getLocationString() { return annotation.getLocationString();}
	public String getChromosomeString() { return annotation.getChromosomeString(); }
	
	public SequenceRegion extractRegionBasedOnGC(float targetGC, int size, int buffer) {
		SequenceRegion theRegion = super.extractRegionBasedOnGC(targetGC, size, buffer);
		if(theRegion != null) {
			theRegion.setRegionStart(getRegionStart() + theRegion.getRegionStart());
			theRegion.setRegionEnd(getRegionStart() + theRegion.getRegionEnd());
			theRegion.setChromosome(annotation.getChromosome());
		}
		return theRegion;
	}
	public boolean contains(LightweightGenomicAnnotation other) {
		return annotation.contains(other);
	}
	public int getDistanceTo(LightweightGenomicAnnotation other) {
		return annotation.getDistanceTo(other);
	}
	public void setOrientation(String orientation) {
		annotation.setOrientation(orientation);
		
	}
	public int compareTo(GenomicAnnotation arg0) {
		return annotation.compareTo(arg0);
	}
	public List<Annotation> disect(Annotation a) {
		return annotation.disect(a);
	}
	public List<Annotation> disect(List<? extends Annotation> disectors) {
		return annotation.disect(disectors);
	}
	
	public SequenceRegion getRegion(int start, int end) {
		SequenceRegion region = super.getRegion(start, end);
		if(annotation.getChromosome() != null) {
			region.setChromosome(annotation.getChromosome());
		}
		return region;
	}
	public int getOrientedEnd() {
		return annotation.getOrientedEnd();
	}
	public int getOrientedStart() {
		return annotation.getOrientedStart();
	}
	
	public boolean isUnoriented() {
		return annotation.isUnoriented();
	}
	
	public void addBlock(String name, int start, int end) {
		//DO nothing as there is nothing to do here.
		
	}
	public List<? extends Annotation> getBlocks() {
		return annotation.getBlocks();
	}
	
	public List<? extends Annotation> getBlocks(boolean oriented) {
		return annotation.getBlocks(oriented);
	}
	
	public boolean mayHaveBlocks() {
		return false;
	}
	
	public String toUCSC() {
		return annotation.toUCSC();
	}
	
	public int getLength() {
		return length();
	}


	public boolean overlaps(Collection<? extends Annotation> others, int buffer) {
		return annotation.overlaps(others, buffer);
	}


	public boolean overlaps(Collection<? extends Annotation> others) {
		return annotation.overlaps(others);
	}

	
	public void removeExtraScores() {
		annotation.removeExtraScores();
	}

	public int compareTo(LightweightGenomicAnnotation b) {
		return annotation.compareTo(b);
	}

	public LightweightGenomicAnnotation intersect(LightweightGenomicAnnotation other) {

		return annotation.intersect(other);
	}

	public List<Annotation> intersect(List<? extends Annotation> annotations) { 
		return annotation.intersect(annotations); 
	}
	public double percentGC() {
		char[] seqChar=this.getSequenceBases().toCharArray();
		double GCCount=0;
		for(int i=0; i<seqChar.length; i++){
			if(seqChar[i]=='G' || seqChar[i]=='C' || seqChar[i]=='g' || seqChar[i]=='c'){GCCount++;}
		}
		double percentGC=GCCount/seqChar.length;
		return percentGC;
	}

	public void addExtraScore(double score) {
		annotation.addExtraScore(score);
		
	}

	public List<Double> getExtraScores() {
		return annotation.getExtraScores();
	}
	
	/**
	 * Checks whether two annotations differ by a small (fudge) factor
	 * @param fudge Maximum difference at either end or start to consider similar
	 * @return
	 */
	public boolean almostEqual(LightweightGenomicAnnotation other, int fudge) {
		return annotation.almostEqual(other, fudge);
	}
	
	@Override
	public int getSAMStart() {
		return annotation.getSAMStart();
	}
	@Override
	public int getSAMEnd() {
		return annotation.getSAMEnd();
	}
	@Override
	public String getReferenceName() {
		return annotation.getReferenceName();
	}
	@Override
	public String getChr() {
		return annotation.getChr();
	}
	@Override
	public Strand getStrand() {
		return annotation.getStrand();
	}
	@Override
	public boolean hasOrientation() {
		return annotation.hasOrientation();
	}
	@Override
	public int numBlocks() {
		return annotation.numBlocks();
	}
	@Override
	public int getSize() {
		return annotation.getSize();
	}
	@Override
	public int size() {
		return annotation.size();
	}
	@Override
	public int getLengthOnReference() {
		return annotation.getLengthOnReference();
	}
	@Override
	public int getReferenceCoordinateAtPosition(int positionInAnnotation) {
		return annotation.getReferenceCoordinateAtPosition(positionInAnnotation);
	}
	@Override
	public int getPositionAtReferenceCoordinate(int referenceCoordinate) {
		return annotation.getPositionAtReferenceCoordinate(referenceCoordinate);
	}
	@Override
	public int getReferenceCoordinateAtPosition(int positionInAnnotation,
			boolean ignoreOrientation) {
		return annotation.getReferenceCoordinateAtPosition(positionInAnnotation, ignoreOrientation);
	}
	@Override
	public int getPositionAtReferenceCoordinate(int referenceCoordinate,
			boolean ignoreOrientation) {
		return annotation.getPositionAtReferenceCoordinate(referenceCoordinate, ignoreOrientation);
	}
	@Override
	public void setOrientation(char orientation) {
		annotation.setOrientation(orientation);
		
	}
	@Override
	public void setOrientation(Strand orientation) {
		annotation.setOrientation(orientation);
		
	}
	@Override
	public void setOrientedStart(int orientedStart) {
		annotation.setOrientedStart(orientedStart);
	}
	@Override
	public void setOrientedEnd(int orientedEnd) {
		annotation.setOrientedEnd(orientedEnd);
	}
	@Override
	public void setReferenceName(String refName) {
		annotation.setReferenceName(refName);
	}
	@Override
	public boolean equals(Annotation other) {
		return annotation.equals(other);
	}
	@Override
	public void expand(int deltaStart, int deltaEnd) {
		annotation.expand(deltaStart, deltaEnd);
	}
	@Override
	public Annotation trim(int deltaStart, int deltaEnd) {
		return annotation.trim(deltaStart, deltaEnd);
	}
	@Override
	public Annotation copy() {
		return annotation.copy();
	}
	@Override
	public int getDistanceTo(Annotation other) {
		return annotation.getDistanceTo(other);
	}
	@Override
	public String toBED() {
		return annotation.toBED();
	}
	@Override
	public String toShortBED() {
		return annotation.toShortBED();
	}
	@Override
	public String toBEDGraph() {
		return annotation.toBEDGraph();
	}
	@Override
	public boolean overlaps(Annotation other, int buffer) {
		return annotation.overlaps(other, buffer);
	}
	@Override
	public boolean overlaps(Annotation other) {
		return annotation.overlaps(other);
	}
	@Override
	public boolean overlaps(Annotation other, int buffer, boolean considerStrand) {
		return annotation.overlaps(other, buffer, considerStrand);
	}
	@Override
	public boolean overlapsStranded(Annotation other) {
		return annotation.overlapsStranded(other);
	}
	@Override
	public int getOverlap(Annotation other) {
		return annotation.getOverlap(other);
	}
	@Override
	public boolean contains(Annotation other) {
		return annotation.contains(other);
	}
	@Override
	public Annotation union(Annotation other) {
		return annotation.union(other);
	}
	@Override
	public Annotation intersect(Annotation other) {
		return annotation.intersect(other);
	}
	@Override
	public int compareToAnnotation(Annotation b) {
		return annotation.compareToAnnotation(b);
	}
	@Override
	public void stitchTo(Annotation next) {
		annotation.stitchTo(next);
	}
	@Override
	public String toCigar() {
		return annotation.toCigar();
	}
	@Override
	public boolean fullyContains(Annotation annotation) {
		return annotation.fullyContains(annotation);
	}
	@Override
	public int compareTo(Annotation arg0) {
		return annotation.compareTo(arg0);
	}
	@Override
	public void moveToCoordinate(int i) {
		throw new UnsupportedOperationException();
	}
	@Override
	public void shift(int delta) {
		throw new UnsupportedOperationException();
	}
	@Override
	public boolean isNegativeStrand() {
		return annotation.isNegativeStrand();
	}

	@Override
	public Annotation complement() {
		return annotation.complement();
	}

	@Override
	public Annotation minus(Collection<? extends Annotation> others) {
		return annotation.minus(others);
	}

	@Override
	public boolean equals(Annotation other, boolean useOrientation) {
		return annotation.equals(other, useOrientation);
	}

	@Override
	public Collection<? extends Annotation> getSpliceConnections() {
		return annotation.getSpliceConnections();
	}

	@Override
	public String toBED(int r, int g, int b) {
		throw new UnsupportedOperationException("TODO");
	}

	@Override
	public boolean overlaps(Annotation other, boolean considerOrientation) {
		throw new UnsupportedOperationException("TODO");
	}

	@Override
	public int getMidpoint() {
		//get midpoint
		int mid=length()/2;
		//convert to reference space
		return getReferenceCoordinateAtPosition(mid);
	}

	@Override
	public String getFullInfoString() {
		// TODO Auto-generated method stub
		return null;
	}
}
