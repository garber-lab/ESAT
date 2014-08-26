package broad.core.annotation;

import java.util.List;

import nextgen.core.annotation.Annotation;

import broad.core.sequence.Sequence;

public interface GenomicAnnotation extends Cloneable ,Feature, LightweightGenomicAnnotation {
	//String getDisplayName();
	void setReversedOrientation(boolean isInReversedOrientation);
	
	/**
	 * @deprecated Use length() instead
	 * @return
	 */
	int getLength();
	
	Sequence getSequence();
	void setSequence(Sequence seq);
	public void setOrientation(String orientationString);
	public Strand getOrientation();
	
	/**
	 * Genomic annotation integer identification if such exists.
	 */
	public String getId();
	public void setId(String id);
	
	/**
	 * @return int - the five prime buffer (i.e. bases not belonging to the annotation proper).
	 */
	public int getFivePrimeBases();

	/**
	 * @return int - the three prime buffer (i.e. bases not belonging to the annotation proper).
	 */
	public int getThreePrimeBases();
	
	/**
	 * Checks whether a twoSubjectAnnotation flanks the instance
	 */
	boolean isFlankedBy(TwoSubjectAnnotation twoSubjectAnnotation, int buffer);
	
	
	
	
	
	
	
	/**
	 * Gets the start of the annotation considering its orientation 
	 * @return getStart() if the feature is in direct orientation or getEnd() otherwise
	 */
	int getOrientedStart();
	
	/**
	 * Gets the end of the annotation considering its orientation 
	 * @return getEnd() if the feature is in direct orientation or getStart() otherwise
	 */
	int getOrientedEnd();
	
	/**
	 * Returns a string of the form chr<chromosome>:start-end
	 * @return
	 */
	public String getLocationString();
	
	/**
	 * Returns true if the annotation (like a full BED or a RefSeq) has sub annotations like exons or blocks
	 */
	public boolean mayHaveBlocks();
	
	/**
	 * Returns a list of blocks if the annotations has any (@see containsBlocks)
	 */
	public List<? extends Annotation> getBlocks();
	
	/**
	 * Adds a block to the annotation if it supports blocks
	 */
	public void addBlock(String name, int start, int end);
	
}
