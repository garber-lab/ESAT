package umms.core.alignment;

import java.util.Collection;

import umms.core.coordinatesystem.CoordinateSpace;
import umms.core.feature.Window;
import umms.core.alignment.AbstractPairedEndAlignment.TranscriptionRead;
import umms.core.annotation.Annotation;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;

/**
 * Based on the Alignment interface in IGV. This interface represents any read or pair of reads aligned to the genome
 * @author skadri
 *
 */
public interface Alignment extends Annotation {

	//TODO Implement all of the SAM functions
	
	/**
	 * This method returns the name of the alignment
	 * @return
	 */
	public String getReadName();

    /**
     * Returns the start position of the alignment in genomic space.
     * @return
     */
    public int getFragmentStart();

	/**
	 * Returns the end position of the alignment in genomic space.
	 * @return
	 */
    public int getFragmentEnd();
	
    /**
     * Returns the fragment size (insert+lengths of the reads) 
     * @return
     */
    public Collection<Integer> getFragmentSize(CoordinateSpace C);

    /**
     * Returns the mapping quality
     * @return
     */
    public int getMappingQuality();

    /**
     * Returns true if the alignment is true.
     * @return
     */
    public boolean isPaired();

    /**
     * For paired end, returns true if the read in direction of transcription is negative stranded
     * @return
     */
    public boolean isNegativeStrand();
    
    /**
     * For single end, returns false.
     * For paired end, returns true if reads successfully map to different chromosomes.
     * @return
     */
    public boolean isChimera();
    
    /**
     * Returns true if the alignment is flagged as a duplicate
     * @return
     */
    public boolean isDuplicate();

    /**
     * Sets the isDuplicate Flag
     * @param duplicateFlag
     */
    public void setDuplicateFlag(boolean duplicateFlag);
    
    /**
     * Returns the strand for the fragment
     * @return
     */
    public Strand getFragmentStrand();
    
	/**
	 * Get the ending position of the read considering strand
	 * @return Lowest position if negative strand, highest position otherwise
	 */
	public int getLastFragmentPositionStranded();

	/**
	 * Get the ending position of the read considering strand
	 * @return Lowest position if negative strand, highest position otherwise
	 */
	public int getFirstFragmentPositionStranded();

	/**
	 * Get midpoint of fragment with respect to coordinate space
	 * @param annot Parent annotation
	 * @return Fragment midpoint in coordinate space
	 */
	public int getFragmentMidpoint(Annotation annot);
	
    /**
     * This method returns the fragment formed by the alignment in the coordinate space C, specified.
     * @param C
     * @return
     */
    public Collection<? extends Window> getFragment(CoordinateSpace C);
	
	
	 /**
     * This method returns the alignment blocks that the read actually aligns to.
     * @param C
     * @return
     */
	public Annotation getReadAlignmentBlocks(CoordinateSpace C);

	/**
	 * Get start of the fragment
	 * @return
	 */
	public int getStart();

	/**
	 * Get end of the fragment
	 * @return
	 */
	public int getEnd();

	/**
	 * Get chromosome
	 * @return
	 */
	public String getChr();

	/**
	 * Returns the start position of the alignment
	 * @return
	 */
	public int getAlignmentStart();

	/**
	 * Returns the end position of the alignment
	 * @return
	 */
	public int getAlignmentEnd();

	/**
	 * Returns true if the alignment is mapped
	 * @return
	 */
	public boolean isMapped();

	/**
	 * Returns the value for the specified attribute
	 * @param string
	 * @return
	 */
	public Object getAttribute(String string);

	public String getReadSequence();

	/**
	 * The weight to give this alignment in the scoring.
	 * For example, if you want to scale by number of hits, the weight would be 1/NH
	 * @return
	 */
	public double getWeight();

	/**
	 * Returns whether the alignment represented a proper pair defined by the aligner
	 * @return
	 */
	public boolean isProperPair();

	/**
	 * Set whether the alignment represented a proper pair defined by the aligner
	 * @return
	 */
	public void setProperPairFlag(boolean properPairFlag);
	
	public SAMRecord toSAMRecord();

	/**
	 * This will return the actual read alignments
	 * @param space Coordinate Space to return the blocks in
	 */
	public Collection<Annotation> getReadAlignments(CoordinateSpace space);
	
	/**
	 * This will return the actual read alignment objects NOT BLOCKS
	 */
	public Collection<Alignment> getReadMates();

	/**
	 * Return the splice connections contained within this read
	 * @return
	 */
	public Collection<? extends Annotation> getSpliceConnections();
	
	/**
	 * Returns whether the read contains an insertion or deletion
	 * @return
	 */
	public boolean hasIndel();
	
	public void setFragmentStrand(TranscriptionRead strand);
	
	/**
	 * Get the coordinates of the first and last position of the interval between the reads
	 * @return The coordinates of the space between the reads or null if single read or reads overlap
	 */
	public int[] getIntervalBetweenReads();
	
	
	/**
	 * Sets the SAM File header for this alignment
	 * @param header
	 */
	public void setHeader(SAMFileHeader header);
	
	/**
	 * Returns the SAM file header for this alignment
	 * @return
	 */
	public SAMFileHeader getHeader();
	
	
}
