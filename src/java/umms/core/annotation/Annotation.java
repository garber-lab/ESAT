package umms.core.annotation;

import java.util.Collection;
import java.util.List;

public interface Annotation extends Comparable<Annotation> {

	/**
	 * @author engreitz
	 * Use this enumeration whenever possible to represent strand (rather than using a string or character)
	 */
	public enum Strand {
		POSITIVE('+'), NEGATIVE('-'), UNKNOWN('*');
		private char value;
		
		private Strand(char value) {
			this.value = value;
		}
		
		private Strand(String value) {
			if (value.length() > 1) throw new IllegalArgumentException("Illegal strand string");
			this.value = value.charAt(0);
		}

		public String toString() {
			return "" + value;
		}
		
		public Strand getReverseStrand() {
			if(this.equals(POSITIVE)){
				return NEGATIVE;
			}
			else if(this.equals(NEGATIVE)){
				return POSITIVE;
			}
			else{
				return UNKNOWN;
			}
		}
		
		public static Strand fromString(String value) {
			if (value.equals("+")) return POSITIVE;
			if (value.equals("-")) return NEGATIVE;
			else return UNKNOWN;
		}
	}
	
	/**
	 * @return start of the annotation (inclusive)
	 */
	public int getStart();
	
	
	/**
	 * Returns the start position of this annotation into SAM coordinate space
	 * SAM coordinates are 1-based and inclusive whereas all of our objects are 0-based exclusive
	 * @return
	 */
	public int getSAMStart();
	
	/**
	 * @return end of the annotation (exclusive)
	 */
	public int getEnd();
	
	/**
	 * Returns the start position of this window into SAM coordinate space
	 * SAM coordinates are 1-based and inclusive whereas all of our objects are 0-based exclusive
	 * @return
	 */
	public int getSAMEnd();
	
	/**
	 * @return name of the reference sequence
	 */
	public String getReferenceName();
	
	
	/**
	 * @return name of the reference sequence
	 */
	public String getChr();
	
	/**
	 * @return the name of the annotation
	 */
	public String getName();
	
	/**
	 * @return orientation or strand
	 */
	public Strand getOrientation();
	
	/**
	 * @return orientation or strand
	 */
	public Strand getStrand();
	
	/**
	 * @return true if strand is known
	 */
	public boolean hasOrientation();
	
	public boolean isNegativeStrand();
	
	/**
	 * @return number of blocks in the Annotation
	 */
	public int numBlocks();
	
	
	/**
	 * @return list of blocks in coordinate order
	 */
	public List<? extends Annotation> getBlocks();
	

	/**
	 * 
	 * @param oriented set to true if the list should be ordered according to the annotation orientation and by position
	 * @return
	 */
	public List<? extends Annotation> getBlocks(boolean oriented);
	
	/**
	 * @return length of the Annotation; if it has blocks, the gaps between blocks are not included
	 */
	public int length();
	
	/**
	 * @return same as size();
	 */
	public int getSize();
	
	/**
	 * @return same as length()
	 */
	public int size();
	
	/**
	 * @return The distance between the start of the first block and end of the last on the reference sequence (including spaces between blocks)
	 */
	public int getLengthOnReference();
	
	/**
	 * @return the 5-prime end of the annotation, considering its orientation or strand
	 */
	public int getOrientedStart();
	
	/**
	 * @return the 3-prime end of the annotation, considering its orientation or strand
	 */
	public int getOrientedEnd();
	
	public double getScore();
	
	/**
	 * @return false if the annotation is not in either negative or positive orientation
	 */
	public boolean isUnoriented();
	/**
	 * @param positionInAnnotation 0-based position in the Annotation
	 * @return 0-based reference coordinate at the specified position, considering orientation.  
	 * NOTE:  positionInAnnotation=0 for a negative-strand annotation will match getEnd()
	 */
	public int getReferenceCoordinateAtPosition(int positionInAnnotation);
	
	
	/**
	 * @param referenceCoordinate 0-based reference coordinate
	 * @return 0-based position in the Annotation, considering orientation
	 */
	public int getPositionAtReferenceCoordinate(int referenceCoordinate);
	public int getReferenceCoordinateAtPosition(int positionInAnnotation, boolean ignoreOrientation);
	public int getPositionAtReferenceCoordinate(int referenceCoordinate, boolean ignoreOrientation);
	
	
	/**
	 * Sets a new start position.  Will delete and/or truncate blocks if the new start position 
	 * is greater than the old one.  Will extend the first block if the new start position is less than
	 * old one.
	 * @param new start position (ignoring orientation) of the annotation
	 */
	public void setStart(int start);
	
	/**
	 * Sets a new end position.  Will delete and/or truncate blocks if the new end position 
	 * is less than the old one.  Will extend the last block if the new end position is greater than
	 * old one.
	 * @param new end position (ignoring orientation) of the annotation
	 */
	public void setEnd(int end);
	
	public void setOrientation(char orientation);
	public void setOrientation(Strand orientation);
	
	public void setOrientedStart(int orientedStart);
	public void setOrientedEnd(int orientedEnd);
	
	/**
	 * @param refName name of the reference sequence / chromosome
	 */
	public void setReferenceName(String refName);
	
	
	/**
	 * @param name name of the annotation
	 */
	public void setName(String name);
	
	public void setScore(double score);
	
	
	public boolean equals(Annotation other);

	/**
	 * Expand annotation by a set number of bases, taking into account blocks
	 * @param deltaStart  positive extends the gene.  
	 * @param deltaEnd  positive extends the gene.
	 */
	public void expand(int deltaStart, int deltaEnd);
	
	/**
	 * Trim annotation by a set number of bases, taking into account blocks. 
	 * @param deltaStart  (Positive) The number of bases to remove from the beginning of gene.  e.g. deltaStart=1 removes one base from the front
	 * @param deltaEnd  (Positive) The number of bases to remove from the end of gene.   e.g. deltaStart=1 removes one base fromt he end
	 * @return 
	 */
	public Annotation trim(int deltaStart, int deltaEnd);
	
	
	/**
	 * Shift the positioning of the annotation by a given number of base pairs.
	 * @param delta Positive or negative shift to apply to all blocks in the annotation. {@code delta=0} leaves the annotation unchanged.
	 */
	public void shift(int delta);
	
	
	/**
	 * Move an annotation, preserving the relationships between its blocks, to a new coordinate.
	 * @param coordinateInReference
	 */
	public void moveToCoordinate(int coordinateInReference);
	
	
	/**
	 * @return A deep clone of this Annotation object
	 */
	public Annotation copy();
		
	
	/**
	 * Fragments the current annotation if it overlaps the provided one.
	 * it returns an empty list if this annotation does not overlap the
	 * one passed in.
	 * 
	 * @param Annotation annotation to disect the current one
	 * @return List<Annotation> of the disected annotations, an empty list is returned if 
	 * 			the annotations do not overlap.
	 */
	List<Annotation> disect(Annotation a);
	
	
	/**
	 * Fragments the current annotation if it overlaps the provided ones.
	 * It returns a list with one component (this annotation) if no annotation 
	 * in the provided list overlaps the discted annotaion.
	 * 
	 * @param List<GenomicAnnotation> <b>sorted</b> annotations with which to disect the current one.
	 * @return List<GenomicAnnotation> of the disected annotations, a list just this annotation is returned if 
	 * 			the annotations do not overlap.
	 */
	List<Annotation> disect(List<? extends Annotation> disectors);
	
	
	/**
	 * Returns the difference (all regions in this genomic annotation not in the given other annotation)
	 * between this genomic annotation and the one
	 * 
	 * @param other - the annotations to take the difference against
	 * @return A blocked annotation of the resulting difference.
	 */
	public Annotation minus(Annotation other);
	
	/**
	 * Returns the difference (all regions in this genomic annotation not in the given list)
	 * between this genomic annotation and the given list
	 * 
	 * @param others - the annotations to take the difference against
	 * @return A blocked annotation of the resulting difference.
	 */
	public Annotation minus(Collection<? extends Annotation> others); 
	
	/**
	 * Calculates the distance to the another genomic annotation.
	 * @return 0 if the annotations overlap or the minimum between the edges otherwise.
	 */
	public int getDistanceTo(Annotation other);

	
	/**
	 * @return chrX:start-end
	 */
	public String toUCSC();
	public String toBED();
	public String toBED(int r, int g, int b);
	public String toShortBED();
	public String toBEDGraph();
	
	/**
	 * @return A string with no whitespace that includes all information about the annotation
	 */
	public String getFullInfoString();
	
	/**
	 * 
	 * @param other - other genomic annotation
	 * @param buffer if the overlap is within buffer they will be considered overlapping 
	 *        even if they do not overlap within their original boundaries.
	 * @return true if they overlap in this extended definition, not considering orientation
	 */
	public boolean overlaps(Annotation other, int buffer);

	/**
	 * 
	 * @param other - other genomic annotation
	 * @param buffer if the overlap is within buffer they will be considered overlapping 
	 *        even if they do not overlap within their original boundaries.
	 * @return true if they overlap in this extended definition, not considering orientation
	 */
	public boolean overlaps(Collection<? extends Annotation> others, int buffer);

	/**
	 * 
	 * @param other GenomicAnnotation to check of overlap
	 * @return true if the current instance overlaps with the other one, not considering orientation
	 */
	public boolean overlaps(Annotation other);

	/**
	 * 
	 * @param other GenomicAnnotations to check of overlap
	 * @return true if the current instance overlaps with the other one, not considering orientation
	 */
	public boolean overlaps(Collection<? extends Annotation> others);
	
	public boolean overlaps(Annotation other, int buffer, boolean considerOrientation);
	public boolean overlaps(Annotation other, boolean considerOrientation);
	public boolean overlapsStranded(Annotation other);
	
	/**
	 * Returns the number of bases that overlap between the two annotations.
	 * @param other
	 * @return
	 */
	public int getOverlap(Annotation other) ;

	/**
	 * @param other
	 * @return true if the blocks of other are fully contained within the blocks of this object
	 */
	public boolean contains(Annotation other);
	
	/**
	 * @param other
	 * @return union of current instance and other object.  If strands agree, then strand is set to that.  If the strands conflict, strand is set to UNKNOWN
	 */
	public Annotation union(Annotation other);

	/**
	 * Returns the result of intersecting this instance with another
	 * @param other
	 * @return returns null if there is no intersection
	 */
	public Annotation intersect(Annotation other);


	/**
	 * Returns the result of intersecting this instance with a list of annotations,
	 * with one (possibly blocked) intersection Annotation for each provided annotation
	 * @param other annotations
	 * @return List of intersections or empty list if none
	 */
	public List<Annotation> intersect(List<? extends Annotation> annotations);

	
	/**
	 * Should behave like compareTo for a basic annotation
	 * @param b
	 * @return
	 */
	public int compareToAnnotation(Annotation b);


	public void stitchTo(Annotation next);  // TODO what does this do?

	public String toCigar();
	

	/**
	 * Test whether the annotation is fully contained in this
	 * @param annotation
	 * @return
	 */
	public boolean fullyContains(Annotation annotation);
	
	/**
	 * @return a new CompoundInterval containing the gaps between the blocks in this annotation
	 */
	public Annotation complement();

	/**
	 * Test whether two annotations are equal regardless of whether the orientation are the same
	 * @param other The annotation to test equality to
	 * @param useOrientation true=ensure strands match false=OK if strands dont match
	 * @return
	 */
	public boolean equals(Annotation other, boolean useOrientation);

	/**
	 * This will return a collection of splice junctions
	 * This is not the same as all gaps, it will only return splice junctions
	 * @return A collection of splice junctions
	 */
	public Collection<? extends Annotation> getSpliceConnections();

	/**
	 * This will return the position in the center of this annotation taking into account the blocked structure
	 * @return
	 */
	public int getMidpoint();

}
