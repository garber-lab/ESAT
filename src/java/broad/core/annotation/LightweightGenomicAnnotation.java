package broad.core.annotation;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Annotation.Strand;

public interface LightweightGenomicAnnotation extends Annotation{

	public int getStart();

	public int getEnd();
	
	long getMiddle();
	
	public double getScore();
	
	public String getName();

	public void setScore(double score);

	public void setName(String name);

	public Strand getOrientation();
	
	public void setOrientation(String orientation);

	public int length();

	public void setStart(int start);

	public void setEnd(int end);

	public String getChromosome();

	public void setChromosome(String chr);
	
	boolean inReversedOrientation();
	
	/**
	 * Calculates the distance to the another genomic annotation.
	 * @return 0 if the annotations overlap or the minimum between the edges otherwise.
	 */
	public int getDistanceTo(LightweightGenomicAnnotation other) ;
	
	/**
	 * 
	 * Similar to getLocationString() @see getLocationString and provided only  for backwards 
	 * compatibility.
	 */
	public String toUCSC();
	
	
	/**
	 * 
	 * @param other GenomicAnnotation to check of overlap
	 * @return true if the current instance overlaps with the other one.
	 */
	//public boolean overlaps(LightweightGenomicAnnotation other);
	
	
	/**
	 * Returns the number of bases that overlap between the two annotations.
	 * @param other
	 * @return
	 */
	public int getOverlap(LightweightGenomicAnnotation other) ;
	
	public boolean contains(LightweightGenomicAnnotation other) ;
	
	/**
	 * Change the current instance to represent its intersection
	 * with the provided annotation 
	 * @param other
	 */
	public void takeIntersection(LightweightGenomicAnnotation other);
	
	/**
	 * Change the current instance to represent its union
	 * with the provided annotation 
	 * @param GenomicAnnotation other
	 */
	public void takeUnion(LightweightGenomicAnnotation other);
	
	/**
	 * Enlarges this annotation by 
	 * @param other
	 */
	public void stitchTo(LightweightGenomicAnnotation other);
	
	/**
	 * Returns the usual way to write a chromosome: chrSymbol if a chromosome or the full name of the scaffold 
	 * 
	 */
	public String getChromosomeString();
	

	/**
	 * Returns the ith extra score
	 */
	public double getExtraScore(int i);
	
	/**
	 * Adds a new extra score
	 * @param score
	 */
	public void addExtraScore (double score);
	
	/**
	 * 
	 * @return all extra scores
	 */
	public List<Double> getExtraScores();
	
	/**
	 * removes extra scores
	 */
	public void removeExtraScores();
	
	/**
	 * Checks whether two annotations differ by a small (fudge) factor
	 * @param other The other annotation
	 * @param fudge Maximum difference at either end or start to consider similar
	 * @return
	 */
	public boolean almostEqual(LightweightGenomicAnnotation other, int fudge) ;

}