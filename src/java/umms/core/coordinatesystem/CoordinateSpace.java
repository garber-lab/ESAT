package umms.core.coordinatesystem;

import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.io.File;
import java.io.IOException;

import umms.core.annotation.Annotation;
import umms.core.annotation.BasicAnnotation;
import umms.core.annotation.Gene;
import umms.core.feature.GeneWindow;
import umms.core.feature.Window;
import broad.core.annotation.ShortBEDReader;
import broad.core.datastructures.IntervalTree;
import broad.core.error.ParseException;

/**
 * This interface represents the coordinate space in which we analyze the models
 * @author skadri
 *	@author mguttman
 */
public interface CoordinateSpace {

	/**
	 * Returns a fragment of type Window between start and end on chromosome
	 * @param chr
	 * @param start
	 * @param end
	 * @return
	 */
	Collection<? extends Window> getFragment(String chr, int start, int end);
	
	/**
	 * Returns a fragment of type Window between start and end on chromosome
	 * @param Annotation
	 * @return
	 */
	Collection<? extends Window> getFragment(Annotation annotation);
	
	/**
	 * Get the size of the annotation with respect to the coordinate space
	 * @param region
	 * @return
	 */
	public int getSize(Annotation region);
	
	/**
	 * Returns a collection of regions overlapping start and end on chromosome
	 * @param chr
	 * @param start
	 * @param end
	 * @return
	 */
	Collection<? extends Window> getOverlappingRegion(String chr,int start,int end);
	
	/**
	 * Iterate over windows within the defined region
	 * @param windowSize
	 * @param chr
	 * @param start
	 * @param end
	 * @return returns an iterator that iterates over the region
	 */
	Iterator<? extends Window> getWindowIterator(int windowSize,String chr,int start,int end,int overlap);
	
	/**
	 * 
	 * @return returns an iterator that goes from first window to last window on chromosome chr sequentially
	 */
	Iterator<? extends Window> getWindowIterator(String chr, int windowSize, int overlap);
		
	/**
	 * 
	 * @return returns an iterator that goes from first window to last window in the coordinate space sequentially
	 */
	Iterator<? extends Window> getWindowIterator(int windowSize, int overlap);
	
	/**
	 * Get a window iterator over a collection of GeneWindows
	 * @param baseGenes The regions over which to iterate
	 * @param windowSize Window size
	 * @param overlap Overlap size
	 * @return Iterator over windows on the collection of regions
	 */
	public Iterator<Window> getWindowIterator(Collection<Gene> baseGenes, int windowSize, int overlap);

	
	/**
	 * @return returns an iterator over an annotation sequentially
	 */
	Iterator<? extends Window> getWindowIterator(Annotation window, int windowSize, int overlap);
	
	/**
	 * Returns the total number of positions in the coordinate space
	 * For non-contiguous coordinate spaces, this will return a collapsed count such that overlapping blocks are not double counted
	 * @return number of positions in space
	 */
	long getLength();
	
	/**
	 * Returns the total number of positions in the coordinate space
	 * For non-contiguous coordinate spaces, this will return a collapsed count such that overlapping blocks are not double counted
	 * @param only positions on a given chromosome
	 * @return number of positions in space
	 */
	long getLength(String chr);
	
	
	/**
	 * Permutes an annotation to a new location on the same reference in the given coordinate space, preserving the 
	 * length of the annotation and maintaining a consistent block structure.  Does not check whether the given annotation 
	 * was valid in the coordinate space in the first place.
	 * @param a Annotation to permute.  Will be modified.
	 * @return The permuted Annotation
	 * @throws PermutationNotFoundException if an annotation cannot be found
	 */
	public Annotation permuteAnnotation(Annotation a);
	
	/**
	 * Permutes an annotation to a new location within the given bounds, preserving the 
	 * length of the annotation and maintaining a consistent block structure.  Does not 
	 * check whether the given annotation was valid in the coordinate space in the first place.
	 * @param a Annotation to permute.  Will be modified.
	 * @param bounds Boundaries within which the annotation will be permuted
	 * @return The permuted Annotation
	 * @throws PermutationNotFoundException if a permuted location cannot be found
	 */
	public Annotation permuteAnnotation(Annotation a, Annotation bounds);
	
	
	/**
	 * @return A list of reference names contained in the coordinate space (e.g. chr1, chr2, etc.)
	 */
	public Collection<String> getReferenceNames();
	
	/**
	 * @return A list of chromosome names contained in the coordinate space
	 */
	public Collection<String> getChromosomeNames();
	
	/**
	 * @return The reference annotations contained in this coordinate space.
	 */
	public Collection<? extends Annotation> getReferenceAnnotations();  // TODO should eventually be an AnnotationReader
	
	/** 
	 * @param name
	 * @return The named reference as an annotation
	 */
	public Annotation getReferenceAnnotation(String name);
	
	/**
	 * Returns true if the coordinate space contains any annotations on the queried chromosome
	 * @param chr
	 * @return
	 */
	public boolean hasChromosome(String chr);
	
	/**
	 * Get an annotation spanning an entire chromosome
	 * @param chrName The chromosome name
	 * @return Annotation consisting of the entire chromosome
	 */
	public Annotation getEntireChromosome(String chrName);
	
	/**
	 * Check if window in coordinate space and not masked
	 * @param window
	 * @return
	 */
	public boolean isValidWindow(Annotation window);
}

