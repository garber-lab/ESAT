package umms.core.annotation;

import java.io.IOException;
import java.util.Collection;
import java.util.List;

import net.sf.samtools.util.CloseableIterator;
import umms.core.coordinatesystem.CoordinateSpace;
import umms.core.model.score.CountScore;
import umms.core.model.score.WindowProcessor;
import umms.core.model.score.WindowScore;
import umms.core.model.score.WindowScoreIterator;

import org.apache.commons.collections15.Predicate;

/**
 * @author engreitz
 * Interface to allow for counting over both AnnotationReader and ReadCollection. 
 */
public interface AnnotationCollection<T extends Annotation> {
	
	/**
	 * @return Total number of annotations.  Considers filters.
	 */
	public int size();
	
	
	/**
	 * @return Global count, possibly weighted.  Might be different from size if only counting
	 * reads that map, etc.
	 */
	public double getGlobalCount();
	
	/**
	 * This method will return the number of reads overlapping a given region
	 * @param region The reference region to count overlapping elements
	 * @param fullyContained whether to only count overlapping elements if they are fully contained within the reference region
	 * @return returns a double to allow implementing methods to utilize a weight object
	 * @throws IOException
	 */
	public double getCount(Annotation region, boolean fullyContained);
	
	/**
	 * This method will return the number of reads overlapping a given region (not necessarily fully contained)
	 * @param region The reference region to count overlapping elements
	 * @return returns a double to allow implementing methods to utilize a weight object
	 * @throws IOException
	 */
	public double getCount(Annotation region);
	
	
	public double getCount(AnnotationList<? extends Annotation> set);
	public double getCount(AnnotationList<? extends Annotation> set, boolean fullyContained);
	
	/**
	 * This method will return the total number of bases that are covered by alignments in the collection
	 * If discontinuous region, it will only count blocks
	 * The maximum possible value is defined by region.getTranscriptSize();
	 * @param region The region to determine overlap
	 * @param fullyContained Whether to count only fully contained reads
	 * @return
	 */
	public int getBasesCovered(Annotation region, boolean fullyContained);
	
	/**
	 * This method will return the total number of bases that are covered by alignments in the collection
	 * If discontinuous region, it will only count blocks.
	 * The maximum possible value is defined by region.getTranscriptSize();
	 * Annotations counted here do not need to be fully contained in the region.
	 * @param region The region to determine overlap
	 * @return
	 */
	public int getBasesCovered(Annotation region);

	
	/**
	 * Returns the total base-level coverage for the given region.
	 * @param region
	 * @return
	 */
	public double getTotalCoverage(Annotation region);
	
	/**
	 * Returns the average base-level coverage for the given region.
	 * @param region
	 * @return
	 */
	public double getAverageCoverage(Annotation region);
	
	/**
	 * Collapse all reads within the defined region into the overlapping set
	 * @param region The region to determine overlap
	 * @param fullyContained Whether to count only fully contained reads
	 * @return null if there are no Annotations that overlap
	 */
	public Annotation collapse(Annotation region, boolean fullyContained);
	

	
	/**
	 * This method will return a collection of the annotations that overlap a given region (not necessarily fully contained)
	 * @param region The region to overlap
	 * @return
	 */
	public CloseableIterator<T> getOverlappingAnnotations(Annotation region);
	
	
	/**
	 * This method will return a collection of the annotations that overlap a given region
	 * @param region The region to determine overlap
	 * @param fullyContained Whether to count only fully contained reads
	 * @return
	 */
	public CloseableIterator<T> getOverlappingAnnotations(Annotation region, boolean fullyContained);
	

	/**
	 * @return A closeable iterator through all annotations that pass the filters. 
	 */
	public CloseableIterator<T> iterator();
	
	
	/**
	 * Adds a filter to use to decide whether an alignment is valid
	 * @param filter the filter to use
	 */
	public void addFilter(Predicate<T> filter);
	
	
	/**
	 * Add a set of filters to decide whether an alignment is valid
	 * @param filters Collection of filters to add
	 */
	public void addFilters(Collection<Predicate<T>> filters);
	
	
	/**
	 * Shuffle the locations of annotations within the given region and return an iterator
	 * over the permuted annotations.  Permuted annotations will not be necessarily sorted.
	 * NOTE: This will alter the annotations in the reader, so re-reading or reloading 
	 * might be necessary.
	 * @param region reads reads fully contained in this region will be shuffled
	 */
	public CloseableIterator<T> getPermutedAnnotations(Annotation region);

	
	/**
	 * Compute the reads fully contained in region that are NOT overlapping excluded
	 * @param region
	 * @param excluded
	 * @return Count contained in region that do not overlap the excluded region
	 */
	public double getCountExcludingRegion(Annotation region, Annotation excluded);
	
	public double getCountStrandedExcludingRegion(Annotation region, Annotation excluded);
	
	/**
	 * @return The coordinate space currently used by the AnnotationCollection
	 */
	public CoordinateSpace getCoordinateSpace();
	
	
	public WindowScoreIterator<CountScore> scan(int windowSize, int overlap);
	public <W extends WindowScore> WindowScoreIterator<W> scan(int windowSize, int overlap, WindowProcessor<W> processor);
	public <W extends WindowScore> WindowScoreIterator<W> scan(Annotation region, int windowSize, int overlap, WindowProcessor<W> processor);
	public WindowScoreIterator<CountScore> scan(Annotation region, int windowSize, int overlap);


	public double getRefSequenceLambda(String refname);
	
}
