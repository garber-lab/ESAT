package broad.pda.seq.segmentation;

import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;

import org.apache.commons.collections15.Predicate;

import net.sf.samtools.util.CloseableIterator;
import nextgen.core.alignment.Alignment;
import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Gene;

import broad.core.annotation.LightweightGenomicAnnotation;
import broad.core.datastructures.IntervalTree;
import broad.core.sequence.Sequence;

//TODO Add updateInterval tree to interface

public interface AlignmentDataModel{

	//Get all Annotation overlapping a specified region
	public CloseableIterator<Alignment> getAnnotationOverlappingRegion(Annotation align) throws IOException;
	
	//Get all represented chromosomes and their lengths in the genome
	public Map<String, Integer> getChromosomeLengths();

	/**
	 * Get length of a chromosome
	 * @param chr
	 * @return
	 */
	public int getChromosomeLength(String chr);

	/**
	 * Counts the number of reads for the specified chromosome
	 * @return
	 * @throws IOException 
	 */
	public double getCounts(String chr) throws IOException;
	
	//get the number of Annotation overlapping a given region
	public double getCountsPerAlignment(Annotation first, int EF) throws IOException;

	//get the number of Annotation overlapping a given region with a cached interval tree
	public double getCountsPerAlignment(Annotation align, IntervalTree<Alignment> tree, int EF);
	
	//get the number of Annotation overlapping a given region with a cached interval tree
	public double getCountsPerAlignment(Gene gene, IntervalTree<Alignment> tree, int EF);

	//get the number of Annotation overlapping a given region with a cached interval tree
	public double getCountsPerAnnotationtranded(Gene gene, IntervalTree<Alignment> tree, int EF, String orientation);

	//given an array of Annotation generate the counts over the whole region
	public double getCountsPerAlignment(Annotation[] Annotation, IntervalTree<Alignment> tree, int EF);
	
	//get the number of reads that are fully contained within the region (doesn't start before or extend past)
	public double getCountsPerAlignmentFullyContained(Gene gene, IntervalTree<Alignment> tree, int EF);
	
	/**
	 * Similar to @see getCountsPerAlignment(Annotation[] Annotation, IntervalTree<Alignment> tree, int EF)
	 * but delegates the management of the alignment interval tree to the implementing class rather than the caller
	 * @param Annotation
	 * @param EF
	 * @return
	 * @throws IOException 
	 */
	public double getCountsPerAlignment(Annotation[] Annotation,  int EF) throws IOException;

	public double getCountsPerAlignment(Annotation align,
	Map<String, IntervalTree<Annotation>> goodExonTree,
	IntervalTree<Alignment> tree, int extensionFactor);

	public int getCountsWithinExons(Annotation align, Iterator<Alignment> iter, int EF);

	public int getCountsWithinExons(Annotation align, IntervalTree<Alignment> tree, int EF);
	
	public double getCountsOfUniqueOverlappers(Annotation target, Annotation exclude, IntervalTree<Alignment> tree, int EF);
	
	public double getCountsOfUniqueOverlappers(Annotation target, Gene exclude, IntervalTree<Alignment> tree, int EF);
	
	public double getCountsOnBase(String chr, int index) throws IOException;
	
	public int getBasesCovered(Annotation record, int EF) throws IOException;
	
	//get the number of spliced reads overlapping a given region
	public double getNumberOfSplicedReads(Annotation align) throws IOException;
	
	public GeneCounts getGeneCounts(Gene gene, int extensionFactor) throws IOException;

	//get the actual number of bases covered by Annotation overlapping a given region with a cached interval tree (if a single 76 base read spans the gene, return 76)
	public double getBasesCoveredPerAlignment(Gene gene, IntervalTree<Alignment> tree, int EF);
			
	//get the exonic regions overlapping a given region
	public Iterator<Annotation> getExonAnnotationOverlappingRegion(Annotation align) throws IOException;
	
	public Gene getPeak(Gene gene, Gene startCodon, IntervalTree<Alignment> tree, int EF);
	
	public Collection<Annotation> getOverlappingRegions(String chr) throws IOException;
	
	//compute a cached interval tree for a given region
	//Replace with counter per node
	public IntervalTree<Alignment> getIntervalTree(String chr, int start, int end) throws IOException;
	//public IntervalTree<Annotation> getIntervalTreeCachedAnnotation(String chr, int start, int end);
	
	
	public IntervalTree<Alignment> getIntervalTreeCached(String chr, int start, int end) throws IOException;

	public Gene updateGeneByFirstCounts(Gene gene, IntervalTree<Alignment> tree, int EF);
	
	//get all spliced reads overlapping a region
	//break into intronic regions and store them as introns
	public Map<Annotation, Integer> getSplicedReads(Annotation align, int minIntronSize, int maxIntronSize) throws IOException;
	
	public Map<Annotation, Integer> getSplicedReads(Annotation region, Collection<Predicate<Annotation>> filters) throws IOException;
	public Map<Annotation, Integer> getSplicedReads(Annotation region) throws IOException;
	public Map<Annotation, Integer> getSplicedReads(Annotation align, final int minIntronSize, final int maxIntronSize , int minNumIntrons, Sequence chrSeq) throws IOException;
	
	public Map<Annotation, Integer> getSplicedReads(Annotation region, Collection<Predicate<Annotation>> filters, int minNumIntrons) throws IOException;

	public Map<Annotation, Integer> getSplicedReads(Annotation region, Collection<Predicate<Annotation>> filters, int minNumIntrons, Sequence chrSeq) throws IOException;

	public Map<Annotation, Integer> getSplicedReads(Annotation region, Predicate<Annotation> filter) throws IOException;

	public Map<Annotation, Integer> getSplicedReads(int minIntronSize, int maxIntronSize) throws IOException;

	public double getSpliceWeightFactor();

	public Collection<? extends Annotation> getSplicedReadExons(String chr) throws IOException;

	public CloseableIterator<Alignment> getReadIterator();
	public CloseableIterator<Alignment> getReadIterator(Annotation region) throws IOException;
	
	
	
	//public int getFirstReads(String chr);
	
	/**
	 * @param the chromosome over which to iterate
	 * @return CloseableIterator over all reads mapped to the given Chromosome
	 */
	public CloseableIterator<Alignment> getChromosomeReadIterator(String chromosome) throws IOException;

	public double getTotalNumberOfStrandedReads();

	public double getScorePerAlignmentFromCache(Gene window, IntervalTree<Alignment> chunkAlignmentTree, int extensionFactor);
	public void resetTreeCache();

	public int getChunkStart();

	public int getChunkEnd();
	
	public IntervalTree<Alignment> getIntervalTreeTruncatedCached(String chromosome, int start, int end) throws IOException;

	public IntervalTree<Annotation> getFullIntervalTreeAsAnnotation(String chr) throws IOException;

	
	public double getMinimumMappingQuality();

	public void setChunkSize(int chunkSize);

	public void setNegativeStranded() ;

	public void setPositiveStranded() ;
	
	public void setSecondRead() ;
	
	public void setFirstRead() ;

	void setNormalizationStrategy(ReadCountNormalizer normalizer);

	void setMinimumMappingQuality(double minimumMapQuality);

	//get the number of Annotation overlapping a given region with a cached interval tree
	public boolean passes(Gene gene, IntervalTree<Alignment> tree, int EF, double cutoff);

	//get the number of Annotation overlapping a given region with a cached interval tree
	public long passesCutoff(Annotation align, IntervalTree<Alignment> tree, int EF, double threshold);

	public boolean hasNextExon(String chr) throws IOException;
	public void restartIterator();

	public void resetTreeCache(String chr);
	
	public void clearFullIntervalTreeAsAlignmentsCached(String chr);
	
	public boolean isStranded();
	
	public void unsetStranded();
	public boolean isPositiveStranded();
	public boolean isNegativeStranded() ;
	public String getModelFilePath();
	void setExtensionFactor(int factor);
}
