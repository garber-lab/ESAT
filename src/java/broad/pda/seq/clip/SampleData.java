/**
 * CLIP/Protect seq dataset for one sample
 */
package broad.pda.seq.clip;

import broad.core.parser.StringParser;

import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

import org.apache.log4j.Logger;


import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Gene;
import nextgen.core.coordinatesystem.TranscriptomeSpace;
import nextgen.core.model.TranscriptomeSpaceAlignmentModel;
import nextgen.core.model.score.ScanStatisticScore;
import nextgen.core.readFilters.FragmentLengthFilter;
import nextgen.core.readFilters.GenomicSpanFilter;
import nextgen.core.readFilters.MappingQualityFilter;
import nextgen.core.readFilters.NumHitsFilter;


/**
 * @author shari
 *
 */
public class SampleData {

	protected String sampleName;
	protected TranscriptomeSpaceAlignmentModel data;
	protected TranscriptomeSpaceAlignmentModel maxFragmentLengthData;
	protected Map<Gene, ScanStatisticScore> maxFragmentLengthGeneScores;
	protected Map<Gene, ScanStatisticScore> geneScores;
	protected Map<Gene, Double> geneAvgCoverage;
	public static Logger logger = Logger.getLogger(SampleData.class.getName());
	protected Map<String, Collection<Gene>> genesByChr;
	protected Map<String, Gene> genesByName;
	protected double expressionCutoffValue;
	private static int DEFAULT_MAX_GENOMIC_SPAN = 100000;
	private static int DEFAULT_MAX_FRAGMENT_LENGTH = 150;
	protected boolean expressionByScanPval;
	private String originalBamFile;
	private boolean read1TranscriptionStrand;
	protected boolean fullyContainedReads;
	
	/**
	 * @param bamFile Bam file
	 * @param firstReadTranscriptionStrand Whether read 1 is in direction of transcription
	 * @param genes Genes by chromosome
	 * @param window Window size
	 * @param step Step size
	 * @param expressionCutoff P value cutoff for scan test of gene expression
	 * @param expByScanPval Expression is assessed by scan P value. If false, uses average read depth
	 * @param fullyContained Count fully contained reads only
	 * @throws IOException 
	 */
	public SampleData(String bamFile, boolean firstReadTranscriptionStrand, Map<String, Collection<Gene>> genes, double expressionCutoff, boolean expByScanPval, boolean fullyContained) throws IOException {
		fullyContainedReads = fullyContained;
		geneAvgCoverage = new TreeMap<Gene, Double>();
		originalBamFile = bamFile;
		read1TranscriptionStrand = firstReadTranscriptionStrand;
		StringParser p = new StringParser();
		p.parse(bamFile, "\\.");
		String withoutSuffix = p.asString(0);
		p.parse(withoutSuffix, "/");
		sampleName = p.asString(p.getFieldCount() - 1);
		for(int i = 1 ; i < p.getFieldCount() - 1; i++) {
			sampleName += "." + p.asString(i);
		}
		expressionCutoffValue = expressionCutoff;
		expressionByScanPval = expByScanPval;
		genesByChr = genes;
		data = new TranscriptomeSpaceAlignmentModel(bamFile, new TranscriptomeSpace(genes));
		maxFragmentLengthData = new TranscriptomeSpaceAlignmentModel(bamFile, new TranscriptomeSpace(genes));
		
		// Read filters
		data.addFilter(new GenomicSpanFilter(DEFAULT_MAX_GENOMIC_SPAN));
		data.addFilter(new MappingQualityFilter(5,10));
		data.addFilter(new NumHitsFilter(1));
		
		maxFragmentLengthData.addFilter(new GenomicSpanFilter(DEFAULT_MAX_GENOMIC_SPAN));
		maxFragmentLengthData.addFilter(new MappingQualityFilter(5,10));
		maxFragmentLengthData.addFilter(new NumHitsFilter(1));
		maxFragmentLengthData.addFilter(new FragmentLengthFilter(maxFragmentLengthData.getCoordinateSpace(), DEFAULT_MAX_FRAGMENT_LENGTH));
		
		genesByName = new TreeMap<String, Gene>();
		for(String chr : genesByChr.keySet()) {
			for(Gene gene : genesByChr.get(chr)) {
				genesByName.put(gene.getName(), gene);
			}
		}
		geneScores = new TreeMap<Gene, ScanStatisticScore>();
		maxFragmentLengthGeneScores = new TreeMap<Gene, ScanStatisticScore>();
		logger.info("Instantiated sample data object. Name = " + sampleName);
	}
	
	@Override
	public int hashCode() {
		return sampleName.hashCode();
	}
	
	/**
	 * Set genome wide scan P value cutoff for expression of transcript
	 * @param expressionScanPvalCutoff P value cutoff for transcript expression against genomic background
	 */
	public void setExpressionScanPvalueCutoff(double expressionScanPvalCutoff) {
		expressionCutoffValue = expressionScanPvalCutoff;
	}

	/**
	 * Get genome wide scan P value cutoff for expression of transcript
	 * @return P value cutoff for transcript expression against genomic background
	 */
	public double getExpressionScanPvalueCutoff() {
		return expressionCutoffValue;
	}
	
	/**
	 * Get number of fragments mapping to gene
	 * Get score from cache or calculate and cache score
	 * @param gene The gene
	 * @return The number of fragments mapping to the gene
	 */
	public double getGeneCount(Gene gene) {
		if(geneScores.containsKey(gene)) {
			return geneScores.get(gene).getCount();
		}
		ScanStatisticScore score = new ScanStatisticScore(data, gene, fullyContainedReads);
		return score.getCount();
	}
	
	/**
	 * Get number of fragments mapping to gene after fragment length filter
	 * Get score from cache or calculate and cache score
	 * @param gene The gene
	 * @return The number of fragments mapping to the gene
	 */
	public double getFragmentLengthGeneCount(Gene gene) {
		if(maxFragmentLengthGeneScores.containsKey(gene)) {
			return maxFragmentLengthGeneScores.get(gene).getCount();
		}
		ScanStatisticScore score = new ScanStatisticScore(maxFragmentLengthData, gene, fullyContainedReads);
		return score.getCount();
	}
	
	/**
	 * Get coordinate space wide scan P value of number of fragments mapping to the gene
	 * Get score from cache or calculate and cache score
	 * @param gene The gene
	 * @return The scan P value of the number of fragments mapping to the gene with respect to teh coordinate space
	 */
	public double getGeneScanPval(Gene gene) {
		if(geneScores.containsKey(gene)) {
			return geneScores.get(gene).getScanPvalue();
		}
		ScanStatisticScore score = new ScanStatisticScore(data, gene, fullyContainedReads);
		geneScores.put(gene, score);
		return score.getScanPvalue();
	}
	
	
	/**
	 * Get average coverage of gene
	 * Get score from cache or calculate and cache score
	 * @param gene The gene
	 * @return The average coverage of the gene
	 */
	public double getGeneAverageCoverage(Gene gene) {
		if(geneAvgCoverage.containsKey(gene)) {
			return geneAvgCoverage.get(gene).doubleValue();
		}
		ScanStatisticScore score = new ScanStatisticScore(data, gene, fullyContainedReads);
		geneScores.put(gene, score);
		double avgCoverage = score.getAverageCoverage(data);
		geneAvgCoverage.put(gene, Double.valueOf(avgCoverage));
		return avgCoverage;
	}
	
	/**
	 * Whether the gene is significantly expressed in the sample
	 * @param gene The gene
	 * @return Whether the gene is expressed at the given significance level
	 */
	public boolean isExpressed(Gene gene) {
		if(!geneScores.containsKey(gene)) {
			ScanStatisticScore score = new ScanStatisticScore(data, gene, fullyContainedReads);
			logger.debug("CHECK_GENE_EXPRESSION\t" + gene.getName());
			geneScores.put(gene, score);			
		}
		ScanStatisticScore score = geneScores.get(gene);
		if(expressionByScanPval) {
			return score.getScanPvalue() <= expressionCutoffValue;
		}
		double avgDepth = getGeneAverageCoverage(gene);
		logger.debug("cutoff=" + expressionCutoffValue + "\tcount=" + score.getCount() + "\tsize=" + score.getCoordinateSpace().getSize(score.getAnnotation()) + "\tscore=" + avgDepth);
		return avgDepth >= expressionCutoffValue;
	}
	
	
	/**
	 * Whether read 1 is in direction of transcription
	 * @return True iff first read transcription strand is true
	 */
	public boolean firstReadTranscriptionStrand() {
		return read1TranscriptionStrand;
	}
	
	/**
	 * Get enrichment of a window over a gene
	 * @param gene The gene
	 * @param window Window contained in the gene
	 * @return Enrichment of window over gene background
	 */
	public double getEnrichmentOverGene(Gene gene, Annotation window) {
		if(!gene.contains(window)) {
			throw new IllegalArgumentException("Gene must contain window.");
		}
		double geneAvgCov = getGeneAverageCoverage(gene);
		ScanStatisticScore windowScore = new ScanStatisticScore(data, window, data.getCount(gene), gene.getSize());
		double windowAvgCoverage = windowScore.getAverageCoverage(data);
		double enrichment = windowAvgCoverage / geneAvgCov;
		return enrichment;
	}

	/**
	 * Get the logger
	 * @return The logger
	 */
	@SuppressWarnings("static-method")
	public Logger getLogger() {
		return logger;
	}
	
	/**
	 * Get the sample name
	 * @return Sample name
	 */
	public String getSampleName() {
		return sampleName;
	}
	
	/**
	 * Get name of original bam file
	 * @return Bam file name
	 */
	public String getOriginalBamFile() {
		return originalBamFile;
	}
	
	/**
	 * Set the sample name
	 * @param name Sample name
	 */
	public void setSampleName(String name) {
		sampleName = name;
	}
	
/*	*//**
	 * Write the window scores to a file for future use
	 * @param m Multi sample binding site caller that this belongs to
	 * @throws IOException
	 *//*
	public void writeWindowScoresToFile(MultiSampleBindingSiteCaller m) throws IOException {
		windowScoreFile.writeWindowScoresToFile(m);
	}
*/	
	/**
	 * Get the alignment data
	 * @return Alignment model
	 */
	public TranscriptomeSpaceAlignmentModel getData() {
		return data;
	}
	
	/**
	 * Get the alignment data filtered by max fragment length
	 * @return Alignment model
	 */
	public TranscriptomeSpaceAlignmentModel getFragmentLengthFilterData() {
		return maxFragmentLengthData;
	}

	
}
