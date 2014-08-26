/**
 * 
 */
package broad.pda.seq.clip;

import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

import nextgen.core.annotation.Gene;
import nextgen.core.coordinatesystem.GenomicSpace;
import nextgen.core.model.AlignmentModel;
import nextgen.core.model.score.MultiScore;
import nextgen.core.model.score.ScanStatisticScore;
import nextgen.core.readFilters.MappingQualityFilter;
import nextgen.core.readFilters.NumHitsFilter;

/**
 * @author prussell
 *
 */
public class GenomeSpaceSampleData extends SampleData {

	private AlignmentModel genomeData;
	private Map<Gene, ScanStatisticScore> genomeScores;
	private static String SCAN_STAT_KEY = "scan_statistic";

	/**
	 * @param bamFile Bam alignment file
	 * @param chrSizeFile Chromosome size file
	 * @param genes Genes by chromosome
	 * @param window Window size
	 * @param step Step size
	 * @param expressionCutoff Genome wide scan P value cutoff for expression
	 * @param fullyContained Count fully contained reads only
	 * @throws IOException
	 */
	public GenomeSpaceSampleData(String bamFile, String chrSizeFile, Map<String, Collection<Gene>> genes, double expressionCutoff, boolean fullyContained) throws IOException {
		super(bamFile, false, genes, expressionCutoff, true, fullyContained);
		genomeData = new AlignmentModel(bamFile, new GenomicSpace(chrSizeFile));
		genomeData.addFilter(new MappingQualityFilter(5,10));
		genomeData.addFilter(new NumHitsFilter(1));
		genomeScores = new TreeMap<Gene, ScanStatisticScore>();
	}

	/**
	 * Whether the gene is significantly expressed genome wide
	 * @param gene The gene
	 * @return Whether the gene is expressed at the given significance level
	 */
	@Override
	public boolean isExpressed(Gene gene) {
		if(!genomeScores.containsKey(gene)) {
			ScanStatisticScore score = new ScanStatisticScore(genomeData, gene, fullyContainedReads);
			logger.debug("CHECK_GENE_EXPRESSION\t" + gene.getName());
			logger.debug("CHECK_GENE_EXPRESSION\t" + gene.getChr() + ":" + gene.getStart() + "-" + gene.getEnd());
			logger.debug("CHECK_GENE_EXPRESSION\t" + "global_length=" + score.getGlobalLength());
			logger.debug("CHECK_GENE_EXPRESSION\t" + "global_count=" + score.getTotal());
			logger.debug("CHECK_GENE_EXPRESSION\t" + "global_lambda=" + score.getGlobalLambda());
			logger.debug("CHECK_GENE_EXPRESSION\t" + "window_size=" + score.getCoordinateSpace().getSize(gene));
			logger.debug("CHECK_GENE_EXPRESSION\t" + "window_count=" + score.getCount());
			logger.debug("CHECK_GENE_EXPRESSION\t" + "pval=" + score.getPvalue());
			genomeScores.put(gene, score);			
		}
		ScanStatisticScore score = genomeScores.get(gene);
		return score.getPvalue() <= expressionCutoffValue;
	}


	/**
	 * Get genome wide scan P value of number of fragments mapping to the gene
	 * Get score from cache or calculate and cache score
	 * @param gene The gene
	 * @return The scan P value of the number of fragments mapping to the gene
	 */
	public double getGenePval(Gene gene) {
		if(!genomeScores.containsKey(gene)) {
			ScanStatisticScore score = new ScanStatisticScore(genomeData, gene, fullyContainedReads);
			logger.debug("GET_GENE_SCAN_PVAL\t" + gene.getName());
			logger.debug("GET_GENE_SCAN_PVAL\t" + gene.getChr() + ":" + gene.getStart() + "-" + gene.getEnd());
			logger.debug("GET_GENE_SCAN_PVAL\t" + "global_length=" + score.getGlobalLength());
			logger.debug("GET_GENE_SCAN_PVAL\t" + "global_count=" + score.getTotal());
			logger.debug("GET_GENE_SCAN_PVAL\t" + "global_lambda=" + score.getGlobalLambda());
			logger.debug("GET_GENE_SCAN_PVAL\t" + "window_size=" + score.getCoordinateSpace().getSize(gene));
			logger.debug("GET_GENE_SCAN_PVAL\t" + "window_count=" + score.getCount());
			logger.debug("GET_GENE_SCAN_PVAL\t" + "pval=" + score.getPvalue());
			genomeScores.put(gene, score);			
		}
		ScanStatisticScore score = genomeScores.get(gene);
		return score.getPvalue();
	}

	
}
