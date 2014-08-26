/**
 * 
 */
package broad.pda.countreads;

import broad.core.parser.CommandLineParser;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;

import org.apache.commons.collections15.Predicate;
import org.apache.log4j.Logger;

import broad.core.math.Statistics;
import broad.pda.annotation.BEDFileParser;

import nextgen.core.alignment.Alignment;
import nextgen.core.alignment.AbstractPairedEndAlignment.TranscriptionRead;
import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.BasicAnnotation;
import nextgen.core.annotation.Gene;
import nextgen.core.coordinatesystem.GenomicSpace;
import nextgen.core.coordinatesystem.TranscriptomeSpace;
import nextgen.core.model.AlignmentModel;
import nextgen.core.model.score.CountScore;

/**
 * @author prussell
 *
 */
public class AnnotationCounts {

	private Map<String, Collection<Gene>> genes;
	private AlignmentModel transcriptomeData;
	private AlignmentModel genomeData;
	private static Logger logger = Logger.getLogger(AnnotationCounts.class.getName());
	private boolean includePositionLevelInfo;
	private int minCountForAvgWithMin;
	
	
	private AnnotationCounts(String bamFile, String bedFile, String chrSizeFile, boolean positionLevelInfo, int minForAvgWithMin) throws IOException {
		includePositionLevelInfo = positionLevelInfo;
		minCountForAvgWithMin = minForAvgWithMin;
		logger.info("Loading genes...");
		genes = BEDFileParser.loadDataByChr(new File(bedFile));
		logger.info("Creating transcriptome space...");
		transcriptomeData = new AlignmentModel(bamFile, new TranscriptomeSpace(genes), new ArrayList<Predicate<Alignment>>(), true, TranscriptionRead.UNSTRANDED, true, null);
		logger.info("Creating genomic space...");
		genomeData = new AlignmentModel(bamFile, new GenomicSpace(chrSizeFile));
		logger.info("Done loading data.");
	}
	
	@SuppressWarnings("unused")
	private void writeCounts(String outFile) throws IOException {
		writeCounts(outFile, null);
	}
	
	private static double getAverageCountPerPositionWithThreshold(double[] counts, int minCount) {
		List<Double> nonzeroCounts = new ArrayList<Double>();
		for(int i=0; i<counts.length; i++) {
			double count = counts[i];
			if(count > minCount) {
				nonzeroCounts.add(Double.valueOf(count));
			}
		}
		return Statistics.mean(nonzeroCounts);
	}
	
	private static double getAverageCountPerPosition(double[] counts) {
		return Statistics.mean(counts);
	}
	
	private static double getMaxPositionCount(double[] counts) {
		double max = 0;
		for(int i=0; i<counts.length; i++) {
			double count = counts[i];
			if(count > max) {
				max = count;
			}
		}
		return max;
	}

	
	private void writeCounts(String outFile, String singleChr) throws IOException {
		logger.info("Writing counts to file " + outFile);
		FileWriter w = new FileWriter(outFile);
		String header = "gene_name\t";
		header += "count_full_span\t";
		header += "count_annotation\t";
		header += "count_introns_fully_contained\t";
		header += "fpkm\t";
		if(includePositionLevelInfo) {
			header += "average_position_coverage\t";
			header += "average_coverage_positions_above_" + minCountForAvgWithMin + "_only\t";
			header += "max_position_coverage\t";
		}
		w.write(header + "\n");
		Collection<String> chrs = new TreeSet<String>();
		if(singleChr != null) {
			chrs.add(singleChr);
		} else {
			chrs.addAll(genes.keySet());
		}
		for(String chr : chrs) {
			if(!genes.containsKey(chr)) {
				continue;
			}
			for(Gene gene : genes.get(chr)) {
				logger.info(gene.getName() + "\t" + gene.toUCSC());
				int start = gene.getStart();
				int end = gene.getEnd();
				double fullSpanCount = genomeData.getCount(new BasicAnnotation(chr, start, end));
				double annotationCount = transcriptomeData.getCount(gene);
				double intronCount = 0;
				for(Annotation intron : gene.getIntronSet()) {
					intronCount += genomeData.getCount(intron, true);
				}
				CountScore countScore = new CountScore(transcriptomeData, gene);
				double rpkm = countScore.getRPKM();
				String write = gene.getName() + "\t";
				write += fullSpanCount + "\t";
				write += annotationCount + "\t";
				write += intronCount + "\t";
				write += rpkm + "\t";
				if(includePositionLevelInfo) {
					double[] positionCounts = transcriptomeData.getCountsPerPosition(gene);
					double averageCoverage = getAverageCountPerPosition(positionCounts);
					write += averageCoverage + "\t";
					double averageCoverageWithThreshold = getAverageCountPerPositionWithThreshold(positionCounts, minCountForAvgWithMin);
					write += averageCoverageWithThreshold + "\t";
					double maxCoverage = getMaxPositionCount(positionCounts);
					write += maxCoverage + "\t";
					
				}
				w.write(write + "\n");
			}
		}
		w.close();
		logger.info("Done writing counts.");
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {

		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-b", "Bam file", true);
		p.addStringArg("-g", "Bed gene annotation", true);
		p.addStringArg("-s", "Chromsome size file", true);
		p.addStringArg("-c", "Single chromosome to write", false, null);
		p.addStringArg("-o", "Output file", true);
		p.addBooleanArg("-p", "Include position level information for each annotation (slow)", false, false);
		p.addIntArg("-m", "Min count to include a position for -p option", false, 20);
		p.parse(args);
		String bamFile = p.getStringArg("-b");
		String bedFile = p.getStringArg("-g");
		String sizeFile = p.getStringArg("-s");
		String outFile = p.getStringArg("-o");
		String singleChr = p.getStringArg("-c");
		boolean positionInfo = p.getBooleanArg("-p");
		int minCount = p.getIntArg("-m");
		
		AnnotationCounts a = new AnnotationCounts(bamFile, bedFile, sizeFile, positionInfo, minCount);
		a.writeCounts(outFile, singleChr);
		
	}

}
