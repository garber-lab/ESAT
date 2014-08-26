/**
 * 
 */
package broad.core.math;

import broad.core.parser.CommandLineParser;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import broad.pda.annotation.BEDFileParser;

import net.sf.samtools.util.CloseableIterator;
import nextgen.core.alignment.Alignment;
import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.AnnotationList;
import nextgen.core.annotation.Gene;
import nextgen.core.coordinatesystem.CoordinateSpace;
import nextgen.core.coordinatesystem.TranscriptomeSpace;
import nextgen.core.feature.GeneWindow;
import nextgen.core.model.AlignmentModel;
import nextgen.core.model.score.CountScore;
import nextgen.core.model.score.WindowScoreIterator;
import nextgen.core.readFilters.GenomicSpanFilter;
import nextgen.core.writers.WigWriter;

/**
 * @author prussell
 *
 */
public class PermutationDistribution {
	
	private Annotation region;
	private int regionSize;
	private AlignmentModel data;
	private int numPerms;
	private Map<Integer, EmpiricalDistribution> scanDistributions;
	private CoordinateSpace coord;
	private double totalCount;
	private ArrayList<AnnotationList<Alignment>> permutations;
	private static int DEFAULT_MAX_GENOMIC_SPAN = 300000;
	
	private static Logger logger = Logger.getLogger(PermutationDistribution.class.getName());
	
	/**
	 * @param parentRegion Parent region
	 * @param bamFile Bam file
	 * @param bedFile Bed file for transcriptome space
	 * @param numPermutations Number of permutations
	 * @throws IOException
	 */
	public PermutationDistribution(Annotation parentRegion, String bamFile, String bedFile, int numPermutations) throws IOException {
		this(parentRegion, bamFile, bedFile, numPermutations, DEFAULT_MAX_GENOMIC_SPAN);
	}
	
	/**
	 * @param parentRegion Parent region
	 * @param bamFile Bam file
	 * @param bedFile Bed file for transcriptome space
	 * @param numPermutations Number of permutations
	 * @param maxGenomicSpan Max genomic span for paired reads
	 * @throws IOException
	 */
	public PermutationDistribution(Annotation parentRegion, String bamFile, String bedFile, int numPermutations, int maxGenomicSpan) throws IOException {
		this(parentRegion, new AlignmentModel(bamFile, TranscriptomeSpace.getTranscriptomeSpace(bedFile)), TranscriptomeSpace.getTranscriptomeSpace(bedFile), numPermutations, maxGenomicSpan);
	}
	
	/**
	 * @param parentRegion Parent annotation
	 * @param alignments Collection of annotations to count e.g. read mappings
	 * @param coordinateSpace Coordinate space
	 * @param numPermutations Number of permutations
	 */
	public PermutationDistribution(Annotation parentRegion, AlignmentModel alignments, CoordinateSpace coordinateSpace, int numPermutations) {
		this(parentRegion, alignments, coordinateSpace, numPermutations, DEFAULT_MAX_GENOMIC_SPAN);
	}

	/**
	 * @param parentRegion Parent annotation
	 * @param alignments Collection of annotations to count e.g. read mappings
	 * @param coordinateSpace Coordinate space
	 * @param numPermutations Number of permutations
	 * @param maxGenomicSpan Max genomic span for paired reads
	 */
	public PermutationDistribution(Annotation parentRegion, AlignmentModel alignments, CoordinateSpace coordinateSpace, int numPermutations, int maxGenomicSpan) {
		logger.debug("Instantiating permutation distribution for parent region " + parentRegion.getName());
		region = parentRegion;
		regionSize = region.getSize();
		data = alignments;
		data.addFilter(new GenomicSpanFilter(maxGenomicSpan));
		totalCount = data.getCount(region);
		logger.debug("Region size is " + regionSize + ". Total count is " + totalCount + ".");
		numPerms = numPermutations;
		coord = coordinateSpace;
		scanDistributions = new TreeMap<Integer, EmpiricalDistribution>();
		generatePermutations();
	}
	
	/**
	 * Get scan P value of window
	 * @param window The window
	 * @return The percentage of permutations with a window of the same size that has more overlapping fragments
	 */
	public double getScanPvalue(Annotation window) {
		if(!region.contains(window)) {
			throw new IllegalArgumentException("Parent region " + region.getName() + " does not contain annotation " + window.getName() + ".");
		}
		int windowSize = window.getSize();
		double windowCount = data.getCount(window);
		double rtrn = 1 - getScanDistribution(windowSize).getCumulativeProbability(windowCount);
		logger.debug("Window " + window.getName() + " size " + windowSize + " count " + windowCount + " scan P value " + rtrn);
		return rtrn;
	}
	
	/**
	 * Get the max count of all count scores returned by the iterator
	 * @param iter Iterator
	 * @return The max count over all scores
	 */
	private static double maxCount(WindowScoreIterator<CountScore> iter) {
		double rtrn = Double.MIN_VALUE;
		while(iter.hasNext()) {
			CountScore score = iter.next();
			double count = score.getCount();
			score = null;
			if(count > rtrn) {
				rtrn = count;
			}
		}
		return rtrn;
	}
	
	/**
	 * Get scan distribution for the window size
	 * @param windowSize Window size
	 * @return Scan distribution
	 */
	public EmpiricalDistribution getScanDistribution(int windowSize) {
		if(!scanDistributions.containsKey(Integer.valueOf(windowSize))) {
			calculateDistribution(windowSize);
		}
		return scanDistributions.get(Integer.valueOf(windowSize));
	}
	
	/**
	 * Calculate and store scan distribution for the window size
	 * @param windowSize Window size
	 */
	private void calculateDistribution(int windowSize) {
		logger.debug("Calculating scan distribution for windows of size " + windowSize + "...");
		EmpiricalDistribution dist = new EmpiricalDistribution(regionSize, 0, totalCount);
		int printStep = Math.max(1, permutations.size() / 100);
		int numDone = 0;
		for(AnnotationList<Alignment> permutation : permutations) {
			WindowScoreIterator<CountScore> iter = permutation.scan(region, windowSize, windowSize-1);
			double maxCount = maxCount(iter);
			dist.add(maxCount);
			numDone++;
			if(numDone % printStep == 0) {
				logger.debug("Finished " + numDone + "/" + permutations.size() + " permutations.");
			}
		}
		scanDistributions.put(Integer.valueOf(windowSize), dist);
	}
	
	/**
	 * Generate all permutations
	 */
	private void generatePermutations() {
		logger.debug("");
		logger.debug("Generating permutations...");
		int printStep = Math.max(1, numPerms / 100);
		int numDone = 0;
		permutations = new ArrayList<AnnotationList<Alignment>>();
		for(int i = 0; i < numPerms; i++) {
			// Get iterator over randomized annotations
			CloseableIterator<Alignment> permutedIter = data.getPermutedAnnotations(region);
			// Instantiate annotation list with randomized annotations
			ArrayList<Alignment> randomizedAnnotations = new ArrayList<Alignment>();
			while(permutedIter.hasNext()) {
				randomizedAnnotations.add(permutedIter.next());
			}
			permutedIter.close();
			AnnotationList<Alignment> randomizedAnnotationList = new AnnotationList<Alignment>(coord, randomizedAnnotations);			
			permutations.add(randomizedAnnotationList);
			numDone++;
			if(numDone % printStep == 0) {
				logger.debug("Made " + numDone + "/" + numPerms + " permutations.");
			}
		}
		logger.debug("Done generating permutations.");
	}
	
	/**
	 * Write a wig file of -log10(Pvalue) at each position of annotation
	 * @param outWig Output wig file
	 * @throws IOException
	 */
	public void writeWigEachPosition(String outWig) throws IOException {
		logger.info("");
		logger.info("Writing -log10(Pvalue) at each position to file " + outWig + "...");
		Gene gene = new Gene(region);
		Collection<GeneWindow> windows = gene.getWindows(1);
		int numWindows = windows.size();
		int printStep = Math.max(1, numWindows / 100);
		int numDone = 0;
		FileWriter w = new FileWriter(outWig);
		TreeMap<Integer, Double> counts = new TreeMap<Integer, Double>();
		for(GeneWindow window : windows) {
			int start = window.getStart();
			double pval = getScanPvalue(window);
			double score = -1 * Math.log10(pval);
			counts.put(Integer.valueOf(start), Double.valueOf(score));
			numDone++;
			if(numDone % printStep == 0) {
				logger.debug("Finished " + numDone + "/" + numWindows + " windows.");
			}
		}
		WigWriter.write(w, region.getChr(), counts, false);
		w.close();
		logger.info("Done writing file.");		
	}
	
	/**
	 * Write windows to bed file with scan P value in score field
	 * @param windowSize Window size
	 * @param outBed Output bed file
	 * @throws IOException
	 */
	public void writeWindowsWithPvalue(int windowSize, String outBed) throws IOException {
		logger.info("");
		logger.info("Writing windows of size " + windowSize + " with P-values to file " + outBed + "...");
		Gene gene = new Gene(region);
		Collection<GeneWindow> windows = gene.getWindows(windowSize);
		int numWindows = windows.size();
		int printStep = Math.max(1, numWindows / 100);
		int numDone = 0;
		FileWriter w = new FileWriter(outBed);
		for(GeneWindow window : windows) {
			window.setName(region.getName() + "_window_" + numDone);
			double pval = getScanPvalue(window);
			window.setBedScore(pval);
			w.write(window + "\n");
			numDone++;
			if(numDone % printStep == 0) {
				logger.debug("Finished " + numDone + "/" + numWindows + " windows.");
			}
		}
		w.close();
		logger.info("Done writing file.");
	}
	
	/**
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-b", "Bed file for transcriptome space", true);
		p.addStringArg("-a", "Bam file of alignments", true);
		p.addStringArg("-g", "Name of parent gene in bed file", true);
		p.addIntArg("-p", "Number of permutations", true);
		p.addIntArg("-w", "Window size to write to bed file", false, 1);
		p.addStringArg("-ob", "Output bed file of windows", false, null);
		p.addStringArg("-ow", "Output wig file of -log10(Pvalue) at each position", false, null);
		p.addBooleanArg("-d", "Debug logging", false, false);
		p.addIntArg("-mg", "Max genomic span for paired reads", false, DEFAULT_MAX_GENOMIC_SPAN);
		p.parse(args);
		String bedFile = p.getStringArg("-b");
		String bamFile = p.getStringArg("-a");
		String geneName = p.getStringArg("-g");
		int numPerms = p.getIntArg("-p");
		Gene parentGene = BEDFileParser.getGeneByName(geneName, bedFile);
		int windowSize = p.getIntArg("-w");
		String outBed = p.getStringArg("-ob");
		String outWig = p.getStringArg("-ow");
		boolean debug = p.getBooleanArg("-d");
		int maxGenomicSpan = p.getIntArg("-mg");
		
		if(debug) {
			logger.setLevel(Level.DEBUG);
		}
		
		PermutationDistribution pd = new PermutationDistribution(parentGene, bamFile, bedFile, numPerms, maxGenomicSpan);
		
		if(outBed != null) {
			pd.writeWindowsWithPvalue(windowSize, outBed);
		}
		
		if(outWig != null) {
			pd.writeWigEachPosition(outWig);
		}
		
		logger.info("");
		logger.info("All done.");
		
	}
	
	
}
