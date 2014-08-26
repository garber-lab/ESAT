package broad.pda.seq.rap;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

import broad.core.annotation.MaximumContiguousSubsequence;
import broad.core.annotation.ShortBEDReader;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.util.Log;
import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.BasicAnnotation;
import nextgen.core.model.AlignmentModel;
import nextgen.core.model.score.RatioScore;
import nextgen.core.model.score.WindowProcessor;
import nextgen.core.model.score.WindowScoreIterator;

import broad.core.math.EmpiricalDistribution;
import org.apache.commons.math3.stat.descriptive.rank.*;

public class RatioPermutationPeakCaller extends GenomeCommandLineProgram {
    private static final Log log = Log.getInstance(RatioPermutationPeakCaller.class);
    
    @Usage
    public String USAGE = "Call peaks based on a permutation model of ratios.  Run broad.pda.rap.BuildRatioNullDistribution first.";
    
	@Option(doc="Window size")
	public int WINDOW;
	
	@Option(doc="Overlap between windows")
	public int OVERLAP;

	@Option(doc="Input SAM or BAM file")
	public File TARGET;

	@Option(doc="Control SAM or BAM file for normalization.")
	public File CONTROL;
	
	@Option(doc="Output file basename")
	public File OUTPUT;
	
	@Option(doc="File containing empirical null distribution of ratios", optional=true)
	public File NULL_DIST = null;
    
	@Option(doc="Percentile cut-off (e.g. 0.001)", optional=true)
	public Double CUTOFF = 0.001;
	
	@Option(doc="Whether to look for regions BELOW the given cutoff, rather than above.")
	public Boolean DEPLETED=false;
	
	@Option(doc="Percentile cut-off for peak trimming", optional=true)
	public Double TRIM = null;

	
	
	/**
	 * Stock main method.
	 *
	 * @param args main arguments
	 */
	public static void main(final String[] args) {
		System.exit(new RatioPermutationPeakCaller().instanceMain(args));
	}
	

	@Override
	protected int doWork() {
		
		try {
			if (OUTPUT.exists()) {
				OUTPUT.delete();
			}
			
			List<Annotation> regions = getRegions();

			AlignmentModel target = loadAlignmentModel(TARGET);
			AlignmentModel control = loadAlignmentModel(CONTROL);

			double ratioCutoff = getRatioCutoff();
			log.info("Cutoff = " + ratioCutoff);
			
			
			BufferedWriter bw = new BufferedWriter(new FileWriter(OUTPUT,true));
			BufferedWriter trimWriter = null;
			if (TRIM != null) trimWriter = new BufferedWriter(new FileWriter(new File(OUTPUT.getAbsolutePath()+".trim.bed"),true));
			
			for (Annotation region : regions) {
				log.info("Starting: " + region.toUCSC());
				RatioScore.Processor processor = new RatioScore.Processor(target, control);
				
				List<Annotation> sigWindows = getSignificantWindows(target, control, ratioCutoff, region, processor);
				List<RatioScore> scoredWindows = scoreWindows(processor, sigWindows);
				writeResults(scoredWindows, bw);
				
				if (TRIM != null) {
					trimPeaks(target, control, sigWindows, WINDOW/10, TRIM, processor.getNumeratorRegionTotal(), processor.getDenominatorRegionTotal());
					scoredWindows = scoreWindows(processor, sigWindows);
					writeResults(scoredWindows, bw);
				}
			}
			
			bw.close();
			if (TRIM != null) trimWriter.close();
			
		} catch (Exception e) {
			log.error(e);
		}
		
		return 0;
	}

	
	private double getRatioCutoff() throws IOException {
		EmpiricalDistribution ed = new EmpiricalDistribution(NULL_DIST);
		if (!DEPLETED) CUTOFF = 1 - CUTOFF;
		return ed.getQuantileOfBins(CUTOFF);
	}

	
	private List<Annotation> getSignificantWindows(AlignmentModel target, AlignmentModel control, double ratioCutoff, Annotation region, WindowProcessor<RatioScore> processor) throws IOException {
		// TODO: write or use existing SlideAndCount file
		LinkedList<Annotation> sigWindows = new LinkedList<Annotation>();
		EmpiricalDistribution allWindows = BuildRatioNullDistribution.getEmptyEmpiricalDistribution();
		
		int count = 0;
		WindowScoreIterator<RatioScore> itr = target.scan(region, WINDOW, OVERLAP, processor);
		while (itr.hasNext()) {
			RatioScore curr = itr.next();
			allWindows.add(curr.getLog2RegionRatio());
			
			if ((DEPLETED && curr.getLog2RegionRatio() <= ratioCutoff) ||
			    (!DEPLETED && curr.getLog2RegionRatio() >= ratioCutoff)) {
				
				count++;
				
				// If this window overlaps the previous one, combine into one window
				boolean overlapping = false;
				if (sigWindows.size() > 0) {
					Annotation previous = sigWindows.getLast();
					if (curr.getAnnotation().overlaps(previous)) {
						sigWindows.removeLast();
						sigWindows.addLast(previous.union(curr.getAnnotation()));
						overlapping = true;
					}
				} 

				if (!overlapping) sigWindows.addLast(curr.getAnnotation());
			}
		}
		itr.close();
		
		// Write empirical distribution of all window scores
		BufferedWriter bw = new BufferedWriter(new FileWriter(new File(OUTPUT.getAbsolutePath() + ".allWindowDistribution.txt")));
		allWindows.write(bw);
		bw.close();
		
		log.info("Found " + count + " significant windows.");
		log.info("Collapsed to " + sigWindows.size() + " non-overlapping regions.");
		return sigWindows;
	}
	
	
	private List<RatioScore> scoreWindows(WindowProcessor<RatioScore> processor, List<? extends Annotation> windows) {
		// Using the same processor initialized above preserves the region counts
		
		LinkedList<RatioScore> scored = new LinkedList<RatioScore>();
		for (Annotation window : windows) {
			scored.addLast(processor.processWindow(window));
		}
		return scored;
	}
	
	
	private void writeResults(List<RatioScore> scored, BufferedWriter bw) throws IOException {
		for (RatioScore score : scored) {
			bw.write(score.toString());
			bw.newLine();
		}
	}
	
	
	private void trimPeaks(AlignmentModel target, AlignmentModel control, List<? extends Annotation> peaks, int step, double quantileCutoff, double numeratorRegionTotal, double denominatorRegionTotal) throws IOException {
		for (Annotation peak : peaks) {
			trimPeak(target, control, peak, step, quantileCutoff, numeratorRegionTotal, denominatorRegionTotal);
		}
	}
	
	private void trimPeak(AlignmentModel target, AlignmentModel control, Annotation peak, int step, double quantileCutoff, double numeratorRegionTotal, double denominatorRegionTotal) throws IOException {
		// total hack code JE 1/8/13
		
		log.info("Trimming " + peak.toUCSC() + " " + peak.numBlocks()); 
		
		// Calculate ratios for all subwindows 
		WindowProcessor<RatioScore> processor = new RatioScore.Processor(target, control);
		WindowScoreIterator<RatioScore> itr = target.scan(peak, step, 0, processor);
		List<RatioScore> scores = new LinkedList<RatioScore>();
		while (itr.hasNext()) scores.add(itr.next());
		double[] ratios = new double[scores.size()];
		for (int i = 0; i < scores.size(); i++) {
			RatioScore score = scores.get(i);
			score.setNumeratorRegionTotal(numeratorRegionTotal);
			score.setDenominatorRegionTotal(denominatorRegionTotal);
			ratios[i] = score.getLog2RegionRatio();
			if (DEPLETED) ratios[i] = -1 * ratios[i];
		}
		
		// First use the maximum subsequence algorithm.  All windows that pass the global cutoff get a positive score,
		// and those that don't get a negative score.
		double zero = getRatioCutoff();
		if (DEPLETED) zero = -1 * zero;
		for (int i = 0; i < ratios.length; i++) {
			ratios[i] -= zero;
		}
		
		double[] result = MaximumContiguousSubsequence.maxSubSum3(ratios); // requires that some values are positive

		// Next try the quantile trim method
		// RESULTS:  it looks like this part of the method is mostly unnecessary
		/*
		if (DEPLETED) quantileCutoff = 1 - quantileCutoff;
		double cutoff = new Percentile().evaluate(ratios, quantileCutoff*100);
		System.out.println("Cutoff = " + cutoff);
		if (!DEPLETED) cutoff = -1 * cutoff - zero;
		int midIndex = scores.size() / 2;
		int endIndex = MaximumContiguousSubsequence.contiguousStartSubSequenceOverMin(Arrays.copyOfRange(ratios, midIndex, ratios.length), cutoff) + midIndex - 1;
		System.out.println("endIndex = " + endIndex);
		int startIndex = MaximumContiguousSubsequence.contiguousEndSubSequenceOverMin(Arrays.copyOfRange(ratios, 0, midIndex), cutoff) + 1;
		System.out.println("startIndex = " + startIndex); 
		
		// Take the largest possible segment that passes either cutoff
		startIndex = Math.min(startIndex, (int) result[1]);
		endIndex = Math.max(endIndex, (int) result[2]);
		*/
		int startIndex = (int) result[1];
		int endIndex = (int) result[2];
		// Move the peak to the new coordinates
		peak.moveToCoordinate(scores.get(startIndex).getAnnotation().getStart());
		peak.setEnd(scores.get(endIndex).getAnnotation().getEnd());		
		log.info("Trimmed to: " + peak.toUCSC());
	}
}
