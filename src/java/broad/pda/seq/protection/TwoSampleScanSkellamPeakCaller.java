/**
 * 
 */
package broad.pda.seq.protection;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;

import org.apache.log4j.Logger;

import broad.core.annotation.MaximumContiguousSubsequence;
import broad.core.math.Distribution;
import broad.core.math.Statistics;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Gene;
import nextgen.core.annotation.Annotation.Strand;
import nextgen.core.coordinatesystem.CoordinateSpace;
import nextgen.core.coordinatesystem.TranscriptomeSpace;
import nextgen.core.feature.GeneWindow;
import nextgen.core.model.ScanStatisticDataAlignmentModel;
import nextgen.core.model.TranscriptomeSpaceAlignmentModel;
import nextgen.core.model.score.CountScore;
import nextgen.core.model.score.ScanStatisticScore;
import nextgen.core.model.score.WindowProcessor;

/**
 * @author prussell
 *
 */
public class TwoSampleScanSkellamPeakCaller {
	
	protected Map<String, Map<Gene, Double>> backgroundScanPvalues;
	protected ScanStatisticDataAlignmentModel backgroundData;
	protected String backgroundName;
	protected TranscriptomeSpaceAlignmentModel signalData;
	protected String signalName;
	protected Map<String, Collection<Gene>> genes;
	
	/**
	 * Transcriptome wide scan P value cutoff for expression of gene
	 */
	public static double EXPRESSION_PVALUE_CUTOFF = 0.01;
	protected static int MIN_WINDOW_COUNT = 10;
	protected int windowSize;
	protected int stepSize;
	protected double pValueCutoffSkellam;
	protected double pValueCutoffScan;
	protected double trimPeakByQuantile;
	protected static Logger logger = Logger.getLogger(TwoSampleScanSkellamPeakCaller.class.getName());
	protected TranscriptomeSpace transcriptomeSpace;
	protected WindowProcessor<CountScore> backgroundProcessor;
	protected WindowProcessor<ScanStatisticScore> signalProcessor;

	public Collection<Annotation> filterWindows(Collection<Annotation> windows) throws IOException {
		TreeSet<Annotation> rtrn = new TreeSet<Annotation>();
		for(Annotation window : windows) {
			double scanPval = signalData.scoreWindow(new GeneWindow(window)).getScanPvalue();
			if(scanPval < pValueCutoffScan) rtrn.add(window);
		}
		return rtrn;
	}

	
	public void writeResults(Collection<Annotation> windows, String out) throws IOException {
		FileWriter w = new FileWriter(out);
		for(Annotation window : windows)	{
			w.write(window.toBED() + "\n");
		}
		w.close();
	}
	
	public void writeResult(Collection<Annotation> windows, FileWriter writer) throws IOException {
		for(Annotation window : windows)	{
			writer.write(window.toBED() + "\n");
		}
	}


	/**
	 * Given a total number of reads, a region length, and a window size, get the parameter of the Poisson distribution modeling the number of reads in a window of this size
	 * @param regionLength The total region length
	 * @param d The number of reads
	 * @param windowSize The window size
	 * @return Lambda, the expected number of reads in the window under the null Poisson model
	 */
	public static double getPoissonLambda(int regionLength, double d, int windowSize) {
		if(d <= 0) {
			throw new IllegalArgumentException("Can't get poisson expected number of reads for " + d + " in window.");
		}
		return windowSize * d / regionLength;
	}


	/**
	 * Compute Skellam P-value of read counts in a region given the two parameters
	 * @param backgroundLambda Poisson lambda for background sample
	 * @param signalLambda Poisson lambda for signal sample
	 * @param backgroundCount The count for background sample
	 * @param signalCount The count for signal sample
	 * @return The probability under the null hypothesis of observing a greater difference
	 * @throws IOException
	 */
	protected static double getSkellamPvalue(double backgroundLambda, double signalLambda, int backgroundCount, int signalCount) throws IOException {
		return Distribution.skellamRightTail(signalCount - backgroundCount, signalLambda, backgroundLambda);
	}

	/**
	 * Trim the peak to max contiguous subsequence of counts above a certain quantile
	 * @param window The peak region
	 * @return The trimmed peak
	 * @throws IOException
	 */
	protected Annotation trimPeak(Gene window) throws IOException {

		List<Double> coverageData = signalData.getPositionCountList(window);
		Annotation rtrn = SampleData.trimMaxContiguous(window, coverageData, trimPeakByQuantile);
		int trimmedStart = rtrn.getStart();
		int trimmedEnd = rtrn.getEnd();
		rtrn.setName(rtrn.getReferenceName() + ":" + trimmedStart + "-" + trimmedEnd);
		rtrn.setOrientation(window.getOrientation());	
		return rtrn;

	}


	
}
