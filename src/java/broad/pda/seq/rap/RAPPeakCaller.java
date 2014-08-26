package broad.pda.seq.rap;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;

import broad.core.math.Distribution;
import broad.core.math.EmpiricalDistribution;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.annotation.BEDFileParser;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.AnnotationList;
import nextgen.core.annotation.BasicAnnotation;
import nextgen.core.annotation.Gene;
import nextgen.core.coordinatesystem.CoordinateSpace;
import nextgen.core.coordinatesystem.GenomicSpace;
import nextgen.core.coordinatesystem.MaskedGenomicSpace;
import nextgen.core.exception.IncompleteMethodImplementationException;
import nextgen.core.feature.Window;
import nextgen.core.model.ScanStatisticDataAlignmentModel;

/**
 * 
 * @author skadri
 *
 */

public class RAPPeakCaller {


	private ScanStatisticDataAlignmentModel backgroundSample;
	private ScanStatisticDataAlignmentModel rapSample;
	
	private CoordinateSpace maskedGenome; 
	
	private int windowSize;
	private int step;
	private double skellamCutoff;
	private double klCutoff;
	private double scanCutoff;
	
	Map<String, Collection<Gene>> maskedRegions;
	
	protected static Logger logger = Logger.getLogger(RAPPeakCaller.class.getName());
	
	Map<String,Integer> chrSizes;
	
	
	static final String usage = "Usage: RAPPeakCaller -task <task name> "+
			"\n\tscore: Computes significant peaks in RAP data from the input sample " + 
			"\n\t\t-background <Alignment (mapped to genome) of input/AS sample> "+
			"\n\t\t-rap <Alignment (mapped to genome) of RAP sample> "+
			"\n\t\t-sizes <Chromosome sizes>"+
			"\n\t\t-maskedRegions <masked regions>"+
			"\n\t\t-window <windowSize>"+
			"\n\t\t-step <Step size>"+
			"\n\t\t-skellam <max skellam p-value>"+
			"\n\t\t-scan <max scan p-value>"+
			"\n\t\t-out <Output file>";

	/**
	 * 
	 * @param backgroundFile Alignment file for the input sample
	 * @param rapFile Alignment file for the rap file
	 * @param sizes Chromosome sizes
	 * @param windowSize Window Size
	 * @param step Step Size / Overlap between consecutive windows
	 * @param maxSkellam
	 * @param maxScanPvalue
	 * @throws IOException
	 */
	public RAPPeakCaller(String backgroundFile, String rapFile, String sizes, String maskedregions, int windowSize, int step) throws IOException {
		
		chrSizes = BEDFileParser.loadChrSizes(sizes);
		logger.info("Chromosome sizes are loaded");
		//maskedGenome = new MaskedGenomicSpace(chrSizes, BEDFileParser.loadAlignmentDataByChr(new File(maskedRegions)), 0);
		maskedGenome = new GenomicSpace(sizes, maskedregions,50);
		logger.info("Masked Genome space is created");
		backgroundSample = new ScanStatisticDataAlignmentModel(backgroundFile, maskedGenome);
		logger.info("Background sample:"+backgroundFile+" is read");
		rapSample = new ScanStatisticDataAlignmentModel(rapFile, maskedGenome);
		logger.info("RAP sample:"+rapFile+" is read");
		
		this.windowSize = windowSize;
		this.step = step;
		
		maskedRegions = BEDFileParser.loadDataByChr(new File(maskedregions));
	}
	
	/**
	 * Call peaks using skellam
	 * @param outFile
	 * @throws IOException
	 */
	public void callSkellamPeaks(String outFile,double maxSkellam, double maxScanPvalue) throws IOException{
		
		this.skellamCutoff = maxSkellam;
		this.scanCutoff = maxScanPvalue;
		BufferedWriter bw = new BufferedWriter(new FileWriter(outFile));
		BufferedWriter bwV = new BufferedWriter(new FileWriter(outFile+".values"));
		
		TreeSet<Window> sigWindows = new TreeSet<Window>();
		
		for(String chr:chrSizes.keySet()){
			
			double backgroundLambda = backgroundSample.getRefSequenceLambda(chr);
			double length = (double)chrSizes.get(chr) - (double)getMaskedLength(chr);
			backgroundLambda = (backgroundLambda*backgroundSample.getRefSequenceLength(chr))/length;
			double rapLambda = rapSample.getRefSequenceLambda(chr);
			rapLambda = (rapLambda*rapSample.getRefSequenceLength(chr))/length;
			
			logger.info("Background lambda for "+chr+": "+backgroundLambda);
			logger.info("RAP lambda for "+chr+": "+rapLambda);
			
			if(rapLambda ==0.0){
				logger.error(chr +" is not expressed in RAP file");
			}
			else{
				//logger.info("Get the window iterator:");
				//Get an iterator over the masked chr
				//Step is overlap
				Iterator<? extends Window> iter = maskedGenome.getWindowIterator(chr,this.windowSize,this.step);
				//logger.info("Got the window iterator");
				while(iter.hasNext()){
					//logger.info("Enter the loop.");
					Window window = iter.next();
					//logger.info(window.toUCSC());
					int backgroundCount = (int) Math.round(backgroundSample.getCount(window));
					
			/*		if(backgroundCount <= 0) 
						continue;*/
					
					int rapCount = (int) Math.round(rapSample.getCount(window));
					//TODO: If signal is low, do not consider
					if(rapCount < rapLambda) {
						//logger.info("Count for this window: "+rapCount+" is less than lambda: "+rapLambda);
						continue; 
					}
					
					double skellamPval = getSkellamPValue(backgroundLambda, rapLambda, backgroundCount, rapCount);
					//Calculate scan statistic
					double scanPval = rapSample.scoreWindow(window).getScanPvalue();
					
					bwV.write(window.toUCSC()+"\t"+skellamPval+"\t"+scanPval+"\n");
					
					if(skellamPval <= this.skellamCutoff) {
						logger.info("Skellam P-value "+skellamPval+" for "+window.toUCSC()+" is less than cutoff: "+skellamCutoff);
						// filter for scan P value
						if(scanPval < scanCutoff) {
							logger.info("Scan statistic P-value "+scanPval+" for "+window.toUCSC()+" is less than cutoff: "+scanCutoff);
							window.setScore(skellamPval);
							sigWindows.add(window);
							bw.write(window.toBED() + "\n");
						}
						else{
							logger.info("Scan statistic P-value "+scanPval+" for "+window.toUCSC()+" is more than cutoff: "+scanCutoff);
						}
					}
					else{
						if(skellamPval<1.0)
							logger.info("Skellam P-value "+skellamPval+" for "+window.toUCSC()+" is more than cutoff: "+skellamCutoff);
					}
				}
			}
		}
		
		bw.close();
		bwV.close();
	}
	
	
	private int getMaskedLength(String chr){
		
		int len = 0;
		if(maskedRegions.containsKey(chr)){
			for(Gene g:maskedRegions.get(chr)){
				len = g.getSize();
			}
		}
		return len;
	}
	/**
	 * For each window, calculate the poisson lambda for the given window over the entire chromosome
	 * @param model
	 * @param window
	 * @param windowSize
	 * @return
	 * @throws IOException 
	 */
/*	private double getPoissonLambda(DataAlignmentModel model, Window window,int windowSize) throws IOException{
		double counts = model.getCount(chr);
		if(counts < 0) {
			throw new IllegalArgumentException("Can't get poisson expected number of reads for " + counts + " in window.");
		}
		return windowSize * counts / maskedGenome.getLength(window.getChr());
	}*/
	
	public void callKLPeaks(String outFile,double klCutOff) throws IOException{
		
		this.klCutoff = klCutOff;
		BufferedWriter bw = new BufferedWriter(new FileWriter(outFile));
		BufferedWriter bwV = new BufferedWriter(new FileWriter(outFile+".values"));
		
		TreeSet<Gene> sigWindows = new TreeSet<Gene>();
		
		for(String chr:chrSizes.keySet()){
			
			//double backgroundLambda = backgroundSample.getRefSequenceLambda(chr);
			double rapLambda = rapSample.getRefSequenceLambda(chr);
			
			//logger.info("Background lambda for "+chr+": "+backgroundLambda);
			logger.info("RAP lambda for "+chr+": "+rapLambda);
			
			if(rapLambda ==0.0){
				logger.error(chr +" is not expressed in RAP file");
			}
			else{
				Iterator<? extends Window> iter = maskedGenome.getWindowIterator(chr,this.windowSize,this.step);
				//logger.info("Got the window iterator");
				while(iter.hasNext()){
					//logger.info("Enter the loop.");
					Gene window = new Gene(iter.next());
					
					int rapCount = (int) Math.round(rapSample.getCount(window));
					//TODO: If signal is low, do not consider
					if(rapCount < rapLambda) {
						//logger.info("Count for this window: "+rapCount+" is less than lambda: "+rapLambda);
						continue; 
					}

					int backgroundCount = (int) Math.round(backgroundSample.getCount(window));
					/*if(backgroundCount <= 0) 
						continue;
					*/
					double[] qArr = backgroundSample.getCountsPerPosition(window);
					double[] pArr = rapSample.getCountsPerPosition(window);
					double[] minMax=minMax(pArr,qArr);
					
					EmpiricalDistribution pDist = new EmpiricalDistribution(pArr,100,0.0,minMax[1]);
					EmpiricalDistribution qDist = new EmpiricalDistribution(qArr,100,0.0,minMax[1]);
					
					// ADD PSUEDOCOUNTS
					pDist.addPsuedocounts(1.0);
					qDist.addPsuedocounts(1.0);
					
					//KL DIVERGENCE
					double KLdiv = pDist.KLDivergence(qDist);
					
					bwV.write(window.toUCSC()+"\t"+KLdiv+"\n");
					
					if(KLdiv >= klCutoff) {
						logger.info("KL divergence "+KLdiv+" for "+window.toUCSC()+" is less than cutoff: "+klCutoff);
						sigWindows.add(window);
						bw.write(window.toBED() + "\n");
					}
				}
			}
		}
		
		bw.close();
		bwV.close();
	}
	/**
	 * Based on prussell's function in PairedSampleCoverage
	 * Compute Skellam P-value of read counts in a region given the two parameters
	 * @param backgroundLambda Poisson lambda for background sample
	 * @param signalLambda Poisson lambda for signal sample
	 * @param backgroundCount The count for background sample
	 * @param signalCount The count for signal sample
	 * @return The probability under the null hypothesis of observing a greater difference
	 * @throws IOException
	 */
	private double getSkellamPValue(double backgroundLambda, double signalLambda, int backgroundCount, int signalCount){
		
		return Distribution.skellamRightTail(signalCount - backgroundCount, signalLambda, backgroundLambda);
	
	}
	
	

	/**
	 * Returns the coordinate space for this peak caller
	 */
	public CoordinateSpace getCoordinateSpace() {
		return maskedGenome;
	}
	

	/**
	 * This function returns the min and max values (combined) in two arrays
	 * @param arr1
	 * @param arr2
	 * @return
	 */
	private double[] minMax(double[] arr1, double[] arr2){
		double min=Double.MAX_VALUE;
  		double max=-Double.MAX_VALUE;
  		
  		for(int i=0; i<arr1.length; i++){
  			min=Math.min(arr1[i], min);
  			max=Math.max(arr1[i], max);
  		}
  		for(int i=0; i<arr2.length; i++){
  			min=Math.min(arr2[i], min);
  			max=Math.max(arr2[i], max);
  		}
  		
  		double[] minMax={min, max};
  		return minMax;
	}
	/**
	 * 
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException{
		
		/*
		 * Gives a log4j error. Check later.
		 */
		Globals.setHeadless(true);

		/*
		 * @param for ArgumentMap - size, usage, default task
		 * argMap maps the command line arguments to the respective parameters
		 */
		ArgumentMap argMap = CLUtil.getParameters(args,usage,"full");
		RAPPeakCaller dummy = new RAPPeakCaller(argMap.get("background"),argMap.get("rap"),argMap.get("sizes"),argMap.get("maskedRegions"),argMap.getInteger("window"),argMap.getInteger("step"));
		//dummy.callSkellamPeaks(argMap.getOutput(),argMap.getDouble("skellam"),argMap.getDouble("scan"));
		dummy.callKLPeaks(argMap.getOutput(),argMap.getDouble("kl"));
	}


}
