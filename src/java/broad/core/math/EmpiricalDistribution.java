package broad.core.math;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import umontreal.iro.lecuyer.probdist.KolmogorovSmirnovDist;

import org.apache.commons.math3.stat.inference.ChiSquareTest;
import org.apache.log4j.Logger;

import Jama.Matrix;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;

public class EmpiricalDistribution {
	private double max;
	private double min;
	private double[][] intervals;
	private long [] intervalDataNumber;
	ChiSquareTest chi = new ChiSquareTest();
	boolean includeOutOfRange = false;
	
	/**
	 * All the data values that contribute to the bins
	 * Subject to same filters like min and max
	 */
	private Collection<Double> allDataValues;
	
	static Logger logger = Logger.getLogger(EmpiricalDistribution.class.getName());
	
	public EmpiricalDistribution(int numberOfBins, double min, double max) {
		this.min = min;
		this.max = max;
		this.allDataValues = new ArrayList<Double>();
		init(numberOfBins);
	}
	
	public EmpiricalDistribution(int numberOfBins, double min, double max, boolean includeOutOfRange) {
		this(numberOfBins, min, max);
		this.allDataValues = new ArrayList<Double>();
		this.includeOutOfRange = includeOutOfRange;
	}
	
	public EmpiricalDistribution(Collection<Double> vals, int numberOfBins) {
 	   double[] minMax=minMax(vals);
 	  this.allDataValues = new ArrayList<Double>();
        this.min = minMax[0];
        this.max = minMax[1];

        double binSize = ( max - min)/(double)numberOfBins ;
        //System.out.println("binSize: " + binSize);
        intervals = new double[numberOfBins][2];
        intervalDataNumber = new long[numberOfBins];

        //System.err.println("binsize " + binSize);
        for(int i = 0; i < numberOfBins; i++) {
                intervals[i][0] = min + i*binSize;
                intervals[i][1] = min + (i+1)*binSize;
                //System.out.println("Interval " + intervals[i][0] + "-" + intervals[i][1]);
                intervalDataNumber[i] = 0;
        }
		this.allDataValues = new ArrayList<Double>();
        addAll(vals);
	}
	
	public double getMax() {return max;}
	public double getMin() {return min;}
	public int getBinNumber() {return intervalDataNumber.length;}
	
	public EmpiricalDistribution(Collection<Double> vals) {
		this(vals, 200);
	}
	
	public EmpiricalDistribution(Collection<Double> vals, int numberOfBins, boolean absVal) {
		
		if(absVal){
			vals=abs(vals);
		}
		
		double[] minMax=minMax(vals);
		
		this.min = minMax[0];
		this.max = minMax[1];
		
		double binSize = ( max - min)/(double)numberOfBins ;
		//System.out.println("binSize: " + binSize);
		intervals = new double[numberOfBins][2];
		intervalDataNumber = new long[numberOfBins];
		
		//System.err.println("binsize " + binSize);
		for(int i = 0; i < numberOfBins; i++) {
			intervals[i][0] = min + i*binSize;
			intervals[i][1] = min + (i+1)*binSize;
			//System.out.println("Interval " + intervals[i][0] + "-" + intervals[i][1]);
			intervalDataNumber[i] = 0;
		}
		this.allDataValues = new ArrayList<Double>();
		addAll(vals);	
	}
	
	private Collection<Double> abs(Collection<Double> vals){
		Collection<Double> rtrn=new ArrayList<Double>();
		
		for(Double val: vals){
			rtrn.add(Math.abs(val));
		}
		
		return rtrn;
	}
	
    public EmpiricalDistribution(double[] vals, int numberOfBins) {
 	   double[] minMax=minMax(vals);
 	   
        this.min = minMax[0];
        this.max = minMax[1];
        this.allDataValues = new ArrayList<Double>();
        init(numberOfBins);
        addAll(vals);
    }
    
    public EmpiricalDistribution(List<Double> vals, int numberOfBins) {
  	   double[] minMax=minMax(vals);
  	   
         this.min = minMax[0];
         this.max = minMax[1];
         this.allDataValues = new ArrayList<Double>();
         init(numberOfBins);
         addAll(vals);
     }
    
    public EmpiricalDistribution(double[] vals, int numberOfBins, double min, double max) {
   	   
    	this.min = min;
		this.max = max;
		this.allDataValues = new ArrayList<Double>();
		init(numberOfBins);
        addAll(vals);
      }
    
    /**
     * This function will calculate the significance of the chi square statistic between this and other empirical distribution.
     * @param other
     * @return
     * @throws IllegalArgumentException
     * @throws MathException
     */
    public double testGoodnessOfFit(EmpiricalDistribution other) throws IllegalArgumentException {
    	if(other.intervalDataNumber.length != intervalDataNumber.length) {
    		throw new IllegalArgumentException("Distributions have different number of bins: " +intervalDataNumber.length +" vs " +other.intervalDataNumber.length);
    	}
    	
    	return chi.chiSquareTestDataSetsComparison(intervalDataNumber, other.intervalDataNumber);
    }

    /**
     * This function will calculate the chi square distance between this and other empirical distribution.
     * A psuedocount of 1 is added to both distributions before calculating the statistic otherwise zeros will return a MathException
     * @author skadri
     * @param other
     * @return
     * @throws IllegalArgumentException
     * @throws MathException
     */
    public double chiSquareDistance(EmpiricalDistribution other) throws IllegalArgumentException, IOException {
    	if(other.intervalDataNumber.length != intervalDataNumber.length) {
    		throw new IllegalArgumentException("Distributions have different number of bins: " +intervalDataNumber.length +" vs " +other.intervalDataNumber.length);
    	}
    	double sum = 0.0;
    	for(int i = 0; i < intervals.length; i++) {
    		sum += (double)(Math.pow((double)(intervalDataNumber[i]-other.intervalDataNumber[i]),2.0))/(double)(intervalDataNumber[i]+other.intervalDataNumber[i]);
    	}
    	sum = sum/2.0;
    	//return chi.chiSquareTestDataSetsComparison(intervalDataNumber, other.intervalDataNumber);
    	return sum;
    }
    
    /**
     * This function will calculate the KULLBACK-LIEBLER DIVERGENCE between this and other empirical distribution.
     * MAKE SURE THAT GIBB'S INEQUALITY HOLDS. FOR THIS, ADD 1.0 AS PSUEDOCOUNTS (USE addPsuedocounts() FUNCTION)
     * @param other
     * @return
     * @throws IllegalArgumentException
     * @throws MathException
     */
    public double KLDivergence(EmpiricalDistribution other) throws IllegalArgumentException {
    	/*
    	 * All inputs must have same dimension.
    	 */
    	if(other.intervalDataNumber.length != intervalDataNumber.length) {
    		throw new IllegalArgumentException("Distributions have different number of bins: " +intervalDataNumber.length +" vs " +other.intervalDataNumber.length);
    	}

    	/*
    	 * Check probabilities sum to 1 +/- 0.0001
    	 */
    	double thisSum = 0.0;
    	double otherSum = 0.0;
    	for(int i=0;i<intervals.length;i++){
    		thisSum += getRawDensity(i);
    		otherSum += getRawDensity(i);
    	}
    	if((thisSum-1)>0.0001 || (otherSum-1)>0.0001){
    		System.err.println("Probabilities don't sum to 1.0: "+thisSum+" "+otherSum);
    		return -1.0;
    	}
    	else{
	    	double sum = 0.0;
	    	double sum1 = 0.0;
	    	double sum2 = 0.0;
		    for(int i = 0; i < this.intervals.length; i++) {
		   		//Q(x)>0 for all x whenever P(x) >0. 0ln0 is interpreted as 0
		    	if(this.getRawDensity(i)==0.0){
		    		continue;
		    	}
		    	if(other.getRawDensity(i)==0.0){
		    		continue;
		    	}
		    	sum1 += ((double)this.getRawDensity(i))*(Math.log((double)this.getRawDensity(i)))/Math.log(2.0);
		    	sum2 += ((double)this.getRawDensity(i))*(Math.log((double)other.getRawDensity(i)))/Math.log(2.0);
		    	
    			/*double ratio = (Math.log((double)this.getRawDensity(i))) - (Math.log((double)other.getRawDensity(i)));
	    		if(ratio==0.0){
	 //   			System.err.println("Log Ratio is 0.0"+ (Math.log((double)getRawDensity(i)))+" "+(Math.log((double)other.getRawDensity(i))));
	    			sum += 0.0;
	    		}
	    		else{
	    			sum += (((double)this.getRawDensity(i))*ratio)/Math.log(2);	
	    		}*/
	    	}
	//	   	if(sum1<sum2){
	/*	   	if(writeEmpirical){
		   		bw.write("p(x)\t");
		   		for(int i = 0; i < this.intervals.length; i++) {
		   			bw.write(new Double(this.getRawDensity(i)).toString()+"\t");
		   		}
		   		bw.newLine();
		   		bw.write("q(x)\t");
		   		for(int i = 0; i < other.intervals.length; i++) {
		   			bw.write(new Double(other.getRawDensity(i)).toString()+"\t");
		   		}
		   		bw.newLine();
		   //		System.err.println("KLD is negative");
		   	}*/
		   	sum = sum1-sum2;
		   	return sum;
	    }
    }
    
    public double KLDivergenceSym(EmpiricalDistribution other) throws IllegalArgumentException, IOException {
    	
    	//addPsuedocounts(1.0);
    	//other.addPsuedocounts(1.0);
    	double sum = this.KLDivergence(other);
    	sum += other.KLDivergence(this);
    	
    	return sum;
    }
    
	private void init(int numberOfBins) {
		double binSize = ( max - min)/(double)numberOfBins ;
        //System.out.println("binSize: " + binSize);
        intervals = new double[numberOfBins][2];
        intervalDataNumber = new long[numberOfBins];

        //System.err.println("binsize " + binSize);
        for(int i = 0; i < numberOfBins; i++) {
                intervals[i][0] = min + i*binSize;
                intervals[i][1] = min + (i+1)*binSize;
                //System.out.println("Interval " + intervals[i][0] + "-" + intervals[i][1]);
                intervalDataNumber[i] = 0;
        }
	}
	
	public EmpiricalDistribution(InputStream is) throws IOException {
		this.allDataValues = new ArrayList<Double>();
		BufferedReader br = new BufferedReader(new InputStreamReader(is));
		String line = null;
		LinkedHashMap<Double,Integer> observations = new LinkedHashMap<Double, Integer>();
		while((line = br.readLine() ) != null) {
			if(line.startsWith("#") || line.trim().length() == 0) {
				continue;
			}
			String [] lineInfo = line.split("\t");
			observations.put(Double.parseDouble(lineInfo[0]), Integer.parseInt(lineInfo[1]));
		}
		
		intervals = new double[observations.size()][2];
		intervalDataNumber = new long[observations.size()];
		
		Iterator<Double> midPointIt = observations.keySet().iterator();
		int i = 0;
		double midPoint = 0;
		while(midPointIt.hasNext()) {
			midPoint = midPointIt.next();
			intervalDataNumber[i] = observations.get(midPoint);
			intervals[i][0] = midPoint;
			
			if(i > 0) {
				double dist = midPoint - intervals[i-1][0];
				intervals[i-1][1] =  intervals[i-1][0] + dist/2d;
				intervals[i-1][0] = intervals[i-1][1] - dist;
			}
			i++;
		}
		
		intervals[intervals.length - 1][0] = intervals[intervals.length - 2][1];
		intervals[intervals.length - 1][1] = (midPoint - intervals[intervals.length - 1][0] ) +  midPoint;
	}
	
	public EmpiricalDistribution(File source) throws IOException {
		this(new FileInputStream(source));
	}
	
	public EmpiricalDistribution(Matrix permutations, int numberOfBins) {
		double[] minMax=minMax(permutations);
	  	   
         this.min = minMax[0];
         this.max = minMax[1];
         this.allDataValues = new ArrayList<Double>();
         init(numberOfBins);
         addAll(permutations);
	}

	public void add(double observation) {
		this.allDataValues.add(Double.valueOf(observation));
		if(observation > max) {
			if (includeOutOfRange) {
				intervalDataNumber[intervals.length-1]++;
			} else {
				logger.trace("Observation "+observation + " is too large (max="+max+")");
				return;
			}
		} else if (observation < min) {
			if (includeOutOfRange) {
				intervalDataNumber[0]++;
			} else {
				//System.err.println("Observation "+observation + " is too small (min="+min+")");
				return;
			}
		}
		
		for(int i = 0; i < intervals.length; i++) {
			if((intervals[i][0] <= observation && intervals[i][1] > observation) || (intervals[i][1] == this.max && observation == this.max)) {
				intervalDataNumber[i]++;
				//logger.trace("Added obs " + observation + " so far we have " + intervalDataNumber[i]);
				break;
			}
		}

	}
	
	public double leftIntersect(EmpiricalDistribution other) {
		if(intervalDataNumber.length != other.getIntervalDataNumber().length) {
			throw new IllegalArgumentException("Distributions must have the same number of bins");
		}
		
		if(intervals[0][0] != other.intervals[0][0] || intervals[intervals.length - 1][1] != other.intervals[intervals.length - 1][1]) {
			throw new IllegalArgumentException("Distributions must have the same bins");
		}

		double diff = getDensity(0) - other.getDensity(0);
		
		double intersect = (intervals[0][1] + intervals[0][0])/(double)2;
		
		for(int i = 1; i < intervalDataNumber.length; i++) {
			
			double intervalDiff = getDensity(i) - other.getDensity(i);
			if(diff == 0) {
				diff = intervalDiff;
			} else 	if(diff * intervalDiff < 0) {
				intersect =  (intervals[i][1] + intervals[i][0])/(double)2;
				break;
			}
		}
		
		return intersect;
	}
	
	public int getBin(double value) {
		return (int) Math.floor( (value - intervals[0][0]) * intervalDataNumber.length/(intervals[intervals.length - 1][1] - intervals[0][0]));
	}
	
	public double getBinMidPoint(int bin) {
		return (intervals[bin][0] + intervals[bin][1])/(double)2;
	}
	
	public double getBinStart(int bin) {
		return intervals[bin][0];
	}
	
	public double getBinEnd(int bin) {
		return intervals[bin][1];
	}

	
	public double rightIntersect(EmpiricalDistribution other) {
		if(intervalDataNumber.length != other.getIntervalDataNumber().length) {
			throw new IllegalArgumentException("Distributions must have the same number of bins");
		}
		
		if(intervals[0][0] != other.intervals[0][0] || intervals[intervals.length - 1][1] != other.intervals[intervals.length - 1][1]) {
			throw new IllegalArgumentException("Distributions must have the same bins");
		}

		double diff = getDensity(intervals.length - 1) - other.getDensity(intervals.length - 1);
		
		double intersect = (intervals[intervals.length - 1][1] + intervals[intervals.length - 1][0])/(double)2;
		
		for(int i = intervals.length - 2; i >= 0; i--) {
			double intervalDiff = getDensity(i) - other.getDensity(i);
			if(diff == 0) {
				diff = intervalDiff;
			} else	if(diff * intervalDiff < 0) {
				intersect =  (intervals[i][1] + intervals[i][0])/(double)2;
				break;
			}
		}
		
		return intersect;
	}
	
	public double getFirstIntersectLeftFromMedian(EmpiricalDistribution other) {
		if(intervalDataNumber.length != other.getIntervalDataNumber().length) {
			throw new IllegalArgumentException("Distributions must have the same number of bins");
		}
		
		if(intervals[0][0] != other.intervals[0][0] || intervals[intervals.length - 1][1] != other.intervals[intervals.length - 1][1]) {
			throw new IllegalArgumentException("Distributions must have the same bins");
		}

		int medianBin = getMedianBin();
		double diff = getDensity(medianBin) - other.getDensity(medianBin);
		
		double intersect = (intervals[0][1] + intervals[0][0])/(double)2;
		
		for(int i = medianBin; i >= 0; i--) {
			
			double intervalDiff = getDensity(i) - other.getDensity(i);
			if(diff * intervalDiff < 0) {
				intersect =  (intervals[i][1] + intervals[i][0])/(double)2;
				break;
			}
		}
		
		return intersect;
	}
	
	public double getFirstIntersectRightFromMedian(EmpiricalDistribution other) {
		if(intervalDataNumber.length != other.getIntervalDataNumber().length) {
			throw new IllegalArgumentException("Distributions must have the same number of bins");
		}
		
		if(intervals[0][0] != other.intervals[0][0] || intervals[intervals.length - 1][1] != other.intervals[intervals.length - 1][1]) {
			throw new IllegalArgumentException("Distributions must have the same bins");
		}

		int medianBin = getMedianBin();
		double diff = getDensity(medianBin) - other.getDensity(medianBin);
		
		double intersect = (intervals[intervals.length - 1][1] + intervals[intervals.length - 1][0])/(double)2;
		
		for(int i = medianBin; i < intervals.length; i++) {
			double intervalDiff = getDensity(i) - other.getDensity(i);
			if(diff * intervalDiff < 0) {
				intersect =  (intervals[i][1] + intervals[i][0])/(double)2;
				break;
			}
		}
		
		return intersect;
	}
	
	public int getTotalObservations() {
		int total = 0;
		for(int i = 0; i < intervalDataNumber.length; i++) {
			total += intervalDataNumber[i];
		}
		return total;
	}
	
	public void ensureNoEmptyBins() {
		for(int i = 0; i < intervalDataNumber.length; i++) {
			if(intervalDataNumber[i] == 0) {
				intervalDataNumber[i] = 1;
			}
		}
	}
	
	public void addPsuedocounts(double psuedocount) {
		for(int i = 0; i < intervalDataNumber.length; i++) {
				intervalDataNumber[i] += psuedocount;
		}
	}
	public double getDensity(int bin) {
		int total = getTotalObservations(); //must store this as a class field.
		return intervalDataNumber[bin]/((double)total * (intervals[bin][1] - intervals[bin][0]));
	}
	
	/**
	 * Get the total number of observations in the bin
	 * @param bin Bin number
	 * @return
	 */
	public double getHistogram(int bin) {
		return intervalDataNumber[bin];
	}
	
	/**
	 * Get distribution of number of observations per bin
	 * @return Distribution of number of observations per bin
	 */
	public EmpiricalDistribution getDistributionOfHistogramValues() {
		Collection<Double> counts = new ArrayList<Double>();
		for(int i = 0; i < intervalDataNumber.length; i++) {
			counts.add(Double.valueOf(getHistogram(i)));
		}
		return new EmpiricalDistribution(counts);
	}
	
	public double getRawDensity(int bin) {
		int total = getTotalObservations(); //must store this as a class field.
		return intervalDataNumber[bin]/((double)total);
	}
	
	public double getPercentOfValuesBetween(double leftBound, double rightBound) {
		int total = getTotalObservations();
		int totalInRange = 0;
		for(int i = 0; i < intervalDataNumber.length; i++) {
			if(intervals[i][0] > leftBound && intervals[i][1] < rightBound) {
				totalInRange += intervalDataNumber[i];
			}
		}
		
		return totalInRange/(double)total;
	}
	
	public double getPercentOfValuesLargerThan(double value) {
		int total = getTotalObservations();
		if(total == 0) {
			return 0;
		}
		int totalInRange = 0;
		for(int i = 0; i < intervalDataNumber.length; i++) {
			if(intervals[i][0] > value) {
				totalInRange += intervalDataNumber[i];
			}
		}
		
		return totalInRange/(double)total;		
	}
	
	public double getMean() {
		
		double mean = 0;
		double total = getTotalObservations();		
		for(int i = 0; i < intervalDataNumber.length; i++) {
			mean += (intervals[i][1] + intervals[i][0])/(double)2 * intervalDataNumber[i]/total;
		}
		return mean;
	}
		
	public double getStandardDeviation() {
		double mean = getMean();
		double total = getTotalObservations();
		double sd = 0;
		for(int i = 0; i < intervalDataNumber.length; i++) {
			sd += Math.pow((intervals[i][1] + intervals[i][0])/(double)2-mean,2) * intervalDataNumber[i];
		}
		return Math.sqrt(sd/total);
	}
	
	private int getMedianBin() {
		int bin = 0;
		long total = getTotalObservations();
		long totalSoFar = 0;
		//System.out.println("total " + total + " so far " + totalSoFar);
		while(totalSoFar < total/2) {
			totalSoFar += intervalDataNumber[bin++];
			//System.out.println("total " + total + " so far " + totalSoFar);
		}
		return bin;
	}
	
	/**
	 * Get median of all data values including out of range
	 * @return Median of all individual values including out of range or -99 if no data
	 */
	public double getMedianOfAllDataValues() {
		if(allDataValues.isEmpty()) {
			throw new IllegalStateException("No data");
		}
		return Statistics.median(allDataValues);
	}
	
	//This is not efficient, probably should cache this data?
	public double getQuantileOfBins(double probability) {
		int total = getTotalObservations();
		int totalSoFar = 0;
		int bin = 0;
		while(bin < intervalDataNumber.length) {
			totalSoFar += intervalDataNumber[bin];
			double pctSoFar = totalSoFar/(double)total;
			if(probability < pctSoFar) {
				break;
			}
			bin++;
		}
		
		return  (intervals[bin][1] + intervals[bin][0])/(double)2;
		
	}
	
	public double getCumulativeProbability(double value) {
		int total = getTotalObservations();
		if(total == 0) {
			return 0;
		}
		int totalInRange = 0;
		for(int i = 0; i < intervalDataNumber.length; i++) {
			if(intervals[i][0] < value) {
				totalInRange += intervalDataNumber[i];
			}
		}
		
		return totalInRange/(double)total;		
	}
	
	public int getNumberOfValuesLargerThan(double value) {
		int totalInRange = 0;
		for(int i = 0; i < intervalDataNumber.length; i++) {
			if(intervals[i][0] >= value) {
				totalInRange += intervalDataNumber[i];
			}
		}
		
		return totalInRange;		
	}
	
    
    /**
     * Compared a list of observed with this empirical distribution using the Kolmogorov Smirnov Test.
     * Will be approximate since we don't store every value in this distribution.
     * 
     * Get signifiance with: 1 - new umontreal.iro.lecuyer.probdist.KolmogorovSmirnovDist().cdf(D, observed.size());
     * @param observed	List of observations
     * @return			D statistic
     * 
     * NOTE: UNTESTED
     */
    public double getKSStatistic(List<Double> observed) {
    	double maxD = 0;
    	double nObs = observed.size();
    	double nExp = getTotalObservations();
    	
    	Collections.sort(observed);
    	
    	double currD = maxD;
    	int expBin = 0, expCount = 0;
    	for (int i = 0; i < nObs; i++) {
    		while (intervals[expBin][0] < observed.get(i)) {
    			expBin++;
    			expCount += intervalDataNumber[expBin];
    		}
    		currD = Math.abs((i+1.0)/nObs - (expCount / nExp));
    		if (currD > maxD)
    			maxD = currD;
    	}
    	return maxD;
    }
	
    
    public double getKSStatistic(double[] observed) {
    	List<Double> o = new ArrayList<Double>();
    	for (int i=0; i<observed.length; i++) 
    		o.add(observed[i]);
    	return getKSStatistic(o);
    }
   
    
    
	public void write(BufferedWriter bw) throws IOException {
		DecimalFormat format = new DecimalFormat("###0.#####");
		int total = getTotalObservations();
		int totalSoFar = 0;
		//System.err.println("total " + total);
		
		String header = "interval_start_inclusive\t";
		header += "interval_end_exclusive\t";
		header += "count\t";
		header += "density\t";
		header += "cumulative_count\t";
		header += "cumulative_density";
		bw.write(header + "\n");
		
		for(int i = 0; i < intervalDataNumber.length; i++) {
			totalSoFar += intervalDataNumber[i];
			//double midPoint = (intervals[i][0] + intervals[i][1])/(double)2;
			//System.out.println("midpoint: " +midPoint);
			//bw.write(format.format(midPoint));
			bw.write(Double.valueOf(intervals[i][0]).toString());
			bw.write("\t");
			bw.write(Double.valueOf(intervals[i][1]).toString());
			bw.write("\t");
			bw.write(String.valueOf(intervalDataNumber[i]));
			bw.write("\t");
			bw.write(format.format(intervalDataNumber[i]/((double)total * (intervals[i][1] - intervals[i][0]))));
			bw.write("\t");
			bw.write(String.valueOf(totalSoFar));
			bw.write("\t");
			bw.write(format.format(totalSoFar/(double)total));
			//System.out.println(buf);
			bw.newLine();
		}

		
	}
	
	public void write(String save) throws IOException {
		BufferedWriter bw=new BufferedWriter(new FileWriter(save));
		write(bw);
		bw.close();
	}

	public void addDistribution(EmpiricalDistribution other) {
		long    []   otherData = other.getIntervalDataNumber();
		
		if(otherData.length != intervalDataNumber.length) {
			throw new IllegalArgumentException("Distributions need to have the same number of bins this has " + intervalDataNumber.length + " other has " + otherData.length);
		}
		
		for(int i = 0; i < otherData.length; i++) {
			intervalDataNumber[i] += otherData[i];
			allDataValues.addAll(other.getAllDataValues());
		}
	}
	
	/**
	 * Get list of all the data values that contributed to the bins
	 * @return List of all data values subject to same filters as the binned counts
	 */
	public Collection<Double> getAllDataValues() {
		return allDataValues;
	}

	private long[] getIntervalDataNumber() {
		return intervalDataNumber;
	}
	
    private double[] minMax(Collection<Double> vals){
   		double min=Double.MAX_VALUE;
   		double max=-Double.MAX_VALUE;
   		
   		for(Double val: vals){
   			min=Math.min(val, min);
   			max=Math.max(val, max);
   		}
   		
   		double[] minMax={min, max};
   		return minMax;
   	}
       
       
       private double[] minMax(double[] vals){
      		double min=Double.MAX_VALUE;
      		double max=-Double.MAX_VALUE;
      		
      		for(int i=0; i<vals.length; i++){
      			min=Math.min(vals[i], min);
      			max=Math.max(vals[i], max);
      		}
      		
      		double[] minMax={min, max};
      		return minMax;
      	}
       
       private double[] minMax(Matrix vals){
     		double min=Double.MAX_VALUE;
     		double max=-Double.MAX_VALUE;
     		
     		for(int i=0; i<vals.getRowDimension(); i++){
     			for(int j=0; j<vals.getColumnDimension(); j++){
     				double val=vals.get(i, j);
     				min=Math.min(val, min);
         			max=Math.max(val, max);
     			}
     		}
     		
     		double[] minMax={min, max};
     		return minMax;
     	}
       
       public void addAll(Collection<Double> values){
    	   for(Double val: values){
    		   add(val);
    	   }
       }
       
       public void addAll(double[] values){
    	   for(int i=0; i<values.length; i++){
    		   add(values[i]);
    	   }
       }
       
       public void addAll(Matrix vals){
    	   for(int i=0; i<vals.getRowDimension(); i++){
    		   for(int j=0; j<vals.getColumnDimension(); j++){
    			   add(vals.get(i,j));
    		   }
    	   }
       }

	public static final String USAGE = "Usage: EmpiricDistribution TASK=<task_num> <task_args>\n" +
	"\tTasks:\n" +
	"\t\t1. Get empirical distribution -col <column to use 0 for the first one> [-in <Input file standard input is default> -out <output file, standard output by default> -separator <tab (\t) is default>]" +
	"\t\t\t -binNum <Number of bins, default is 200> -maxValue <Maximum value to consider> -minValue <Minimum value to consider>" +
	"\n\t2. Combine distributions -indir <directory> -suffix <suffix of distribution files to combine> -out <output file name defalut is standard out>" +
	"\n\t3. Get quantile -probability -in <Input distribution standard input is default>" +
	"\n";
	
	public static void main(String[] args) throws Exception {
		
		ArgumentMap argMap = CLUtil.getParameters(args, USAGE);
		
		if ("1".equals(argMap.getTask())) {	
			BufferedReader br = argMap.getInputReader();
			BufferedWriter bw = argMap.getOutputWriter(); 
			String separator = argMap.containsKey("separator") ? argMap.get("separator") : "\t";
			int col = argMap.getInteger("col");
			int binNum = argMap.containsKey("binNum") ? argMap.getInteger("binNum") : 200;
			float maxVal = argMap.getFloat("maxValue");
			float minVal = argMap.getFloat("minValue");
			
			EmpiricalDistribution dist = new EmpiricalDistribution(binNum, minVal, maxVal);
			
			String line = null;
			while((line = br.readLine())  != null) {
				if(line.startsWith("#") || line.trim().length() == 0) {
					continue;
				}
				String [] lineInfo = line.split(separator);
				double val = Double.parseDouble(lineInfo[col]);
				dist.add(val);
			}
			br.close();
			dist.write(bw);
			bw.close();
		} else if ("2".equals(argMap.getTask())) {
			File dir = new File(argMap.getInputDir());
			final String suffix = argMap.getMandatory("suffix");
			
			
			String [] distributions = dir.list(new FilenameFilter() {

				public boolean accept(File arg0, String arg1) {
					return arg1.endsWith(suffix);
				}
				
			});
			
			if(distributions.length == 0) {
				return;
			}
			
			EmpiricalDistribution combinedDist = new EmpiricalDistribution(new File(dir + "/" + distributions[0]));
			
			for(int i = 1; i < distributions.length; i++) {
				combinedDist.addDistribution(new EmpiricalDistribution(new File(dir + "/" + distributions[i])));
			}
			
			BufferedWriter bw = argMap.getOutputWriter();
			combinedDist.write(bw);
			bw.close();
		} else if ("3".equals(argMap.getTask())) {
			double probability = argMap.getDouble("probability");
			
			InputStream is = argMap.getInputStream();
			EmpiricalDistribution distribution = new EmpiricalDistribution(is);
			is.close();
			
			System.out.println(distribution.getQuantileOfBins(probability));
		}
	}

	/**
	 * Get the distribution as a map
	 * @return Map associating bin lower bound with number of values in the bin
	 */
	public Map<Double, Double> getCountByBinStart() {
		Map<Double, Double> rtrn = new TreeMap<Double, Double>();
		for(int i=0; i < intervals.length; i++) {
			double binStart = intervals[i][0];
			double numValues = intervalDataNumber[i];
			rtrn.put(Double.valueOf(binStart), Double.valueOf(numValues));
		}
		return rtrn;
	}
	
	//This method gets the probability between size and -size (centered around 0)
	public double getProbability(double size) {
		//need to find the mirror of this number
		double avg=this.getMean();
		
		double flipSize=avg-(size-avg);
		
		double lower=Math.min(size, flipSize);
		double upper=Math.max(size, flipSize);
		
		double prob=this.getCumulativeProbability(upper)-this.getCumulativeProbability(lower);
		
		//System.err.println(avg+" "+size+" "+flipSize+" "+prob);
		
		return prob;
	}

	public double getPValue(double score) {
		return getPValue(score, 0);
	}

	
	public double getPValue(double score, double min) {
		double prob= 1 - this.getCumulativeProbability(score);
		if(prob < min){return min;}
		return prob;
	}
	
	/**
	 * Get the number of bins
	 * @return The number of bins
	 */
	public int numBins() {
		return intervalDataNumber.length;
	}
	
	/**
	 * P(x>Score)
	 * @param score
	 * @return
	 */
	public double getPValue2(double score) {
		
		int total = getTotalObservations();
		int totalInRange = 0;
		for(int i = 0; i < intervalDataNumber.length; i++) {
			if(intervals[i][0] >= score) {
				break;
			}
			totalInRange += intervalDataNumber[i];
		}
		
		return (double)1.0-((double)totalInRange/(double)total);
	}
	
	public double getZscore(double value){
		
		double mean = this.getMean();
		double controlVariance = this.getStandardDeviation();
		return (value - mean)*Math.sqrt((double)this.intervalDataNumber.length)/(controlVariance);
	}
	
	
	public static EmpiricalDistribution getEmpiricalMaxStatisticDistribution(Collection<EmpiricalDistribution> distributions) {
		
		Collection<Double> maxVals = new ArrayList<Double>();
		for(EmpiricalDistribution dist : distributions) {
			maxVals.add(Double.valueOf(dist.getMax()));
		}
		
		return new EmpiricalDistribution(maxVals, 1000);
		
	}

	
}
