package broad.core.math;

import java.io.BufferedReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.commons.collections15.iterators.ArrayIterator;
import org.apache.commons.math3.distribution.HypergeometricDistribution;
import org.apache.commons.math3.special.Erf;


import jsc.independentsamples.SmirnovTest;

import Jama.Matrix;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;

public class Statistics {
	
	/*public static double median(List<? extends Number> list) {
		return quantile(list, 0.5);
	}*/
	
	
	public static double median(Collection<Double> list){
		Object[] array=list.toArray();
		Arrays.sort(array);
		if(array.length%2==0){
			double val1=(Double)array[array.length/2];
			double val2=(Double)array[(array.length/2)-1];
			double median=(val1+val2)/2.0;
			return median;
		}
		else{return (Double)array[array.length/2];}
	}
	
	public static double median(double[] list){
		if(list.length == 1) {
			return list[0];
		}
		ArrayList<Double> withoutNaN = new ArrayList<Double>();
		for(int i=0; i<list.length; i++) {
			if(!Double.isNaN(list[i])) {
				withoutNaN.add(Double.valueOf(list[i]));
			}
		}
		double[] listWithoutNaN = new double[withoutNaN.size()];
		for(int i = 0; i < listWithoutNaN.length; i++) {
			listWithoutNaN[i] = withoutNaN.get(i).doubleValue();
		}
		Arrays.sort(listWithoutNaN);
		if(listWithoutNaN.length%2==0){
			double val1=listWithoutNaN[listWithoutNaN.length/2];
			double val2=listWithoutNaN[(listWithoutNaN.length/2)-1];
			double median=(val1+val2)/2.0;
			return median;
		}
		return listWithoutNaN[listWithoutNaN.length/2];
	}
	
	public static double median(int[] list){
		Arrays.sort(list);
		if(list.length%2==0){
			double val1=list[list.length/2];
			double val2=list[(list.length/2)-1];
			double median=(val1+val2)/2.0;
			return median;
		}
		return list[list.length/2];
	}

	
	/**
	 * Uses the Weighted average method.
	 * @param list - Ordered list of numbers (i.e. short, float, double, etc)
	 * @param pct  - Desired quantile
	 * @return the estimated quantile requested: x s.t. P(X <= x) >= pct
	 */
	public static double quantile(List<? extends Number> list, double pct) {
		if(list.size() == 1) { return list.get(0).doubleValue();}
		if(list.size() == 0) { return 0;}
		if(pct==0.0){
			return (list.get(0).doubleValue());
		}
		if(pct==1.0){
			return list.get(list.size()-1).doubleValue();
		}
		int idx = (int)  Math.floor( pct * (list.size() - 1));
		double reminder = pct * (list.size() - 1) - idx;
		double idxthTerm = list.get(idx).doubleValue();
		double idxNextthTerm = list.get(idx + 1).doubleValue();
		//System.out.println("pct " + pct + " # " + list.size() + " reminder " + reminder + " idx " + idx + " idxthTerm " + idxthTerm + " idxNextthTerm " + idxNextthTerm);
		return  idxthTerm + reminder*(idxNextthTerm - idxthTerm);
	}
	
	public static double quantile(double [] vals, double pct) {
		if(vals.length == 1) { return vals[0];}
		if(vals.length == 0) { return 0;}
		int idx = (int)  Math.floor( pct * (vals.length - 1));
		double reminder = pct * (vals.length - 1) - idx;
		double idxthTerm = vals[idx];
		double idxNextthTerm = vals[idx + 1];
		//System.out.println("pct " + pct + " # " + list.size() + " reminder " + reminder + " idx " + idx + " idxthTerm " + idxthTerm + " idxNextthTerm " + idxNextthTerm);
		return  idxthTerm + reminder*(idxNextthTerm - idxthTerm);
	}
	
	public static double mean(double [] values) {
		double total = 0;
		int counter=0;
		for (int i = 0; i < values.length; i++) {
			if(!new Double(values[i]).equals(Double.NaN)){
				total = total + values[i];
				counter++;
			}
		}
		
		return total/counter;
	}
	
	public static double mean(Matrix data) {
		double total = 0;
		for (int i = 0; i < data.getRowDimension(); i++) {
			for(int j = 0; j < data.getColumnDimension(); j++) {
				total = total + data.get(i,j);	
			}
		}
		
		return total/(double)(data.getRowDimension() * data.getColumnDimension());
	}
	
	public static double mean(Collection<? extends Number>  values) {
		double total = 0;
		int size = values.size();
		for (Number n: values) {
			total = total + n.doubleValue();
		}
		
		return total/(double)size;
	}
	
	
	
	public static double average(Collection<? extends Number> values ) {
		return mean(values);
	}
	
	public static double average(double [] values) {
		return mean(values);
	}
	
	public static double variance(double [] values) {
		return variance(values, mean(values));
	}
	
	public static double variance(double [] values, double mean) {
		double sumOfDeviationSquares = 0;
		for (int i = 0; i < values.length; i++) {
			double deviation = values[i] - mean;
			sumOfDeviationSquares = sumOfDeviationSquares +  deviation * deviation;
		}
		
		return Math.sqrt(sumOfDeviationSquares/(double)(values.length - 1));
	}
	
	public static double variance(List<? extends Number> values) {
		return variance(values, mean(values));
	}
	
	public static double variance(List<? extends Number> values, double mean) {
		double sumOfDeviationSquares = 0;
		int size = values.size();
		for (int i = 0; i < size; i++) {
			double deviation = values.get(i).doubleValue() - mean;
			sumOfDeviationSquares = sumOfDeviationSquares +  deviation * deviation;
		}
		
		return Math.sqrt(sumOfDeviationSquares/(double)(size - 1));
	}
	
	public static double covariance(List<? extends Number> values1, List<? extends Number> values2, double mean1, double mean2) {
		if(values1.size() != values2.size()) {
			throw new IllegalArgumentException("Estimating covariance requries both samples to be of the same size, but sample 1 had length " + values1.size()+", while sample 2 had length " + values2.size());
		}
		
		double sum = 0;
		for(int i = 0; i < values1.size(); i++) {
			sum += (values1.get(i).doubleValue() - mean1)*(values2.get(i).doubleValue() - mean2); //BAD for large numbers .... but good enough for now.
		}
		
		return sum /(double)(values1.size() - 1);
	}
	
	public static double covariance(List<? extends Number> values1, List<? extends Number> values2) {
		return covariance(values1, values2, mean(values1), mean(values2));
	}
	
	public static double pearsonCorrelation(List<? extends Number> values1, List<? extends Number> values2) {
		double cov = covariance(values1, values2);
		double stdv1 = stdev(values1);
		double stdv2 = stdev(values2);
		return cov/(double)(stdv1 * stdv2);
	}
	
	public static double min(Collection<Double> vals){
		double min=Double.MAX_VALUE;
		for(Double val: vals){
			min=Math.min(min, val);
		}
		return min;
	}
	
	public static double max(Collection<Double> vals){
		double max=-Double.MAX_VALUE;
		for(Double val: vals){
			max=Math.max(max, val);
		}
		return max;
	}
	
	public static double min(double[] vals){
		double min=Double.MAX_VALUE;
		for(int i=0; i<vals.length; i++){
			min=Math.min(min, vals[i]);
		}
		return min;
	}
	
	public static double max(double[] vals){
		double max=-Double.MAX_VALUE;
		for(int i=0; i<vals.length; i++){
			max=Math.max(max, vals[i]);
		}
		return max;
	}
	
	public static int minInt(Collection<Integer> vals){
		int min=Integer.MAX_VALUE;
		for(int val: vals){
			min=Math.min(min, val);
		}
		return min;
	}
	
	public static int maxInt(Collection<Integer> vals){
		int max=-Integer.MAX_VALUE;
		for(int val: vals){
			max=Math.max(max, val);
		}
		return max;
	}
	
	public static int  min(int[] vals){
		int min=Integer.MAX_VALUE;
		for(int i=0; i<vals.length; i++){
			min=Math.min(min, vals[i]);
		}
		return min;
	}
	
	public static int max(int[] vals){
		int max=-Integer.MAX_VALUE;
		for(int i=0; i<vals.length; i++){
			max=Math.max(max, vals[i]);
		}
		return max;
	}
	
	/*public static double metaAvg(Collection<List<Double>> data){
		double rtrn=0;
		int count=0;
		for(List<Double> vals: data){
			rtrn+=Statistics.sum(vals);
			count+=vals.size();
		}
		return rtrn/count;
	}*/

	
	public static double metaAvg(Collection<double[]> data){
		double rtrn=0;
		int count=0;
		for(double[] vals: data){
			rtrn+=sum(vals);
			count+=vals.length;
		}
		return rtrn/count;
	}
	
	public static double sum(List<Double> vals){
		double rtrn=0;
		for(int i=0; i<vals.size(); i++){
			rtrn+=(Double)vals.get(i);
		}
		return rtrn;
	}
	
	public static double sum(double[] vals){
		double rtrn=0;
		for(int i=0; i<vals.length; i++){
			rtrn+=vals[i];
		}
		return rtrn;
	}
	
	
	public static double tstat(List<Double> g1, List<Double> g2, double tstat_tuning_param){
		double avg1=average(g1);
		double avg2=average(g2);
		
		double sd1=stdev(g1);
		double sd2=stdev(g2);
		
		double n1=g1.size();
		double n2=g2.size();
		
		
		double S = Math.sqrt((Math.pow(sd1,2)*(n1-1) + (Math.pow(sd2, 2)*(n2-1)))/(n1+n2-2));
		double t = (avg1 - avg2)*Math.sqrt(n1*n2)/((tstat_tuning_param+S)*Math.sqrt(n1+n2));
		return t;
	}
	
	public static double tstat(double[] g1, double[] g2, double tstat_tuning_param){
		double avg1=average(g1);
		double avg2=average(g2);
		
		double sd1=stdev(g1);
		double sd2=stdev(g2);
		
		double n1=g1.length;
		double n2=g2.length;
		
		
		double S = Math.sqrt((Math.pow(sd1,2)*(n1-1) + (Math.pow(sd2, 2)*(n2-1)))/(n1+n2-2));
		double t = (avg1 - avg2)*Math.sqrt(n1*n2)/((tstat_tuning_param+S)*Math.sqrt(n1+n2));
		return t;
	}
	
	public static double diff(double[] g1, double[] g2){
		double avg1=average(g1);
		double avg2=average(g2);
		
		
		return avg2-avg1;
	}
	
	
	
	//computes the unnormalized maxmean as defined by Tibshirani
	public static double maxmean(double[] vals){
		//first compute the sum of all positive values and negative values
		double positiveSum=0;
		double negativeSum=0;
		for(int i=0; i<vals.length; i++){
			if(vals[i]>0){positiveSum+=vals[i];}
			else{negativeSum+=vals[i];}
		}
		
		double positiveMean=positiveSum/vals.length;
		double negativeMean=negativeSum/vals.length;
		
		if(positiveMean > Math.abs(negativeMean)){return positiveMean;}
		else{return negativeMean;}
				
	}
	
	
	public static double pairedTStat(List<Double> g1, List<Double> g2){
		List<Double> dVector=new ArrayList<Double>();
		
		for(int i=0; i<g1.size(); i++){
			double d=g1.get(i)-g2.get(i);
			dVector.add(d);
		}
		
		double avg=average(dVector);
		double stdev=stdev(dVector);
		int n=dVector.size();
		
		double t=(avg*Math.sqrt(n))/stdev;
		return t;
	}
	
	
	public static double pairedTStat(double[] gr1, double[] gr2) {
		List<Double> dVector=new ArrayList<Double>();
		
		for(int i=0; i<gr1.length; i++){
			double d=gr1[i]-gr2[i];
			dVector.add(d);
		}
		
		double avg=average(dVector);
		double stdev=stdev(dVector);
		int n=dVector.size();
		
		double t=(avg*Math.sqrt(n))/stdev;
		return t;
	}	
	
	public static double stdev(Collection<? extends Number> values){
		double avg=average(values);
		double stdev=0;
		for(Number val: values){
			stdev+=Math.pow(val.doubleValue()-avg,2);
		}
		stdev=Math.sqrt(stdev/(values.size()-1));
		return stdev;
	}
	
	public static double stdev(double[] values){
		double avg=average(values);
		double stdev=0;
		for(double val: values){
			stdev+=Math.pow(val-avg,2);
		}
		stdev=Math.sqrt(stdev/(values.length-1));
		return stdev;
	}
	
	public static double pearsonDistance(double[] x, double[] y){
		double sum_sq_x = 0;
		double sum_sq_y = 0;
		double sum_coproduct = 0;
		double mean_x = average(x);
		double mean_y = average(y);
		int N=x.length;
		
		for (int i=0; i<x.length; i++){
		    double sweep = (i + 1.0) / (i+2);
		    double delta_x = x[i] - mean_x;
		    double delta_y = y[i] - mean_y;
		    sum_sq_x += delta_x * delta_x * sweep;
		    sum_sq_y += delta_y * delta_y * sweep;
		    sum_coproduct += delta_x * delta_y * sweep;
		    mean_x += delta_x / (i+2);
		    mean_y += delta_y / (i+2) ;
		}
		double pop_sd_x = Math.sqrt( sum_sq_x / N );
		double pop_sd_y = Math.sqrt( sum_sq_y / N );
		double cov_x_y = sum_coproduct / N;
		double correlation = cov_x_y / (pop_sd_x * pop_sd_y);		
		return correlation;
	}
	
	public static double pearsonDistance(List<Double> xA, List<Double> yA){
		double[] x=new double[xA.size()];
		double[] y=new double[yA.size()];
		
		for(int i=0; i<xA.size(); i++){
			x[i]=xA.get(i);
			y[i]=yA.get(i);
		}
		
		return pearsonDistance(x, y);
	}
	
	public static double spearmanCorrelation(List<Double> xA, List<Double> yA){
		double[] x=new double[xA.size()];
		double[] y=new double[yA.size()];
		
		for(int i=0; i<xA.size(); i++){
			x[i]=xA.get(i);
			y[i]=yA.get(i);
		}
		double[] rX=rank(x, true);
		double[] rY=rank(y, true);
		
		return pearsonDistance(rX,rY);
		
	}
	
	public static double spearmanCorrelation(double [] x, double [] y){

		double[] rX=rank(x, true);
		double[] rY=rank(y, true);
		
		return pearsonDistance(rX,rY);
		
	}
	
	public static double euclideanDistance(double[] x, double[] y){
		if (isNaN(x) || isNaN(y))
			return Double.NaN;
		double sum=0;
		for(int i=0; i<x.length; i++){
			sum+=Math.pow(x[i]-y[i],2);
		}
		return Math.sqrt(sum);
	}
	
	public static double euclideanDistance(List<Double> x, List<Double> y){
		if (isNaN(x) || isNaN(y))
			return Double.NaN;
		double sum=0;
		for(int i=0; i<x.size(); i++){
			sum+=Math.pow(x.get(i)-y.get(i),2);
		}
		return Math.sqrt(sum);
	}
	
	public double manhattanDistance(double[] x, double[] y){
		double sum=0; 
		for(int i=0; i<x.length; i++){
			sum+=Math.abs(x[i]-y[i]);
		}
		return sum;
	}
	
	public static double FDRCorrect(double[] pvals, double alpha){
		Arrays.sort(pvals);
		double g=pvals.length;
		int k=0;
		for(int i=0; i<pvals.length; i++){
			double qi=(i/g)*alpha;
			if(pvals[i]<=qi){k=i;}
		}
		return pvals[k];
	}
	
	public static double FDRCorrect(List<Double> pvals, double alpha){
		Collections.sort(pvals);
		double g=pvals.size();
		int k=0;
		for(int i=0; i<pvals.size(); i++){
			double qi=(i/g)*alpha;
			if(pvals.get(i)<=qi){k=i;}
		}
		return pvals.get(k);
	}
	
	public static double FDR(EmpiricalDistribution observedDist, EmpiricalDistribution[] randomDist, double x){
		double Rk=observedDist.getNumberOfValuesLargerThan(Math.abs(x));
		//double Rk2=countNumGreaterThanThresh(tstats, k);
		
		//System.err.println("vals: "+dist.getPercentOfValuesLargerThan(k)+" "+(1-dist.getCummulativeProbability(k))+" "+Rk1+" "+Rk2);
		
		double[] Vk=new double[randomDist.length];
		for(int i=0; i<randomDist.length; i++){
			Vk[i]=randomDist[i].getNumberOfValuesLargerThan(Math.abs(x));
		}
		
		double EVk=mean(Vk);
		double FDRk=0;
		if(Rk>0){FDRk=EVk/Rk;}
		//System.err.println(FDRk+" "+EVk+" "+Rk);
		return FDRk;		
	}
	

	/**
	 * Adjust a set of P values for multiple testing with the Hochberg method
	 * @param pvals The P values
	 * @return Array of adjusted P values in same order
	 */
	public static double[] hochbergPvalueAdjust(double[] pvals) {
		
		int m = pvals.length;
		if(m == 0) throw new IllegalArgumentException("Array of P-values is empty");
		if(m == 1) return pvals;
		double[] origOrder = Arrays.copyOf(pvals, m);
		Arrays.sort(pvals);
		double[] corrected = new double[m];
		double[] rtrn = new double[m];
		
		TreeSet<Double> pvalsDone = new TreeSet<Double>();
		
		corrected[m-1] = pvals[m-1];
		for(int j=0; j<m; j++) {
			if(origOrder[j] == pvals[m-1]) {
				pvalsDone.add(Double.valueOf(pvals[m-1]));
				rtrn[j] = corrected[m-1];
				//System.err.println("Corrected P-value for " + origOrder[j] + " is " + rtrn[j]);
			}
		}


		for(int i=2; i <= m; i++) {
			corrected[m-i] = Math.min(corrected[m-i+1], i*pvals[m-i]);
			if(pvalsDone.contains(Double.valueOf(pvals[m-i]))) continue;
			for(int j=0; j<m; j++) {
				if(origOrder[j] == pvals[m-i]) {
					rtrn[j] = corrected[m-i];
					if(!pvalsDone.contains(Double.valueOf(pvals[m-i]))) {
						pvalsDone.add(Double.valueOf(pvals[m-i]));
						//System.err.println("Corrected P-value for " + origOrder[j] + " is " + rtrn[j]);
					}
				}
			}
		}

		//System.err.println("Returning:");
		//for(int i=0; i<m; i++) System.err.print(rtrn[i] + "\t");
		//System.err.println();
		
		return rtrn;
		
	}
	
	/**
	 * Computes what is P(X > value) using list of values
	 * @param data 
	 * @param value
	 * @param isSorted true if the data list is sorted, it will avoid sorting here.
	 * @return
	 */
	public static double pvalue (List<Double> data, double value, boolean isSorted) {
		if(!isSorted) {
			Collections.sort(data);
		}
		int idx = 0;
		int size = data.size();
		while(idx < size) {
			if(data.get(idx) >= value) {
				break;
			}
			idx++;
		}
		
		return 1.0 - idx /(double) size;
	}
	
	/**
	 * Computes what is P(X > value) using list of values
	 * @param data 
	 * @param value
	 * @param isSorted true if the data list is sorted, it will avoid sorting here.
	 * @return
	 */
	public static double pvalue (double [] data, double value, boolean isSorted) {
		if(!isSorted) {
			Arrays.sort(data);
		}
		int idx = 0;
		int size = data.length;
		while(idx < size) {
			if(data[idx] >= value) {
				break;
			}
			idx++;
		}
		
		return 1.0 - idx /(double) size;
	}
	
	/**
	 * Computes what is P(X > value) using list of values
	 * @param data 
	 * @param value
	 * @return
	 */
	public static double pvalue (List<Double> data, double value) {
		return pvalue(data, value, false);
	}
	
	public static final String USAGE = "Usage: Statistics TASK=<task_num> <task_args>\n" +
	"\tTasks:\n" +
	"\t\t1. Compute Benjamini-Hochberg FDR adjusted pvalue of a given pvalue column in a file " +
	"\n\t\t\t-in <file (standard input is default)> -fdr <Maximum false discovery rate> -col <Data column, 0 based> [-separator <character used to separate columns, default is tab '\\t>'] " +
	"\n";
	
	/*public static void main(String[] args) throws Exception {
		
		ArgumentMap argMap = CLUtil.getParameters(args, USAGE);
		
		if ("1".equals(argMap.getTask())) {	
			double fdr = argMap.getDouble("fdr");
			int col = argMap.getInteger("col");
			String separator = argMap.containsKey("separator") ? argMap.get("separator") : "\t";
			BufferedReader br = argMap.getInputReader();
			
			String line = null;
			ArrayList<Double> pvals = new ArrayList<Double>();
			System.out.println("reading file");
			while((line = br.readLine()) != null) {
				//System.out.println(line);
				String [] lineInfo = line.split(separator);
				//System.out.println("lineInfo size: " +  lineInfo.length + " col value <"+lineInfo[col] +">");
				pvals.add(Double.parseDouble(lineInfo[col]));
			}
			//System.out.println(" .. done ... starting array");
			double [] valArray = new double [pvals.size()];
			for(int i = 0; i < pvals.size(); i++) {
				valArray[i] = pvals.get(i);
			}
			System.out.println("done. calling benjamini hochberg");
			System.out.println("Desired FDR: " + fdr + " BH adjusted pvalue: " + FDRCorrect(valArray, fdr));
		} else {
			System.err.println(USAGE);
		}
	}*/
	

	public static void main(String[] args) throws Exception {
		
		ArgumentMap argMap = CLUtil.getParameters(args, USAGE);
		
		if ("1".equals(argMap.getTask())) {	
			double fdr = argMap.getDouble("fdr");
			int col = argMap.getInteger("col");
			String separator = argMap.containsKey("separator") ? argMap.get("separator") : "\t";
			BufferedReader br = argMap.getInputReader();
			
			String line = null;
			ArrayList<Double> pvals = new ArrayList<Double>();
			System.out.println("reading file");
			while((line = br.readLine()) != null) {
				//System.out.println(line);
				String [] lineInfo = line.split(separator);
				//System.out.println("lineInfo size: " +  lineInfo.length + " col value <"+lineInfo[col] +">");
				pvals.add(Double.parseDouble(lineInfo[col]));
			}
			//System.out.println(" .. done ... starting array");
			double [] valArray = new double [pvals.size()];
			for(int i = 0; i < pvals.size(); i++) {
				valArray[i] = pvals.get(i);
			}
			System.out.println("done. calling benjamini hochberg");
			System.out.println("Desired FDR: " + fdr + " BH adjusted pvalue: " + FDRCorrect(valArray, fdr));
		} else {
			System.err.println(USAGE);
		}
	}

	public static double zScore(double observation, double[] controlObservations, double fudge_factor) {
		double mean = mean(controlObservations);
		double controlVariance = variance(controlObservations, mean);
		return (observation - mean)*Math.sqrt(controlObservations.length)/(controlVariance + fudge_factor);
	}
	
	public static double zScore(double observation, double mean, double variance, int length) {
		return (observation - mean)*Math.sqrt(length)/(variance);
	}
	
	/**
	 * Two sample T statistic
	 * @param g1 Measurements from group 1
	 * @param g2 Measurements from group 2
	 * @return T statistic (positive iff mean of group 1 is greater than mean of group 2)
	 */
	public static double tstat(List<Double> g1, List<Double> g2){
		return tstat(g1,g2, 3.0);
	}
	
	public static double zScore(double observation, double[] controlObservations) {
		return zScore(observation, controlObservations, 0);
	}
	
	public static double []  computePermutedZScores(double observation, double [] controlObservations, double fudge_factor) {
		double [] randomizedControls = new double[controlObservations.length];
		for(int i = 0; i < controlObservations.length; i++) {
			randomizedControls[i] = controlObservations[i];
		}
		
		double [] randomizedZScores = new double[controlObservations.length];
		for(int i = 0; i < controlObservations.length; i++) {
			double testVal = controlObservations[i];
			randomizedControls[i] = observation;
			randomizedZScores[i] =  zScore(testVal, randomizedControls);
			randomizedControls[i] = testVal;
		}
		
		return randomizedZScores;
	}
	
	public static double []  computePermutedZScores(double observation, double [] controlObservations) {
		return computePermutedZScores(observation, controlObservations, 0);
	}


	//given an array of values covert to a vector of relative ranks within vector
	public static double[] rank(double[] vals, boolean ascending){
		double[] rtrn=new double[vals.length];
		
		//make count map
		Map<Double, Integer> map=new TreeMap<Double, Integer>();
		for(int i=0; i<vals.length; i++){
			int counter=0;
			if(map.containsKey(vals[i])){counter=map.get(vals[i]);}
			counter++;
			map.put(vals[i], counter);
		}
		
		//assign ranks to each number
		Map<Double, Double> rankMap=new TreeMap<Double, Double>(); //original value and rank
		
		if(ascending){
			int rank=0;
			for(Double val: map.keySet()){
				int count=map.get(val);
				ArrayList<Integer> list=new ArrayList<Integer>();
				for(int i=0; i<count; i++){list.add(rank+i);}
				rank+=count;
				rankMap.put(val, average(list));
			}
		}
		else{
			int rank=vals.length-1;
			for(Double val: map.keySet()){
				int count=map.get(val);
				ArrayList<Integer> list=new ArrayList<Integer>();
				for(int i=0; i<count; i++){list.add(rank+i);}
				rank-=count;
				rankMap.put(val, average(list));
			}
		}
				
		for(int i=0; i<rtrn.length; i++){
			rtrn[i]=rankMap.get(vals[i]);
		}
		
		return rtrn;
	}
	
	public static double absFold(Collection<Double> g1, Collection<Double> g2, boolean log){
		double avg1=average(g1);
		double avg2=average(g2);
		
		double fold=Math.max(avg1,avg2)/Math.min(avg1, avg2);
		
		
		if(log){
			double sum=0;
			for(Double val: g1){
				sum+=Math.pow(2, val);
			}
			avg1=sum/g1.size();
			
			sum=0;
			for(Double val: g2){
				sum+=Math.pow(2, val);
			}
			avg2=sum/g2.size();
			
			fold=Math.max(avg1,avg2)/Math.min(avg1, avg2);
		}
		
		
		
		return fold;
	}
	
	public static double absFold(double[] g1, double[] g2, boolean log){
		double avg1=average(g1);
		double avg2=average(g2);
		
		double fold=Math.max(avg1,avg2)/Math.min(avg1, avg2);
		
		if(log){
			double sum=0;
			for(int i=0; i<g1.length; i++){
				sum+=Math.pow(2, g1[i]);
			}
			avg1=sum/g1.length;
			
			sum=0;
			for(int i=0; i<g2.length; i++){
				sum+=Math.pow(2, g2[i]);
			}
			avg2=sum/g2.length;
			fold=Math.max(avg1,avg2)/Math.min(avg1, avg2);
		}
		
		
		return fold;
	}
	
	public static double fold(Collection<Double> g1, Collection<Double> g2, boolean log){
		double avg1=average(g1);
		double avg2=average(g2);
		
		double fold=avg1 / avg2;
		
		
		if(log){
			double sum=0;
			for(Double val: g1){
				sum+=Math.pow(2, val);
			}
			avg1=sum/g1.size();
			
			sum=0;
			for(Double val: g2){
				sum+=Math.pow(2, val);
			}
			avg2=sum/g2.size();
			
			fold=avg1/avg2;
		}
		
		
		
		return fold;
	}
	
	public static double fold(double[] g1, double[] g2, boolean log){
		double avg1=average(g1);
		double avg2=average(g2);
		
		double fold=avg1/ avg2;
		
		if(log){
			double sum=0;
			for(int i=0; i<g1.length; i++){
				sum+=Math.pow(2, g1[i]);
			}
			avg1=sum/g1.length;
			
			sum=0;
			for(int i=0; i<g2.length; i++){
				sum+=Math.pow(2, g2[i]);
			}
			avg2=sum/g2.length;
			fold=avg1/avg2;
		}
		
		
		return fold;
	}
	
	public static double fold(double[] g1, boolean log){
		double max=max(g1);
		double min=min(g1);
		
		double fold=max/min;
		
		if(log){
			fold=Math.pow(2, max-min);
		}
		
		
		return fold;
	}

	public static double medianCollection(Collection<Double> values) {
		return median(new ArrayList(values));
	}

	
	public static double hypergeometric(int numberOfGenesTotal, int diffExpressedGenes, int genesInSet, int diffExpressedGenesInSet) {
		if(genesInSet==0 || diffExpressedGenesInSet==0){return 1.0;}
		
		double sum=0;
		cern.jet.random.HyperGeometric hyper=new cern.jet.random.HyperGeometric(numberOfGenesTotal, diffExpressedGenes, genesInSet, new cern.jet.random.engine.DRand() );
		
		HypergeometricDistribution hyper2=new HypergeometricDistribution(numberOfGenesTotal, diffExpressedGenes, genesInSet); 
		
		//System.err.println(1-hyper2.cumulativeProbability(diffExpressedGenesInSet));
		
		for(int i=diffExpressedGenesInSet; i<Math.min(genesInSet, diffExpressedGenes); i++){
			//System.err.println(i+" "+hyper.pdf(i));
			sum+=hyper.pdf(i);
		}
		
		return sum;	
	}
	
	public static double MAD(double[] array){
		double median=median(array);
		double[] rtrn=new double[array.length];
		for(int i=0; i<array.length; i++){
			rtrn[i]=Math.abs(array[i]-median);
		}
		return median(rtrn);
	}

	public static double maxStdev(Matrix data) {
		double maxStdev=-Double.MAX_VALUE;
		for(int i=0; i<data.getColumnDimension(); i++){
			double[] vals=data.getColumn(i);
			double mad=stdev(vals);
			maxStdev=Math.max(mad, maxStdev);
		}
		return maxStdev;
	}

	public static double maxMAD(Matrix data) {
		double maxMAD=-Double.MAX_VALUE;
		for(int i=0; i<data.getColumnDimension(); i++){
			double[] vals=data.getColumn(i);
			double mad=MAD(vals);
			maxMAD=Math.max(mad, maxMAD);
		}
		return maxMAD;
	}
	
	public static double geometricMean(double[] vals){
		double mean=1;
		
		for(int i=0; i<vals.length; i++){
			mean*=vals[i];
		}
		
		mean=Math.pow(mean, (double)1/vals.length);
		return mean;
	}

	public static double geometricMeanMAD(Matrix data) {
		double mean=1;
		for(int i=0; i<data.getColumnDimension(); i++){
			double[] vals=data.getColumn(i);
			double mad=MAD(vals);
			//System.err.println(mad);
			mean=mean*mad;
		}
		mean=Math.pow(mean, (double)1/data.getColumnDimension());
		return mean;
	}

	public static double anovaFStat(Collection<double[]> lists) {
		ANOVA anova=new ANOVA(lists);
		return anova.F();
	}

	//returns an array with a random permutation of the numbers 1-n
	public static int[] randomPermutation(int n){
	
		int i, cnt, a, b, t;

	    Random rng=new Random(System.currentTimeMillis());
		   int arr[]=new int[n];
		   //set swap times:
		    cnt = Math.abs(rng.nextInt())% n + n;

		   // initialize array
		   for(i=0; i<n; i++)arr[i]=i+1;

		  //swap a,b positions in loop:
		   for(i=0; i<cnt; i++){   
		     a = Math.abs(rng.nextInt()%n);
		     b = Math.abs(rng.nextInt()%n);
		     //swap:
		     t= arr[a];
		     arr[a]=arr[b];
		     arr[b]=t;
		   }
		   return(arr);

		
	}

	public static double [] toDoubleArray(List<Double> lst) {
		
		double [] rtrn=new double[lst.size()];
		for(int i=0; i<lst.size();i++) {rtrn[i]=lst.get(i);}
		return rtrn;
	}

	public static ArrayList<Double> toDoubleArrayList(double[] intervalScores) {
		
		ArrayList<Double> arr= new ArrayList<Double>();
		for (int i=0; i<intervalScores.length; i++) arr.add(new Double(intervalScores[i]));
		return arr;
		
	}

	public static ArrayList<Double> toUniqueSortedList( ArrayList<Double> List) {
		Collections.sort(List);
		ArrayList<Double> unq= new ArrayList<Double>();
		for (int i=0; i<List.size(); i++)
			if (! unq.contains(List.get(i)))unq.add(List.get(i));
		Collections.sort(unq);
		return unq;
	}

	public static Matrix absoluteValue(Matrix permutations) {
		Matrix rtrn=permutations.copy();
		
		for(int i=0; i<rtrn.getRowDimension(); i++){
			for(int j=0; j<rtrn.getColumnDimension(); j++){
				rtrn.set(i, j, Math.abs(permutations.get(i, j)));
			}
		}
		return rtrn;
	}

	public static double[] absoluteValue(double[] observed) {
		double[] rtrn=new double[observed.length];
		
		for(int i=0; i<observed.length; i++){
			rtrn[i]=Math.abs(observed[i]);
		}
		
		return rtrn;
	}

	public static double entropy(double[] p_vec) {
		double res=0;
		for (int i=0; i<p_vec.length; i++){
			double p= p_vec[i];
			if (p !=0)
				res+= p* (Math.log(p)/Math.log(2));
		}
		res=-1*res;
		return res;
	}
	
	/**
	 * computes the hamming distance between two vectors, given a threshold.
	 * The threshold is used to reduce the vectors to a trinary form were, for example
	 * x[i] is set to 1 when x[i]>threshold, x[i] is set to -1 when x[i]< threshold and 0 otherwise
	 * @param x
	 * @param y
	 * @param threshold to use in order to reduce the vectors
	 * @return
	 */
	public static double hamming(double [] x, double [] y, double threshold) {
		short [] tx = threshold(x, threshold);
		short [] ty = threshold(y, threshold);
		
		return hamming(tx, ty);
	}
	
	
	/**
	 * computes the hamming distance between two vectors, given a threshold.
	 * The threshold is used to reduce the vectors to a trinary form were, for example
	 * x.get(i) is set to 1 when x.get(i)>threshold, x.get(i)is set to -1 when x.get(i)< threshold and 0 otherwise
	 * @param x
	 * @param y
	 * @param threshold to use in order to reduce the vectors
	 * @return
	 */
	public static double hamming(List<Double> x, List<Double> y, double threshold) {
		short [] tx = threshold(x, threshold);
		short [] ty = threshold(y, threshold);
		
		return hamming(tx, ty);
	}
	
	private static short [] threshold(double [] x, double threshold) {
		short [] tx = new short[x.length];
		for (int i = 0; i < x.length; i++) {
			if(x[i] > threshold ) {
				tx[i] = 1;
			} else if ( x[i] < threshold) {
				tx[i] = -1;
			} else {
				tx[i] = 0;
			}
		}
		return tx;
	}
	
	private static double hamming (short [] x, short [] y) {
		double d = 0;
		if(x.length != y.length) {
			throw new IllegalArgumentException("Length of vectors is not the same, can't compute Hamming distance ");
		}
		for (int i = 0; i < x.length; i++) {
			if(x[i] != y[i]) {
				d++;
			}
		}
		
		return d;
	}
	
	private static short [] threshold(List<Double> x, double threshold) {
		short [] tx = new short[x.size()];
		for (int i = 0; i < x.size(); i++) {
			if(x.get(i) > threshold ) {
				tx[i] = 1;
			} else if ( x.get(i) < threshold) {
				tx[i] = -1;
			} else {
				tx[i] = 0;
			}
		}
		return tx;
	}

	public static Double JSDist(double[] p1, double[] p2) {
		int len=p1.length;
		double res=Double.NaN;
		if (isNaN(p1) || isNaN(p2))
			return res;
		double [] pm= new double[len];
		for (int i=0; i<len;i++) {pm[i]=(p1[i]+p2[i])/2;}
		double hpm=Statistics.entropy(pm);
		double hp1=Statistics.entropy(p1);
		double hp2=Statistics.entropy(p2);
		res= Math.sqrt(hpm - ((hp1+hp2)/2));
		return res;
	}

	
	public static double[] normalizeToRelativeAbundance(double[] vals) {
		int len=vals.length;
		double[] res= new double[len];
		double sum= 0;
		for (int i=0; i<len;i++) {sum+=vals[i];}
		for (int i=0; i<len;i++) {res[i]=(vals[i]/sum);}
		return res;
	}
	
	private static boolean isNaN(double[] p1) {
		boolean res=false;
		for (int i=0; i<p1.length;i++){
			if (Double.isNaN(p1[i])){
				res=true; break;}
		}
		return res;
	}
	
	private static boolean isNaN(List<Double> p1) {
		boolean res=false;
		for (int i=0; i<p1.size();i++){
			if (Double.isNaN(p1.get(i))){
				res=true; break;}
		}
		return res;
	}
	
	/**
     * Calculate p-value given a Z-score
     *
     * @return double p-value
     */
    public static double zscoreToPvalue(double zscore) {
        zscore = zscore / Math.sqrt(2.0);
        double lPvalue = 0.0;
        lPvalue = Erf.erf(zscore);
        return lPvalue;
    }

    public static double zScore2(double observation, double[] controlObservations) {
		double mean = mean(controlObservations);
		double controlVariance = variance(controlObservations, mean);
		return (observation - mean)/(controlVariance );
	}
    
    
    public static double getKSStatistic(double[] observed, double[] expected) {
    	SmirnovTest test = new SmirnovTest(observed, expected);
    	return test.getStatistic();
    }
    
    /**
     * Jesse Engreitz
     * September 28, 2012
     * @param observed
     * @param expected
     * @return
     * This works ... but better to use the existing library in the above method
    public static double getKSStatistic(List<Double> observed, List<Double> expected) {

    	double maxD = 0;
    	double nObs = observed.size();
    	double nExp = expected.size();
    	
    	Collections.sort(observed);
    	Collections.sort(expected);
    	
    	double currD = maxD;
    	int expCounter = 0;
    	for (int i = 0; i < nObs; i++) {
    		
    		while (expected.get(expCounter) < observed.get(i)) expCounter++;
    		currD = Math.abs(((i+1.0) / nObs) - ((expCounter) / nExp));
    		if (currD > maxD)
    			maxD = currD;
    	}
    	
    	return maxD;
    }
    */

    
	
}
