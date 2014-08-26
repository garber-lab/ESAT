package broad.core.math;

import java.util.*;

/**
 * @author prussell
 * Mann Whitney U test
 */
public class MannWhitney {

	double z;
	double pValue;
	private static String SAMPLE_1_IDENTIFIER = "sample1";
	private static String SAMPLE_2_IDENTIFIER = "sample2";
	
	/**
	 * Instantiate with arrays of measurements
	 * @param sample1measurements Sample 1 measurements
	 * @param sample2measurements Sample 2 measurements
	 */
	public MannWhitney(double[] sample1measurements, double[] sample2measurements){
		TreeMap<Double, List<String>> map=new TreeMap<Double, List<String>>();
		double[] all = new double[sample1measurements.length + sample2measurements.length];
		for(int i=0; i < sample1measurements.length; i++) {
			all[i]=sample1measurements[i];
			if(map.containsKey(Double.valueOf(all[i]))) {
				List<String> list=map.get(Double.valueOf(all[i]));
				list.add(SAMPLE_1_IDENTIFIER);
			}
			else{
				List<String> list=new ArrayList<String>();
				list.add(SAMPLE_1_IDENTIFIER);
				map.put(Double.valueOf(all[i]), list);
			}
		}
		for(int i = sample1measurements.length; i < all.length; i++){
			all[i] = sample2measurements[i-sample1measurements.length];
			if(map.containsKey(Double.valueOf(all[i]))){
				List<String> list=map.get(Double.valueOf(all[i]));
				list.add(SAMPLE_2_IDENTIFIER);
			}
			else{
				List<String> list=new ArrayList<String>();
				list.add(SAMPLE_2_IDENTIFIER);
				map.put(Double.valueOf(all[i]), list);
			}
		}
		List<Double>[] measurementsRanked=rankOrder(map);
		z=calculateZ(measurementsRanked, sample1measurements, sample2measurements);
		pValue=calculatePvalue(z);
	}
	
	/**
	 * Instantiate with lists of measurements
	 * @param sample1measurements Sample 1 measurements
	 * @param sample2measurements Sample 2 measurements
	 */
	public MannWhitney(ArrayList<Double> sample1measurements, ArrayList<Double> sample2measurements) {
		TreeMap<Double, List<String>> map = new TreeMap<Double, List<String>>();
		double[] all = new double[sample1measurements.size() + sample2measurements.size()];
		for(int i=0; i < sample1measurements.size(); i++){
			all[i] = sample1measurements.get(i).doubleValue();
			if(map.containsKey(Double.valueOf(all[i]))){
				List<String> list= map.get(Double.valueOf(all[i]));
				list.add(SAMPLE_1_IDENTIFIER);
			}
			else{
				List<String> list=new ArrayList<String>();
				list.add(SAMPLE_1_IDENTIFIER);
				map.put(Double.valueOf(all[i]), list);
			}
		}
		for(int i=sample1measurements.size(); i<all.length; i++){
			all[i]=sample2measurements.get(i-sample1measurements.size()).doubleValue();
			if(map.containsKey(Double.valueOf(all[i]))){
				List<String> list= map.get(Double.valueOf(all[i]));
				list.add(SAMPLE_2_IDENTIFIER);
			}
			else{
				List<String> list=new ArrayList<String>();
				list.add(SAMPLE_2_IDENTIFIER);
				map.put(Double.valueOf(all[i]), list);
			}
		}
		List<Double>[] measurementsRanked=rankOrder(map);
		z=calculateZ(measurementsRanked, sample1measurements, sample2measurements);
		pValue=calculatePvalue(z);
	}

	/**
	 * Calculate P value assuming a standard normal distribution of test statistic under the null hypothesis
	 * @param statistic The test statistic
	 * @return The P value
	 */
	private static double calculatePvalue(double statistic){
		cern.jet.random.Normal norm=new cern.jet.random.Normal(0,1, new cern.jet.random.engine.DRand());
		double cdf=norm.cdf(statistic);
		return Math.min(1, Math.min((1-cdf), cdf)*2);
	}
	
	/**
	 * Get the test statistic
	 * @return The test statistic
	 */
	public double getZ() {
		return z;
	}
	
	/**
	 * Get P value
	 * @return P value
	 */
	public double getPvalue() {
		return pValue;
	}
	
	/**
	 * Calculate test statistic based on arrays of measurements
	 * @param measurementsRanked Combined measurements ranked
	 * @param sample1measurements Sample 1 measurements
	 * @param sample2measurements Sample 2 measurements
	 * @return Z statistic
	 */
	private static double calculateZ(List<Double>[] measurementsRanked, double[] sample1measurements, double[] sample2measurements){
		double U=Statistics.sum(measurementsRanked[0]);
		double mu=(sample1measurements.length*(sample1measurements.length+sample2measurements.length+1))/2.0;
		double var=((sample1measurements.length*sample2measurements.length)*(sample1measurements.length+sample2measurements.length+1))/12.0;
		double sigma=Math.sqrt(var);
		double Z=(U-mu)/sigma;
		return Z;
	}
	
	/**
	 * Calculate test statistic based on lists of measurements
	 * @param measurementsRanked Combined measurements ranked
	 * @param sample1measurements Sample 1 measurements
	 * @param sample2measurements Sample 2 measurements
	 * @return Z statistic
	 */
	private static double calculateZ(List<Double>[] measurementsRanked, List<Double> sample1measurements, List<Double> sample2measurements) {
		double U=Statistics.sum(measurementsRanked[0]);
		double mu=(sample1measurements.size()*(sample1measurements.size()+sample2measurements.size()+1))/2.0;
		double var=((sample1measurements.size()*sample2measurements.size())*(sample1measurements.size()+sample2measurements.size()+1))/12.0;
		double sigma=Math.sqrt(var);
		double Z=(U-mu)/sigma;
		return Z;
	}
	
	/**
	 * Get the overall ranks of measurements from both samples
	 * 1 = low, N = high
	 * @param measurements Map of each measurement to the set of samples it belongs to (one or both samples)
	 * @return List whose first element is the list of sample 1 ranks and whose second element is the list of sample 2 ranks
	 */
	protected static List<Double>[] rankOrder(TreeMap<Double, List<String>> measurements){
		ArrayList<Double> sample1ranks = new ArrayList<Double>();
		ArrayList<Double> sample2ranks = new ArrayList<Double>();
		double count=0;
		for(Double key: measurements.keySet()){
			// List of samples the value belongs to (either one or two samples)
			List<String> list= measurements.get(key);
			double end = count + list.size();
			// Rank accounting for ties
			double rank = rank(count, list);
			count=end;
			for(String str: list){
				if(str.equalsIgnoreCase(SAMPLE_1_IDENTIFIER)){
					sample1ranks.add(Double.valueOf(rank));
				}
				if(str.equalsIgnoreCase(SAMPLE_2_IDENTIFIER)){
					sample2ranks.add(Double.valueOf(rank));
				}
			}
			
		}
		//System.err.println(x);
		List[] rtrn={sample1ranks,sample2ranks};
		return rtrn;
	}
	

	protected static double rank(double count, List<String> list){
		// For tie (duplicated) values, returns the mean of the ranks of the duplicated values: 
		// for first rank r and 4 duplicated values it should return  r+ (r+1) + (r+2) + (r+3) = 4*r + 1+2+3 =(in general) N*r + N(N-1)/2
		// But the rank passed (count) is one behind so we add 1.
		double N = list.size();
		return (N*(count+1) + (N*(N - 1))/2)/N;
	}
	
	
}
