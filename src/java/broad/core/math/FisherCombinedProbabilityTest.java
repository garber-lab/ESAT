/**
 * 
 */
package broad.core.math;

import java.util.Collection;
import java.util.TreeSet;

import org.apache.commons.math3.analysis.function.Log;
import org.apache.commons.math3.distribution.ChiSquaredDistribution;

/**
 * @author prussell
 * Fisher's method to combine results of several independent tests bearing upon the same overall hypothesis
 */
public class FisherCombinedProbabilityTest {
	
	private Collection<Double> individualPvals;
	private int numTests;
	
	/**
	 * Instantiate with array of P values
	 * @param pvals The P values from the independent tests
	 */
	public FisherCombinedProbabilityTest(double[] pvals) {
		individualPvals = new TreeSet<Double>();
		for(int i=0; i<pvals.length; i++) {
			individualPvals.add(Double.valueOf(pvals[i]));
		}
		numTests = individualPvals.size();
	}
	
	/**
	 * Instantiate with collection of P values
	 * @param pvals The P values from the independent tests
	 */
	public FisherCombinedProbabilityTest(Collection<Double> pvals) {
		individualPvals = pvals;
		numTests = individualPvals.size();
	}
	
	/**
	 * Get the Fisher test statistic X^2
	 * @return The test statistic
	 */
	public double getTestStatistic() {
		Log log = new Log();
		double rtrn = 0;
		for(Double pval : individualPvals) {
			rtrn += log.value(pval.doubleValue());
		}
		return -2 * rtrn;
	}
	
	/**
	 * Get the number of degrees of freedom of the null chi square distribution
	 * @return The number of degrees of freedom
	 */
	public double getNullDistDegreesOfFreedom() {
		return (double)2 * numTests;
	}
	
	/**
	 * Get the P value of the test statistic
	 * @return The P value with respect to the null chi square distribution
	 */
	public double getCombinedPvalue() {
		
		String message = "Pvals: ";
		for(Double p : individualPvals) {
			message += p.toString() + "\t";
		}
		message += "combined=";
		
		// If one of the P values is exactly zero return zero
		for(Double p : individualPvals) {
			if(p.doubleValue() == 0) {
				message += "0";
				System.err.println(message);
				return 0;
			}
		}
		
		double testStatistic = getTestStatistic();
		ChiSquaredDistribution dist = new ChiSquaredDistribution(getNullDistDegreesOfFreedom());
		double rtrn = 1 - dist.cumulativeProbability(testStatistic);
		message += rtrn;
		System.err.println(message);
		return rtrn;
	}
	
	
}
