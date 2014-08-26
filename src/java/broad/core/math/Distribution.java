/**
 * 
 */
package broad.core.math;

import umontreal.iro.lecuyer.probdist.ChiSquareNoncentralDist;


/**
 * @author prussell
 * Statistical distributions
 */
public class Distribution {
	
	
	/**
	 * Compute the pdf of a noncentral chi square distribution
	 * @param x Value of random variable
	 * @param nu Degrees of freedom
	 * @param lambda Noncentrality parameter
	 * @return The pdf of the noncentral chi square distribution with the specified parameters, evaluated at the specified value
	 */
	public static double noncentralChiSquareDensity(double x, double nu, double lambda) {
		
		if(nu <= 0) {throw new IllegalArgumentException("Degrees of freedom must be > 0.");}
		if(lambda <= 0) {throw new IllegalArgumentException("Non-centrality parameter must be > 0.");}
		
		return ChiSquareNoncentralDist.density(nu, lambda, x);
		
	}
		
	/**
	 * Compute the cdf of a noncentral chi square distribution
	 * @param x Value of random variable
	 * @param nu Degrees of freedom
	 * @param lambda Noncentrality parameter
	 * @return The cdf of the noncentral chi square distribution with the specified parameters, evaluated at the specified value
	 */	
	public static double noncentralChiSquareCdf(double x, double nu, double lambda) {
		
		if(nu <= 0) {throw new IllegalArgumentException("Degrees of freedom must be > 0.");}
		if(lambda <= 0) {throw new IllegalArgumentException("Non-centrality parameter must be > 0.");}
		
		return ChiSquareNoncentralDist.cdf(nu, lambda, x);		
		
	}

	/**
	 * Compute the pdf of a Skellam distribution, the distribution of the difference of two Poisson random variables
	 * @param k Value of random variable
	 * @param lambda1 Poisson parameter for distribution 1
	 * @param lambda2 Poisson parameter for distribution 2
	 * @return Probability that the difference (variable1 minus variable2) is equal to the specified value
	 */
	public static double skellamDensity(int k, double lambda1, double lambda2) {
		
		if(lambda1 < 0 || lambda2 < 0) {throw new IllegalArgumentException("Both Poisson parameters must be >= 0.");}
		
		if(k >= 0) return 2 * noncentralChiSquareDensity( 2*lambda1, 2*(k+1), 2*lambda2 );
		
		return 2 * noncentralChiSquareDensity( 2*lambda2, 2*(1-k), 2*lambda1 );
		
	}
	
	/**
	 * Compute the cdf of a Skellam distribution, the distribution of the difference of two Poisson random variables
	 * @param k Value of random variable
	 * @param lambda1 Poisson parameter for distribution 1
	 * @param lambda2 Poisson parameter for distribution 2
	 * @return Probability that the difference (variable1 minus variable2) is less than or equal to the specified value
	 */
	public static double skellamCdf(int k, double lambda1, double lambda2) {
		
		if(lambda1 <= 0 || lambda2 <= 0) {throw new IllegalArgumentException("Both Poisson parameters must be > 0.");}
		
		if(k >= 0) return 1 - noncentralChiSquareCdf( 2*lambda1, 2*(k+1), 2*lambda2 );
		
		return noncentralChiSquareCdf( 2*lambda2, -2*k, 2*lambda1 );
		
	}
	
	/**
	 * Compute the one-tailed P-value of a Skellam distribution, the distribution of the difference of two Poisson random variables
	 * @param k Value of random variable
	 * @param lambda1 Poisson parameter for distribution 1
	 * @param lambda2 Poisson parameter for distribution 2
	 * @return The probability of observing at least as extreme a value under the null hypothesis of a Skellam-distributed random variable (whichever tail is smaller)
	 */
	public static double skellamPvalue(int k, double lambda1, double lambda2) {
		
		if(lambda1 <= 0 || lambda2 <= 0) {throw new IllegalArgumentException("Both Poisson parameters must be > 0.");}
		return Math.min(skellamLeftTail(k, lambda1, lambda2), skellamRightTail(k, lambda1, lambda2));
		
	}
	
	/**
	 * Compute the left tail of a Skellam distribution, the distribution of the difference of two Poisson random variables
	 * @param k Value of random variable
	 * @param lambda1 Poisson parameter for distribution 1
	 * @param lambda2 Poisson parameter for distribution 2
	 * @return The probability of observing at most this value under the null hypothesis of a Skellam-distributed random variable
	 */
	public static double skellamLeftTail(int k, double lambda1, double lambda2) {
		
		if(lambda1 <= 0 || lambda2 <= 0) {throw new IllegalArgumentException("Both Poisson parameters must be > 0.");}
		return skellamCdf(k, lambda1, lambda2);
		
	}
	
	
	/**
	 * Compute the right tail of a Skellam distribution, the distribution of the difference of two Poisson random variables
	 * @param k Value of random variable
	 * @param lambda1 Poisson parameter for distribution 1
	 * @param lambda2 Poisson parameter for distribution 2
	 * @return The probability of observing at least this value under the null hypothesis of a Skellam-distributed random variable
	 */
	public static double skellamRightTail(int k, double lambda1, double lambda2) {
		
		if(lambda1 <= 0 || lambda2 <= 0) {throw new IllegalArgumentException("Both Poisson parameters must be > 0.");}
		return 1 - skellamCdf(k-1, lambda1, lambda2);
		
	}

	
		
}
