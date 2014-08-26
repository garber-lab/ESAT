package broad.core.math;

import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Random;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;

public class MathUtil {
	private static final double LOG10_2 = Math.log10(2);
	/**
	 * Returns the roots of a polynomial with given coefficients. 
	 * 
	 * @param coefs - Polynomial coefficients in increasing order 
	 * (e.g. the polynomial is of the form coefs[0] + coefs[1]*X + coefs[2]*X^2 + ...
	 * @return
	 */
	public static double [][] roots (double [] coefs) {
		Matrix m = new Matrix(coefs.length -1, coefs.length - 1);
		
		m.set(0,0,-coefs[coefs.length - 2]);
		for(int i = 1; i < coefs.length - 1; i++) {
			m.set(i,i-1,1); 
			m.set(0,i,-coefs[coefs.length - i - 2]);
		}
		
		m.timesEquals(1d/coefs[coefs.length - 1]);
		//m.print(3, 5);
		EigenvalueDecomposition ed = m.eig();
		
		double [] real = ed.getRealEigenvalues();
		double [] imaginary = ed.getImagEigenvalues();
		
		double [][] roots = new double[real.length][2];
		for(int i = 0; i < real.length; i++) {
			roots[i][0] = real[i];
			roots[i][1] = imaginary[i];
		}
		
		return roots;
	}
	
	public static double log2(double val) { return Math.log10(val)/LOG10_2;}
	
	public static double max(Collection<? extends Number> list) {
		double max = Double.NEGATIVE_INFINITY;
		if(list != null && !list.isEmpty()) {
			Iterator<? extends Number> it = list.iterator();
			while(it.hasNext()) {
				Number n = it.next();
				double dn = n.doubleValue();
				if(dn > max) {
					max = dn;
				}
			}
		}
		return max;
	}
	
	public static double [] shuffle(double [] original) {
		double [] shuffle = original;
		int N = shuffle.length;
		Random r = new Random();
		
		for (int i = 0; i < N; i++) {
		    int rIdx = i + r.nextInt(N-i);     // between i and N-1
		    double temp = shuffle[i];
		    shuffle[i] = shuffle[rIdx];
		    shuffle[rIdx] = temp;
	    }
		
		return shuffle;
	}
	
	public static<T> List<T> shuffle(List<T> original) {
		List<T> shuffle = original;
		int N = shuffle.size();
		Random r = new Random();
		
		for (int i = 0; i < N; i++) {
		    int rIdx = i + r.nextInt(N-i);     // between i and N-1
		    T temp = shuffle.get(i);
		    shuffle.set(i, shuffle.get(rIdx));
		    shuffle.set(rIdx,temp);
	    }
		
		return shuffle;
	}
	
	   /** sqrt(a^2 + b^2) without under/overflow. **/

	   public static double hypot(double a, double b) {
		   return Jama.util.Maths.hypot(a, b);
	   }
	   
	   /**
	    * The binomial coefficient n choose k
	    * @param n n
	    * @param k k
	    * @return n choose k
	    */
	   public static long binomialCoefficient(int n, int k) {
		   return factorial(n) / (factorial(k) * factorial(n-k));
	   }
	   
	   /**
	    * Factorial function
	    * @param n n
	    * @return n!
	    */
	   public static long factorial(int n) {
	        if      (n <  0) throw new RuntimeException("Underflow error in factorial");
	        else if (n > 20) throw new RuntimeException("Overflow error in factorial");
	        else if (n == 0) return 1;
	        else             return n * factorial(n-1);
	   }

	public static boolean closeTo1(double d) {
		return 0.99999999 < d && d < 1.00000001;
	}

	   
	   
	   
}
