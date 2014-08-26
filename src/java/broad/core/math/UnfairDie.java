package broad.core.math;

import java.util.Random;

import org.apache.log4j.Logger;

/**
 * 
 * @author prussell
 *
 */
public class UnfairDie {
	
	@SuppressWarnings("unused")
	private static Logger logger = Logger.getLogger(UnfairDie.class.getName());
	
	// Probabilities per position
	private double[] probs;
	private int numPositions;
	// Helper list of boundaries of intervals within a uniform distribution
	private double[] uniformBounds;
	private Random rand;
	
	/**
	 * @param probabilities Vector of probabilities where the ith probability corresponds to the number i
	 */
	public UnfairDie(double[] probabilities) {
		
		checkAddsToOne(probabilities);
		if(probabilities.length < 1) {
			throw new IllegalArgumentException("Array of probabilities is empty");
		}
		
		probs = probabilities;
		numPositions = probs.length;
		rand = new Random(System.currentTimeMillis());
		
		// Construct the list of bounds for uniform distribution
		uniformBounds = new double[numPositions - 1];
		uniformBounds[0] = probs[0];
		for(int i = 1; i < numPositions - 1; i++) {
			uniformBounds[i] = uniformBounds[i-1] + probs[i];
		}
		
	}
	
	/**
	 * Roll the unfair die and get a random integer
	 * @return Random integer generated from unfair die
	 */
	public int roll() {
		double d = rand.nextDouble();
		if(d < uniformBounds[0]) {
			return 0;
		}
		if(d > uniformBounds[uniformBounds.length - 1]) {
			return uniformBounds.length;
		}
		for(int i = 0; i < uniformBounds.length - 1; i++) {
			if(d >= uniformBounds[i] && d < uniformBounds[i + 1]) {
				return i + 1;
			}
		}
		throw new IllegalStateException("Didn't find interval for random number " + d + ".");
	}
	
	private void checkAddsToOne(double[] list) {
		double total = 0;
		for(int i = 0; i < list.length; i++) {
			total += list[i];
		}
		if(!MathUtil.closeTo1(total)) {
			throw new IllegalArgumentException("Probabilities add up to " + total);
		}
	}
	
}
