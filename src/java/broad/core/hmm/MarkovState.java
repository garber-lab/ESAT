package broad.core.hmm;

public interface MarkovState<T> {
	
	/**
	 * 
	 * @param observation The roll of a dice, an alignment column, a nucleotide, etc.
	 * @return P(observation | this state)
	 */
	double getEmissionProbability(T observation);
	
	/**
	 * @return The name of this Markov State
	 */
	String getName();
	
	/**
	 * In many cases the implementing class may be able to implement this method so that the logarithm
	 * is not on every call.
	 * @param observation The roll of a dice, an alignment column, a nucleotide, etc.
	 * @return natural log of P(observation | this state)
	 */
	double getEmissionLogProbability(T observation);
	
	/**
	 * Uses its background distribution to emit a sampled observation.
	 * @return
	 */
	T emitObservation();
	
	

}
