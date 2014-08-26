package broad.core.hmm;

import java.util.ArrayList;
import java.util.Date;
import java.util.Iterator;
import java.util.List;
import java.util.Random;
import java.util.Stack;

import Jama.Matrix;

public class MarkovModel<T> {
	
	private Matrix stateTransitionMatrix;
	private List<MarkovState<T>> states;
	private double [] initialStateTransitionProbabilities;
	private double [] endStateTransitionProbabilities;
	private static final double roundingError = 0.00001;
	private static final double LOG_ZERO_PROB = -1000000;
	
	
	/**
	 * Initialize the Model by specifying the number of states.
	 * The probabilities of starting with a given state will be initialized as uniform: All states are equally likely at start
	 * @param numberOfStates
	 */
	public MarkovModel (int numberOfStates) {
		stateTransitionMatrix = new Matrix(numberOfStates, numberOfStates);
		states = new ArrayList<MarkovState<T>>(numberOfStates);
		initialStateTransitionProbabilities = new double[numberOfStates];
		for(int i = 0; i < numberOfStates; i++) {
			initialStateTransitionProbabilities[i] = 1/(double)numberOfStates;
		}
		
		endStateTransitionProbabilities = new double[numberOfStates];
		for(int i = 0; i < numberOfStates; i++) {
			endStateTransitionProbabilities[i] = 1/(double)numberOfStates;
		}
	}
	
	public void addState(MarkovState<T> state) {
		states.add(state);
	}
	
	public void addState(MarkovState<T> state, int stateIdx) {
		states.set(stateIdx,state);
	}
	
	public void setInitialStateTransitionProbability(int stateIdx, double probability) {
		initialStateTransitionProbabilities[stateIdx] = probability;
	}
	
	public void setEndStateTransitionProbability(int stateIdx, double probability) {
		endStateTransitionProbabilities[stateIdx] = probability;
	}
	
	public MarkovState<T> getState(String stateName) {
		MarkovState<T> theState = null;
		Iterator<MarkovState<T>> it = states.iterator();
		
		while(it.hasNext() && theState == null) {
			MarkovState<T> state = it.next();
			if(state.getName().equals(stateName)) {
				theState = state;
			}
		}
		
		return theState;
	}
	
	public List<Integer> emitPath(int numOfSteps) {
		Stack<Integer> stateStack = new Stack<Integer>();
		stateStack.push(drawInitialState());
		
		for(int i = 0; i < numOfSteps; i++) {
			int prior = stateStack.peek();
			stateStack.push(drawNextState(prior));
		}
		return stateStack;
	}


	public void setStateTransitionProbability(int stateIdx1, int stateIdx2, double probability) {
		stateTransitionMatrix.set(stateIdx1, stateIdx2, probability);
	}
	
	public double computePathLogLikelihood(List<T> observationSequence, short [] path, int initialStateIdx, int endStateIdx) {
		if(observationSequence.size() != path.length) {
			throw new IllegalArgumentException("Got " + observationSequence.size() + " but a path of " + path.length + " states, both should be of the same length");
		}
		Matrix transitionLogP   = logTransformTransitionMatrix();
		
		double logLikelihood = transitionLogP.get(initialStateIdx, path[0]) + states.get(path[0]).getEmissionLogProbability(observationSequence.get(0));
		
		for(int i = 1; i < path.length; i++) {
			logLikelihood = logLikelihood + states.get(path[i]).getEmissionLogProbability(observationSequence.get(i)) + transitionLogP.get(path[i-1], path[i]);
		}
		
		return logLikelihood + transitionLogP.get(path[path.length - 1], endStateIdx);
	}
	
	public short[] viterbiMostLikelyEstimation(List<T> observationSequence) throws BadModelException {
		short[] stateSequence = new short[observationSequence.size()];
		
		if(observationSequence.size() == 0) {
			return stateSequence;
		}
		
		double [] initStateLogP = logTransformInitialStateTransitionProbabilities();
		double [] endStateLogP = logTransformEndStateTransitionProbabilities();
		Matrix transitionLogP   = logTransformTransitionMatrix();
		
		double [] llsofar = new double [states.size()]; //log likelihoods for each state computed so far as they are computed.
		double [] llprevious = new double [states.size()]; //log likelihoods as they were computed in prior iteration, sort of a cache.
		short [][] pointerMatrix = new short[states.size()][observationSequence.size()];
		
		if(getConsistencyProblems().size() >0){
			throw new BadModelException(this);
		}
		
		// Initialize the first observation: Use initialStateProbabilities for this.
		// Also initialize the first observation pointers to -1
		for(int i = 0; i < states.size(); i++) {
			MarkovState<T> state = states.get(i);
			llsofar[i] = state.getEmissionLogProbability(observationSequence.get(0)) + initStateLogP[i];
			//System.err.println("init, state " + i+" observation " + observationSequence.get(0) + " state prob: " + llsofar[i]);
			pointerMatrix[i][0] = -1;
		}
		
		// Dynamic programming step.
		for(int i = 1; i < observationSequence.size(); i++) {
			//System.err.println("position " + i + " obs " + observationSequence.get(i));
			for(short stateIdx = 0; stateIdx < states.size(); stateIdx++) {
				llprevious[stateIdx] = llsofar[stateIdx];
			}
			for(short stateIdx = 0; stateIdx < states.size(); stateIdx++) {
				short backPtr = -1;
				double maxStateLL = Double.NEGATIVE_INFINITY;
				for(short priorStateIdx = 0; priorStateIdx < states.size(); priorStateIdx++) {
					double logLikelihoodOfTransitioningHereFromThere =  llprevious[priorStateIdx] + transitionLogP.get(priorStateIdx, stateIdx);
					if(logLikelihoodOfTransitioningHereFromThere == Double.NEGATIVE_INFINITY) { logLikelihoodOfTransitioningHereFromThere = LOG_ZERO_PROB;} // AVOID Negative infinity.
					//System.err.print("\t\ttransfromHereToThere " + logLikelihoodOfTransitioningHereFromThere + " maxStateLL " + maxStateLL);
					if(maxStateLL < logLikelihoodOfTransitioningHereFromThere) {
						maxStateLL = logLikelihoodOfTransitioningHereFromThere;
						backPtr = priorStateIdx;
						//System.err.print(" Was higher than maxStateLL!");
					} 
					//System.err.println(" finished with state " + priorStateIdx + " backPtr: " +  backPtr);
				}
				double obsStateEmissionProb = states.get(stateIdx).getEmissionLogProbability(observationSequence.get(i));
				if(obsStateEmissionProb == Double.NEGATIVE_INFINITY) { obsStateEmissionProb = LOG_ZERO_PROB;} // AVOID Negative infinity.
				llsofar[stateIdx] =  maxStateLL + obsStateEmissionProb;
				//System.err.println("\tstate " + stateIdx + " emission prob: " + states.get(stateIdx).getEmissionLogProbability(observationSequence.get(i)) + " corrected prob: " + obsStateEmissionProb);
				pointerMatrix[stateIdx][i] = backPtr;
			}
		}
		
		//Retrace step, likelihood is the max of llsofar
		
		//first find last pointer and total likelihood
		short lastState = -1;
		double logLikelihood = Double.NEGATIVE_INFINITY;
		for(short stateIdx = 0; stateIdx < states.size(); stateIdx++) {
			//System.err.print("llsofar["+stateIdx+"] "+llsofar[stateIdx] + "---- endStateLogP["+stateIdx+"] " + endStateLogP[stateIdx] + " sum " + (llsofar[stateIdx] + endStateLogP[stateIdx]) + "  loglikel " + logLikelihood);
			if(llsofar[stateIdx] + endStateLogP[stateIdx] > logLikelihood) {
				lastState = stateIdx;
				logLikelihood = llsofar[stateIdx] + endStateLogP[stateIdx];
			}
			//System.err.println(" last state for  was " + stateIdx + " likelihood " + logLikelihood);
		}
		
		//System.out.println("loglikelihood: " + logLikelihood);
		// traceback
		stateSequence[observationSequence.size() - 1]= lastState;
		for(int i = observationSequence.size() - 1; i > 0; i--) {
			//System.err.print("i " + i  +" stateSequence ");
			//System.err.println(  stateSequence[i]);
			stateSequence[i - 1]= pointerMatrix[stateSequence[i]][i];
		}
		
		return stateSequence;
	}
	
	public ForwardResult<T> runForwardAlgorithm(List<T> observationSequence) throws BadModelException {
		ForwardResult<T> result = new ForwardResult<T>(this, observationSequence);
		
		if(getConsistencyProblems().size() >0){
			throw new BadModelException(this);
		}
		
		//Initialization step.
		for(int i = 0; i < states.size(); i++) {
			MarkovState<T> state = states.get(i);
			result.set(i, 0, state.getEmissionProbability(observationSequence.get(0)) * initialStateTransitionProbabilities[i]);
		}
		result.normalizeColumn(0);
		// Dynamic programming step.
		for(int i = 1; i < observationSequence.size(); i++) {
			for(short stateIdx = 0; stateIdx < states.size(); stateIdx++) {
				double probOfSiteAndState = 0;
				for(short priorStateIdx = 0; priorStateIdx < states.size(); priorStateIdx++) {
					probOfSiteAndState =  probOfSiteAndState + result.get(priorStateIdx, i - 1)* stateTransitionMatrix.get(priorStateIdx, stateIdx);
				}
				
				result.set(stateIdx, i, probOfSiteAndState * states.get(stateIdx).getEmissionProbability(observationSequence.get(i)));
			}
			
			result.normalizeColumn(i);
			
		}
		
		
		return result;
	}
	
	
	public BackwardResult<T> runBackwardAlgorithm(List<T> observationSequence, ForwardResult<T> forwardData) throws BadModelException {
		//long start = (new Date()).getTime();
		BackwardResult<T> result = new BackwardResult<T>( forwardData);
		if(getConsistencyProblems().size() >0){
			throw new BadModelException(this);
		}
		//System.err.println("\tcreated backwardResult: " + ((new Date()).getTime() - start) );
		
		//Initialization
		for(int i = 0; i < states.size(); i++) {
			result.set(i, observationSequence.size() - 1, 1/*endStateTransitionProbabilities[i]*/); //TODO: FIX this
		}
		//System.err.println("\tInitialization: " + ((new Date()).getTime() - start) );
		// Recursion
		for(int i = observationSequence.size() - 2; i >=0; i--) {
			// hack to avoid computing emmision  probabilities numOfStates^2 times.
			double [] nextStateEmissionProbs = new double[states.size()];
			for(int stateIdx = 0; stateIdx < states.size(); stateIdx++) {
				nextStateEmissionProbs[stateIdx] = states.get(stateIdx).getEmissionProbability(observationSequence.get(i+1));
			}
			for(int stateIdx = 0; stateIdx < states.size(); stateIdx++) {
				double jointProbOfSeqAndStateAtPos = 0;
				for(int nextStateIdx = 0; nextStateIdx < states.size(); nextStateIdx++) {
					jointProbOfSeqAndStateAtPos = jointProbOfSeqAndStateAtPos + 
							result.get(nextStateIdx, i + 1)*stateTransitionMatrix.get(stateIdx, nextStateIdx)*nextStateEmissionProbs[nextStateIdx];
				}
				
				result.set(stateIdx, i, jointProbOfSeqAndStateAtPos / forwardData.getScaling(i+1));
				//System.err.println("\tNormilzed : " + ((new Date()).getTime() - start) );
			}
		}
		//System.err.println("\trecursion: " + ((new Date()).getTime() - start) );
		return result;
	}
	
	public ForwardResult<T> runForwardAlgorithmNoRescaling(List<T> observationSequence) {
		ForwardResult<T> result = new ForwardResult<T>(this, observationSequence);
		
		double [] initStateLogP = logTransformInitialStateTransitionProbabilities();
		
		//Initialization step.
		for(int i = 0; i < states.size(); i++) {
			MarkovState<T> state = states.get(i);
			result.set(i, 0, state.getEmissionLogProbability(observationSequence.get(0)) + initStateLogP[i]);
		}
		
		// Dynamic programming step.
		for(int i = 1; i < observationSequence.size(); i++) {
			for(short stateIdx = 0; stateIdx < states.size(); stateIdx++) {
				double probOfSiteAndState = 0;
				for(short priorStateIdx = 0; priorStateIdx < states.size(); priorStateIdx++) {
					probOfSiteAndState =  probOfSiteAndState + Math.exp(result.get(priorStateIdx, i - 1))* stateTransitionMatrix.get(priorStateIdx, stateIdx);
				}
				
				result.set(stateIdx, i, Math.log(probOfSiteAndState) + states.get(stateIdx).getEmissionLogProbability(observationSequence.get(i)));
			}
		}
		
		
		return result;
	}
	private BackwardResult<T> runBackwardAlgorithmNoScaling(List<T> observationSequence, ForwardResult<T> forwardData) {
		BackwardResult<T> result = new BackwardResult<T>( forwardData);
		double [] initStateLogP = logTransformInitialStateTransitionProbabilities(); // will use initial state probs to initialize the algorithm.
		
		//Initialization
		for(int i = 0; i < states.size(); i++) {
			result.set(i, observationSequence.size() - 1, initStateLogP[i]);
		}
		
		// Recursion
		for(int i = observationSequence.size() - 2; i >=0; i--) {

			for(int stateIdx = 0; stateIdx < states.size(); stateIdx++) {
				double jointProbOfSeqAndStateAtPos = 0;
				for(int nextStateIdx = 0; nextStateIdx < states.size(); nextStateIdx++) {
					jointProbOfSeqAndStateAtPos = jointProbOfSeqAndStateAtPos + 
							Math.exp(result.get(nextStateIdx, i + 1))*stateTransitionMatrix.get(stateIdx, nextStateIdx)*states.get(nextStateIdx).getEmissionProbability(observationSequence.get(i+1));
				}
				result.set(stateIdx, i, Math.log(jointProbOfSeqAndStateAtPos));
			}
		}
		
		return result;
	}

	protected Matrix getStateTransitionMatrix() {
		return stateTransitionMatrix;
	}

	protected List<MarkovState<T>> getStates() {
		return states;
	}

	protected double[] getInitialStateTransitionProbabilities() {
		return initialStateTransitionProbabilities;
	}
	
	protected List<String> getConsistencyProblems() {
		ArrayList<String> problemList = new ArrayList<String>();
		
		if(states == null || states.isEmpty()) {
			problemList.add("No states were specified");
		}
		
		for(int i = 0; i < stateTransitionMatrix.getRowDimension(); i++) {
			double  total = 0;
			for(int j = 0; j < stateTransitionMatrix.getColumnDimension(); j++) {
				total = total + stateTransitionMatrix.get(i, j);
			}
			if(Math.abs( 1 - total) > roundingError) {
				problemList.add("Row " + i + " of state transition matrix does not add up to " + total + " it must add up to 1");
			}
		}
		
		double total = 0;
		for(int i = 0; i < initialStateTransitionProbabilities.length; i++) {
			total = total + initialStateTransitionProbabilities[i];
		}
		if(Math.abs(total - 1) > roundingError) {
			problemList.add("Initial state probabilities add up to " + total + " they must add up to 1");
		}
		
		return problemList;
	}
	
	protected int drawInitialState() {
		double draw  =  new Random().nextDouble();
		int state = 0;
		double cummulativeProb = initialStateTransitionProbabilities[0];
		while(state < initialStateTransitionProbabilities.length && draw > cummulativeProb ) {
			cummulativeProb += initialStateTransitionProbabilities[++state];
		}
		return state;
	}
	
	protected int drawNextState(int prior) {
		double draw  =  new Random().nextDouble();
		int nextStateIdx = 0;
		double cummulativeProb = stateTransitionMatrix.get(prior, 0);
		while(nextStateIdx < stateTransitionMatrix.getColumnDimension() && draw > cummulativeProb ) {
			cummulativeProb += stateTransitionMatrix.get(prior, ++nextStateIdx);
		}
		return nextStateIdx;
	}
	
	private double [] logTransformInitialStateTransitionProbabilities() {
		double [] t = new double[initialStateTransitionProbabilities.length];
		for(int i = 0; i < initialStateTransitionProbabilities.length; i++) {
			t[i] = Math.log(initialStateTransitionProbabilities[i]);
		}
		
		return t;
	}
	
	private double [] logTransformEndStateTransitionProbabilities() {
		double [] t = new double[endStateTransitionProbabilities.length];
		for(int i = 0; i < endStateTransitionProbabilities.length; i++) {
			t[i] = Math.log(endStateTransitionProbabilities[i]);
		}
		
		return t;
	}
	
	private Matrix logTransformTransitionMatrix() {
		Matrix t = new Matrix(stateTransitionMatrix.getRowDimension(), stateTransitionMatrix.getColumnDimension());
		for(int i = 0; i < stateTransitionMatrix.getRowDimension(); i++) {
			for(int j = 0; j < stateTransitionMatrix.getColumnDimension(); j++) {
				t.set(i,j,Math.log(stateTransitionMatrix.get(i,j)));
			}
		}
		
		return t;
	}
	
	public static class ForwardResult<T> {

		private static final long serialVersionUID = -1295929611441551168L;
		MarkovModel<T> sourceModel;
		List<T> observations;
		double[] scalings;
		double logProbability;
		double[][] data;

		protected ForwardResult(MarkovModel<T> sourceModel, List<T> observationSequence) {
			//super(sourceModel.getStates().size(), observationSequence.size());
			data = new double[sourceModel.getStates().size()][observationSequence.size()];
			this.sourceModel = sourceModel;
			this.observations = observationSequence;
			this.scalings = new double[observationSequence.size()];
		}
		
		public double getLogProbability() {
			return logProbability;	
		}
		
		public void setScaling(int idx, double val) {
			scalings[idx] = val;
			logProbability = logProbability + Math.log(val);
		}
		
		public double getScaling(int idx) { return scalings[idx]; }
		
		public double getProbability() {
			return Math.exp(logProbability);
		}
		
		public BackwardResult<T> runBackwardAlgorithm() throws BadModelException {
			return sourceModel.runBackwardAlgorithm(observations, this);
		}
		
		void normalizeColumn(int observationIdx) {
			double normConstant = 0;
			for(int i = 0; i < sourceModel.getStates().size(); i++) {
				normConstant = normConstant + data[i][observationIdx];
			}
			setScaling(observationIdx,normConstant);
			for(int i = 0; i < sourceModel.getStates().size(); i++) {
				data[i][observationIdx]= data[i][observationIdx]/normConstant;
			}
		}
		
		public int getRowDimension() {return data.length;}
		public int getColumnDimension() {return data.length == 0 ? 0 : data[0].length;}
		public double get(int row, int col) { return data[row][col];}
		public void set(int row, int col, double val) { data[row][col] = val;}
	}
	
	public static class BackwardResult<T>  {

		private static final long serialVersionUID = -2641581777198536941L;
		private ForwardResult<T> forwardData;
		private double observationLogProb;
		double [][] data;
		protected BackwardResult(ForwardResult<T> forwardData) {
			//super(forwardData.getRowDimension(), forwardData.getColumnDimension());
			data = new double[forwardData.getRowDimension()][forwardData.getColumnDimension()];
			this.forwardData = forwardData;
			observationLogProb = forwardData.getLogProbability();
		}
		
		/**
		 * Computes the posterior probability of a given observation being emitted by a state
		 * @param stateIdx The state index
		 * @param atObservation the observation number
		 * @return The posterior log probability
		 */
		public double getPosteriorLogProbability(int stateIdx, int atObservation) {
			return Math.log(getPosteriorProbability(stateIdx, atObservation));
		}
		
		public double getPosteriorProbability(int stateIdx, int atObservation) {
			return forwardData.get(stateIdx, atObservation) * get(stateIdx, atObservation) ;
		}
		
		public int getRowDimension() {return data.length;}
		public int getColumnDimension() {return data.length == 0 ? 0 : data[0].length;}
		public double get(int row, int col) { return data[row][col];}
		public void set(int row, int col, double val) { data[row][col] = val;}
	}

}
