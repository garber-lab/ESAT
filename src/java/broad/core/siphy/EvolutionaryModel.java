package broad.core.siphy;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Stack;

import nextgen.core.annotation.Annotation;

import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.apache.commons.math3.stat.Frequency;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import broad.core.annotation.LightweightGenomicAnnotation;
import broad.core.math.MathUtil;
import broad.core.sequence.Sequence;

/**
 * Represents an Markovian Evolutionary Model
 * @author mgarber
 *
 */
public class EvolutionaryModel {
	private EvolutionaryModelParameters parameters;
	private HashMap<Double, Matrix> transitionMatrixCache = new HashMap<Double, Matrix>();
	static final int MAX_ITERATIONS = 10;
	static final double SMALL_DIFF  = 0.0001;
	static final double TINY_DIFF   = 0.00001;
	private static DecimalFormat numberFormat = new  DecimalFormat("##0.####");
	int alphabetSize;
    ///test for version 
	Matrix R;  //Rate matrix after extracting frequency;
	Matrix Q; //Full Rate matrix copy (R*diagRho).
	
	Matrix D;  //Eigenvalue diagonal matrix 
	Matrix V;  //Eigenvector matrix
	Matrix iV; //inverse eigenvector matrix
	Matrix pi; // Fitted stationary distribution
	private double omega; //Fitted proportionality constraint.
	
	LinkedHashMap<Double, Frequency> binStatistics;
	
	Map<Integer, NodeLikelihoodParameters> nodeFittingParamMap;
	
	public EvolutionaryModel(EvolutionaryModelParameters parameters) {
		super();
		setModelParameters(parameters);
		resetPI();
		omega = 1;
		nodeFittingParamMap = new HashMap<Integer, NodeLikelihoodParameters>();
	}
	
	public void setOmega(double omega) {
		this.omega = omega;
		clearCaches();
	}
	
	public void RPIDecomposition() {
		R = extractEquilibriumFromRateMatrix();
	}
	
	/**
	 * Computes current model's rate matrix. If the model is Q = R*rho 
	 * then it computes Q, otherwise it returns the rate matrix in the 
	 * Model Parameters
	 * @return An update rate matrix.
	 */
	protected Matrix getRateMatrix() {
		Matrix Q = parameters.getRateMatrix();
		if(R != null && pi != null) {
			//Q.print(4, 4);
			//pi.print(4, 4);
			//R.print(4,4);
			
			//Jan 26/09 changed to test Or Zuk's suggestion
			Q = R.times(pi); //Original decomposition
			//addPseudoCounts(pi);
			//Q = pi.inverse().times(R); //Or's suggestion
			
			//Q.print(4, 4);
			//System.out.println("Printed getRateMatrix Q,pi,R,Q");
			
			//pi.print(15, 10);
			
			for(int i = 0; i < Q.getRowDimension(); i++) {
				double t = 0;
				for(int j = 0; j < Q.getColumnDimension(); j++) {
					if(j!=i) {
						t += Q.get(i,j);
					}
					Q.set(i, i, -t);
				}
				//System.out.println("row i, diag + off diag: " + t);			
			}
			
		}

		return Q;
	}
	
	public Map<String, Short> sample() {
		double [] bgDistribution = parameters.getBackgroundNucleotideFreqs();
		short rootState = sampleFromAlphabet(bgDistribution);
		
		Phylogeny tree = parameters.getTree();
		Map<String, Short> column = new LinkedHashMap<String, Short>(parameters.getTree().getNumberOfExternalNodes());
		column.putAll(sample(tree.getRoot().getChildNode1(), rootState));
		column.putAll(sample(tree.getRoot().getChildNode2(), rootState));
		return column;
	}


	
	public Map<String,Short> sample(PhylogenyNode node, short parentState) {
		Map<String, Short> column = new HashMap<String, Short>(parameters.getTree().getNumberOfExternalNodes());
		
		Matrix edgeTransitionMatrix = computeTransitions(node.getDistanceToParent());
		double [] transitionProbs = new double [edgeTransitionMatrix.getRowDimension()];
		for(int i = 0; i < transitionProbs.length; i++) {
			transitionProbs[i] = edgeTransitionMatrix.get(parentState, i);
		}
		
		short nodeState = sampleFromAlphabet(transitionProbs);
		if(node.isExternal()) {
			column.put(node.getSeqName(), nodeState);
		} else {
			column.putAll(sample(node.getChildNode1(), nodeState));
			column.putAll(sample(node.getChildNode2(), nodeState));
		}
		
		return column;
	}
	
	/**
	 * Implementation of the EM algorithm to fit the base distribution rho at a site
	 * assuming omega is constant.
	 * @throws UnableToFitException 
	 */
	
	public PiFit fitPI(Map<String, Matrix> alignmentColumn, Phylogeny tree)  {
		PiFit fit = null;
		PhylogenyNode root = tree.getRoot();
		
		double newOldDist = 10;
		short iteration = 0;
		clearCaches();
		resetPI();
		R = extractEquilibriumFromRateMatrix();
		prepareRateMatrix();
		double initialLogLikelihood = Math.log(pruneAndPeel(alignmentColumn, root, 0));
		
		//System.out.println("\tNeutral likelihood: " + initialLogLikelihood);
		
		do {
			//System.out.println("\tStart of iteration " + (iteration));
			Matrix preIterationPI = pi.copy();
			
			try {
				fit = piEMIteration(alignmentColumn, root, preIterationPI);
			} catch (UnableToFitException e) {
				e.printStackTrace();
				break;
			}
			
			newOldDist = 0;
			for(int i = 0; i < pi.getRowDimension(); i++) {
				newOldDist += Math.abs(preIterationPI.get(i, i) - fit.fittedPi.get(i, i));
			}
			
			/*
			System.out.println("\n\tEnding Iteration, Old PI ");
			preIterationPI.print(4, 6);
			System.out.println("\tFitted likelihood " + fit.fittedLikelihood + " PI " );
			fit.fittedPi.print(4, 6);
			*/
			this.pi = fit.fittedPi;
			//R = extractRhoFromRateMatrix(diagRho);
			//prepareRateMatrix();
			//if(iteration == 0) {
			//	initialLogLikelihood = fit.unfittedLikelihood;
			//}
		} while(newOldDist > SMALL_DIFF && iteration++ < MAX_ITERATIONS);
		if(fit != null) {
			fit.numOfIterations = iteration;
			fit.unfittedLikelihood = initialLogLikelihood;
		}
		nodeFittingParamMap = new HashMap<Integer, NodeLikelihoodParameters>(root.getNumberOfChildNodes());
		clearCaches();
		return fit;
	}

	/**
	 * Implementation of Felsenstein's prunning algorithm to compute the 
	 * likelihood given this model.
	 * @param observedData
	 * @return
	 * @throws MathException 
	 */
	public OmegaFit fitOmega(Map<String, Matrix> leafValues, Phylogeny tree, int window) {
		PhylogenyNode root = tree.getRoot();
		double originalOmega = omega;
		ChiSquaredDistribution chiSq = new ChiSquaredDistribution(1);
		/*
		System.out.println("Q");
		Q.print(6, 4);
		RPIDecomposition();
		System.out.println("R");
		R.print(6,4);
		System.out.println("Q again");
		Q.print(6, 4);		
		*/
		//System.out.println("Using tree " + tree.toNewHampshire(true));
		double newOmega = omega;
		int iteration = 0;
		OmegaFit fit = new OmegaFit();
		double [] data = null;
		do {
			//System.out.println("Iteration " + iteration);
			data = omegaEMIteration(leafValues, root, newOmega, window);
			newOmega = data[0];
			if(iteration == 0) {
				fit.initialLogLikelihood = data[3];
			}
			//System.out.println("\t results: omega " + data[0] + " transitions " + data[1] + " totalTime " + data[2] + " log likelihood " + data[3]);
		}while(Math.abs(omega - newOmega) > SMALL_DIFF && iteration++ < MAX_ITERATIONS);

		fit.fittedLogLikelihood = data[3];
		fit.omega = data[0];
		fit.numOfIterations = iteration;
		fit.transitions = data[1];
		fit.totalTime = data[2];
		fit.pVal = 1- chiSq.cumulativeProbability(fit.getLogOddsScore());
		//System.out.println("\tdone, logodds " + fit.getLogOddsScore() + " pval " + fit.pVal);		
		
		setOmega(originalOmega); // reset omega since fitting parameters should not change omega. This is dangerous and really not thread safe.
		nodeFittingParamMap = new HashMap<Integer, NodeLikelihoodParameters>(root.getNumberOfChildNodes());
		clearCaches();
		//System.out.println("Final omega for site: " + data[0] + " log odds " + (fit.fittedLogLikelihood - fit.initialLogLikelihood));
		return fit;
	}
	
	public PiFit piEMIteration(Map<String, Matrix> column, PhylogenyNode root, Matrix newPI) throws UnableToFitException {		
		PiFit fit = new PiFit();
		nodeFittingParamMap = new HashMap<Integer, NodeLikelihoodParameters>(root.getNumberOfChildNodes());
		
		//fit.fittedPi = newPI;
		clearCaches();
		prepareRateMatrix();
		//System.out.println("\t\tQ:\n");
		//Q.print(6, 4);
		
		fit.fittedPi = newPI;
		
		double startLikelihood = pruneAndPeel(column, root, 0);
		double initialLogLikelihood = Math.log(startLikelihood);
		//NodeLikelihoodParameters rootPrms = nodeFittingParamMap.get(root.getID());
		//Matrix rootProbs = rootPrms.computeProbabilityOfLetterAtNode();
		//System.out.println("\t\tStart Likelihood " + startLikelihood + ", precise logLikelihood " + initialLogLikelihood+ " root probs (" + rootProbs.get(0, 0) + "," + 
		//		rootProbs.get(1, 0) + "," + rootProbs.get(2, 0) + "," + rootProbs.get(3, 0) + ")");
		//fit.fittedLikelihood = initialLikelihood;
		Matrix N = computeSufficientStatistics(root, startLikelihood);
		
		// Interim computations, the following quantities are repeatedly used in computing the coefficients for the Lagrange multipliers.
		// 1 Compute T (total time matrix) from N (# of transition matrix).
		double [] T = new double[N.getRowDimension()];
		for(int i = 0; i < N.getRowDimension(); i++) {
			T[i] = N.get(i, i) /(omega * Q.get(i, i));  //Recall that T(i) = N(i,i)/Q_ii = N(i,i)/(omega * R_ii * rho_i) but our Q is just R*rho so we need to adjust for omega.
		}
		/*
		System.out.println("\t\tDuration " + T[0] +"," + T[1] + "," + T[2] + "," + T[3]);
		System.out.println("\t\tTransitions A to C,G,T " + N.get(0, 1) +"," + N.get(0,2) + "," + N.get(0,3) );
		System.out.println("\t\tTransitions C to A,G,T " + N.get(1, 0) +"," + N.get(1,2) + "," + N.get(1,3) );
		System.out.println("\t\tTransitions G to C,A,T " + N.get(2, 0) +"," + N.get(2,1) + "," + N.get(2,3) );
		*/
		// We need the posterior probability of a given base at the root.
		NodeLikelihoodParameters rootLikelihoods = nodeFittingParamMap.get(root.getID());
		Matrix posteriorRoot = rootLikelihoods.computeProbabilityOfLetterAtNode();
		
		//System.out.print("posterior Root probs: ");
		//posteriorRoot.transpose().print(4,6);
		
		//double previousLogLikelihood = computeLogLikelihood(N, T, posteriorRoot);
		fit.unfittedLikelihood = initialLogLikelihood;//previousLogLikelihood;
		//System.out.println("\t\tStarting with loglikelihood " + previousLogLikelihood);
		
		double [] A = new double[N.getRowDimension()];
		//System.out.print("\t\tA: ");
		for(int i = 0; i < N.getRowDimension(); i++) {
			A[i] = posteriorRoot.get(i, 0);
			for(int j = 0; j < N.getRowDimension(); j++) {
				if(j != i) {
					A[i] += N.get(j, i);
				}
			}
			//System.out.print(A[i] + ", ");
		}
		//System.out.print("\n\t\tB: ");
		double []B = new double[N.getRowDimension()];
		for(int i = 0; i < N.getRowDimension(); i++) {
			for(int j = 0; j < N.getRowDimension(); j++) {
				if(j != i) {
					B[i] += T[j] * R.get(j, i) ;
				}
			}
			B[i] *= omega;
			//System.out.print(B[i] +",");
		}
		//System.out.println("");
		//To avoid multiplying the same quantities over and over
		double [][] Bduos = new double[N.getRowDimension() - 1][N.getRowDimension()];
		for(int i = 0; i < N.getRowDimension() - 1; i++) {
			for (int j = i+1; j < N.getRowDimension(); j++) {
				Bduos[i][j] = B[i]*B[j];
			}
		}
		
		double [][][] Btrios = new double[N.getRowDimension() - 2][N.getRowDimension()- 1][N.getRowDimension()];
		for(int i = 0; i < N.getRowDimension() - 2; i++) {
			for (int j = i+1; j < N.getRowDimension() - 1; j++) {
				for(int k = j+1; k < N.getRowDimension(); k++) {
					Btrios[i][j][k] = Bduos[i][j]*B[k];
				}
			}
		}
		
		double [] mlRhoPolynomialCoeffs = new double[5];
		
		
		mlRhoPolynomialCoeffs[0] = B[0]*B[1]*B[2]*B[3] - 
			( Btrios[1][2][3]*A[0] + Btrios[0][2][3]*A[1] +
			  Btrios[0][1][3]*A[2] + Btrios[0][1][2]*A[3]);
			 
		mlRhoPolynomialCoeffs[1] = Btrios[0][1][2] + Btrios[0][2][3] + Btrios[1][2][3] + Btrios[0][1][3] - 
			(A[0]*(Bduos[1][2] + Bduos[1][3] + Bduos[2][3]) + 
			 A[1]*(Bduos[0][2] + Bduos[0][3] + Bduos[2][3]) +
			 A[2]*(Bduos[0][1] + Bduos[0][3] + Bduos[1][3]) +
			 A[3]*(Bduos[0][1] + Bduos[0][2] + Bduos[1][2])
			);
		
		mlRhoPolynomialCoeffs[2] = Bduos[0][1] + Bduos[0][2] + Bduos[0][3] + Bduos[1][2] + Bduos[1][3] + Bduos[2][3] -
			( A[0]*(B[1]+B[2]+B[3]) + A[1]*(B[0]+B[2]+B[3]) + A[2]*(B[0]+B[1]+B[3]) + A[3]*(B[0]+B[1]+B[2]));
		
		mlRhoPolynomialCoeffs[3] = B[0]+B[1]+B[2]+B[3] - (A[0]+A[1]+A[2]+A[3]);
		
		mlRhoPolynomialCoeffs[4] = 1;
		
		//System.out.println("\t\tp = x^4 + " + mlRhoPolynomialCoeffs[3] +"x^3 " + mlRhoPolynomialCoeffs[2] + "x^2 " + mlRhoPolynomialCoeffs[1] +"x " + mlRhoPolynomialCoeffs[0]);
		double [][] roots = MathUtil.roots(mlRhoPolynomialCoeffs);
		List<Double> realRoots = new ArrayList<Double>();
		// we now check for real roots and ones for which the denominator (B_i + root) >= 0
		for (int i = 0; i < roots.length; i++) {
			boolean good = false;
			//System.out.println("\t\troot " + i + " " + roots[i][0]  + " " + roots[i][1] + "i");
			if(Math.abs(roots[i][1]) < TINY_DIFF) {
				good = true;
				for(int j = 0; j < B.length; j++) {
					//System.out.println("\tB["+j+"] + root["+i+"] = "+(B[j] + roots[i][0]));
					if(B[j] + roots[i][0] <= 0) {
						good = false;
					}
				}
			} 
			if(good) {
				realRoots.add(roots[i][0]);
			}
		}
		
		if(realRoots.size() < 1) {
			System.out.println("Could not fit: ");
			System.out.println("\t\tp = x^4 + " + mlRhoPolynomialCoeffs[3] +"x^3 " + mlRhoPolynomialCoeffs[2] + "x^2 " + mlRhoPolynomialCoeffs[1] +"x " + mlRhoPolynomialCoeffs[0]);
			throw new UnableToFitException("Oh Oh no roots where found for rho estimator lagrange multiplier ");
		}

		Iterator<Double> lagrangeGoodRootsIt = realRoots.iterator();
		double lastRootsLikelihood = initialLogLikelihood;//previousLogLikelihood;
		double rootLogLikelihood = 0;
		int rootNum = 1;
		//System.out.println("\t\tTrying all roots:");
		while(lagrangeGoodRootsIt.hasNext()) {
			double lMultiplier = lagrangeGoodRootsIt.next();
			Matrix resultingPI  = new Matrix(newPI.getRowDimension(), newPI.getColumnDimension());
			double total = 0;
			for(int i = 0; i < resultingPI.getColumnDimension(); i++) {
				double numerator = A[i];
				double denominator = B[i] + lMultiplier;
				double pi_i = numerator/denominator;
				resultingPI.set(i, i, pi_i);
				total += pi_i;
			}
			//System.out.println("\t\tUsed root " + rootNum++ + " - " + lMultiplier + ", Trying pi ");
			resultingPI.timesEquals(1d/total);
			//resultingPI.print(4, 6);
			pi = resultingPI;
			clearCaches();
			prepareRateMatrix();
			//rootLikelihood = computeLogLikelihood(N, T, posteriorRoot);
			
			double rootLikelihood = pruneAndPeel(column, root, 0);
			if(rootLikelihood > 0) {
				rootLogLikelihood = Math.log(rootLikelihood);
			} else {
				//System.out.println("\t\troot likelihood was 0");
				rootLogLikelihood = -100000;
			}
			//Compute log likelihood function
			

			//System.out.println("\t\tLikelihood of previous "+lastRootsLikelihood+" new " + rootLogLikelihood);

			if(rootLikelihood > 0 && lastRootsLikelihood< rootLogLikelihood) {
				fit.fittedPi = resultingPI;
				fit.fittedLikelihood = rootLogLikelihood;
				lastRootsLikelihood = rootLogLikelihood;
			} else {
				pi = fit.fittedPi;
			}
			
		}
		return fit;
	}

	private double computeLogLikelihood(Matrix N, double[] T, Matrix posteriorRoot) {
		double likelihood = 0;
		Matrix Q = getRateMatrix().times(omega);
		for(int a = 0; a < Q.getRowDimension(); a++) {
			likelihood += T[a] * Q.get(a, a) + posteriorRoot.get(a,0) * Math.log(pi.get(a, a));
			for(int b = 0; b < Q.getRowDimension(); b++) {
				if(b != a) {
					likelihood += N.get(a, b)*Math.log(Q.get(a,b));
				}
			}
		}
		return likelihood;
	}

	void clearCaches() {
		clearComputedTransitionsCache();
		clearNodeLikelihoods();
	}

	public double[] omegaEMIteration(Map<String, Matrix> leafValues, PhylogenyNode root, double newOmega, int window) {
		setOmega(newOmega); // update global omega to the value computed in prior iteration
		// E step (compute quantities required to get the gradient of the log likelihood bound function)
		nodeFittingParamMap = new HashMap<Integer, NodeLikelihoodParameters>(root.getNumberOfChildNodes());
		clearComputedTransitionsCache();
		double numOfTransitions = 0d;
		double totalTime = 0d; //it is close to expected time it stays on a transition
		double logLikelihood = 0;
		for(int site = 0; site < window; site++) {
			double likelihood = pruneAndPeel(leafValues, root, site);
			//System.out.println("Likelihood " + likelihood);
			// This matrix aids in the decomposition of the derivative of exp(Qt) w.r.t our parameters
			Matrix N = computeSufficientStatistics(root, likelihood);		

			//System.out.println("N: ");
			//N.print(alphabetSize, 10);

			// M-step, set omega.

			for(int i = 0; i < alphabetSize; i++) {
				for(int j = 0; j < alphabetSize; j++) {
					if(i == j) {
						totalTime -= N.get(i, i);
					} else {
						numOfTransitions += N.get(i, j);
					}
				}
			}
			logLikelihood += Math.log(likelihood);
		}
		//System.out.println("\ttransitions: " + numOfTransitions + ", totalTime: " + totalTime);
		newOmega = numOfTransitions * omega / totalTime;  
		//System.out.println("new omega " + newOmega);
		double [] data = {newOmega, numOfTransitions, totalTime, logLikelihood};
		return data;
	}

	private Matrix computeSufficientStatistics(PhylogenyNode root, double likelihood) {
		computeJMatrix(root);
		Iterator<NodeLikelihoodParameters> fitIt = nodeFittingParamMap.values().iterator();
		Matrix E = new Matrix(alphabetSize, alphabetSize, 0);
		while(fitIt.hasNext()) {
			//System.out.println("Iteration " + count++);
			NodeLikelihoodParameters fit = fitIt.next();
			if(fit.node.isRoot()) {
				continue;
			}
			Matrix projectedAlpha = iV.times(fit.alpha);
			Matrix projectedBeta  = V.transpose().times(fit.beta);
			/*
			System.out.println("\nnode " + fit.node.getID() + "-("+fit.node.getSeqName() + ")");
			System.out.println("\talpha (" + fit.alpha.get(0,0) +"," + fit.alpha.get(1,0) +"," + fit.alpha.get(2,0) +"," + fit.alpha.get(3,0) +")");
			System.out.println("\tbeta  (" + fit.beta.get(0,0) +"," + fit.beta.get(1,0) +"," + fit.beta.get(2,0) +"," + fit.beta.get(3,0) +")");
			System.out.println("\tprojected beta (" + projectedBeta.get(0,0) +"," + projectedBeta.get(1,0) +"," + projectedBeta.get(2,0) +"," + projectedBeta.get(3,0) +")");
			System.out.println("\tprojected alpha (" + projectedAlpha.get(0,0) +"," + projectedAlpha.get(1,0) +"," + projectedAlpha.get(2,0) +"," + projectedAlpha.get(3,0) +")");
			*/
			Matrix nodeE = projectedBeta.times(projectedAlpha.transpose()); // E = t(pBeta) * pAlpha)
			
			//System.out.println("J:");
			//fit.J.print(6, 4);
			nodeE.arrayTimesEquals(fit.J);
			E.plusEquals(nodeE);
		}

		E = E.times(1d/likelihood);
		//System.out.println("E/L:" );
		//E.print(alphabetSize, 12);

		Matrix N = iV.transpose().times(E).times(V.transpose());
		N.arrayTimesEquals(Q);

		return N;
	}
	
	public double computeLikelihood(Map<String, Matrix> leafValues,  PhylogenyNode root, int site) {
		double likelihood = 0;
		updateNodeTransition(root);
		
		peelModel(leafValues, root, site);
		
		NodeLikelihoodParameters rootFit = nodeFittingParamMap.get(root.getID());
		Matrix rootAlpha = rootFit.alpha;

		for(int i = 0; i < rootAlpha.getRowDimension(); i++) {
			likelihood += pi.get(i,i) * rootAlpha.get(i, 0); 
		}
		return likelihood;
	}

	public double pruneAndPeel(Map<String, Matrix> leafValues, PhylogenyNode root, int site) {
		double likelihood = computeLikelihood(leafValues, root, site);	
		pruneModel(root);
		return likelihood;
	}

	public void computeJMatrix(PhylogenyNode node) {
		//System.out.println("Called computeJMatrix node " + node.getID() + "-" + node.getSeqName() +"(" + node.getDistanceToParent() +")" );
		NodeLikelihoodParameters nodeFit = nodeFittingParamMap.get(node.getID());
		if(!node.isRoot()) {
			double dist = node.getDistanceToParent() * omega; 
			Matrix J = new Matrix(alphabetSize, alphabetSize);
			for(int i = 0; i < alphabetSize; i++) {
				for(int j = 0; j < alphabetSize; j++) {
					double di = D.get(i, i);
					double dj = D.get(j, j);
					if( Math.abs(di - dj) < 0.0001) { //If they are close, assume they are the same
						J.set(i, j, dist * Math.exp(di * dist));
					} else {
						J.set(i,j, (Math.exp(di * dist) - (Math.exp(dj * dist))) / (di - dj));
					}
				}
			}
			nodeFit.J = J;
			//System.out.println("node " + node.getID() + "-" + node.getSeqName() +" J");
			//J.print(alphabetSize, 10);
		} else {
			nodeFit.J = Matrix.identity(alphabetSize, alphabetSize);
		}
		
		if(node.getNumberOfChildNodes() > 0) {
			computeJMatrix(node.getChildNode1());
			computeJMatrix(node.getChildNode2());
		}
	}
	/**
	 * Implementation of Felsenstein's (dynamic programming, belief propagation) algorithm to
	 * estimate probabilities of seeing each base at the parent of the given node in the complement
	 * tree, obained by removing the node and its subtree.
	 * <b>Note:<\b> You must call peelModel before calling this method as it uses the peeling data 
	 * @param leafValues
	 * @param atNode
	 */
	
	private void pruneModel(PhylogenyNode atNode) {
		Matrix beta = new Matrix(alphabetSize, 1);			
		if(atNode.isRoot()) {

			for(int i = 0; i < alphabetSize; i++) {
				beta.set(i, 0, pi.get(i, i));
				//beta.set(i, 0, 1);
			}
			NodeLikelihoodParameters rootFit = nodeFittingParamMap.get(atNode.getID());
			//rootFit.beta = beta.times(1d/likelihood); // avoid small numbers.
			rootFit.beta = beta;
		} else {
			PhylogenyNode parent = atNode.getParent();
			PhylogenyNode sybling = parent.getChildNode1().getID() == atNode.getID() 
				? parent.getChildNode2() 
				: parent.getChildNode1();
			
			NodeLikelihoodParameters parentFit  = nodeFittingParamMap.get(parent.getID());
			
			NodeLikelihoodParameters syblingFit = nodeFittingParamMap.get(sybling.getID());	
			Matrix syblingSubtreeLikelihoods = syblingFit.transition.times(syblingFit.alpha);	
			//System.out.println("sybling alpha ");
			//syblingFit.alpha.print(alphabetSize, 5);
			Matrix parentComplementTreeLikelihoods = parentFit.transposedTransition.times(parentFit.beta);//parentFit.transition.times(parentFit.beta);
			//System.out.println("parent beta ");
			//parentFit.beta.print(alphabetSize, 5);
			beta = syblingSubtreeLikelihoods.arrayTimes(parentComplementTreeLikelihoods);
			NodeLikelihoodParameters atNodeFit = nodeFittingParamMap.get(atNode.getID());
			atNodeFit.beta = beta;
			
		}
		
		//System.out.println("Node " + atNode.getID() + "-" + atNode.getSeqName() + "(" + atNode.getDistanceToParent() + ") beta:");
		//beta.print(alphabetSize, 4);
		//System.out.print(" and alpha:");
		//nodeFittingParamMap.get(atNode.getID()).alpha.print(alphabetSize, 4);
		
		if(atNode.getNumberOfChildNodes() > 0) {
			pruneModel(atNode.getChildNode1());
			pruneModel(atNode.getChildNode2());
		}
	}
	/**
	 * Implementation of Felsenstein's peeling (dynamic programming, belief propagation) algorithm to 
	 * estimate the probabilities of seeing each of the letters of the alphabet at the 'atNode' given
	 * the evidence (leafValues).
	 * @param leafValues - observed values, the aligned (ungapped) column
	 * @param atNode - node whose subtree it roots will be used.
	 */
	private void peelModel(Map<String, Matrix> leafValues, PhylogenyNode atNode, int site) {
		Matrix alpha = new Matrix(alphabetSize, 1);
		
		NodeLikelihoodParameters fittingParams = nodeFittingParamMap.get(atNode.getID());
		
		if(atNode.isExternal()) {
			
			/*
			 * For leaves we observed the value, or the sampled base probabilities if missing data was sampled.!
			 */
			Matrix seqSeqMatrix = leafValues.get(atNode.getSeqName());
			if(seqSeqMatrix == null) {
				throw new IllegalStateException ("Alignment data does not contain a row for sequence " + atNode.getSeqName() + ". The aligned sequences must inlcude all sequences in the tree");
			}
			for(int i = 0; i < alphabetSize; i++) {
				alpha.set(i, 0, seqSeqMatrix.get(i,site));
			}
			
			//System.out.println("leaf node " + atNode.getSeqName() + " nucleotide prob dist:");
			//seqSeqMatrix.print(10, 8);

		} else {
			for(int i = 0; i < alphabetSize; i++)  {
				alpha.set(i, 0, 1);
			}

			// compute subtree likelihoods
			for(int child = 0; child < 2; child++) {
				PhylogenyNode childNode = atNode.getChildNode(child);
				peelModel(leafValues, childNode, site);
				NodeLikelihoodParameters childFit = nodeFittingParamMap.get(childNode.getID());
				Matrix likelihoodsForChild = childFit.transition.times(childFit.alpha);
				for(int j = 0; j < alphabetSize; j++) {
					alpha.set(j, 0, likelihoodsForChild.get(j, 0) * alpha.get(j,0)); 
				}
			}			
		}	
		//System.out.println("atNode " + atNode.getID() + "(" + atNode.getDistanceToParent() + ") alpha: ");
		fittingParams.alpha = alpha;
		//fittingParams.transition.print(alphabetSize, 10);
		//alpha.print(alphabetSize, 10);
	}
	
	public Map<Integer,NodeLikelihoodParameters> getNodeLikelihoodParameterMap() {
		return nodeFittingParamMap;
	}
	
	/**
	 * Computes the vector of transition probabilities for the given branch length and value of the end state.
	 * @param branchLenght distance to the end state.
	 * @return P the transition probabilities matrix for the given branch length.
	 */
	public Matrix computeTransitions(double branchLength) {
		Matrix result = transitionMatrixCache.get(branchLength);
		
		if(result == null) {
			Matrix DExp = D.copy();
			for(int i = 0; i < DExp.getRowDimension(); i++) {
				DExp.set(i, i, Math.exp(DExp.get(i, i) * omega * branchLength));
			}
		
			result = V.times(DExp);
		
			result = result.times(iV);
			
			transitionMatrixCache.put(branchLength, result);
		} 
		
		return result;
	}
	
	void clearComputedTransitionsCache() {
		transitionMatrixCache = new HashMap<Double, Matrix>();
	}
	
	public void setModelParameters(EvolutionaryModelParameters parameters) throws IllegalArgumentException{
		this.parameters = parameters;
		prepareRateMatrix();
		alphabetSize = parameters.getBackgroundNucleotideFreqs().length;
	}

	public EvolutionaryModelParameters getModelParameters() {
		return parameters;
	}
	
	private void updateNodeTransition (PhylogenyNode atNode) {
		NodeLikelihoodParameters nodeFittingParameters = new NodeLikelihoodParameters(atNode);
		nodeFittingParamMap.put(atNode.getID(), nodeFittingParameters);
		
		
		Matrix branchLengthTransitions = 
			computeTransitions(atNode.isRoot() ? 0 : atNode.getDistanceToParent());

		nodeFittingParameters.transition = branchLengthTransitions;
		nodeFittingParameters.transposedTransition = branchLengthTransitions.transpose();
		
		if(atNode.getNumberOfChildNodes() > 0) {
			updateNodeTransition(atNode.getChildNode1());
			updateNodeTransition(atNode.getChildNode2());
		}
	}
	
	private void prepareRateMatrix() {
		Q = getRateMatrix() ; // use temporary Q matrix if available
		
		EigenvalueDecomposition ed = Q.eig();
		
		if(!ed.getD().isDiagonal()) {
			new IllegalArgumentException("rate matrix  has non real eigenvalues: " + ed.getD() + " this is not supported in this version");
		}
		
		V = ed.getV();
		D = ed.getD();
		/*
		ArrayList<Double> eigens = new ArrayList<Double>(4);
		eigens.add(D.get(0,0));
		eigens.add(D.get(1,1));
		eigens.add(D.get(2,2));
		eigens.add(D.get(3,3));
		
		Collections.sort(eigens);
		Matrix nV = new Matrix(V.getRowDimension(), V.getColumnDimension());
		for(int i = 0; i < V.getColumnDimension(); i++) {
			for (int j = 0; j < V.getColumnDimension(); j++) {
				if(eigens.get(i) == D.get(j,j)) {
					int sign = V.get(0,j) < 0 ? -1 : 1;
					for(int k = 0; k < V.getRowDimension(); k++) {
						nV.set(k, i, V.get(k,j) * sign);
					}
				}
			}
		}
		V = nV;
		D.set(0, 0,eigens.get(0));
		D.set(1, 1,eigens.get(1));
		D.set(2, 2,eigens.get(2));
		D.set(3, 3,eigens.get(3));
		*/

		
		iV = V.inverse();

		
		/*
		System.out.println("Q or R");
		Q.print(Q.getRowDimension(), 10);
		System.out.println("V");
		V.print(V.getRowDimension(), 10);
		System.out.println("iV");
		iV.print(iV.getRowDimension(), 10);
		System.out.println("D");
		D.print(D.getRowDimension(), 10);
		
		System.out.println("Should be the Q (V*D*iV):");
		V.times(D).times(iV).print(iV.getRowDimension(), 7);
		*/
		
	}
	
	public Phylogeny getTree() {
		return parameters.getTree();
	}
	
	public static class NodeLikelihoodParameters {
		public Matrix transposedTransition;
		PhylogenyNode node;
		Matrix transition;
		Matrix alpha; // Stores probability of finding each letter at node given leaf values
		Matrix beta;  // Stores probability of finding each base at given node given the complement of the subtree that this node roots.
		Matrix J;     // The J matrix obtained to simplify derivative of exp(Qt) w.r.t model hidden parameters
		
		public NodeLikelihoodParameters(PhylogenyNode node) {
			this.node = node;
		}
		
		/**
		 * Estimates base probabilities at the node. Originally it was implemented
		 * so that the caller would provide the tree likelihood (to avoid the extra addition)
		 * but small rounding differences had a net impact of 0.05 the worse cases (it was in general
		 * in the order of 10 ^ -20). To avoid this we recompute the total likelihood by adding up 
		 * the likelihoods of each of the nucleotides.
		 * @param totalTreeLikelihood
		 * @return
		 */
		public Matrix computeProbabilityOfLetterAtNode(/*double totalTreeLikelihood*/) {
			Matrix Qa = new Matrix(transition.getRowDimension(), 1);
			double likelihood = 0;
			for (int i = 0; i < transition.getRowDimension(); i++) {
				
				for (int j = 0; j < transition.getColumnDimension(); j++) {
					Qa.set(i, 0, Qa.get(i, 0) + beta.get(j, 0) * transition.get(j, i));
				}
				Qa.set(i, 0, Qa.get(i,0)  * alpha.get(i,0));
				likelihood += Qa.get(i, 0);
			}
			//System.out.println( "Tree likelihood " + totalTreeLikelihood + " node computed likelihood " + likelihood + " Difference in likelihoods " + (totalTreeLikelihood - likelihood));
			return Qa.times(1d/likelihood); 
		}
	}
	
	private short sampleFromAlphabet(double[] frequencies) {
		Random r = new Random();
		double draw = r.nextDouble();
		//System.out.println("draw " + draw);
		short base = -1;
		short i = 0;
		double cummulativeProb = 0;
		while( base < 0) {
			cummulativeProb += frequencies[i];
			//System.out.println("\tcumm dist: " + cummulativeProb);
			if(draw <= cummulativeProb) {
				//System.out.println("base " + i);
				base = i;
			}
			i++;
		}
		return base;
	}

	//TODO: read Korbel et al. to formalize binning of data. For now using 10 bins 
	public void setOmegaDistByTreeLength(String neutralOmegaDistributionByTreeLengthFile) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(neutralOmegaDistributionByTreeLengthFile));
		String line = null;
		long valueCount = 0;
		LinkedHashMap<Double, List<Double>> treeLengthOmegasMap = new LinkedHashMap<Double, List<Double>>();
		while((line = br.readLine()) != null) {
			String [] lineInfo = line.split("\t");
			double treeLength = Double.parseDouble(lineInfo[0]);
			String [] omegaStrs = lineInfo[1].split(",");
			valueCount += omegaStrs.length;
			List<Double> omegaList = new ArrayList<Double>(omegaStrs.length);
			for(int i = 0; i < omegaStrs.length; i++) {
				omegaList.add(Double.parseDouble(omegaStrs[i]));
			}
			
			treeLengthOmegasMap.put(treeLength, omegaList);
		}		
		br.close();
		binStatistics = new LinkedHashMap<Double, Frequency>();
		double binDataSize = Math.ceil(valueCount/10);
		System.out.println("Loaded " + valueCount + " omegas estimated from neutral sequence bin sizes will be of " + binDataSize);
		Iterator<Double> treeLengthIt = treeLengthOmegasMap.keySet().iterator();
		double currentBinSize = 0;
		Stack<Frequency> frequencyStack = new Stack<Frequency>();
		Stack<Double> upperLimmitFrequencyStack = new Stack<Double>();
		double priorTreeLength = 0;
		while(treeLengthIt.hasNext()) {
			double treeLength = treeLengthIt.next();

			if(frequencyStack.isEmpty() || currentBinSize > binDataSize) {
				if(!frequencyStack.isEmpty()) {
					upperLimmitFrequencyStack.push(priorTreeLength);
					System.out.println("Pushing new limming: " + priorTreeLength);
				}
				Frequency freq = new Frequency();
				frequencyStack.push(freq);
				currentBinSize = 0;
				//System.out.println("Pushing new Frequency ");
			}
			
			Iterator<Double> omegaIt = treeLengthOmegasMap.get(treeLength).iterator();
			while(omegaIt.hasNext()) {
				currentBinSize++;
				frequencyStack.peek().addValue(omegaIt.next());
			}
			priorTreeLength = treeLength;
		}
		
		upperLimmitFrequencyStack.push(priorTreeLength);
		
		assert(upperLimmitFrequencyStack.size() == frequencyStack.size());
		
		for(int i = 0; i < upperLimmitFrequencyStack.size(); i++) {
			binStatistics.put(upperLimmitFrequencyStack.get(i), frequencyStack.get(i));
		}
		/*
		System.out.println("Example, p values for and omega 0.1:5 \n");
		Iterator<Double> upperLimitIt = binStatistics.keySet().iterator();
		while(upperLimitIt.hasNext()) {
			double length = upperLimitIt.next();
			Frequency f = binStatistics.get(length);
			BufferedWriter bw = new BufferedWriter(new FileWriter("omega_dist_lengths_up_to"+String.valueOf(length).substring(0,4) +".txt"));
			Iterator<Double> valIt = f.valuesIterator();
			while(valIt.hasNext()) {
				bw.write(String.valueOf(valIt.next()));
				bw.newLine();
			}
			
			bw.close();
			
			System.out.println("\tlength " + length + " p-val " + f.getCumPct(new Double(0.1d)));
		}
		*/
	}

	public double getPVal(OmegaFit fit, boolean useLogOddsLikelihood) {
		double pVal = -1;
		
		Iterator<Double> treeLengthIt = binStatistics.keySet().iterator();
		while(treeLengthIt.hasNext()) {
			double length = treeLengthIt.next();
			if(length >= fit.getTreeLength()) {
				Frequency f = binStatistics.get(length);
				//pVal = useLogOddsLikelihood ? f.getCumPct(new Double(fit.getOmega()));
				break;
			}
		}
		return pVal;
	}

	private void clearNodeLikelihoods() {
		nodeFittingParamMap.clear();		
	}

	public void setNeutralOmegaDistribution(String neutralOmegaFile) throws IOException {
		Frequency omegaDist = new Frequency();
		
		BufferedReader br = new BufferedReader(new FileReader(neutralOmegaFile));
		String line = null;
		double lastTreeLength = -1;
		int lineNum = 0;
		while((line = br.readLine()) != null) {
			lineNum++;
			if(line.startsWith("#")) {
				continue;
			}
			
			String [] lineInfo = line.split("\\s");
			
			double omega = Double.parseDouble(lineInfo[1]);
			double treeLength = Double.parseDouble(lineInfo[2]);
			if(lineNum >1 && lastTreeLength != treeLength) {
				throw new IllegalArgumentException("The neutral omega file must all be computed on the same tree, found two different tree lentgths in lines " + (lineNum - 1) + " and line " + lineNum);
			}
			
			lastTreeLength = treeLength;
			omegaDist.addValue(omega);
		}
		br.close();
		binStatistics = new LinkedHashMap<Double, Frequency>(1);
		binStatistics.put(lastTreeLength, omegaDist);
	}

	private Matrix extractEquilibriumFromRateMatrix() {
		/*
		Matrix invDiagPi = new Matrix(pi.getRowDimension(), pi.getColumnDimension());
		for(int i = 0; i < invDiagPi.getRowDimension(); i++) {
			invDiagPi.set(i, i, 1d/pi.get(i, i));
		}
		*/
		Matrix rateMatrix = parameters.getRateMatrix();
		
		//Jan 22/09 changed to test Or Zuk's suggestion

		Matrix R = rateMatrix.times(pi.inverse()); //Original decomposition.
		//Matrix R = pi.times(rateMatrix); //This is Or's suggestion.
		
		//System.out.println("R should be symetric, is it? " + R.isSymmetric());
		if(!R.isSymmetric()) {
			//R.print(4, 6);
			R.plusEquals(R.transpose()).timesEquals(1d/2d);
			//System.out.println("but now it is:");
			//R.print(4,6);
		}
		
		return R;
	}
	
	public static class OmegaFit {
		private double omega;
		private double initialLogLikelihood;
		private double fittedLogLikelihood;
		private int numOfIterations;
		private double treeLength;
		private double transitions;
		private double totalTime;
		private double pVal;
		private Annotation region;
		private double LOD;
		private boolean LODSet=false;
		
		public OmegaFit(){}
		
		public OmegaFit(LightweightGenomicAnnotation region, double omega, double LOD, double p, double treeLength){
			this.region=region;
			this.omega=omega;
			this.LOD=LOD;
			this.LODSet=true;
			this.pVal=p;
			this.treeLength=treeLength;
		}
		
		public void setRegion(Annotation region){this.region=region;}
		
		public double getFittedLogLikelihood() {
			return fittedLogLikelihood;
		}
		public double getLogOddsScore() {
			// TODO Auto-generated method stub
			if(LODSet){return LOD;}
			return omega < 1 ? 2*(fittedLogLikelihood - initialLogLikelihood) : -2*(fittedLogLikelihood - initialLogLikelihood);
		}
		public double getInitialLogLikelihood() {
			return initialLogLikelihood;
		}
		public int getNumOfIterations() {
			return numOfIterations;
		}
		public double getOmega() {
			return omega;
		}
		public double getTreeLength() {
			return treeLength;
		}
		public Annotation getRegion(){
			return region;
		}		
		public void setTreeLength(double treeLength) {
			this.treeLength = treeLength;
		}
		
		public double getPVal() { return pVal;}
		
		public double getTotalTime() {
			return totalTime;
		}
		public double getTransitions() {
			return transitions;
		}
		protected void setFittedLogLikelihood(double fittedLogLikelihood) {
			this.fittedLogLikelihood = fittedLogLikelihood;
		}
		protected void setInitialLogLikelihood(double initialLogLikelihood) {
			this.initialLogLikelihood = initialLogLikelihood;
		}
		protected void setNumOfIterations(int numOfIterations) {
			this.numOfIterations = numOfIterations;
		}
		protected void setOmega(double omega) {
			this.omega = omega;
		}
		protected void setTotalTime(double totalTime) {
			this.totalTime = totalTime;
		}
		protected void setTransitions(double transitions) {
			this.transitions = transitions;
		}
		
		
	}
	
	public static class PiFit implements Fit{
		
		int position;
		private Matrix fittedPi;
		int numOfIterations;
		private double fittedLikelihood;
		private double unfittedLikelihood;
		private double logLikelihoodRatio;
		private double pValue;
		private double treeLength;
		private int treeBit;
		
		public PiFit() {
			super();
		}
		
		public PiFit(String[] lineInfo) {
			super();
			int data = 0;
			position = Integer.parseInt(lineInfo[data++]);
			fittedPi = new Matrix(4,1);
			for(int i = 0; i < 4; i++) {
				fittedPi.set(i,0,Double.parseDouble(lineInfo[data++]));
			}
			logLikelihoodRatio = Double.parseDouble(lineInfo[data++]);
			treeLength = Double.parseDouble(lineInfo[data++]);
			
			if(lineInfo.length > data) {
				treeBit  = Integer.parseInt(lineInfo[data++]);
			}
			
		}

		public Matrix getPI() {
			Matrix colPI = new Matrix(fittedPi.getRowDimension(), 1);
			for(int i = 0; i < colPI.getRowDimension(); i++) {
				colPI.set(i, 0, fittedPi.get(i,i));
			}
			return colPI;
		}
		
		public int getNumberOfIterations() { return numOfIterations;}
		
		public double getLogLikelihoodRatio() {
			return logLikelihoodRatio == 0 ? 2*(fittedLikelihood - unfittedLikelihood) : logLikelihoodRatio;
		}

		public int getPosition() {
			return position;
		}

		public void setPosition(int position) {
			this.position = position;
		}

		public Matrix getFittedPi() {
			return fittedPi;
		}

		public double getPValue() {
			return pValue;
		}

		public void setTreeLength(double treeDist) {
			this.treeLength = treeDist;
			
		}
		
		public double getInformationContent() {
			double ic = 2;
			for (int i = 0; i < 4; i++) {
				if(fittedPi.get(i, 0) > 0) {
					ic = ic + fittedPi.get(i, 0)*MathUtil.log2(fittedPi.get(i, 0));
				}
			}
			return ic;
		}
		
		public double getTreeLength() { return this.treeLength;}

		public void setTreeBit(int treeBitVal) {
			this.treeBit = treeBitVal;
			
		}

		public int getTreeBit() {
			return treeBit;
		}
		
		public String toString() {
			StringBuilder sb = new StringBuilder();
			
			sb.append(String.valueOf(getPosition()));
			
			for(int i = 0; i < fittedPi.getRowDimension(); i++) {
				sb.append("\t").append(numberFormat.format(fittedPi.get(i, 0)));
			}
			
			sb.append("\t")
				.append((numberFormat.format(getLogLikelihoodRatio())))
				.append("\t")
				.append(numberFormat.format(getTreeLength()))
				.append("\t")
				.append(String.valueOf(getTreeBit()));
			
			return sb.toString();
		}
		
	}

	public void resetPI() {
		pi = new Matrix(parameters.getBackgroundNucleotideFreqs().length, parameters.getBackgroundNucleotideFreqs().length);
		for(int i = 0; i < parameters.getBackgroundNucleotideFreqs().length; i++) {
			pi.set(i, i, parameters.getBackgroundNucleotideFreqs()[i]);
		}		
	}

	public EvolutionaryModel copy() {
		EvolutionaryModel copy = new EvolutionaryModel(parameters.copy());
		copy.alphabetSize = alphabetSize;
		copy.binStatistics = binStatistics;
		copy.pi = pi;
		
		return copy;
	}

	public EvolutionaryModelParameters getParameters() {
		return parameters;
	}

	public void setParameters(EvolutionaryModelParameters parameters) {
		this.parameters = parameters;
	}

	public int getAlphabetSize() {
		return alphabetSize;
	}

	public void setAlphabetSize(int alphabetSize) {
		this.alphabetSize = alphabetSize;
	}

	public Matrix getPi() {
		return pi;
	}

	public void setPi(Matrix pi) {
		if(pi.getColumnDimension() == 1) {
			this.pi = new Matrix(pi.getRowDimension(), pi.getRowDimension());
			for(int i = 0; i < pi.getRowDimension(); i++) {
				this.pi.set(i, i, pi.get(i,0));
			}
		} else {
			this.pi = pi;
		}
		prepareRateMatrix();
	}

	private void addPseudoCounts(Matrix pi) {
		if(!pi.isDiagonal()) {
			throw new IllegalArgumentException("Expected diagonal matrix representing the stationary distribution but it was not a symmetric matrix");
		}
		
		int numOfZeroes = 0;
		for(int i = 0; i < pi.getRowDimension(); i++) {
			if(pi.get(i, i) == 0 ) {
				numOfZeroes++;
			}
		}
		
		if(numOfZeroes > 0) {
			for(int i = 0; i < pi.getRowDimension(); i++) {
				if(pi.get(i, i) == 0 ) {
					pi.set(i, i, TINY_DIFF);  
				} else {
					pi.set(i, i, (TINY_DIFF * numOfZeroes)/(double)(pi.getRowDimension() - numOfZeroes));
				}
			}			
		}
	}

	public void pruneTree(List<String> ignoreList) {
		Phylogeny tree = parameters.getTree();
		Phylogeny prunnedTree = ConservationUtils.pruneTree(ignoreList, tree);		
		parameters.setTree(prunnedTree);
	}

	/**
	 * Assumes the given sequence as ancestral, the runs the markov process forward to 
	 * generate sequences at ach node.
	 * @param sequence - ancestral sequence
	 * @return Phylogeny - with the same topology as the model's but its node names have been replaced with the evolved sequence.
	 */
	public Phylogeny evolve(String sequence) {
		Sequence seq = new Sequence(sequence);
		seq.setSequenceBases(sequence);
		Matrix encodedSeq = seq.encodeSequenceAsVector();
		Phylogeny tree = parameters.getTree().copy();
		PhylogenyNode root = tree.getRoot();
		
		root.setSeqName(sequence);
		evolveForwardFromNode(root, encodedSeq);
		
		return tree;
	}

	/**
	 * Starting a a given node it generates "evolved" sequence for all its children, setting the corresponding 
	 * children's node name to the resulting decoded sequence
	 * TODO: Create a new Phylogeny object that can hold the encoded sequence as data. This would avoid having to 
	 *       Decode the sequence and loosing all information.
	 * @param root
	 * @param encodedSeq
	 */
	private void evolveForwardFromNode(PhylogenyNode node, Matrix encodedSeq) {
		PhylogenyNode child1 = node.getChildNode1();
		Matrix evolvedSequenceAtChild1 = runMarkovProcessForward(child1, encodedSeq);
		String sequenceAtChild1Str = sampleSequence(evolvedSequenceAtChild1);
		child1.setSeqName(sequenceAtChild1Str);
		if(child1.getNumberOfChildNodes() > 0) {
			Sequence child1Seq = new Sequence("child1");
			child1Seq.setSequenceBases(sequenceAtChild1Str);
			evolveForwardFromNode(child1, child1Seq.encodeSequenceAsVector());
		}
		
		PhylogenyNode child2 = node.getChildNode2();
		Matrix evolvedSequenceAtChild2 = runMarkovProcessForward(child2, encodedSeq);
		String sequenceAtChild2Str = sampleSequence(evolvedSequenceAtChild2);
		child2.setSeqName(sequenceAtChild2Str);		
		if(child2.getNumberOfChildNodes() > 0) {
			Sequence child2Seq = new Sequence("child2");
			child2Seq.setSequenceBases(sequenceAtChild2Str);
			evolveForwardFromNode(child2, child2Seq.encodeSequenceAsVector());
		}
		
	}

	private String sampleSequence(Matrix encodedSeq) {
		short [] sampledSequence = new short[encodedSeq.getColumnDimension()];
		Random r = new Random();		
		for(int j = 0 ; j < encodedSeq.getColumnDimension(); j++) {
			double prob = r.nextDouble();
			double totalProb = 0;
			short base = -1;
			for(short i = 0; i < encodedSeq.getRowDimension(); i++) {
				totalProb = totalProb + encodedSeq.get(i, j);
				if(totalProb >= prob) {
					base =  i;
					break;
				}
			}
			if(base == -1) {
				throw new RuntimeException ("Shit!, did not find a letter with enough probability at column " + j + " totalProb " + totalProb + " prob " + prob);
			}
			sampledSequence[j] = base;
		}
		//System.out.println("sampled seq " + Sequence.decode(sampledSequence) );
		return Sequence.decode(sampledSequence);
	}

	private Matrix runMarkovProcessForward(PhylogenyNode child, Matrix parentEncodedSeq) {
		Matrix branchLengthTransitions = 
			computeTransitions(child.isRoot() ? 0 : child.getDistanceToParent());
		/*
		System.out.print("Transition ");
		branchLengthTransitions.print(3, 6);
		System.out.print("parent sequence ");
		parentEncodedSeq.print(3, 6);
		System.out.print("child seq ");
		branchLengthTransitions.transpose().times(parentEncodedSeq).print(3, 6);
		*/
		return branchLengthTransitions.transpose().times(parentEncodedSeq);
	}

	public void changeBackground(double pA, double pC, double pG, double pT) {
		Matrix newPi = new Matrix(4,1);
		newPi.set(0,0,pA);
		newPi.set(1,0,pC);
		newPi.set(2,0,pG);
		newPi.set(3,0,pT);
		changeBackground(newPi);
	}
	
	public void changeBackground(Matrix newPi) {
		RPIDecomposition();
		setPi(newPi);
		parameters.setBackgroundNucleotideFreqs(newPi); 
		parameters.setRateMatrix(getRateMatrix());
	}
}
