package broad.core.math;

import java.util.Arrays;
import java.util.Map;
import java.util.TreeMap;

import Jama.Matrix;
import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;


public class FDRDistribution {

	
	double[] sortedAbsObservations;
	Matrix sortedAbsPermutations;
	double[] absFDRs;
	
	//double[] averagePermutations;
	Matrix sortedPermutations;
	Matrix negativeSortedPermutations;
	double[] sortedObservations;
	double[] negativeSortedObservations;
	double[] posFDRs;
	double[] negFDRs;
	double alpha;
	private Map<Double, Integer> indexMap;
	private Map<Double, Integer> absoluteIndexMap;
	
	IntervalTree<Integer>[] posTrees=null;
	IntervalTree<Integer>[] negTrees=null;
	IntervalTree<Integer>[] absTrees=null;
	
	//TODO Merge with Moran's constructor with normalize false
	public FDRDistribution(double[] observed, Matrix permutations, double alpha){
		Matrix absPermutations=Statistics.absoluteValue(permutations);
		this.sortedAbsPermutations=absPermutations.sortColumns();
		this.sortedAbsObservations=Statistics.absoluteValue(observed);
		Arrays.sort(this.sortedAbsObservations);
		
		this.sortedPermutations=permutations.sortColumns();
		Arrays.sort(observed);
		this.sortedObservations=observed;
		this.negativeSortedObservations=flipSort(observed);
		this.negativeSortedPermutations=flipSort(this.sortedPermutations);
		this.alpha=alpha;
		
		posTrees=makeTrees(sortedPermutations);
		negTrees=makeTrees(negativeSortedPermutations);
		absTrees=makeTrees(sortedAbsPermutations);
		
		
		//System.err.println("started computing FDRs");
		this.posFDRs=this.computeAllFDRs(sortedObservations, posTrees);
		this.negFDRs=this.computeAllNegFDRs(negativeSortedObservations, negTrees);
		this.absFDRs=this.computeAllFDRs(sortedAbsObservations, absTrees);
		indexMap=createIndexMap(this.sortedObservations);
		absoluteIndexMap=createIndexMap(this.sortedAbsObservations);
		//average matrix
		//this.averagePermutations=average(sortedPermutations);
	}
	
	public FDRDistribution(double[] observed, FDRDistribution dist){
		this.sortedAbsObservations=Statistics.absoluteValue(observed);
		Arrays.sort(this.sortedAbsObservations);
		Arrays.sort(observed);
		this.sortedObservations=observed;
		this.negativeSortedObservations=flipSort(observed);
		
		this.posTrees=dist.getPosTrees();
		this.negTrees=dist.getNegTrees();
		this.absTrees=dist.getAbsTrees();
		
		this.posFDRs=this.computeAllFDRs(sortedObservations, posTrees);
		this.negFDRs=this.computeAllNegFDRs(negativeSortedObservations, negTrees);
		this.absFDRs=this.computeAllFDRs(sortedAbsObservations, absTrees);
		
		indexMap=createIndexMap(this.sortedObservations);
		absoluteIndexMap=createIndexMap(this.sortedAbsObservations);
		
	}
	
	
	private IntervalTree<Integer>[] getAbsTrees() {
		if(this.absTrees==null){
			this.absTrees=makeTrees(this.sortedAbsPermutations);
		}
		return this.absTrees;
	}

	private IntervalTree<Integer>[] getNegTrees() {
		if(this.negTrees==null){
			this.negTrees=makeTrees(this.negativeSortedPermutations);
		}
		return this.negTrees;
	}

	private IntervalTree<Integer>[] getPosTrees() {
		if(this.posTrees==null){
			this.posTrees=makeTrees(this.sortedPermutations);
		}
		return this.posTrees;
	}

	//TODO Moran: Make sure this actually does some thing
	//Normalize scores in each row according to the permuted values in that row.
	public FDRDistribution(double[] observed, Matrix permutations, double alpha, boolean normalizeRows){
		
		if (normalizeRows==false ) {this.initialize(observed, permutations, alpha);}
		else{
			double[] normObs=normalizeScore(observed,permutations);
			Matrix normPerm= new Matrix(permutations.getNumRows(),permutations.getNumColumns());
			for (int i=0; i<permutations.getNumColumns(); i++){
				normPerm.setColumn(i, normalizeScore(permutations.getColumn(i),permutations));
			}
			//this.initialize(normObs, normPerm, alpha);
			this.sortedPermutations=permutations.sortColumns();
			Arrays.sort(observed);
			this.sortedObservations=observed;
			this.negativeSortedObservations=flipSort(observed);
			this.negativeSortedPermutations=flipSort(this.sortedPermutations);
			//computeCriticalValue(alpha);
			this.alpha=alpha;
			//System.err.println("started computing FDRs");
			this.posFDRs=this.computeAllFDRs(sortedObservations, sortedPermutations);
			this.negFDRs=this.computeAllNegFDRs(negativeSortedObservations, negativeSortedPermutations);
			indexMap=createIndexMap(this.sortedObservations);
		
		}
	}
	
	

	private void initialize (double[] observed, Matrix permutations, double alpha){
		this.sortedPermutations=permutations.sortColumns();
		Arrays.sort(observed);
		this.sortedObservations=observed;
		this.negativeSortedObservations=flipSort(observed);
		this.negativeSortedPermutations=flipSort(this.sortedPermutations);
		//computeCriticalValue(alpha);
		this.alpha=alpha;
		//System.err.println("started computing FDRs");
		this.posFDRs=this.computeAllFDRs(sortedObservations, sortedPermutations);
		this.negFDRs=this.computeAllNegFDRs(negativeSortedObservations, negativeSortedPermutations);
		indexMap=createIndexMap(this.sortedObservations);
		
	}
	
	private double[] average(Matrix sortedPermutations) {
		double[] rtrn=new double[sortedPermutations.getRowDimension()];
		
		for(int i=0; i<rtrn.length; i++){
			double[] vals=sortedPermutations.getRow(i);
			rtrn[i]=Statistics.average(vals);
		}
		
		return rtrn;
	}

	private Map<Double, Integer> createIndexMap(double[] sortedObservations) {
		Map<Double, Integer> rtrn=new TreeMap<Double, Integer>();
		
		for(int i=0; i<sortedObservations.length; i++){
			rtrn.put(sortedObservations[i], i);
		}
		
		return rtrn;
	}

	private double[] computeAllNegFDRs(double[] negativeSortedObservations2, Matrix negativeSortedPermutations2) {
		double[] fdr=this.computeAllFDRs(negativeSortedObservations2, negativeSortedPermutations2);
		double[] rtrn=new double[fdr.length];
		
		int counter=0;
		for(int i=fdr.length-1; i>=0; i--){
			rtrn[counter++]=fdr[i];
		}
		
		return rtrn;
	}
	
	private double[] computeAllNegFDRs(double[] negativeSortedObservations2, IntervalTree<Integer>[] trees) {
		double[] fdr=this.computeAllFDRs(negativeSortedObservations2, trees);
		double[] rtrn=new double[fdr.length];
		
		int counter=0;
		for(int i=fdr.length-1; i>=0; i--){
			rtrn[counter++]=fdr[i];
		}
		
		return rtrn;
	}

	private double[] flipSort(double[] observed) {
		double[] rtrn=new double[observed.length];
		
		int counter=0;
		for(int i=observed.length-1; i>=0; i--){rtrn[counter++]=-observed[i];}
		
		//Arrays.sort(rtrn);
		return rtrn;
	}
	
	private Matrix flipSort(Matrix observed) {
		Matrix rtrn=new Matrix(observed.getNumRows(), observed.getNumColumns());
		for(int i=0; i<observed.getNumColumns(); i++){
			double[] column=observed.getColumn(i);
			rtrn.setColumn(i, flipSort(column));
		}
		return rtrn;
	}

	/*private void computeCriticalValue(double alpha) {
		//System.err.println(sortedObservations.length);
		
		for(int k=0; k<sortedObservations.length; k++){
			if(sortedObservations[k]>0){
				double fdr=computeExactPositiveFDR(this.sortedObservations[k]);
				System.err.println("positive "+sortedObservations[k]+" "+fdr);
				if(fdr<=alpha){
					this.posCriticalValue=this.sortedObservations[k];
					this.posCriticalIndex=k;
					break;
				}
			}
		}
		
		for(int k=sortedObservations.length-1; k>=0; k--){
			if(sortedObservations[k]<0){
				double fdr=computeExactNegativeFDR(this.sortedObservations[k]);
				System.err.println("negative "+sortedObservations[k]+" "+fdr);
				if(fdr<=alpha){
					this.negCriticalValue=this.sortedObservations[k];
					this.negCriticalIndex=k;
					break;
				}
			}
		}
	}*/

	

	public double getFDR(double k){
		int index=indexMap.get(k);
		return Math.min(this.posFDRs[index], this.negFDRs[index]);
	}
	
	public double getAbsFDR(double k){
		int index=this.absoluteIndexMap.get(Math.abs(k));
		return this.absFDRs[index];
	}
	
	private double[] computeAllFDRs(double[] sortedObservations, Matrix sortedPermutations){
		double[] exactFDRs=assignFDR(sortedObservations, sortedPermutations);
		double[] fdr=adjustFDRs(exactFDRs);
		return fdr;
	}
	
	private double[] computeAllFDRs(double[] sortedObservations, IntervalTree<Integer>[] trees){
		double[] exactFDRs=assignFDR(sortedObservations, trees);
		double[] fdr=adjustFDRs(exactFDRs);
		return fdr;
	}
	
	
	private double[] adjustFDRs(double[] exactFDRs) {
		double[] rtrn=new double[exactFDRs.length];
		
		//go through each cell and replace with min value
		//TODO: It might not be enough to back track since there might be values that we never observe for which the FDR might be very low
		double minSoFar=1.0;
		for(int i=0; i<exactFDRs.length; i++){
			//goes from lower scores to higher
			minSoFar=Math.min(minSoFar, exactFDRs[i]);
			rtrn[i]=minSoFar;
		}
		
		return rtrn;
	}

	/*private double[] assignFDR(double[] sortedObservations, Matrix sortedPermutations) {
		//This matrix is for each observed value the number of random greater than it (each column is a permutation)
		
		//Step 1: For each observed value compute the number of permutations greater than it
		Matrix greaterThanNums=new Matrix(sortedObservations.length, sortedPermutations.getColumnDimension());
		for(int i=0; i<sortedPermutations.getColumnDimension(); i++){
			//System.err.println("Permutation "+i);
			greaterThanNums.setColumn(i, computeLargerForArray(sortedObservations, sortedPermutations.getColumn(i)));
		}
		
		double[] fdrs=new double[sortedObservations.length];
		
		double[] observedGreater=computeLargerForArray(sortedObservations, sortedObservations);
		
		//Step 2: Compute the exact FDR
		for(int i=0; i<observedGreater.length; i++){
			double[] perms=greaterThanNums.getRow(i);
			double expectedNumberGreaterThanK=Statistics.average(perms);
			
			double estimateOfMeanGreaterThanK=expectedNumberGreaterThanK;
			//double estimateOfMeanGreaterThanK=estimateMeanGreaterThanK(expectedNumberGreaterThanK, sortedObservations.length, observedGreater[i]);
			if(estimateOfMeanGreaterThanK==0){estimateOfMeanGreaterThanK=((double)1/greaterThanNums.getColumnDimension());}
			
			double fdr=estimateOfMeanGreaterThanK/observedGreater[i];
			fdrs[i]= Math.min(fdr, 1);
		}
		
		return fdrs;
	}*/
	
	//Needed to make this MUCH faster
	private double[] assignFDR(double[] sortedObservations, Matrix sortedPermutations) {
		IntervalTree<Integer>[] trees=makeTrees(sortedPermutations);
		return assignFDR(sortedObservations, trees);
	}
	
	//Needed to make this MUCH faster
	private double[] assignFDR(double[] sortedObservations, IntervalTree<Integer>[] trees) {
		//This matrix is for each observed value the number of random greater than it (each column is a permutation)
		
		//Step 1: For each observed value compute the number of permutations greater than it
		Matrix greaterThanNums=new Matrix(sortedObservations.length, trees.length);
		for(int i=0; i<trees.length; i++){
			greaterThanNums.setColumn(i, computeLargerForArray(sortedObservations, trees[i]));
		}
		
		double[] fdrs=new double[sortedObservations.length];
		
		
		//Step 2: Compute the exact FDR
		for(int i=0; i<sortedObservations.length; i++){
			double observedGreater=sortedObservations.length-i;
			
			double[] perms=greaterThanNums.getRow(i);
			double expectedNumberGreaterThanK=Statistics.average(perms);
			
			double estimateOfMeanGreaterThanK=expectedNumberGreaterThanK;
			//double estimateOfMeanGreaterThanK=estimateMeanGreaterThanK(expectedNumberGreaterThanK, sortedObservations.length, observedGreater[i]);
			if(estimateOfMeanGreaterThanK==0){estimateOfMeanGreaterThanK=((double)1/greaterThanNums.getColumnDimension());}
			
			double fdr=estimateOfMeanGreaterThanK/observedGreater;
			//if(fdr<.2){System.err.println(fdr);}
			fdrs[i]= Math.min(fdr, 1);
		}
		
		return fdrs;
	}

	private double[] computeLargerForArray(double[] sortedObservations,IntervalTree<Integer> tree) {
		double[] rtrn=new double[sortedObservations.length];
		
		for(int i=0; i<sortedObservations.length; i++){
			//TODO Should get a point and get max/min
			Node<Integer> val=tree.min(new Double(sortedObservations[i]*1000).intValue(), new Double(sortedObservations[i]*1000).intValue());
			int j=-1;
			if(val==null){j=sortedObservations.length;}
			else{j=val.getValue();}
			rtrn[i]=sortedObservations.length-j;
		}
		
		return rtrn;
	}


	private IntervalTree<Integer>[] makeTrees(Matrix sortedPermutations) {
		IntervalTree<Integer>[] trees=new IntervalTree[sortedPermutations.getColumnDimension()];
		
		for(int i=0; i<trees.length; i++){
			trees[i]=new IntervalTree<Integer>();
			double[] vals=sortedPermutations.getColumn(i);
			for(int j=0; j<vals.length; j++){
				double current=vals[j];
				trees[i].put(new Double(current*1000).intValue(), new Double(current*1000).intValue(), j);
			}
		}
		
		return trees;
	}


	//return the number based on the index of the observations
	private static double[] computeLargerForArray(double[] sortedObservations, double[] permutationValues){
		double[] rtrn=new double[sortedObservations.length];
		
		for(int i=0; i<sortedObservations.length; i++){
			double val=sortedObservations[i];
			for(int j=0; j<permutationValues.length; j++){
				double permVal=permutationValues[j];
				if(permVal>=val){rtrn[i]=(permutationValues.length-j); break;}
			}
		}
				
		return rtrn;
	}
	
	
	private double computeExactFDR(double k, double[] sortedObservations){
		double numberGreaterThanK=(double)countObservationsLargerThan(k, sortedObservations);
		double[] permutationGreaterThanK=new double[sortedPermutations.getColumnDimension()];
		
		for(int i=0; i<permutationGreaterThanK.length; i++){
			permutationGreaterThanK[i]=(double)countObservationsLargerThan(k, sortedPermutations.getColumn(i));
		}
		
		double expectedNumberGreaterThanK=Statistics.average(permutationGreaterThanK);
		
		double estimateOfMeanGreaterThanK=estimateMeanGreaterThanK(expectedNumberGreaterThanK, sortedObservations.length, numberGreaterThanK);
		//if(estimateOfMeanGreaterThanK==0){estimateOfMeanGreaterThanK=((double)1/sortedPermutations.getColumnDimension());}
		
		
		double fdr=estimateOfMeanGreaterThanK/numberGreaterThanK;
		//if(numberGreaterThanK==0){System.err.println(k+" "+(estimateOfMeanGreaterThanK/numberGreaterThanK)+" "+(expectedNumberGreaterThanK/numberGreaterThanK)+" "+expectedNumberGreaterThanK+" "+estimateOfMeanGreaterThanK+" "+numberGreaterThanK);}
		return Math.min(fdr, 1);
	}
	
		
	//This method implements the PaGE convergence estimate of the expectation of the mean
	private double estimateMeanGreaterThanK(double expectedNumberGreaterThanK, int numGenes, double Rk) {
		double mean=expectedNumberGreaterThanK;
		for(int i=0; i<6; i++){
			mean=(mean/numGenes)*(numGenes-(Rk-mean));
		}
		return mean;
	}

	/*private double computeExactNegativeFDR(double k){
		int numberGreaterThanK=countObservationsLessThan(k, sortedObservations);
		double[] permutationGreaterThanK=new double[sortedPermutations.getColumnDimension()];
		
		for(int i=0; i<permutationGreaterThanK.length; i++){
			permutationGreaterThanK[i]=countObservationsLessThan(k, sortedPermutations.getColumn(i));
		}
		
		double expectedNumberGreaterThanK=Statistics.average(permutationGreaterThanK);
		double estimateOfMeanGreaterThanK=estimateMeanGreaterThanK(expectedNumberGreaterThanK, sortedObservations.length, numberGreaterThanK);
		
		//return expectedNumberGreaterThanK/numberGreaterThanK;
		
		//System.err.println(k+" "+(estimateOfMeanGreaterThanK/numberGreaterThanK)+" "+(expectedNumberGreaterThanK/numberGreaterThanK)+" "+expectedNumberGreaterThanK+" "+estimateOfMeanGreaterThanK+" "+numberGreaterThanK);
		
		double fdr=expectedNumberGreaterThanK/numberGreaterThanK;
		return Math.min(fdr, 1);
	}*/
	
	/**
	 * Assumes sorted list!
	 * @param score
	 * @param list
	 * @return
	 */
	/*private static int countObservationsLessThan(double score, double[] vals) {
		int i;
		for(i = 0; i < vals.length; i++) {
			if(vals[i]>=score) {break;}
		}
		return i;
	}*/
	
	/**
	 * Assumes sorted list!
	 * @param score
	 * @param list
	 * @return
	 */
	private static int countObservationsLargerThan(double score, double[] vals) {
		int i;
		for(i = 0; i<vals.length; i++) {
			if(score <= vals[i]) {break;}
		}
		return vals.length-i;
	}
	

	
	private double[] normalizeScore(double[] observed, Matrix permutations) {
		double[] normObs =observed.clone();
		for (int i=0; i<observed.length; i++) {normObs[i]=normalize (observed[i], permutations.getRow(i)); }
		return normObs;
	}
	
	private double normalize(double score, double[] randomScores){
		double normScore=(score-Statistics.average(randomScores))/Statistics.stdev(randomScores);
		return normScore;
	}
	
	public static void main(String[] args){
		double[] rtrn={0,1,2,3,4,5,6,7,8,9,10};
		double[] vals={0,0,0,6,7,8,9,10,11,12,13};
		double[] greater=computeLargerForArray(rtrn, rtrn);
		
		for(int i=0; i<greater.length; i++){
			System.err.println(i+" "+rtrn[i]+" "+greater[i]);
		}
	
	}
	
}
