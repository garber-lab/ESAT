package broad.pda.differentialExpression;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.log4j.Logger;

import Jama.Matrix;
import broad.core.datastructures.MatrixWithHeaders;
import broad.core.math.CombinationGenerator;
import broad.core.math.FDRDistribution;


public class Permutations {
	static Logger logger = Logger.getLogger(Permutations.class.getName());
	//TODO This should call permuteData for multiple groups
	public static MatrixWithHeaders[] permuteData(MatrixWithHeaders data, Collection<String> group1, Collection<String> group2, int numberPermutations, double[] fudgeFactors, boolean useFold, boolean paired) {
		MatrixWithHeaders[] rtrn;
		CombinationGenerator comb=new CombinationGenerator((group1.size()+group2.size()), group1.size());
		ArrayList<String> samples=new ArrayList<String>();
		samples.addAll(group1);
		samples.addAll(group2);
		
		int totalNumberPossible=comb.getTotal();
		if(totalNumberPossible<0){totalNumberPossible=1000; }//to account for overflow
		logger.debug("Size of group1: "+group1.size()+" group2: "+group2.size()+" number of possible permutations: "+ totalNumberPossible);
		
		//TODO: Fix this once its working
		if(totalNumberPossible<numberPermutations){
			logger.debug("Doing all "+totalNumberPossible +" permutations");
			//do all
			rtrn=new MatrixWithHeaders[totalNumberPossible];
			int i=0;
			while(comb.hasMore()){
				System.err.print(" "+i);
				int[] next=comb.getNext();
				rtrn[i]=generatePermutation(data, next, samples, fudgeFactors, useFold, paired);
				i++;
			}
			System.err.println();
		}
		
		else{
			logger.debug("Doing a random "+numberPermutations+" out of "+totalNumberPossible +" possible permutations");
			//do numPermutations randomly
			rtrn=new MatrixWithHeaders[numberPermutations];
			for(int i=0; i<numberPermutations; i++){
				System.err.print(" "+i);
				int[] next=comb.getNextRandom();
				rtrn[i]=generatePermutation(data, next, samples, fudgeFactors, useFold, paired);
			}
			System.err.println();
		}
		
		return rtrn;
	}
	
	private static MatrixWithHeaders generatePermutation(MatrixWithHeaders data, int[] permutedIndeces, ArrayList<String> samples, double[] fudgeFactors, boolean useFold, boolean paired) {
		Collection<String> group1=new TreeSet<String>();
		Collection<String> group2=new TreeSet<String>();
		for(int i=0; i<permutedIndeces.length; i++){
			//System.err.println(permutedIndeces[i]);
			group1.add(samples.get(permutedIndeces[i]));
		}
		for(String sample: samples){if(!group1.contains(sample)){group2.add(sample);}}
		logger.debug("permuted group 1: " + group1 + " perm g2: " + group2);
		return DifferentialScoring.computeTestStatistics(data, group1, group2, fudgeFactors, useFold, paired);
		
	}

	//MG: This is a new generic constructor that supports as many groups as provided
	public static MatrixWithHeaders[] permuteData(MatrixWithHeaders data, Map<String, Collection<String>> groups, int numberPermutations, double[] fudgeFactors, boolean useFold, boolean paired) {
		if(groups.size()==2){
			Iterator<Collection<String>> iter=groups.values().iterator();
			Collection<String> group1=iter.next();
			Collection<String> group2=iter.next();
			return permuteData(data, group1, group2, numberPermutations, fudgeFactors, useFold, paired);
		}
		
		MatrixWithHeaders[] rtrn=new MatrixWithHeaders[numberPermutations];;
		
		//TODO How do we do this permutations systematically
		ArrayList<String> samples=new ArrayList<String>();
		for(String group: groups.keySet()){samples.addAll(groups.get(group));}
		
		
		for(int i=0; i<numberPermutations; i++){
			//Generate permuted groups
			rtrn[i]=permutedGroups(data, fudgeFactors, useFold, paired, groups, samples);
		}
				
		return rtrn;
	}

	private static MatrixWithHeaders permutedGroups(MatrixWithHeaders data, double[] fudgeFactors, boolean useFold, boolean paired, Map<String, Collection<String>> groups, ArrayList<String> samples) {
		Map<String, Collection<String>> permutedGroups=new TreeMap<String, Collection<String>>();
				
		ArrayList<String> remainingSamples=new ArrayList<String>();
		remainingSamples.addAll(samples);
		
		for(String group: groups.keySet()){
			Collection<String> list=new TreeSet<String>();
			int size=groups.get(group).size();
			for(int i=0; i<size; i++){
				int index =new Double(Math.random()*remainingSamples.size()).intValue();
				String sample=remainingSamples.remove(index);
				list.add(sample);
			}
			permutedGroups.put(group, list);
		}		
		
		return DifferentialScoring.computeTestStatistics(data, permutedGroups, fudgeFactors, useFold, paired);
	}

	
	//This method generate permutations of groups from all the data but keeping the sizes of samples and controls the same
	//This is way too conservative
	//Keep all controls and a random sample of other samples from the remainder
	public static MatrixWithHeaders[] permuteDataByGroupSize(MatrixWithHeaders data, int groupSize, Map<String, Collection<String>> groups, int numPerm, double[] fudgeFactors, boolean useFold) {
		Collection<String> controls=groups.get("Control");
		
		//System.err.println(controls);
		
		Collection<String> remainder=new TreeSet<String>();
		
		for(String group: groups.keySet()){
			if(!group.equalsIgnoreCase("Control")){
				remainder.addAll(groups.get(group));
			}
		}
		
		MatrixWithHeaders[] rtrn=new MatrixWithHeaders[numPerm];
		
		//TODO Make part of combination generator
		for(int i=0; i<numPerm; i++){
			System.err.print(" "+i);
			//First pick the samples
			ArrayList<String> all=new ArrayList<String>();
			all.addAll(remainder);
			
			Collection<String> samples=new TreeSet();
			
			for(int j=0; j<groupSize; j++){
				int random=new Double(Math.random()*all.size()).intValue();
				String sample=all.remove(random);
				samples.add(sample);
			}
			
			ArrayList<String> temp=new ArrayList();
			temp.addAll(samples);
			temp.addAll(controls);
			
			//pick random
			CombinationGenerator gen=new CombinationGenerator(samples.size()+controls.size(), samples.size());
			int[] next=gen.getNextRandom();
			rtrn[i]=generatePermutation(data, next, temp, fudgeFactors, useFold, false);
		}
		System.err.println();
		
		return rtrn;
	}

	public static MatrixWithHeaders[] permuteAndAssignFDR(MatrixWithHeaders testStats, MatrixWithHeaders data, Collection<String> group1,Collection<String> group2, int numberPermutations,double[] fudgeFactors, boolean useFold, boolean paired, double alpha) {
		MatrixWithHeaders fdr=new MatrixWithHeaders(testStats.getRowNames(), testStats.getColumnNames()); 
		MatrixWithHeaders absFDR=new MatrixWithHeaders(testStats.getRowNames(), testStats.getColumnNames()); 
		
		
		//MatrixWithHeaders[] permMatrix; //Might not want to save this
		CombinationGenerator comb=new CombinationGenerator((group1.size()+group2.size()), group1.size());
		ArrayList<String> samples=new ArrayList<String>();
		samples.addAll(group1);
		samples.addAll(group2);
		
		int totalNumberPossible=comb.getTotal();
		if(totalNumberPossible<0){totalNumberPossible=1000; }//to account for overflow
		System.err.println("Size of group1: "+group1.size()+" group2: "+group2.size()+" number of possible permutations: "+ totalNumberPossible);
		
		Map<String, Matrix> map=new TreeMap<String, Matrix>();
		for(String column: testStats.getColumnNames()){
			//make new matrix consisting of all genes by all permutations
			Matrix matrix=new Matrix(testStats.rowDimension(), numberPermutations);
			map.put(column, matrix);
		}
		
		if(totalNumberPossible<numberPermutations){
			System.err.println("Doing all "+totalNumberPossible +" permutations");
			//do all
			//permMatrix=new MatrixWithHeaders[totalNumberPossible];
			int i=0;
			while(comb.hasMore()){
				//for each permutation
				System.err.print(" "+i);
				int[] next=comb.getNext();
				MatrixWithHeaders perm=generatePermutation(data, next, samples, fudgeFactors, useFold, paired); //test stat scores for perm
				for(String column: perm.getColumnNames()){
					Matrix matrix=map.get(column);
					matrix.setColumn(i, perm.getColumn(column));
					map.put(column, matrix);
				}
				i++;
			}
			System.err.println();
		}
		
		else{
			System.err.println("Doing a random "+numberPermutations+" out of "+totalNumberPossible +" possible permutations");
			//do numPermutations randomly
			//permMatrix=new MatrixWithHeaders[numberPermutations];
			for(int i=0; i<numberPermutations; i++){
				System.err.print(" "+i);
				int[] next=comb.getNextRandom();
				MatrixWithHeaders perm=generatePermutation(data, next, samples, fudgeFactors, useFold, paired); //test stat scores for perm
				for(String column: perm.getColumnNames()){
					Matrix matrix=map.get(column);
					matrix.setColumn(i, perm.getColumn(column));
					map.put(column, matrix);
				}				
			}
			System.err.println();
		}
		
		
		for(String column: testStats.getColumnNames()){
			double[] observed=testStats.getColumn(column);
			Matrix matrix=map.get(column);
			FDRDistribution dist=new FDRDistribution(observed, matrix, alpha);
			for(String row: testStats.getRowNames()){
				//get test-stat value
				double val=testStats.get(row, column);
				//get FDR value to gene
				double fdrVal=dist.getFDR(val); //TODO: Consider whether to do this or not
				//get absolute FDR value for gene
				double absFDRVal=dist.getAbsFDR(val);
				absFDR.set(row, column, absFDRVal);
				fdr.set(row, column, fdrVal);
			}
		}
		
		
		MatrixWithHeaders[] rtrn={fdr, absFDR};
		
		return rtrn;

		
	}
	
}
