package broad.projection.gsa;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import Jama.Matrix;
import broad.core.datastructures.MatrixWithHeaders;
import broad.core.math.EmpiricalDistribution;
import broad.core.math.FDRDistribution;
import broad.core.math.KSTest;
import broad.core.math.Statistics;
import broad.pda.differentialExpression.DifferentialScoring;
import broad.pda.differentialExpression.Permutations;

public class GeneSetEnrichment {

	int numPerm=100; //TODO Must change!!
	double[] fudgeFactors={0};
	double alpha=.05;
	
	MatrixWithHeaders rankedList;
	MatrixWithHeaders normalizedKSEnrichments;
	MatrixWithHeaders normalizedMaxMeanEnrichments;
	MatrixWithHeaders ksFDR;
	MatrixWithHeaders maxMeanFDR;
	MatrixWithHeaders KSEnrichments;
	MatrixWithHeaders maxMeanEnrichments;
	private boolean useFold=false;
	private boolean paired=false;
	private boolean normalizeScores=true;
	private boolean permutePhenotype=true;
	
	public GeneSetEnrichment(MatrixWithHeaders data, Collection<String> group1, Collection<String> group2, Map<String, Collection<String>> geneSets, boolean normalizeScores, boolean permutePhenotype){
		this.normalizeScores=normalizeScores;
		this.permutePhenotype=permutePhenotype;
		
		//Step 1: Make a ranked list based on the difference between classes
		rankedList=DifferentialScoring.computeTestStatistics(data, group1, group2, fudgeFactors, useFold, paired);
		//TODO Doesnt the ranked list need to be sorted??
		System.err.println("Completed ranking the list");
		
		
		//Step 2: Determine the enrichment score for the genes in the gene set
		KSEnrichments=scoreKSEnrichment(rankedList, geneSets);
		System.err.println("Completed KS Scoring");
		
		maxMeanEnrichments=scoreMaxMeanEnrichment(rankedList, geneSets);
		System.err.println("Completed max mean scoring");
		
		if(this.normalizeScores){
			//Step 3: Normalize enrichment score
			//TODO Assign FDR values based on the gene set permutations
			MatrixWithHeaders[] norm=this.normalizeEnrichmentScores(rankedList, geneSets, numPerm, KSEnrichments, maxMeanEnrichments);
			normalizedKSEnrichments=norm[0];
			normalizedMaxMeanEnrichments=norm[1];
			ksFDR=norm[2];
			maxMeanFDR=norm[3];
			System.err.println("Completed normalization");	
		}
	
		/*if(permutePhenotype){
			//Step 4: Permute the data and compute an FDR for the enrichment scores
			MatrixWithHeaders[] fdrs=computeFDRForGeneSet(data, rankedList, group1, group2, geneSets, fudgeFactors, KSEnrichments, maxMeanEnrichments, alpha);
			ksFDR=fdrs[0];
			maxMeanFDR=fdrs[1];
		}*/	
		
	}
	
	private MatrixWithHeaders[] normalizeEnrichmentScores(MatrixWithHeaders rankedList2,Map<String, Collection<String>> geneSets, int numPerm,	MatrixWithHeaders KS, MatrixWithHeaders maxMeanEnrichments) {
		MatrixWithHeaders rtrnKS=new MatrixWithHeaders(KS.getRowNames(), KS.getColumnNames());
		MatrixWithHeaders rtrnMM=new MatrixWithHeaders(KS.getRowNames(), KS.getColumnNames());
		
		
		int counter=0;
		System.err.print("Scoring genes "+geneSets.size()+"....");
		Matrix permsKS=new Matrix(geneSets.size(), numPerm);
		Matrix permsMM=new Matrix(geneSets.size(), numPerm);
		
		int geneSetIndex=0;
		for(String geneSet: geneSets.keySet()){
			System.err.print(" "+counter);
			for(int i=0; i<KS.columnDimension(); i++){
				double[][] perms=permutations(rankedList, numPerm, geneSets.get(geneSet));
				double[] norm=this.normalize(rankedList, geneSets.get(geneSet), KS.get(geneSet, i), maxMeanEnrichments.get(geneSet, i), numPerm, perms);
				rtrnKS.set(geneSet, i, norm[0]);
				rtrnMM.set(geneSet, i, norm[1]);
				
				permsKS.setRow(geneSetIndex, perms[2]);
				permsMM.setRow(geneSetIndex, perms[3]);
			}
			counter++;
			geneSetIndex++;
		}
		System.err.println();
		
		
		MatrixWithHeaders fdrKS=assignFDR(rtrnKS, permsKS);
		MatrixWithHeaders fdrMM=assignFDR(rtrnMM, permsMM);
				
		MatrixWithHeaders[] rtrn={rtrnKS, rtrnMM, fdrKS, fdrMM};
		return rtrn;
	}

	private MatrixWithHeaders assignFDR(MatrixWithHeaders rtrnKS, Matrix permsKS) {
		MatrixWithHeaders rtrn=new MatrixWithHeaders(rtrnKS.getRowNames(), rtrnKS.getColumnNames());
		FDRDistribution distKS=new FDRDistribution(rtrnKS.getColumn(0), permsKS, alpha);
		
		for(String gene: rtrnKS.getRowNames()){
			for(String experiment: rtrnKS.getColumnNames()){
				double val=rtrnKS.get(gene, experiment);
				double fdr=Math.min(distKS.getFDR(val), distKS.getAbsFDR(val));
				rtrn.set(gene, experiment, fdr);
			}
		}
		return rtrn;
	}

	private double[][] permutations(MatrixWithHeaders rankedList, int numPermutations, Collection<String> geneSet) {
		double[] permutationsKS=new double[numPermutations];
		double[] permutationsMM=new double[numPermutations];
		
		double[] normKS=new double[numPermutations];
		double[] normMM=new double[numPermutations];
		
		for(int i=0; i<numPermutations; i++){
			//System.err.print(" "+i);
			Collection<String> randomSet=permute(geneSet, rankedList.getRowNames());
			permutationsKS[i]=KSTest.KSScores(rankedList, randomSet, 0)[0];
			permutationsMM[i]=maxMeanGeneSet(rankedList, randomSet, 0);
		}
		
		for(int i=0; i<normKS.length; i++){
			double KSNorm=normalize(permutationsKS[i], permutationsKS);
			normKS[i]=KSNorm;
		}
		
		for(int i=0; i<normMM.length; i++){
			double MMNorm=normalize(permutationsMM[i], permutationsMM);
			normMM[i]=MMNorm;
		}
		
		double[][] rtrn={permutationsKS, permutationsMM, normKS, normMM};
		return rtrn;
	}

	public GeneSetEnrichment(MatrixWithHeaders data, Collection<String> group1, Collection<String> group2, Map<String, Collection<String>> geneSets){
		this(data, group1, group2, geneSets, true, true);
	}
	
	//Supports continuous data. i.e: ranking genes according to their correlation with a specific pattern
	//patternData should be a one line matrix
	//Time complexity: n+ 2mn +2*m*n*numPerm*k+  (k=average gene set size, m=num of gene sets, n=num of probes)
	public GeneSetEnrichment(MatrixWithHeaders patternData,MatrixWithHeaders geneData , Map<String, Collection<String>> geneSets, boolean normalizeScores, boolean permutePhenotype){
		this.normalizeScores=normalizeScores;
		this.permutePhenotype=permutePhenotype;
		
		//Step 1: Make a ranked list based on correlation between the gene data and the pattern data
		rankedList=computeContinuousTestStatistics(patternData, geneData);
		System.err.println("Completed ranking the list");
		
		//Step 2: Determine the enrichment score for the genes in the gene set
		KSEnrichments=scoreKSEnrichment(rankedList, geneSets);
		System.err.println("Completed KS Scoring");
		
		maxMeanEnrichments=scoreMaxMeanEnrichment(rankedList, geneSets);
		System.err.println("Completed max mean scoring");
		
		if(this.normalizeScores){
			//Step 3: Normalize enrichment score
			normalizedKSEnrichments=this.normalizeKSEnrichmentScore(rankedList, geneSets, numPerm, KSEnrichments);
			System.err.println("Completed KS normalization");
			
			normalizedMaxMeanEnrichments=this.normalizeMaxMeanEnrichmentScore(rankedList, geneSets, numPerm, maxMeanEnrichments);
			System.err.println("Completed max mean normalization");		
		}
		
		if(permutePhenotype){
			//TODO : support permutation of phenotype
		//Step 4: Permute the data and compute an FDR for the enrichment scores
		//MatrixWithHeaders[] fdrs=computeFDRForGeneSet(data, rankedList, group1, group2, geneSets, fudgeFactors, KSEnrichments, maxMeanEnrichments, alpha);
		//ksFDR=fdrs[0];
		//maxMeanFDR=fdrs[1];
		}
		else // If we are not using phenotype permutation we calculate FDR based on the 
			//normalized KS and MaxMean enrichments
		{
			MatrixWithHeaders[] fdrs=this.computeFDRForGeneSet_permuteGene(rankedList, geneSets, KSEnrichments, maxMeanEnrichments,alpha);
			ksFDR=fdrs[0];
			maxMeanFDR=fdrs[1];
		}
	}
	
	

	public MatrixWithHeaders getRankedList(){return this.rankedList;}
	public MatrixWithHeaders getKSFDR(){return this.ksFDR;}
	public MatrixWithHeaders getMaxMeanFDR(){return this.maxMeanFDR;}
	public MatrixWithHeaders getMaxMeanEnrichments(){return this.maxMeanEnrichments;}
	public MatrixWithHeaders getKSEnrichments(){return this.KSEnrichments;}
	public MatrixWithHeaders getNormalizedKSEnrichments(){return this.normalizedKSEnrichments;}
	public MatrixWithHeaders getNormalizedMaxMeanEnrichments(){return this.normalizedMaxMeanEnrichments;}

	private MatrixWithHeaders[] computeFDRForGeneSet(MatrixWithHeaders data,MatrixWithHeaders rankedList, Collection<String> group1, Collection<String> group2, Map<String, Collection<String>> geneSets, double[] fudgeFactors, MatrixWithHeaders kSEnrichments, MatrixWithHeaders maxMeanEnrichments, double alpha) {
		MatrixWithHeaders[] ks=new MatrixWithHeaders[numPerm];
		MatrixWithHeaders[] maxMean=new MatrixWithHeaders[numPerm];
		
		//Step 1: Generate permuted ranked lists
		MatrixWithHeaders[] permutedRankedLists=Permutations.permuteData(data, group1, group2, numPerm, fudgeFactors, useFold, paired);
		
		//Step 2: Score each permuted ranked list by gene set
		for(int i=0; i<permutedRankedLists.length; i++){
			System.err.println("Permutation "+i);
			ks[i]=this.scoreKSEnrichment(permutedRankedLists[i], geneSets);
			maxMean[i]=this.scoreMaxMeanEnrichment(permutedRankedLists[i], geneSets);
		}
		
		//Step 3: Adding FDR value
		MatrixWithHeaders ksFDR=assignFDRValues(kSEnrichments, ks, alpha, false); //MITCH: if you want FDR to control for gene set size , change the last value to be true
		MatrixWithHeaders maxMeanFDR=assignFDRValues(maxMeanEnrichments, maxMean, alpha,false);//MITCH: if you want FDR to control for gene set size , change the last value to be true
		
		MatrixWithHeaders[] rtrn={ksFDR, maxMeanFDR};
		return rtrn;
	}

	
	//Compute FDR based on gene permutations rather than phenotype permutation
	
	private MatrixWithHeaders[] computeFDRForGeneSet_permuteGene(MatrixWithHeaders rankedList,  Map<String, Collection<String>> geneSets,  MatrixWithHeaders kSEnrichments, MatrixWithHeaders maxMeanEnrichments,double alpha) {
		MatrixWithHeaders[] ks=new MatrixWithHeaders[numPerm];
		MatrixWithHeaders[] maxMean=new MatrixWithHeaders[numPerm];
		
		//Step 1: Generate permuted ranked lists: permute the list of Pearson correlations
		MatrixWithHeaders[] permutedRankedLists=new MatrixWithHeaders[numPerm];
		for (int i=0; i < numPerm; i++){
			//randomly permute the values in each col of rankList
			MatrixWithHeaders tmp=rankedList.copy();
			tmp.randomPermuteColumns();
			permutedRankedLists[i]=tmp;
		}
		//Step 2: Score each permuted ranked list by gene set
		for(int i=0; i<permutedRankedLists.length; i++){
			System.err.println("Permutation "+i);
			ks[i]=this.scoreKSEnrichment(permutedRankedLists[i], geneSets);
			maxMean[i]=this.scoreMaxMeanEnrichment(permutedRankedLists[i], geneSets);
		}
		
		//Step 3: Adding FDR value
		MatrixWithHeaders ksFDR=assignFDRValues(kSEnrichments, ks, alpha,true);
		MatrixWithHeaders maxMeanFDR=assignFDRValues(maxMeanEnrichments, maxMean, alpha, true);
		
		MatrixWithHeaders[] rtrn={ksFDR, maxMeanFDR};
		return rtrn;
	}

	
	private MatrixWithHeaders assignFDRValues(MatrixWithHeaders data, MatrixWithHeaders[] permutations, double alpha, boolean normalizeScr) {
		MatrixWithHeaders fdr=new MatrixWithHeaders(data.getRowNames(), data.getColumnNames()); 
		
		for(String column: data.getColumnNames()){
			Matrix matrix=new Matrix(data.rowDimension(), permutations.length);
			for(int i=0; i<permutations.length; i++){matrix.setColumn(i, permutations[i].getColumn(column));}
			double[] observed=data.getColumn(column);
			FDRDistribution dist=new FDRDistribution(observed, matrix, alpha);
			for(String row: data.getRowNames()){
				double fdrVal=dist.getFDR(data.get(row, column));
				fdr.set(row, column, fdrVal);
			}
		}
		
		return fdr;
	}

	private MatrixWithHeaders scoreKSEnrichment(MatrixWithHeaders rankedList, Map<String, Collection<String>> geneSets) {
		MatrixWithHeaders KS=new MatrixWithHeaders(new ArrayList(geneSets.keySet()), rankedList.getColumnNames());
				
		for(String geneSet: geneSets.keySet()){
			//System.err.println("GENE SET SIZE: "+geneSet+"\t"+geneSets.get(geneSet).size());
			//Step 1: score the gene set 
			for(int i=0; i<rankedList.columnDimension(); i++){
				//TODO: Convert the coordinates if needed
				//TODO We should do this conversion only once
				//Collection<String> geneList=convertIDs(rankedList.getPIDToName(), geneSets.get(geneSet));
				Collection<String> geneList=geneSets.get(geneSet);
				KS.set(geneSet, i, KSTest.KSScores(rankedList, geneList, i)[0]); //This is not normalized.
			}
		}
		
		return KS;
	}
	
	//The gene list is either PIDs or gene names but we want o make it the same as the UID as the ranked list
	/*private Collection<String> convertIDs(Map<String, String> map, Collection<String> geneSet) {
		Collection<String> rtrn=new TreeSet<String>();
		//check if any of the genes arein the PIDs
		int pidCounter=0;
		for(String gene: geneSet){
			if(map.containsKey(gene)){pidCounter++;}
		}
		if(pidCounter>0){return geneSet;}
		
		for(String pid: map.keySet()){
			String name=map.get(pid);
			if(geneSet.contains(name)){rtrn.add(pid);}
		}
		return rtrn;
	}*/

	private MatrixWithHeaders scoreMaxMeanEnrichment(MatrixWithHeaders rankedList, Map<String, Collection<String>> geneSets) {
		MatrixWithHeaders KS=new MatrixWithHeaders(new ArrayList(geneSets.keySet()), rankedList.getColumnNames());
				
		for(String geneSet: geneSets.keySet()){
			//System.err.println(geneSet+"\t"+geneSets.get(geneSet).size());
			//Step 1: score the gene set 
			for(int i=0; i<rankedList.columnDimension(); i++){
				KS.set(geneSet, i, maxMeanGeneSet(rankedList, geneSets.get(geneSet), i));
			}
		}
		
		return KS;
	}
	
	private MatrixWithHeaders normalizeKSEnrichmentScore(MatrixWithHeaders rankedList, Map<String, Collection<String>> geneSets, int numPerm, MatrixWithHeaders KS){
		MatrixWithHeaders rtrn=new MatrixWithHeaders(KS.getRowNames(), KS.getColumnNames());
		
		for(String geneSet: geneSets.keySet()){
			//System.err.println(geneSet+"\t"+geneSets.get(geneSet).size());
			for(int i=0; i<KS.columnDimension(); i++){
				double norm=this.normalizeKS(rankedList, geneSets.get(geneSet), KS.get(geneSet, i), numPerm);
				rtrn.set(geneSet, i, norm);
			}
		}
		
		return rtrn;
	}
	
	private MatrixWithHeaders normalizeMaxMeanEnrichmentScore(MatrixWithHeaders rankedList, Map<String, Collection<String>> geneSets, int numPerm, MatrixWithHeaders maxMean){
		MatrixWithHeaders rtrn=new MatrixWithHeaders(maxMean.getRowNames(), maxMean.getColumnNames());
		
		for(String geneSet: geneSets.keySet()){
			//System.err.println(geneSet+"\t"+geneSets.get(geneSet).size());
			for(int i=0; i<maxMean.columnDimension(); i++){
				double norm=this.normalizeMaxMean(rankedList, geneSets.get(geneSet), maxMean.get(geneSet, i), numPerm);
				rtrn.set(geneSet, i, norm);
			}
		}
		
		return rtrn;
	}
	
	//When calculating permutation of gene sets, we are just storing the permutation for ONE gene set 
	//and not all gene sets. 
	private double normalizeKS(MatrixWithHeaders rankedList,Collection<String> geneSet, double KS, int numPermutations) {
		//Generate random permutations of the gene set
		double[] permutations=new double[numPermutations];
		for(int i=0; i<numPermutations; i++){
			Collection<String> randomSet=permute(geneSet, rankedList.getRowNames());
			permutations[i]=KSTest.KSScores(rankedList, randomSet, 0)[0];
		}
		return normalize(KS, permutations);
	}
	
	//When calculating permutation of gene sets, we are just storing the permutation for ONE gene set 
	//and not all gene sets. 
	private double[] normalize(MatrixWithHeaders rankedList,Collection<String> geneSet, double KS, double maxMean, int numPermutations, double[][] perms) {
		//Generate random permutations of the gene set
		double[] permutationsKS=new double[numPermutations];
		double[] permutationsMM=new double[numPermutations];
		if(perms==null){
			//System.err.print("Permutations");
			for(int i=0; i<numPermutations; i++){
				//System.err.print(" "+i);
				Collection<String> randomSet=permute(geneSet, rankedList.getRowNames());
				permutationsKS[i]=KSTest.KSScores(rankedList, randomSet, 0)[0];
				permutationsMM[i]=maxMeanGeneSet(rankedList, randomSet, 0);
			}
		}
		else{
			permutationsKS=perms[0];
			permutationsMM=perms[1];
		}
		
		//System.err.println("\n");
		double KSNorm=normalize(KS, permutationsKS);
		double MMNorm=normalize(maxMean, permutationsMM);
		
		double KSFDR=pvalue(KS, permutationsKS);
		double MMFDR=pvalue(maxMean, permutationsMM);
		
		double[] rtrn={KSNorm, MMNorm, KSFDR, MMFDR};
		return rtrn;
	}
	
	private double pvalue(double kS, double[] permutationsKS) {
		EmpiricalDistribution dist=new EmpiricalDistribution(permutationsKS, 200);
		return dist.getCumulativeProbability(kS);
	}

	private double normalizeMaxMean(MatrixWithHeaders rankedList,Collection<String> geneSet, double maxMean, int numPermutations) {
		//Generate random permutations of the gene set
		double[] permutations=new double[numPermutations];
		for(int i=0; i<numPermutations; i++){
			Collection<String> randomSet=permute(geneSet, rankedList.getRowNames());
			permutations[i]=this.maxMeanGeneSet(rankedList, randomSet, 0);
		}
		return normalize(maxMean, permutations);
	}
	
	
	private double normalize(double score, double[] randomScores){
		double normScore=(score-Statistics.average(randomScores))/Statistics.stdev(randomScores);
		return normScore;
	}

	private Collection<String> permute(Collection<String> geneSet, List<String> genes){
		Collection<String> rtrn=new HashSet<String>();
						
		for(int i=0; i<geneSet.size(); i++){
			int index=new Double(Math.random()*genes.size()).intValue();
			String randomGene=(String)genes.remove(index);
			rtrn.add(randomGene);
		}
		
		return rtrn;
	}
	
	private double maxMeanGeneSet(MatrixWithHeaders rankedList, Collection<String> geneSet, int index){
		MatrixWithHeaders matrix=rankedList.submatrixByRowNames(geneSet);
		//genes are the rows and columns are samples
		//want to compute a statisitc on the rows by columns
		
		double[] rtrn=new double[matrix.getNumberColumns()];
		
		for(int i=0; i<rtrn.length; i++){
			double[] vals=matrix.getColumn(i);
			//simple statistic for testing
			rtrn[i]=Statistics.maxmean(vals);
		}
		
		return rtrn[index]; 
	}
	
	
	// Calculate the correlation between the geneData and the patternData and returns a rank list
	private MatrixWithHeaders computeContinuousTestStatistics(
			MatrixWithHeaders patternData, MatrixWithHeaders geneData) {
		
		List<String> columns=new ArrayList<String>();
		columns.add("Pearson");
		MatrixWithHeaders statisticMatrix=new MatrixWithHeaders(geneData.getRowNames(),columns );
		
		double[] gr1=patternData.getRow(patternData.getRowNames().get(0));
		
		for(String gene: geneData.getRowNames()){
			double[] gr2=geneData.getRow(gene);
			statisticMatrix.set(gene,"Pearson",Statistics.pearsonDistance(gr1,gr2));
		}
		
		return statisticMatrix;
	}
}
