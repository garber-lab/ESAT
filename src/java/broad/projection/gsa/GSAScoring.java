package broad.projection.gsa;

import java.io.*;
import java.util.*;

import broad.core.datastructures.MatrixWithHeaders;
import broad.core.math.KSTest;
import broad.core.math.Statistics;
import broad.core.util.GMTParser;
import broad.core.util.ParseGCTFile;
import broad.pda.geneexpression.PaGEJava;
import broad.pda.geneexpression.ScoreDiffGenes;

public class GSAScoring {
	
	int minNumber = 0;//TODO: check what the minNUmber should be
	int numPerms=10;
	double fudgeFactor=.1; 
	MatrixWithHeaders[] restandardizedScores;
	double alpha=1;
	double minFold=2.0;

	//allow for fixing control class rather than all enumeration
	public GSAScoring(File expressionFile, File clsFile, File geneSetFile, String save)throws Exception{
		this.restandardizedScores=this.computeEnrichmentScores(expressionFile, clsFile, geneSetFile);
		this.restandardizedScores[0].writeGCT(save+".KS.gct");
		this.restandardizedScores[1].writeGCT(save+".MaxMean.gct");
	}
	
	public GSAScoring(File expressionFile, File clsFile, String save)throws Exception{
		this.restandardizedScores=this.computeEnrichmentScores(expressionFile, clsFile);
		this.restandardizedScores[0].writeGCT(save+".KS.gct");
		this.restandardizedScores[1].writeGCT(save+".MaxMean.gct");
	}
	
	public MatrixWithHeaders getKSScores(){return this.restandardizedScores[0];}
	public MatrixWithHeaders getMaxMeanScores(){return this.restandardizedScores[1];}
	
	
	private MatrixWithHeaders[] computeEnrichmentScores(File expressionFile, File clsFile, File geneSetFile)throws Exception{
		MatrixWithHeaders expression=new MatrixWithHeaders(expressionFile.getAbsolutePath());
		Map<String, ArrayList> groupIndexes=ParseGCTFile.getGroupIndexes(clsFile);
		Map<String, Collection<String>> geneSets=GMTParser.ParseGMTFile(geneSetFile, minNumber);
		
		//TODO come up with convention for controlClass indication
		String controlClass=ParseGCTFile.getControlClass(clsFile);
		
		return this.computeEnrichmentScores(expression, groupIndexes, geneSets, null);
	}
	
	private MatrixWithHeaders[] computeEnrichmentScores(File expressionFile, File clsFile)throws Exception{
		MatrixWithHeaders expression=new MatrixWithHeaders(expressionFile.getAbsolutePath());
		Map<String, ArrayList> groupIndexes=ParseGCTFile.getGroupIndexes(clsFile);
		//TODO come up with convention for controlClass indication
		String controlClass=ParseGCTFile.getControlClass(clsFile);
		
		PaGEJava page=new PaGEJava(expression, groupIndexes, controlClass);
		Map<String, Collection<String>> geneSets=page.getSignificantGeneSets(alpha, minFold);
		
		
		return this.computeEnrichmentScores(expression, groupIndexes, geneSets, controlClass);
	}
	
	
	private MatrixWithHeaders[] computeEnrichmentScores(MatrixWithHeaders expression, Map<String, ArrayList> groupIndexes, Map<String, Collection<String>> geneSets, String controlClass){
		//Score each list by t-statistic
		MatrixWithHeaders scoreList=ScoreDiffGenes.scoreListByTStatistic(expression, groupIndexes, controlClass, fudgeFactor);
		
		//KSScore each geneSet
		MatrixWithHeaders KSScores=this.KSScoreGeneSets(scoreList, geneSets);
		
		//maxmean matrix by genesets
		MatrixWithHeaders maxMeanGeneSets=this.scoreGeneSets(scoreList, geneSets);
		
		//restandardize scores by permutations
		//KS, MM
		
		MatrixWithHeaders[] restandardized=this.restandardize(scoreList, geneSets, expression.getRowNames(), KSScores, maxMeanGeneSets);
		
		return restandardized;
		
	} 
	
		
	//normalize by expectation and variance of permuted geneSets
	private MatrixWithHeaders[] restandardize(MatrixWithHeaders scoreList, Map<String, Collection<String>> geneSets, List<String> genes, MatrixWithHeaders KSScores, MatrixWithHeaders maxMeanScores){
		
		MatrixWithHeaders[] randomizationsKS=new MatrixWithHeaders[numPerms];
		MatrixWithHeaders[] randomizationsMM=new MatrixWithHeaders[numPerms];
		
		for(int i=0; i<numPerms; i++){
			//permute geneSet
			Map<String, Collection<String>> permutedGeneSet=permuteGeneSets(geneSets, genes);
			
			//KSScore each geneSet
			randomizationsKS[i]=this.KSScoreGeneSets(scoreList, permutedGeneSet);
				
			//maxmean matrix by genesets
			randomizationsMM[i]=this.scoreGeneSets(scoreList, permutedGeneSet);
		}
		
		MatrixWithHeaders[] rtrn={restandardize(KSScores, randomizationsKS),	restandardize(maxMeanScores, randomizationsMM)};
		
		
		return rtrn;
	}
	
	private MatrixWithHeaders restandardize(MatrixWithHeaders scores, MatrixWithHeaders[] randomizations){
		MatrixWithHeaders rtrn=new MatrixWithHeaders(scores.getRowNames(), scores.getColumnNames());
		
		for(int i=0; i<scores.rowDimension(); i++){
			for(int j=0; j<scores.columnDimension(); j++){
				double score=(scores.get(i, j)-average(randomizations,i,j))/stdev(randomizations, i,j);
				rtrn.set(i, j, score);
				System.err.println(scores.getRowNames().get(i)+" "+scores.getColumnNames().get(j)+" "+score);
			}
		}
		return rtrn;
	}
	
	private double average(MatrixWithHeaders[] random, int i, int j){
		double[] vals=new double[random.length];
		for(int k=0; k<vals.length; k++){vals[k]=random[k].get(i, j);}
		return Statistics.average(vals);
	}
	
	private double stdev(MatrixWithHeaders[] random, int i, int j){
		double[] vals=new double[random.length];
		for(int k=0; k<vals.length; k++){vals[k]=random[k].get(i, j);}
		double stdev=Statistics.stdev(vals);
		System.err.println(Statistics.average(vals)+" "+stdev);
		if(Double.isNaN(stdev)){for(int k=0; k<vals.length; k++){System.err.println(vals[k]);}}
		
		return stdev;
	}
	
	private Map<String, Collection<String>> permuteGeneSets(Map<String, Collection<String>> geneSets, List<String> genes){
		Map rtrn=new TreeMap();
	
		for(String name: geneSets.keySet()){
			Collection<String> permuted=permute(geneSets.get(name).size(), genes);
			rtrn.put(name, permuted);
		}
		
		return rtrn;
	}
	
	private Collection<String> permute(int size, Collection<String> geneUniverse){
		Collection<String> rtrn=new HashSet();
		
		ArrayList list=new ArrayList();
		for(String str: geneUniverse){list.add(str);}
		
		for(int i=0; i<size; i++){
			int index=new Double(Math.random()*list.size()).intValue();
			String randomGene=(String)list.remove(index);
			rtrn.add(randomGene);
		}
		
		return rtrn;
	}
	
	//KSScore each geneSet
	private MatrixWithHeaders KSScoreGeneSets(MatrixWithHeaders scoreList, Map<String, Collection<String>> geneSets){
		MatrixWithHeaders rtrn=new MatrixWithHeaders(new ArrayList(geneSets.keySet()), scoreList.getColumnNames());
		
		//go through each scored experiments
		for(int i=0; i<scoreList.getNumberColumns(); i++){
			//sort it
			MatrixWithHeaders rankedList=scoreList.sortList(i);
			//score each gene set for that experiment
			for(String geneSet: geneSets.keySet()){
				double[] KSScore=KSTest.KSScores(rankedList, geneSets.get(geneSet), i);
				rtrn.set(geneSet, i, KSScore[0]);
			}
		}
		
		return rtrn;
	}
	
	private MatrixWithHeaders scoreGeneSets(MatrixWithHeaders rankedList, Map<String, Collection<String>> geneSets){
		//score gene sets
		//extract gene set values
		MatrixWithHeaders geneSetScoreMatrix=new MatrixWithHeaders(new ArrayList(geneSets.keySet()), rankedList.getColumnNames());
		for(String geneSet: geneSets.keySet()){
			MatrixWithHeaders submatrix=rankedList.submatrixByRowNames(geneSets.get(geneSet));
			//Now score each gene set either by KS or maxmean
			double[] scores=scoreGeneSet(submatrix);
			geneSetScoreMatrix.setRow(geneSet, scores);
		}
		
		return geneSetScoreMatrix;
	}
	
	private double[] scoreGeneSet(MatrixWithHeaders matrix){
		//genes are the rows and columns are samples
		//want to compute a statisitc on the rows by columns
		
		double[] rtrn=new double[matrix.getNumberColumns()];
		
		for(int i=0; i<rtrn.length; i++){
			double[] vals=matrix.getColumn(i);
			//simple statistic for testing
			rtrn[i]=Statistics.maxmean(vals);
		}
		
		return rtrn;
	}
	
	
	
	
	public static void main(String[] args)throws Exception{
		if(args.length>3){
			File expressionFile=new File(args[0]);
			File clsFile=new File(args[1]);
			String save=args[2];
			String gmtFile=args[3];
			new GSAScoring(expressionFile, clsFile, new File(gmtFile), save);
			
		}
		else{System.err.println(usage);}
	}

	static String usage="V2 args[0]=expression file \n args[1]=cls file \n args[2]=save file \n args[3]=gmtFile";
}
