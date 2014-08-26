package broad.pda.differentialExpression;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import nextgen.core.annotation.Gene;


import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.MatrixWithHeaders;
import broad.core.datastructures.Pair;
import broad.core.error.ParseException;
import broad.core.math.EmpiricalDistribution;
import broad.core.math.FDRDistribution;
import broad.core.math.Statistics;
import broad.core.util.GMTParser;
import broad.pda.annotation.BEDFileParser;
import broad.pda.gene.GeneWithIsoforms;
import broad.pda.geneexpression.ExpressionExperimentInfo;
import broad.pda.geneexpression.agilent.AgilentUtils;
import broad.projection.gsa.GeneSetEnrichment;

public class CisVsTrans {

	//int[] geneSizes={5,10,20};
	int numRandom=1000;
	double nominalAlpha=.01;
	
	//Find consecutive regions of enrichment
	public CisVsTrans(MatrixWithHeaders data, File neighborFile, File gmtFile, String save, File experimentInfoFile, File chipFile) throws IOException{
		//Initialize
		Map<String, Collection<String>> experimentInfo=AgilentUtils.parseExperimentInfoFileToGroups(experimentInfoFile);
		Map<String, Collection<String>> chipInfo=GMTParser.parseCHIPFileByName(chipFile);
		Map<String, Collection<String>> geneSets=GMTParser.ParseGeneSetFile(gmtFile, 0);
				
		//Get the neighbors from the file
		Map<String, Pair<String>> neighbors=getNeighbors(neighborFile, experimentInfo);
		
		//make gene sets
		Map<String, Map<String, Collection<String>>> geneSetsByLinc=makeGeneSets(neighbors, geneSets);
		geneSetsByLinc=convertGeneSets(geneSetsByLinc, data, chipInfo, 0);
		
		//make random gene sets
		Map<String, Collection<String>> randomGeneSets=makeRandomGeneSets(geneSets, numRandom);
		randomGeneSets=convert(randomGeneSets, data, chipInfo, 0);
		
		//compute normalized scores for geneSet
		//runGeneSetComparisons(data, experimentInfo, geneSetsByLinc, randomGeneSets, save);
		
		//Get best of the neighbors
		getBestNeighbor(data, geneSetsByLinc, randomGeneSets, save);
		
	}
	
	private void getBestNeighbor(MatrixWithHeaders data, Map<String, Map<String, Collection<String>>> geneSetsByLinc, Map<String, Collection<String>> randomGeneSets, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String linc: data.getColumnNames()){
			System.err.println(linc);
			Map<String, Collection<String>> geneSets=geneSetsByLinc.get(linc);
			if(geneSets!=null){
			double best=0;
			for(String geneSet: geneSets.keySet()){
				Collection<String> genes=geneSets.get(geneSet);
				double val=getBest(data, genes, linc);
				if(val!=0){best=val;}
			}
			writer.write(linc+"\t"+best+"\n");
		}else{System.err.println("Skipping "+linc);}
		}
		
		writer.close();
	}

	private double getBest(MatrixWithHeaders data, Collection<String> genes, String linc) {
		double[] vals=new double[genes.size()];
		
		int i=0;
		for(String gene: genes){
			vals[i]=data.get(gene, linc);
			i++;
		}
		
		return computeConsensusScore(vals);
	}

	private double computeConsensusScore(double[] vals) {
		//Return max if positive
		double max=max(vals);
		
		//Return min if negative
		double min=min(vals);
		
		return largestAbsValue(max, min);
		
	}
	private double max(double[] vals) {
		double max=-Double.MAX_VALUE;
		for(int i=0; i<vals.length; i++){
			max=Math.max(max, vals[i]);
		}
		return max;
	}
	
	private double min(double[] vals) {
		double max=Double.MAX_VALUE;
		for(int i=0; i<vals.length; i++){
			max=Math.min(max, vals[i]);
		}
		return max;
	}

	private double largestAbsValue(double max, double min) {
		double maxAbs=Math.abs(max);
		double minAbs=Math.abs(min);
		if(maxAbs>minAbs){return max;}
		return min;
	}
	
	private Map<String, Collection<String>> makeRandomGeneSets(Map<String, Collection<String>> geneSets, int numRandom) {
		Map<String, Collection<String>> rtrn=new TreeMap<String, Collection<String>>();
		
		ArrayList<String> list=new ArrayList<String>(geneSets.keySet());
		
		for(int i=0; i<Math.min(numRandom, list.size()); i++){
			int random=new Double(Math.random()*list.size()).intValue();
			String geneSet=list.remove(random);
			rtrn.put(geneSet, geneSets.get(geneSet));
		}
	
		return rtrn;
	}

	private Map<String, Map<String, Collection<String>>> makeGeneSets(Map<String, Pair<String>> neighbors, Map<String, Collection<String>> geneSets) throws IOException {
		Map<String, Map<String, Collection<String>>> rtrn=new TreeMap<String, Map<String, Collection<String>>>();
		
		for(String linc: neighbors.keySet()){
			Pair<String> genes=neighbors.get(linc);
			Map<String, Collection<String>> allSets=new TreeMap<String, Collection<String>>();
			allSets.putAll(geneSetsContainingGene(genes.getValue1(), geneSets));
			allSets.putAll(geneSetsContainingGene(genes.getValue2(), geneSets));
			rtrn.put(linc, allSets);
			System.err.println("Name "+linc);
		}
				
		return rtrn;
	}

	private Map<String, Collection<String>> geneSetsContainingGene(String goi, Map<String, Collection<String>> geneSets) {
		Map<String, Collection<String>> rtrn=new TreeMap<String, Collection<String>>();
		
		for(String geneSet: geneSets.keySet()){
			Collection<String> genes=geneSets.get(geneSet);
			boolean validSet=false;
			for(String gene: genes){
				if(gene.equalsIgnoreCase(goi)){validSet=true;}
			}
			if(validSet){rtrn.put(geneSet, genes);}
		}
				
		//System.err.println(goi+" "+rtrn);
		
		return rtrn;
	}

	private Map<String, Pair<String>> getNeighbors(File neighborFile, Map<String, Collection<String>> experimentInfo) throws IOException {
		Collection<String> lines=BEDFileParser.loadList(neighborFile.getAbsolutePath(), true);
		Map<String, Pair<String>> rtrn=new TreeMap<String, Pair<String>>();
		
		for(String line: lines){
			String[] tokens=line.split("\t");
			String name=tokens[1];
			Pair<String> genes=new Pair<String>(tokens[2], tokens[3]);
			if(experimentInfo.containsKey(name)){rtrn.put(name, genes);}
			else{System.err.println("Skipping "+name+" because not in truncated info file");}
		}
		
		return rtrn;
	}

	private void runGeneSetComparisons(MatrixWithHeaders data, Map<String, Collection<String>> samplesByName, Map<String, Map<String, Collection<String>>>geneSetsByLinc, Map<String, Collection<String>> randomGeneSets, String save) throws IOException {
		Collection<String> negatives=samplesByName.get("Control");
				
		ArrayList<String> geneSets1=new ArrayList<String>();
		for(String linc: geneSetsByLinc.keySet()){
			for(String geneSet: geneSetsByLinc.get(linc).keySet()){
				String name=linc+"_"+geneSet;
				geneSets1.add(name);
			}
		}
		
		ArrayList<String> list=new ArrayList<String>(samplesByName.keySet());
		list.remove("Control");
		MatrixWithHeaders maxMeanMatrix=new MatrixWithHeaders(geneSets1, list);
		MatrixWithHeaders KSMatrix=new MatrixWithHeaders(geneSets1, list);
		MatrixWithHeaders MMPValue=new MatrixWithHeaders(geneSets1, list);
		MatrixWithHeaders KSPValue=new MatrixWithHeaders(geneSets1, list);
		
		for(String name: samplesByName.keySet()){
			System.err.println(name);
			if(!name.equalsIgnoreCase("Control")){
				GeneSetEnrichment random=new GeneSetEnrichment(data, negatives, samplesByName.get(name), randomGeneSets, true, true);
				for(String linc: geneSetsByLinc.keySet()){
					System.err.println(linc);
					Map<String, Collection<String>> geneSets=geneSetsByLinc.get(linc);
					GeneSetEnrichment diff=new GeneSetEnrichment(data, negatives, samplesByName.get(name), geneSets, true, true);
									
					//diff has gene sets by scores
					MatrixWithHeaders maxMean=diff.getNormalizedMaxMeanEnrichments();
					MatrixWithHeaders randomMaxMean=random.getNormalizedMaxMeanEnrichments();
					MatrixWithHeaders ks=diff.getNormalizedKSEnrichments();
					MatrixWithHeaders randomKS=random.getNormalizedKSEnrichments();
					
										
					EmpiricalDistribution mmDist=makeDist(randomMaxMean);
					EmpiricalDistribution ksDist=makeDist(randomKS);
					
					for(String geneSet: maxMean.getRowNames()){
						double score=maxMean.get(geneSet, 0);
						double[] randomScores=randomMaxMean.getColumn(0);
						double normScore=(score-Statistics.average(randomScores))/Statistics.stdev(randomScores);
						String name1=linc+"_"+geneSet;
						maxMeanMatrix.set(name1, name, normScore);
						double p=mmDist.getPValue(score, 1.0/randomScores.length);	
						MMPValue.set(name1, name, p);
						
						score=ks.get(geneSet, 0);
						randomScores=randomKS.getColumn(0);
						normScore=(score-Statistics.average(randomScores))/Statistics.stdev(randomScores);
						KSMatrix.set(name1, name, normScore);
						p=ksDist.getPValue(score, 1.0/randomScores.length);
						KSPValue.set(name1, name, p);
					}
					
					//TODO
					//get max/min normalized KS, maxmean as the score for the linc and assign it
					//assign normalized score (maxMean/average(random)))
					//assign p-value/fdr
				}
		}
		}
		
			
			
		maxMeanMatrix.writeGCT(save+".MaxMean.gct");
		KSMatrix.writeGCT(save+".KS.gct");
		KSPValue.writeGCT(save+".ksFDR.gct");
		MMPValue.writeGCT(save+".mmFDR.gct");
		
		
		filter(maxMeanMatrix, MMPValue);
		filter(KSMatrix, KSPValue);
		
		maxMeanMatrix.writeGCT(save+".mmFiltered.gct");
		KSMatrix.writeGCT(save+".ksFiltered.gct");
		
		//KS.writeGCT(save+".KS.gct");
		//maxMeanFDR.writeGCT(save+".maxMeanFDR.gct");
		//KSFDR.writeGCT(save+".KSFDR.gct");
	}
	
	private void filter(MatrixWithHeaders maxMeanMatrix, MatrixWithHeaders mMPValue) {
		for(String gene: maxMeanMatrix.getRowNames()){
			for(String experiment: maxMeanMatrix.getColumnNames()){
				double pvalue=mMPValue.get(gene, experiment);
				if(pvalue>this.nominalAlpha){maxMeanMatrix.set(gene, experiment, 0);}
			}
		}
	}

	private EmpiricalDistribution makeDist(MatrixWithHeaders randomMaxMean) {
		return new EmpiricalDistribution(randomMaxMean.getColumn(0), 200);
	}

	private Map<String, Map<String, Collection<String>>> convertGeneSets(Map<String, Map<String, Collection<String>>> geneSets, MatrixWithHeaders data, Map<String, Collection<String>> chipInfo, int minNum) {
		Map<String, Map<String, Collection<String>>> rtrn=new TreeMap<String, Map<String, Collection<String>>>();
		
		for(String linc: geneSets.keySet()){
			Map<String, Collection<String>> map=geneSets.get(linc);
			rtrn.put(linc, convert(map, data, chipInfo, minNum));
		}
		
		return rtrn;
	}
	
	private Map<String, Collection<String>> convert(Map<String, Collection<String>> geneSets, MatrixWithHeaders data, Map<String, Collection<String>> chipInfo, int minNum) {
		Map<String, Collection<String>> rtrn=new TreeMap<String, Collection<String>>();
		
		for(String geneSet: geneSets.keySet()){
			Collection<String> genes=geneSets.get(geneSet);
			Collection<String> pids=convert(genes, data, chipInfo);
			if(pids.size()>minNum){rtrn.put(geneSet, pids);}
		}
		
		return rtrn;
	}

	private Collection<String> convert(Collection<String> genes, MatrixWithHeaders data, Map<String, Collection<String>> chipInfo) {
		//First check if the current gene names are in the row names
		//if so return itself
		//Else convert
		Collection<String> rtrn=new TreeSet<String>();
		
		int found=0;
		int notFound=0;
		
		for(String gene: genes){
			if(data.hasRow(gene)){rtrn.add(gene); found++;}
			else{
				Collection<Integer> indexes=data.getIndecesForRowDescription(gene);
				if(indexes!=null){
					for(Integer index: indexes){rtrn.add(data.getRowName(index));}
					found++;
				}
				else if(chipInfo.get(gene)!=null){
					Collection<String> pids=chipInfo.get(gene);
					for(String pid: pids){
						if(data.containsRow(pid)){rtrn.add(pid);} //Only if the gene is actually in the data matrix do we want to add it
						found++;
					}
				}
				else{notFound++;}
			}
			
		}
		//System.err.println("Found "+(((double)found/(found+notFound)))*100+"%");
		
		return rtrn;
	}
	
	
	public static void main(String[] args)throws IOException, ParseException{
		if(args.length>5){
			MatrixWithHeaders data=new MatrixWithHeaders(args[0]);
			File neighborFile=new File(args[1]);
			File gmtFile=new File(args[2]);
			String save=args[3];
			File experimentInfoFile=new File(args[4]);
			File chipFile=new File(args[5]);
			
			new CisVsTrans(data, neighborFile, gmtFile, save, experimentInfoFile, chipFile);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=data \n args[1]=neighbor file \n args[2]=gmtFile \n args[3]=save \n args[4]=experiment Info \n args[5]=chip file";
	
}
