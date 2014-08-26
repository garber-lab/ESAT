package broad.pda.geneexpression;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import nextgen.core.annotation.Gene;



import Jama.Matrix;
import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.MatrixWithHeaders;
import broad.core.datastructures.IntervalTree.Node;
import broad.core.error.ParseException;
import broad.core.math.CombinationGenerator;
import broad.core.math.Statistics;
import broad.core.util.CLUtil;
import broad.core.util.ParseGCTFile;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.GeneTools;
import broad.pda.gene.GeneWithIsoforms;
import broad.pda.geneexpression.clustering.Kmeans;

public class geneExpUtils {

	static String usage="Usage: GeneExpUtils -task <task name> "+
	"\n\tNeighborCorrelation:\n calculate the correlation in expression level of a gene and its neighbors: \n\t\t -set1 <set of genes(BED format supported only)> \n\t\t -set2 <all genes that are potential neighbors> \n\t\t -set1ExpTab <GCT file of gene set> \n\t\t -set2ExpTab <GCT file of neighbor genes> \n\t\t -set1NameMap <map of names between set1's BED file and GCT file> \n\t\t -set2NameMap <map of names between set2's BED file and GCT file> -out <Output file>" +
	"\n\tMakeRpkmMat:\n generate an expression matrix for a ref set of genomic loci from a set f RNASeq data sets: \n\t\t -ref <reference BED file> \n\t\t -lst <file with a lst of bed file with score=expMeasure; lst structure <file path>\t<name>> \n\t\t -exactLoci <optional flag if exp vale is infered by an isoform that is identical to the reference isoform>  -extraField <optional, extra field column that has the score> \n\t\t -pvalField <optional, extra field column with enrichment score pval> \n"+
	"\n\tAllPairsCorr:\n calculate the pearson correlation between all pairs of genes from 2 sets: \n\t\t -gct1<gct file for set1 >\n\t\t  -gct2<gct file for set2>\n\t\t  -outPrefix<prefix for output files>\n\t\t  -numPerm <int> \n"+
	"\n\tGetSignificantCorrSet\n: get all genes that are significantly correlated with another set: \n\t\t -set1 <list of genes>\n\t\t -corrMat<all pairs correlation mat between set1 to set2 >\n\t\t -pvalMat<all pairs correlation p-value mat between set1 to set2 >\n\t\t -outPrefix<set2 out fname>\n"+
	"\n\tCalcGlobalRPKM\n: calculate a global RPKM value based on Scripture's score output: \n\t\t -set1 <output bed of scripture's score task> \n\t\t -out <outfile>"+
	"\n\tclusterByJS\n: partition exp matrix to the closest centroid using the Jensen-Shannon divergence metric: \n\t\t -exp <exp matrix (def: GCT format)> \n\t\t -numHeader<num header rows> \n\t\t -numHeaderCol <num header column> \n\t\t -centroids <same format GCT file> \n\t\t -reportJSSpecificity <optional, report all specificity scores> \n\t\t -outprefix <string> \n" +
	"\n\tkmeans\n: apply kmeans clustering on a gct file: \n\t\t -exp <exp matrix (def: GCT format)> \n\t\t  -centroids <optional; same format GCT file with predefined centroids> \n\t\t -refCentroids <optional; same format GCT file with predefined centroids; will calculate distance from these patterns>  \n\t\t -metric <JS/pearson/euclidean>\n\t\t -k <number of clusters> \n\t\t -chooseK <OPTIONAL \"x,y,z...\" comma seperated list of Ks to run and report Sillouette> \n\t\t -outprefix <string> \n" +
	"\n";
	
	static int permCashSize=10000;
	static int permNum=1000;
	
	public static void main(String [] args) throws Exception  {
		ArgumentMap argmap = CLUtil.getParameters(args, usage, "orient");
		
		if("NeighborCorrelation".equalsIgnoreCase(argmap.getTask())) {
			String set1bed = argmap.getMandatory("set1");
			String set2bed = argmap.getMandatory("set2");
			String set1gct = argmap.getMandatory("set1ExpTab");
			String set2gct = argmap.getMandatory("set2ExpTab");
			String set1map = argmap.getMandatory("set1NameMap");
			String set2map = argmap.getMandatory("set2NameMap");
			boolean revmap1 = argmap.containsKey("revmap1") ? true : false;
			boolean revmap2 = argmap.containsKey("revmap2") ? true : false;
			
			BufferedWriter bw = argmap.getOutputWriter();
			CalcNeighborsCorrelation(set1bed,set2bed,set1gct,set2gct,set1map,set2map,bw,revmap1,revmap2);
			bw.close();
		}
		else if ("collapseByGene".equalsIgnoreCase(argmap.getTask())){
			String setgct = argmap.getMandatory("setExpTab");
			String setmap = argmap.getMandatory("setNameMap");
			boolean revmap = argmap.containsKey("revmap") ? true : false;
			BufferedWriter bw = argmap.getOutputWriter();
			collapseExpMatByGene(setgct,setmap,bw,revmap);
			bw.close();
		}
		else if ("MakeRpkmMat".equalsIgnoreCase(argmap.getTask())){
			String refLociBed= argmap.getMandatory("ref");
			String bedLst=argmap.getMandatory("lst");
			boolean byExactLoci = argmap.containsKey("exactLoci")? true: false;
			int extraField =  (argmap.containsKey("extraField")? Integer.valueOf(argmap.get("extraField")): -1);
			int pvalField =  (argmap.containsKey("pvalField")? Integer.valueOf(argmap.get("pvalField")): -1);
			BufferedWriter bw = new BufferedWriter( new FileWriter (argmap.getMandatory("outf"))); 
			makeRPKMmat(refLociBed,bedLst,byExactLoci,bw,extraField,pvalField);
			bw.close();
		}
		else if("AllPairsCorr".equalsIgnoreCase(argmap.getTask())) {
			String set1gct = argmap.getMandatory("gct1");
			String set2gct = argmap.getMandatory("gct2");
			String outPrefix = argmap.getMandatory("outPrefix");
			Integer permNum =new Integer(argmap.getMandatory("numPerm"));
			BufferedWriter bw1 = new BufferedWriter( new FileWriter (outPrefix+".Corr"));
			BufferedWriter bw2 = new BufferedWriter( new FileWriter (outPrefix+".Pval"));
			CalcAllPairsCorrelation(set1gct,set2gct,bw1,bw2,permNum);
			bw1.close();
			bw2.close();
		}
		else if("GetSignificantCorrSet".equalsIgnoreCase(argmap.getTask())){
			String set1 = argmap.getMandatory("set1");
			String corrMatF = argmap.getMandatory("corrMat");
			String pvalMatF = argmap.getMandatory("pvalMat");
			String outPrefix = argmap.getMandatory("outPrefix");
			BufferedWriter bw1 = new BufferedWriter( new FileWriter (outPrefix+".corrGeneLst"));
			BufferedWriter bw2 = new BufferedWriter( new FileWriter (outPrefix+".corrGeneMat"));
			BufferedWriter bw3 = new BufferedWriter( new FileWriter (outPrefix+".pvalGeneMat"));
			
			BufferedWriter bw = argmap.getOutputWriter();
			GetSignificantCorrSet(set1,corrMatF,pvalMatF,bw1,bw2,bw3);
			bw1.close();
			bw2.close();
			bw3.close();
			
		}
		else if("CalcGlobalRPKM".equalsIgnoreCase(argmap.getTask())){
			String set1 = argmap.getMandatory("set1");	
			BufferedWriter bw = argmap.getOutputWriter();
			CalcGlobalRpkmFromScriptureOutput(set1,bw);
			bw.close();
			
		}
		else if ("clusterByJS".equalsIgnoreCase(argmap.getTask())){
			String expMat_F = argmap.getMandatory("exp");
			String clusterCentroids_F =argmap.getMandatory("centroids");
			Double numHeader=argmap.getDouble("numHeader");
			Double numHeaderCol=argmap.getDouble("numHeaderCol");
			String outprefix = argmap.getMandatory("outprefix");
			boolean reportJSSpecificity= argmap.containsKey("reportJSSpecificity")? true: false;
			clusterByJS(expMat_F,clusterCentroids_F,numHeader,outprefix,reportJSSpecificity);
		
		}
		else if ("kmeans".equalsIgnoreCase(argmap.getTask())){
			String expMat_F = argmap.getMandatory("exp");
			int k= argmap.getInteger("k");
			String outprefix = argmap.getMandatory("outprefix");
			String metric = argmap.getMandatory("metric");
			String clusterCentroids_F= argmap.containsKey("centroids")? argmap.getMandatory("centroids"):"" ;
			String refCentroids_F= argmap.containsKey("refCentroids")? argmap.getMandatory("refCentroids"):"" ;
			String chooseK = argmap.containsKey("chooseK")? argmap.getMandatory("chooseK"):"" ;
			if (chooseK.equals(""))
				runKmeans(expMat_F,clusterCentroids_F,metric,k,outprefix,refCentroids_F);
			else
				runKmeansWithSilhouette(expMat_F,clusterCentroids_F,metric,k,outprefix,refCentroids_F,chooseK);
			
		}
		
		
		else{System.err.println(usage);}
	}



	


	
	private static void collapseExpMatByGene(String setgct, String setmapF,
			BufferedWriter bw, boolean revmap) throws IOException, ParseException {
		MatrixWithHeaders setexp= new MatrixWithHeaders(setgct);
		Map <String, List<String>> setmap= loadNameMap(setmapF,revmap);
		ArrayList<String> nameLst= new ArrayList<String>();
		nameLst.addAll((setmap.keySet()));
		MatrixWithHeaders res = new MatrixWithHeaders (new Matrix(setmap.keySet().size(),setexp.getColumnNames().size()),nameLst,setexp.getColumnNames());
		
		for (String name:setmap.keySet()){
			
			List<String> probeLst= setmap.get(name);
			MatrixWithHeaders currMat = new MatrixWithHeaders (new Matrix(probeLst.size(),setexp.getColumnNames().size()),probeLst,setexp.getColumnNames());
			for (String probe:probeLst){
				String probeU=probe.toUpperCase();
				if (setexp.containsRow (probeU) ){
					double [] a= setexp.getRow(probeU);
					currMat.setRow(probe,setexp.getRow(probeU));
				}
			}
			double [] newRow= currMat.getMedianOverAllRows();
			res.setRow(name,  newRow);
		}
		
		res.writeGCT(bw);
		
	}


	private static void CalcNeighborsCorrelation(String set1bed,
			String set2bed, String set1gct, String set2gct, String set1mapF,
			String set2mapF, BufferedWriter bw,boolean revmap1,boolean revmap2) throws IOException, ParseException {
		
		BEDFileParser set1=new BEDFileParser(set1bed);
		BEDFileParser set2=new BEDFileParser(set2bed);
		MatrixWithHeaders set1exp= new MatrixWithHeaders(set1gct);
		MatrixWithHeaders set2exp= new MatrixWithHeaders(set2gct);
		Map <String, List<String>> set1map= loadNameMap(set1mapF,revmap1);
		Map <String, List<String>> set2map= loadNameMap(set2mapF,revmap2);
		
		
		//Find all neighbors
		List <Gene> genes= new LinkedList();
		List <Gene> neighbors= new LinkedList();
		
		List <Gene> set1Lst =set1.GetGenes();
		for (Gene g: set1Lst){
			getRightNeighbors(g,set2,genes,neighbors);
			getLeftNeighbors(g,set2,genes,neighbors);
		}
		
		/*  This code can be re-activated once I find a random combination generator
		Matrix permCash=new Matrix(permCashSize, set2exp.columnDimension());
		//pre-compute col permutations
		for (int i=0; i<permCashSize; i++){
			permCash.setRow(i,Statistics.randomPermutation(set2exp.columnDimension()));
		}
		CombinationGenerator combGenerator= new CombinationGenerator (permCashSize,permNum);
		*/
		
		// Generate a random distribution of correlation values by looking at all random pairs
		ArrayList <Double> nullDistribution = getPairsCorrNullDistribution(set1exp, set2exp);
		Collections.sort(nullDistribution);
		BufferedWriter bw2 = new BufferedWriter(new FileWriter("pairCorrNullDistribution.txt"));
		bw2.write(nullDistribution.toString());
		bw2.close();
		
		
		//Calc Pearson correlation
		for (int i=0;i<genes.size();i++){
			Gene g=genes.get(i);
			Gene n=neighbors.get(i);
			if (g.getName() != null  & n.getName() != null){
				List<String> gprobes=set1map.get(g.getName());
				List<String> nprobes=set2map.get(n.getName());
				if (gprobes != null  & nprobes != null){
				for (String gprob: gprobes){
					String clngprob=gprob.toUpperCase().replaceAll(" ", "");
					if (set1exp.containsRow(clngprob)){
						double[] g_exp=set1exp.getRow(clngprob);
						for (String nprob:nprobes){
							String clnNprob=nprob.toUpperCase().replaceAll(" ", "");
							if (set2exp.containsRow(clnNprob)){
								double[] n_exp=set2exp.getRow(clnNprob);
								double corr = Statistics.pearsonDistance(g_exp,n_exp);
								//int[] selectedPerms= combGenerator.getNext();
								//double pval = calcEmpiricalPval(g_exp,n_exp,permCash,selectedPerms);
								//for (int p=0; p<permNum; p++){
								//	permCash.setRow(p,Statistics.randomPermutation(set2exp.columnDimension()));
								//}
								//double pval = calcEmpiricalPval(g_exp,n_exp,permCash);
								
								double pval = calcEmpiricalPval_randomPairs(g_exp,n_exp,nullDistribution);
								
								String str=g.getName()+"\t"+gprob+"\t"+n.getName()+"\t"+nprob+"\t"+corr+"\t"+pval+"\n";
								bw.write(str);
							}
						}
					}
				}
			}
			}
			
		}
		
		
	}


	
	
	private static ArrayList<Double> getPairsCorrNullDistribution(
			MatrixWithHeaders set1exp, MatrixWithHeaders set2exp) {
		
		ArrayList<Double> allCorr=new ArrayList<Double>();
		for (int i=0; i<set1exp.rowDimension(); i++){
			for (int j=0; j<set2exp.rowDimension();j++)
				allCorr.add(Statistics.pearsonDistance(set1exp.getRow(i),set2exp.getRow(j)));
		}
		
		return allCorr;
	}


	private static double calcEmpiricalPval_permutePattern(double[] gExp, double[] nExp,
			Matrix permCash ) {
		
		double corr = Statistics.pearsonDistance(gExp,nExp);
		
		double larger=0.0;
		ArrayList <Double> randomCor= new ArrayList <Double>() ;
		for (int i=0; i < permNum; i++){
			//double [] newNexp=reoderArr(nExp,permCash.getRow(selectedPerms[i]));
			double [] newNexp=reorderArr(nExp,permCash.getRow(i));
			if (newNexp != null) randomCor.add(i,Statistics.pearsonDistance(gExp,newNexp));
			if ((corr > 0 & randomCor.get(i)>corr) | (corr < 0 & randomCor.get(i)<corr)) larger++;
				
		}
		
		return (larger/randomCor.size());
	}


	private static double[] reorderArr(double[] arr, double[] ixs) {
		
		double[] narr = new double[arr.length];
		for (int i=0; i<arr.length; i++)
		{
			int ix=(int)ixs[i]-1;
			if (ix < arr.length) 
			narr[i]=arr[ix];
			else
				return null;
		}
		return narr;
	}
	
	private static ArrayList<Double> reorderArr(ArrayList<Double> arr,ArrayList<Double> ixs) {
		ArrayList<Double> narr = new ArrayList<Double>();
		for (int i=0; i<arr.size(); i++)
		{
			int ix=(int) (ixs.get(i)-1);
			if (ix < arr.size()) 
				narr.add(i,arr.get(ix));
			else
				return null;
		}
		return narr;
	}


	//Assumes that the nullDistribution is sorted
	private static double calcEmpiricalPval_randomPairs(double[] gExp, double[] nExp,
			ArrayList <Double>  sortedNullDistribution ) {
		
		double corr = Statistics.pearsonDistance(gExp,nExp);
		return (calcEmpiricalPval(sortedNullDistribution,corr));
	}

	private static double calcEmpiricalPval(ArrayList <Double>  sortedNullDistribution , double val){
		double larger=0.0;
		int i;
		if (val > 0) {
			i=sortedNullDistribution.size()-1;
			while (sortedNullDistribution.get(i) >= val) {larger++; i--;}
		}
		else{
			i=0;
			while (sortedNullDistribution.get(i) <= val) {larger++; i++;}
		}
		
		return (larger/sortedNullDistribution.size());
	}

	private static Map<String, List<String>> loadNameMap(String mapName,boolean reverseMap) throws IOException {
		Map<String, List<String>> rtrn= new HashMap<String, List<String>>();
		String line;
		BufferedReader br = new BufferedReader(new FileReader(mapName));
		
		while((line = br.readLine()) != null){
			String [] lineSplit = line.split("\t");
			
			String n1=lineSplit[0];
			String n2=lineSplit[1];
			if (reverseMap){
				n1=lineSplit[1];
				n2=lineSplit[0];
			}
			if (! rtrn.containsKey(n1)){
				List<String> lst = new LinkedList<String>();
				rtrn.put(n1, lst);
			}
			rtrn.get(n1).add(n2);
		}
		br.close();
		return rtrn;
	}


	private static void getLeftNeighbors(Gene g, BEDFileParser set2,
			List<Gene> genes, List<Gene> neighbors) {

		IntervalTree<GeneWithIsoforms> tree=set2.getChrTree(g.getChr());  
		Node<GeneWithIsoforms> LeftNode=tree.max(g.getStart(), g.getEnd()-1);
		if (LeftNode != null){
			genes.add(g);
			neighbors.add(LeftNode.getValue());
			int size=g.getStart()-LeftNode.getValue().getEnd();
			System.err.println(g.getName()+"\tLeft: "+LeftNode.getValue().getName()+"\t distance " + size);
		}
		
	}


	private static void getRightNeighbors(Gene g, BEDFileParser set2,
			List<Gene> genes, List<Gene> neighbors) {
		IntervalTree<GeneWithIsoforms> tree=set2.getChrTree(g.getChr());  
		Node<GeneWithIsoforms> RightNode=tree.min(g.getStart()+1, g.getEnd());
		if (RightNode != null){
			genes.add(g);
			neighbors.add(RightNode.getValue());
			int size=RightNode.getValue().getStart()-g.getEnd();
			System.err.println(g.getName()+"\tRight: "+RightNode.getValue().getName()+"\t distance " + size);
		}
		
	}
			
	
	//mem efficient version
	private static void makeRPKMmat(String refLociBed, String bedLst,
			boolean byExactLoci, BufferedWriter bw, int extraField, int pvalField) throws IOException {
		
		BEDFileParser ref=new BEDFileParser(refLociBed);
		
			//Mem efficient
		BufferedReader br = new BufferedReader(new FileReader(bedLst));
		String line;
		Map<String,String> setsPath= new HashMap<String,String>();
		 while((line = br.readLine()) != null) {
				line = line.trim();
				String [] lineSplit = line.split("\t");
				setsPath.put(lineSplit[1], lineSplit[0]);	
		 }		
		
		LinkedList<String> setsLst=new LinkedList(setsPath.keySet());
		Collections.sort(setsLst);
		
		Map <String, ArrayList<Double>> resMap= new HashMap<String, ArrayList<Double>>();
		int i=0;
		for (String setName: setsLst){
			BEDFileParser bed= new BEDFileParser(setsPath.get(setName));
			Iterator<String> chrIt=ref.getChromosomeIterator();
			while (chrIt.hasNext()){
				String chr=chrIt.next();
				Iterator<GeneWithIsoforms> geneIt= ref.getChrTree(chr).valueIterator();
				while (geneIt.hasNext()){
					GeneWithIsoforms loci =geneIt.next();
					if (!resMap.containsKey(loci.getName())){
						ArrayList<Double> arr=new ArrayList<Double>();
						resMap.put(loci.getName(),arr);
						
					}
					Gene  isoform;
					if (byExactLoci){
						isoform=bed.getExactIsoform(loci);
					}
					else
						isoform=GeneTools.selectHighScoringIsoform(loci,bed);
					double val=0;
					if (isoform != null){
							if (extraField>=0){
								if (pvalField!=-1 && Double.valueOf(isoform.getExtraFields(pvalField))>0.05)
									val=0;
								else
									val=Double.valueOf(isoform.getExtraFields(extraField));
							}
							else
								val=isoform.getBedScore();
					}
					ArrayList<Double> arr=resMap.get(loci.getName());
					arr.add(i,val);
					resMap.put(loci.getName(),arr);
				}
				
			}
			
			//Go to the next set 
			i++;
		}
		System.err.println("The number of genes included is: " + resMap.size());
		bw.write("#1.2\n");
		bw.write(resMap.size()+"\t"+setsLst.size()+"\n");
		bw.write("name\tdescription");
		for(String columnName : setsLst) {
			bw.write("\t"+columnName);
		}
		bw.write("\n");
		for (String name: resMap.keySet()){
			bw.write(name+"\t"+name);
			ArrayList<Double> a=resMap.get(name);
			for (int j=0;  j< setsLst.size();j++ ){
				bw.write("\t"+a.get(j));
			}
			bw.write("\n");
		}
		
		
		
		
	}
	
	
	
	
	
	@SuppressWarnings("unchecked")
	private static void CalcAllPairsCorrelation(String set1gct, String set2gct,BufferedWriter bwCorr, BufferedWriter bwPval,int permNum) throws IOException {
		
		Map<String, ArrayList<Double>> lincExp=ParseGCTFile.loadData(new File(set1gct));
		Map<String, ArrayList<Double>> geneExp=ParseGCTFile.loadData(new File(set2gct));
		CalcAllPairsCorrelation(lincExp, geneExp, bwCorr, bwPval, permNum) ;
	}
	
    @SuppressWarnings("unchecked")
	private static void CalcAllPairsCorrelation(Map<String, ArrayList<Double>> lincExp, Map<String, ArrayList<Double>> geneExp,BufferedWriter bwCorr, BufferedWriter bwPval,int permNum) throws IOException {
			
		
		int colNum=lincExp.get("header").size();
		lincExp.remove("header");
		geneExp.remove("header");
		
		LinkedList<String> geneLst=new LinkedList(geneExp.keySet());
		
		
		
		bwCorr.write("name");
		bwPval.write("name");
		for (int i=0; i<geneLst.size(); i++){
			bwCorr.write("\t"+geneLst.get(i));
			bwPval.write("\t"+geneLst.get(i));
		}
		bwCorr.write("\n");
		bwPval.write("\n");
		ArrayList<ArrayList<Double>> permOrder= generateRandomOrders(permNum,colNum);
		for (String linc: lincExp.keySet()){
			ArrayList<Double> lincVal=lincExp.get(linc);
			//ArrayList<ArrayList<Double>> lincsRandomPerms= generateRandomPerms(lincVal,permNum);
			ArrayList<ArrayList<Double>> lincsRandomPerms= generateRandomPerms(lincVal,permNum,permOrder);
			ArrayList<Double> sortedRandomCorr= generateRandomCorr(geneExp,lincsRandomPerms,geneLst);
			ArrayList<Double> resCorr=oneToAllCorr(lincVal,geneExp,geneLst);
			ArrayList<Double> resPval=calcEmpiricalPval(resCorr,sortedRandomCorr);
			bwCorr.write(linc);
			bwPval.write(linc);
		
			for (int i=0; i<geneLst.size(); i++){
				bwCorr.write("\t"+resCorr.get(i));
				bwPval.write("\t"+resPval.get(i));
			}
			bwCorr.write("\n");
			bwPval.write("\n");
		}
		
	}




	





	private static ArrayList<ArrayList<Double>> generateRandomOrders(int permNum2, int colNum) {
		ArrayList<Double> vec=new ArrayList<Double>();
		for (double i=0; i<colNum; i++)
			vec.add(i);
		return (generateRandomPerms(vec,permNum2));
	}





	private static ArrayList<Double> calcEmpiricalPval(ArrayList<Double> resCorr, ArrayList<Double> randomCorr) {

		ArrayList<Double> res = new ArrayList<Double>();
		for (int i=0; i<resCorr.size(); i++)
			res.add(i,calcEmpiricalPval(randomCorr,resCorr.get(i)));
		return res;
	}

	private static ArrayList<ArrayList<Double>> generateRandomPerms(ArrayList<Double> lincVal,int permNum) {
		ArrayList<ArrayList<Double>> res=new ArrayList<ArrayList<Double>>();
		ArrayList<Double> arr=new ArrayList<Double> (lincVal);
		for (int i=0;i<permNum; i++ ){
			Collections.shuffle(arr);
			ArrayList<Double> a=new ArrayList<Double>(arr);
			res.add(a);
		}
		
		return res;
	}

	private static ArrayList<ArrayList<Double>> generateRandomPerms(ArrayList<Double> lincVal, int permNum2,ArrayList<ArrayList<Double>> permOrder) {
		ArrayList<ArrayList<Double>> res=new ArrayList<ArrayList<Double>>();
		ArrayList<Double> arr=new ArrayList<Double> (lincVal);
		for (int i=0;i<permNum; i++ ){
			ArrayList<Double> a=reorderArr(arr,permOrder.get(i));
			res.add(arr);
		}
		
		return res;
	}
	
	




	private static ArrayList<Double> generateRandomCorr( Map<String, ArrayList<Double>> geneExp,
			ArrayList<ArrayList<Double>> lincsRandomPerms,LinkedList<String> geneLst) {
		ArrayList<Double> res=new ArrayList<Double>();
		for (int i=0; i<lincsRandomPerms.size(); i++)
			res.addAll(oneToAllCorr(lincsRandomPerms.get(i),geneExp,geneLst));
		Collections.sort(res);
		return res;
	}


	private static ArrayList<Double> oneToAllCorr(ArrayList<Double> lincVal,Map<String, ArrayList<Double>> geneExp, LinkedList<String> geneLst) {
		ArrayList<Double> res= new ArrayList<Double>();
		for (int i=0; i< geneLst.size();i++){
			ArrayList<Double> a=geneExp.get(geneLst.get(i));
			res.add(Statistics.pearsonDistance(lincVal,a));
		}
			
		return res;
	}



	private static void GetSignificantCorrSet(String set1F, String corrMatF,	String pvalMatF, BufferedWriter bw, BufferedWriter bw2, BufferedWriter bw3) throws IOException {
		
		HashMap<String,Integer> set1= new HashMap<String,Integer>();
		BufferedReader br= new BufferedReader(new FileReader(set1F));
		String line;
		while ((line = br.readLine())  != null){
			line.trim();
			set1.put(line,0);			
		}
		br.close();
		
		BufferedReader br_corr= new BufferedReader(new FileReader(corrMatF));
		BufferedReader br_pval= new BufferedReader(new FileReader(pvalMatF));
		
		// first line is set2 universe
		line=br_corr.readLine();
		String[] colNames=line.split("\t");
		line=br_pval.readLine();
		String[] pval_colNames=line.split("\t");
		for(int i=0; i<colNames.length ; i++){
			if (!colNames[i].equalsIgnoreCase(pval_colNames[i]))
				System.err.printf("Col %d differes between corr and pval mat\n", i);
		}
		bw2.write(line+"\n");
		bw3.write(line+"\n");
		
		int totalFound=0;
		int totalDefined=set1.keySet().size();
		String line1;
		String line2;
		HashMap<String,String> set2 =new HashMap<String,String> ();
		
		while((line1 = br_corr.readLine())  != null){
			String[] cur_corr= line1.split("\t");
			line2 = br_pval.readLine();
			String[] cur_pval= line2.split("\t");
			
			if (set1.containsKey(cur_corr[0])){
				if (set1.get(cur_corr[0])==0){ 
					set1.put(cur_corr[0], 1);
					totalFound++;
					bw2.write(line1+"\n");
					bw3.write(line2+"\n");
				}
												
				for (int i=1 ; i < cur_pval.length; i++){
					Double val= new Double(cur_pval[i]);
					Double cor= new Double(cur_corr[i]);
					if (cor>0.5 & val<=0.05){
						set2.put(pval_colNames[i],"");
					}
					//if (! cor.isNaN())
						//System.err.printf("%s\t%s\n",cur_corr[i],cur_pval[i] );
						
				}
			}
			if (totalFound==totalDefined)
				break;
			
		}
		
		br_corr.close();
		br_pval.close();
		
		for (String s: set2.keySet()){
			bw.write(s+"\n");
		}
		bw.close();
		
	}



	private static void CalcGlobalRpkmFromScriptureOutput(String set1,BufferedWriter bw) throws IOException {
		
		BEDFileParser bed= new BEDFileParser(set1);
		
		Iterator<String> chrIt=bed.getChromosomeIterator();
		double totalLibCount=0;
		while (chrIt.hasNext()){
			String chr=chrIt.next();
			Iterator<GeneWithIsoforms> geneIt= bed.getChrTree(chr).valueIterator();
			while (geneIt.hasNext()){
				GeneWithIsoforms loci =geneIt.next();
				double localRpkm= Double.valueOf(loci.getExtraFields(5));
				double avCov=Double.valueOf(loci.getExtraFields(4));
				if (localRpkm !=0 && avCov!=0 ){
					double chrTot=(avCov /localRpkm) * Math.pow(10,9);
					totalLibCount+=chrTot;
					System.err.println(chr+"\tTotal read: "+chrTot);
					break;
				}
			}
		}
		chrIt=bed.getChromosomeIterator();
		System.err.println("Total reds in lib "+totalLibCount);
		while (chrIt.hasNext()){
			String chr=chrIt.next();
			Iterator<GeneWithIsoforms> geneIt= bed.getChrTree(chr).valueIterator();
			while (geneIt.hasNext()){
				GeneWithIsoforms loci =geneIt.next();
				for (Gene g:loci.getAllIsoforms()){
					double globalRPKM= (Double.valueOf(g.getExtraFields(4))/totalLibCount) * Math.pow(10,9);
					g.setBedScore(globalRPKM);
					g.setExtraFields(globalRPKM,5);
					bw.write(g.toBEDwithBedScore()+"\n");
					
				}
			}
		}
		
		
	}


	private static void clusterByJS(String expMatF, String clusterCentroidsF,
			Double numHeader, String outprefix, boolean reportJSSpecificity) throws IOException, ParseException {
		
		BufferedWriter bw2=null;
		BufferedWriter bw3=null;
		
				
		MatrixWithHeaders centroids= new MatrixWithHeaders(clusterCentroidsF);
		int extraCluster=centroids.rowDimension()+1;
		normalizeToRelativeAbundance(centroids);
		double maxDimDist=calcNormJSFactor(centroids);
		
		if (reportJSSpecificity==false)
			bw2= new BufferedWriter (new FileWriter (outprefix + ".clusterMap.txt"));
		else{
			bw3= new BufferedWriter (new FileWriter (outprefix + ".JSspecificity.txt"));
			bw3.write("Name\tDescription");
			for (int i=0; i<centroids.rowDimension();i++)
				bw3.write("\t"+centroids.getRowName(i));
			bw3.write("\n");
		}
		
		BufferedReader br = new BufferedReader(new FileReader(expMatF));
		
		//read header 
		String line = br.readLine();
		
		if(line.startsWith("#1.2")) {
			line = br.readLine();
		}
	
		String [] dimensionsStr = line.split("\t");
		int expectedRowDimension = Integer.parseInt(dimensionsStr[0]);
		int expectedColumnDimension = Integer.parseInt(dimensionsStr[1]);
		
		String header = br.readLine();
		
		
	
		String [] columnNames = header.split("\t");
		List<String> columnNameList = new ArrayList<String> (columnNames.length - 2);
		if(columnNames.length - 2 != expectedColumnDimension) {
			System.err.println("WARNING: expected "+expectedColumnDimension+ " columns but read " +  (columnNames.length - 2) );
		}
		int lineNum=4;
		// read line by line and classify to centroid
		while( (line = br.readLine()) != null) {
			String [] info = line.split("\t");
			if(info.length != columnNames.length) {
				throw new ParseException("Line " + lineNum + " has " + info.length + " columns but header had " + columnNames.length + " columns");
			}
			String rowName=info[0];
			String rowDesc=info[1];
			double[] lineData = new double[(info.length - 2)];
			double sum=0.0;
			for(int i = 2 ; i < info.length; i++) {
				lineData[i-2]=Double.parseDouble(info[i]);
				sum+=lineData[i-2];
			}
			lineNum++;
			String[] centroid_res=new String[1];
			ArrayList<Double> distArr= new ArrayList<Double>();
			if (sum == 0)
				centroid_res[0]=String.valueOf(extraCluster);
			else{
				centroid_res[0] = findClosestCenterByJS(lineData,centroids,distArr);
				
			}
			if (reportJSSpecificity==false)
				bw2.write(rowName+"\t"+rowDesc+"\t"+centroid_res[0]+"\n");
			else{
				ArrayList<Double> specificityArr = calcJSspecificity(distArr,maxDimDist);
				bw3.write(rowName+"\t"+rowDesc);
				if (specificityArr.isEmpty()){
					for (int j=0; j<centroids.rowDimension(); j++)
					bw3.write("\t0");
				}
				else{	
				for (int j=0; j<specificityArr.size(); j++)
					bw3.write("\t"+specificityArr.get(j));
				}
				bw3.write("\n");
			}
		}
		
		if (reportJSSpecificity==false)
			bw2.close();
		else
			bw3.close();
		br.close();
			
		
		
		
	}


	private static ArrayList<Double> calcJSspecificity(ArrayList<Double> distArr, double maxDimDist) {

		ArrayList<Double> sp=new ArrayList<Double>();
		for (int i=0; i<distArr.size();i++){
			double d=distArr.get(i);
			double dn=1-d/maxDimDist;
			if(new Double(d).isNaN() || dn<0)
				sp.add(new Double(0));
			else
				sp.add(dn);
			
		}
		
		return sp;
	}







	private static double calcNormJSFactor(MatrixWithHeaders centroids) {

		double [] a= new double[centroids.columnDimension()];
		double [] b= new double[centroids.columnDimension()];
		for (int i=0; i<a.length; i++)
			{a[i]=0; b[i]=0;}
			
		a[0]=1;
		b[1]=1;
		double res = Statistics.JSDist(a,b);
		return res;
	}







	private static String findClosestCenterByJS(double [] lineData,
			MatrixWithHeaders centroids, ArrayList<Double> JS) {
		String res ="";
		double[] nLineData=normalizeToRelativeAbundance(lineData);
		int i=0;
		for (i=0; i< centroids.rowDimension(); i++){
			double[] currRow= centroids.getRow(i);
			JS.add(i,getJSDist(nLineData,currRow));
		}
		double m_val= Statistics.min(JS);
		for ( i=0;i<centroids.rowDimension(); i++){
			if (JS.get(i)==m_val)
				break;
		}
		res=centroids.getRowName(i);
		
		return res;
	}

	private static void normalizeToRelativeAbundance(MatrixWithHeaders mat) {
		for (int i=0; i<mat.rowDimension();i++){
			String row=mat.getRowName(i);
			double []vals=mat.getRow(i);
			mat.setRow(row, normalizeToRelativeAbundance(vals));
		}
	}


	private static double[] normalizeToRelativeAbundance(double[] vals) {
		int len=vals.length;
		double[] res= new double[len];
		double sum= 0;
		for (int i=0; i<len;i++) {sum+=vals[i];}
		for (int i=0; i<len;i++) {res[i]=(vals[i]/sum);}
		return res;
	}


	private static Double getJSDist(double[] p1, double[] p2) {
		int len=p1.length;
		double res=100;
		double [] pm= new double[len];
		for (int i=0; i<len;i++) {pm[i]=(p1[i]+p2[i])/2;}
		double hpm=Statistics.entropy(pm);
		double hp1=Statistics.entropy(p1);
		double hp2=Statistics.entropy(p2);
		res= hpm - ((hp1+hp2)/2);
		return res;
	}


	
	
	private static void runKmeans(String expMatF, String clusterCentroidsF,
			String metric, int k, String outprefix, String refCentroidsF) throws IOException, ParseException {

		MatrixWithHeaders expMat=new MatrixWithHeaders(expMatF);
		MatrixWithHeaders centroids=null;
		if (!(metric.equalsIgnoreCase("JS")||metric.equalsIgnoreCase("pearson")|| metric.equalsIgnoreCase("euclidean") ))
		{	System.err.println("only supported distance metrics are : JS, PEARSON, EUCLIDEAN");	
			return;	}
		
		if (!clusterCentroidsF.equalsIgnoreCase("")){
			 centroids= new MatrixWithHeaders(clusterCentroidsF);
			 if (metric.equalsIgnoreCase("JS"))
				 normalizeToRelativeAbundance(centroids);
		}
		Kmeans kmeans_res=new Kmeans(expMat,centroids,k,metric);
		
		BufferedWriter bw= new BufferedWriter (new FileWriter (outprefix + ".clusterMap.txt"));
		kmeans_res.writeClusters(bw);
		bw.close();
		
		bw= new BufferedWriter (new FileWriter (outprefix + ".clusterCentroids.gct"));
		kmeans_res.writeClustersCentroids(bw);
		bw.close();
		
		bw= new BufferedWriter (new FileWriter (outprefix + ".distFromCluster.gct"));
		kmeans_res.writeDistanceFromCentroids(bw);
		bw.close();
		
		if (! refCentroidsF.equalsIgnoreCase("")){
			
			centroids= new MatrixWithHeaders(refCentroidsF);
			 if (metric.equalsIgnoreCase("JS"))
				 normalizeToRelativeAbundance(centroids);
			 
			bw= new BufferedWriter (new FileWriter (outprefix + ".distFromInputCenteroids.gct"));
			kmeans_res.writeDistanceFromInputCentroids(centroids,bw);
			bw.close();
		}
		
		
	}
	
	
	private static void runKmeansWithSilhouette(String expMatF,
			String clusterCentroidsF, String metric, int k, String outprefix,
			String refCentroidsF, String chooseK) throws IOException, ParseException {


		MatrixWithHeaders expMat=new MatrixWithHeaders(expMatF);
		MatrixWithHeaders centroids=null;
		MatrixWithHeaders refcentroids=null;
		if (!(metric.equalsIgnoreCase("JS")||metric.equalsIgnoreCase("pearson")|| metric.equalsIgnoreCase("euclidean") ))
		{	System.err.println("only supported distance metrics are : JS, PEARSON, EUCLIDEAN");	
			return;	}
		
		if (!clusterCentroidsF.equalsIgnoreCase("")){
			 centroids= new MatrixWithHeaders(clusterCentroidsF);
			 if (metric.equalsIgnoreCase("JS"))
				 normalizeToRelativeAbundance(centroids);
		}
		
		if (! refCentroidsF.equalsIgnoreCase("")){
			refcentroids= new MatrixWithHeaders(refCentroidsF);
			 if (metric.equalsIgnoreCase("JS"))
				 normalizeToRelativeAbundance(refcentroids);
		}
		
		//extract K to run on 
		String[] kStrs=chooseK.split(",",-1);
		int[] Ks = new int[kStrs.length]; 
		for (int i=0; i<kStrs.length;i++)
			Ks[i]=Integer.valueOf(kStrs[i]);
		
		//Make output files
		String range=Ks[0]+"_"+Ks[Ks.length-1];
		BufferedWriter bw_geneS=new BufferedWriter (new FileWriter (outprefix + ".allGeneSilhouetteVals."+range+".txt"));
		BufferedWriter bw_allS=new BufferedWriter (new FileWriter (outprefix + ".KAvSilhouette."+range+".txt"));
		BufferedWriter bw_clusterS=new BufferedWriter (new FileWriter (outprefix + ".ClusterAvSilhouette."+range+".txt"));
	

		
		//Run on different Ks:
		for (int i=0; i<Ks.length; i++){
			k=Ks[i];
			System.err.println("Start runing k=" +k);
			Kmeans kmeans_res=new Kmeans(expMat,centroids,k,metric);
			kmeans_res.calcLightSilhouette();
			
			//Write cluster res
			BufferedWriter bw= new BufferedWriter (new FileWriter (outprefix + ".clusterMap."+k+".txt"));
			kmeans_res.writeClusters(bw);
			bw.close();
			
			bw= new BufferedWriter (new FileWriter (outprefix + ".clusterCentroids."+k+".gct"));
			kmeans_res.writeClustersCentroids(bw);
			bw.close();
			
			bw= new BufferedWriter (new FileWriter (outprefix + ".distFromCluster."+k+".gct"));
			kmeans_res.writeDistanceFromCentroids(bw);
			bw.close();
			
			if (refcentroids != null){
			bw= new BufferedWriter (new FileWriter (outprefix + ".distFromInputCenteroids."+k+".gct"));
			kmeans_res.writeDistanceFromInputCentroids(refcentroids,bw);
			bw.close();
			}
			
			//Write silhouette res
			kmeans_res.writeSilhouetteRes(k, bw_geneS, bw_allS, bw_clusterS);	
			System.err.println("K="+k+" Silhouette=" + kmeans_res.getSilhouette());
			
		}
		
		bw_geneS.close();
		bw_allS.close();
		bw_clusterS.close();
		
	}

	
	
}
