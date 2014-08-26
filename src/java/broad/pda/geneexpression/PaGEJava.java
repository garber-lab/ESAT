package broad.pda.geneexpression;

import java.io.*;
import java.util.*;

import broad.core.datastructures.MatrixWithHeaders;
import broad.core.math.ComputeFDR;
import broad.core.math.EmpiricalDistribution;
import broad.core.math.Statistics;
import broad.core.util.ParseGCTFile;


public class PaGEJava {
	
	double fudgeFactor=.1;
	int numBins=1000;
	int numPerm=100;
	boolean isLog=true;

	MatrixWithHeaders FDRs;
	MatrixWithHeaders folds;
	MatrixWithHeaders expression;
	
	public PaGEJava(MatrixWithHeaders expression, Map<String, ArrayList> groupIndexes, String controlGroup){
		init(expression, groupIndexes, controlGroup);
	}
	
	public PaGEJava(File gctFile, File clsFile)throws Exception{
		MatrixWithHeaders expression=new MatrixWithHeaders(gctFile.getAbsolutePath());
		Map<String, ArrayList> groupIndexes=ParseGCTFile.getGroupIndexes(clsFile);
		String controlClass=ParseGCTFile.getControlClass(clsFile);
		init(expression, groupIndexes, controlClass);
	}
	
	private void init(MatrixWithHeaders expression, Map<String, ArrayList> groupIndexes, String controlGroup){
		this.expression=expression;
		MatrixWithHeaders[] diffExpression=this.computeDifferentialStatsForAllGroups(expression, groupIndexes, controlGroup);
		this.FDRs=diffExpression[0];
		this.folds=diffExpression[1];
	}
	
	public Map<String, Collection<String>> getSignificantGeneSets(double alpha, double minFold){
		Map<String, Collection<String>> rtrn=new TreeMap<String, Collection<String>>();
		
		for(String geneSet: FDRs.getColumnNames()){
			Set<String> set=new TreeSet<String>();
			for(String gene: FDRs.getRowNames()){
				double FDR=FDRs.get(gene, geneSet);
				double fold=folds.get(gene, geneSet);
				if(FDR<alpha && fold>minFold){set.add(gene);}
			}
			rtrn.put(geneSet, set);
		}
		
		return rtrn;
	}
	
	
	//TODO Remove duplicate computations ES vs MEF and MEF vs ES
	private MatrixWithHeaders[] computeDifferentialStatsForAllGroups(MatrixWithHeaders expression, Map<String, ArrayList> groupIndexes, String controlGroup){
		MatrixWithHeaders FDRs=null;
		MatrixWithHeaders folds=null;
		
		boolean useControlClass=false;
		if(controlGroup!=null && groupIndexes.containsKey(controlGroup)){
			useControlClass=true;
		}
		
		for(int i=0; i<groupIndexes.size(); i++){
			for(int j=(i+1); j<groupIndexes.size(); j++){
				String group=groupIndexes.keySet().toArray()[i].toString();
				String group2=groupIndexes.keySet().toArray()[j].toString();
				System.err.println(group+" "+group2);
			if(!useControlClass || group.equalsIgnoreCase(controlGroup)){
				String name=group+"_vs_"+group2;
				
				
				Map<String, double[]> g1=expression.submatrixByColumnIndex(groupIndexes.get(group)).toMap();
				Map<String, double[]> g2=expression.submatrixByColumnIndex(groupIndexes.get(group2)).toMap();
				int n1=groupIndexes.get(group).size();
				int n2=groupIndexes.get(group2).size();
				
				if(group.equalsIgnoreCase(group2)){
					name=group+"_vs_REST";
					ArrayList restOfGenes = ScoreDiffGenes.getREST(group, groupIndexes);
					g2=expression.submatrixByColumnIndex(restOfGenes).toMap();
					n2=restOfGenes.size();
				}
				
				Map<String, Double> fdr=computeDifferentiallyExpressedGenes(g1, g2, numPerm, n1, n2);
				Map<String, Double> fold=computeFoldChange(g1, g2, isLog);
				MatrixWithHeaders fdrMatrix=convertToMatrix(fdr, name);
				MatrixWithHeaders foldsMatrix=convertToMatrix(fold, name);
				if(FDRs==null){FDRs=fdrMatrix;}
				else{FDRs.appendColumns(fdrMatrix);}
				if(folds==null){folds=foldsMatrix;}
				else{folds.appendColumns(foldsMatrix);}
			}
			}
		}
		MatrixWithHeaders[] rtrn={FDRs, folds};
		return rtrn;
	}
	
	private MatrixWithHeaders convertToMatrix(Map<String, Double> map, String name){
		List<String> columns=new ArrayList();
		columns.add(name);
		MatrixWithHeaders rtrn=new MatrixWithHeaders(new ArrayList(map.keySet()), columns);
		
		for(String gene: map.keySet()){
			rtrn.set(gene, 0, map.get(gene));
		}
		
		return rtrn;
	}
	
	private Map<String, Double> computeDifferentiallyExpressedGenes(Map<String, double[]> g1, Map<String, double[]> g2, int numPerm, int n1, int n2){
		Map<String, Double> tstats=tstat(g1, g2);
		Map[] randomTStats=this.randomTFoldStats(g1, g2, numPerm, n1, n2);
		Map<String, Double> fdr=this.computeFDR(tstats, randomTStats);
		return fdr;
	}
	
	private Map<String, Double> computeFoldChange(Map<String, double[]> g1, Map<String, double[]> g2, boolean isLog){
		Map<String, Double> fold=fold(g1, g2, isLog);
		return fold;
	}
	
	private Map<String, Double> fold(Map<String, double[]> g1, Map<String, double[]> g2, boolean isLog){
		Map<String, Double> rtrn=new TreeMap();
		for(String gene: g1.keySet()){
			if(!gene.equalsIgnoreCase("header")){
				double[] gr1=g1.get(gene);
				double[] gr2=g2.get(gene);
				double fold=Statistics.absFold(gr1, gr2, isLog);
				//System.err.println(gene+" "+Statistics.average(gr1)+" "+Statistics.average(gr2)+" "+fold+" "+gr1.length+" "+gr2.length);
				rtrn.put(gene, fold);
			}
		}
		return rtrn;
	}
	
	private Map[] randomTFoldStats(Map<String, double[]> g1, Map<String, double[]> g2, int numPerm, int n1, int n2){
		Map<String, Double> tstats=tstat(g1, g2);
		Map<String, double[]> geneExpression=merge(g1,g2);
		Map[] randomTStats=new Map[numPerm];
		for(int i=0; i<numPerm; i++){
			if(i%10 ==0){System.err.println("Permutation: "+i);}
			Map[] rand=permute(geneExpression, n1, (n1+n2));
			randomTStats[i]=tstat(rand[0], rand[1]);
		}
		return randomTStats;		
	}
	
	
	
	private Map merge(Map<String, double[]> g1, Map<String, double[]> g2){
		Map rtrn=new TreeMap();
		for(String gene: g1.keySet()){
			double[] list1=g1.get(gene);
			double[] list2=g2.get(gene);
			double[] list=new double[list1.length+list2.length];
			for(int i=0; i<list1.length; i++){list[i]=list1[i];}
			for(int i=0; i<list2.length; i++){list[i+list1.length]=list2[i];}
			rtrn.put(gene, list);
		}
		return rtrn;
	}
	
	private Map[] permute(Map<String, double[]> g1, int n1, int n){
		Set<Integer> positions=new TreeSet();
		for(int i=0; i<n1; i++){
			positions.add(new Double(Math.random()*n).intValue());
		}
		
		Map rand1=new TreeMap();
		Map rand2=new TreeMap();
		for(String gene: g1.keySet()){
			double[] list=g1.get(gene);
			ArrayList l1=new ArrayList();
			ArrayList l2=new ArrayList();
			for(int i=0; i<n; i++){
				if(positions.contains(i)){l1.add(list[i]);}
				else{l2.add(list[i]);}
			}
			rand1.put(gene, l2a(l1));
			rand2.put(gene, l2a(l2));
		}
		Map[] rtrn={rand1, rand2};
		return rtrn;
	}
	
	
	private double[] l2a(ArrayList<Double> list){
		double[] rtrn=new double[list.size()];
		
		for(int i=0; i<rtrn.length; i++){rtrn[i]=list.get(i);}
		
		return rtrn;
	}
	

	private Map tstat(Map<String, double[]> g1, Map<String, double[]> g2){
		Map rtrn=new TreeMap();
		for(String gene: g1.keySet()){
			if(!gene.equalsIgnoreCase("header")){
				double[] gr1=g1.get(gene);
				double[] gr2=g2.get(gene);
				double tstat=Statistics.tstat(gr1, gr2, fudgeFactor);
				rtrn.put(gene, tstat);
			}
		}
		return rtrn;
	}
	
	private Map computeFDR(Map<String, Double> tstats, Map[] randomTStats){
		Map fdr=new TreeMap();
		EmpiricalDistribution[] dists=this.makeDistributions(tstats, randomTStats);
		
		for(String gene: tstats.keySet()){
			double t=tstats.get(gene);
			double FDR=ComputeFDR.FDR(dists[0], dists, t);
			fdr.put(gene, FDR);
		}
		return fdr;
	}

	public void write(String save, File chipfile, double alpha, double fold)throws IOException{
		FileWriter writer=new FileWriter(save);
		
		Map<String, String> mapping=ParseGCTFile.parseChipFile(chipfile);
		Map<String, Collection<String>> geneSets=this.getSignificantGeneSets(alpha, fold);
		
		Set<String> genes=new TreeSet();
		for(String geneSet: geneSets.keySet()){genes.addAll(geneSets.get(geneSet));}
		
		writer.write("PID\tName");
		for(String header: this.expression.getColumnNames()){
			writer.write("\t"+header);
		}
		writer.write("\n");
		
		for(String gene: genes){
			double[] expression=this.expression.getRow(gene);
			writer.write(gene+"\t"+mapping.get(gene));
			for(int i=0; i<expression.length; i++){writer.write("\t"+expression[i]);}
			writer.write("\n");
		}
		
		writer.close();
	}
	
	private EmpiricalDistribution[] makeDistributions(Map<String, Double> tstats, Map[] randomTStats){
		EmpiricalDistribution tDist=new EmpiricalDistribution(tstats.values(), numBins, true);
		EmpiricalDistribution[] randomTDist=new EmpiricalDistribution[randomTStats.length+1];
		randomTDist[0]=tDist;
		for(int i=0; i<randomTStats.length; i++){
			randomTDist[i+1]=new EmpiricalDistribution(randomTStats[i].values(), numBins, true);
		}
		return randomTDist;
	}
	
	public static void main(String[] args)throws Exception{
		if(args.length>3){
			File expressionFile=new File(args[0]);
			File clsFile=new File(args[1]);
			String save=args[2];
			File chipFile=new File(args[3]);
			PaGEJava page=new PaGEJava(expressionFile, clsFile);
			page.write(save, chipFile, .05, 2.0);
			page.FDRs.writeGCT(save+".FDRs.gct", chipFile);page.folds.writeGCT(save+".fold.gct", chipFile);
		}
		else{System.err.println(usage);}
	}

	static String usage="V3 args[0]=expression file \n args[1]=cls file \n args[2]=save file \n args[3]=chip file";
	
}
