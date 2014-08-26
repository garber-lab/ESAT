package broad.core.siphy.tools.conservation;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.*;


import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Gene;

import org.apache.commons.math3.stat.descriptive.rank.Percentile;

import broad.core.siphy.EvolutionaryModel.OmegaFit;
import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;

public class ComputePercentageOfKmers {
	
	
	public ComputePercentageOfKmers(File BEDScores, String save, double alpha)throws IOException{
		Map<Alignments, List<OmegaFit>> omega=parseScores(BEDScores);
		Map<Alignments, Double> percentageMap=getPercentageOfKmers(omega, alpha);
		write(save, percentageMap);
		//write(save+".pval", omega, alpha);
	}
	
	public ComputePercentageOfKmers(File BEDScores, Collection<Gene> genes, String save, double alpha)throws IOException{
		Map<Alignments, List<OmegaFit>> omega=parseScores(BEDScores);
		Map<Gene, Double> percentageMap=getPercentageOfKmers(omega, genes, alpha);
		write(save, percentageMap);
	}
	
	
	
	private void write(String save, Map<? extends Object, Double> omega)throws IOException{
		FileWriter writer=new FileWriter(save);
		
		for(Object align: omega.keySet()){
			Double p=omega.get(align);
			writer.write(align+"\t"+p+"\n");
			
		}
		
		writer.close();
	}
	
	
	private Map<Alignments, List<OmegaFit>> parseScores(File file)throws IOException{
		Map<Alignments, List<OmegaFit>> rtrn=new TreeMap();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine;
		int counter=0;
        while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
        	String[] tokens=nextLine.split("\t");
        	
        	Alignments align=new Alignments(tokens[0], tokens[1], tokens[2]);
        	Alignments region=new Alignments(tokens[3]);
        	double omega=new Double(tokens[4]);
        	double LOD=new Double(tokens[5]);
        	double p=new Double(tokens[6]);
        	double treeLength=new Double(tokens[7]);
        	
        	OmegaFit fit=new OmegaFit(region, omega, LOD, p, treeLength);
        	
        	List<OmegaFit> list=new ArrayList();
        	if(rtrn.containsKey(align)){list=rtrn.get(align);}
        	list.add(fit);
        	rtrn.put(align, list);
        	if(counter%10000 ==0){System.err.println(counter+" "+align);}
        	
        	counter++;
        }
		
        reader.close();
		return rtrn;
	}
	

	private Map getPercentageOfKmers(Map<Alignments, List<OmegaFit>> map, double alpha){
		Map rtrn=new TreeMap();
		
		for(Alignments align: map.keySet()){
			List<OmegaFit> list=map.get(align);
			double p=percent(list, alpha);
			rtrn.put(align, p);
		}
		
		return rtrn;
	}
	
	private Map<Gene, Double> getPercentageOfKmers(Map<Alignments, List<OmegaFit>> map, Collection<Gene> genes, double alpha){
		Map<Gene, Double> rtrn=new TreeMap<Gene, Double>();
		
		for(Gene gene: genes){
			Collection<? extends Annotation> exons=gene.getSortedAndUniqueExons();
			List<OmegaFit> list=new ArrayList<OmegaFit>();
			for(Annotation exon: exons){
				List temp=map.get(exon);
				if(temp!=null){list.addAll(temp);}
			}
			double p=percent(list, alpha);
			rtrn.put(gene, p);
		}
				
		return rtrn;
	}
	
	
	private double percent(List<OmegaFit> list, double alpha){
		double sum=0;
		int count=0;
		
		for(OmegaFit fit: list){
			if(fit.getPVal()<alpha){sum++;}
			count++;
		}
		
		return sum/count;
	}
	
	private Map<Alignments, List<OmegaFit>> getTopPercentileElements(Map<Alignments, List<OmegaFit>> omega, double percentile){
		Map rtrn=new TreeMap();
		
		for(Alignments align: omega.keySet()){
			System.err.println(align);
			List<OmegaFit> list=omega.get(align);
			double crit=percentileCutoff(list, percentile);
			List<OmegaFit> filtered=filter(list, crit);
			rtrn.put(align, filtered);
		}
		
		return rtrn;
	}
	
	private List<OmegaFit> filter(List<OmegaFit> list, double crit){
		List<OmegaFit> rtrn=new ArrayList();
		
		for(OmegaFit fit: list){
			double score=fit.getLogOddsScore();
			if(score>crit){rtrn.add(fit);}
		}
		
		return rtrn;
	}
	
	private double percentileCutoff(List<OmegaFit> list, double percentile){
		double[] array=makeDoubleArray(list);
		Percentile per=new Percentile(percentile);
		double crit=per.evaluate(array, percentile);
		return crit;
	}
	
	private double[] makeDoubleArray(List<OmegaFit> list){
		double[] rtrn=new double[list.size()];
		
		int i=0;
		for(OmegaFit fit: list){
			rtrn[i++]=fit.getLogOddsScore();
		}
		
		return rtrn;
	}
	
	
	
	private Map<Alignments, List<OmegaFit>> computeOmegaPerExon(Set<Alignments> set, File modelFile, String alnDir, String alnFileFormat, int windowSize)throws Exception{
		Map rtrn=new TreeMap();
		for(Alignments align: set){
			System.err.println(align);
			String alnFile=alnDir+"/"+align.getChr()+".maf";
			List<OmegaFit> fit=EstimateOmegaPerExon.slideWindowComputeOmega(align, modelFile, alnFile, alnFileFormat, windowSize, windowSize-1);
			rtrn.put(align, fit);
		}
		return rtrn;
	}
	
	public static void main(String[] args)throws Exception{
		if(args.length>3){
			File bedScores=new File(args[0]);
			String save=args[1];
			double percentile=new Double(args[2]);
			Collection<Gene> genes=BEDFileParser.loadData(new File(args[3]));
			//double alpha=new Double(args[3]);
			new ComputePercentageOfKmers(bedScores, genes, save, percentile);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=BED Scores \n args[1]=save file \n args[2]=alpha \n args[3]=genes";
	
	
}
