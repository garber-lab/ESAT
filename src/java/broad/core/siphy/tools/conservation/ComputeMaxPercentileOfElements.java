package broad.core.siphy.tools.conservation;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.*;


import nextgen.core.annotation.Annotation;

import org.apache.commons.math3.stat.descriptive.rank.Percentile;

import broad.core.siphy.EvolutionaryModel.OmegaFit;
import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;

public class ComputeMaxPercentileOfElements {
	
	
	public ComputeMaxPercentileOfElements(File BEDScores, String save, double percentile)throws IOException{
		Map<Annotation, List<OmegaFit>> omega=parseScores(BEDScores);
		Map<Annotation, List<OmegaFit>> percentileMap=getTopPercentileElements(omega, percentile);
		write(save, percentileMap);
		//write(save+".pval", omega, alpha);
	}
	
	private void write(String save, Map<Alignments, List<OmegaFit>> omega, double alpha)throws IOException{
		FileWriter writer=new FileWriter(save);
		
		for(Alignments align: omega.keySet()){
			List<OmegaFit> list=omega.get(align);
			for(OmegaFit fit: list){
				if(fit.getPVal()<alpha){writer.write(align+"\t"+fit.getRegion()+"\t"+fit.getOmega()+"\t"+fit.getLogOddsScore()+"\t"+fit.getPVal()+"\n");}
			}
		}
		
		writer.close();
	}
	
	
	private Map<Annotation, List<OmegaFit>> parseScores(File file)throws IOException{
		Map<Annotation, List<OmegaFit>> rtrn=new TreeMap();
		
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
	

	public ComputeMaxPercentileOfElements(File BEDFile, File modelFile, String alnDir, String alnFileFormat, String save, int windowSize, double percentile)throws Exception{
		Set<Annotation> set=BEDFileParser.loadAlignmentData(BEDFile);
		Map<Annotation, List<OmegaFit>> omega=this.computeOmegaPerExon(set, modelFile, alnDir, alnFileFormat, windowSize);
		Map<Annotation, List<OmegaFit>> percentileMap=getTopPercentileElements(omega, percentile);
		write(save, percentileMap);
	}
	
	private Map<Annotation, List<OmegaFit>> getTopPercentileElements(Map<Annotation, List<OmegaFit>> omega, double percentile){
		Map rtrn=new TreeMap();
		
		for(Annotation align: omega.keySet()){
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
	
	private void write(String save, Map<Annotation, List<OmegaFit>> omega)throws IOException{
		FileWriter writer=new FileWriter(save);
		
		for(Annotation align: omega.keySet()){
			List<OmegaFit> fitList=omega.get(align);
			for(OmegaFit fit: fitList){
				writer.write(align+"\t"+fit.getRegion()+"\t"+fit.getOmega()+"\t"+fit.getLogOddsScore()+"\t"+fit.getPVal()+"\n");
			}
		}
		
		writer.close();
	}
	
	private Map<Annotation, List<OmegaFit>> computeOmegaPerExon(Set<Annotation> set, File modelFile, String alnDir, String alnFileFormat, int windowSize)throws Exception{
		Map rtrn=new TreeMap();
		for(Annotation align: set){
			System.err.println(align);
			String alnFile=alnDir+"/"+align.getChr()+".maf";
			List<OmegaFit> fit=EstimateOmegaPerExon.slideWindowComputeOmega(align, modelFile, alnFile, alnFileFormat, windowSize, windowSize-1);
			rtrn.put(align, fit);
		}
		return rtrn;
	}
	
	public static void main(String[] args)throws Exception{
		if(args.length>2){
			File bedScores=new File(args[0]);
			String save=args[1];
			double percentile=new Double(args[2]);
			//double alpha=new Double(args[3]);
			new ComputeMaxPercentileOfElements(bedScores, save, percentile);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=BED Scores \n args[1]=save file \n args[2]=percentile";
	
	
}
