     package broad.core.siphy.tools.conservation;


import java.io.*;
import java.util.*;

import nextgen.core.annotation.Annotation;

import broad.core.siphy.EvolutionaryModel.OmegaFit;
import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;


public class ComputeKmersForRegion {
	
	public ComputeKmersForRegion(File BEDFile, File modelFile, String alnDir, String alnFileFormat, String save, int windowSize)throws Exception{
		Set<Annotation> set=BEDFileParser.loadAlignmentData(BEDFile);
		set=filter(set, windowSize);
		Map<Annotation, List<OmegaFit>> omega=computeOmegaPerExon(set, modelFile, alnDir, alnFileFormat, windowSize);
		write(save, omega);
	}
	
	private Set filter(Set<Annotation> set, int windowSize){
		Set rtrn=new TreeSet();
		
		for(Annotation align: set){
			if(align.length()>windowSize){rtrn.add(align);}
		}
		
		return rtrn;
	}
	
	private void write(String save, Map<Annotation, List<OmegaFit>> omega)throws IOException{
		FileWriter writer=new FileWriter(save);
		
		for(Annotation align: omega.keySet()){
			List<OmegaFit> fitList=omega.get(align);
			for(OmegaFit fit: fitList){
				writer.write(align+"\t"+fit.getRegion().toUCSC()+"\t"+fit.getOmega()+"\t"+fit.getLogOddsScore()+"\t"+fit.getPVal()+"\t"+fit.getTreeLength()+"\n");
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
		if(args.length>5){
			File bedFile=new File(args[0]);
			File modelFile=new File(args[1]);
			String alnDir=args[2];
			String alnFormat=args[3];
			String save=args[4];
			int windowSize=new Integer(args[5]);
			new ComputeKmersForRegion(bedFile, modelFile, alnDir, alnFormat, save, windowSize);
		}
		else{System.err.println(usage);}
	}
	
	static String usage="v1.3 args[0]=BED file \n args[1]=model file \n args[2]=alignment directory \n args[3]=alignment file format \n args[4]=save file \n args[5]=window size";
	
}
