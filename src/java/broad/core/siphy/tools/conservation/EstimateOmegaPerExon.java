package broad.core.siphy.tools.conservation;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Gene;

import broad.core.siphy.TreeScaler;
import broad.core.siphy.EvolutionaryModel.OmegaFit;
import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;

public class EstimateOmegaPerExon {

	public EstimateOmegaPerExon(File BEDFile, File modelFile, String alnDir, String alnFileFormat, String save)throws Exception{
		Set<Gene> genes=BEDFileParser.loadData(BEDFile);
		writeOmegaPerExon(genes, modelFile, alnDir, alnFileFormat, save);
	}
	
	public EstimateOmegaPerExon(Collection<Gene> genes, File modelFile, String alnDir, String alnFileFormat, String save)throws Exception{
		writeOmegaPerExon(genes, modelFile, alnDir, alnFileFormat, save);
	}
	
	private void write(Annotation align, OmegaFit fit, FileWriter writer) throws IOException {
		if(fit!=null){writer.write(align+"\t"+fit.getOmega()+"\t"+fit.getLogOddsScore()+"\t"+fit.getPVal()+"\t"+fit.getTreeLength()+"\n");}
		else{System.err.println("NULL: "+align);}
	}
	
	private void write(String save, Map<Alignments, OmegaFit> omega)throws IOException{
		FileWriter writer=new FileWriter(save);
		
		for(Alignments align: omega.keySet()){
			OmegaFit fit=omega.get(align);
			if(fit!=null){writer.write(align+"\t"+fit.getOmega()+"\t"+fit.getLogOddsScore()+"\t"+fit.getPVal()+"\t"+fit.getTreeLength()+"\n");}
			else{System.err.println("NULL: "+align);}
		}
		
		writer.close();
	}
	
	private void writeOmegaPerExon(Collection<Gene> genes, File modelFile, String alnDir, String alnFileFormat, String save) throws Exception {
		FileWriter writer = new FileWriter(save);
		for(Gene gene: genes){
			//System.err.println(gene);
			String alnFile=alnDir+"/"+gene.getChr()+".maf";
			Collection<? extends Annotation> exons=gene.getSortedAndUniqueExons();
			for(Annotation align: exons){
				OmegaFit fit=fitOmega(align, modelFile, alnFile, alnFileFormat);
				if(fit!=null){System.out.println(align+"\t"+fit.getOmega()+"\t"+fit.getLogOddsScore()+"\t"+fit.getPVal()+"\t"+fit.getTreeLength());}
				else{System.err.println("NULL: "+align);}
				write(align, fit, writer);
			}
		}
		writer.close();
	}
	
	private Map<Annotation, OmegaFit> computeOmegaPerExon(Collection<Gene> genes, File modelFile, String alnDir, String alnFileFormat)throws Exception{
		Map<Annotation, OmegaFit> rtrn=new TreeMap<Annotation, OmegaFit>();
		for(Gene gene: genes){
			//System.err.println(gene);
			String alnFile=alnDir+"/"+gene.getChr()+".maf";
			Collection<? extends Annotation> exons=gene.getSortedAndUniqueExons();
			for(Annotation align: exons){
				OmegaFit fit=fitOmega(align, modelFile, alnFile, alnFileFormat);
				rtrn.put(align, fit);
			}
		}
		return rtrn;
	}
	
	private void writeOmegaPerExon(Set<Alignments> set, File modelFile, String alnDir, String alnFileFormat, String save)throws Exception{
		FileWriter writer = new FileWriter(save);
		int i=0;
		for(Alignments align: set){
			//System.err.println(align +" "+i+" "+set.size());	
			String alnFile=alnDir+"/"+align.getChr()+".maf";
			OmegaFit fit=fitOmega(align, modelFile, alnFile, alnFileFormat);
			//if(fit!=null){
			//	System.out.println(align+"\t"+fit.getOmega()+"\t"+fit.getLogOddsScore()+"\t"+fit.getPVal()+"\t"+fit.getTreeLength());
			//}
			if(fit!=null){System.out.println(align+"\t"+fit.getOmega()+"\t"+fit.getLogOddsScore()+"\t"+fit.getPVal()+"\t"+fit.getTreeLength());}
			else{System.err.println("NULL: "+align);}
			write(align, fit, writer);
			i++;
		}
		writer.close();
	}

	
	private static Map computeOmegaPerExon(Set<Alignments> set, File modelFile, String alnDir, String alnFileFormat)throws Exception{
		Map rtrn=new TreeMap();
		int i=0;
		for(Alignments align: set){
			//System.err.println(align +" "+i+" "+set.size());	
			String alnFile=alnDir+"/"+align.getChr()+".maf";
			OmegaFit fit=fitOmega(align, modelFile, alnFile, alnFileFormat);
			//if(fit!=null){
			//	System.out.println(align+"\t"+fit.getOmega()+"\t"+fit.getLogOddsScore()+"\t"+fit.getPVal()+"\t"+fit.getTreeLength());
			//}
			rtrn.put(align, fit);
			i++;
		}
		return rtrn;
	}
	
	public static OmegaFit fitOmega(Annotation align, File modelFile, String alnFile, String alnFileFormat)throws Exception{
		TreeScaler scaler=new TreeScaler(align, modelFile, alnFile, alnFileFormat);
		OmegaFit regionFit = scaler.scaleRegion(align);
		return regionFit;
	}
	
	public static ArrayList<OmegaFit> slideWindowComputeOmega(Annotation align, File modelFile, String alnFile, String alnFileFormat, int windowSize, int overlap)throws Exception{
		TreeScaler scaler=new TreeScaler(align, modelFile, alnFile, alnFileFormat);
		ArrayList<OmegaFit> regions=scaler.scaleTree(windowSize, new ArrayList(), overlap);
		return regions;
	}
	
	
	public static OmegaFit fitOmega(Alignments align, TreeScaler scaler)throws Exception{
		OmegaFit regionFit = scaler.scaleRegion(align);
		return regionFit;
	}
	
	public static ArrayList<OmegaFit> slideWindowComputeOmega(Alignments align, TreeScaler scaler, int windowSize, int overlap)throws Exception{
		ArrayList<OmegaFit> regions=scaler.scaleTree(windowSize, new ArrayList(), overlap);
		return regions;
	}
	
	public static void main(String[] args)throws Exception{
		
		
		if(args.length>4){
			File bedFile=new File(args[0]);
			File modelFile=new File(args[1]);
			String alnDir=args[2];
			String alnFormat=args[3];
			String save=args[4];
			EstimateOmegaPerExon estimateOmegaPerExon = new EstimateOmegaPerExon(bedFile, modelFile, alnDir, alnFormat, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage="v1.2 args[0]=BED file \n args[1]=model file \n args[2]=alignment directory \n args[3]=alignment file format \n args[4]=save file";
	
}
