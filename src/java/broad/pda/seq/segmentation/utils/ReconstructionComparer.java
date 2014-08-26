package broad.pda.seq.segmentation.utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.*;

import nextgen.core.annotation.Gene;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;

import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.annotation.BEDFileParser;
import broad.pda.annotation.GTFFileParser;
import broad.pda.gene.GeneWithIsoforms;
/**
 * @author skadri
 *
 */
public class ReconstructionComparer{

	
	static final String usage = "Usage: ReconstructionComparer -task <task name> "+
			"\n\tcompare: Compares the reference with another annotation file and outputs a table with results" +
			"\n\t\t-refSet <RefSeq bed file for species 1 [BED by default]> "+
			"\n\t\t-compSet <RefSeq bed file for species 2 [BED by default]> "+
			"\n\t\t-minOverlap MinimumOverlap required for comparison. [Default is 75%] "+
			"\n\t\t-out <Output file [Defaults to stdout]> "
			;
	
	BEDFileParser reference;
	BEDFileParser annotations;
	static Logger logger = Logger.getLogger(ReconstructionComparer.class.getName());
	double minOverlap;
	
	public ReconstructionComparer (String[] args)throws IOException{
	
		Globals.setHeadless(true);
		/*
		 * @param for ArgumentMap - size, usage, default task
		 * argMap maps the command line arguments to the respective parameters
		 */
		ArgumentMap argMap = CLUtil.getParameters(args,usage,"compare");
		
		/*
		 * Read the annotation file for first species
		 */
		String annotations1File = argMap.getMandatory("refSet");
		/*
		 * Read the annotation file for second species
		 */
		String annotations2File = argMap.getMandatory("compSet");
		/*
		 * Check the format of the annotation files and call the GTF or BED parser accordingly
		 */
		reference =  annotations1File.endsWith(".gtf") || annotations1File.endsWith(".GTF")? new GTFFileParser(annotations1File) : new BEDFileParser(annotations1File);
		annotations =  annotations2File.endsWith(".gtf") || annotations2File.endsWith(".GTF")? new GTFFileParser(annotations2File) : new BEDFileParser(annotations2File);

		/*  
		 * This is actually a bidirectional map but not implemented as such.
		 */
		Map<String,String> orthologs = Collections.synchronizedMap(new HashMap<String,String>());
		
		minOverlap = argMap.isPresent("minOverlap")? argMap.getDouble("minOverlap") : 0.75;
		
		/*
		 * Output Filename
		 */
		BufferedWriter bw = argMap.getOutputWriter();
		
		//BufferedWriter bwBed = new BufferedWriter(new FileWriter(argMap.getOutput()+".bed"));
		BufferedWriter bwSum = new BufferedWriter(new FileWriter(argMap.getOutput()+".summary.txt"));
		bwSum.write("geneName"+"\t"+"StartMisannotated"+"\t"+"EndMisannotated"+"\n");

		List<Gene> refGenes = reference.GetGenes();
		List<Gene> compareGenes = annotations.GetGenes();
		
		//For every reference gene
		for(Gene gene1:refGenes){
			for(Gene gene2:compareGenes){
				double orientedStartDiff;
				double orientedEndDiff;
				//if(passesOverlap(gene1,gene2)){
				if(gene2.getName().contains(gene1.getName()) && gene1.getChr().equals(gene2.getChr()) && gene1.getNumExons()==gene2.getNumExons()){
					//compare starts and ends
					if(gene1.getOrientation().equals("+")){
						/*
						 * POSITIVE VALUES = EXTENDED
						 * NEGATIVE VALUES = TRIMMED
						 */
						orientedStartDiff = gene1.getStart()-gene2.getStart();
						orientedEndDiff = gene2.getEnd()-gene1.getEnd();
					}
					else{
						orientedEndDiff = gene1.getStart()-gene2.getStart();
						orientedStartDiff = gene2.getEnd()-gene1.getEnd();
					}
					
					bw.write(gene1.getName()+"\t");
					bw.write(gene2.toUCSC()+"\t");
					bw.write(orientedStartDiff+"\t");
					bw.write(orientedEndDiff+"\n");
					
					if(orientedEndDiff==0 && orientedStartDiff==0){
						//nothing has changed
					}
					else{
						bwSum.write(gene1.getName()+"\t");
						if(orientedStartDiff>0){
							bwSum.write("extend"+"\t");
						}
						else if(orientedStartDiff<0){
							bwSum.write("trim"+"\t");
						}
						else{
							bwSum.write("unchanged"+"\t");
						}
						
						if(orientedEndDiff>0){
							bwSum.write("extend"+"\n");
						}
						else if(orientedEndDiff<0){
							bwSum.write("trim"+"\n");
						}
						else{
							bwSum.write("unchanged"+"\n");
						}
					}
				}
			}
		}
		bw.close();
		//bwBed.close();
		bwSum.close();
	}
	
	private boolean passesOverlap(Gene gene1, Gene gene2){
		if(gene1.getOrientation().equals(gene2.getOrientation())){
			
			if(gene1.percentOverlapping(gene2)>minOverlap){
				return true;
			}
			else{
				return false;
			}
		}
		else{
			return false;
		}
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException{

		ReconstructionComparer dummy = new ReconstructionComparer(args);
	}

}
