package broad.pda.seq.chip;

import java.io.*;
import java.util.*;

import nextgen.core.annotation.Annotation;

import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;
import broad.pda.seq.segmentation.AlignmentDataModel;
import broad.pda.seq.segmentation.ContinuousDataAlignmentModel;
import broad.pda.seq.segmentation.GenericAlignmentDataModel;

public class ComputeEnrichments {

	public ComputeEnrichments(File[] files, Collection<Annotation> regions, String save, String sizes)throws IOException{
		Map<Annotation, double[]>[] scoreArray=new Map[files.length];
		ArrayList<String> names=new ArrayList<String>();
		
		for(int i=0; i<files.length; i++){
			if(!files[i].getName().endsWith("sai")){
				System.err.println(files[i]);
				AlignmentDataModel data=new GenericAlignmentDataModel(files[i].getAbsolutePath(), sizes);
				ContinuousDataAlignmentModel model=new ContinuousDataAlignmentModel(data);
				Map<Annotation, double[]> scores=model.scoreSegments(regions);
				scoreArray[i]=scores;
				names.add(files[i].getName());
			}
		}
		write(save, scoreArray, names);
	}
	
	private void write(String save, Map<Annotation, double[]> [] scores, ArrayList<String> names)throws IOException{
		FileWriter writerP=new FileWriter(save+".pvalue");
		FileWriter writerE=new FileWriter(save+".enrichment");
		
		writerE.write("#1.2\n");
		writerE.write(scores[0].keySet().size()+"\t"+names.size()+"\n");
		writerE.write("Region\tRegion");
		for(String name: names){writerE.write("\t"+name);}
		writerE.write("\n");
		
		writerP.write("Region");
		for(String name: names){writerP.write("\t"+name);}
		writerP.write("\n");
		
		for(Annotation align: scores[0].keySet()){
			writerE.write(align.toUCSC()+"\t"+align.toUCSC());
			writerP.write(align.toUCSC());
			for(int i=0; i<scores.length; i++){
				try{
				writerP.write("\t"+scores[i].get(align)[0]);
				writerE.write("\t"+scores[i].get(align)[1]);
				}catch(NullPointerException ex){}
			}
			writerP.write("\n");
			writerE.write("\n");
		}
		
		writerE.close();
		writerP.close();
	}
	
	public static void main(String[] args)throws IOException{
		if(args.length>3){
			File[] files=new File(args[0]).listFiles();
			Collection<Annotation> regions=BEDFileParser.loadAlignmentData(new File(args[1]));
			String save=args[2];
			String sizes=args[3];
			new ComputeEnrichments(files, regions, save, sizes);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=alignment files (dir) \n args[1]=regions \n args[2]=save file \n args[3]=sizes";
}
