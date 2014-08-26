package broad.pda.seq.segmentation;

import java.io.*;
import java.util.*;

import nextgen.core.annotation.AbstractAnnotation;
import nextgen.core.annotation.Annotation;

import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.GeneTools;

public class GlobalSpliceFilter {

	//use the rate of spliced reads globally to identify outliers on the negative tail of the distribution
	//lambda estimation should come from the number of spliced reads dividied by the number of spliced locations
	//also estimate locally, if passes both cutoffs then get rid of it
	
	//only filter if there are alternative isoforms in the graph
	//if only one path exists then keep it
	
	public GlobalSpliceFilter(AlignmentDataModel data, String save, String genomeDirectory, String chr)throws IOException{
		Map<Annotation, double[]> map=filterSpliceSites(data, genomeDirectory, chr);
		
		write(save, map);
	}

	private void write(String save, Map<Annotation, double[]> map) throws IOException {
		FileWriter writer=new FileWriter(save);
				
		for(Annotation intron: map.keySet()){
			double[] vals=map.get(intron);
			writer.write(intron+"\t"+intron.getOrientation()+"\t"+vals[0]+"\t"+vals[1]+"\n");
		}
		
		writer.close();
	}
	
	private Sequence getChrSeq(String genomeDir, String chr) throws IOException{
		String seqFile=genomeDir+"/"+chr.replaceAll("chr", "")+"/"+chr+".fa";
		FastaSequenceIO fsio = new FastaSequenceIO(seqFile);
		List<Sequence> seqs = fsio.extractRecordsWithIDLike(chr, false);
		if(!seqs.isEmpty()) {
			return seqs.get(0);
		}
		else{return null;}
	}

	private Map<Annotation, double[]> filterSpliceSites(AlignmentDataModel data, String genomeDirectory, String chr) throws IOException {
		Map<Annotation, double[]> rtrn=new TreeMap<Annotation, double[]>();
		
		//for(String chr: data.getChromosomeLengths().keySet()){
			System.err.println(chr);
			Alignments chrRegion=new Alignments(chr, 0, data.getChromosomeLengths().get(chr));
			Map<Annotation, Integer> spliced=data.getSplicedReads(chrRegion);
			if(!spliced.isEmpty()){
				double lambda=count(spliced)/spliced.size();
				System.err.println(chr+" "+lambda);
				rtrn.putAll(score(spliced, lambda, genomeDirectory, chr));
			}
		//}
		
		return rtrn;
	}
	
	
	
	private double count(Map<Annotation, Integer> spliced) {
		double num=0;
		
		for(Annotation intron: spliced.keySet()){
			num+=spliced.get(intron);
		}
		
		return num;
	}

	private Map<Annotation, double[]> score(Map<Annotation, Integer> spliced, double lambda, String genomeDir, String chr) throws IOException {
		Map<Annotation, double[]> rtrn=new TreeMap();
		
		for(Annotation intron: spliced.keySet()){
			Sequence chrSeq=getChrSeq(genomeDir, chr);
			String orientation=GeneTools.orientationFromSpliceSites(intron, chrSeq);
			intron.setOrientation(AbstractAnnotation.getStrand(orientation));
			//double cdf=poissonCDF(spliced.get(intron), lambda);
			double[] array={1, spliced.get(intron)};
			rtrn.put(intron, array);
		}
	
		return rtrn;
	}
	
	private double poissonCDF(int k, double lambda){
		cern.jet.random.Poisson poiss=new cern.jet.random.Poisson(lambda, new cern.jet.random.engine.DRand());
		return poiss.cdf(k);
	}

	public static void main(String[] args)throws Exception{
		if(args.length>4){
			AlignmentDataModel data=new GenericAlignmentDataModel(args[0], args[1]);
			String save=args[2];
			String genomeDirectory=args[3];
			String chr=args[4];
			new GlobalSpliceFilter(data, save, genomeDirectory, chr);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=alignment file (SAM) \n args[1]=sizes \n args[2]=save \n args[3]=genomeDirectory \n args[4]=chr";
	
	
}
