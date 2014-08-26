package umms.core.scripture.statistics;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;

import net.sf.samtools.util.CloseableIterator;
import umms.core.alignment.AbstractPairedEndAlignment.TranscriptionRead;
import umms.core.alignment.Alignment;
import umms.core.annotation.Annotation;
import umms.core.annotation.BasicAnnotation;
import umms.core.annotation.Gene;
import umms.core.annotation.Annotation.Strand;
import umms.core.coordinatesystem.TranscriptomeSpace;
import umms.core.general.CloseableFilterIterator;
import umms.core.model.AlignmentModel;
import umms.core.readFilters.PairedAndProperFilter;
import umms.core.readFilters.SameOrientationFilter;
import umms.core.readFilters.SplicedReadFilter;
import umms.core.scripture.BuildScriptureCoordinateSpace;
import org.apache.commons.collections15.Predicate;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;

import broad.core.annotation.MaximumContiguousSubsequence;
import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.Pair;
import broad.core.error.ParseException;
import broad.core.math.ScanStatistics;
import broad.core.math.Statistics;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.annotation.BEDFileParser;
import broad.pda.seq.segmentation.AlignmentDataModelStats;

public class ConnectDisconnectedTranscripts {
	
	static Logger logger = Logger.getLogger(ConnectDisconnectedTranscripts.class.getName());
	Map<String,Collection<Gene>> annotations;
	private AlignmentModel model;
	private int constant = 10000;
	private double medianInsertSize=300;
	private static int DEFAULT_INSERT_SIZE = 500;
	private static int DEFAULT_NUM_BINS = 100;
	private String outName;
	static final String usage = "Usage: ConnectDisconnectedTranscripts -task doWork "+
			"\n**************************************************************"+
			"\n\t\tArguments"+
			"\n**************************************************************"+
			"\n\t\t-out <Output file [Defaults to stdout]> "+
			"\n\t\t-in <Reconstruction bed file. [BED by default]> "+
			"\n\t\t-alignment <Alignment bam file with data on which reconstructions were calculated> "+
			"\n\t\t-strand <first/second : mate in the direction of transcription> "+
			"\n\t\t-maxSize <Maximum insert size. Default=500bp> "+
			"\n";
	
	public ConnectDisconnectedTranscripts(String annotationFile,String outputName,File bamFile,TranscriptionRead strand,int maxSize) throws IOException{
		
		annotations= BEDFileParser.loadDataByChr(new File(annotationFile));
		model=new AlignmentModel(bamFile.getAbsolutePath(), null, new ArrayList<Predicate<Alignment>>(),true,strand,false);
		model.addFilter(new PairedAndProperFilter());
		outName = annotationFile;
		//Compute insert size distribution
		logger.info("Compute insert size distribution");
		//Default number of bins =10
		medianInsertSize = model.getReadSizeDistribution(new TranscriptomeSpace(annotations), maxSize, DEFAULT_NUM_BINS).getMedianOfAllDataValues();
		//medianInsertSize = 600;
		logger.info("Median size = "+medianInsertSize);
		doWork(outputName);
	}
	
	
	public void doWork(String outputName) throws IOException{
	
		boolean somethingWasConnected = true;

		Map<String,Collection<Gene>> conn = null;
		int loop=0;
		int counter=0;
		while(somethingWasConnected && loop<10){
			somethingWasConnected =false;
			loop++;
			//Set counter for each loop
			counter=0;
			logger.info("Connected disconnected transcripts: Loop "+loop);
			conn = new HashMap<String,Collection<Gene>>();
			for(String chr:annotations.keySet()){
				conn.put(chr, new TreeSet<Gene>());			
			}
			
			for(String chr:annotations.keySet()){
				//For all genes on this chromosome
				logger.debug("Connecting, Processing "+chr);
				Collection<Gene> newGenes = new TreeSet<Gene>();
				//MAKE AN INTERVAL TREE OF THE GENES on this chr
				IntervalTree<Gene> tree = new IntervalTree<Gene>();
				for(Gene g:annotations.get(chr)){
					conn.get(chr).add(g);
					tree.put(g.getStart(), g.getEnd(), g);
					counter++;
				}
				//For each transcript
				//Iterate over all reconstructions
				Iterator<Gene> iter=tree.toCollection().iterator();
				while(iter.hasNext()){
					Gene gene=iter.next();
					//For all assemblies downstream of this assembly in 10kB regions
					Iterator<Gene> overlappers=tree.overlappingValueIterator(gene.getEnd(), gene.getEnd()+constant);
					
					newGenes.add(gene);
					while(overlappers.hasNext()){
						Gene other = overlappers.next();
						if(isCandidate(gene,other)){
							if(pairedEndReadSpansTranscripts(gene, other)){ 
								//if(secondTranscriptIsSingleExon(gene,other)){
									logger.debug("Attempt to connect "+gene.getName()+" "+gene.toUCSC()+" and "+other.getName()+" "+other.toUCSC());
									
									//Connect the genes
									Annotation connected = getConnectedTranscript(gene,other);
									if(connected!=null){
										somethingWasConnected = true;
										Gene newConnected = new Gene(connected);
										double[] scores = getScores(newConnected);
										double[] fields = new double[4];
										newConnected.setName(gene.getName()+"_"+other.getName());
										//[0] : sum
										fields[0] = scores[0];
										//[1] : p-value
										fields[1] = scores[1];
										//[2] : FPK
										fields[2] = (scores[0]*1000.0)/newConnected.getSize();
										//[3] : FPKM
										//Calculate FPKM
										fields[3] = fields[2]*((double)1000000.0)/model.getGlobalPairedFragments();
										logger.debug("For isoform : "+newConnected.getName()+"\tNum of exons: "+newConnected.getSpliceConnections().size()+"\t"+fields[0]+"\t"+fields[1]);
										newConnected.setBedScore(fields[3]);
										logger.debug(newConnected.toBED());
										newConnected.setExtraFields(fields);
										newGenes.remove(gene);
										newGenes.add(newConnected);
										boolean there = conn.get(chr).remove(gene);
										if(there){counter--;}
										there = conn.get(chr).remove(other);
										if(there){counter--;}
										conn.get(chr).add(newConnected);
										counter++;
									}
								//}
							}
							else{
								//logger.info("The genes "+gene.getName()+" "+gene.toUCSC()+" and "+other.getName()+" "+other.toUCSC()+" do not have paired reads");
							}
						}
					}				
				}
				logger.info(chr+"\t"+counter);
			}
			annotations = conn;
		}
		
		//FileWriter bw = new FileWriter(outName+".12connected.bed");
		FileWriter bw = new FileWriter(outName+".connected.bed");
		for(String name:conn.keySet()){
			Iterator<Gene> ter = conn.get(name).iterator();
			while(ter.hasNext()){
				Gene isoform = ter.next();
				Gene iso = new Gene(trimEnds(isoform,0.1));
				bw.write(iso.toBED()+"\n");
			}
		}
		bw.close();
	}
	
	/**
	 * This function will return true if gene and other have at least 1 paired end read in common
	 * @param gene
	 * @param other
	 * @return
	 */
	private boolean pairedEndReadSpansTranscripts(Gene gene,Gene other){
		
		boolean rtrn = false;
		//Get all overlapping paired end reads in same orientation as the gene
		CloseableFilterIterator<Alignment> iter = new CloseableFilterIterator<Alignment>(model.getOverlappingReads(new BasicAnnotation(gene.getChr(),gene.getStart(),other.getEnd()),true), new SameOrientationFilter(gene));
		while(iter.hasNext()){
			Alignment read = iter.next();
			List<Annotation> mates = (List<Annotation>) read.getReadAlignments(model.getCoordinateSpace());
			if(mates.get(0).getOrientation().equals(gene.getOrientation()) && mates.get(1).getOrientation().equals(gene.getOrientation())){
				if((gene.contains(mates.get(0)) && other.contains(mates.get(1)))
						||
						gene.contains(mates.get(1)) && other.contains(mates.get(0))){
					rtrn = true;
					break;
				}
			}
			
		}
		iter.close();
		return rtrn;
	}
	
	/**
	 * Returns true if the transcript at the oriented 3' end is a single exon transcript.
	 * @param gene
	 * @param other
	 * @return
	 */
	private boolean secondTranscriptIsSingleExon(Gene gene,Gene other){
		
		Pair<Gene> orderedGenes = ConnectDisconnectedTranscripts.getOrderedAssembly(gene, other);
		if(gene.isNegativeStrand()){
			if(orderedGenes.getValue1().getBlocks().size()==1){
				return true;
			}
		}
		else{
			if(orderedGenes.getValue2().getBlocks().size()==1){
				return true;
			}
		}
		return false;
	}
	
	/**
	 * Connect the two genes by fusing the last exon of the first and the first exon of the last gene
	 * @param gene
	 * @param other
	 * @return
	 * @throws IOException 
	 */
	private Annotation getConnectedTranscript(Gene gene,Gene other) throws IOException{
		
		Pair<Gene> orderedGenes = getOrderedAssembly(gene, other);
		Annotation connected = null;
		/*
		 * If the distance between the transcripts is less than the insert size,
		 * paste through
		 */
		if(orderedGenes.getValue2().getStart()-orderedGenes.getValue1().getEnd() < medianInsertSize){
			Gene firstGene = orderedGenes.getValue1().copy();
			firstGene.setEnd(orderedGenes.getValue2().getStart());
			connected = firstGene.union(orderedGenes.getValue2());
		}
		/*
		 * If the distance between the transcripts is more than the insert size,
		 * 		if there is at least 1 splice read spanning the junction,
		 * 			union it
		 */
		else{
			boolean flag= false;
			Annotation junction = orderedGenes.getValue1().getLastExon().union(orderedGenes.getValue2().getFirstExon());
			CloseableFilterIterator<Alignment> splicedIter=new CloseableFilterIterator<Alignment>(model.getOverlappingReads(junction,false), new SplicedReadFilter());
			while(splicedIter.hasNext()){
				Alignment read = splicedIter.next();
				//if there is at least 1 splice read spanning the junction,
				if(BuildScriptureCoordinateSpace.compatible(read,junction)){
					flag=true;
					break;
				}
			}
			splicedIter.close();
			if(flag){
				connected = gene.union(other);				
			}
			/*
			 * Else NO CONNECTION
			 */
		}
		return connected;		
	}
	
	public static Pair<Gene> getOrderedAssembly(Gene gene1,Gene gene2) {
		Pair<Gene> rtrn=new Pair<Gene>();
		//Order by CompareTo
		if(gene1.compareTo(gene2)<0){
			rtrn.setValue1(gene1);
			rtrn.setValue2(gene2);
		}
		else{
			rtrn.setValue1(gene2);
			rtrn.setValue2(gene1);
		}
		return rtrn;
	}
	/**
	 * This function returns true if
	 * 		1. gene and other do not overlap
	 * 	    2. gene and other are in the same orientation
	 * @param gene
	 * @param other
	 * @return
	 */
	private boolean isCandidate(Gene gene,Gene other){ 
		
		if(gene.overlaps(other)){
			//logger.info(gene.getName()+" overlaps "+gene.toUCSC());
			return false;
		}
		else{
			//if transcript is in an intron of the other transcript
			if(gene.getEnd()>other.getStart() && other.getEnd()>gene.getStart()){
				//logger.info("One is inside the intron of the other");
				return false;
			}
			if(gene.getOrientation().equals(other.getOrientation())){
				//logger.info("The orientation is same");
				return true;
			}
		}
		//logger.info("Something else");
		return false;
	}
	
	/**
	 * Returns paired end counts and scan p-value for the specified gene
	 * @param gene
	 * @return
	 */
	private double[] getScores(Annotation gene){
		double[] scores = new double[2];
		scores[0] = 0.0;
		double globalPairedLambda = model.getGlobalPairedFragments()/model.getGlobalLength();
		//Get all reads overlapping the transcript
		CloseableIterator<Alignment> iter = new CloseableFilterIterator<Alignment>(model.getOverlappingReads(gene,true), new PairedAndProperFilter());
		//For each read,
		while(iter.hasNext()){					
			Alignment read = iter.next();
			boolean countRead = true;
			for(Annotation mate:read.getReadAlignments(model.getCoordinateSpace())){
				if(!BuildScriptureCoordinateSpace.compatible(gene,mate)){
					//logger.debug("Read "+mate.toUCSC()+" is not compatible with isoform with "+isoform.getExons().length);
					countRead=false;
					break;
				}
			}
			//For the assembly, we need to treat each read separately	
			if(countRead){
				scores[0] += read.getWeight();
			}
		}
		iter.close();
		//logger.debug("Count = "+scores[0]+" Int version "+new Double(scores[0]).intValue()+" global paired lambda = "+globalPairedLambda+" gene size = "+model.getCoordinateSpace().getSize(gene)+ " or "+gene.size()+" global length = "+model.getGlobalLength()+" global lambda = "+model.getGlobalLambda());
		scores[1] = ScanStatistics.calculatePVal(new Double(scores[0]).intValue(), globalPairedLambda,gene.size(), model.getGlobalLength());
		
//		scores[1] = AlignmentDataModelStats.calculatePVal(new Double(scores[0]).intValue(), model.getGlobalLambda(), model.getCoordinateSpace().getSize(gene), model.getGlobalLength());
		
		return scores;
	}
	public static void main (String [] args) throws ParseException, IOException {
		
		/*
		 * Gives a log4j error. Check later.
		 */
		Globals.setHeadless(true);
		/*
		 * @param for ArgumentMap - size, usage, default task
		 * argMap maps the command line arguments to the respective parameters
		 */
		ArgumentMap argMap = CLUtil.getParameters(args,usage,"doWork");
		TranscriptionRead strand = TranscriptionRead.UNSTRANDED;
		if(argMap.get("strand").equalsIgnoreCase("first")){
			//System.out.println("First read");
			strand = TranscriptionRead.FIRST_OF_PAIR;
		}
		else if(argMap.get("strand").equalsIgnoreCase("second")){
			//System.out.println("Second read");
			strand = TranscriptionRead.SECOND_OF_PAIR;
		}
		else
			System.out.println("no strand");
		
		
		new ConnectDisconnectedTranscripts(argMap.getInput(),argMap.getOutput(),new File(argMap.getMandatory("alignment")),strand,argMap.getInteger("maxSize", DEFAULT_INSERT_SIZE));
		
	}
	
	/**
	 * Returns the gene trimmed at both ends
	 * @param gene
	 * @return
	 */
	private Gene trimEnds(Gene gene,double pct){
		
		List<Double> counts = model.getCountsStrandedPerPosition(gene);
		double[] cntArr = BuildScriptureCoordinateSpace.l2a(counts);

		Collections.sort(counts);		
		double cutoff = Math.max(2, Statistics.quantile(counts, pct));
		int trimStart = MaximumContiguousSubsequence.contiguousStartSubSequenceOverMin(cntArr, cutoff);
		int trimEnd   =  gene.size() - MaximumContiguousSubsequence.contiguousEndSubSequenceOverMin(cntArr, cutoff);
		logger.debug(gene.getName()+" trimStart " + trimStart + " trimEnd " + trimEnd + ", transcript length " + gene.size()+ " with cutoff "+cutoff);

		if(trimStart>trimEnd){
			return gene;
		}
		if(gene.getOrientation().equals(Strand.NEGATIVE)){
			int temp = trimStart;
			trimStart = trimEnd;
			trimEnd = temp;
			logger.debug("Reset trimStart and TrimEnd to  " + trimStart + " - " + trimEnd);
		}	
			
		Gene newGene = gene.copy();
			
		if(((trimEnd+trimStart) < gene.size()) && trimStart<trimEnd && trimStart<gene.size() && trimEnd<gene.size()){
			
			double score = gene.getBedScore();
			String[] extras = gene.getExtraFields();
			newGene = new Gene(gene.trim(trimStart, trimEnd));		
			newGene.setBedScore(score);
			if(extras!=null)
				newGene.setExtraFields(extras);
			logger.debug("trimming ("+trimStart +" - "+ trimEnd+") gene was: " + gene.toBED() + " and now is: \n" +newGene.toBED());
		}
		return newGene;
	}
}
