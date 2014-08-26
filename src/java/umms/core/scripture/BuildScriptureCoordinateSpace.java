package umms.core.scripture;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.commons.collections15.Predicate;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.jgrapht.GraphPath;

import broad.core.annotation.MaximumContiguousSubsequence;
import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.core.datastructures.Pair;
import broad.core.math.ScanStatistics;
import broad.core.math.Statistics;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.datastructures.Alignments;

import net.sf.samtools.util.CloseableIterator;
import umms.core.alignment.Alignment;
import umms.core.alignment.AbstractPairedEndAlignment.TranscriptionRead;
import umms.core.annotation.Annotation;
import umms.core.annotation.Annotation.Strand;
import umms.core.annotation.BasicAnnotation;
import umms.core.annotation.Gene;
import umms.core.coordinatesystem.CoordinateSpace;
import umms.core.coordinatesystem.TranscriptInGenomicSpace;
import umms.core.coordinatesystem.TranscriptomeSpace;
import umms.core.general.CloseableFilterIterator;
import umms.core.model.JCSAlignmentModel;
import umms.core.readFilters.GenomicSpanFilter;
import umms.core.readFilters.PairedEndFilter;
import umms.core.readFilters.ReadsToReconstructFilter;
import umms.core.readFilters.SplicedReadFilter;
import umms.core.readFilters.UniqueMappedReadsFilter;
import umms.core.readers.PairedEndReader.AlignmentType;
import umms.core.scripture.OrientedChromosomeTranscriptGraph.TranscriptGraphEdge;
import umms.core.scripture.statistics.ConnectDisconnectedTranscripts;

public class BuildScriptureCoordinateSpace {

	private static final TranscriptionRead DEFAULT_TXN_READ =  TranscriptionRead.UNSTRANDED;
	static Logger logger = Logger.getLogger(BuildScriptureCoordinateSpace.class.getName());
	private JCSAlignmentModel model;	
	String genomeSeq = null;
	int windowSize=20000000;
	private CoordinateSpace space;
	private Map<String, ChromosomeTranscriptGraph> graphs;
	private static boolean forceStrandSpecificity=true; //TODO This should be passed or at least determined from data
	private static double DEFAULT_MIN_COV_THRESHOLD = 0.2;
	private static double MIN_SPLICE_PERCENT = 0.05;
	private double coveragePercentThreshold = DEFAULT_MIN_COV_THRESHOLD;
	//private static TranscriptionRead DEFAULT_TXN_READ = TranscriptionRead.UNSTRANDED;
	String outName = null;
	private double DEFAULT_ALPHA = 0.01;
	private double alpha=DEFAULT_ALPHA;
	private static double MIN_SPLICE_READS = 3.0;
	private double THRESHOLD_SPURIOUS = 0.95;
	private double minSpliceReads = MIN_SPLICE_READS;
	private double minSplicePercent = MIN_SPLICE_PERCENT;
	int counter = 1000;
	int globalCounter = 1000;
	File bamFileName;
	private int constant = 10000;
//	private double globalPairedLambda=0.0;
	File bamfile;
	double medianInsertSize;
	TranscriptionRead strand;
	Map<Double,Double> fragmentSizeDistribution;
	boolean isSingleEnd;
	
	public BuildScriptureCoordinateSpace(File bamFile,String genomeDir){
		this(bamFile,genomeDir,bamFile.getName()+".reconstructions",true,DEFAULT_TXN_READ,null);
	}
	/**
	 * 
	 * @param bamFile
	 * @param threshold
	 * @param genomeDir
	 * @param outputName
	 * @param forceStrandedness
	 */
	public BuildScriptureCoordinateSpace(File bamFile,String genomeDir,String outputName,boolean forceStrandedness,TranscriptionRead st,ArgumentMap argMap){
			
		isSingleEnd = argMap.isPresent("singleEnd");
		if(isSingleEnd)
			logger.info("Data is singleEnd");
		else
			logger.info("Data is paired end");
		strand = st;
		bamfile=bamFile;
		this.graphs=new TreeMap<String, ChromosomeTranscriptGraph>();
		genomeSeq = genomeDir;
		bamFileName = bamFile;
		forceStrandSpecificity = forceStrandedness;
		outName = outputName;
		//TODO: CHECK IF DATA IS SINGLE END USING THE SAME HEADER
		JCSAlignmentModel libmodel=new JCSAlignmentModel(bamfile.getAbsolutePath(), null, new ArrayList<Predicate<Alignment>>(),!isSingleEnd,strand,false);
		
		model=new JCSAlignmentModel(bamfile.getAbsolutePath(), new TranscriptInGenomicSpace(libmodel.getRefSequenceLengths()), new ArrayList<Predicate<Alignment>>(),!isSingleEnd,strand,false);
		space=model.getCoordinateSpace();

		setThresholds(argMap);
//		globalFragments = calculateGlobalFragments();
		logger.info("Parameters used: " +
				"\nIntron intention Filter: "+THRESHOLD_SPURIOUS+
				"\nPremature Assembly Filter: "+coveragePercentThreshold+
				"\nSplice junction Filter : "+
				"\n\tNumber of spliced reads : "+minSpliceReads+
				"\n\tPercentage of total spliced reads: "+minSplicePercent+
				"\nAlpha for single exon assemblies : "+alpha
				);
		
		assemble(strand);
		
		Map<String, Collection<Gene>> rtrn=getPaths();
		try {
			FileWriter writer=new FileWriter(outName+".08graph.all.paths.bed");
			for(String chr: rtrn.keySet()){
				for(Gene g:rtrn.get(chr)){
					writer.write(g.toBED()+"\n");
				}
			}			
			writer.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		try {
//			FileWriter writer=new FileWriter(outName+".pairedGenes.bed");
//			FileWriter writer2=new FileWriter(outName+".pairedCounts.txt");
			Map<String, Collection<Gene>> genes = setFPKMScores(rtrn);
//			writer.close();
//			writer2.close();

//			write(outName+".10final.paths.bed",genes);
			
			postProcess(genes);
		} catch (IOException e) {
			e.printStackTrace();
		}
	} 
	
	/**
	 * 
	 * @param bamFile
	 * @param threshold
	 * @param genomeDir
	 * @param outputName
	 * @param forceStrandedness
	 * @param chr
	 */
	public BuildScriptureCoordinateSpace(File bamFile,String genomeDir,String outputName,boolean forceStrandedness,TranscriptionRead strand,String chr,ArgumentMap argMap){
			
		isSingleEnd = argMap.isPresent("singleEnd");
		if(isSingleEnd)
			logger.info("Data is singleEnd");
		else
			logger.info("Data is paired end");
		bamfile=bamFile;
		this.graphs=new TreeMap<String, ChromosomeTranscriptGraph>();
		genomeSeq = genomeDir;
		bamFileName = bamFile;
		forceStrandSpecificity = forceStrandedness;
		//this.strand = strand;
		outName = outputName;
		
		setThresholds(argMap);
		//assemble(strand);
		logger.info("Parameters used: " +
				"\nIntron Retention Filter: "+THRESHOLD_SPURIOUS+
				"\nPremature Assembly Filter: "+coveragePercentThreshold+
				"\nSplice junction Filter : "+
				"\n\tNumber of spliced reads : "+minSpliceReads+
				"\n\tPercentage of total spliced reads: "+minSplicePercent+
				"\nAlpha for single exon assemblies : "+alpha
				);
		
		ChromosomeTranscriptGraph graph=assemble(chr,strand);
		graphs.put(chr, graph);
		
		Map<String, Collection<Gene>> rtrn=getPaths();
		try {
			FileWriter writer=new FileWriter(outName+".08graph.all.paths.bed");
			for(Gene g:rtrn.get(chr)){
				writer.write(g.toBED()+"\n");
			}			
			writer.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		try {
//			FileWriter writer=new FileWriter(outName+".pairedGenes.bed");
//			FileWriter writer2=new FileWriter(outName+".pairedCounts.txt");
			Map<String, Collection<Gene>> genes = setFPKMScores(rtrn);
//			writer.close();
//			writer2.close();

			write(outName+".10final.paths.bed",genes);
			
			postProcess(genes);
		} catch (IOException e) {
			e.printStackTrace();
		}
	} 
	
	private void setThresholds(ArgumentMap argMap){ 
		
		if(argMap!=null){
			coveragePercentThreshold = argMap.getDouble("coverage", DEFAULT_MIN_COV_THRESHOLD);
			alpha = argMap.getDouble("alpha", DEFAULT_ALPHA);
			minSpliceReads = argMap.getDouble("minSpliceReads", MIN_SPLICE_READS);
			minSplicePercent = argMap.getDouble("percentSpliceReads", MIN_SPLICE_PERCENT);
		}
	}

	private void assemble(TranscriptionRead strand) {
		//Iterate over all chromosomes
		for(String chr: space.getReferenceNames()){
			logger.info("Reference name: "+chr);
			if(model.getRefSequenceLambda(chr)==0.0){
				logger.info(chr+" is not expressed in the alignment file");
			}
			else{
				ChromosomeTranscriptGraph graph=assemble(chr,strand);
				this.graphs.put(chr, graph);
			}
		}
	}

	private void postProcess(Map<String, Collection<Gene>> oldGenes) throws IOException{
		
		Map<String, Collection<Gene>> newGenes = new HashMap<String, Collection<Gene>>();
		for(String chr:oldGenes.keySet()){
			//For all paths on this chromosome
			
			//MAKE AN INTERVAL TREE OF THE GENES on this chr
			IntervalTree<Gene> tree = new IntervalTree<Gene>();
			for(Gene g:oldGenes.get(chr)){
				tree.put(g.getStart(), g.getEnd(), g);
			}
			/*
			 * MERGE THE GENES
			 */
			//Step 1: Iterate through the working assembly and find paths that overlap
			//if overlaps but is incompatible, see if we can branch
			//iterate
			Iterator<Gene> iter=tree.toCollection().iterator();			

			Collection<Gene> considered = new HashSet<Gene>();
	 		logger.info("Enter merging");
			while(iter.hasNext()){
				Gene branch1=iter.next();
				considered.add(branch1);
				//get overlapping branches
				Iterator<Gene> overlappers=tree.overlappingValueIterator(branch1.getStart(), branch1.getEnd());
				//Collection<Assembly> toRemove=new TreeSet<Assembly>();
				while(overlappers.hasNext()){
					Gene branch2=overlappers.next();
					if(!considered.contains(branch2)){
						if(!branch1.equals(branch2) && compatible(branch1, branch2)){
							logger.debug("Merging: "+ branch1.getName()+ " and "+branch2.getName());
							Collection<Annotation> rtrn=new TreeSet<Annotation>();
							rtrn.addAll(branch1.getBlocks());
							rtrn.addAll(branch2.getBlocks());
							Gene merged=new Gene(rtrn);
							merged.setName(branch1.getName());
							//SET THE SCORE TO THE MAX FPKM
//							merged.setBedScore(Math.max(branch1.getBedScore(),branch2.getBedScore()));
							double[] scores;
							if(!isSingleEnd)
								scores = getScores(merged,fragmentSizeDistribution);
							else
								scores = getScores(merged);
							double[] fields = new double[4];
							//[0] : sum
							fields[0] = scores[0];
							//[1] : p-value
							fields[1] = scores[1];
							//[2] : FPK
							fields[2] = (scores[0]*1000.0)/merged.getSize();
							//[3] : FPKM
							//Calculate FPKM
							fields[3] = fields[2]*((double)1000000.0)/model.getGlobalPairedFragments();
							merged.setBedScore(fields[2]);
							logger.debug(merged.toBED());
							merged.setExtraFields(fields);
							//remove annotation1 and annotation2
							tree.remove(branch1.getStart(), branch1.getEnd(), branch1);
							tree.remove(branch2.getStart(), branch2.getEnd(), branch2);
							//add merged
							tree.put(merged.getStart(), merged.getEnd(), merged);
						}
					}
				}
			}
			
			/*
			 * FPKM THRESHOLD
			 */
			//MAKE A MAP OF GENE TO ISOFORMS
			Map<Gene,Set<Gene>> isoformMap = getIsoformMap(tree.toCollection());		
			Set<Gene> finalSet = new HashSet<Gene>();
			//For each gene
			for(Gene gene:isoformMap.keySet()){
				//For each transcript
				double max = Double.MIN_VALUE;
				for(Gene isoform:isoformMap.get(gene)){					
					double score = isoform.getBedScore();
					if(score>max){
						max = score;
					}
				}
				//CHeck if the FPKM of each is at least 10% of the mac FPKM
				for(Gene isoform:isoformMap.get(gene)){		
					double score = isoform.getBedScore();
					if(score>0.1*max && score>0.0){
						//filtered.remove(isoform.getStart(), isoform.getEnd(), isoform);
						finalSet.add(isoform);
					}
				}
			}
			newGenes.put(chr, finalSet);
		}
	
		//write(outName+".11postprocessed.paths.bed",newGenes);
		write(outName+".scripture.paths.bed",newGenes);
		if(!isSingleEnd){
			connectDisconnectedTranscripts(newGenes);
		}
		else
			logger.info("Data is single end so disconnected transcripts due to coverage drops cannot be connected");
	}
	
	
	private ChromosomeTranscriptGraph assemble(String chr,TranscriptionRead strand) {
		
		//Option 1: Scan through the space and collapse into compatible edges and nodes
		ChromosomeTranscriptGraph graph=assembleDirectly(chr,strand);
		
		return graph;
	}
	
	/**
	 * 
	 * @param fragSizeDist
	 * @param k
	 * @param w
	 * @param d
	 * @return
	 */
	private double calculatePvalue(Map<Double,Double> fragSizeDist,int k,double w,double d){
		//Return 1 if less than 2 counts OR
		//if minimum insert size > annotation.size
		if(k<=2 || w< Statistics.min(fragSizeDist.keySet())){
			return 1.0;
		}
		
		double lambdaTw = 0.0;
		double lambdaW = 0.0; // parameter for Poisson distribution
		for(double f:fragSizeDist.keySet()){
			if(f<=w){
				lambdaW += fragSizeDist.get(f)*(w-f);
				lambdaTw += fragSizeDist.get(f)*(d-w+f);
			}
		}
		
		double a=((k-lambdaW)/k)*(lambdaTw*ScanStatistics.poisson(k-1, lambdaW));     // poisson function = Poisson PDF
		double result=ScanStatistics.Fp(k-1, lambdaW)*Math.exp(-a);					   // Fp = Poisson CDF
		double p=1-result;
		p=Math.abs(p);
		p=Math.min(1, p);
		//p=Math.max(0, p);
		//logger.info("Count = "+k+" Window = "+w+"LambdaW = "+lambdaW+" Pvalue = "+p);
		return p;
	}
	
	/**
	 * Computes the exact distribution of fragment sizes for the specified reconstructions (transcriptome space)
	 * @param cs
	 * @return
	 */
	private Map<Double,Double> calculateReadSizeDistribution(Map<String, Collection<Gene>> annotations) {
		
		Predicate<Alignment> f= new PairedEndFilter();
		model.addFilter(f);
		Collection<Double> allDataValues = new ArrayList<Double>();
		
		//Use the current reconstructions
		Map<String,Collection<Gene>> temp = new HashMap<String,Collection<Gene>>();
		for(String chr:annotations.keySet()){
			if(!annotations.get(chr).isEmpty()){
				temp.put(chr, annotations.get(chr));
			}
		}
		temp = filterSingleIsoforms(temp);
		//Compute the fragment size distribution
		CoordinateSpace cs = new TranscriptomeSpace(temp);
		
		CloseableIterator<Alignment> iter = model.getReadIterator();
		
		Map<Double,Double> dist = new TreeMap<Double,Double>();
		int done = 0;
		while(iter.hasNext()) {
			Alignment align = iter.next();
			done++;
			if(done % 1000000 == 0) logger.info("Got " + done + " records.");
			try {
				for(Integer size : align.getFragmentSize(cs)) {
					double doublesize = size.doubleValue();
					if(dist.containsKey(doublesize))
						dist.put(doublesize, (dist.get(doublesize)+1));
					else
						dist.put(doublesize, 1.0);
					
					allDataValues.add(doublesize);
				}
			} catch (NullPointerException e) {
				// Catch NullPointerException if the fragment is from a chromosome that is not present in the TranscriptomeSpace
				continue;
			}
		}
		iter.close();
		model.removeFilter(f);
		
		for(Double fragsize:dist.keySet())
			dist.put(fragsize, dist.get(fragsize)/model.getGlobalLength());
		
		medianInsertSize = Statistics.median(allDataValues);
		return dist;
	}
	
	/**
	 * Filters out all annotations that have overlapping isoforms and returns only the single isoforms
	 * For genes with multiple isoforms, only keeps one of the isoforms
	 */
	private Map<String,Collection<Gene>> filterSingleIsoforms(Map<String,Collection<Gene>> annotations){
		
		Map<String,Collection<Gene>> rtrn = new HashMap<String,Collection<Gene>>();
		
		for(String chr:annotations.keySet()){
			Collection<Gene> considered = new HashSet<Gene>();
			//MAKE AN INTERVAL TREE OF THE GENES on this chr
			IntervalTree<Gene> tree = new IntervalTree<Gene>();
			for(Gene g:annotations.get(chr)){
				tree.put(g.getStart(), g.getEnd(), g);
			}
			
			for(Gene gene:annotations.get(chr)){
				
				//If gene has not already been processed
				if(!considered.contains(gene)){
					considered.add(gene);
					
					if(!hasAnOverlappingGene(gene, tree, 0.2)){
						tree.remove(gene.getStart(), gene.getEnd(), gene);
					}
				}			
			}
			rtrn.put(chr, tree.toCollection());
		}
		return rtrn;
	}
		
	/**
	 * This function connects disconnected reconstructions using paired end reads.
	 * @param annotations
	 * @throws IOException 
	 */
	private void connectDisconnectedTranscripts(Map<String, Collection<Gene>> annotations) throws IOException{
		
		// Do not attempt connect if the data is not paired end
		if(!(model.getAlignmentType()==AlignmentType.PAIRED_END)){
			return;
		}
//		model=new JCSAlignmentModel(bamfile.getAbsolutePath(), null, new ArrayList<Predicate<Alignment>>(),true,strand,true);
		model.addFilter(new PairedEndFilter());
//		model.addFilter(new IndelFilter());
//		model.addFilter(new GenomicSpanFilter(20000000));
		
		logger.info("Paired end data used to connect disconnected transcripts missed due to drop in coverage");
		//Use only paired end reads
		
		int loop=0;
		boolean somethingWasConnected = true;
		
//		medianInsertSize += model.getReadSizeDistribution(new TranscriptomeSpace(temp), 800, 100).getMedianOfAllDataValues();
		logger.info("Median size = "+medianInsertSize);

		//Will contain the final set of annotations
		Map<String,Collection<Gene>> conn = new HashMap<String,Collection<Gene>>();
		//Set of transcripts just connected
		//Will contain everything initially then only those changed
		Map<String,Collection<Gene>> justConnected = new HashMap<String,Collection<Gene>>();
		//Initiate
		for(String chr:annotations.keySet()){
			conn.put(chr, new TreeSet<Gene>());	
			justConnected.put(chr, new TreeSet<Gene>());
			for(Gene g:annotations.get(chr)){
				conn.get(chr).add(g);
			}
		}
		
		while(somethingWasConnected && loop<10){
			somethingWasConnected =false;
			loop++;
 
			logger.info("Connected disconnected transcripts: Loop "+loop);
			
			//FIRST TIME GO OVER ALL GENES
			if(loop==1){
				//For each chromosome
				for(String chr:annotations.keySet()){
					
					logger.debug("Connecting, Processing "+chr);
//					Collection<Gene> newGenes = new TreeSet<Gene>();
					//MAKE AN INTERVAL TREE OF THE GENES on this chr
					//TODO: This can be optimized
					IntervalTree<Gene> tree = new IntervalTree<Gene>();
					for(Gene g:annotations.get(chr)){
						tree.put(g.getStart(), g.getEnd(), g);
					}
					//For each transcript
					//Iterate over all reconstructions
					Iterator<Gene> iter=tree.toCollection().iterator();
					while(iter.hasNext()){
						Gene gene=iter.next();
						//For all assemblies downstream of this assembly in 10kB regions
						Iterator<Node<Gene>> overlappers=tree.overlappers(gene.getEnd()+1, gene.getEnd()+constant);
						
						while(overlappers.hasNext()){
							Gene other = overlappers.next().getValue();
							if(isCandidate(gene,other)){
								if(pairedEndReadSpansTranscripts(gene, other)){ 
									//if(secondTranscriptIsSingleExon(gene,other)){
										logger.debug("Attempt to connect "+gene.getName()+" "+gene.toUCSC()+" and "+other.getName()+" "+other.toUCSC());
										
										//Connect the genes
										Annotation connected = getConnectedTranscript(gene,other,medianInsertSize);
										if(connected!=null){
											somethingWasConnected = true;
											Gene newConnected = new Gene(connected);
											double[] scores = getScores(newConnected,fragmentSizeDistribution);
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
											newConnected.setBedScore(fields[2]);
											logger.debug(newConnected.toBED());
											newConnected.setExtraFields(fields);
//											newGenes.add(newConnected);
											justConnected.get(chr).add(newConnected);
											conn.get(chr).remove(gene);
											conn.get(chr).remove(other);
											conn.get(chr).add(newConnected);
										}
									//}
								}
								else{
									//logger.info("The genes "+gene.getName()+" "+gene.toUCSC()+" and "+other.getName()+" "+other.toUCSC()+" do not have paired reads");
								}
							}
						}				
					}
					logger.info("Connected "+justConnected.get(chr).size());
				}
				annotations = conn;
			}
			else{
				//For each chromosome
				for(String chr:annotations.keySet()){
					//ONLY PROCESS IF JUST CONNECTED CONTAINS SOMETHING FROM THIS CHROMOSOME
					if(!justConnected.get(chr).isEmpty()){
						//For all just connected on this chromosome
						logger.debug("Connecting, Processing "+chr);
						Collection<Gene> newGenes = new TreeSet<Gene>();
						//MAKE AN INTERVAL TREE OF THE GENES on this chr
						//TODO: This can be optimized
						IntervalTree<Gene> tree = new IntervalTree<Gene>();
						for(Gene g:conn.get(chr)){
							tree.put(g.getStart(), g.getEnd(), g);
						}
						//For each transcript
						//Iterate over all reconstructions
						for(Gene gene:justConnected.get(chr)){
							//For all assemblies downstream of this assembly in 10kB regions
							Iterator<Gene> overlappers=tree.overlappingValueIterator(gene.getEnd(), gene.getEnd()+constant);
							
							while(overlappers.hasNext()){
								Gene other = overlappers.next();
								if(isCandidate(gene,other)){
									if(pairedEndReadSpansTranscripts(gene, other)){ 
										//if(secondTranscriptIsSingleExon(gene,other)){
											logger.debug("Attempt to connect "+gene.getName()+" "+gene.toUCSC()+" and "+other.getName()+" "+other.toUCSC());
											
											//Connect the genes
											Annotation connected = getConnectedTranscript(gene,other,medianInsertSize);
											if(connected!=null){
												somethingWasConnected = true;
												Gene newConnected = new Gene(connected);
												double[] scores = getScores(newConnected,fragmentSizeDistribution);
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
												newGenes.add(newConnected);
												conn.get(chr).remove(gene);
												conn.get(chr).remove(other);
												conn.get(chr).add(newConnected);
											}
										//}
									}
									else{
										//logger.info("The genes "+gene.getName()+" "+gene.toUCSC()+" and "+other.getName()+" "+other.toUCSC()+" do not have paired reads");
									}
								}
							}				
						}
						justConnected.put(chr, newGenes);
					}
					logger.info("Connected "+justConnected.get(chr).size());
				}
				annotations = conn;
			}				
		}
		
		try{write(outName+"."+"connected.bed",conn);}catch(IOException ex){}
		//FileWriter bw = new FileWriter(outName+".12connected.bed");
/*		FileWriter bw = new FileWriter(outName+".connected.bed");
		for(String name:conn.keySet()){
			Iterator<Gene> ter = conn.get(name).iterator();
			while(ter.hasNext()){
				Gene isoform = ter.next();
				Gene iso = new Gene(trimEnds(isoform,0.1));
				bw.write(isoform.toBED()+"\n");
			}
		}
		bw.close();*/
	}
	
	
	/**
	 * Connect the two genes by fusing the last exon of the first and the first exon of the last gene
	 * @param gene
	 * @param other
	 * @return
	 * @throws IOException 
	 */
	private Annotation getConnectedTranscript(Gene gene,Gene other,double medianInsertSize) throws IOException{
		
		Pair<Gene> orderedGenes = ConnectDisconnectedTranscripts.getOrderedAssembly(gene, other);
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
				if(compatible(read,junction)){
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
	
	/**
	 * This function returns true if
	 * 		1. gene and other do not overlap
	 * 	    2. gene and other are in the same orientation
	 * @param gene
	 * @param other
	 * @return
	 */
	private boolean isCandidate(Gene gene,Gene other){ 
		
		//if they dont overlap
		//if transcript is NOT in an intron of the other transcript
		//Same orientation
		if(!gene.overlaps(other) && !(gene.getEnd()>other.getStart() && other.getEnd()>gene.getStart()) && gene.getOrientation().equals(other.getOrientation()))
			return true;
		return false;
	}
	
	/**
	 * This function will return true if gene and other have at least 1 paired end read in common
	 * @param gene
	 * @param other
	 * @return
	 */
/*	private boolean pairedEndReadSpansTranscripts(Gene gene,Gene other){
		
		boolean rtrn = false;
		//Get all overlapping paired end reads in same orientation as the gene
		CloseableFilterIterator<Alignment> iter = new CloseableFilterIterator<Alignment>(model.getOverlappingReads(new BasicAnnotation(gene.getChr(),gene.getStart(),other.getEnd()),false), new SameOrientationFilter(gene));
		while(iter.hasNext()){
			Alignment read = iter.next();
			List<Annotation> mates = (List<Annotation>) read.getReadAlignments(model.getCoordinateSpace());
			if(mates.get(0).getOrientation().equals(gene.getOrientation()) && mates.get(1).getOrientation().equals(gene.getOrientation())){
				if((gene.overlaps(mates.get(0)) && other.overlaps(mates.get(1)))
						||
						gene.overlaps(mates.get(1)) && other.overlaps(mates.get(0))){
					logger.info("Yes");
					rtrn = true;
					break;
				}
			}
			
		}
		iter.close();
		return rtrn;
	}*/
	
	/**
	 * This function will return true if gene and other have at least 1 paired end read in common
	 * @param gene
	 * @param other
	 * @return
	 */
	private boolean pairedEndReadSpansTranscripts(Gene gene,Gene other){
		
		boolean rtrn = false;
		//Get all overlapping paired end reads in same orientation as the gene
		CloseableIterator<Alignment> iter = model.getOverlappingReads(new BasicAnnotation(gene.getChr(),gene.getStart(),other.getEnd(),gene.getOrientation()),false);
		while(iter.hasNext()){
			Alignment read = iter.next();
			if(read.overlaps(gene) && read.overlaps(other)){
					//logger.info(read.toUCSC()+" "+read.getName());
					rtrn = true;
					break;
			}			
		}
		iter.close();
		return rtrn;
	}
	

	public static void write(String save, Map<String, Collection<Gene>> rtrn) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String chr: rtrn.keySet()){
			Collection<Gene> genes=rtrn.get(chr);
			for(Gene gene: genes){
				writer.write(gene.toBED()+"\n");
			}
		}
		
		writer.close();
	}
	
	public static void memoryStats(){
		int mb = 1024*1024;
        
        //Getting the runtime reference from system
        Runtime runtime = Runtime.getRuntime();
         
        logger.debug("##### Heap utilization statistics [MB] #####");
         
        //Print used memory
        logger.debug("Used Memory:"
            + (runtime.totalMemory() - runtime.freeMemory()) / mb);
 
        //Print free memory
        logger.debug("Free Memory:"
            + runtime.freeMemory() / mb);
         
        //Print total available memory
        logger.debug("Total Memory:" + runtime.totalMemory() / mb);
 
        //Print Maximum available memory
        logger.debug("Max Memory:" + runtime.maxMemory() / mb);
	}
	
	/**
	 * Makes assemblies using all reads and subjects all assemblies to a list of filters.
	 * @param chr
	 * @param strand
	 * @return
	 */
	private ChromosomeTranscriptGraph assembleDirectly(String chr,TranscriptionRead strand){
		
		JCSAlignmentModel libmodel=new JCSAlignmentModel(bamfile.getAbsolutePath(), null, new ArrayList<Predicate<Alignment>>(),!isSingleEnd,strand,false);
		
		model=new JCSAlignmentModel(bamfile.getAbsolutePath(), new TranscriptInGenomicSpace(libmodel.getRefSequenceLengths()), new ArrayList<Predicate<Alignment>>(),!isSingleEnd,strand,false);
		model.addFilter(new UniqueMappedReadsFilter());
		model.addFilter(new ReadsToReconstructFilter());
		model.addFilter(new GenomicSpanFilter(20000000));
		this.space=model.getCoordinateSpace();

		long S = System.currentTimeMillis();	
//		logger.info("Assembling spliced reads");
		long start = System.currentTimeMillis();		
		//SPLICED READS
		//CloseableFilterIterator<Alignment> splicedIter=new CloseableFilterIterator<Alignment>(model.getOverlappingReads(chr), new CanonicalSpliceFilter(genomeSeq));
		CloseableFilterIterator<Alignment> splicedIter=new CloseableFilterIterator<Alignment>(model.getOverlappingReads(chr), new SplicedReadFilter());
		
		IntervalTree<Assembly> splicedAssemblies=assembleDirectly(splicedIter,strand);
		long end = System.currentTimeMillis();
		logger.debug("TIME: ASSEMBLE SPLICED: "+(end-start));
		//try{write(splicedAssemblies, outName+"."+chr+"."+"01splicedAssemblies.bed");}catch(IOException ex){}		
		logger.info("Size of spliced assemblies: "+splicedAssemblies.size());
		
		// ALL READS
		CloseableIterator<Alignment> iter=model.getOverlappingReads(chr);
		IntervalTree<Assembly> workingAssemblies=assembleDirectly(iter, splicedAssemblies,strand);
		end = System.currentTimeMillis();
/*		try{write(workingAssemblies, outName+"."+chr+"."+"01directAsseblies.bed");}catch(IOException ex){}
		CloseableIterator<Alignment> iter=model.getOverlappingReads(chr);
		IntervalTree<Assembly> workingAssemblies=assembleAllDirectly(iter, strand);
		try{write(workingAssemblies, outName+"."+chr+"."+"01assemblies.bed");}catch(IOException ex){}	
*/		
		
		/*
		 * ONLY ALL READS
		 */ 
/*		long start = System.currentTimeMillis();
		logger.info("Assembling all reads");
		CloseableIterator<Alignment> iter=model.getOverlappingReads(chr);
		IntervalTree<Assembly> workingAssemblies=assembleDirectly(iter, strand);
		long end = System.currentTimeMillis();
		logger.debug("TIME: ASSEMBLE NON SPLICED: "+(end-start));*/
		logger.info("Size of direct assemblies: "+workingAssemblies.size());		
				
		logger.info("Intron retention filter");
		//REMOVE SPURIOUS
		start = System.currentTimeMillis();
		IntervalTree<Assembly> unspuriousAssemblies = intronRetentionFilter(workingAssemblies,chr);
		end = System.currentTimeMillis();
		logger.debug("TIME: REMOVE SPURIOUS: "+(end-start));
		try{write(unspuriousAssemblies, outName+"."+chr+"."+"03intronRetentionAssemblies.bed");}catch(IOException ex){}
									
		logger.info("Merge assemblies");
		start = System.currentTimeMillis();
		mergeAssembly(unspuriousAssemblies);
		end = System.currentTimeMillis();
		logger.debug("TIME: MERGE: "+(end-start));
		try{write(unspuriousAssemblies, outName+"."+chr+"."+"04mergedAssemblies.bed");}catch(IOException ex){}

		logger.info("Extend compatible assemblies");
		start = System.currentTimeMillis();
		extendAssembly(unspuriousAssemblies);
		end = System.currentTimeMillis();
		logger.debug("TIME: BRANCHING: "+(end-start));
		try{write(unspuriousAssemblies, outName+"."+chr+"."+"05extendedAssemblies.bed");}catch(IOException ex){}
				
		logger.info("Remove premature assemblies");
		start = System.currentTimeMillis();
		// Flag and remove premature
		IntervalTree<Assembly> filteredAssemblies = removePrematureAssemblies(unspuriousAssemblies);
		end = System.currentTimeMillis();
		logger.debug("TIME: FILTER PREMATURE: "+(end-start));
		//IntervalTree<Assembly> filteredAssemblies = removePrematureAssembliesUsingPairedEnds(unspuriousAssemblies);
		try{write(filteredAssemblies, outName+"."+chr+"."+"06filteredAssemblies.bed");}catch(IOException ex){}

		start = System.currentTimeMillis();

		IntervalTree<Assembly> highSplicedAssemblies = removeLowSpliceJunctionAssemblies(filteredAssemblies);
		end = System.currentTimeMillis();
		logger.debug("TIME: REMOVE LOW SPLICED JUNCTIONS: "+(end-start));
		try{write(highSplicedAssemblies, outName+"."+chr+"."+"07highSplicedAssemblies.bed");}catch(IOException ex){}
		
		//Lets split potential preprocessed and mature transcripts and test whether to include certain nodes
		
		start = System.currentTimeMillis();
		//Make graph
		ChromosomeTranscriptGraph graph=makeGraph(highSplicedAssemblies, chr); //TODO Should use the merged set
		end = System.currentTimeMillis();
		logger.debug("TIME: MAKE GRAPH: "+(end-start));
		long E = System.currentTimeMillis();	
		logger.debug("TIME: TOTAL "+(E-S));
		return graph;
	}

	private void extendAssembly(IntervalTree<Assembly> tree) {
		//We have a set of assemblies that are all incompatible
		//We want to link up parts
		logger.debug("Enter extend assembly");
		//iterate through and get overlapping assemblies
		Iterator<Assembly> iter=tree.toCollection().iterator();
		Collection<Assembly> considered = new HashSet<Assembly>();
		
		//IntervalTree<Assembly> currentAssemblies = new IntervalTree<Assembly>();
		
		while(iter.hasNext()){
			//try and merge the non-overlapping portions
			Assembly assembly1=iter.next();
			//currentAssemblies.put(assembly1.getStart(), assembly1.getEnd(), assembly1);
			//boolean a1Changed = false;
			considered.add(assembly1);
			Iterator<Assembly> overlappers=tree.overlappingValueIterator(assembly1.getStart(), assembly1.getEnd());
			while(overlappers.hasNext()){
				Assembly assembly2=overlappers.next();
				//logger.debug("Assembly1: "+assembly1.getName()+" Assembly2: "+assembly2.getName());
				//logger.debug("overlaps "+assembly2.toUCSC());
				if(!considered.contains(assembly2)){
					//try and merge the non-overlapping portions
					Collection<Assembly> merged=branchAssemblies(assembly1, assembly2);
					if(!merged.isEmpty()){
						logger.debug("Branching "+assembly1.getName()+" and "+assembly2.getName());
						
						//remove assembly1
						tree.remove(assembly1.getStart(), assembly1.getEnd(), assembly1);
						//remove assembly2
						tree.remove(assembly2.getStart(), assembly2.getEnd(), assembly2);						
						considered.remove(assembly1);
						considered.remove(assembly2);
						//currentAssemblies.remove(assembly1.getStart(), assembly1.getEnd(), assembly1);
						//currentAssemblies.remove(assembly2.getStart(), assembly2.getEnd(), assembly2);
						//add merged
						//logger.debug("No. of assemblies "+ merged.size());
						for(Assembly merge: merged){
							tree.put(merge.getStart(), merge.getEnd(), merge);
							considered.add(merge);
							//currentAssemblies.put(merge.getStart(), merge.getEnd(), merge);
						}
					}
				}
			}
		}
		//return currentAssemblies;
	}

	private Collection<Assembly> branchAssemblies(Assembly assembly1, Assembly assembly2) {
		Collection<Assembly> rtrn=new TreeSet<Assembly>();
		
		//order the two assemblies by which starts first
		Pair<Assembly> orderedAssembly=getOrderedAssembly(assembly1, assembly2);
		
		//go from first.getFirstExon() to second.getFirstExon
		Annotation portionToConsiderAdding=orderedAssembly.getValue1().intersect(new Alignments(assembly1.getChr(), orderedAssembly.getValue1().getStart(), orderedAssembly.getValue2().getBlocks().iterator().next().getEnd()));
//		Annotation portionToTest = portionToConsiderAdding.minus(orderedAssembly.getValue2());
		//ONLY IF THE REGION OF THE PORTIONTOCONSIDERADDING THAT DOES NOT OVERLAP THE SECOND ASSEMBLY IS NOT SPLICED
//		if(portionToTest.getSpliceConnections().isEmpty()){
			//if this is compatible with second assemebly
			if(compatible(orderedAssembly.getValue2(), portionToConsiderAdding)){
				Assembly merged=mergeToAssembly(portionToConsiderAdding, orderedAssembly.getValue2());
				merged.setName(orderedAssembly.getValue2().getName());
				rtrn.add(merged);
				rtrn.add(orderedAssembly.getValue1());
				//System.out.println(merged);
			}
//		}
		//get first exon of later start site exon
		//trim assembly 2 from start till this point
		//ask if compatible
		//if so, merge trim with rest of later assembly
		return rtrn;
	}
	
	private Assembly branchIntrons(Assembly exon,Assembly intron) {
		Assembly rtrn = null;
		
		//order the two assemblies by which starts first
		Pair<Assembly> orderedAssembly=getOrderedAssembly(exon, intron);
		
		Annotation portionToConsiderAdding;
		//if exon is first
		if(orderedAssembly.getValue1().equals(exon)){
			portionToConsiderAdding=orderedAssembly.getValue1().intersect(new Alignments(exon.getChr(), orderedAssembly.getValue1().getStart(), orderedAssembly.getValue2().getBlocks().iterator().next().getEnd()));
		}
		else{
			//Last exon
			int index = orderedAssembly.getValue1().getBlocks().size();
			Annotation exon2 = orderedAssembly.getValue1().getBlocks().get(index-1);
			portionToConsiderAdding=orderedAssembly.getValue2().intersect(new Alignments(exon.getChr(), exon2.getStart(), orderedAssembly.getValue2().getEnd()));
			if(orderedAssembly.getValue1().contains(portionToConsiderAdding))
				return rtrn;
		}
		if(compatible(portionToConsiderAdding,intron)){
			rtrn = mergeToAssembly(portionToConsiderAdding, intron);
			rtrn.setName("gene_v_"+globalCounter);
			globalCounter++;
			//System.out.println(merged);
		}
		return rtrn;
	}

	/**
	 * Helper function to assembleAllDirectly
	 * @param assembly1
	 * @param assembly2
	 * @return
	 */
/*	private Assembly branchExonWithAssembly(Assembly exon, Assembly assembly) {
		Assembly rtrn = null;
		
		//order the two assemblies by which starts first
		Pair<Assembly> orderedAssembly=getOrderedAssembly(exon, assembly);
		
		Annotation portionToConsiderAdding;
		//if exon is first
		if(orderedAssembly.getValue1().equals(exon)){
			portionToConsiderAdding=orderedAssembly.getValue1().intersect(new Alignments(exon.getChr(), orderedAssembly.getValue1().getStart(), orderedAssembly.getValue2().getBlocks().iterator().next().getEnd()));
		}
		else{
			//Last exon
			int index = orderedAssembly.getValue1().getBlocks().size();
			Annotation exon2 = orderedAssembly.getValue1().getBlocks().get(index-1);
			portionToConsiderAdding=orderedAssembly.getValue2().intersect(new Alignments(exon.getChr(), exon2.getStart(), orderedAssembly.getValue2().getEnd()));
			if(orderedAssembly.getValue1().contains(portionToConsiderAdding))
				return rtrn;
		}
		if(compatible(portionToConsiderAdding,assembly)){
			rtrn = mergeToAssembly(portionToConsiderAdding, assembly);
			rtrn.setName(assembly.getName());
		}
		return rtrn;
	}*/

	private Pair<Assembly> getOrderedAssembly(Assembly assembly1, Assembly assembly2) {
		Pair<Assembly> rtrn=new Pair<Assembly>();
		//Order by CompareTo
		if(assembly1.compareTo(assembly2)<0){
			rtrn.setValue1(assembly1);
			rtrn.setValue2(assembly2);
		}
		else{
			rtrn.setValue1(assembly2);
			rtrn.setValue2(assembly1);
		}
		return rtrn;
	}

	private void write(IntervalTree<Assembly> tree, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		Collection<Assembly> set=tree.toCollection();
		
		for(Assembly align: set){
			writer.write(align+"\n");
		}
		
		writer.close();
	}
	
	private void writeIntrons(IntervalTree<SpliceJunction> tree, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		Collection<SpliceJunction> set=tree.toCollection();
		
		for(SpliceJunction align: set){
			writer.write(align.toAssembly()+"\n");
		}
		
		writer.close();
	}
	
	private void writeExons(IntervalTree<Exon> tree, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		Collection<Exon> set=tree.toCollection();
		
		for(Exon align: set){
			writer.write(align+"\n");
		}
		
		writer.close();
	}
		
	private IntervalTree<Assembly> assembleDirectly(CloseableIterator<Alignment> iter,TranscriptionRead strand) {
		IntervalTree<Assembly> workingAssemblies=new IntervalTree<Assembly>();
		return assembleDirectly(iter, workingAssemblies,strand);
	}


	/**
	 * Makes a graph using the specified assemblies for the specified chrosomome
	 * This function will go through each gene and only add those that pass the p-value threshold
	 * @param workingAssemblies
	 * @param chr
	 * @return
	 */
	private ChromosomeTranscriptGraph makeGraph(IntervalTree<Assembly> workingAssemblies, String chr) {
		ChromosomeTranscriptGraph graph=new ChromosomeTranscriptGraph(chr);
		Iterator<Node<Assembly>> iter=workingAssemblies.iterator();
		//For each assembly node
		while(iter.hasNext()){
			Collection<Assembly> genes=iter.next().getContainedValues();
			//For each gene(assembly)
			for(Annotation gene: genes){
				//Get all the exons of the gene
				List<? extends Annotation> blocks=gene.getBlocks();
				//if the gene has 1 exon
				if(blocks.size()==1){
					double pval = getScores(gene)[1];
					//System.err.println(gene.toUCSC()+" Count: "+score.getCount()+" pval:"+pval+"new count: "+s+"  new p-val "+pval2);
					if(pval<alpha){
						//System.err.println("passes");
						graph.connectVertexToGraph(blocks.get(0));
					}
					else{
						logger.debug(gene.getName()+" is not added to the graph because it does not pass the p-value threshold.");
					}
				}
				else{
			/*		double pval = getScores(gene)[1];
					if(pval<DEFAULT_ALPHA){*/
						graph.addAnnotationToGraph(gene);
			/*		}
					else{
						logger.debug("Gene "+gene.toUCSC()+" is filtered out because it does not meet the significance threshold : "+pval);
					}*/
				}
			}
		}
		return graph;
	}	

	/**
	 * Returns paired end counts and scan p-value for the specified gene
	 * @param gene
	 * @return
	 */
	private double[] getScores(Annotation gene){
		double[] scores = new double[2];
		scores[0] = 0.0;
		//Get all reads overlapping the transcript
//		CloseableIterator<Alignment> iter = new CloseableFilterIterator<Alignment>(model.getOverlappingReads(gene,true), new PairedEndFilter());
		CloseableIterator<Alignment> iter = model.getOverlappingReads(gene,true);
		//For each read,
		while(iter.hasNext()){					
			Alignment read = iter.next();
			boolean countRead = true;
			for(Annotation mate:read.getReadAlignments(space)){
				if(!compatible(gene,mate)){
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
		//logger.info("Count = "+scores[0]+" Int version "+new Double(scores[0]).intValue()+" global paired lambda = "+globalPairedLambda+" gene size = "+model.getCoordinateSpace().getSize(gene)+ " or "+gene.size()+" global length = "+model.getGlobalLength()+" global lambda = "+model.getGlobalLambda());
		//TODO: REPLACE GLOBAL LAMBDA
		scores[1] = ScanStatistics.calculatePVal(new Double(scores[0]).intValue(), model.getGlobalLambda(),gene.size(), model.getGlobalLength());
		//logger.info("Pvalue = "+scores[1]);
//		scores[1] = ScanStatistics.calculatePVal(new Double(scores[0]).intValue(), model.getGlobalLambda(), model.getCoordinateSpace().getSize(gene), model.getGlobalLength());
		
		return scores;
	}
	
	
	/**
	 * Returns paired end counts and scan p-value BASED ON CONTAINMENT for the specified gene
	 * @param gene
	 * @return
	 */
	private double[] getScores(Annotation gene,Map<Double,Double> fragSizeDist){
		double[] scores = new double[2];
		scores[0] = 0.0;
		//Get all reads overlapping the transcript
		CloseableIterator<Alignment> iter = new CloseableFilterIterator<Alignment>(model.getOverlappingReads(gene,true), new PairedEndFilter());
		
		//For each read,
		while(iter.hasNext()){					
			Alignment read = iter.next();
			boolean countRead = true;
			for(Annotation mate:read.getReadAlignments(space)){
				if(!compatible(gene,mate)){
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
		//logger.info("Count = "+scores[0]+" Int version "+new Double(scores[0]).intValue()+" global paired lambda = "+globalPairedLambda+" gene size = "+model.getCoordinateSpace().getSize(gene)+ " or "+gene.size()+" global length = "+model.getGlobalLength()+" global lambda = "+model.getGlobalLambda());
		scores[1] = calculatePvalue(fragSizeDist, new Double(scores[0]).intValue(), gene.size(), model.getGlobalLength());
		
		return scores;
	}
	
	
	private void mergeAssembly(IntervalTree<Assembly> tree) {
		//Step 1: Iterate through the working assembly and find paths that overlap
		Iterator<Assembly> iter=tree.toCollection().iterator();
		//if overlaps but is incompatible, see if we can branch
		//iterate
		Collection<Assembly> considered = new HashSet<Assembly>();
 		logger.debug("Enter merging");
		while(iter.hasNext()){
			Assembly branch1=iter.next();
			considered.add(branch1);
			//get overlapping branches
			Iterator<Assembly> overlappers=tree.overlappingValueIterator(branch1.getStart(), branch1.getEnd());
			//Collection<Assembly> toRemove=new TreeSet<Assembly>();
			while(overlappers.hasNext()){
				Assembly branch2=overlappers.next();
				if(!considered.contains(branch2)){
					if(!branch1.equals(branch2) && compatible(branch1, branch2)){
						logger.debug("Merging: "+ branch1.getName()+ " and "+branch2.getName());
						Assembly merged=mergeToAssembly(branch1, branch2);
						merged.setName(branch1.getName());
						//remove annotation1 and annotation2
						tree.remove(branch1.getStart(), branch1.getEnd(), branch1);
						tree.remove(branch2.getStart(), branch2.getEnd(), branch2);
						//add merged
						tree.put(merged.getStart(), merged.getEnd(), merged);
					}
				}
			}
		}
	}
	
	/**
	 * This function removes any overlapping bad isoforms
	 * @param assemblies
	 * @return
	 */
	private IntervalTree<Assembly> intronRetentionFilter(IntervalTree<Assembly> assemblies,String chr) {
		//Iterate over all assemblies
		Iterator<Assembly> iter=assemblies.toCollection().iterator();
		
/*		double normalizationFactor = (double)spliceCount/(double)(spliceCount+nonSpliceCount);
		logger.error("R = "+normalizationFactor);*/
		
		//Check for coverage and add only those that pass the coverage mark
		IntervalTree<Assembly> currentAssemblies = new IntervalTree<Assembly>();
		//For each assembly
		while(iter.hasNext()){
			Assembly assembly1=iter.next();
			//If assembly is not already marked as spurious
			if(!assembly1.isSpurious()){

				boolean	toRemoveAssembly1 = false;
				//For all overlapping assemblies
				Iterator<Assembly> overlappers=assemblies.overlappingValueIterator(assembly1.getStart(), assembly1.getEnd());
				while(overlappers.hasNext() && !toRemoveAssembly1){
					Assembly assembly2 = overlappers.next();
					if(assembly1.getOrientation().equals(assembly2.getOrientation())){
						//For all exons of this assembly
						for(Annotation intron1: assembly1.getSpliceConnections()){
							if(!toRemoveAssembly1){
								for(Annotation exon2:assembly2.getBlocks()){
									//If intron in assembly1 is contained in exon of assembly2
									if(exon2.contains(intron1)){
										//CHECK FOR COVERAGE
										double intronCount = model.getIntronCounts(intron1);
										double exonCount = (model.getCount(new BasicAnnotation(exon2.getChr(), intron1.getStart(), intron1.getEnd(), exon2.getOrientation()), false)
																	/(double)(intron1.getEnd()-intron1.getStart()));
										double ratioI = intronCount/(intronCount+exonCount);
										double ratioE = exonCount/(intronCount+exonCount);
										//CASE 1 : Assembly with intron is mature and assembly with exon is immature
										if(ratioI>THRESHOLD_SPURIOUS){
											logger.debug("Comparing exon "+ exon2.toUCSC()+" in assembly "+assembly2.getName()+ " to intron "+intron1.toUCSC()+" in assembly "+assembly1.getName());
											assembly2.setSpurious(true);
											logger.debug("CASE 1: assembly "+ assembly2.getName()+" is spurious. Intron count "+intronCount+" > "+exonCount+" Ratio: "+ratioI);
										}
										else{
											//CASE 2: Assembly with exon is real and the other 
											if(ratioE>THRESHOLD_SPURIOUS){
												//logger.warn(assembly.toUCSC()+" is filtered as pre-mature because "+exon1.toUCSC()+" overlaps "+intron2.toUCSC());
												logger.debug("Comparing exon "+ exon2.toUCSC()+" in assembly "+assembly2.getName()+ " to intron "+intron1.toUCSC()+" in assembly "+assembly1.getName());
												assembly1.setSpurious(true);
												toRemoveAssembly1 = true;
												logger.debug("CASE 2: assembly "+ assembly1.getName()+" is spurious.Exon count "+exonCount+" > "+intronCount+" Ratio: "+(ratioE));
												break;
											}
											//else{
												//Both are real
											//}
										}
									}
								}
							}
							else{
								break;
							}
						}
					}
				}
				if(!toRemoveAssembly1){
					currentAssemblies.put(assembly1.getStart(), assembly1.getEnd(), assembly1);
					//for this assembly, get all the 
					//writer.addAlignment(alignment)
				}
			}
		}
		
/*		iter=currentAssemblies.toCollection().iterator();
		Set<Alignment> reads = new HashSet<Alignment>();
		while(iter.hasNext()){
			Assembly assembly = iter.next();
			CloseableIterator<Alignment> riter = model.getOverlappingReads(assembly, true);
			while(riter.hasNext()){
				Alignment read = riter.next();
				//System.out.println(read.getName());
				if(!reads.contains(read)){
					//writer.addAlignment(read.toSAMRecord());
					reads.add(read);
				}
			}
			riter.close();
		}*/
		//writer.close();
		return currentAssemblies;
	}
	
	/**
	 * This function removes any premature assemblies that have coverage conflicts in the exons
	 * @param assemblies
	 * @return
	 */
	private IntervalTree<Assembly> removePrematureAssemblies(IntervalTree<Assembly> assemblies) {
		//Iterate over all assemblies
		Iterator<Assembly> iter=assemblies.toCollection().iterator();
		
		//Check for coverage and add only those that pass the coverage mark
		IntervalTree<Assembly> currentAssemblies = new IntervalTree<Assembly>();

		while(iter.hasNext()){
			Assembly assembly=iter.next();
			boolean	toRemove = false;
			//For all overlapping assemblies
			Iterator<Assembly> overlappers=assemblies.overlappingValueIterator(assembly.getStart(), assembly.getEnd());
			while(overlappers.hasNext()){
				Assembly overlapper = overlappers.next();
				if(assembly.getOrientation().equals(overlapper.getOrientation())){
					
					//Check if the intron"ed" assembly is confident so it can be used
					if(!overlapper.isConfidentIsSet()){
						setConfidence(overlapper);
					}
					if(overlapper.isConfident()){
						//logger.warn("Compare " +assembly.getName()+" and "+overlapper.getName());
						//CHECK IF THE PARTS OF THE ASSEMBLY OTHER THAN THIS EXON OF THE ASSEMBLY) ARE COMPATIBLE
						//For all exons of this assembly
						for(Annotation exon1: assembly.getBlocks()){
							if(!toRemove){
								//If exon overlaps introns of overlapping assembly
								for(Annotation intron2:overlapper.getSpliceConnections()){
									//Flag as premature				
									if(exon1.overlaps(intron2)){
										Annotation p = assembly.minus(exon1);
										Annotation q = overlapper.minus(exon1);
										if(compatible(p, q) || p.length()==0){
											assembly.setPossiblePremature(true);
											//CHECK FOR COVERAGE
											if(coveragePassesCheck(exon1,intron2)){
												//toRemove = false;
											}
											else{
												//logger.warn(assembly.getName()+" is filtered as pre-mature because "+exon1.toUCSC()+" overlaps "+intron2.toUCSC()+" of "+overlapper.getName());
												toRemove = true;
												break;
											}
										}
									}
								}
							}
							else{
								break;
							}
						}
					}
				}
			}
			if(!toRemove){
				currentAssemblies.put(assembly.getStart(), assembly.getEnd(), assembly);
			}
		}
		
		return currentAssemblies;
	}
	
	private void setConfidence(Assembly overlapper){
		//for each intron
		for(Annotation intron2:overlapper.getSpliceConnections()){
			Annotation[] exons = overlapper.getFlankingBlocks(intron2);
			//exon[0] is left of intron
			//coverage to left of left exon
			//Since it is one base, we dont need coverage
			double leftScore = model.getCount(new BasicAnnotation(exons[0].getChr(),exons[0].getEnd()-2,exons[0].getEnd()-1,overlapper.getOrientation()),false);
			//coverage to right of right exon
			double rightScore = model.getCount(new BasicAnnotation(exons[1].getChr(),exons[1].getStart()+1,exons[1].getStart()+2,overlapper.getOrientation()),false);
			
			if(leftScore<rightScore){
				if(leftScore<rightScore*coveragePercentThreshold){ 
					logger.debug(overlapper.getName()+" does not pass confidence test because "+exons[0].toUCSC()+" has score "+leftScore+" compared to "+
								exons[1].toUCSC()+" which has "+rightScore);
					overlapper.setConfident(false);
				}
			}else{
				if(rightScore<leftScore*coveragePercentThreshold){
					logger.debug(overlapper.getName()+" does not pass confidence test because "+exons[1].toUCSC()+" has score "+rightScore+" compared to "+
							exons[0].toUCSC()+" which has "+leftScore);
					overlapper.setConfident(false);
				}
			}
		}
		if(!overlapper.isConfidentIsSet())
			overlapper.setConfident(true);
	}
	/**
	 * This function removes any premature assemblies that have coverage conflicts in the exons
	 * @param assemblies
	 * @return
	 */
/*	private IntervalTree<Assembly> removePrematureAssembliesUsingPairedEnds(IntervalTree<Assembly> assemblies) {
		//Iterate over all assemblies
		Iterator<Assembly> iter=assemblies.toCollection().iterator();
		
		//Check for coverage and add only those that pass the coverage mark
		IntervalTree<Assembly> currentAssemblies = new IntervalTree<Assembly>();

		while(iter.hasNext()){
			Assembly assembly=iter.next();
			boolean	toRemove = false;
			//For all overlapping assemblies
			Iterator<Assembly> overlappers=assemblies.overlappingValueIterator(assembly.getStart(), assembly.getEnd());
			while(overlappers.hasNext()){
				Assembly overlapper = overlappers.next();
				//For all exons of this assembly
				for(Annotation exon1: assembly.getBlocks()){
					if(!toRemove){
						//If exon overlaps introns of overlapping assembly
						for(Annotation intron2:overlapper.getSpliceConnections()){
							//Flag as premature				
							if(exon1.overlaps(intron2)){
								assembly.setPossiblePremature(true);
								//CHECK FOR COVERAGE
								if(exonPassesPairedEndTest(exon1,overlapper,intron2)){
									//toRemove = false;
								}
								else{
									logger.warn(assembly.toUCSC()+" is filtered as pre-mature because "+exon1.toUCSC()+" overlaps "+intron2.toUCSC());
									toRemove = true;
									break;
								}
							}
						}
					}
					else{
						break;
					}
				}
			}
			if(!toRemove){
				currentAssemblies.put(assembly.getStart(), assembly.getEnd(), assembly);
			}
		}
		
		return currentAssemblies;
	}*/
	
	
	/**
	 * Returns an array of size 2 where [0] is the first mate and [1] is the second mate
	 * @param read
	 * @return
	 */
	private Annotation[] getPairedEndReads(Alignment read){
		
		//We implement this specific to our purpose with paired ends
		Annotation[] annArr = new Annotation[2];
		int i=0;
		for(Annotation a: read.getReadAlignments(space)){
			if(i>1){
				logger.error("More than 2 mates are returned");
			}
			annArr[i] = a;
			i++;
		}
		return annArr;
	}
	
	/**
	 * This function returns true is the first mate is in the first exon and the second mate is in the second exon
	 * @param exon1
	 * @param exon2
	 * @param mates
	 * @return
	 */
	private boolean boundaryIsCompatibleWithPairedEnd(Annotation exon1,Annotation exon2,Alignment read){
		
		Annotation[] mates = getPairedEndReads(read);
		if((compatible(exon1,mates[0]) && 
				compatible(exon2,mates[1]))||
					(compatible(exon2,mates[0]) 
							&& compatible(exon1,mates[1])))
			return true;
		else
			return false;
	}
	/**
	 * This function will remove any assemblies that have introns with very low number of splice junctions
	 * @param assemblies
	 * @return
	 */
	private IntervalTree<Assembly> removeLowSpliceJunctionAssemblies(IntervalTree<Assembly> assemblies){
	
		logger.info("Remove low splice junction assemblies");
		//Iterate over all assemblies
		Iterator<Assembly> iter=assemblies.toCollection().iterator();
				
		//Check for coverage and add only those that pass the coverage mark
		IntervalTree<Assembly> currentAssemblies = new IntervalTree<Assembly>();//assemblies;

		//For each assembly
		while(iter.hasNext()){
			Assembly assembly=iter.next();
			
			boolean toRemove=false;
			//If the assembly does not pass the confidence test, remove it
			if(!assembly.isConfidentIsSet()){
				logger.debug("Confidence is set in the splice filter");
				setConfidence(assembly);
			}
/*			if(!assembly.isConfident()){
				logger.debug(assembly.getName()+" removed because does not pass the confidence test.");
			} else{*/
			if(assembly.isConfident()){
				//IF THE GENE HAS ONE INTRON, check min #splice reads
				if(assembly.getSpliceConnections().size()==1){
					for(Annotation intron:assembly.getSpliceConnections()){
						double count=model.getIntronCounts(intron);
						if(count<minSpliceReads){
							logger.debug(assembly.getName()+" has "+count+" spliced reads. Removing.");
							toRemove = true;
						}
					}
					if(!toRemove)
						currentAssemblies.put(assembly.getStart(), assembly.getEnd(), assembly);
				}
				else{
					Map<Annotation,Double> intronToSplicedCountMap = new TreeMap<Annotation,Double>();
					double avgCount = 0.0;
					//boolean toRemove = false;
					//for each intron
					for(Annotation intron:assembly.getSpliceConnections()){
						//find the number of supporting splice junctions
						double count=model.getIntronCounts(intron);											
						avgCount +=count;
						intronToSplicedCountMap.put(intron, count);
					}
					
					//Calculate average
					avgCount = (avgCount / (double)intronToSplicedCountMap.keySet().size());
					toRemove = false;
					//Check for all introns
/*					int cnt=0;
					//First check the first and last intron. If they do not pass the threshold
					//then trim and then go over this again.
					for(Annotation intron:intronToSplicedCountMap.keySet()){
						if(cnt==0 || cnt==intronToSplicedCountMap.keySet().size()-1){
							if(intronToSplicedCountMap.get(intron)<=avgCount*minSplicePercent){
								assembly = new Assembly(trimEnds(assembly,0.25));
							}
						}
						cnt++;
					}
*/					
					for(Annotation intron:assembly.getSpliceConnections()){
						if(intronToSplicedCountMap.get(intron)<=avgCount*minSplicePercent){
							//If  assembly is flagged to be removed because of an intron,
							logger.debug(assembly.getName()+" removed because intron "+intron.toUCSC()+" has coverage "+intronToSplicedCountMap.get(intron)+" compared to "+avgCount);
							toRemove = true; 
							//return true;
						}
					}
					if(!toRemove){
						//Removed by this filter but confidence was set to remove
						/*if(!assembly.isConfident())
							logger.debug(assembly.getName()+" is NOT confident but RETAINED.");*/
						currentAssemblies.put(assembly.getStart(), assembly.getEnd(), assembly);
					}
					else{
						//Assembly flagged to remove
						logger.debug(assembly.getName()+" is confident but flagged to be REMOVED. Check trimming.");
						//trimEnds(assembly);						
					}
				}
			}
			else{
				logger.debug(assembly.getName()+" is not confident.");
			}
		}
		return currentAssemblies;
	}

	/**
	 * Returns the assembly trimmed at both ends
	 * @param assembly
	 * @return
	 */
	private Annotation trimEnds(Assembly assembly,double pct){
		
		List<Double> counts = model.getCountsStrandedPerPosition(assembly);
		double[] cntArr = l2a(counts);

		Collections.sort(counts);		
		double cutoff = Math.max(2, Statistics.quantile(counts, pct));
		int trimStart = MaximumContiguousSubsequence.contiguousStartSubSequenceOverMin(cntArr, cutoff);
		int trimEnd   =  assembly.size() - MaximumContiguousSubsequence.contiguousEndSubSequenceOverMin(cntArr, cutoff);
		logger.debug(assembly.getName()+" trimStart " + trimStart + " trimEnd " + trimEnd + ", transcript length " + assembly.size()+ " with cutoff "+cutoff);

		if(trimStart>trimEnd){
			return assembly;
		}
		if(assembly.getOrientation().equals(Strand.NEGATIVE)){
			int temp = trimStart;
			trimStart = trimEnd;
			trimEnd = temp;
			logger.debug("Reset trimStart and TrimEnd to  " + trimStart + " - " + trimEnd);
		}	
			
		Annotation newAssembly = new Assembly(assembly);
			
		if(((trimEnd+trimStart) < assembly.size()) && trimStart<trimEnd && trimStart<assembly.size() && trimEnd<assembly.size()){
			
			newAssembly = assembly.trim(trimStart, trimEnd);			
			logger.debug("trimming ("+trimStart +" - "+ trimEnd+") gene was: " + assembly.toBED() + " and now is: " +newAssembly.toBED());
		}
		return newAssembly;
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
	/**
	 * Helper function to removePrematureAssemblies
	 * @param exon1
	 * @param intron2
	 * @return
	 */
	private boolean coveragePassesCheck(Annotation exon1, Annotation intron2){
		Annotation overlap = exon1.intersect(intron2);
		//Annotation nonoverlap = exon1.minus(intron2);
		
		if(overlap==null){
			return true;
		}
		//If fully contained
		//TODO: A minimum coverage threshold?
		if(overlap.equals(exon1)){
			//logger.error("Fully Contained");
			return true;
		}
		
		if(overlap.getLengthOnReference()<=3){
			//logger.warn(exon1.toUCSC()+" overlaps "+intron2.toUCSC()+" with <=3");
			return false;
		}
		
		Annotation exonBoundary = null;
		Annotation intronBoundary = null;
		if(intron2.getEnd()>exon1.getEnd()){
			exonBoundary = new BasicAnnotation(exon1.getChr(),intron2.getStart()-2,intron2.getStart()-1,exon1.getOrientation());
			intronBoundary = new BasicAnnotation(intron2.getChr(),intron2.getStart()+1,intron2.getStart()+2,exon1.getOrientation());
		}
		else{
			exonBoundary = new BasicAnnotation(exon1.getChr(),intron2.getEnd()+1,intron2.getEnd()+2,exon1.getOrientation());
			intronBoundary = new BasicAnnotation(intron2.getChr(),intron2.getEnd()-2,intron2.getEnd()-1,exon1.getOrientation());
		}
		double overlapScore = model.getCount(intronBoundary,false);
		double nonoverlapScore = model.getCount(exonBoundary,false);
//		logger.error(overlap.toUCSC()+" overlap coverage: "+overlapScore+" against "+nonoverlapScore);
		if(overlapScore>=nonoverlapScore*coveragePercentThreshold){
			//logger.error(exon1.toUCSC()+" passes coverage test with intron "+ intron2.toUCSC());
			//logger.error(overlap.toUCSC()+" overlap coverage: "+overlapScore+" against "+nonoverlapScore);
			//logger.debug("Intron boundary score for "+intronBoundary.toUCSC()+ " : "+overlapScore+" Exon boundary score for "+exonBoundary.toUCSC()+ " : "+nonoverlapScore);
			return true;
		}
		else{
			//logger.error(exon1.toUCSC()+" does not pass coverage test with intron "+ intron2.toUCSC());
			logger.debug("REMOVED: Intron boundary score for "+intronBoundary.toUCSC()+ " : "+overlapScore+" Exon boundary score for "+exonBoundary.toUCSC()+ " : "+nonoverlapScore);
			return false;
		}		
	}

	/**
	 * 
	 * @param iter
	 * @param workingAssemblies
	 * @param strand
	 * @return
	 */
	private IntervalTree<Assembly> assembleDirectly(CloseableIterator<Alignment> iter, IntervalTree<Assembly> workingAssemblies,TranscriptionRead strand) {
		
		String linc="gene_v2_";
		boolean flagPremature=!workingAssemblies.isEmpty();
		while(iter.hasNext()){
			
			Alignment reads=iter.next();	
			//For the assembly, we need to treat each read separately
			for(Annotation read: reads.getReadAlignments(space)){			
				//Find all compatible assemblies
				Collection<Assembly> compatibleAssemblies = new ArrayList<Assembly>();

				//EACH READ HAS THE FRAGMENT STRAND
				//for each read, get overlapping assemblies
				Iterator<Node<Assembly>> overlappers=workingAssemblies.overlappers(read.getStart(), read.getEnd());
				//if no overlappers add read as assembly
				if(!overlappers.hasNext()){
					//add the read as an annotation
					//Flag this as likely premature
					Assembly readAssembly=new Assembly(read, false);
					readAssembly.setName(linc+globalCounter);
					globalCounter++;
					if(flagPremature){
						readAssembly.setPossiblePremature(true);
					}
					workingAssemblies.put(readAssembly.getStart(), readAssembly.getEnd(), readAssembly);
				}
				else{
					boolean hasCompatible=false;
					while(overlappers.hasNext()){
						Collection<Assembly> assemblies=new TreeSet<Assembly>(overlappers.next().getContainedValues());
						//Find all compatible assemblies
						for(Assembly assembly: assemblies){
							if(compatible(assembly, read)){
								compatibleAssemblies.add(assembly);
								Assembly merged=mergeToAssembly(assembly, read);
								merged.setName(assembly.getName());
								//remove assembly
								workingAssemblies.remove(assembly.getStart(), assembly.getEnd(), assembly);
								//add merged
								workingAssemblies.put(merged.getStart(), merged.getEnd(), merged);
								hasCompatible=true;
							}
						}
						
					}
					
					if(!hasCompatible){
							Assembly readAssembly=new Assembly(read, false);
							readAssembly.setName(linc+globalCounter);
							globalCounter++;
							if(flagPremature){
								readAssembly.setPossiblePremature(true);
							}
							workingAssemblies.put(readAssembly.getStart(), readAssembly.getEnd(), readAssembly);
					}
				}
				
			}
		}		
		iter.close(); //close the iterator
		
		return workingAssemblies;
	}

	
	/**
	 * 
	 * @param iter
	 * @param workingAssemblies
	 * @param strand
	 * @return
	 */
/*	private IntervalTree<Assembly> assembleAllDirectlyOLD(CloseableIterator<Alignment> iter, TranscriptionRead strand) {

		String linc="gene_v_";
		IntervalTree<Assembly> workingAssemblies=new IntervalTree<Assembly>();
		IntervalTree<SpliceJunction> introns = new IntervalTree<SpliceJunction>();
		IntervalTree<Exon> exons = new IntervalTree<Exon>();
		int intronCounter=0;
		int exonCounter=0;
		//Iterate over all reads
		while(iter.hasNext()){
			
			Alignment reads=iter.next();	
			//For the assembly, we need to treat each mate separately
			for(Alignment read: reads.getReadMates()){
				
				//Non-splcied reads
				if(read.getSpliceConnections().isEmpty()){
					
					//Find all compatible EXONS - No assemblies are spliced right now. So only blocks
					Collection<Exon> compatibleExons = new ArrayList<Exon>();
					//for each read, get overlapping assemblies
					Iterator<Exon> overlappers=exons.overlappingValueIterator(read.getStart(), read.getEnd());
					//if no overlappers add read as assembly
					if(!overlappers.hasNext()){
						//add the read as an annotation
						//Flag this as likely premature
						Exon readExon=new Exon(read,exonCounter);
						exonCounter++;
						exons.put(readExon.getStart(), readExon.getEnd(), readExon);
					}
					else{
						boolean hasCompatible=false;
						//Collection<Assembly> assemblies=new TreeSet<Assembly>(overlappers.next().getContainedValues());
						//For all compatible - overlaps in same orientation because NOT spliced reads
						//Find all compatible assemblies
						//for(Assembly assembly: assemblies){
						while(overlappers.hasNext()){
							Exon overlapper = overlappers.next();
							if(overlap(overlapper, read)){
								Assembly merged=mergeToAssembly(overlapper, read);
								Exon mergedExon = new Exon(merged);
								merged.setName(overlapper.getName());
								compatibleExons.add(mergedExon);
								//remove assembly
								exons.remove(overlapper.getStart(), overlapper.getEnd(), overlapper);
								//add merged
								exons.put(mergedExon.getStart(), mergedExon.getEnd(), mergedExon);
								hasCompatible=true;
							}
						}
						if(!hasCompatible){
							Exon readExon = new Exon(read,exonCounter);
							exonCounter++;
							exons.put(readExon.getStart(), readExon.getEnd(), readExon);
						}
						else //if more than 1 compatible assemblies, merge them.
						if(compatibleExons.size()>1){
							Assembly merged = null;
							//Merge them
							for(Exon overlapper:compatibleExons){
								if(merged==null){
									merged = new Assembly(overlapper);
									exons.remove(overlapper.getStart(), overlapper.getEnd(), overlapper);
								}
								else{
									merged = mergeToAssembly(merged,overlapper);
									merged.setName(overlapper.getName());
									exons.remove(overlapper.getStart(), overlapper.getEnd(), overlapper);
								}
							}
							Exon mergedExon = new Exon(merged);
							mergedExon.setName(merged.getName());
							exons.put(mergedExon.getStart(), mergedExon.getEnd(), mergedExon);
						}
					}
				}
				//If spliced, 
				else{
					//get all splice junctions
					Assembly readStructure = new Assembly(read);
					for(Annotation junction:readStructure.getSpliceConnections()){
						junction.setOrientation(reads.getOrientation());
						Iterator<SpliceJunction> overlappers=introns.overlappingValueIterator(junction.getStart(), junction.getEnd());
						Annotation[] fexons = readStructure.getFlankingBlocks(junction);
						//No junctions in this region
						if(!overlappers.hasNext()){
//							logger.info("New junction");
							junction.setName("intron_"+intronCounter);
							intronCounter++;
							SpliceJunction intron = new SpliceJunction(junction,fexons[0].getStart(),fexons[1].getEnd()); 
							introns.put(junction.getStart()-1, junction.getEnd()+1, intron);
						}
						else{
							boolean hasCompatible=false;
							
							//For all compatible - overlaps in same orientation because NOT spliced reads
							//Find all compatible assemblies
							while(overlappers.hasNext()){
								SpliceJunction spliced = overlappers.next();
								if(spliced.isSameAs(junction)){
//									logger.info("Update existing junction");
									//remove assembly
									introns.remove(spliced.getStart(), spliced.getEnd(), spliced);	
									spliced.update(junction,fexons[0].getStart(),fexons[1].getEnd());
									//add merged
									introns.put(spliced.getStart(), spliced.getEnd(), spliced);
									hasCompatible=true;
								}
							}
							//None of the introns are compatible. New Junction
							if(!hasCompatible){
//								logger.info("No other junctions are compatible. New junction");
								junction.setName("intron_"+intronCounter);
								intronCounter++;
//								logger.info(exons[0].toBED());
//								logger.info(exons[1].toBED());
								SpliceJunction intron = new SpliceJunction(junction,fexons[0].getStart(),fexons[1].getEnd()); 
								introns.put(junction.getStart(), junction.getEnd(), intron);
							}
						}
					}
				}
			}
		}
		//close the iterator
		iter.close();		
		
		try{writeExons(exons, outName+"."+"01exons.bed");}catch(IOException ex){}
		try{writeIntrons(introns, outName+"."+"01introns.bed");}catch(IOException ex){}
		
		for(Exon exon:exons.toCollection()){
			workingAssemblies.put(exon.getStart(), exon.getEnd(), exon);
		}
		// NOW MAKE ASSEMBLIES
		//for each splice junction
		Iterator<SpliceJunction> intronIter = introns.toCollection().iterator();
		while(intronIter.hasNext()){
			SpliceJunction intron = intronIter.next();
			if(intron.getCount()>=minSpliceReads){
			logger.info("INTRON: "+intron.toBED());
			Collection<Assembly> considered = new TreeSet<Assembly>();
			Collection<Assembly> intronAssemblies = new TreeSet<Assembly>();
			
			boolean hasOverlap = false;
			Assembly intronAssembly = intron.toAssembly();
			boolean firstExon = true;
			for(Annotation exon:intronAssembly.getBlocks()){
				if(firstExon){
					logger.info("Exon 1: "+exon.toBED());
					firstExon = false;
				}
				else
					logger.info("Exon 2: "+exon.toBED());
				//Get All OVERLAPPING ASSEMBLIES for this exon
				logger.info("Overlaps "+workingAssemblies.numOverlappers(exon.getStart(), exon.getEnd())+" assemblies.");
				Iterator<Assembly> iter2=workingAssemblies.overlappingValueIterator(exon.getStart(),exon.getEnd());
			
				while(iter2.hasNext()){
					hasOverlap=true;
					
					Assembly assembly = iter2.next();
	//				Node<Assembly> assemblies=iter2.next();
	//				for(Assembly assembly:assemblies.getContainedValues()){
						//For all compatible - overlaps in same orientation because NOT spliced reads
						logger.info("OVERLAPS ASSEMBLY: "+assembly.toBED());
						if(!considered.contains(assembly)){
							//1. If assembly contains the exon, then split exon and add new assembly - branch
							//If the exon is equal to one of the exons of the assembly then don't branch - extend as below
							if(exonIsSubsetOfAssembly(assembly,exon)){// && (assembly.getEnd()!=exon.getEnd() && assembly.getStart()!=exon.getStart())){
								
								Assembly merged = branchIntrons(assembly,intronAssembly);
								//TODO: To optimize, done remove assembly. Only add the extended one.
								if(merged!=null){
									logger.info("BRANCH : "+merged.toBED());
									considered.add(merged);
									boolean notMerge = true;
									Collection<Assembly> result = new TreeSet(intronAssemblies);
									for(Assembly ia: result){
										if(compatible(ia,merged)){
											notMerge = false;
											workingAssemblies.remove(ia.getStart(), ia.getEnd(), ia);
											logger.info("Merge "+ia.getName()+" and "+merged.getName());
											Assembly merged2 = mergeToAssembly(ia,merged);
											intronAssemblies.add(merged2);
											intronAssemblies.remove(ia);
											logger.info("BRANCH MERGED: "+merged2.toBED());
											workingAssemblies.put(merged2.getStart(), merged2.getEnd(), merged2);
										}
									}
									if(notMerge){
										workingAssemblies.put(merged.getStart(), merged.getEnd(), merged);
										intronAssemblies.add(merged);
									}
								}
								//Considered only contains the assemblies resulting from the branching.
								considered.add(assembly);
							}
							//2. If assembly partially contains or ends at junction boundary, extend assembly
							else{
								//TODO: To optimize, remove this check
								if(overlap(assembly,exon)){
									Assembly merged=mergeToAssembly(assembly, intronAssembly);
									merged.setName(assembly.getName());
									intronAssemblies.add(merged);
									logger.info("EXTEND: "+merged.toBED());
									//remove assembly
									workingAssemblies.remove(assembly.getStart(), assembly.getEnd(), assembly);
									//add merged
									workingAssemblies.put(merged.getStart(), merged.getEnd(), merged);
									considered.add(merged);
								}
								else{
									logger.info("Non-overlapping assemblies are tested.");
								}
							}
						}
			//		}
				}
			}
			//Making sure both exons don't overlap anything, then make a new assembly, else dont
			if(!hasOverlap){
				logger.info("NEW INTRON ASSEMBLY : "+intronAssembly.toBED());
				intronAssembly.setName(linc+globalCounter);
				intronAssemblies.add(intronAssembly);
				globalCounter++;
				workingAssemblies.put(intronAssembly.getStart(), intronAssembly.getEnd(), intronAssembly);
			}
		}
		}
		return workingAssemblies;
	}
*/	
	
	/**
	 * 
	 * @param iter
	 * @param workingAssemblies
	 * @param strand
	 * @return
	 */
/*	private IntervalTree<Assembly> assembleAllDirectly(CloseableIterator<Alignment> iter, TranscriptionRead strand) {

		String linc="gene_v_";
		IntervalTree<Assembly> workingAssemblies=new IntervalTree<Assembly>();
		IntervalTree<SpliceJunction> introns = new IntervalTree<SpliceJunction>();
		IntervalTree<Exon> exons = new IntervalTree<Exon>();
		int intronCounter=0;
		int exonCounter=0;
		//Iterate over all reads
		while(iter.hasNext()){
			
			Alignment reads=iter.next();	
			//For the assembly, we need to treat each mate separately
			for(Alignment read: reads.getReadMates()){
				
				//Non-splcied reads
				if(read.getSpliceConnections().isEmpty()){
					
					//Find all compatible EXONS - No assemblies are spliced right now. So only blocks
					Collection<Exon> compatibleExons = new ArrayList<Exon>();
					//for each read, get overlapping assemblies
					Iterator<Exon> overlappers=exons.overlappingValueIterator(read.getStart(), read.getEnd());
					//if no overlappers add read as assembly
					if(!overlappers.hasNext()){
						//add the read as an annotation
						//Flag this as likely premature
						Exon readExon=new Exon(read,exonCounter);
						exonCounter++;
						exons.put(readExon.getStart(), readExon.getEnd(), readExon);
					}
					else{
						boolean hasCompatible=false;
						//Collection<Assembly> assemblies=new TreeSet<Assembly>(overlappers.next().getContainedValues());
						//For all compatible - overlaps in same orientation because NOT spliced reads
						//Find all compatible assemblies
						//for(Assembly assembly: assemblies){
						while(overlappers.hasNext()){
							Exon overlapper = overlappers.next();
							if(overlap(overlapper, read)){
								Assembly merged=mergeToAssembly(overlapper, read);
								Exon mergedExon = new Exon(merged);
								merged.setName(overlapper.getName());
								compatibleExons.add(mergedExon);
								//remove assembly
								exons.remove(overlapper.getStart(), overlapper.getEnd(), overlapper);
								//add merged
								exons.put(mergedExon.getStart(), mergedExon.getEnd(), mergedExon);
								hasCompatible=true;
							}
						}
						if(!hasCompatible){
							Exon readExon = new Exon(read,exonCounter);
							exonCounter++;
							exons.put(readExon.getStart(), readExon.getEnd(), readExon);
						}
						else //if more than 1 compatible assemblies, merge them.
						if(compatibleExons.size()>1){
							Assembly merged = null;
							//Merge them
							for(Exon overlapper:compatibleExons){
								if(merged==null){
									merged = new Assembly(overlapper);
									exons.remove(overlapper.getStart(), overlapper.getEnd(), overlapper);
								}
								else{
									merged = mergeToAssembly(merged,overlapper);
									merged.setName(overlapper.getName());
									exons.remove(overlapper.getStart(), overlapper.getEnd(), overlapper);
								}
							}
							Exon mergedExon = new Exon(merged);
							mergedExon.setName(merged.getName());
							exons.put(mergedExon.getStart(), mergedExon.getEnd(), mergedExon);
						}
					}
				}
				//If spliced, 
				else{
					//get all splice junctions
					Assembly readStructure = new Assembly(read);
					for(Annotation junction:readStructure.getSpliceConnections()){
						junction.setOrientation(reads.getOrientation());
						Iterator<SpliceJunction> overlappers=introns.overlappingValueIterator(junction.getStart(), junction.getEnd());
						Annotation[] fexons = readStructure.getFlankingBlocks(junction);
						//No junctions in this region
						if(!overlappers.hasNext()){
//							logger.info("New junction");
							junction.setName("intron_"+intronCounter);
							intronCounter++;
							SpliceJunction intron = new SpliceJunction(junction,fexons[0].getStart(),fexons[1].getEnd()); 
							introns.put(junction.getStart()-1, junction.getEnd()+1, intron);
						}
						else{
							boolean hasCompatible=false;
							
							//For all compatible - overlaps in same orientation because NOT spliced reads
							//Find all compatible assemblies
							while(overlappers.hasNext()){
								SpliceJunction spliced = overlappers.next();
								if(spliced.isSameAs(junction)){
//									logger.info("Update existing junction");
									//remove assembly
									introns.remove(spliced.getStart(), spliced.getEnd(), spliced);	
									spliced.update(junction,fexons[0].getStart(),fexons[1].getEnd());
									//add merged
									introns.put(spliced.getStart(), spliced.getEnd(), spliced);
									hasCompatible=true;
								}
							}
							//None of the introns are compatible. New Junction
							if(!hasCompatible){
//								logger.info("No other junctions are compatible. New junction");
								junction.setName("intron_"+intronCounter);
								intronCounter++;
//								logger.info(exons[0].toBED());
//								logger.info(exons[1].toBED());
								SpliceJunction intron = new SpliceJunction(junction,fexons[0].getStart(),fexons[1].getEnd()); 
								introns.put(junction.getStart(), junction.getEnd(), intron);
							}
						}
					}
				}
			}
		}
		//close the iterator
		iter.close();		
		
		try{writeExons(exons, outName+"."+"01exons.bed");}catch(IOException ex){}
		try{writeIntrons(introns, outName+"."+"01introns.bed");}catch(IOException ex){}
		// NOW MAKE ASSEMBLIES
		//for each splice junction
		Iterator<SpliceJunction> intronIter = introns.toCollection().iterator();
		while(intronIter.hasNext()){
			SpliceJunction intron = intronIter.next();
			if(intron.getCount()>=minSpliceReads){
			logger.info("\nINTRON: "+intron.toBED());
			Collection<Assembly> considered = new TreeSet<Assembly>();
			Collection<Assembly> intronAssemblies = new TreeSet<Assembly>();
			
			boolean hasOverlap = false;
			Assembly intronAssembly = intron.toAssembly();
			boolean firstExon = true;
			for(Annotation exon:intronAssembly.getBlocks()){
				if(firstExon)
					logger.info("Exon 1: "+exon.toBED());					
				else
					logger.info("Exon 2: "+exon.toBED());
				//STEP 1: All OVERLAPPING ASSEMBLIES for this exon
				Iterator<Assembly> iter2=workingAssemblies.overlappingValueIterator(exon.getStart(),exon.getEnd());
			
				while(iter2.hasNext()){
					Assembly assembly = iter2.next();
	//				Node<Assembly> assemblies=iter2.next();
	//				for(Assembly assembly:assemblies.getContainedValues()){
						//For all compatible - overlaps in same orientation because NOT spliced reads
						logger.info("Overlaps ASSEMBLY: "+assembly.getName());
						if(!considered.contains(assembly)){
							//1. If assembly contains the exon, then split exon and add new assembly - branch
							//If the exon is equal to one of the exons of the assembly then don't branch - extend as below
							if(exonIsSubsetOfAssembly(assembly,exon)){// && (assembly.getEnd()!=exon.getEnd() && assembly.getStart()!=exon.getStart())){
								hasOverlap=true;
								Assembly merged = branchIntrons(assembly,intronAssembly);
								//TODO: To optimize, done remove assembly. Only add the extended one.
								if(merged!=null){
									logger.info("BRANCH : "+merged.toBED());
									considered.add(merged);
									boolean notMerge = true;
									Collection<Assembly> result = new TreeSet(intronAssemblies);
									for(Assembly ia: result){
										if(compatible(ia,merged)){
											notMerge = false;
											workingAssemblies.remove(ia.getStart(), ia.getEnd(), ia);
											logger.info("Merge "+ia.getName()+" and "+merged.getName());
											Assembly merged2 = mergeToAssembly(ia,merged);
											intronAssemblies.add(merged2);
											intronAssemblies.remove(ia);
											logger.info("BRANCH MERGED: "+merged2.toBED());
											workingAssemblies.put(merged2.getStart(), merged2.getEnd(), merged2);
										}
									}
									if(notMerge){
										workingAssemblies.put(merged.getStart(), merged.getEnd(), merged);
										intronAssemblies.add(merged);
									}
								}
								//Considered only contains the assemblies resulting from the branching.
								considered.add(assembly);
							}
							//2. If assembly partially contains or ends at junction boundary, extend assembly
							else{
								//TODO: To optimize, remove this check
								if(compatible(assembly, intronAssembly)){
									hasOverlap=true;
									Assembly merged=mergeToAssembly(assembly, intronAssembly);
									merged.setName(assembly.getName());
									intronAssemblies.add(merged);
									logger.info("EXTEND: "+merged.toBED());
									//remove assembly
									workingAssemblies.remove(assembly.getStart(), assembly.getEnd(), assembly);
									//add merged
									workingAssemblies.put(merged.getStart(), merged.getEnd(), merged);
									considered.add(merged);
								}
								else{
									logger.info("Non-overlapping assemblies are tested.");
								}
							}
						}
			//		}
				}
				
				//STEP 2: GET ALL OVERLAPPING EXONS
				Iterator<Exon> iterE= exons.overlappingValueIterator(exon.getStart(),exon.getEnd());
				Collection<Exon> consideredExons = new TreeSet<Exon>();
				
				while(iterE.hasNext()){
					
					Exon unsplicedExon = iterE.next();
					logger.info("Overlaps EXON: "+unsplicedExon.getName());	
					if(!consideredExons.contains(unsplicedExon)){
						
						if(firstExon){
							Assembly intronAss = new Assembly(intronAssembly);
							if(intronAssemblies.isEmpty()|| !exonContainedInAllAssemblies(unsplicedExon,intronAssemblies)){
							if(exonIsSubsetOfAssembly(unsplicedExon,exon)){
								hasOverlap=true;
								Assembly merged = branchIntrons(unsplicedExon,intronAss);
								//TODO: To optimize, done remove assembly. Only add the extended one.
								if(merged!=null){									
									
									logger.info("EXON BRANCH : "+merged.toBED());
									considered.add(merged);
									boolean notMerge = true;
									Collection<Assembly> result = new TreeSet(intronAssemblies);
									for(Assembly ia: result){
										if(compatible(ia,merged)){
											logger.error("MERGING WITH EXON ASSEMBLIES");
											notMerge = false;
											workingAssemblies.remove(ia.getStart(), ia.getEnd(), ia);
											logger.info("Merge "+ia.getName()+" and "+merged.getName());
											Assembly merged2 = mergeToAssembly(ia,merged);
											intronAssemblies.remove(ia);
											intronAssemblies.add(merged2);
											logger.info("EXON BRANCH MERGED: "+merged2.toBED());
											workingAssemblies.put(merged2.getStart(), merged2.getEnd(), merged2);
										}
									}
									if(notMerge){
										workingAssemblies.put(merged.getStart(), merged.getEnd(), merged);
										intronAssemblies.add(merged);
									}
								}
								//Considered only contains the assemblies resulting from the branching.
								consideredExons.add(unsplicedExon);
							}
							else{
								//TODO: To optimize, remove this check
								if(compatible(unsplicedExon, intronAss)){
									hasOverlap=true;
									Assembly merged=mergeToAssembly(unsplicedExon, intronAss);
									merged.setName("gene_v_"+globalCounter);
									globalCounter++;
									intronAssemblies.add(merged);
									logger.info("EXON EXTEND: "+merged.toBED());
									workingAssemblies.remove(intronAss.getStart(),intronAss.getEnd(),intronAss);
									//add merged
									workingAssemblies.put(merged.getStart(), merged.getEnd(), merged);
									considered.add(merged);
									unsplicedExon.setConsideredFlag(true);
								}
								else{
									logger.info("Non-overlapping assemblies are tested.");
								}
							}
							}
							else{
								unsplicedExon.setConsideredFlag(true);
							}
						}
						//SECOND EXON - instead of using the intron exon - use the assemblies put together thus far
						else{
							Collection<Assembly> assemblies = new TreeSet<Assembly>(intronAssemblies);
							for(Assembly intronAss:assemblies){
							
							if(exonIsSubsetOfExon(unsplicedExon,exon)){								
								Assembly merged = branchExonWithAssembly(unsplicedExon,intronAss);
								//TODO: To optimize, done remove assembly. Only add the extended one.
								if(merged!=null){									
									logger.info("EXON BRANCH : "+merged.toBED());
									considered.add(merged);
									boolean notMerge = true;
									Collection<Assembly> result = new TreeSet(intronAssemblies);
									for(Assembly ia: result){
										if(compatible(ia,merged)){
											logger.error("MERGING WITH EXON ASSEMBLIES");
											notMerge = false;
											workingAssemblies.remove(ia.getStart(), ia.getEnd(), ia);
											logger.info("Merge "+ia.getName()+" and "+merged.getName());
											Assembly merged2 = mergeToAssembly(ia,merged);
											intronAssemblies.remove(ia);
											intronAssemblies.add(merged2);
											logger.info("EXON BRANCH MERGED: "+merged2.toBED());
											workingAssemblies.put(merged2.getStart(), merged2.getEnd(), merged2);
										}
									}
									if(notMerge){
										workingAssemblies.put(merged.getStart(), merged.getEnd(), merged);
										intronAssemblies.add(merged);
									}
								}
								//Considered only contains the assemblies resulting from the branching.
								consideredExons.add(unsplicedExon);
							}
							else{
								if(compatible(unsplicedExon,intronAss)){
									Assembly merged=mergeToAssembly(unsplicedExon, intronAss);
									merged.setName(intronAss.getName());
									intronAssemblies.add(merged);
									logger.info("EXON EXTEND: "+merged.toBED());
									workingAssemblies.remove(intronAss.getStart(),intronAss.getEnd(),intronAss);
									//add merged
									workingAssemblies.put(merged.getStart(), merged.getEnd(), merged);
									considered.add(merged);
									unsplicedExon.setConsideredFlag(true);
								}
								else{
									logger.info("Non-overlapping assemblies are tested.");
								}
							}
						}
						}
					}
					// THE ONLY TIME THE EXON WILL BE IN THE "CONSIDERED" LIST IS WHEN THE EXON OVERLAPS THE ENTIRE JUNCTION
					else{
						
					}
				}
				if(firstExon)
					firstExon = false;
			}
			//Making sure both exons don't overlap anything, then make a new assembly, else dont
			if(!hasOverlap){
				logger.info("NEW INTRON ASSEMBLY : "+intronAssembly.toBED());
				intronAssembly.setName(linc+globalCounter);
				intronAssemblies.add(intronAssembly);
				globalCounter++;
				workingAssemblies.put(intronAssembly.getStart(), intronAssembly.getEnd(), intronAssembly);
			}
		}
		}
		
		for(Exon exon:exons.toCollection()){
			if(!exon.isConsidered()){
				workingAssemblies.placeInTree(exon.getStart(), exon.getEnd(), exon);
				logger.info("ADDING EXON "+exon.getName()+" TO ASSEMBLIES");
			}
		}
		return workingAssemblies;
	}
	
	
	private boolean exonIsSubsetOfExon(Exon unsplicedExon, Annotation exon) {
		if(unsplicedExon.contains(exon)){
			if(unsplicedExon.getStart()==exon.getStart())
				return false;
			else
				return true;
		}
		return false;
	}
	private boolean exonContainedInAllAssemblies(Exon unsplicedExon,
			Collection<Assembly> intronAssemblies) {
		for(Assembly assembly:intronAssemblies)
			if(!assembly.contains(unsplicedExon))
				return false;
		return true;
	}*/
	/** 
	 * Returns true is the exon is a subset of the assembly
	 * AND 
	 * the exon does not match an exon of the assembly
	 * @param assembly
	 * @param exon
	 * @return
	 */
/*	private boolean exonIsSubsetOfAssembly(Annotation assembly, Annotation exon) {
		if(assembly.contains(exon)){
			for(Annotation exon1:assembly.getBlocks()){
				if(exon1.overlaps(exon)){
					if(exon1.getEnd()==exon.getEnd())
						return false;
					else
						return true;
				}
			}
		}
		return false;
	}*/

	/*private boolean partiallyCompatible(Annotation assembly, Annotation read) {
		//The annotations are not compatible, so we will check if they are partially compatible
		Collection<? extends Annotation> assemblyIntrons=assembly.getSpliceConnections();
		Collection<? extends Annotation> readIntrons=read.getSpliceConnections();
		
		//If they dont overlap they cant be partially compatible
		if(!assembly.overlaps(read)){return false;}
		
		//Otherwise, they are partially compatible if:
		//(i) some of the introns are compatible
		if(!assemblyIntrons.isEmpty() && !readIntrons.isEmpty()){
			//if there are some introns that
		}
		
		//(ii) or if exons of one overlap the other but dont overlap the intron of the other
		
	}*/

	
	
	//MG: This was working well
	public static boolean compatible(Annotation assembly, Annotation read) {
		//Two alignments will be defined as compatible if:
		//(i) the intronic locations are exactly same
		//if both have introns, ensure that there are no non-overlapping introns
		//logger.debug(read.getName());
		if(!assembly.getSpliceConnections().isEmpty() 
				|| !read.getSpliceConnections().isEmpty()){
			//check if introns are compatible
			if(areIntronsCompatible(assembly, read)){
				//logger.debug("Introns are compatible");
				return true;
			}
			//logger.debug("Introns are NOT compatible");
			//TODO In the case of partial compatibility we should split and make new path
			return false;
		}
		//(ii) one is spliced, the other is not 
		/*if(assembly.getBlocks().size()>1 || !read.getSpliceConnections().isEmpty()){
			//and the non-spliced does not overlap the intron at all, but overlaps the exon
			boolean overlapsWithoutCrossingIntron=overlapsWithoutCrossingIntron(assembly, read);
			if(overlapsWithoutCrossingIntron){return true;}
			return false;
		}*/
		//(iii) both are unspliced and overlap
		boolean overlap=overlap(assembly, read);
		
		if(overlap){return true;}
		return false;
	}
	

	private static boolean areIntronsCompatible(Annotation assembly, Annotation read) {
		Collection<? extends Annotation> assemblyIntrons=assembly.getSpliceConnections();
		Collection<? extends Annotation> readIntrons=read.getSpliceConnections();
		
		//introns are compatible if:
		//(i) all overlapping introns are identical
		for(Annotation intron1: assemblyIntrons){
			for(Annotation intron2: readIntrons){
				//if overlaps but not identical
				if(intron1.overlaps(intron2,forceStrandSpecificity)){
					if(!intron1.equals(intron2, forceStrandSpecificity)){
					//	logger.debug("Case1");
						return false;
					} //TODO: This should use strand info or not based on flag
				}
			}
		}
		
		//(ii) if none of the exons overlap an intron in the other
		for(Annotation exon1: assembly.getBlocks()){
			for(Annotation intron2: readIntrons){
				if(exon1.overlaps(intron2,forceStrandSpecificity)){
					//logger.debug("Case2a");
					return false;
				}
			}
		}
		
		for(Annotation exon2: read.getBlocks()){
			for(Annotation intron1: assemblyIntrons){
				if(exon2.overlaps(intron1,forceStrandSpecificity)){
					//logger.debug("Case2b");
					return false;
				}
			}
		}
		
		//(ii) the introns dont overlap but the exons do
		//just need to test that any exons overlap
		for(Annotation exon1: assembly.getBlocks()){
			for(Annotation exon2: read.getBlocks()){
				if(exon1.overlaps(exon2,forceStrandSpecificity)){
					return true;
				}
			}
		}
		//logger.debug("Case3");
		return false;
	}
	

	private static boolean overlap(Annotation assembly, Annotation read) {
		return assembly.overlaps(read,forceStrandSpecificity);
	}

	private Assembly mergeToAssembly(Annotation assembly, Annotation read) {
		//since these are compatable, we will simply add the read to the assembly by merging exons
		Collection<Annotation> rtrn=new TreeSet<Annotation>();
		rtrn.addAll(read.getBlocks());
		rtrn.addAll(assembly.getBlocks());
		Assembly a = new Assembly(rtrn);
		a.setName(assembly.getName());
		return a;

	}
		


	public Map<String, Collection<Gene>> getPaths() {
		Map<String, Collection<Gene>> rtrn=new TreeMap<String, Collection<Gene>>();
		
		//For each chromosome
		for(String chr: this.graphs.keySet()){
			List<GraphPath<Annotation, TranscriptGraphEdge>> paths=this.graphs.get(chr).getPaths();
			Collection<Gene> genes=new TreeSet<Gene>();
			for(GraphPath<Annotation, TranscriptGraphEdge> path: paths){
				//logger.debug(path.toString());
				Gene gene=this.graphs.get(chr).pathToGene(path);
				genes.add(gene);
			}
			rtrn.put(chr, genes);
		}
		//Add orphan genes
		for(String chr: this.graphs.keySet()){
			for(Gene g:this.graphs.get(chr).getOrphanGenes()){
				rtrn.get(chr).add(g);
			}			
		}
		return rtrn;
	}
	
	/**
	 * Converts a list to an array
	 * @param list
	 * @return
	 */
	public static double[] l2a(List<Double> list){
		double[] rtrn=new double[list.size()];
	
		int i=0;
		for(Double val: list){rtrn[i++]=val;}
	
		return rtrn;
	}
	

	/**
	 * This function will go through each gene and calculate the number of paired end reads 
	 * fully contained within each isoform for the gene
	 * @return
	 * @throws IOException 
	 */
	private Map<String,Collection<Gene>> setFPKMScores(Map<String,Collection<Gene>> geneMap){//,FileWriter writer,FileWriter writer2){
		
		//FIRST CALCULATE THE FRAGMENT SIZE DISTRIBUTION
		//TODO: For genome level calculate once
		if(!isSingleEnd){
			logger.info("Calculating the fragment size distribution");
			fragmentSizeDistribution  = calculateReadSizeDistribution(geneMap);
		}
		Map<String,Collection<Gene>> filteredGenes = new HashMap<String,Collection<Gene>>();
		String name = "gene.v2.1_";
		for(String chr:geneMap.keySet()){
			//logger.debug("For chromosome "+chr+" graph gave "+geneMap.get(chr).size()+" genes");
			Collection<Gene> filtered = new ArrayList<Gene>();
//			try{
			//MAKE A MAP OF GENE TO ISOFORMS
			Map<Gene,Set<Gene>> isoformMap = getIsoformMap(geneMap.get(chr));
			Map<Gene,Double> geneToPairedCount = new HashMap<Gene,Double>();
			
			//For each gene
			for(Gene gene:isoformMap.keySet()){
//				writer2.write("Gene:\n");
				//logger.debug("Starting new gene");
				//For each transcript
				for(Gene isoform:isoformMap.get(gene)){
										
					double[] scores;
					if(!isSingleEnd)
						scores = getScores(isoform,fragmentSizeDistribution);
					else
						scores = getScores(isoform);
					double[] fields = new double[4];
					logger.debug(gene.toUCSC()+" "+scores[0]+"\t"+scores[1]);
					if(scores[1]<alpha){
						isoform.setName(name+new Double(counter).toString()+"_"+isoform.getChr());
						//[0] : sum
						fields[0] = scores[0];
						//[1] : p-value
						fields[1] = scores[1];
						//[2] : FPK
						fields[2] = (scores[0]*1000.0)/isoform.getSize();
						//[3] : FPKM
						//Calculate FPKM
						fields[3] = fields[2]*((double)1000000.0)/model.getGlobalPairedFragments();
						//logger.debug("For isoform : "+isoform.getName()+"\tNum of exons: "+isoform.getSpliceConnections().size()+"\t"+fields[0]+"\t"+fields[1]);
						isoform.setBedScore(fields[2]);
						isoform.setExtraFields(fields);
						counter++;
//						writer.write(isoform+"\n");
//						writer2.write(isoform.toBED()+"\t"+fields[0]+"\n");
						geneToPairedCount.put(isoform, fields[0]);
						filtered.add(isoform);
					}
					else{
						logger.debug("Gene "+gene.toUCSC()+" is filtered out because it does not meet the significance threshold : "+scores[1]);
					}
				}
			}
/*			} catch (IOException e) {
				e.printStackTrace();
			}*/
			logger.info("After significance "+filtered.size()+" genes");
			filteredGenes.put(chr, filtered);
		}
		return filteredGenes;
		
	}

	/**
	 * 
	 * @param genes
	 * @return
	 */
	public static Map<Gene,Set<Gene>> getIsoformMap(Collection<Gene> genes){
				
		Map<Gene,Set<Gene>> isoformMap = new HashMap<Gene,Set<Gene>>();
		
		for(Gene g:genes){
			//If gene has not already been processed
			if(!isoformMapContains(isoformMap,g)){
				//Add gene to its own isoform map
				Set<Gene> set = new HashSet<Gene>();
				set.add(g);
				Collection<Gene> potentialIsoforms = getOverlappingGenes(g,genes,0.5);
				for(Gene p:potentialIsoforms){
					if(!isoformMapContains(isoformMap,p)){
						set.add(p);
					}
				}
				isoformMap.put(g, set);
			}			
		}
		return isoformMap;
	}
	
	public static boolean isoformMapContains(Map<Gene,Set<Gene>> isoformMap,Gene gene){
	
		for(Set<Gene> set:isoformMap.values()){
			for(Gene g:set)
				if(g.equals(gene))
					return true;
		}
		return false;
	}
	/**
	 * 
	 * @param gene
	 * @param allGenes
	 * @param minPctOverlap
	 * @return
	 * TODO: optimize using interval tree
	 */
	public static Collection<Gene> getOverlappingGenes(Gene gene,Collection<Gene> allGenes, double minPctOverlap){
		
		Collection<Gene> overlappers = new HashSet<Gene>();
		for(Gene g:allGenes){
			if(gene.overlaps(g, minPctOverlap) && g.getOrientation().equals(gene.getOrientation()))
				overlappers.add(g);
		}
		return overlappers;
	}
	
	/**
	 * returns a collection of genes in the interval tree that overlap this gene with atleast minPctOverlap and in the same orientation
	 * @param gene
	 * @param allGenesTree
	 * @param minPctOverlap
	 * @return
	 */
	public static Collection<Gene> getOverlappingGenes(Gene gene,IntervalTree<Gene> allGenesTree, double minPctOverlap){
		
		Collection<Gene> overlappers = new HashSet<Gene>();
		Iterator<Gene> OL=allGenesTree.overlappingValueIterator(gene.getStart(), gene.getEnd());
		while(OL.hasNext()){
			Gene gene2=OL.next();
			if(!gene.equals(gene2) && gene.overlaps(gene2, minPctOverlap) && gene2.getOrientation().equals(gene.getOrientation())){
				overlappers.add(gene2);
			}		
		}
		return overlappers;
	}
	
	/**
	 * Returns true if has an overlapping gene with atleast minPctOverlap and in the same orientationin the interval tree
	 * @param gene
	 * @param allGenesTree
	 * @param minPctOverlap
	 * @return
	 */
	public static boolean hasAnOverlappingGene(Gene gene,IntervalTree<Gene> allGenesTree, double minPctOverlap){
		
		Iterator<Gene> OL=allGenesTree.overlappingValueIterator(gene.getStart(), gene.getEnd());
		while(OL.hasNext()){
			Gene gene2=OL.next();
			if(gene.overlaps(gene2, minPctOverlap) && gene2.getOrientation().equals(gene.getOrientation())){
				return true;
			}		
		}
		return false;
	}
	
	public static void main(String[] args)throws IOException{
		
		Globals.setHeadless(true);
		/*
		 * @param for ArgumentMap - size, usage, default task
		 * argMap maps the command line arguments to the respective parameters
		 */
		ArgumentMap argMap = CLUtil.getParameters(args,usage,"reconstruct");
		
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
			logger.info("no strand");
		
		if(argMap.containsKey("chr")){
			new BuildScriptureCoordinateSpace(new File(argMap.getMandatory("alignment")),argMap.getMandatory("genome"),argMap.getOutput(),true, strand,argMap.getMandatory("chr"),argMap);
		}
		else{
			new BuildScriptureCoordinateSpace(new File(argMap.getMandatory("alignment")),argMap.getMandatory("genome"),argMap.getOutput(),true, strand,argMap);
		}
/*		if(args.length>1){
			File bamFile=new File(args[0]);
			String genomeSeqFile = null;
			double threshold = new Double(args[1]);
			genomeSeqFile = args[2];
			if(args.length==3)
				new BuildScriptureCoordinateSpace(bamFile,threshold,genomeSeqFile,args[3]);
			if(args.length==3)
				new BuildScriptureCoordinateSpace(bamFile,threshold,null,args[2],true, TranscriptionRead.SECOND_OF_PAIR);
			if(args.length>=5){
				TranscriptionRead strand = TranscriptionRead.UNSTRANDED;
				if(args[4].equalsIgnoreCase("first")){
					//System.out.println("First read");
					strand = TranscriptionRead.FIRST_OF_PAIR;
				}
				else if(args[4].equalsIgnoreCase("second")){
					//System.out.println("Second read");
					strand = TranscriptionRead.SECOND_OF_PAIR;
				}
				else
					System.out.println("no strand");
				if(args.length==5){
					new BuildScriptureCoordinateSpace(bamFile,threshold,genomeSeqFile,args[3],true, strand);
				}
				else{
					new BuildScriptureCoordinateSpace(bamFile,threshold,genomeSeqFile,args[3],true, strand,args[5]);
				}
			}
		}
		else{System.err.println(usage);}*/
	}
	
	static final String usage = "Usage: BuildScriptureCoordinateSpace -task reconstruct "+
			"\n**************************************************************"+
			"\n\t\tArguments"+
			"\n**************************************************************"+
			"\n\n\t\t-alignment <Alignment file to be used for reconstruction.> "+
			"\n\n\t\t-genome <Fasta file with Genome sequence to be used as reference.> "+
			"\n\t\t-out <Output file [Defaults to stdout]> "+		
			"\n\t\t-strand <VALUES: first, second, unstranded. Specifies the mate that is in the direction of transcription DEFAULT: Unstranded> "+
			
			"\n\n**************************************************************"+
			"\n\t\tOptional Arguments"+
			"\n**************************************************************"+
			"\n\t\t-chr <If specified, Scripture will be run for this chromosome only.> "+
			"\n\t\t-coverage <Specifies the minimum percentage of drop in coverage allowed for an exon. DEFAULT: 0.2> "+
			"\n\t\t-minSpliceReads <The minimum number of splice reads allowed to support a single intron transcript. DEFAULT: 3> "+
			"\n\t\t-percentSpliceReads <The minimum percentage of the average splice counts for a transcript, that an intron can be supported by. DEFAULT: 0.05> "+
			"\n\t\t-alpha <The significance p-value threshold for reconstructions. DEFAULT: 0.01> "+
			"\n";
	
	//static String usage=" args[0]=bam file \n\t args[1]=minimum percentage threshold for coverage \n\t args[2]: Fasta file with Genome sequence"
	//					+"\n\t args[3]= outputName \n\t args[4] transcription strand \n\targs[5] If specified only this chromosome";
	
}
