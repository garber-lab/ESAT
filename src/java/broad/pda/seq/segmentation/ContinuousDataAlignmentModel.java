package broad.pda.seq.segmentation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Stack;
import java.util.TreeMap;
import java.util.TreeSet;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;
import nextgen.core.alignment.Alignment;
import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Gene;
import nextgen.core.annotation.Annotation.Strand;
import nextgen.core.feature.Window;

import org.apache.commons.collections15.Predicate;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;
import org.broad.igv.sam.ReadMate;
import org.jgrapht.util.VertexPair;

import broad.core.annotation.BED;
import broad.core.annotation.BasicLightweightAnnotation;
import broad.core.annotation.LightweightGenomicAnnotation;
import broad.core.annotation.MaximumContiguousSubsequence;
import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.core.math.EmpiricalDistribution;
import broad.core.math.ScanStatistics;
import broad.core.math.Statistics;
import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.core.util.CollapseByIntersection;
import broad.pda.annotation.BEDFileParser;
import broad.pda.annotation.GTFFileParser;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.GeneTools;
import broad.pda.gene.GeneWithIsoforms;
import broad.pda.seq.alignment.AlignmentUtils;
import broad.pda.seq.graph.ChromosomeWithBubblesJGraphT.BubbleEdge;
import broad.pda.seq.graph.ChromosomeWithBubblesJGraphT.EdgeSourceType;
import broad.pda.seq.graph.ChromosomeWithBubblesJGraphT.WindowIterator;
import broad.pda.seq.graph.ChromosomeWithBubblesJGraphT;
import broad.pda.seq.graph.Path;
import edu.uci.ics.jung.graph.util.Pair;


//rewrite of my classic method making use of the direct alignment counts to avoid any voodoo counts coming from overlapping regions and unneccessary 25bp indexing
//Now support RNA-Seq by incorporating a graph structure for segmentation
public class ContinuousDataAlignmentModel{

	static Logger logger = Logger.getLogger(ContinuousDataAlignmentModel.class.getName());
	AlignmentDataModelStats data;
	Map<String, Integer> chromosomeLengths;
	int extensionFactor=0; //ChIP extension factor
	Collection<SAMRecord> introns;
	int chunkSize=DEFAULT_CHUNK_SIZE;

	double minNumberOfReadsAtEnd=1.01;//can also filter by  median of segment and see what happens
	int minNumberOfSplices;

	Map<String, Integer> maskedRegions;
	double alpha=.05;
	private int minIntronSize = 10;
	//private int maxIntronSize=1000000;
	private static final int DEFAULT_CHUNK_SIZE = 1000000;
	AlignmentDataModelStats pairedData;
	AlignmentDataModelStats strandSpecificReads;
	boolean isStrandSpecific=false;
	boolean upWeightSplices;
	private boolean findMaxContiguous;
	private boolean trimEnds;
	private double trimQuantile = 0.25;
	private int minAnnotationSize = 0;
	public static int DEFAULT_MIN_MAPPING_QUALITY = 5;
	public static int DEFAULT_INSERT_SIZE_FUDGE = 20;
	public static double DEFAULT_INS_SIZE_PVAL = 0.05;
	public static double DEFAULT_MAKE_GENE_OVERLAP = 0.3;

	/**
	 * Constructor to build the data alignment model
	 * @param data the alignment data model that stores the read data
	 * @throws IOException
	 */
	public ContinuousDataAlignmentModel(AlignmentDataModel data) throws IOException{
		this(new AlignmentDataModelStats(data), null, 0);
	}

	/**
	 * Constructor to build the data alignment model
	 * @param data the alignment data model that computes the statistics and stores the read data
	 * @throws IOException
	 */
	public ContinuousDataAlignmentModel(AlignmentDataModelStats data) throws IOException{
		this(data, null, 0);
	}

	/**
	 * Constructor to build the data alignment model
	 * @param data the alignment data model that computes the statistics and stores the read data
	 * @param EF an extension factor used to extend reads
	 * @param minSpliceSupport the minimum number of spliced reads needed to add a connection in the graph
	 * @throws IOException
	 */
	public ContinuousDataAlignmentModel(AlignmentDataModelStats data, int EF, int minSpliceSupport ) throws IOException{
		this(data, null, EF, minSpliceSupport);
	}

	/**
	 * Constructor to build the data alignment model
	 * @param data the alignment data model that computes the statistics and stores the read data
	 * @param maskFiles the mask files
	 * @param minSpliceSupport the minimum number of spliced reads needed to add a connection in the graph
	 * @throws IOException
	 */
	public ContinuousDataAlignmentModel(AlignmentDataModelStats data, File[] maskFiles, int minSpliceSupport )throws IOException{
		this(data, maskFiles, 0, minSpliceSupport);
	}


	/**
	 * Constructor to build the data alignment model
	 * @param data the alignment data model that computes the statistics and stores the read data
	 * @param maskFiles the mask files
	 * @param EF an extension factor used to extend reads
	 * @param minSpliceSupport the minimum number of spliced reads needed to add a connection in the graph
	 * @param pairedData the alignment data model that stores the stats for paired end reads and the paired end read data
	 * @param strandSpecificReads an alignment data model for the strand specific reads
	 * @throws IOException
	 */
	public ContinuousDataAlignmentModel(AlignmentDataModelStats data, File [] maskFiles, int EF, int minSpliceSupport, AlignmentDataModelStats pairedData, AlignmentDataModelStats strandSpecificReads) throws IOException{
		this(data, parseMaskFiles(maskFiles), EF, minSpliceSupport, pairedData, strandSpecificReads, DEFAULT_CHUNK_SIZE, false);		
	}

	/**
	 * Constructor to build the data alignment model
	 * @param data the alignment data model that computes the statistics and stores the read data
	 * @param maskFiles the mask files
	 * @param EF an extension factor used to extend reads
	 * @param minSpliceSupport the minimum number of spliced reads needed to add a connection in the graph
	 * @param pairedData the alignment data model that stores the stats for paired end reads and the paired end read data
	 * @param strandSpecificReads an alignment data model for the strand specific reads
	 * @param chunkSize the size of the cache chunk used
	 * @throws IOException
	 */
	public ContinuousDataAlignmentModel(AlignmentDataModelStats data, File [] maskFiles, int EF, int minSpliceSupport, AlignmentDataModelStats pairedData, AlignmentDataModelStats strandSpecificReads, boolean upweight) throws IOException{
		this(data, parseMaskFiles(maskFiles), EF, minSpliceSupport, pairedData, strandSpecificReads, DEFAULT_CHUNK_SIZE, upweight);		
	}


	/**
	 * Constructor to build the data alignment model
	 * @param data the alignment data model that computes the statistics and stores the read data
	 * @param maskFileData a map of all unalignable regions of the genome
	 * @param EF an extension factor used to extend reads
	 * @param minSpliceSupport the minimum number of spliced reads needed to add a connection in the graph
	 * @param pairedData the alignment data model that stores the stats for paired end reads and the paired end read data
	 * @param strandSpecificReads an alignment data model for the strand specific reads
	 * @throws IOException
	 */
	public ContinuousDataAlignmentModel(AlignmentDataModelStats data, Map<String, Integer> maskFileData, int EF, int minSpliceSupport, AlignmentDataModelStats pairedData, AlignmentDataModelStats strandSpecificReads) throws IOException{
		this(data, maskFileData, EF, minSpliceSupport, pairedData, strandSpecificReads, DEFAULT_CHUNK_SIZE, false);
	}

	/**
	 * Constructor to build the data alignment model
	 * @param data the alignment data model that computes the statistics and stores the read data
	 * @param maskFileData a map of all unalignable regions of the genome
	 * @param EF an extension factor used to extend reads
	 * @param minSpliceSupport the minimum number of spliced reads needed to add a connection in the graph
	 * @param pairedData the alignment data model that stores the stats for paired end reads and the paired end read data
	 * @param strandSpecificReads an alignment data model for the strand specific reads
	 * @param upWeightSplices whether to weight spliced reads based on their proportioncoverage scores
	 * @throws IOException
	 */
	public ContinuousDataAlignmentModel(AlignmentDataModelStats data, Map<String, Integer> maskFileData, int EF, int minSpliceSupport, AlignmentDataModelStats pairedData, AlignmentDataModelStats strandSpecificReads, boolean upWeight) throws IOException{
		this(data, maskFileData, EF, minSpliceSupport, pairedData, strandSpecificReads, DEFAULT_CHUNK_SIZE, upWeight);
	}

	/**
	 * Constructor to build the data alignment model
	 * @param data the alignment data model that computes the statistics and stores the read data
	 * @param maskFileData a map of all unalignable regions of the genome
	 * @param EF an extension factor used to extend reads
	 * @param minSpliceSupport the minimum number of spliced reads needed to add a connection in the graph
	 * @param pairedData the alignment data model that stores the stats for paired end reads and the paired end read data
	 * @param strandSpecificReads an alignment data model for the strand specific reads
	 * @param chunkSize the size of the cache chunk used
	 * @param upWeightSplices whether to weight spliced reads based on their proportion when computing coverage scores
	 * @throws IOException
	 */
	public ContinuousDataAlignmentModel(AlignmentDataModelStats data, Map<String, Integer> maskFileData, int EF, int minSpliceSupport, AlignmentDataModelStats pairedData, AlignmentDataModelStats strandSpecificReads, int chunkSize, boolean upWeightSplices) throws IOException{
		this.strandSpecificReads=strandSpecificReads;
		this.pairedData=pairedData;
		this.data=data;
		this.extensionFactor=EF;
		this.chromosomeLengths=data.getChromosomeLengths();
		this.maskedRegions=maskFileData;
		this.minNumberOfSplices = minSpliceSupport;	
		this.chunkSize = chunkSize;
		data.setChunkSize(this.chunkSize);
		this.upWeightSplices=upWeightSplices;
	}

	/**
	 * Constructor to build the data alignment model
	 * @param data the alignment data model that computes the statistics and stores the read data
	 * @param maskFiles the mask files
	 * @param EF an extension factor used to extend reads
	 * @param minSpliceSupport the minimum number of spliced reads needed to add a connection in the graph
	 * @param pairedData the alignment data model that stores the stats for paired end reads and the paired end read data
	 * @throws IOException
	 */
	public ContinuousDataAlignmentModel(AlignmentDataModelStats data, File[] maskFiles, int EF, int minSpliceSupport, AlignmentDataModelStats pairedData) throws IOException{
		this(data, maskFiles, EF, minSpliceSupport, pairedData, null);
	}

	/**
	 * Constructor to build the data alignment model
	 * @param data the alignment data model that computes the statistics and stores the read data
	 * @param maskFiles the mask files
	 * @param EF an extension factor used to extend reads
	 * @param minSpliceSupport the minimum number of spliced reads needed to add a connection in the graph
	 * @throws IOException
	 */
	public ContinuousDataAlignmentModel(AlignmentDataModelStats data, File[] maskFiles, int EF, int minSpliceSupport ) throws IOException{
		this(data, maskFiles, EF, minSpliceSupport, null);
	}



	public void setFindMaxContiguous(boolean findMaxContiguous) {
		this.findMaxContiguous = findMaxContiguous;

	}

	public void setTrimEnds(boolean trimEnds) {
		this.trimEnds = trimEnds;

	}	

	public void setMinContguousSegmentSize(int minSize) {
		this.minAnnotationSize = minSize;

	}
	public void setAlpha(double alpha) {
		this.alpha = alpha;
	
	}
	
	public void setUpWeightSplices(boolean setFlag) {
		this.upWeightSplices = setFlag; 
	}
	
	public void setMinNumberOfSplices(int minSpliceSupport) {
		this.minNumberOfSplices = minSpliceSupport;
	}
	
	public void setExtensionFactor(int extensionFactor) {
		this.extensionFactor = extensionFactor;
		data.setExtensionFactor(extensionFactor);
	}
	
	public void setMaskFileData(Map<String, Integer> maskedRegionMap) {
		this.maskedRegions = maskedRegionMap;
	}

	public void setMinimumMappingQuality(double minMapQual) { this.data.setMinimumMappingQuality(minMapQual);}

	public void setTrimQuantile(double quantile) {this.trimQuantile = quantile;}

	/**
	 * Get the underlying AlignmentDataModelStats object
	 * @return the AlignmentDataModelStats object
	 */
	public AlignmentDataModelStats getData() {return this.data;}
	
	/**
	 * Get an estimate of the expected read coverage per base
	 * @param chr chromosome to retrieve
	 * @return lambda
	 * @throws IOException 
	 */
	public double getLambda(String chr) throws IOException{
		return data.getLambda(chr);
	}

	/**
	 * Gets the total number of reads per chr
	 * @param chr chromosome to retrieve
	 * @return total number of reads
	 */
	public double getSum(String chr)throws IOException{
		return data.getSum(chr);
	}

	/**
	 * Gets the total number of reads per chr
	 * @param chr chromosome to retrieve
	 * @return total number of reads
	 */
	public double getNumberOfReads(String chr)throws IOException{
		return getSum(chr);
	}

	/**
	 * Gets the number of allowable bases per chr
	 * @param chr chromosome to retrieve
	 * @return total number of allowable bases
	 */
	public  double getNumberMarkers(String chr)throws IOException{
		return data.getNumberMarkers(chr);
	}

	/**
	 * Get total number of bases
	 * @return total number of bases
	 */
	public double getNumberOfBPs(){
		double rtrn=0;
		for(String chr: this.chromosomeLengths.keySet()){
			rtrn+=this.chromosomeLengths.get(chr);
		}
		return rtrn;
	}

	public int getCount(Collection<Annotation> Annotation) throws IOException {		
		return data.getCounts(Annotation, this.extensionFactor);
	
	}

	
	public int[] getGeneCounts(Collection<Gene> genes)throws IOException{
		return data.getGeneCounts(genes, this.extensionFactor);
	}
	
	public int[] getGeneCountsWithoutDatasetTotals(Collection<Gene> genes)throws IOException{
		return data.getGeneCountsWithoutDatasetTotals(genes, this.extensionFactor);
	}
	
	
	

	public double[] getCountsPerBp(Annotation align, IntervalTree<Alignment> tree){
		double[] vals=new double[align.getSize()];
	
		for(int i=align.getStart(); i<align.getEnd(); i++){
			double counter=data.getCountsPerAlignment(new Alignments(align.getChr(), i, i), tree, this.extensionFactor); //consider removing extension factor
			vals[i-align.getStart()]=counter;
		}
	
		return vals;
	}

	public List<Double> getCountsPerBp(Gene gene, IntervalTree<Alignment> tree){
		ArrayList<Double> rtrn=new ArrayList<Double>();
		//double[] vals=new double[align.getSize()];
		//String chr=align.getChr();
	
		for(int i=0; i<gene.getExons().length; i++){
			double[] vals=this.getCountsPerBp(gene.getExons()[i], tree);
			for(int j=0; j<vals.length; j++){rtrn.add(vals[j]);}
		}
	
		return rtrn;
	}

	private List<Double> getNormalizedCountsPerBp(Gene gene, IntervalTree<Alignment> tree) throws IOException{
		ArrayList<Double> rtrn=new ArrayList<Double>();
		//double[] vals=new double[align.getSize()];
		//String chr=align.getChr();
	
		for(int i=0; i<gene.getExons().length; i++){
			double[] vals=getCountsPerBp(gene.getExons()[i], tree);
			for(int j=0; j<vals.length; j++){rtrn.add(vals[j]/getLambda(gene.getChr()));}
		}
	
		return rtrn;
	}

	public Map<Gene, double[]> getGeneExpression(Collection<Gene> genes)throws IOException{
		return data.getGeneExpression(genes, extensionFactor);
	}

	public int getCounts(Map<String, IntervalTree<Annotation>> geneTree, Collection<String> multimappers)throws IOException{
	
		int count=0;
	
		for(String chr: geneTree.keySet()){
			CloseableIterator<Alignment> iter=data.getData().getAnnotationOverlappingRegion(new Alignments(chr, 0, this.getChromosomeLength(chr)));
			while(iter.hasNext()){
				Alignment read=iter.next();
				if(multimappers!=null && multimappers.contains(read.getReadName())){}
				else{
					Iterator<Node<Annotation>> genes=geneTree.get(chr).overlappers(read.getStart(), read.getEnd());
					if(genes.hasNext()){count++;}
				}
			}
			iter.close();
		}
	
		return count;
	}

	private Gene get3PrimeBases(Gene gene, int size) {
	
		if(gene.getOrientation().equals(Strand.POSITIVE)){
			int relativeStart=gene.getTranscriptLength()-size;
			int relativeEnd=gene.getTranscriptLength();
			Gene align=gene;
			if(relativeStart>0){align=gene.trimGene(relativeStart, relativeEnd);}
			return align;
		}
		else if(gene.getOrientation().equals(Strand.NEGATIVE)){
			int relativeStart=0;
			int relativeEnd=size;
			Gene align=gene;
			if(relativeEnd<gene.getTranscriptLength()){align=gene.trimGene(relativeStart, relativeEnd);}
			return align;
		}
		else{
			//System.err.println(gene.getAlignment().toUCSC()+" "+gene.getOrientation());
			int relativeStart=gene.getTranscriptLength()-size;
			int relativeEnd=gene.getTranscriptLength();
			Gene align=gene;
			if(relativeStart>0){align=gene.trimGene(relativeStart, relativeEnd);}
			return align;
		}
	}

	public int getChromosomeLength(String chr){
		if(this.chromosomeLengths.containsKey(chr)){return this.chromosomeLengths.get(chr);}
		return 0;
	}

	public Map<String, Integer> getChromosomeLengths() {
		return data.getChromosomeLengths();
	}

	public static int[] getWidths(String str){
		String[] vals=str.split(",");
		int[] rtrn=new int[vals.length];
	
		for(int i=0; i<vals.length; i++){
			rtrn[i]=new Integer(vals[i]);
		}
	
		return rtrn;
	}

	private Collection<? extends Annotation> getIntrons(Alignment record){
		return getIntronsFromExons(toExons(record));
	}
	
	private Collection<? extends Annotation> getIntronsFromExons(Collection<? extends Annotation> exons){
		Gene gene=new Gene(exons);
		Collection<? extends Annotation> rtrn=gene.getIntronSet();
		return rtrn;
	}
	

	private Collection<Annotation> getIntrons(Collection<Gene> allPaths) {
		Collection<Annotation> introns=new TreeSet<Annotation>();
	
		for(Gene gene: allPaths){
			Collection<? extends Annotation> s=gene.getIntronSet();
			for(Annotation intron: s){
				intron.setOrientation(gene.getOrientation());
				introns.add(intron);
			}
		}
	
		return introns;
	}

	private Collection<Annotation> getExons(Collection<Gene> allPaths) {
		Collection<Annotation> exons=new TreeSet<Annotation>();
	
		for(Gene gene: allPaths){
			exons.addAll(gene.getExonSet());
		}
	
		return exons;
	}

	

	private int getPathCoverageData(Path path, int numIntrons, TreeMap<Annotation, Double> exonMap, TreeMap<Annotation, Double> intronMap, Path overlappingPath) {
		Gene gene=overlappingPath.toGene();
		Set<? extends Annotation> exons = gene.getExonSet();
		for(Annotation exon : exons) {
			exonMap.put(exon, path.getNodeCount(exon));
		}
		Collection<BubbleEdge> introns = overlappingPath.getEdges();
		for(BubbleEdge intron : introns) {
			intronMap.put(new Alignments(intron.getConnection()), intron.getSplicedCounts());
		}			
		numIntrons=Math.max(numIntrons, gene.getNumExons()-1);
		return numIntrons;
	}

	private Map<Annotation, Double> getNodeCounts(Collection<Annotation> decollapsed, Map<Annotation, Integer> exons) throws IOException{
		Map<Annotation, Double> rtrn=new TreeMap<Annotation, Double>();
	
		Map<String, IntervalTree<Annotation>> trees=CollapseByIntersection.makeIntervalTree(exons.keySet());
	
		Collection<Annotation> remainder=new TreeSet<Annotation>();
	
		for(Annotation exon: decollapsed){
			IntervalTree<Annotation> exonChrTree = trees.get(exon.getChr());
			if(exonChrTree != null) {
				Iterator<Node<Annotation>> iter= exonChrTree.overlappers(exon.getStart(), exon.getEnd());
				int i=0;
				double count=0;
				while(iter.hasNext()){
					Annotation align=iter.next().getValue();
					count=exons.get(align);
					i++;
				}
				if(i>1){
					//compute from scratch
					remainder.add(exon);
					//System.out.println(exon);
				}
				rtrn.put(exon, count);
			}
		}
	
		Map<Annotation, Double> counter=countExons(remainder);
		rtrn.putAll(counter);
	
		return rtrn;
	}

	private Map<Annotation, double[]> getDataForAnnotation(Collection<Annotation> Annotation) throws IOException{
		Map<Annotation, double[]> rtrn=new TreeMap<Annotation, double[]>();
	
		//cache in memory by chunks
		IntervalTree<Alignment> tree=null;
		int sizeContains=0;
		for(Annotation align: Annotation){
	
			if(tree==null){tree=data.getIntervalTree(align.getChr(), align.getStart(), align.getStart()+this.chunkSize); sizeContains=align.getStart()+this.chunkSize;}
			else if(align.getStart()>=sizeContains || align.getEnd()>=sizeContains){tree=data.getIntervalTree(align.getChr(), align.getStart(), align.getEnd()+this.chunkSize); sizeContains=align.getEnd()+this.chunkSize;}
			double[] vals=getCountsPerBp(align, tree);
			rtrn.put(align, vals);
	
		}
	
	
		return rtrn;
	}

	/**
	 * This function returns an array of type double with counts for each base pair across the genomic region of the gene
	 * @author skadri
	 * @param align
	 * @return
	 * @throws IOException
	 */
	public double[] getDataForAlignment(Annotation align) throws IOException{	
		
		IntervalTree<Alignment> tree=data.getIntervalTree(align.getChr(), align.getStart(), align.getEnd()); //sizeContains=align.getStart()+this.chunkSize;
		
		double[] vals=getCountsPerBp(align, tree);
		
		return vals;
	}
	
	public List<Double> getDataForGene(Gene gene) throws IOException{
		
		IntervalTree<Alignment> tree=data.getIntervalTreeCached(gene.getChr(), gene.getStart(), gene.getEnd());
				
		List<Double> vals=getCountsPerBp(gene, tree);
		
		return vals;
	}	
	
	public Map<Gene, List<Double>> getDataForGene(Collection<Gene> Annotation) throws IOException{
		Map<Gene, List<Double>> rtrn=new TreeMap<Gene, List<Double>>();
	
		//TODO Implement a caching scheme for the tree
	
		long treeTime=0;
		//long sortedTime=0;
		long countTime=0;
	
		//cache in memory by chunks
		for(Gene align: Annotation){
			long start=System.currentTimeMillis();
			IntervalTree<Alignment> tree=data.getIntervalTreeCached(align.getChr(), align.getStart(), align.getEnd());
			long end=System.currentTimeMillis();
			treeTime+=(end-start);
	
			/*start=System.currentTimeMillis();
			Collection<Annotation> exons=align.getSortedAndUniqueExons();
			end=System.currentTimeMillis();
			sortedTime+=(end-start);
	*/
			start=System.currentTimeMillis();
			List<Double> vals=getCountsPerBp(align, tree);
			end=System.currentTimeMillis();
			countTime+=(end-start);
	
			rtrn.put(align, vals);
		}
	
		System.err.println("Tree Time: "+treeTime);
		//System.err.println("Sorted Time: "+sortedTime);
		System.err.println("Count Time: "+ countTime);
	
		//data.resetTreeCache();
	
		return rtrn;
	}

	

	public AlignmentDataModelStats getPairedData() {
		return this.pairedData;
	}

	private Collection<Annotation> scan(int windowSize, double alpha,  String chr)throws IOException{
		double T=getNumberMarkers(chr);
	
		long start=System.currentTimeMillis();
		int criticalValue=calculateCriticalValue(new Double(T).intValue(), windowSize, T, alpha, getLambda(chr));
		long end=System.currentTimeMillis();
	
		logger.info("Computing critical values  for window size "+ windowSize+" took: "+(end-start)/1000.0 + " sec. " +criticalValue);
	
		start=System.currentTimeMillis();
		IntervalTree<Annotation> significantWindows = scanGenome(windowSize, criticalValue, chr);
		end=System.currentTimeMillis();
	
		logger.info("Scanning window NEW took: "+(end-start)/1000.0 +  " sec.");
	
	
		start=System.currentTimeMillis();
		Collection<Annotation> windows = dedup(significantWindows);
	
		logger.info("Going to findMaxContiguous segments " + findMaxContiguous + " trim ends? " + trimEnds);
	
		if(findMaxContiguous) {
			logger.info("Finding max contiguous regions within windows ... ");
			windows =findMaxContiguous(windows);
			logger.info("Done");
		}
	
		if(trimEnds) {
			logger.info("Trimming window ends ... ");
			windows = trimEnds(windows);
			end=System.currentTimeMillis();
			logger.info(" took: "+(end-start)/1000.0);
			logger.info(" Done");
		}
	
	
		return windows;
	}

	public Collection<Annotation> scan(int[] windowSizes, double alpha) throws IOException{
		Collection<Annotation> rtrn=new TreeSet<Annotation>();
	
		for(String chr: this.chromosomeLengths.keySet()){
			Collection<Annotation> set=this.scan(windowSizes, alpha, chr);
			rtrn.addAll(set);
		}
	
		return rtrn;
	}

	public Collection<Annotation> scan(int[] windowSizes, double alpha,  String chr) throws IOException{
		Collection<Annotation> rtrn=new TreeSet<Annotation>();
	
		double lambdaVal=getLambda(chr);
		if(lambdaVal>0){
			for(int i=0; i<windowSizes.length; i++){
				Collection<Annotation> list=scan(windowSizes[i], alpha, chr);
				rtrn.addAll(list);
			}
		}
	
		System.gc();
	
	
		return rtrn;
	}

	
	public void scanGenome(ContinuousDataAlignmentModel data2, int windowSize, String chr, FileWriter writer, boolean filterSignificance) throws IOException{
		int counter=0;
	
		double numMarkers=getNumberMarkers(chr);
	
		
		//load chunks into memory: first chunk
		int chunkNumber=1;
	
		IntervalTree<Alignment> chunkAlignmentTree=data.getIntervalTree(chr, 0, chunkNumber*chunkSize);
		IntervalTree<Alignment> chunkAlignmentTree2=data2.getData().getIntervalTree(chr, 0, chunkNumber*chunkSize);
	
	
	
		double score1=0;
		double score2=0;
		
		Annotation previous=null;
		boolean cached=false;
	
		for(int i=0; i<data.getChromosomeLengths().get(chr); i++){
	
			
			//if score is 0 jump to end of fixed width
			int start=i;
			int end=start+windowSize;
			Annotation current=new Alignments(chr, start, end);
			//System.err.println("Current: " + current.toUCSC());
			double sum=0;
			double sum2=0;
	
			if(score1==0 && score2==0){
				sum=data.getCountsPerAlignment(current, chunkAlignmentTree, extensionFactor); 
				sum2=data2.getData().getCountsPerAlignment(current, chunkAlignmentTree2, extensionFactor);
					
				if(sum ==0 && sum2==0) {
					i=start+windowSize;
					//System.err.println("\tScore was 0 and so was sum, advancing from " + start + " to " +i);
				} else {
					//System.err.println("\tScore was 0, used reads to compute sum: " + sum);
				}
				cached=false;
			}else{
				/*
				if have part of the interval scored and cached then to get next score just need to compute piece not contained
				get score for base not covered and add
				get score for base covered in previous that no longer contained and subtract
				 */		
				cached = true;
				Annotation startPosition = new Alignments(chr, previous.getStart(), current.getStart());
				double subtractVal=data.getCountsOfUniqueOverlappers(startPosition, current, chunkAlignmentTree, extensionFactor);
				double subtractVal2=data2.getData().getCountsOfUniqueOverlappers(startPosition, current, chunkAlignmentTree2, extensionFactor);
				Annotation endPosition = new Alignments(chr, previous.getEnd(), current.getEnd());
				double addVal=data.getCountsOfUniqueOverlappers(endPosition, previous, chunkAlignmentTree, extensionFactor);
				double addVal2=data2.getData().getCountsOfUniqueOverlappers(endPosition, previous, chunkAlignmentTree2, extensionFactor);
				sum=(score1-subtractVal)+addVal;
				sum2=(score2-subtractVal2)+addVal2;
				//System.err.println("\tp: " + previous.toUCSC() + " c: " +current.toUCSC()+ " " +data.getCountsPerAlignment(current, chunkAlignmentTree, extensionFactor)+" "+subtractVal+" "+addVal+" "+data.getCountsPerAlignment(previous, chunkAlignmentTree, extensionFactor)+" score: "+score + " sum: " + sum);
			}
	
			double percent=(i/numMarkers)*100;
			if(counter% 1000000 ==1 || sum < 0){
				double memoryPercent=Runtime.getRuntime().freeMemory()/(double)Runtime.getRuntime().totalMemory();
				//System.err.println(counter+" % of markers done "+percent+", % memory used  "+(memoryPercent*100)+", window: "+window.toUCSC()+", sum: "+sum);
			}
			//assert(sum > 0);
			score1=sum;
			score2=sum2;
			previous=current;
	
			int midPoint=(current.getEnd()-current.getStart())/2;
			
			double ratio=(sum+1)/(sum2+1);
			
			double scaledRatio=ratio;
			if(ratio<1){scaledRatio=(-1.0/ratio);}
			
			//double logratio=Math.log(ratio)/Math.log(2);
			if((sum==0 && sum2==0)){scaledRatio=0;}
			if(sum<0 || sum2<0){}
			//else{writer.write(current.getMidPoint()+"\t"+scaledRatio+"\n");}
			else if(midPoint<data.getChromosomeLengths().get(chr)){
				double p=calculatePVal(new Double(sum2).intValue(), getLambda(chr), windowSize, getNumberMarkers(chr));
				if(!filterSignificance || p<alpha){
					writer.write(current.getChr()+"\t"+midPoint+"\t"+(midPoint+1)+"\t"+scaledRatio+"\n");
				}
			}
			
			
			if(end>=(chunkNumber*chunkSize-1)){
				chunkNumber++; 
				chunkAlignmentTree=data.getIntervalTree(chr, i, Math.max(chunkNumber*this.chunkSize, end));
				chunkAlignmentTree2=data2.getData().getIntervalTree(chr, i, Math.max(chunkNumber*this.chunkSize, end));
			}
	
	
			counter++;
		}
	
	}
	
	// to compute score take current score and subtract base pair before and add base pair after ensuring that none of the added is
	// already overlapping the segment and none of the subtracted is in segment
	// iterate from 0 to numMarkers-w
	// reverse iterate from numMarkers to w
	private IntervalTree<Annotation> scanGenome(int fixedWidth, int critVal, String chr) throws IOException{
		int counter=0;
	
		// chr size minus masked regions
		double numMarkers=getNumberMarkers(chr);
	
		IntervalTree<Annotation> rtrnTree=new IntervalTree<Annotation>();
	
		//load chunks into memory: first chunk
		int chunkNumber=1;
	
		// get interval tree for the first chunk
		// key is mapped coordinates of read
		// value is number of reads mapped to position
		IntervalTree<Alignment> chunkAlignmentTree=data.getIntervalTree(chr, 0, chunkNumber*chunkSize);
		
		double score=0;
		Annotation previous=null;
		boolean cached=false;
	
		// for every position in chromosome
		for(int i=0; i<data.getChromosomeLengths().get(chr); i++) {
				
			long startTime=System.currentTimeMillis();
			
			//if score is 0 jump to end of fixed width
			int start=i;
			int end=start+fixedWidth;
			Annotation current=new Alignments(chr, start, end); // the region of length fixedWidth beginning at start
			//System.err.println("Current: " + current.toUCSC());
			double sum=0;
	
			if(score==0) {
				
				// set sum to the number of mappings in current interval
				//sum=data.getCountsPerAlignment(current, chunkAlignmentTree, extensionFactor); 
				sum=data.getCountsPerAlignment(current, chunkAlignmentTree, 0); //9/9/12 re-did handling of extension factors at the get tree level 
	
				// if there are no reads in the current interval, skip to next interval
				if(sum ==0) {
					i=start+fixedWidth;
					//System.err.println("\tScore was 0 and so was sum, advancing from " + start + " to " +i);
				} else {
					//System.err.println("\tScore was 0, used reads to compute sum: " + sum);
				}
				cached=false; // no part of the interval has been scored
			} else {
				/*
				if have part of the interval scored and cached then to get next score just need to compute piece not contained
				get score for base not covered and add
				get score for base covered in previous that no longer contained and subtract
				 */		
				cached = true; // part of the interval has been scored
				// interval between start position of previous interval and start position of current interval - the part that has been covered
				Annotation startPosition = new Alignments(chr, previous.getStart(), current.getStart());
				// the number of mappings in current that have been counted
				double subtractVal=data.getCountsOfUniqueOverlappers(startPosition, current, chunkAlignmentTree, 0); //9/9/12 re-did handling of extension factors at the get tree level
				// interval between end position of previous interval and end position of current interval - the part that has not been covered
				Annotation endPosition = new Alignments(chr, previous.getEnd(), current.getEnd());
				// the number of mappings in current that have not been counted
				double addVal=data.getCountsOfUniqueOverlappers(endPosition, previous, chunkAlignmentTree, 0);//9/9/12 re-did handling of extension factors at the get tree level
				
				// score is the sum of previous interval
				// sum is the score of current interval
				sum=(score-subtractVal)+addVal;
				//System.err.println("\tp: " + previous.toUCSC() + " c: " +current.toUCSC()+ " " +data.getCountsPerAlignment(current, chunkAlignmentTree, extensionFactor)+" "+subtractVal+" "+addVal+" "+data.getCountsPerAlignment(previous, chunkAlignmentTree, extensionFactor)+" score: "+score + " sum: " + sum);
			}
	
			//double percent=(i/numMarkers)*100;
			//if(counter% 1000000 ==1 || sum < 0) {
				//double memoryPercent=Runtime.getRuntime().freeMemory()/(double)Runtime.getRuntime().totalMemory();
				//System.err.println(counter+" % of markers done "+percent+", % memory used  "+(memoryPercent*100)+", window: "+window.toUCSC()+", sum: "+sum);
			//}
			//assert(sum > 0);
			score=sum;
			previous=current;
	
			/***********BAD For Testing*********************/
			//if(sum>critVal){rtrn.add(align);}
			if(sum>critVal){
				// merge current interval into rtrnTree 
				Iterator<Node<Annotation>> iter=rtrnTree.overlappers(current.getStart(), current.getEnd());
				rtrnTree=mergeAndRemove(iter, current, rtrnTree);
			}
			/**********************************************/
	
			// move on to next chunk
			if(end>=(chunkNumber*chunkSize-1)){chunkNumber++; chunkAlignmentTree=data.getIntervalTree(chr, i, Math.max(chunkNumber*this.chunkSize, end));}
	
	
			//System.out.println("Free memory after iteration: " + i + " " + Runtime.getRuntime().freeMemory());
			counter++;
		}
	
		return rtrnTree;
	}




	private double[] scanPRate(Annotation first)throws IOException{
		return data.scanPRate(new Gene(first), 0);////9/9/12 re-did handling of extension factors at the get tree level
	}

	private double[] scanPRate(Annotation align, IntervalTree<Alignment> tree)throws IOException{
		String chr=align.getChr();
		double sum=data.getCountsPerAlignment(align, tree, 0);//consider getting number of unique reads as well //9/9/12 re-did handling of extension factors at the get tree level
		int count=align.getSize();
		double enrich=(sum/count)/getLambda(chr);
		//System.err.println("sum " + sum + "  - regionsize " + count + "  - lamgda Chr " + getLambda(chr) + " -  avg " + (sum/count));
		double[] rtrn={calculatePVal(new Double(sum).intValue(), getLambda(chr), count, getNumberMarkers(chr)), enrich, sum, (sum/count)};
		return rtrn;
	}

	public double[] scanPRate(Gene gene)throws IOException{
		return data.scanPRate(gene, 0);//9/9/12 re-did handling of extension factors at the get tree level
	}

	public double[] scanPRate(Gene gene, IntervalTree<Alignment> tree)throws IOException{
		return data.scanPRate(gene, tree, 0); //9/9/12 re-did handling of extension factors at the get tree level
	}

	
	public double[] scanPRate(Gene gene, IntervalTree<Alignment> tree, double localLambda) throws IOException{
		return data.scanPRate(gene, tree, localLambda, 0);//9/9/12 re-did handling of extension factors at the get tree level
	
	}

	public double[] scanPRate(Gene gene, IntervalTree<Alignment> tree, double localLambda, int numMarkers) throws IOException{
		return data.scanPRate(gene, tree, localLambda, 0, numMarkers); //9/9/12 re-did handling of extension factors at the get tree level
	}

	private double[] scanPRate(Annotation align,Map<String, IntervalTree<Annotation>> goodExonTree,	IntervalTree<Alignment> tree, double localLambda) throws IOException {
	
	
		String chr=align.getChr();
		//Get scores for align that dont overlap goodExonTree
		double sum=data.getCountsPerAlignment(align, goodExonTree, tree, 0);//9/9/12 re-did handling of extension factors at the get tree level//consider getting number of unique reads as well
		int count=align.getSize();
		double enrich=(sum/count)/localLambda;
	
		double[] rtrn={calculatePVal(new Double(sum).intValue(), localLambda, count, getNumberMarkers(chr)), enrich, sum, (sum/count)};
		return rtrn;
	}

	

	public double[] scanRegion(int fixedWidth, LightweightGenomicAnnotation region) throws IOException{
		return data.scanRegion(fixedWidth, region);	
	}

	public double[] scoreSegment(Annotation align) throws IOException{
		//System.err.println("Scoring segment: " + align.toUCSC());
		IntervalTree<Alignment> tree=data.getIntervalTreeCached(align.getChr(), align.getStart(), align.getEnd());
		double[] pRate=scanPRate(align, tree);
	
		return pRate;
	}
	

	public Map<Annotation, double[]> scoreSegments(Collection<Annotation> set) throws IOException{
		Map<Annotation, double[]> rtrn=new TreeMap<Annotation, double[]>();
	
		//Keep interval tree to speed it up
		IntervalTree<Alignment> tree=null;
		String chr="";
	
		int i=0;
		for(Annotation align: set){
			if(!chr.equalsIgnoreCase(align.getChr())){chr=align.getChr(); tree=data.getIntervalTree(chr, 0, this.chromosomeLengths.get(chr));}
			//if(i%10000 ==0){System.err.println(i+" "+align.toUCSC());}
			double[] pRate=scanPRate(align, tree);
			rtrn.put(align, pRate);
			i++;
		}
		return rtrn;
	}

	public Map<Annotation, double[]> scoreSegments(Collection<Annotation> set, String chrToUse) throws IOException { return scoreSegments(set, chrToUse, false); }
	public Map<Annotation, double[]> scoreSegments(Collection<Annotation> set, String chrToUse, boolean ignoreAlpha) throws IOException{
		Map<Annotation, double[]> rtrn=new TreeMap<Annotation, double[]>();
	
		int i=0;
		int j=0;
		for(Annotation align: set){
			j++;
			if(chrToUse==null || chrToUse.equalsIgnoreCase(align.getChr())){
				if(this.chromosomeLengths.containsKey(align.getChr())){
					IntervalTree<Alignment> tree=data.getIntervalTreeCached(align.getChr(), align.getStart()-1, align.getEnd()+1);
					double[] pRate=scanPRate(align, tree);
					//align.setCountScore(pRate[4]);//setting rpkm as gene scores
					if (ignoreAlpha || pRate[0] < alpha) {
						rtrn.put(align, pRate);
						i++;
					}
				}
			}
		}
		//System.err.println("Processed " + j + " but only " + i + " passed alpha " + alpha);
		data.resetTreeCache();
	
		return rtrn;
	}

	/**
	 * 
	 * @param align
	 * @return {calculatePVal(new Double(sum).intValue(), localLambda, geneLength, getNumberMarkers(chr)), enrich, sum, avgCoverage, rpkm, localLambda, geneLength,nominalP, fullyContained}
	 * @throws IOException
	 */
	public double[] scoreGene(Gene align) throws IOException {
		double [] pRate = null;
		if(this.chromosomeLengths.containsKey(align.getChr())){
			IntervalTree<Alignment> tree=data.getIntervalTreeCached(align.getChr(), align.getStart()-1, align.getEnd()+1);
			pRate=scanPRate(align, tree);
			align.setBedScore(pRate[4]);//setting rpkm as gene scores
		}
		return pRate;
	}
	
	public double[] scoreGene(Gene align, double localLambda) throws IOException {
		double [] pRate = null;
		if(this.chromosomeLengths.containsKey(align.getChr())){
			IntervalTree<Alignment> tree=data.getIntervalTreeCached(align.getChr(), align.getStart()-1, align.getEnd()+1);
			pRate=scanPRate(align, tree, localLambda);
			align.setBedScore(pRate[4]);//setting rpkm as gene scores
		}
		return pRate;
	}
	
	public GeneScore scoreGene(Gene align, IntervalTree<Alignment> tree, double localLambda) throws IOException {
		double [] pRate = null;
		if(this.chromosomeLengths.containsKey(align.getChr())){
			pRate=scanPRate(align, tree, localLambda);
			align.setBedScore(pRate[4]);//setting rpkm as gene scores
		}
		return new GeneScore(align, pRate);
	}
	
	public GeneScore scoreGene(Gene align, IntervalTree<Alignment> tree, double localLambda, int numMarkers) throws IOException {
		double [] pRate = null;
		if(this.chromosomeLengths.containsKey(align.getChr())){
			pRate=scanPRate(align, tree, localLambda, numMarkers);
			align.setBedScore(pRate[4]);//setting rpkm as gene scores
		}
		return new GeneScore(align, pRate);
	}
	
	public Map<Gene, double[]> scoreGenes(Collection<Gene> set, String chrToUse, IntervalTree<Alignment> tree) throws IOException{
		Map<Gene, double[]> rtrn=new TreeMap<Gene, double[]>();

		int i=0;
		for(Gene align: set){
			//System.err.println("Scoring "+ align.getName() +" " + align.toUCSC());
			if(chrToUse==null || chrToUse.equalsIgnoreCase(align.getChr())){
				double[] pRate = scoreGene(align, tree);
				if(pRate != null)
				rtrn.put(align, pRate);
				i++;
			}
		}

		data.resetTreeCache();

		return rtrn;
	}

	public GeneScore getPeak(Gene gene, Gene startCodon, IntervalTree<Alignment> tree) throws IOException {
		//TODO Lets get the width of the peak over the startCodon
		//Get all reads overlapping this region
		//then define width as RefSeqGene going from start of first read to end of last
		//set GeneScore to this count and region
		Gene peak=data.getPeak(gene, startCodon, tree, 0);
		if(peak!=null){
			double[] scores=this.scoreGene(peak, tree);
			//System.err.println("Peak");
			//System.err.println(peak);
			return new GeneScore(peak, scores);
		}
		else{
			double[] scores=this.scoreGene(startCodon, tree);
			return new GeneScore(startCodon, scores);
		}
	}
	
	
	
	
	public Gene getPeak(Gene gene, Gene startCodon) throws IOException {
		IntervalTree<Alignment> tree=data.getIntervalTree(gene.getChr(), gene.getStart(), gene.getEnd());
		return data.getPeak(gene, startCodon, tree, 0);
	}
	
	
	
	public AlignmentDataModelStats getAlignmentDataModelStats(){return this.data;}
	
	public double[] scoreGene(Gene align, IntervalTree<Alignment> tree) throws IOException {
		double [] pRate=scanPRate(align, tree);
		align.setBedScore(pRate[4]);//setting rpkm as gene scores
		return pRate;
	}

	public Map<Gene, double[]> scoreGenes(Collection<Gene> genes) throws IOException{
		return scoreGenes(genes, null);
	}

	public Map<Gene, double[]> scoreGenes(Collection<Gene> set, String chrToUse) throws IOException{
		Map<Gene, double[]> rtrn=new TreeMap<Gene, double[]>();
	
		int i=0;
		for(Gene align: set){
			//logger.info("Scoring "+ align.getName() +" " + align.toUCSC());
			if(chrToUse==null || chrToUse.equalsIgnoreCase(align.getChr())){
				double[] pRate = scoreGene( align);
				if(pRate != null)
					rtrn.put(align, pRate);
				i++;
			}
		}
	
		data.resetTreeCache();
	
		return rtrn;
	}
	
	public Map<Gene, double[]> scoreGenes(Collection<Gene> set, String chrToUse, int size) throws IOException{
		Map<Gene, double[]> rtrn=new TreeMap<Gene, double[]>();
	
		//Keep interval tree to speed it up
	
		int i=0;
		for(Gene gene: set){
			Gene align=get3PrimeBases(gene, size);
			if(chrToUse==null || chrToUse.equalsIgnoreCase(align.getChr())){
				IntervalTree<Alignment> tree=data.getIntervalTreeCached(align.getChr(), align.getStart()-1, align.getEnd()+1);
				double[] pRate=scanPRate(align, tree);
				rtrn.put(align, pRate);
				i++;
			}
		}
	
		return rtrn;
	}

	public void scoreGenes(BEDFileParser genes) throws IOException{
		for (String chr : getChromosomeLengths().keySet()) {
			if(genes.containChr(chr) & genes.getChrTree(chr) != null ) {
				Iterator<GeneWithIsoforms> chrGeneIt = genes.getChrTree(chr).valueIterator();
				while(chrGeneIt.hasNext()) {
					GeneWithIsoforms gene = chrGeneIt.next();
					Iterator<Gene> isoIt = gene.getAllIsoforms().iterator();
					while(isoIt.hasNext()) {
						Gene iso = isoIt.next();
						double [] pRate = scoreGene(iso);
						iso.setBedScore(pRate[4]);
						iso.setExtraFields(pRate);
					}
				}

			} else {
				logger.info("No gene found for chromosome " + chr);
			}
			resetTreeCache();
		}
	}
	
	
	

	public Map<Path, double[]> scorePathsForOnlyRPKM(Collection<Path> paths) {
		Map<Path, double[]> rtrn = new TreeMap<Path, double[]>();
	
		for (Path path : paths) {
			double score = path.getCoverage();
			double size = path.getSize();
			double avgCoverage = score / size;
			double rpkm = avgCoverage * data.getRPKMConstant(path.getChromosome());
	
			double[] array = {0, 0, score, avgCoverage, rpkm, 0, size, 0};
			rtrn.put(path, array);
		}
	
		return rtrn;
	}

	private static Map<Path, double[]> scorePathsForOnlyRPKM(Collection<Path> paths, double rpkmConstant) {
		Map<Path, double[]> rtrn = new TreeMap<Path, double[]>();
	
		for (Path path : paths) {
			double score = path.getCoverage();
			double size = path.getSize();
			double avgCoverage = score / size;
			double rpkm = avgCoverage * rpkmConstant;
	
			double[] array = {0, 0, score, avgCoverage, rpkm, 0, size, 0};
			rtrn.put(path, array);
		}
	
		return rtrn;
	}

	public Map<Path, double[]> scorePaths(Collection<Path> paths, double lambda) throws IOException{
		Map<Path, double[]> rtrn=new TreeMap<Path, double[]>();
	
	
		int i=0;
		for(Path path: paths){
			double score=path.getCoverage();
			int size=path.getSize();
			double scanP=calculateApproximatePVal(new Double(score).intValue(), lambda, size, data.getNumberMarkers(path.getChromosome()), this.alpha);
			double localP=calculateApproximatePVal(new Double(score).intValue(), path.getLocalLambda(), size, data.getNumberMarkers(path.getChromosome()), this.alpha);
	
	
			double avgCoverage = score/size;
	
			double enrich=avgCoverage/lambda;
			double rpkm = avgCoverage * data.getRPKMConstant(path.getChromosome());
			double[] array={scanP, enrich, score, avgCoverage, rpkm, lambda, size, localP};
	
			rtrn.put(path, array);
			i++;
			//if(i%1000 ==0){System.err.println(i+" "+paths.size());}
		}
	
		return rtrn;
	}

	public static Map<Path, double[]> scorePaths(Collection<Path> paths, ChromosomeWithBubblesJGraphT graph, double alpha) {
		return scorePaths(paths, graph.getLambda(), graph.getNumberOfMarkers(), graph.getNumberOfReads(), graph.getRPKMConstant(), graph.getLocalRate(), alpha);
	}

	public static Map<Path, double[]> scorePaths(Collection<Path> paths, double lambda, double numMarkers, double numReads, double RPKMContant, Map<Annotation, Double> localRates, double alpha){
		Map<Path, double[]> rtrn=new TreeMap<Path, double[]>();
	
	
		for(Path path: paths){
	
			double score=path.getCoverage();
			int size=path.getSize();
			double scanP=calculateApproximatePVal(new Double(score).intValue(), lambda, size, numMarkers, alpha);
	
			double localLambda=0;
			if(path.getLocalLambda()==0){
				Collection<? extends Annotation> exons=path.toGene().getExonSet();
				for(Annotation exon: exons){
					if(localRates!=null && localRates.containsKey(exon)){
						localLambda=Math.max(localLambda, localRates.get(exon));
					}
				}
			}
			else{localLambda=path.getLocalLambda();}
	
			double localP=calculateApproximatePVal(new Double(score).intValue(), localLambda, size, numMarkers, alpha);
	
	
			double avgCoverage = score/size;
	
			double enrich=avgCoverage/lambda;
			double rpkm = avgCoverage * RPKMContant;
			double[] array={scanP, enrich, score, avgCoverage, rpkm, lambda, size, localP};
	
			//if(path.toGene().getNumExons()==1){System.err.println(path+"\t"+scanP+"\t"+localP+"\t"+localLambda);}
	
			//if(path.toGene().getNumExons()==1){System.err.println(score);}
	
			rtrn.put(path, array);
		}
	
		return rtrn;
	}

	public double count(Annotation align)throws IOException{
		return data.count(align, 0);//9/9/12 re-did handling of extension factors at the get tree level
	}

	private Map<Annotation, Double> countExons(Collection<Annotation> exons) throws IOException{
		Map<Annotation, Double> rtrn=new TreeMap<Annotation, Double>();
	
		int i=0;
		for(Annotation exon: exons){
			//System.err.println(i+" "+exons.size()+" "+exon);
			double val=data.countWithinExon(exon, 0);
			//int val=new Double(data.count(exon, 0)).intValue();
			rtrn.put(exon, val);
			i++;
		}
	
		return rtrn;
	}

	private Collection<Gene> collapse(Collection<Gene> genes){
		Collection<Gene> rtrn=new TreeSet<Gene>();
	
		Collection<Annotation> exons=new TreeSet<Annotation>();
		for(Gene gene: genes){
			exons.addAll(gene.getExonSet());
		}
	
		exons=CollapseByIntersection.collapseByIntersection(exons, false);
	
		for(Annotation exon: exons){rtrn.add(new Gene(exon));}
	
		return rtrn;
	}

	private Collection<Annotation> collapsePaths(Collection<Path> paths) {
		Collection<Annotation> temp=new TreeSet<Annotation>();
	
		for(Path path: paths){
			temp.add(path.toGene().getAlignment());
		}
	
		return CollapseByIntersection.collapseByIntersection(temp, false);
	}

	private Collection<Gene> collapseAllSingleExons(Collection<Gene> c){
		Collection<Gene> rtrn=new TreeSet<Gene>();
	
		Collection<Annotation> singleExon=new TreeSet<Annotation>();
	
		for(Gene gene: c){
			if(gene.getNumExons()>1){rtrn.add(gene);}
			else{singleExon.add(gene.getAlignment());}
		}
	
		singleExon=CollapseByIntersection.collapseByIntersection(singleExon, false);
	
		for(Annotation exon: singleExon){rtrn.add(new Gene(exon));}
	
		return rtrn;
	}
	
	private Annotation collapseMaxContiguous(Annotation align, double[] array)  throws IOException{
		array=subtract(array, Math.max(getLambda(align.getChr()), this.minNumberOfReadsAtEnd));
	
		//double[] array=this.getDataForAnnotation(align);
		double[] maxSum=MaximumContiguousSubsequence.maxSubSum3(array);
		//System.err.println("region " + align.toUCSC() + " dist" + maxSum[0]+" "+maxSum[1]+"-"+maxSum[2]);
		Annotation newAlign=null;
		if(maxSum[0]>0){
			newAlign=new Alignments(align.getChr(),
					new Double(((align.getStart()+(maxSum[1]-1)))).intValue(),
					new Double(align.getStart()+(maxSum[2]-1)).intValue());
		}
		return newAlign;
	}

	/**
	 * Trim a region to a contiguous subregion with max number of positions having coverage above the quantile
	 * @param region The original region
	 * @param quantile The coverage quantile of the original region
	 * @return The contiguous subregion with maximal number of positions having coverage above the quantile of the original region or null if can't trim
	 * @throws IOException 
	 */
	public Gene trimMaxContiguous(Gene region, double quantile) throws IOException {
		
		int treeStartPos = Math.min(0,region.getStart());
		int treeEndPos = region.getEnd();
		IntervalTree<Alignment> tree = getData().getData().getIntervalTree(region.getChr(), treeStartPos, treeEndPos);

		List<Double> coverageData = getCountsPerBp(region,tree);
		Gene rtrn = collapseMaxContiguous(region, coverageData, quantile);
		int trimmedStart = rtrn.getStart();
		int trimmedEnd = rtrn.getEnd();
		rtrn.setName(rtrn.getChr() + ":" + trimmedStart + "-" + trimmedEnd);
		rtrn.setOrientation(region.getOrientation());
		rtrn.setCDSRegion(Math.max(trimmedStart, region.getStart()), Math.min(trimmedEnd, region.getEnd()));
		
		return rtrn;
		
	}
	
	private Gene collapseMaxContiguous(Gene align, List<Double> data, double quantile) throws IOException{
		double[] array=l2a(data);
		Collections.sort(data);
		double cutoff = Math.max(minNumberOfReadsAtEnd, Statistics.quantile(data, quantile));
		array=subtract(array, Math.max(getLambda(align.getChr()), cutoff));
	
		//double[] array=this.getDataForAnnotation(align);
		double[] maxSum=MaximumContiguousSubsequence.maxSubSum3(array);
	
		Gene newAlign=null;
		if(maxSum[0]>0){
			if(maxSum[1]-1<0 && maxSum[2]-1<0){}
			else{
				newAlign=align.trimGene(new Double(maxSum[1]-1).intValue(), new Double(maxSum[2]-1).intValue());
			}
		}
		return newAlign;
	}

	public static double calculatePVal(int k, double lambda, double w, double T){
		return ScanStatistics.calculatePVal(k, lambda, w, T);
	}

	private int calculateCriticalValue(int maxLength, double w, double T, double alpha, double lambda){
		int num=0;
		for(int k=1; k<maxLength; k++){
			num++;
			double p=calculatePVal(k, lambda, w, T); //if this gets too slow can consider precomputing
			//System.err.println(k+" "+p+" "+lambda+" "+w+" "+T);
			if(p<alpha){break;}
		}
		return num;
	}

	private static double calculateApproximatePVal(int k, double lambda, double w, double T, double alpha){
		return ScanStatistics.calculateApproximatePVal(k, lambda, w, T, alpha);
	}

	/**
	 * Segments locally, computing a local lambda using the overlapper genes
	 * lambda may be different than chromosome lambda only if an overlapping genes is multiexonic.
	 * @param genes
	 * @param tree
	 * @param chr
	 * @return
	 */
	private double computeLambda(Iterator<Node<Gene>> genes, IntervalTree<Alignment> tree, String chr) throws IOException{
	
		int numIntrons=0;
		Collection<Annotation> geneCollection=new TreeSet<Annotation>();
		while(genes.hasNext()){
			Gene gene=genes.next().getValue();
			numIntrons=Math.max(numIntrons, gene.getNumExons()-1);
			geneCollection.add(gene.getAlignment());
		}
	
		Collection<Annotation> collapsed=CollapseByIntersection.collapseByIntersection(geneCollection, false);
	
		double reads=0;
		int counts=0;
		//Go through each gene and get reads from start to end
		for(Annotation align: collapsed){
			reads+=data.getCountsPerAlignment(align, tree, 0);
			counts+=align.getSize();
		}
	
		double lambda=getLambda(chr);
		if(numIntrons>0) {
			lambda=Math.max(reads/counts, lambda);
		} 
		return lambda;
	}

	private double computeLambdaFromPath(Iterator<Node<Path>> genes, IntervalTree<Alignment> tree, Path path, IntervalTree<LightweightGenomicAnnotation> cache) throws IOException{
	
		Iterator<Node<LightweightGenomicAnnotation>> pathComputedOverlappers = cache.overlappers(path.getStart(), path.getEnd());
		if( pathComputedOverlappers.hasNext() ) {
			LightweightGenomicAnnotation containigRegion = pathComputedOverlappers.next().getValue();
			if(containigRegion.contains(new BasicLightweightAnnotation(path.getChromosome(), path.getStart(), path.getEnd()))){
				return containigRegion.getScore();
			}
		}
	
		// If the chaced tree did not contain a region that contained the current path, go ahead and compute its local lambda.
	
		//Collection<Annotation> geneCollection=new TreeSet<Annotation>();
		TreeMap<Annotation, Double> exonMap = new TreeMap<Annotation, Double>();
		TreeMap<Annotation, Double> intronMap = new TreeMap<Annotation, Double>();
	
		int numIntrons=getPathCoverageData(path, 0, exonMap, intronMap, path);
		while(genes.hasNext()){
			Path overlappingPath=genes.next().getValue();
			numIntrons = getPathCoverageData(path, numIntrons, exonMap,intronMap, overlappingPath);
			//geneCollection.addAl(gene.getAlignment());
		}
	
		Collection<Annotation> collapsedExons =CollapseByIntersection.collapseByIntersection(exonMap.keySet(), false);
	
		double reads=0;
		int counts=0;
		//Go through each gene and get reads from start to end
		for(Annotation exon: collapsedExons){
			if (exonMap.containsKey(exon)) {
				reads+= exonMap.get(exon); 
			} else {
				reads += data.getCountsPerAlignment(exon, tree, 0);
			}
			counts+=exon.getSize();
		}
	
		//Now that we have all reads within exons, lets count spliced reads.
		Set<Annotation> intronSet = intronMap.keySet();
		for(Annotation intron: intronSet){			
			reads+= intronMap.get(intron); 
		}
	
		double lambda=getLambda(path.getChromosome());
		if(numIntrons>0) {
			lambda=Math.max(reads/counts, lambda);
		} 
	
		LightweightGenomicAnnotation mergedLoci = new BasicLightweightAnnotation(path.getChromosome(), exonMap.firstKey().getStart(), exonMap.lastKey().getEnd());
		mergedLoci.setScore(lambda);
		cache.put(mergedLoci.getStart(), mergedLoci.getEnd(), mergedLoci);
		return lambda;
	}

	private Collection<Gene> acrossGaps(IntervalTree<Alignment> chunkAlignmentTree, int critVal, int fixedWidth, GenomeWithGaps2 gwg){
		long startTime=System.currentTimeMillis();
		long endTime=0;
	
		Collection<Gene> rtrn=new TreeSet<Gene>();
	
		int counter=0;		
		for(int i=0; i<gwg.getRelativeGenomeLength(); i++){
			Gene window=gwg.getRelativeWindow(i, i+fixedWidth);
			double sum=data.getCountsPerAlignment(window, chunkAlignmentTree, 0); //9/9/12 re-did handling of extension factors at the get tree level
			if(sum>critVal){rtrn.add(window);}
	
			/*******Iteration controller***************/
			if(sum==0){i=i+fixedWidth;}
			/******************************************/
	
			counter++;
			if(counter % 10000 ==0){endTime=System.currentTimeMillis(); logger.debug(window.getAlignment().toUCSC()+" "+window.getExons().length+" "+counter+" "+gwg.getRelativeGenomeLength()+" "+(counter/(double)gwg.getRelativeGenomeLength())+" "+(endTime-startTime)); startTime=System.currentTimeMillis();}
			//TODO Collapse on the fly to avoid huge memory footprint
		}
		return rtrn;
	}

	private Collection<Gene> mergeIntoTranscripts(Collection<Gene> genes, String chr) throws IOException{
	
		Collection<Annotation> exons=new TreeSet<Annotation>();
		Collection<Annotation> introns=new TreeSet<Annotation>();
	
		for(Gene gene: genes){
			exons.addAll(gene.getExonSet());
			introns.addAll(gene.getIntronSet());
		}		
	
		long start=System.currentTimeMillis();
		exons=CollapseByIntersection.collapseByIntersection(exons, false);
	
		long end=System.currentTimeMillis();
		logger.debug("Collapse: "+(end-start));
	
		//for(Annotation exon: exons){System.out.println(exon);}
	
		start=System.currentTimeMillis();
		if(!introns.isEmpty()){
			exons=CollapseByIntersection.DecollapseByIntronLocation(exons, introns);
		}
		end=System.currentTimeMillis();
		logger.info("Decollapse: "+(end-start));
	
	
	
		//writeTest(exons, introns);
	
		logger.debug(exons.size()+" "+introns.size());
		ChromosomeWithBubblesJGraphT bubbles=new ChromosomeWithBubblesJGraphT(chr, exons, introns, data.getLambda(chr), data.getNumberMarkers(chr), data.getNumberOfReads(chr), minNumberOfSplices);
	
		Collection<Gene> rtrn=new TreeSet<Gene>();
		rtrn.addAll(bubbles.getGenePaths(0));
		rtrn.addAll(bubbles.getOrphanNodes());
	
		return rtrn;
	}

	// remove intervals in iter
	// replace with one interval consisting of the union of the span of all intervals in iter and align
	private IntervalTree<Annotation> mergeAndRemove(Iterator<Node<Annotation>> iter, Annotation align, IntervalTree<Annotation> tree){
		IntervalTree<Annotation> rtrn=tree;
	
		int start=align.getStart();
		int end=align.getEnd();
	
		while(iter.hasNext()){
			Annotation next=iter.next().getValue();
			start=Math.min(next.getStart(), start);
			end=Math.max(next.getEnd(), end);
			rtrn.remove(next.getStart(), next.getEnd());
		}
	
		Annotation newAlign=new Alignments(align.getChr(), start, end);
	
		rtrn.put(start, end, newAlign);
	
		return rtrn;
	}

	private IntervalTree<Annotation> mergeAndRemove(IntervalTree<Annotation> tree, Collection<? extends Annotation> exons){
		IntervalTree<Annotation> rtrn=tree;
	
		for(Annotation exon: exons){
			//System.out.println(exon);
			int start=exon.getStart();
			int end=exon.getEnd();
			Iterator<Node<Annotation>> overlappers=rtrn.overlappers(exon.getStart(), exon.getEnd());
			while(overlappers.hasNext()){
				Annotation next=overlappers.next().getValue();
				start=Math.min(next.getStart(), start);
				end=Math.max(next.getEnd(), end);
				rtrn.remove(next.getStart(), next.getEnd());
			}
			Annotation newAlign=new Alignments(exon.getChr(), start, end);
			rtrn.put(start, end, newAlign);
		}
		return rtrn;
	}

	//TODO: Might need to consider an iterative addition
	private Collection<Annotation> makeAdditionalNodes(Collection<Annotation> significantPieces, Map<String, IntervalTree<Annotation>> goodExons) {
		Collection<Annotation> rtrn=new TreeSet<Annotation>();
	
		for(Annotation region: significantPieces){
			Iterator<Node<Annotation>> right=goodExons.get(region.getChr()).overlappers(region.getStart(), region.getEnd()+1);
			while(right.hasNext()){
				Annotation addition=right.next().getValue();
				Annotation align=new Alignments(region.getChr(), Math.min(region.getStart(), addition.getStart()), Math.max(region.getEnd(), addition.getEnd()));
				rtrn.add(align);
			}

			Iterator<Node<Annotation>> left=goodExons.get(region.getChr()).overlappers(region.getStart()-1, region.getEnd());
			while(left.hasNext()){
				Annotation addition=left.next().getValue();
				Annotation align=new Alignments(region.getChr(), Math.min(region.getStart(), addition.getStart()), Math.max(region.getEnd(), addition.getEnd()));
				rtrn.add(align);
			} 

			Iterator<Node<Annotation>> full=goodExons.get(region.getChr()).overlappers(region.getStart()-1, region.getEnd());
			//TODO Should add both to the same exon
			while(full.hasNext()){
				Annotation addition=full.next().getValue();
				Annotation align=new Alignments(region.getChr(), Math.min(region.getStart(), addition.getStart()), Math.max(region.getEnd(), addition.getEnd()));
				rtrn.add(align);
			}
		}
	
		return rtrn;
	}



	private Collection<Annotation> trimEnds(Collection<Annotation> windows) throws IOException{
		Set<Annotation> rtrn = new TreeSet<Annotation>();
	
		//loop through each alignment and compute number of reads for the ends
		Map<Annotation, double[]> counts=this.getDataForAnnotation(windows);
		for(Annotation align: windows){
			Annotation trunc=trimEnds(align, counts.get(align));
			if(trunc!=null){rtrn.add(trunc);}
		}
	
		return rtrn;
	}

	private Collection<Gene> trimEndsForGenes(Collection<Gene> allSegments, double quantile) throws IOException{
		Set<Gene> rtrn=new TreeSet<Gene>();
	
		long start=System.currentTimeMillis();
		//loop through each alignment and compute number of reads for the ends
		logger.debug("Getting reads for annotations");
		//System.err.println("Getting reads for annotations");
		Map<Gene, List<Double>> counts=getDataForGene(allSegments);
		long end=System.currentTimeMillis();
		logger.debug("Counting data: "+(end-start));
		//System.err.println("Counting data: "+(end-start)+"ms.");
		
		for(Gene align: allSegments){
			List<Double> alignData =  counts.get(align);
			if(trimEnds){
			//	start=System.currentTimeMillis();
				Gene trunc = trimEnds(align, alignData, quantile);
	
				if(trunc!=null){rtrn.add(trunc);}
			//	end=System.currentTimeMillis();
			//	System.err.println("Trimming ends for Gene: "+align.getName()+" in "+(end-start)+"ms.");
			}
			else{
			//	start=System.currentTimeMillis();
				Gene trunc=collapseMaxContiguous(align, alignData, quantile);
				if(trunc != null) {
					rtrn.add(align);
				}
			//	end=System.currentTimeMillis();
			//	System.err.println("Collapsing max contiguous for Gene: "+align.getName()+" in "+(end-start)+"ms.");
			}
		}
	
		return rtrn;
	
	}
	
	/**
	 * This function will filter the genes in the input Map by a number of constraints on the total number of reads on the genes
	 * 1. All genes with expression < 5.0 read counts will be removed
	 * 2. All genes with expression in the 75th quantile of the remaining expressed genes will be returned 
	 * @param countsMap: A map of RefSeqGene to a List of counts along the length of the gene 
	 * @return
	 */
	public List<Gene> filterGenesByExpression(Map<Gene, List<Double>> countsMap){
		
		Map<Gene, Double> SumCountsMap=new HashMap<Gene, Double>();
		List<Gene> genes = new ArrayList<Gene>();
		/*
		 * STEP 1 : GET A DISTRIBUTION OF TOTAL COUNTS OF EACH GENE
		 */
		int removed = 0;
		for(Gene gene:countsMap.keySet()){
			/*
			 * STEP 2: DO NOT ADD ANY GENES WITH NUMBER OF COUNTS <5
			 */
			double sum = Statistics.sum(countsMap.get(gene));
			if(sum>=5.0){
				SumCountsMap.put(gene, sum);
			}
			else{
				removed ++;
			}
		}
		
		//TO DO: CHECK the quantile function needs an ordered list
		System.out.println(removed+" genes removed because of expression < 5 reads. Considering remaining "+SumCountsMap.size()+" genes.");
		ArrayList<Double> C = new ArrayList(SumCountsMap.values());
		Collections.sort(C);
		double quant = Statistics.quantile(C, 0.25);
		
		for(Gene gene:SumCountsMap.keySet()){
			if(SumCountsMap.get(gene)>quant){
				genes.add(gene);
			}
		}
		System.out.println((SumCountsMap.size()-genes.size())+" genes removed because of expression < 25th quantile. Considering remaining "+genes.size()+" genes.");
		return genes;
	}
	

	private Annotation trimEnds(Annotation align, double[] array){
		double min = Math.max(minNumberOfReadsAtEnd, Statistics.quantile(array,trimQuantile));
		//double[] array=this.getDataForAnnotation(align);
		int trimStart = align.getStart() + MaximumContiguousSubsequence.contiguousStartSubSequenceOverMin(array, min);
		int trimEnd   = align.getStart() + MaximumContiguousSubsequence.contiguousEndSubSequenceOverMin(array, min);
		Annotation newAlign=null;
		if(trimEnd > trimStart){
			if(trimEnd - trimStart < minAnnotationSize) {
				int toAdd = minAnnotationSize - trimEnd + trimStart + 1;
				trimEnd = trimEnd + toAdd/2;
				trimStart = trimStart - toAdd/2;
			}
			newAlign=new Alignments(align.getChr(),trimStart, trimEnd);
		} 
		//System.err.println("region " + align.toUCSC() + " min coverage " + min + " trimStart " + trimStart + " trimEnd  " + trimEnd);
		return newAlign;
	}

	private Gene trimEnds(Gene align, List<Double> data, double quantile){
		double[] array=l2a(data);
	//	long start=System.currentTimeMillis();
		Collections.sort(data);
	//	long end=System.currentTimeMillis();
	//	System.err.println("Sorting in "+(end-start)+"ms.");
		
		double cutoff = Math.max(minNumberOfReadsAtEnd, Statistics.quantile(data, quantile));
		//System.err.print("Trimming " + align.getName() + " (" + align.getOrientation() + ") cutoffs " + cutoff);
		//double[] array=this.getDataForAnnotation(align);
		int trimStart = MaximumContiguousSubsequence.contiguousStartSubSequenceOverMin(array, cutoff);
		int trimEnd   =  MaximumContiguousSubsequence.contiguousEndSubSequenceOverMin(array, cutoff);
		logger.debug(" trimStart " + trimStart + " trimEdn " + trimEnd + ", transcript length " + align.getTranscriptLength());


		// We only want to trim the last exons not the cut into spliced ones
		if(align.getNumExons() > 1 ) {
			//int genomicTrimStart = align.transcriptToGenomicPosition(trimStart);
			//int genomicTrimEnd   = align.transcriptToGenomicPosition(trimEnd);
			
			int transcriptFirstExonEnd = "-".equals(align.getOrientation())  
					? align.genomicToTranscriptPosition(align.getLastExon().getStart())
					: align.genomicToTranscriptPosition(align.getFirstExon().getEnd() - 1);

			int transcriptLastExonStart = "-".equals(align.getOrientation()) 
					? align.genomicToTranscriptPosition(align.getFirstExon().getEnd() - 1)
					: align.genomicToTranscriptPosition(align.getLastExon().getStart() );
			
			logger.trace("first exon end in transcript " + transcriptFirstExonEnd + " last exon start " + transcriptLastExonStart);

			if(trimStart > transcriptFirstExonEnd) {
				trimStart = Math.max(0, transcriptFirstExonEnd - 50);
			}
			
			if(trimEnd < transcriptLastExonStart) {
				trimEnd = Math.min(align.getTranscriptLength(), transcriptLastExonStart + 50);
			}

			if("-".equals(align.getOrientation()) ){   
				int tmpTrimStart = trimStart;
				trimStart = align.getTranscriptLength() - trimEnd;
				trimEnd = align.getTranscriptLength() - tmpTrimStart;
			}
			
			logger.trace("Reset trimStart and TrimEnd to  " + trimStart + " - " + trimEnd);
			
			//trimStart = genomicTrimStart > align.getFirstExon().getEnd() ? 1 : trimStart;
			//trimEnd = genomicTrimEnd < align.getLastExon().getStart() ? align.getTranscriptLength() : trimEnd;
			//logger.debug("genomic trim end: " + genomicTrimEnd + " genomic trim start " + genomicTrimStart +/* " lastExon start: " + lastExonStart + " firstExonEnd: " + firstExonEnd +*/ " Reset trimStart and TrimEnd to  " + trimStart + " - " + trimEnd);
		}
	
		
		Gene newAlign=align;
		
		if(trimEnd > trimStart && (1 > trimStart || trimEnd < align.getTranscriptLength()-1)){
			newAlign=align.trimGene(trimStart, trimEnd);
			logger.debug("trimming ("+trimStart +" - "+ trimEnd+") gene was: " + align.toBED() + " and now is: " +newAlign.toBED());
		} 
		return newAlign;
	}

	private Collection<Gene> trimGenes(Collection<Gene> set, double quantile) throws IOException{
		return trimEndsForGenes(set, quantile);
	}

	public Gene trimGeneEnds(Gene gene, double quantile) throws IOException{
		long start=System.currentTimeMillis();
		//loop through each alignment and compute number of reads for the ends
		List <Gene> oneGeneList = new ArrayList<Gene>(1);
		oneGeneList.add(gene);
		List<Double> counts=getDataForGene(oneGeneList).get(gene);
		long end=System.currentTimeMillis();
	
		Gene trunc = trimEnds(gene, counts, quantile);
		return trunc;
	
	}

	private static void writeFullBED(String save, Map<Path, double[]> segments)throws IOException{
		FileWriter writer=new FileWriter(save);
	
		for(Path align: segments.keySet()){
			Gene gene=align.toGene();
			double[] vals=segments.get(align);
			for(double val : vals) {
				gene.addExtraField(String.valueOf(val));
			}
			//gene.setName(align.getLocalLambda()+"_"+vals[0]+"_"+vals[7]);
			//gene.setName(""+vals[2]);
			writer.write(gene+"\n");
		}
	
		writer.close();
	}

	public static void writeFullBED(String save, Map<Gene, double[]> map, double alpha)throws IOException{
		FileWriter writer=new FileWriter(save);
	
		for(Gene align: map.keySet()){
			double[] ps=map.get(align);
			if(ps[0]<alpha){
				writer.write(align.toBED());
				for(double p : ps) {
					writer.write("\t"+p);
				}
				writer.write("\n");
			}else {
				logger.debug("Gene " + align.getName() + " in " + align.toUCSC() + " alpha " + ps[0] +" not printed");
			}
	
		}
	
		writer.close();
	}

	public static void writeFullBED(String save, Collection<Gene> segments)throws IOException{
		BufferedWriter writer=new BufferedWriter(new FileWriter(save));
	
		writeFullBED(writer, segments);
	
		writer.close();
	}

	private static void writeFullBED(BufferedWriter writer, Collection<Gene> segments)throws IOException{
		for(Gene align: segments){
			writer.write(align.toBED());
			writer.newLine();
		}
	}

	public void writeConnections(String save, ChromosomeWithBubblesJGraphT graph)throws IOException{
		FileWriter writer=new FileWriter(save);
	
		IntervalTree<BubbleEdge> edges=new IntervalTree<BubbleEdge>();
	
		for(BubbleEdge edge: graph.edgeSet()){
			edges.put(edge.getConnection().getStart(), edge.getConnection().getEnd(), edge);
		}
	
		Collection<Path> paths=graph.getPaths(1);
		Collection<Annotation> transcribedRegions=collapsePaths(paths);
	
		for(Annotation transcript: transcribedRegions){
			writer.write("GENE: "+transcript.toUCSC()+"\n");
			Iterator<Node<BubbleEdge>> edgeIter=edges.overlappers(transcript.getStart(), transcript.getEnd());
			Collection<Annotation> exons=new TreeSet<Annotation>();
			while(edgeIter.hasNext()){
				BubbleEdge edge=edgeIter.next().getValue();
				if(edge.getType().equals(EdgeSourceType.SPLICED)){
					VertexPair<Annotation> nodes=graph.getNodePair(edge);
					Annotation first=nodes.getFirst();
					Annotation second=nodes.getSecond();
					writer.write(first.toUCSC()+"\t"+second.toUCSC()+"\t"+edge.getSplicedCounts()+"\t"+data.getSpliceWeightFactor()+"\n");
					exons.add(first);
					exons.add(second);
				}
			}
			for(Annotation exon: exons){
				writer.write(exon.toUCSC()+"\t"+scanPRate(exon)[2]+"\n");
			}
		}
	
		writer.close();
	}

	private static Map<Gene, double[]> PathsToGenes(Map<Path, double[]> pathsScores) {
	
		Map<Gene, double[]> rtrn=new HashMap<Gene, double[]>();
		for (Path p: pathsScores.keySet()){
			Gene g = p.toGene();
			g.setCountScore(pathsScores.get(p)[4]);
			rtrn.put(g, pathsScores.get(p));
		}
		return rtrn;
	}

	private static Map<Path, Gene> pathsToGenesMap(Map<Path, double[]> pathsScores) {
		Map<Path, Gene> rtrn=new TreeMap<Path, Gene>();
	
		for(Path gene: pathsScores.keySet()){
			double[] p=pathsScores.get(gene);
			//RefSeqGene g=gene.toGene();
			Gene g=gene.toGene(p);
			g.setBedScore(p[4]);
			rtrn.put(gene, g);
	
		}
	
		return rtrn;
	}

	private static Collection<Gene> pathToGenesCollection(Map<Path, double[]> pathsScores) {
		Collection<Gene> rtrn = new TreeSet<Gene>();
	
		for (Path gene : pathsScores.keySet()) {
			double[] p=pathsScores.get(gene);
			//RefSeqGene g=gene.toGene();
			Gene g=gene.toGene(p);
			g.setBedScore(p[4]);
			rtrn.add(g);
		}
	
		return rtrn;
	}

	public boolean hasDataForChromosome(String chr) { return this.chromosomeLengths.containsKey(chr);}

	private boolean overlaps(Alignment read, Annotation current) {
		Annotation t=new Alignments(read.getChr(), read.getStart(), read.getEnd());
		return current.overlaps(t);
	}

	private boolean isPaired() {
		return this.data.isPaired();
	}

	

	/*
	private static Collection<EdgeSourceType> typesToFollow(boolean followPairedEnds) {
		ArrayList<EdgeSourceType> types = new ArrayList<EdgeSourceType>();
		types.add(EdgeSourceType.SPLICED);
		if(followPairedEnds) {
			types.add(EdgeSourceType.PAIRED);
		}
		return types;
	}
	 */
	
	public static void main(String[] args) throws Exception {
		Globals.setHeadless(true);
		System.out.println("Using Version R4.4");
		logger.debug("DEBUG ON");
		ArgumentMap argmap = CLUtil.getParameters(args, usage, "full");
		double lambda = argmap.isPresent("lambda") ? argmap.getDouble("lambda") : 0;
		if ("score".equalsIgnoreCase(argmap.getTask())) {
			String alignmentFile = argmap.getMandatory("alignment");
			boolean useConstituentExons = argmap.containsKey("useConstituentExons");
			boolean useConstituentIntrons = argmap.containsKey("useConstituentIntrons");
			double minMappingQuality = argmap.containsKey("minMappingQuality") ? argmap.getDouble("minMappingQuality") : DEFAULT_MIN_MAPPING_QUALITY;
	
			//Map<String, Collection<RefSeqGene>> annotations = BEDFileParser.loadDataByChr(new File(argmap.getInput()));
			String annotationFile = argmap.getInput();
			BEDFileParser annotationParser =  annotationFile.endsWith(".gtf") || annotationFile.endsWith(".GTF")? new GTFFileParser(annotationFile) : new BEDFileParser(annotationFile);
			Map<String, Collection<Gene>> annotations = null;
			if(useConstituentExons) {
				annotationParser.makeGenes(0.1);
				annotations = annotationParser.toConstituentIsoformMap();
				//annotationParser.writeFullBed("Madegenes.bed");
				BufferedWriter ciw = new BufferedWriter(new FileWriter("constituentExons.bed"));
				for(String chr : annotations.keySet()) {
					Collection<Gene> constituentIsoforms = annotations.get(chr);
					for (Gene g : constituentIsoforms) {
						ciw.write(g.toBED());
						ciw.newLine();
					}
				}
				ciw.close();
			} else if (useConstituentIntrons){
				annotationParser.makeGenes(0.1);
				//annotationParser.writeFullBed("Madegenes.bed");
				annotations = annotationParser.toConstituentIntroformMap();
			}else {
				annotations = annotationParser.toMap();
			}
			Collection<Gene> annotationCollection = new ArrayList<Gene>();
			for(Collection<Gene> chrAnnotations : annotations.values()) {
				annotationCollection.addAll(chrAnnotations);
			}
			//BEDFileParser.writeFullBED("constituentIsoforms.bed", annotationCollection);
			String save = argmap.getOutput();
			String sizes = argmap.get("sizeFile");
			File[] maskFiles  = argmap.isPresent("maskFileDir") ?  new File(argmap.get("maskFileDir")).listFiles(): null ;
			Map<String, Integer> maskFileData = parseMaskFiles(maskFiles);
			logger.info("Minimum mapping quality to count reads: " + minMappingQuality);
			boolean isStranded = argmap.containsKey("stranded");
	
			if(!isStranded) {
				logger.info("Scoring using all reads ");
				AlignmentDataModel Annotation=new GenericAlignmentDataModel(alignmentFile, sizes, false, minMappingQuality);
				Map<Gene, double[]> scores=new TreeMap<Gene, double[]>();				
				runScore(annotations, save, maskFileData, Annotation, scores);
				writeFullBED(save, scores, 1.1);
			} else {
				logger.info("Scoring minus reads");
				AlignmentDataModel Annotation=new GenericAlignmentDataModel(alignmentFile, sizes, false, minMappingQuality);
				Annotation.setNegativeStranded();
				Map<Gene, double[]> scores=new TreeMap<Gene, double[]>();				
				runScore(annotations, save, maskFileData, Annotation, scores);
				writeFullBED(save+".minus", scores, 1.1);
	
				logger.info("Scoring plus reads");
				Annotation.setPositiveStranded();
				scores=new TreeMap<Gene, double[]>();				
				runScore(annotations, save, maskFileData, Annotation, scores);
				writeFullBED(save+".plus", scores, 1.1);
			}
		}	else if ("trim".equalsIgnoreCase(argmap.getTask())) {
			String alignmentFile = argmap.getMandatory("alignment");
			double quantile = argmap.containsKey("quantile") ? argmap.getDouble("quantile") : 0.25;
			Map<String, Collection<Gene>> annotations = BEDFileParser.loadDataByChr(new File(argmap.getInput()));
			String save = argmap.getOutput();
			String sizes = argmap.getMandatory("sizeFile");
			File[] maskFiles  = argmap.isPresent("maskFileDir") ?  new File(argmap.get("maskFileDir")).listFiles(): null ;
			Map<String, Integer> maskFileData = parseMaskFiles(maskFiles);
			AlignmentDataModel Annotation=new GenericAlignmentDataModel(alignmentFile, sizes, false);
			Collection<Gene> trimmed=new TreeSet<Gene>();
			for(String chr : annotations.keySet()) {
				try{
					logger.info("processing " + chr);
					AlignmentDataModelStats alignmentData = new AlignmentDataModelStats(Annotation, maskFileData, chr, false); //Setting lambda greater than one to avoid computing chromosome stats.
					ContinuousDataAlignmentModel data = new ContinuousDataAlignmentModel(alignmentData, maskFileData, 0, 1, null, null);
					data.setTrimEnds(true);
					Collection<Gene> chrAnnotations = annotations.get(chr);
					logger.info("Scanning annotations");
					trimmed.addAll(data.trimEndsForGenes(chrAnnotations, quantile));
				}catch(Exception ex){logger.warn("Skipping "+chr);}
			}
			writeFullBED(save, trimmed);
		} else if ("getIdenticalGappedReadsTranscripts".equalsIgnoreCase(argmap.getTask())) {
			String alignmentFile = argmap.getMandatory("alignment");
			Map<String, Collection<Gene>> annotations = BEDFileParser.loadDataByChr(new File(argmap.getInput()));
			String save = argmap.getOutput();
			String sizes = argmap.getMandatory("sizeFile");
			String name = argmap.getMandatory("name");
			File[] maskFiles  = argmap.isPresent("maskFileDir") ?  new File(argmap.get("maskFileDir")).listFiles(): null ;
			Map<String, Integer> maskFileData = parseMaskFiles(maskFiles);
			GenericAlignmentDataModel Annotation=new GenericAlignmentDataModel(alignmentFile, sizes, false);
			ArrayList<Gene> transcripts =new ArrayList<Gene>();
	
	
			//TODO: support filters and sequence files
			/*
			Sequence chrSequence = null;
			if(chrSequenceFile != null) {
				FastaSequenceIO fsio = new FastaSequenceIO(chrSequenceFile);
	
				List<Sequence> seqs = fsio.loadAll();
	
				if(!seqs.isEmpty()) {
					chrSequence = seqs.get(0);
				} else {
					System.err.println("Sequence for " + chr + " was not found in file " + chrSequenceFile + " continuing without splice site analisys");
				}
				System.err.println("Loaded chromosome Sequence");
			}
			 */
			double total=0; 
			double passed=0;
			double biexon=0;
			for(String chr : annotations.keySet()) {
				//System.err.println("processig " + chr);
				for (Gene gene:annotations.get(chr) )
				{
					Annotation region =new Alignments(gene.getChr(), gene.getStart(), gene.getEnd());
					//boolean pass=Annotation.AreSplicedReadsIdentical(region, Collection<ReadFilter> filters,  Sequence chrSeq);
					transcripts.add(gene); passed++;
					total++;
					if(gene.getNumExons()==2) {biexon++;}
				}
				//System.err.println("done. saving data  to " + save+" for "+chr);
			}	
			writeFullBED(save, transcripts);
			System.out.println(name+"\t"+passed/total+"\t"+passed/biexon);
	
		} else if (argmap.getTask().toUpperCase().contains("PAIREDFILE")){
			String pair1 = argmap.getMandatory("pair1");
			int minMappingQuality = argmap.containsKey("minMappingQuality") ? argmap.getInteger("minMappingQuality") : DEFAULT_MIN_MAPPING_QUALITY;
			boolean forChIP = argmap.containsKey("forChIP");
			String out   = argmap.getOutput();
			boolean isSorted = argmap.containsKey("sorted");
			boolean usePair2Orientation = argmap.containsKey("usePair2Orientation");
	
		}else if ("togff".equalsIgnoreCase(argmap.getTask())) {
			String source = argmap.getMandatory("source");
			boolean toCufflinks = argmap.containsKey("cufflinks");
			String prefix = argmap.get("prefix") ;
			boolean keepids = argmap.containsKey("keepIds");
			String file = argmap.getInput();
			Map<String, Collection<Gene>> data = BEDFileParser.loadDataByChr(new File(file));
	
	
			BufferedWriter bw = argmap.getOutputWriter();
	
			for(String  chr : data.keySet() ) {
				int id = 0;
				for(Gene g : data.get(chr)) {
					if(g.getExtraFields() != null && g.getExtraFields().length > 4) {
						g.addAttribute("RPKM", g.getExtraFields()[4]);
					}
					String name = keepids ? g.getName() : (prefix != null ? prefix  : "") +  "SCRPTR."+g.getChr() + "." + id++;
					g.setName(name);
					bw.write(toCufflinks ? g.toCufflinksGTF(source, name, name, "") :  g.toGTF(source));
				}
			}
			bw.close();
		}else if (argmap.getTask().toUpperCase().contains("CHIP")){
			String out = argmap.getOutput();
			int minMappingQuality = argmap.containsKey("minMappingQuality") ? argmap.getInteger("minMappingQuality") : DEFAULT_MIN_MAPPING_QUALITY;
			String alignmentFile = argmap.getMandatory("alignment");
			int[] windows= getWidths(argmap.getMandatory("windows"));
			int extensionFactor = argmap.containsKey("extensionFactor") ? argmap.getInteger("extensionFactor") : 0;
			String sizes = argmap.get("sizeFile");
			String chr = argmap.get("chr");
			boolean printFullScores = argmap.containsKey("fullScores");
			boolean loadPairsAsFragments = argmap.containsKey("loadPairsAsFragments") || argmap.containsKey("pairedEnd");
	
			//Optional parameters
			boolean findMaxContiguous = argmap.containsKey("findMaxContiguous");
			boolean trimEnds = argmap.containsKey("trim");
			double trimQuantile = argmap.isPresent("trimQuantile") ? argmap.getDouble("trimQuantile") : 0.25;
			int minRemainingLength = argmap.isPresent("minLength") ? argmap.getInteger("minLength") : Statistics.min(windows);
			double alpha = argmap.containsKey("alpha") ? argmap.getDouble("alpha") : .05;
			File[] maskFiles  = argmap.isPresent("maskFileDir") ?  new File(argmap.get("maskFileDir")).listFiles(): null ;
			Map<String, Integer> maskFileData = parseMaskFiles(maskFiles);
	
			ContinuousDataAlignmentModel data = AlignmentUtils.loadAlignmentData(alignmentFile, true, minMappingQuality, true, false, null, loadPairsAsFragments);	
			data.setMaskFileData(maskFileData);
			data.setExtensionFactor(extensionFactor);
			//AlignmentDataModel Annotation=new GenericAlignmentDataModel(alignmentFile, sizes, false, minMappingQuality); //Does this solve the scoring issue?
	
			//logger.info("AlignmentDataModel loaded, initializing model stats");
			//AlignmentDataModelStats alignmentData = new AlignmentDataModelStats(Annotation, maskFileData, chr, lambda, false);
			//ContinuousDataAlignmentModel data = new ContinuousDataAlignmentModel(alignmentData, maskFiles, extensionFactor, 0);
			List<String> chromosomes = new ArrayList<String>();
			if(chr!= null && chr.length() > 0) {
				chromosomes.add(chr);
			} else {
				Map<String, Integer> chrSizes = data.getChromosomeLengths();
				chromosomes = new ArrayList<String>(chrSizes.keySet());
			}
			data.setTrimEnds(trimEnds);
			data.setMinContguousSegmentSize(minRemainingLength);
			data.setTrimQuantile(trimQuantile);
			data.setFindMaxContiguous(findMaxContiguous);
	
			Map<Annotation, double[]> scores = new HashMap<Annotation, double[]>();
			int totalMappedReads = 0;
			for(String workChr : chromosomes) {
				logger.info("Processing chromosome " + workChr + ", Scanning windows");
				totalMappedReads += data.data.getSum(workChr);
				Collection<Annotation> segments = data.scan(windows, alpha, workChr);
				logger.info("Scoring segments");
				scores.putAll(data.scoreSegments(segments, workChr));
				logger.info("Printing results");
			}
			logger.info("Done. Total mapped reads: " + totalMappedReads);
			BEDFileParser.writeSortedBED(out, scores);
			if(printFullScores) {
				BEDFileParser.writeSortedBEDWithScores(out+".scores", scores);
			}
	
		} else if ("scoresegments".equalsIgnoreCase(argmap.getTask())) {
			String alignmentFile = argmap.getMandatory("alignment");
			Map<String, Collection<Annotation>> annotations = BEDFileParser.loadAlignmentDataByChr(new File(argmap.getInput()));
			String save = argmap.getOutput();
			String sizes = argmap.getMandatory("sizeFile");
			File[] maskFiles  = argmap.isPresent("maskFileDir") ?  new File(argmap.get("maskFileDir")).listFiles(): null ;
			double alpha = argmap.isPresent("alpha") ? argmap.getDouble("alpha") : 0.05;
			boolean ignoreAlpha = argmap.isPresent("ignoreAlpha");
			
			if (ignoreAlpha) {
				logger.info("Ignoring alpha - outputting all scores"); 
			} else {
				logger.info("Using alpha = " + alpha);
			}
			Map<String, Integer> maskFileData = parseMaskFiles(maskFiles);
			AlignmentDataModel Annotation=new GenericAlignmentDataModel(alignmentFile, sizes, false);
			Map<Annotation, double[]> scores=new TreeMap<Annotation, double[]>();
			AlignmentDataModelStats alignmentData = new AlignmentDataModelStats(Annotation, maskFileData);
			ContinuousDataAlignmentModel data = new ContinuousDataAlignmentModel(alignmentData, maskFileData, 0, 1, null, null);
			data.alpha = alpha;
			for(String chr : annotations.keySet()) {
				logger.info("processing " + chr);
				Collection<Annotation> chrAnnotations = annotations.get(chr);
				scores.putAll(data.scoreSegments(chrAnnotations, chr, ignoreAlpha));
			}
			BEDFileParser.writeBEDWithScores(save, scores);
		}else if("trimSegments".equalsIgnoreCase(argmap.getTask())) { 
			String alignmentFile = argmap.getMandatory("alignment");
			double trimQuantile = argmap.isPresent("trimQuantile") ? argmap.getDouble("trimQuantile") : 0.25;
			int minRemainingLength = argmap.isPresent("minLength") ? argmap.getInteger("minLength") : 200;
			Map<String, Integer> maskFileData = null;
			if(argmap.isPresent("maskFileDir")) {
				maskFileData=parseMaskFiles(new File(argmap.get("maskFileDir")).listFiles());
			}
			Map<String, Collection<Annotation>> annotations = BEDFileParser.loadAlignmentDataByChr(new File(argmap.getInput()));
			String save = argmap.getOutput();
			String sizes = argmap.get("sizeFile");
			boolean findMaxContiguous = argmap.containsKey("findMaxContiguous");
			AlignmentDataModel Annotation=new GenericAlignmentDataModel(alignmentFile, sizes, false);
			Map<Annotation, double[]> scores=new TreeMap<Annotation, double[]>();
			AlignmentDataModelStats alignmentData = new AlignmentDataModelStats(Annotation, maskFileData);
			ContinuousDataAlignmentModel data = new ContinuousDataAlignmentModel(alignmentData, maskFileData, 0, 1, null, null);
			data.setMinContguousSegmentSize(minRemainingLength);
			data.setTrimQuantile(trimQuantile);
			data.alpha = 1.1;
			for(String chr : annotations.keySet()) {
				logger.info("processig " + chr);
				Collection<Annotation> chrAnnotations = annotations.get(chr);
				if(findMaxContiguous) {
					chrAnnotations = data.findMaxContiguous(chrAnnotations);
				} else {
					chrAnnotations = data.trimEnds(chrAnnotations);
				}
	
				scores.putAll(data.scoreSegments(chrAnnotations, chr));
				BEDFileParser.writeBEDWithScores(save, scores);
			}
		}else if("adjustEnds".equalsIgnoreCase(argmap.getTask())) { 
			String alignmentFile = argmap.getMandatory("alignment");
			double trimQuantile = argmap.isPresent("trimQuantile") ? argmap.getDouble("trimQuantile") : 0.25;
			int minRemainingLength = argmap.isPresent("minLength") ? argmap.getInteger("minLength") : 200;
			Map<String, Integer> maskFileData = null;
			if(argmap.isPresent("maskFileDir")) {
				maskFileData=parseMaskFiles(new File(argmap.get("maskFileDir")).listFiles());
			}
			Map<String, Collection<Annotation>> annotations = BEDFileParser.loadAlignmentDataByChr(new File(argmap.getInput()));
			String save = argmap.getOutput();
			String sizes = argmap.get("sizeFile");
			boolean findMaxContiguous = argmap.containsKey("findMaxContiguous");
			AlignmentDataModel Annotation=new GenericAlignmentDataModel(alignmentFile, sizes, false);
			Map<Annotation, double[]> scores=new TreeMap<Annotation, double[]>();
			AlignmentDataModelStats alignmentData = new AlignmentDataModelStats(Annotation, maskFileData);
			ContinuousDataAlignmentModel data = new ContinuousDataAlignmentModel(alignmentData, maskFileData, 0, 1, null, null);
			data.setMinContguousSegmentSize(minRemainingLength);
			data.setTrimQuantile(trimQuantile);
			data.alpha = 1.1;
			for(String chr : annotations.keySet()) {
				logger.info("processig " + chr);
				Collection<Annotation> chrAnnotations = annotations.get(chr);
				if(findMaxContiguous) {
					chrAnnotations = data.findMaxContiguous(chrAnnotations);
				} else {
					chrAnnotations = data.trimEnds(chrAnnotations);
				}
	
				scores.putAll(data.scoreSegments(chrAnnotations, chr));
				BEDFileParser.writeBEDWithScores(save, scores);
			}
		}
			
		else if("polyAseq".equalsIgnoreCase(argmap.getTask())){
			GenericAlignmentDataModel Annotation=new GenericAlignmentDataModel("temp.bam", "sizes", 5);
			Annotation.setNegativeStranded();
			ContinuousDataAlignmentModel data = new ContinuousDataAlignmentModel(Annotation);
			Annotation al = new Alignments("chrX","1","2600737");
			double XX = data.count(al);
			System.out.println("Negative strand: "+XX);
			Annotation.setPositiveStranded();
			data = new ContinuousDataAlignmentModel(Annotation);
			XX = data.count(al);
			System.out.println("Positive strand: "+XX);
			
		}
		else{
			System.err.println("Invalid task " + argmap.getTask() + "\n" + usage);
		}
	}

	public void resetTreeCache() {
		data.resetTreeCache();
	}

	private static void runScore(Map<String, Collection<Gene>> annotations, String save,
			Map<String, Integer> maskFileData, AlignmentDataModel Annotation,	Map<Gene, double[]> scores) throws IOException {
		AlignmentDataModelStats alignmentData = new AlignmentDataModelStats(Annotation, maskFileData);
		ContinuousDataAlignmentModel data = new ContinuousDataAlignmentModel(alignmentData, maskFileData, 0, 1, null, null);
		for(String chr : annotations.keySet()) {
			logger.info("processig " + chr);
			Collection<Gene> chrAnnotations = annotations.get(chr);
			scores.putAll(data.scoreGenes(chrAnnotations, chr));
			//System.err.print("done. saving data  to " + save+"."+chr);
		}
	}

	private static int sum(Collection<Annotation> maskedRegions){
		int sum=0;
		for(Annotation align: maskedRegions){sum+=align.getSize();}
		return sum;
	}
	
	

	public double CV(Gene gene, IntervalTree<Alignment> tree) throws IOException{
		List<Double> vals=getNormalizedCountsPerBp(gene, tree);
		double avg=Statistics.average(vals);
		double sd=Statistics.stdev(vals);
		return sd/avg;
	}

	private double[] l2a(List<Double> list){
		double[] rtrn=new double[list.size()];
	
		int i=0;
		for(Double val: list){rtrn[i++]=val;}
	
		return rtrn;
	}

	private double[] subtract(double[] array, double factor){
		double[] rtrn=new double[array.length];
	
		for(int i=0; i<array.length; i++){
			rtrn[i]=array[i]-factor;
		}
	
		return rtrn;
	}


	

	private Collection<? extends Annotation> toExons( Alignment record){
		return record.getFragment(null).iterator().next().getBlocks();
	}

	private Collection<Annotation> dedup(IntervalTree<Annotation> allSegments){
		Set<Annotation> set=new TreeSet<Annotation>();
	
		Iterator<Node<Annotation>> iter=allSegments.iterator();
	
		while(iter.hasNext()){
			Annotation align=iter.next().getValue();
			set.add(align);
		}	
	
		return set;
	}

	private Collection<Annotation> findMaxContiguous(Collection<Annotation> windows) throws IOException{
		Set<Annotation> rtrn = new TreeSet<Annotation>();
	
		//loop through each alignment and compute number of reads for the ends
		Map<Annotation, double[]> counts=this.getDataForAnnotation(windows);
		for(Annotation align: windows){
			Annotation trunc=collapseMaxContiguous(align, counts.get(align));
			if(trunc!=null){rtrn.add(trunc);}
		}
	
		return rtrn;
	}


	public Collection<Annotation> permuteRegions(Annotation align) throws IOException{
		Collection<Annotation> rtrn=new TreeSet<Annotation>();
		//get records overlapping this region
		Iterator<Annotation> iter=data.getExonAnnotationOverlappingRegion(align);
	
		int regionSize=align.getSize();
	
		while(iter.hasNext()){
			Annotation record=iter.next();
			//randomly place in the interval
			int start=align.getStart()+new Double(Math.random()*(regionSize)).intValue();
			Annotation r=new Alignments(align.getChr(), start, start+(record.getSize()));
			r.setOrientation(record.getStrand());
			rtrn.add(r);
		}
	
		return rtrn;
	}

	public Collection<Annotation> rawRegions(Annotation align) throws IOException{
		Collection<Annotation>  rtrn=new TreeSet<Annotation> ();
		//get records overlapping this region
		Iterator<Annotation> iter=data.getExonAnnotationOverlappingRegion(align);
	
		while(iter.hasNext()){
			Annotation record=iter.next();
			rtrn.add(record);
		}
	
		return rtrn;
	}

	private static Collection<Gene> filter(Map<Gene, double[]> scores, double alpha2) {
		Collection<Gene> rtrn=new TreeSet<Gene>();
	
		for(Gene gene: scores.keySet()){
			double[] p=scores.get(gene);
			if(p[0]<alpha2){
				rtrn.add(gene);
			}
		}
	
		return rtrn;
	}

	public static Collection<Gene> filterPaths(Map<Path, double[]> scores, double alpha2) {
		Collection<Gene> rtrn=new TreeSet<Gene>();
	
		/*
		 * [0] = P-value
		 * [1] = enrichment 
		 * [2] = score (counts?)
		 * [3] = average coverage
		 * [4] = rpkm
		 * [5] = lambda
		 * [6] = size
		 * [7] = localPvalue
		 */
		for(Path gene: scores.keySet()){
			double[] p=scores.get(gene);
			
			//RefSeqGene g=gene.toGene();
			Gene g=gene.toGene(p);
			//System.out.println("For "+g.toBED()+": ");
			//System.out.println("P-value: "+p[0] +", RPKM: "+p[4]+", Score: "+p[2]);
			g.setBedScore(p[4]);
			if(g.getNumExons()>1){
				if(p[0]<alpha2){rtrn.add(g);}
			}
			else if(p.length>7){
	
				if(p[0]<alpha2 && p[7]<alpha2){rtrn.add(g);}
			}
			else{
	
				if(p[0]<alpha2){rtrn.add(g);}
			}
	
		}
	
		return rtrn;
	}
	
	// For AlignmentCollection interface - should be moved to Generic model eventually. JE
	public int getBasesCovered(Annotation region) throws IOException {
		return data.getData().getBasesCovered(region, extensionFactor); 
	}
	
	// For AlignmentCollection interface - should be moved to Generic model eventually. JE
	public int getBasesCovered(Annotation region, int EF) throws IOException {
		return data.getData().getBasesCovered(region, EF); 
	}

	private static Collection<Path> refSeqGeneToPath( IntervalTree<GeneWithIsoforms> chrTree, ChromosomeWithBubblesJGraphT graph) {
		Collection<Path> rtrn = new LinkedList<Path>();
		Iterator<GeneWithIsoforms> it =chrTree.valueIterator();
		while(it.hasNext()){
			Collection <Gene> allIso=it.next().getAllIsoforms();
			for (Gene iso: allIso){
				Path p=new Path(iso, graph);
				rtrn.add(p);
			}
		}
		return rtrn;
	}

	public static Map<String, Integer> parseMaskFiles(File[] files) throws IOException{
		Map<String, Integer> rtrn=new TreeMap<String, Integer>();
		if(files==null){return rtrn;}
	
		for(int i=0; i<files.length; i++){
			String chr=files[i].getName().split("\\.")[0];
			Collection<Annotation> data = BEDFileParser.loadAlignmentData(files[i]);
			rtrn.put(chr,sum(data));
		}	
		return rtrn;
	}





	//static String usage="version=v6 \n args[0]=BAM File \n args[1]=mask file \n args[2]=save file \n args[3]=sizes \n args[4]=window size \n args[5]=trim ends?";
	static String usage="\nParameters \n -alignment <Alignment file in BAM, SAM or Alignemnt format> \n -maskFileDir <Mask File directory> \n -out <Output file name>"+ 
			"\n -sizeFile <Chromosome size file> \n -chr <Chromsomosome to segment> \n -chrSequence <Necessary to filter spliced reads by splice site information. Notice that this is only compatible with region files that contain regions of only one chromosome> "+
			"\n Optional arguments: \n -windows <Comma separated list of windows to evaluate defaults to contiguous regions of coverage>\n -trim <Include this flag if trimming of the ends of windows based on read coverage  is desired this is expensive> \n -alpha <Desired FDR>" +
			"\n  -dontFilterCanonicalSplice" +
			"\n -start <To segment only a subregion of the chromosome include its start> -end <To segment only a subregion of the chromosome include its end> " + 
			"\n -minSpliceSupport <Minimum count to support splice reads, default is 1> \n-minSpliceFrequency <When there are more than one splice junction, junctions that account for less than the specified portion of junctions are ignored>\n -pairedEnd <Paired end alignment files> -strandSpecificReads <Strand specific alignment file>\n\t-scoreRegions <Full BED to score> -upWeightSplices -lambda <If a prior background expectation for number of reads per base exists> -exons <BED file of exons> -introns <Introns and counts>" +
			"\n\nTask: AddPairs -  Uses a paired end alignment to tune graph \n\t-in <Graph in .dot format. Standard input is assumed> \n\t-pairedEnd <Paired end information (as in previous task), in single line BED format>\n\t -maskFileDir <Directory containing mask files for the genome> \n\t-chr <Chromosome (only a chromosome at a time is supported at this point)> \n\t-sizeFile <Chromosome size file> \n\t-out <Output file name>"+
			"\n\nTask: fastScore -  Computes several expression related scores for a set of annotations using a the graph .dot file \n\t-in <chr.dot file> \n\t-annotations <BED file with annotation to score>\n\t -chr <chr e.g: chrZ>\n\t -alpha <optional>\n\t -out\n" +
			"\n\nTask: score -  Computes several expression related scores for a set of annotations -in <Full BED file with annotations to score> \n\t-alignment <Alignment file in BAM, SAM or Alignemnt format> \n\t-sizeFile <Chromosome size file> \n\t-out <Output file name> \n\t -maskFileDir <Mask File directory>\n\t -useConstituentExons <For each gene a set of contituent exons will be chosen and scored>"+
			"\n\nTask: extractDot - Extracts a graph for the specified region \n\t-in <Dot file from a previous Scripture run \n\t-chr <Chromosome> \n\t-start <Start of region> \n\t-end<End of region> \n\t-out <output file>"+
			"\n\nTask: getIdenticalGappedReadsTranscripts -  Report all transcripts that ALL their introns are spanned by identical gapped reads -in <Full BED file with annotations to score> \n\t-alignment <Alignment file in BAM, SAM or Alignemnt format> \n\t-sizeFile <Chromosome size file> \n\t-out <Output file name> \n\t -maskFileDir <Mask File directory>"+
			"\n\nTask: makePairedFile Makes a paired end alignment file from two sets of independtly aligned left and right ends, ideally the files should be name-sorted \n\t-pair1 <First pair Annotation> \n\t-pair2 <Second pair Annotation> \n\t-out <output consolidated paired end alignment> \n\t-sorted  <Include this flag if the data is already read name sorted, ideally both input files should be sorted by read name using unix sort for example> \n\t-usePair2Orientation <If the second paired rather than the first should be used to orient insert like for dUTP libraries>\n\t-forChIP <If the alignment if for ChIP rather than RNAseq then Ms will be used instead of Ns>" + 
			"\n\nTask: chipScan - Segment the genome assuming contiguous data. Similar to the default task but optimized for contiguous data. \n -alignment <Alignment file in BAM, SAM or Alignemnt format> \n -extensionFactor <Extend reads by this factor (defaults to 0)> \n -maskFileDir <Mask File directory> \n -out <Output file name>"+ 
			"\n -chr <Chromosome to segment>\n -sizeFile <Chromosome size file> \n -windows <Comma separated list of windows to evaluate defaults to contiguous regions of coverage> \n Optional arguments:\n -findMaxContiguous <Each significant window is trimmed by finding the contiguous sub region with coverage over a predefined threshold> -trim <Include this flag if trimming of the ends of windows based on read coverage  is desired this is expensive> \n -alpha <Desired FDR>" +
			"\n\nTask: trim -  Trims end of transcripts by removing all bases whose coverage is below the specified quantile of transcript expression -in <Full BED file with annotations to trim> \n\t-alignment <Alignment file in BAM, SAM or Alignemnt format> \n\t-sizeFile <Chromosome size file> \n\t-out <Output file name> \n\t -maskFileDir <Mask File directory>\n\t-quantile <Coverage quantile below which end bases should be trimmed>"+
			"\n\nTask: trimSegments -  Trims end of continuous segments by removing all bases whose coverage is below the specified quantile of transcript expression -in <Full BED file with annotations to trim> \n\t-alignment <Alignment file in BAM, SAM or Alignemnt format> \n\t[-sizeFile <Chromosome size file>] \n\t-out <Output file name> \n\t -maskFileDir <Mask File directory>\n\t-quantile <Coverage quantile below which end bases should be trimmed> \n\t\t-findMaxContiguous <To break up segments that have peak/valley shapes>"+
			"\n\nTask: adjustEnds -  Takes an annotation set and adjust transcript ends based on the given alignment \n\t-in <Full BED file with annotations to trim> \n\t-alignment <Alignment file in BAM, SAM or Alignemnt format> \n\t[-sizeFile <Chromosome size file>] \n\t-out <Output file name> \n\t -maskFileDir <Mask File directory>\n\t=trimQuantile <Coverage quantile below which end bases should be trimmed> \n\t\t-findMaxContiguous <To break up segments that have peak/valley shapes>"+
			"\n";




}
