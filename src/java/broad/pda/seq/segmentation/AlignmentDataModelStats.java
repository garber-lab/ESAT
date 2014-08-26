package broad.pda.seq.segmentation;

import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import net.sf.samtools.util.CloseableIterator;
import nextgen.core.alignment.Alignment;
import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Gene;

import org.apache.commons.collections15.Predicate;
import org.apache.log4j.Logger;

import broad.core.annotation.LightweightGenomicAnnotation;
import broad.core.datastructures.IntervalTree;
import broad.core.math.ScanStatistics;
import broad.core.sequence.Sequence;
import broad.core.util.CollapseByIntersection;
import broad.pda.datastructures.Alignments;
import broad.pda.seq.alignment.StoredAlignmentStats;


public class AlignmentDataModelStats {
	private static final Logger logger =  Logger.getLogger(AlignmentDataModelStats.class.getName());
	private Map<String, Double> numMarkers; // chr size minus masked regions
	private Map<String, Double> numberOfReads;
	private Map<String, Integer> maskedRegions;
	private Map<String, Double> lambda;
	private Map<String, Double> rpkmConstant;
	private double globalRPKMConstant;
	private AlignmentDataModel data;
	private boolean isPaired;
	
	
	public AlignmentDataModelStats() {
		super();
		initializeVariables();
	}
	
	/*public AlignmentDataModelStats(AlignmentDataModel data, Map<String, Integer> maskedRegions, String chr) throws IOException {
		this(data, maskedRegions, chr, 0, true);		
	}*/
	
	
	public AlignmentDataModelStats(AlignmentDataModel data, Map<String, Integer> maskedRegions, String chr) throws IOException {
		this.data = data;
		this.maskedRegions = maskedRegions;
		initializeVariables();
		//computeGlobalStats(chr);
		Map<String,Integer > chromosomeLengths = data.getChromosomeLengths();
		computeGlobalStats(chr);
		double totalReads=this.getNumberOfReads(chr);
		if(totalReads > 0) {
			globalRPKMConstant = Math.pow(10,9) / totalReads;
		}
	}
	
	public AlignmentDataModelStats(AlignmentDataModel data, Map<String, Integer> maskedRegions, String chr, double lambda, boolean loadChromosomeStats) throws IOException {
		this.data = data;
		this.maskedRegions = maskedRegions;
		initializeVariables();
		//computeGlobalStats(chr);
		Map<String,Integer > chromosomeLengths = data.getChromosomeLengths();
		int totalReads = 0;
		for(String chromosome : chromosomeLengths.keySet()) {
			if(lambda > 0 && chr != null && chr.equals(chromosome)) {
				this.lambda.put(chr, lambda);
				if(maskedRegions!=null && maskedRegions.containsKey(chr)){numMarkers.put(chr, (double) (data.getChromosomeLength(chr) - maskedRegions.get(chr)));}
				else{numMarkers.put(chr, (double) (data.getChromosomeLength(chr)));}
				numberOfReads.put(chr, lambda * numMarkers.get(chr));
				rpkmConstant.put(chr, Math.pow(10,9)/numberOfReads.get(chr));
				
			} else if (loadChromosomeStats){
				computeGlobalStats(chromosome);
			}
			if(numberOfReads.get(chromosome) != null) {
				totalReads += numberOfReads.get(chromosome);
			}
		}
		if(totalReads > 0) {
			globalRPKMConstant = Math.pow(10,9) / totalReads;
		}
		CloseableIterator<Alignment> readIt = data.getReadIterator();
		Alignment alignment = readIt.next();
		isPaired = alignment.isPaired();
		logger.debug("The alignment is paired: " + isPaired);
		readIt.close();
		//System.err.println("Lambda:  " + (chr!=null ? this.lambda.get(chr) : "null chromosome"));
	}
	
	public AlignmentDataModelStats(AlignmentDataModel data, Map<String, Integer> maskedRegions) throws IOException {
		this.data = data;
		init(maskedRegions, 0, true);		
	}
	
	public AlignmentDataModelStats(AlignmentDataModel data, Map<String, Integer> maskedRegions, Map<String, Double> readsPerChromosome) throws IOException {
		this.data = data;
		init(maskedRegions, 0, true, readsPerChromosome);		
	}
	
	public AlignmentDataModelStats(AlignmentDataModel data, Map<String, Integer> maskedRegions, String chr, boolean loadChromosomeStats) throws IOException {
		this(data, maskedRegions, 0, loadChromosomeStats);
	}
	
	public AlignmentDataModelStats(AlignmentDataModel data, Map<String, Integer> maskedRegions, double lambda) throws IOException {
		this.data = data;
		init( maskedRegions, lambda, true);
	}
	
	public AlignmentDataModelStats(AlignmentDataModel data, Map<String, Integer> maskedRegions, double lambda, boolean loadChromosomeStats) throws IOException {
		this.data = data;
		init(maskedRegions, lambda, loadChromosomeStats);
	}
	
	public AlignmentDataModelStats(AlignmentDataModel data, Map<String, Integer> maskedRegions, String chr, double lambda) throws IOException {
		this(data, maskedRegions, chr, lambda, true);
	}
	

	
	public AlignmentDataModelStats(AlignmentDataModel data, String chr) throws IOException {
		this(data, null, chr);
		
	}
	
	public AlignmentDataModelStats(AlignmentDataModel data) throws IOException {
		this.data = data;
		init(null, 0, true);
		
	}

	public void init( Map<String, Integer> maskedRegions, double lambda, boolean loadChromosomeStats) throws IOException{
		this.maskedRegions = maskedRegions;
		initializeVariables();
		//computeGlobalStats(chr);
		Map<String,Integer > chromosomeLengths = data.getChromosomeLengths();
		int totalReads = 0;
		if(lambda > 0) {
			for(String chr : chromosomeLengths.keySet()) {
				this.lambda.put(chr, lambda);
				if(maskedRegions!=null && maskedRegions.containsKey(chr)){
					numMarkers.put(chr, (double) (data.getChromosomeLength(chr) - maskedRegions.get(chr)));
				}else{
					numMarkers.put(chr, (double) (data.getChromosomeLength(chr)));
				}
				numberOfReads.put(chr, lambda * numMarkers.get(chr));
				totalReads += numberOfReads.get(chr);
				rpkmConstant.put(chr, Math.pow(10,9)/numberOfReads.get(chr));
			}
			if(totalReads > 0) {
				globalRPKMConstant = Math.pow(10,9) / totalReads;
			}
		}
		if(loadChromosomeStats) {
			computeGlobalStats();
		}
		CloseableIterator<Alignment> readIt = data.getReadIterator();
		Alignment alignment = readIt.next();
		isPaired = alignment.isPaired();
		logger.info("The alignment is paired: " + isPaired);
		readIt.close();
	}
	
	public void init( Map<String, Integer> maskedRegions, double lambda, boolean loadChromosomeStats, Map<String, Double> readsByChr) throws IOException{
		this.maskedRegions = maskedRegions;
		initializeVariables();
		//computeGlobalStats(chr);
		Map<String,Integer > chromosomeLengths = data.getChromosomeLengths();
		int totalReads = 0;
		if(lambda > 0) {
			for(String chr : chromosomeLengths.keySet()) {
				this.lambda.put(chr, lambda);
				if(maskedRegions!=null && maskedRegions.containsKey(chr)){
					numMarkers.put(chr, (double) (data.getChromosomeLength(chr) - maskedRegions.get(chr)));
				}else{
					numMarkers.put(chr, (double) (data.getChromosomeLength(chr)));
				}
				numberOfReads.put(chr, lambda * numMarkers.get(chr));
				totalReads += numberOfReads.get(chr);
				rpkmConstant.put(chr, Math.pow(10,9)/numberOfReads.get(chr));
			}
			if(totalReads > 0) {
				globalRPKMConstant = Math.pow(10,9) / totalReads;
			}
		}
		if(loadChromosomeStats) {
			computeGlobalStats(readsByChr);
		}
		CloseableIterator<Alignment> readIt = data.getReadIterator();
		Alignment alignment = readIt.next();
		isPaired = alignment.isPaired();
		System.err.println("The alignment is paired: " + isPaired);
		readIt.close();
	}

	private void initializeVariables(){
		numberOfReads=new TreeMap<String, Double>();
		numMarkers=new TreeMap<String, Double>();
		lambda=new TreeMap<String, Double>();
		rpkmConstant = new HashMap<String, Double>();
		maskedRegions = new HashMap<String, Integer>();
	}
	
	private void computeGlobalStats() throws IOException{
		Map<String, Integer> chromosomeLengths = data.getChromosomeLengths();
		int totalReads = 0;
		for(String chromosome : chromosomeLengths.keySet()) {
			computeGlobalStats(chromosome);
			totalReads += numberOfReads.get(chromosome);
		}
		globalRPKMConstant = Math.pow(10,9) / totalReads;
	}
	
	private void computeGlobalStats(Map<String, Double> numberReads) throws IOException{
		Map<String, Integer> chromosomeLengths = data.getChromosomeLengths();
		int totalReads = 0;
		for(String chromosome : chromosomeLengths.keySet()) {
			computeGlobalStats(chromosome, numberReads);
			totalReads += numberOfReads.get(chromosome);
		}
		globalRPKMConstant = Math.pow(10,9) / totalReads;
	}
	
	private void computeGlobalStats(String chr) throws IOException{
		logger.info("Computing aln global stats for " + chr + " - length: " + data.getChromosomeLengths().get(chr));
		if(chr!=null && data.getChromosomeLengths().containsKey(chr)){
			double [] dataPoints = computeDataStats(chr);
			lambda.put(chr, dataPoints[2]);
			numberOfReads.put(chr, dataPoints[0]);
			numMarkers.put(chr, dataPoints[1]);
			rpkmConstant.put(chr, Math.pow(10,9)/dataPoints[0]);
		}
		
	}
	
	private void computeGlobalStats(String chr, Map<String, Double> numberReads) throws IOException{
		System.err.println("Computing aln global stats for " + chr + " - length: " + data.getChromosomeLengths().get(chr));
		if(chr!=null && data.getChromosomeLengths().containsKey(chr)){
			double [] dataPoints = computeDataStats(chr, numberReads);
			lambda.put(chr, dataPoints[2]);
			numberOfReads.put(chr, dataPoints[0]);
			numMarkers.put(chr, dataPoints[1]);
			rpkmConstant.put(chr, Math.pow(10,9)/dataPoints[0]);
		}
		
	}
	
	private double [] computeDataStats(String chr) throws IOException {
		//System.err.print("computing ...  ");
		double numberOfReads = getNumberOfReadsByChr(chr);
		double numMarkers=data.getChromosomeLength(chr);
		
		if(maskedRegions!=null && maskedRegions.containsKey(chr)){
			numMarkers=numMarkers - maskedRegions.get(chr); //need to add support for unalignable bases and remove them from the count
		}
		
		double lambda = numberOfReads/numMarkers;
		logger.info("chromosome "+ chr +" Lambda: " + lambda + " reads " + numberOfReads + " markers " + numMarkers);
		double [] dataPoints = {numberOfReads, numMarkers, lambda};
		return dataPoints;
	}
	
	private double [] computeDataStats(String chr, Map<String, Double> numberReads) throws IOException {
		//System.err.print("computing ...  ");
		double numberOfReads = numberReads.get(chr);
		double numMarkers=data.getChromosomeLength(chr);
		
		if(maskedRegions!=null && maskedRegions.containsKey(chr)){
			numMarkers=numMarkers - maskedRegions.get(chr); //need to add support for unalignable bases and remove them from the count
		}
		
		double lambda = numberOfReads/numMarkers;
		System.err.println("chromosome "+ chr +" Lambda: " + lambda + " reads " + numberOfReads + " markers " + numMarkers);
		double [] dataPoints = {numberOfReads, numMarkers, lambda};
		return dataPoints;
	}
	
	public AlignmentDataModel getData() { return data;}
	
	public double getGlobalRPKMConstant() { return globalRPKMConstant;}
	
	public Map<String, Double> getNumberOfReads(){return this.numberOfReads;}

	public double getLambda(String chr) throws IOException{
		if(!lambda.containsKey(chr) || lambda.get(chr)==null){
			logger.info("computing from scratch "+chr);
			computeGlobalStats(chr);
		}
		return lambda.get(chr);
	}
	public double getSum(String chr) throws IOException{
		if(!numberOfReads.containsKey(chr) || numberOfReads.get(chr)==null){this.computeGlobalStats(chr);}
		return this.numberOfReads.get(chr);
	}	
	public double getNumberOfReads(String chr) throws IOException{return getSum(chr);}
	
	public double getNumberOfReadsByChr(String chr) throws IOException {
		Annotation chrRegion = null;
		if(data.getChromosomeLengths().containsKey(chr)) {
			chrRegion = new Alignments(chr, 0, data.getChromosomeLength(chr));
		}
		
		double count = 0;
		if(chrRegion != null) {
			boolean isStranded = data.isStranded();
			boolean isNegativeStranded = data.isNegativeStranded();
			data.unsetStranded();
			count =  data.getCountsPerAlignment(chrRegion, 0);
			if (isStranded) {
				if(isNegativeStranded) {
					data.setNegativeStranded();
				} else {
					data.setPositiveStranded();
				}
			}
		}
		return count;	
	}

	public int getTotalReads() throws IOException {
		int total=0;
		
		for(String chr: data.getChromosomeLengths().keySet()){
			total+=this.getSum(chr);
		}
		
		return (total);
	}

	private int getTotalSplicedReads() throws IOException{
		int total=0;
		
		for(String chr: data.getChromosomeLengths().keySet()){
			Annotation chrRegion=new Alignments(chr, 0, data.getChromosomeLength(chr));
			total+=data.getNumberOfSplicedReads(chrRegion);
		}
		
		return (total);
	}

	public  double getNumberMarkers(String chr) throws IOException{
		if(!numMarkers.containsKey(chr) || numMarkers.get(chr)==null){this.computeGlobalStats(chr);}
		return this.numMarkers.get(chr);
	}

	public Map<String, Integer> getChromosomeLengths() {
		return data.getChromosomeLengths();
	}
	
	

	public int[] getGeneCounts(Collection<Gene> genes,  int extensionFactor) throws IOException{
		int exon=0;
		int intron=0;
		int UTR5=0;
		int UTR3=0;
		for(Gene gene: genes){
			GeneCounts geneCounts = data.getGeneCounts(gene, extensionFactor);
	
			exon+=geneCounts.getExonicCounts();
			intron+= geneCounts.getdoubleronicCounts();
			UTR3+=geneCounts.getUtr3Counts();
			UTR5+=geneCounts.getUtr5Counts();
		}
		
		int total=getTotalReads();
		int splicedReads=getTotalSplicedReads();
		int intergenic=total-(exon+intron+UTR5+UTR3);
		int[] rtrn={exon, intron, UTR5, UTR3, intergenic, total, splicedReads};
		return rtrn;
	}
	
	/**
	 * Get read counts mapping to exons, introns, 5'UTRs and 3'UTRs in gene set of interest
	 * @param genes
	 * @param extensionFactor
	 * @return An int array consisting of counts for exons, introns, 5'UTRs, 3'UTRs
	 * @throws IOException
	 */
	public int[] getGeneCountsWithoutDatasetTotals(Collection<Gene> genes,  int extensionFactor) throws IOException{
		int exon=0;
		int intron=0;
		int UTR5=0;
		int UTR3=0;
		for(Gene gene: genes){
			GeneCounts geneCounts = data.getGeneCounts(gene, extensionFactor);
	
			exon+=geneCounts.getExonicCounts();
			intron+= geneCounts.getdoubleronicCounts();
			UTR3+=geneCounts.getUtr3Counts();
			UTR5+=geneCounts.getUtr5Counts();
		}
		

		int[] rtrn={exon, intron, UTR5, UTR3};
		return rtrn;
	}

	public Map<Gene, double[]> getGeneExpression(Collection<Gene> genes, int extensionFactor) throws IOException {
		Map<Gene, double[]> rtrn=new TreeMap<Gene, double[]> ();
		
		String chr="";
		IntervalTree<Alignment> tree=null;
		for(Gene gene: genes){
			if(!chr.equalsIgnoreCase(gene.getChr())){
				chr=gene.getChr(); 
				tree=data.getIntervalTree(chr, 0, data.getChromosomeLength(chr));
			}
			double[] scanP= scanPRate(gene, tree, extensionFactor);
			rtrn.put(gene, scanP);
		}
		
		return rtrn;
	}

	public int getCounts(Collection<Annotation> Annotation, int extensionFactor) throws IOException {
		// JE 9/26/12 Use caching in GenericAlignmentDataModel
		
		int rtrn = 0;
		for (Annotation a : Annotation) {
			rtrn += data.getCountsPerAlignment(a, extensionFactor);  // converting double to int here ...
		}
		return rtrn;
	}

	public double getCountsPerAlignment(Annotation align, IntervalTree<Alignment> tree, int extensionFactor) {
		return data.getCountsPerAlignment(align, tree, extensionFactor);
	}

	public double getCountsPerAlignment(Gene window,IntervalTree<Alignment> tree, int extensionFactor) {
		return data.getCountsPerAlignment(window, tree, extensionFactor);
	}

	public double getCountsPerAlignment(Annotation align,Map<String, IntervalTree<Annotation>> exonTree,IntervalTree<Alignment> tree, int extensionFactor) {
		return data.getCountsPerAlignment(align, exonTree, tree, extensionFactor);
	}

	public double getRPKMConstant(String chr) {
		return globalRPKMConstant == 0 ? rpkmConstant.get(chr) : globalRPKMConstant;
	}

	public Map<Annotation, Integer> getSplicedReads(Annotation chrRegion,int i, int j) throws IOException {
		return data.getSplicedReads(chrRegion, i, j);
	}

	public Map<Annotation, Integer> getSplicedReads(Annotation region,List<Predicate<Annotation>> filters, int minNumIntrons) throws IOException {
		return data.getSplicedReads(region, filters, minNumIntrons);
	}

	public Map<Annotation, Integer> getSplicedReads(Annotation region,int minIntronSize, int maxIntronSize, int minNumberOfSplices, Sequence chrSequence) throws IOException {
		return data.getSplicedReads(region, minIntronSize, maxIntronSize, minNumberOfSplices, chrSequence);
	}

	public IntervalTree<Alignment> getIntervalTree(LightweightGenomicAnnotation region) throws IOException {
		return getIntervalTree(region.getChromosome(), region.getStart(), region.getEnd());
	}

	public IntervalTree<Alignment> getIntervalTree(String chr, int start, int end) throws IOException {
		return data.getIntervalTree(chr, start, end);
	}

	public IntervalTree<Alignment> getIntervalTreeCached(String chr, int start,int end) throws IOException {
		return data.getIntervalTreeCached(chr, start, end);
	}

	public IntervalTree<Alignment> getIntervalTreeTruncatedCached(String chromosome, int start, int end) throws IOException {
		return data.getIntervalTreeTruncatedCached(chromosome, start, end);
	}

	public double getScorePerAlignmentFromCache(Gene window,IntervalTree<Alignment> chunkAlignmentTree, int extensionFactor) {
		return data.getScorePerAlignmentFromCache(window, chunkAlignmentTree, extensionFactor);
	}

	public int getChunkStart() {
		return data.getChunkStart();
	}

	public int getChunkEnd() {
		return data.getChunkEnd();
	}

	public double getSpliceWeightFactor() {
		return data.getSpliceWeightFactor();
	}

	public Iterator<Annotation> getExonAnnotationOverlappingRegion(Annotation align) throws IOException {
		return data.getExonAnnotationOverlappingRegion(align);
	}

	public CloseableIterator<Alignment> getAnnotationOverlappingRegion(	Annotation Annotation) throws IOException {
		
		return data.getAnnotationOverlappingRegion(Annotation);
	}

	public double getCountsOfUniqueOverlappers(Annotation startPosition,
			Annotation window, IntervalTree<Alignment> chunkAlignmentTree,
			int extensionFactor) {
		return data.getCountsOfUniqueOverlappers(startPosition, window, chunkAlignmentTree, extensionFactor);
	}

	public CloseableIterator<Alignment> getReadIterator() {
		return data.getReadIterator();
	}

	public CloseableIterator<Alignment> getReadIterator(Annotation region) throws IOException{
		return data.getReadIterator(region);
	}

	public double getMinimumMappingQuality() { return data.getMinimumMappingQuality();}

	public void setAlignmentDataModel(AlignmentDataModel data) {
		this.data = data;
		
	}

	public void setCountNormalizer(ReadCountNormalizer normalizer) {
		data.setNormalizationStrategy(normalizer);
	}

	public void setLambda(String chr, double numberOfReads) {
		double numMarkers=data.getChromosomeLength(chr);
		if(maskedRegions!=null && maskedRegions.containsKey(chr)){
			numMarkers=numMarkers - maskedRegions.get(chr); //need to add support for unalignable bases and remove them from the count
		}
		
		double lambda = numberOfReads/numMarkers;
		this.lambda.put(chr, lambda);
	}

	public void setNumberOfReads(String chr, double numberOfReads2) {
		this.numberOfReads.put(chr, numberOfReads2);
	}

	public void setMinimumMappingQuality(double minMapQual) {this.data.setMinimumMappingQuality(minMapQual);}

	public void setNumberOfMarkers(String chr, double numMarkers2) {
		this.numMarkers.put(chr, numMarkers2);
	}

	public void setChunkSize(int chunkSize) {
		data.setChunkSize(chunkSize);
	}

	public void resetTreeCache() {
		data.resetTreeCache();		
	}

	public void resetTreeCache(String chr) {
		data.resetTreeCache(chr);	
	}

	public double[] scanPRate(Gene gene, int extensionFactor) throws IOException {
		String chr=gene.getChr();
		double sum=0; // number of reads mapping to gene
		int count=0; // gene size (total length of exons)
		
		Set<? extends Annotation> exons=gene.getExonSet();
		for(Annotation exon: exons){
			sum+=data.getCountsPerAlignment(exon, extensionFactor);
			count+=exon.getSize();
		}
		double avgCoverage = sum/(double)count;
		double enrich=(avgCoverage)/getLambda(chr);
		//double[] rtrn={calculatePVal(new Double(sum).intValue(), getLambda(chr), count, getNumberMarkers(chr)), enrich, sum, avgCoverage, avgCoverage *rpkmConstant.get(chr)};
		double[] rtrn={ScanStatistics.calculatePVal(new Double(sum).intValue(), getLambda(chr), count, getNumberMarkers(chr)), enrich, sum, avgCoverage, avgCoverage * getRPKMConstant(chr)};
		return rtrn;
	}

	// to compute score take current score and subtract base pair before and add base pair after ensuring that none of the added is already overlapping the segment and none of the subtracted is in segment
	//iterate from 0 to numMarkers-w
	//reverse iterate from numMarkers to w
	public double [] scanRegion(int fixedWidth,  LightweightGenomicAnnotation region) throws IOException{
		int counter=0;
		double [] rtrn = new double[region.length() - fixedWidth];
		IntervalTree<Alignment> chunkAlignmentTree=getIntervalTree(region);
		String chr = region.getChromosome();
		double score=0;
		Annotation previous=null;
		int until = region.getEnd() - fixedWidth;
		logger.info("Region: "+ region.toUCSC() + ", size " + region.length()+ " array size " + rtrn.length +" until " + until);
		for(int i=region.getStart(); i<until; i++){
	
			long startTime=System.currentTimeMillis();
	
			//if score is 0 jump to end of fixed width
			int start=i;
			int end=start+fixedWidth;
			Annotation current=new Alignments(region.getChromosome(), start, end);
			//System.err.println("Current: " + current.toUCSC());
			double sum=0;
	
			if(score==0){
				sum=data.getCountsPerAlignment(current, chunkAlignmentTree,0); 
				//System.err.println("Computed sum from Annotation: " + sum);
				if(sum ==0) {
					i=start+fixedWidth;
					//System.err.println("\tScore was 0 and so was sum, advancing from " + start + " to " +i);
				} 
				
			}else{
				/*
				if have part of the interval scored and cached then to get next score just need to compute piece not contained
				get score for base not covered and add
				get score for base covered in previous that no longer contained and subtract
				 */		
				Annotation startPosition = new Alignments(chr, previous.getStart(), current.getStart());
				double subtractVal=data.getCountsOfUniqueOverlappers(startPosition, current, chunkAlignmentTree, 0);
				Annotation endPosition = new Alignments(chr, previous.getEnd(), current.getEnd());
				double addVal=data.getCountsOfUniqueOverlappers(endPosition, previous, chunkAlignmentTree, 0);
				sum=(score-subtractVal)+addVal;
				//System.err.println("\tp: " + previous.toUCSC() + " c: " +current.toUCSC()+ " " +data.getCountsPerAlignment(current, chunkAlignmentTree, 0)+" "+subtractVal+" "+addVal+" "+data.getCountsPerAlignment(previous, chunkAlignmentTree, 0)+" score: "+score + " sum: " + sum);
			}
			if(i < until) {
				rtrn[i - region.getStart()] =  sum;
				score=sum;
				previous=current;
	
				counter++;
			}
		}
	
		return rtrn;
	}
	
	public Gene getPeak(Gene gene, Gene startCodon, IntervalTree<Alignment> tree, int EF){
		return this.data.getPeak(gene, startCodon, tree, EF);
	}

	public double[] scanPRate(LightweightGenomicAnnotation annotation, int extensionFactor) throws IOException {
		String chr=annotation.getChromosome();
		double sum=data.getCountsPerAlignment(annotation, extensionFactor); //consider getting number of unique reads as well
		int length=annotation.length();
		double enrich=(sum/(double)length)/getLambda(chr);
		double avgCoverage = sum/(double) length; 
		//System.err.println(chr+" "+sum+" "+count+" "+getLambda(chr));
		//double[] rtrn={calculatePVal(new Double(sum).intValue(), getLambda(chr), length, getNumberMarkers(chr)), enrich, sum, avgCoverage,avgCoverage * rpkmConstant.get(chr)};
		double[] rtrn={ScanStatistics.calculatePVal(new Double(sum).intValue(), getLambda(chr), length, getNumberMarkers(chr)), enrich, sum, avgCoverage,avgCoverage * getRPKMConstant(chr)};
		return rtrn;
	}
	
	public double[] scanPRate(Gene gene, IntervalTree<Alignment> tree,int extensionFactor) throws IOException {
		String chr=gene.getChr();
		return scanPRate(gene, tree, getLambda(chr), extensionFactor);
	}
	
	public double[] scanPRate(Annotation align, IntervalTree<Alignment> tree, double localLambda, int extensionFactor) throws IOException{

		return scanPRate(align, data.getCountsPerAlignment(align, tree, extensionFactor), localLambda, extensionFactor);
	}
	
	public double[] scanPRate(Annotation align, double readCount, int extenstionFactor) throws IOException{
		return scanPRate(align, readCount, getLambda(align.getChr()), extenstionFactor);
	}
	
	public double[] scanPRate(Annotation align, double readCount, double localLambda, int extenstionFactor) throws IOException{
		String chr=align.getChr();
		double sum=readCount;
		int length=align.length();
		double avgCoverage = sum/(double) length;
		double enrich=avgCoverage/localLambda;
		//double rpkm = avgCoverage * rpkmConstant.get(chr) ;
		double rpkm = avgCoverage * getRPKMConstant(chr);
		double nominalP=ScanStatistics.poisson((int) sum,(localLambda*length));
		double[] rtrn={ScanStatistics.calculatePVal(new Double(sum).intValue(), localLambda, length, getNumberMarkers(chr)), enrich, sum, avgCoverage, rpkm, localLambda, length, nominalP};
		return rtrn;
	}
	
	
	public double[] scanPRate(Gene gene, IntervalTree<Alignment> tree,double localLambda, int extensionFactor) throws IOException {
		String chr=gene.getChr();

		double sum=data.getCountsPerAlignment(gene, tree, extensionFactor);
		
		int geneLength=gene.getTranscriptLength();
		double avgCoverage = sum/(double) geneLength;
		double fullyContained=data.getCountsPerAlignmentFullyContained(gene, tree, extensionFactor);
		
		double enrich=avgCoverage/localLambda;
		//double rpkm = avgCoverage * rpkmConstant.get(chr) ;
		double rpkm = avgCoverage * getRPKMConstant(chr) ;
		//Added by Moran 5/11
		double nominalP=ScanStatistics.poisson((int) sum,(localLambda*geneLength));
		double[] rtrn={ScanStatistics.calculatePVal(new Double(sum).intValue(), localLambda, geneLength, getNumberMarkers(chr)), enrich, sum, avgCoverage, rpkm, localLambda, geneLength,nominalP, fullyContained};
		return rtrn;
	}

	public double[] scanPRate(Gene gene, IntervalTree<Alignment> tree,double localLambda, int extensionFactor, double numMarkers) throws IOException {
		String chr=gene.getChr();

		double sum=data.getCountsPerAlignment(gene, tree, extensionFactor);
		
		int geneLength=gene.getTranscriptLength();
		double avgCoverage = sum/(double) geneLength;
		double fullyContained=data.getCountsPerAlignmentFullyContained(gene, tree, extensionFactor);
		
		double enrich=avgCoverage/localLambda;
		//double rpkm = avgCoverage * rpkmConstant.get(chr) ;
		double rpkm = avgCoverage * getRPKMConstant(chr) ;
		//Added by Moran 5/11
		double nominalP=ScanStatistics.poisson((int) sum,(localLambda*geneLength));
		double[] rtrn={ScanStatistics.calculatePVal(new Double(sum).intValue(), localLambda, geneLength, numMarkers), enrich, sum, avgCoverage, rpkm, localLambda, geneLength,nominalP, fullyContained};
		return rtrn;
	}
	
	public double count(Annotation align, int extensionFactor) throws IOException {
		IntervalTree<Alignment> tree=data.getIntervalTreeCached(align.getChr(), align.getStart(), align.getEnd());
		double sum=data.getCountsPerAlignment(align, tree, extensionFactor);//consider getting number of unique reads as well
		return sum;
	}
	
	public int countWithinExon(Annotation align, int extensionFactor) throws IOException {
		//IntervalTree<Alignment> tree=data.getIntervalTreeCached(align.getChr(), align.getStart(), align.getEnd());
		CloseableIterator<Alignment> iter=data.getReadIterator(align);
		int sum=data.getCountsWithinExons(align, iter, extensionFactor); //consider getting number of unique reads as well
		iter.close();
		return sum;
	}

	
	public boolean isPaired() {return isPaired;}

	public int chromosomeLength(String chr) {
		return data.getChromosomeLength(chr);
	}

	public boolean passes(Gene window,IntervalTree<Alignment> chunkAlignmentTree, int i, int critVal) {
		return data.passes(window, chunkAlignmentTree, i, critVal);
	}

	public void normalizeReadCountsToBaseModel(AlignmentDataModelStats baseline) throws IOException {
		setCountNormalizer(new ReadCountNormalizerToBaselineAlignment(this, baseline));
	}
	
	public static class ReadCountNormalizerToBaselineAlignment implements ReadCountNormalizer {
		private double factor;
		
		private ReadCountNormalizerToBaselineAlignment(AlignmentDataModelStats toNormalize,	AlignmentDataModelStats baseline) throws IOException {
			double totalReads = toNormalize.getTotalReads();
			factor = totalReads > 0 ? baseline.getTotalReads()/totalReads : 1;
			logger.info("factor: " + factor);
		}

		@Override
		public double normalize(double count) {
			return factor * count;
		}
		
	}

	public double getMaskedRegionsSize(String chr) {
		return maskedRegions.containsKey(chr) ? maskedRegions.get(chr) : 0;
	}

	public void setStoredStats(StoredAlignmentStats sas) {
		this.globalRPKMConstant = sas.getGlobalRPKMConstant();
		this.isPaired = sas.isAlignmentPaired();
	
		this.numberOfReads = new TreeMap<String, Double>(sas.getNumberOfReads());
		this.maskedRegions = new TreeMap<String, Integer>(sas.getMaskedRegions());
		this.numMarkers    = new TreeMap<String, Double>(sas.getNumMarkers());
		this.rpkmConstant  = new TreeMap<String, Double>(sas.getLocalRpkmConstant());
		this.lambda        = new TreeMap<String, Double>(sas.getLambda());
		
	}

	public void setExtensionFactor(int extensionFactor) {
		data.setExtensionFactor(extensionFactor);
		
	}
}
