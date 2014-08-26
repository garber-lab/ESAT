package umms.esat;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
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
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import umms.core.alignment.Alignment;
import umms.core.alignment.AbstractPairedEndAlignment.TranscriptionRead;
import umms.core.annotation.Annotation;
import umms.core.annotation.BasicAnnotation;
import umms.core.annotation.Gene;
import umms.core.annotation.Annotation.Strand;
import umms.core.coordinatesystem.TranscriptInGenomicSpace;
import umms.core.exception.RuntimeIOException;
import umms.core.model.AlignmentModel;
import umms.core.model.score.BasicScanStatisticScore;
import umms.core.model.score.EndInWindowCountScore;
import umms.core.model.score.ScanStatisticScore;
import umms.core.model.score.StartInWindowCountScore;
import umms.core.model.score.WindowScore.AbstractWindowScore;
import umms.core.readFilters.MappingQualityFilter;
import umms.core.readFilters.PCRDuplicateFilter;
import umms.core.readFilters.UniqueMappedReadsFilter;
//import net.sf.picard.util.IntervalTree;      /* should already be in picard library import, also in IGV library */






import org.apache.commons.collections15.Predicate;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;

//import org.broad.igv.feature.BEDFileParser;
import broad.pda.annotation.BEDFileParser;
import broad.core.datastructures.*;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;

/**
 * @author skadri
 *
 */
public class ESAT {

	static final String usage = "Usage: ESAT -task <task name> "+
			"\n\tTASK 1: score3P: Computes 3' expression of a given annotation set" + 
			"\n\tTASK 2: score5P: Computes 5' expression of a given annotation set " + 
			"\n**************************************************************"+
			"\n\t\tMANDATORY arguments"+
			"\n**************************************************************"+
			"\n\n\t\t-alignment <Alignment (mapped to genome) for which expression must be calculated. Index file MUST be provided.> "+
			"\n\t\t\t OR"+
			"\n\t\t-alignments <Three column text file listing, one line per sample, for which expression must be calculated.Index files MUST be provided." +
			"\n\t\t\t Sample_Name\t Alignment_file\t Grouping_for_sample (Replicates of a condition have the same grouping)> "+
			"\n\n\t\t-annotations <Annotation file for which to calculate expression. [BED by default]> "+

			"\n\n**************************************************************"+
			"\n\t\tOPTIONAL arguments"+
			"\n**************************************************************"+
			"\n\t\t-singleEnd Only if data is single end"+
			"\\nn\t\t-window <Window used to score gene. Default is 500bp> "+
			"\n\t\t-maxIntoGene <Maximum distance from annotated 3' end to find the best scoring window within the gene> "+
			"\n\t\t-maxExtension <It will allow the program to search for the best scoring window maxExtension bases from the 3' end of the gene> "+
			"\n\t\t-minMappingQuality <Only use reads mapping quality greater than the specified>" +
			"\n\t\t-dontWeighReadsFlag <DEPRECATED - NOT AVAILABLE. If provided the read counts will NOT be penalized for multiple mapping >"+
			"\n\t\t-removePCRDuplicatesFlag <If provided, PCR duplicates are removed DEFAULT:false>"+
			"\n\t\t-strand <VALUES: first, second, unstranded. Specifies the mate that is in the direction of transcription DEFAULT: Unstranded> "+
			"\n\t\t-oppositeStrand <If provided indicates that the read in the opposite direction of transcription should be used for quantification. Default is same strand.>"+
			"\n\t\t-reconstructions <If provided a Bed/GTF file with the reconstructed genes is supplied, 3' & 5' ends will be adjusted based on this data.maxExtension & maxIntoGene will be set to 0.>"+
			"\n\t\t-dontCollapseIsoforms <By default, the bed file will be collapse with a 40% overlap between isoforms collapsed to single genes. The output will be on genes and on inidividual isoforms. " +
			"If this flag is provided, this will NOT be done.>"+
			"\n\t\t-out <Output file [Defaults to stdout]> "+
//			"\n\t\t-maskedRegions <Path to file of masked regions in tab delimited format chr start end> "+
			"\n\t\t-filterMultimappers <If provided, multimapped reads will NOT be counted. By default: Multimappers are penalized.> "+
			"\n\t\t-geneMap <If this geneID to geneName tab-delimited map is specified, the annotations are collapsed into overlapping genes and map is used to name the genes. Annotations are then done at the GENE level."+

			"\n\n**************************************************************"+
			"\n\t\tArguments specific to multiple file run"+
			"\n**************************************************************"+
			"\n\n\t\t-normalizedOutput <If present, an additional matrix file with normalized count values; FALSE by default>"+
			"\n\t\t-scoreFullGene <If present, will return the DGE score as the score over the entire gene, after collapsing the annotations. " +
						"If -dontCollapseIsoforms is specified, this will not be calculated. DEFAULT: FALSE>"+
			
			"\n\n\tTASK 3: normalizeMatrix: Takes an input DGE matrix, and outputs a normalized matrix." + 
			"\n\t\t-in <Input matrix> "+
			"\n\t\t-out <Normalized matrix> "+
			"\n\t\t-debugMode <Will output additional results>"+

			"\n";
	/*
	 * Output column names
	 */
	static final String ANNOTATED_END_EXPR_COL = "AnnotatedEndExpression";
	static final String UPSTREAM_EXPR_COL = "BestUpstreamExpression";
	static final String DOWNSTREAM_EXPR_COL = "BestDownstreamExpression";
	static final String BEST_EXPR_COL = "Expression";
	static final String BEST_PVAL_COL = "p-value";
	static final String BEST_DIST_TO_END_COL = "expressionWindowDistToEnd";
	public static String whitespaceDelimiter = "\\s++"; //$NON-NLS-1$
	
	static Logger logger = Logger.getLogger(ESAT.class.getName());/*
	 * Indices of important scores
	 */
	static final int COUNT_SCORE = 0;
	static final int PVAL_SCORE = 6;
	static final int RPKM_SCORE = 1;
	static final int ANNOTATION_LENGTH = 5;
	static final int LOCAL_LAMBDA = 8;
	static final int STEP = 20;
	static final double MIN_OVERLAP = 0.25;
	static final double PVAL_THRESHOLD = 0.05;

	/*
	 * We use the counts as the score to compare window scores
	 */
	static final int USE_SCORE = COUNT_SCORE;
	
	private static String annotationFile;
	private static int maxIntoGene;
	private static int maxExtension;
	private static int minimumMappingQuality;
	private static int window;
	private static boolean weighReadsFlag;
	private static boolean removePCRDuplicatesFlag;
	private static boolean fullGeneScoreFlag;
	private static boolean debugMode;
	//@Deprecated
	//private static boolean isStranded;
	private static TranscriptionRead strand;
	private static boolean collapseIsoforms;
	private static boolean oppositeStrand;
	private static boolean pairedFlag;
	private static boolean filterMultimappers;
	private static HashMap<Gene, String> duplicateNameMap;
	private static Map<Gene, Set<Gene>> collapsedGeneMap;
	private static Map<String, IntervalTree<Gene>> collapsedGenes;
	//each row of the matrix is one Gene name
	static List<String> rows = new ArrayList<String>();
	static String maskedRegionFile;

	public ESAT(String[] args) throws IOException, ParseException {

		/*
		 * Gives a log4j error. Check later.
		 */
		Globals.setHeadless(true);
		/*
		 * @param for ArgumentMap - size, usage, default task
		 * argMap maps the command line arguments to the respective parameters
		 */
		ArgumentMap argMap = CLUtil.getParameters(args,usage,"score3P");

		if ("normalizeMatrix".equalsIgnoreCase(argMap.getTask())) {			
			writeNormalizedMatrix(argMap.getOutput(),new MatrixWithHeaders(argMap.getInput()));
		}
		else{
			/*
			 * Parameters if not provided are set to defaults:
			 * maxIntoGene : default 0
			 * maxExtension : default 0
			 * window: default 500
			 * minMappingQuality : default -1
			 */
			maxIntoGene = argMap.isPresent("maxIntoGene")? argMap.getInteger("maxIntoGene") : 0;
			maxExtension = argMap.isPresent("maxExtension")? argMap.getInteger("maxExtension") : 0;
			minimumMappingQuality = argMap.isPresent("minMappingQuality")? argMap.getInteger("minMappingQuality") : -1;
			window = argMap.isPresent("window")? argMap.getInteger("window"): 500;
			collapseIsoforms = !argMap.isPresent("dontCollapseIsoforms");
			oppositeStrand = argMap.isPresent("oppositeStrand");
			//TODO: Check if this works
			maskedRegionFile = argMap.containsKey("maskedRegions") ? argMap.get("maskedRegions") : null;
			debugMode = argMap.isPresent("debugMode");
			/*
			 * FLAG for WEIGHING READS BY NH FLAG
			 * TRUE by default
			 * Convert string to boolean
			 */
			//@Deprecated
			weighReadsFlag = !argMap.isPresent("dontWeighReadsFlag");
			
			/*
			 * FLAG for REMOVING PCR DUPLICATES 
			 * FALSE by default
			 * Convert string to boolean
	 		 */
			removePCRDuplicatesFlag = argMap.isPresent("removePCRDuplicatesFlag");
			
			filterMultimappers = argMap.isPresent("filterMultimappers");
			/*
			 * FLAG for TO RETURN SCORE OVER FULL GENE 
			 * FALSE by default
			 * Convert string to boolean
			 */
			fullGeneScoreFlag = argMap.isPresent("scoreFullGene");
	
			//@Deprecated
			//isStranded = argMap.isPresent("stranded");
			strand = TranscriptionRead.UNSTRANDED;
			if(argMap.get("strand").equalsIgnoreCase("first")){
				if(oppositeStrand)
					strand = TranscriptionRead.SECOND_OF_PAIR;
				strand = TranscriptionRead.FIRST_OF_PAIR;
			}
			else if(argMap.get("strand").equalsIgnoreCase("second")){
				if(oppositeStrand)
					strand = TranscriptionRead.FIRST_OF_PAIR;
				strand = TranscriptionRead.SECOND_OF_PAIR;
			}
			
			pairedFlag = !argMap.isPresent("singleEnd");
							
			//IF user wants the normalized matrix,
			/*
			 * FLAG to return a normalized matrix. FALSE by default. Convert string to boolean
			 */
			boolean normalizedOutput = argMap.isPresent("normalizedOutput");
	
			/*
			 * Read the annotation file
			 */
			annotationFile = argMap.getMandatory("annotations");		
			/*
			 * Check the format of the annotation file and call the GTF or BED parser accordingly
			 */
			if(annotationFile.endsWith(".gtf") || annotationFile.endsWith(".GTF")){
				logger.error("Please supply an annotation file in the BED format",new RuntimeIOException());
			}
			Map<String,Collection<Gene>> annotations =  BEDFileParser.loadDataByChr(new File(annotationFile));	
			/**
			 * If user wants to reconstruct gene ends.
			 */
			if(argMap.isPresent("reconstructions")) {
				
				String reconstructionFile = argMap.getMandatory("reconstructions");
				BEDFileParser reconstructions =  new BEDFileParser(reconstructionFile);
				logger.info("A reconstructed bed file provided, going to adjust ends of genes and ignore 3' extensions");
	
				BEDFileParser annotationParser = new BEDFileParser(annotationFile);
				annotationParser = reconstructGeneEnds(reconstructions, annotationParser, annotationFile);
				annotations = annotationParser.toMap();
				maxIntoGene = 0;
				maxExtension = 0;
			}
			
			logger.info("Using window: " + window + " max 5' extension: " + maxExtension + " maximum  premature start? "+ maxIntoGene + " minimum mapping quality "+ minimumMappingQuality+" opposite strand "+oppositeStrand);

			//HashMap of gene name to the number of duplicates
			HashMap<String, Integer> duplicateMap = new HashMap<String, Integer>();
			//HashMap of gene to rowName
			duplicateNameMap = new HashMap<Gene, String>();
			int duplicates=0;
			for(String chr:annotations.keySet()){
				for(Gene gene : annotations.get(chr)) {
					if(!rows.contains(gene.getName())) {
						String name = gene.getName();
						if(duplicateNameMap.containsKey(gene)){
							logger.info("Entry for "+name+" already exists");
						}
						else{
							rows.add(name);
							duplicateMap.put(name, 1);
							duplicateNameMap.put(gene, name);
						}
					} 
					// If the gene name has another annotation
					else {
						
						if(duplicateNameMap.containsKey(gene)){
							logger.info("Entry for "+gene.getName()+" already exists in "+duplicateNameMap.get(gene));
						}
						else{
							if(!duplicateMap.containsKey(gene.getName()))
								duplicateMap.put(gene.getName(), 1);
							//Row name is now the geneName appended with the duplicate number
							duplicateMap.put(gene.getName(), (duplicateMap.get(gene.getName())+1));
							String name = (gene.getName()+"_"+duplicateMap.get(gene.getName()));
							rows.add(name);
							duplicateNameMap.put(gene, name);
							//logger.warn("Duplicated annotation : "+ gene.toBED());
							duplicates++;
						}
					}
				}
			}
			logger.info("Found " + duplicates + " duplicates, ignoring them, going to process " + rows.size() + " annotations");
			
//			/* test - AD */
//			Collection<Integer> dupValues = duplicateMap.values();
//			Set<String> dupKeys = duplicateMap.keySet();
//			Integer dupCount = 0;
//			for (String dupKey : dupKeys) {
//				Integer dupVal = duplicateMap.get(dupKey);
//				if (dupVal > 1) {
//					logger.info(dupKey + " has " + dupVal + " duplicates");
//					
//				}
//			}
//		
//			logger.info("duplicateMap has " + dupKeys.size() + " entries");
			
			if(collapseIsoforms)
				setCollapsedGeneMap(annotations);
			/*
			 * If "alignment" is present, it is single file
			 */
			if(argMap.isPresent("alignment")){
				
				scoreSingleSample(annotations,argMap.getMandatory("alignment"),argMap.getOutput(),"score3P".equalsIgnoreCase(argMap.getTask()));			
			}
			/*
			 * if "alignments" is present, it is multiple run
			 */
			 else if(argMap.isPresent("alignments")){
				 /*
				  * Read the names of the alignment file into an array
				  */
				 BufferedReader br = new BufferedReader(new FileReader(argMap.getMandatory("alignments")));
				 Map<String,String> alignmentFiles = new HashMap<String,String>();
				 Map<String,Collection<String>> conditionMaps = new HashMap<String,Collection<String>>();
				 String s;
				 while((s = br.readLine())!= null){
					 //Map of sample name to alignment File
					 alignmentFiles.put(s.split(whitespaceDelimiter)[0],s.split(whitespaceDelimiter)[1]);
					 String cond = s.split(whitespaceDelimiter)[2];
					 if(!conditionMaps.containsKey(cond) ) {
						 conditionMaps.put(cond, new ArrayList<String>());
					 }
					 Collection<String> groupSamples = conditionMaps.get(cond);
					 groupSamples.add(s.split(whitespaceDelimiter)[0]);	
				 }					

				 scoreMultipleSamples(annotations,alignmentFiles,argMap.getOutput(),normalizedOutput,conditionMaps,"score3P".equalsIgnoreCase(argMap.getTask()));
			 }
			 else{
				 throw new IllegalArgumentException("Argument alignment or alignments must be provided\n"+usage);
			 }
		}
	}

	/**
	 * This function collapses the genes
	 * @param annotations
	 * @return
	 */
	private void setCollapsedGeneMap(Map<String,Collection<Gene>> annotations){
		
		collapsedGeneMap = new HashMap<Gene,Set<Gene>>();
		collapsedGenes = new HashMap<String, IntervalTree<Gene>>();
		String geneString = "gene_";
		
		for(String chr:annotations.keySet()){
			IntervalTree<Gene> oldTree = new IntervalTree<Gene>();
			for(Gene g:annotations.get(chr)){
				oldTree.put(g.getStart(), g.getEnd(), g);
			}
			
			IntervalTree<Gene> newTree = new IntervalTree<Gene>();
			collapsedGenes.put(chr, newTree);
			
			//For all genes in old tree
			for(Gene currentElement:oldTree.toCollection()) {
				//Find overlapping genes in new tree
				Iterator<Gene> overlapperIt = newTree.overlappingValueIterator(currentElement.getStart(), currentElement.getEnd());
				Gene overlapper = null;
				int overlapperNum = 0;
				//While there are more overlappers
				while(overlapperIt.hasNext()) {
					Gene overlapperCandidate= overlapperIt.next();
					//If overlapper and gene are compatible
					if(BEDFileParser.isOverlapCompatible(currentElement,overlapperCandidate, MIN_OVERLAP) ) {
						overlapper = overlapperCandidate;
					}				
					if(overlapper != null) {
						overlapperNum++;
						//Remove overlapper from new tree
						newTree.remove(overlapper.getStart(),overlapper.getEnd());
						if(!collapsedGeneMap.containsKey(overlapper)){
							logger.info("Not contains "+overlapper.getName());
						}
						Set<Gene> constituents = collapsedGeneMap.remove(overlapper);
						//Merge gene and overlapper
						Gene mergedElement=new Gene(overlapper.takeUnion(currentElement));
						//Set new name
						mergedElement.setName(overlapper.getName()+"_"+duplicateNameMap.get(currentElement));
						constituents.add(currentElement);
						constituents.add(currentElement);						
						collapsedGeneMap.put(mergedElement, constituents);
						
						//Put merged in new tree
						newTree.put(mergedElement.getStart(), mergedElement.getEnd(), mergedElement);
					} 
					overlapper = null;
				}
				//If no overlappers, add the gene itseld
				if(overlapperNum == 0) {
					//Bug fix : Dec 27 2010 - the function was not merging a refseq gene with isoforms that had several isoforms
					Gene mergedElement=new Gene(currentElement);
					//Set new name
					mergedElement.setName(geneString+duplicateNameMap.get(currentElement));
					Set<Gene> constituents = new TreeSet<Gene>();
					constituents.add(currentElement);
					collapsedGeneMap.put(mergedElement, constituents);
					newTree.put(mergedElement.getStart(), currentElement.getEnd(),  mergedElement);
				}
			}
		}
	}
	/**
	 * This function takes an input count matrix and returns a normalized matrix, where a scaling factor is calculated for each sample,
	 * and all genes for that sample are scaled by this factor.
	 * @param fileName
	 * @param resultMatrix
	 * @throws IOException
	 */
	public static void writeNormalizedMatrix(String fileName, MatrixWithHeaders resultMatrix) throws IOException{
		
		Map<String,Double> factors = calculateNormalizationFactors(resultMatrix);
		MatrixWithHeaders normalizedMatrix = resultMatrix.multiplyColumnsWithConstants(factors);
		String normalizedFileName = fileName+".normalized";
		String normalizedFactorFile = fileName+".coverage";
		BufferedWriter outBw = new BufferedWriter(new FileWriter(normalizedFileName));
		normalizedMatrix.write(outBw);
		outBw.close();
		
		BufferedWriter fBw = new BufferedWriter(new FileWriter(normalizedFactorFile));
		for(String ss:factors.keySet()){
			fBw.write(ss+"\t"+((double)1.0/factors.get(ss))+"\n");
		}
		fBw.close();
	}
	
	/**
	 * Calculates the quantification profile for a single 3p sample
	 * @param annotations
	 * @param alignmentFile
	 * @param outputFile
	 * @throws IOException
	 * @throws ParseException
	 */
	private static void scoreSingleSample(Map<String,Collection<Gene>> annotations, String alignmentFile,String outputFile,boolean is3p) throws IOException, ParseException {
		
		/*
		 * Initialize the lists for output of expression score
		 * These lists names will also be the columns of the matrix
		 */
		List<String> cols = new ArrayList<String>();
		cols.add(BEST_EXPR_COL);
		cols.add(BEST_PVAL_COL);
		cols.add(ANNOTATED_END_EXPR_COL);
		if(maxIntoGene>0 || maxExtension>0){

			if(maxIntoGene>0)
				cols.add(UPSTREAM_EXPR_COL);
			if(maxExtension>0)
				cols.add(DOWNSTREAM_EXPR_COL);

			cols.add(BEST_DIST_TO_END_COL);
		}

		/* 
		 * Initialize the matrix with row and column names
		 */
		MatrixWithHeaders resultMatrix = new MatrixWithHeaders(rows,cols);
		
		//Map of each gene to its peak
		Map<String,Annotation> geneToWindowMap = new HashMap<String,Annotation>();
								
		logger.info("Scoring using all reads:");
		/*
		 * Initialize the data model using the alignment file
		 * @param: <alignment flieName> <load_chromosome_stats> <minMappingQuality> <remove_PCR_duplicates> <weigh reads by NH flag>
		 */
		AlignmentModel model = new AlignmentModel(alignmentFile, null, new ArrayList<Predicate<Alignment>>(), pairedFlag,strand,true,maskedRegionFile);
		AlignmentModel libDataModel = new AlignmentModel(alignmentFile,new TranscriptInGenomicSpace(model.getRefSequenceLengths()), new ArrayList<Predicate<Alignment>>(), pairedFlag,strand,true,maskedRegionFile);
		libDataModel.addFilter(new MappingQualityFilter(minimumMappingQuality,minimumMappingQuality));
		if(removePCRDuplicatesFlag)
			libDataModel.addFilter(new PCRDuplicateFilter());
		if(filterMultimappers)
			libDataModel.addFilter(new UniqueMappedReadsFilter());
		//To report all peaks - not just the best one
		Map<Gene,Set<Gene>> geneToPeaksMap = new HashMap<Gene,Set<Gene>>();
		
		// For each chromosome
		for(String chr:annotations.keySet()) {
			/*
			 * If the alignment data has data from that chromosome
			 */
			if(libDataModel.containsReference(chr)){
				IntervalTree<Gene> tree = new IntervalTree<Gene>();
				/* for each chromosome, create an IntervalTree and add all genes in the reference 
//				 * file in this chromosome to it. - AD */ 	
				Collection<Gene> annotColl = annotations.get(chr);
				for(Gene g : annotations.get(chr)){
					tree.put(g.getStart(), g.getEnd(), g);
				}
				logger.info("Processing " + chr);
				//Parse over all annotations on chromosome
				for(Gene annotation:annotations.get(chr)){
					
					Set<Gene> thisPeaks = new HashSet<Gene>();
					
					Gene annotationEnd = new Gene(annotation);					
					int annotationLength = annotation.getSize();  /* sum of all exons - AD */
					/*
					 * If the length of the annotated transcript > window size being analyzed,
					 * get sub annotation for window length
					 */
					if(annotationLength>window){
						annotationEnd = is3p? getSubAnnotationFromEnd(annotation,window,0) : getSubAnnotationFromStart(annotation,window,0);

						if(annotationEnd == null){
							logger.warn("Annotation end for " + (annotationLength - window) + "-" + annotationLength + " --> " + annotation.toBED() + " was null.");
						}
					}
//					System.out.println(annotationEnd.getSequence());
					/*
					 * Calculate
					 * [0] = count
					 * [1] = RPKM
					 * [2] = region total
					 * [3] = total
					 * [4] = scan p value
					 * [5] = global lambda
					 * [6] = local lambda
					 * [7] = region length  
					 * for the annotated region 	
					 */
					if(oppositeStrand)
						annotationEnd =reverseStrand(annotationEnd);
					double [] bestScores = new ScanStatisticScore(libDataModel, annotationEnd,false).getScores();					 
					/*
					 * Use COUNTS as the score
					 * The best score is the window score for now
					 */
					double annotatedEndExpr = bestScores[USE_SCORE];
					//For all peaks
					if(bestScores[PVAL_SCORE]<PVAL_THRESHOLD){
						annotationEnd.setScore(bestScores[USE_SCORE]); //set RPKM as score
						annotationEnd.setExtraFields(bestScores);
						annotationEnd.setName(duplicateNameMap.get(annotation));
						thisPeaks.add(annotationEnd);
					}

					int bestScoreDistanceFromEnd = 0;
					geneToWindowMap.put(duplicateNameMap.get(annotation), annotationEnd);

					double upstreamExpr = 0;
					double downstreamExpr = 0;

					int intoGene = STEP;
					/*
					 * If upstream extension is allowed, use sliding windows with overlaps of STEP
					 * while retreat region length is smaller than (annotation - window), i.e. it lies within the annotation
					 * and it is at distance less than the max region allowed upstream of end of gene
					 */
					while((intoGene<(annotationLength - window)) && (intoGene<maxIntoGene)){

						/*
						 * get annotation for region of length window, "retreat" length from end of transcript
						 */
						annotationEnd = is3p? getSubAnnotationFromEnd(annotation, window, intoGene): getSubAnnotationFromStart(annotation, window, intoGene);
						
						/*
						 * calculate score for the window and find max score of windows
						 * Uses scanPRate(Gene, tree) which adds extension factor in AlignmentDataModelStats.
						 * Difference between scoregene and scanPRate from above is that this does not include exons.
						 */
						if(annotationEnd !=null){
							if(oppositeStrand)
								annotationEnd =reverseStrand(annotationEnd);
							double[] tmpScores = new ScanStatisticScore(libDataModel, annotationEnd,false).getScores();
							
							//For all peaks
							//If significant
							boolean overlaps = false;
							if(tmpScores[PVAL_SCORE]<PVAL_THRESHOLD){
								//if overlaps any already added peak
								Set<Gene> tempPeaks = new HashSet<Gene>();
								tempPeaks.addAll(thisPeaks);
								for(Gene peak:thisPeaks){
									if(annotationEnd.overlaps(peak)){
										overlaps = true;
										//if better then replace
										if(tmpScores[USE_SCORE]>peak.getBedScore()){
											annotationEnd.setBedScore(tmpScores[USE_SCORE]); //set RPKM as score
											annotationEnd.setExtraFields(bestScores);
											annotationEnd.setName(duplicateNameMap.get(annotation));
											tempPeaks.remove(peak);
											tempPeaks.add(annotationEnd);
										}
									}
								}
								thisPeaks = tempPeaks;
								if(!overlaps){
									annotationEnd.setBedScore(tmpScores[USE_SCORE]); //set RPKM as score
									annotationEnd.setExtraFields(bestScores);
									annotationEnd.setName(duplicateNameMap.get(annotation));
									//does not overlap. Add
									thisPeaks.add(annotationEnd);	
								}
							}
							//For best peak
							if(tmpScores[USE_SCORE]>bestScores[USE_SCORE]){
								for(int i=0;i<bestScores.length;i++){
									bestScores[i] = tmpScores[i];
								}
								bestScoreDistanceFromEnd = is3p? -intoGene : intoGene; 
									
								geneToWindowMap.put(duplicateNameMap.get(annotation), annotationEnd);
							}
							if((is3p && tmpScores[USE_SCORE]>upstreamExpr))
								upstreamExpr = tmpScores[USE_SCORE];
							if(!is3p && tmpScores[USE_SCORE]>downstreamExpr)
								downstreamExpr = tmpScores[USE_SCORE];
						}
						intoGene += STEP;
					}

					int extend = STEP;
					while(extend < maxExtension) {
						Annotation end = null;
						if((is3p && annotation.getOrientation().equals(Strand.NEGATIVE))||(!is3p && annotation.getOrientation().equals(Strand.POSITIVE))){
							if((annotation.getStart() - extend) < 0){
								logger.info(annotation.getStart()+" - "+extend +" is less than 0");
								break;
							}
							end = new BasicAnnotation(annotation.getChr(), annotation.getStart() - extend, annotation.getStart() - (extend-window),annotation.getOrientation());
						}
						else{
							if((annotation.getEnd() + extend)> libDataModel.getRefSequenceLength(chr)){
								logger.info(annotation.getEnd()+" + "+extend +" is more than "+libDataModel.getRefSequenceLength(chr));
								break;
							}
							end = new BasicAnnotation(annotation.getChr(), annotation.getEnd() + (extend - window), annotation.getEnd() + extend,annotation.getOrientation());
						}
						
						/*
						 * Get an interval tree for all/any exons that overlap with the extended region
						 */
						Iterator<Gene> endOverlappersIter = tree.overlappingValueIterator(end.getStart(), end.getEnd());
						/*
						 * While there is an overlap with a gene
						 * and gene is same gene
						 */
						boolean overlapperIsSameGene = true;
						while(endOverlappersIter.hasNext() && overlapperIsSameGene){
							
							Gene overlapper = endOverlappersIter.next();
							//If same orientation
							if(overlapper.getOrientation().equals(annotation.getOrientation())){
								//compare the end coordiantes of the gene
								if(is3p){
									if(!(overlapper.getOrientedEnd() == annotation.getOrientedEnd()))
										if(!BEDFileParser.isOverlapCompatible(annotation,overlapper, MIN_OVERLAP))
											overlapperIsSameGene = false;
								}
								else{
									if(!(overlapper.getOrientedStart() == annotation.getOrientedStart()) && !overlapper.overlaps(annotation))
										if(!BEDFileParser.isOverlapCompatible(annotation,overlapper, MIN_OVERLAP))
											overlapperIsSameGene = false;
								}
							}
						}
						if(!overlapperIsSameGene)
							break;
						annotationEnd = new Gene(end);
						//No overlap so continue with scoring the region
						if(oppositeStrand)
							annotationEnd =reverseStrand(annotationEnd);
						double[] tmpScores = new ScanStatisticScore(libDataModel, annotationEnd,false).getScores();
						//For all peaks
						//If significant
						boolean overlaps = false;
						if(tmpScores[PVAL_SCORE]<PVAL_THRESHOLD){
							//if overlaps any already added peak
							Set<Gene> tempPeaks = new HashSet<Gene>();
							tempPeaks.addAll(thisPeaks);
							for(Gene peak:thisPeaks){
								if(annotationEnd.overlaps(peak)){
									overlaps = true;
									//if better then replace
									if(tmpScores[USE_SCORE]>peak.getBedScore()){
										annotationEnd.setBedScore(tmpScores[USE_SCORE]); //set RPKM as score
										annotationEnd.setExtraFields(bestScores);
										annotationEnd.setName(duplicateNameMap.get(annotation));
										tempPeaks.remove(peak);
										tempPeaks.add(annotationEnd);
									}
								}
							}
							thisPeaks = tempPeaks;
							if(!overlaps){
								annotationEnd.setBedScore(tmpScores[USE_SCORE]); //set RPKM as score
								annotationEnd.setExtraFields(bestScores);
								annotationEnd.setName(duplicateNameMap.get(annotation));
								//does not overlap. Add
								thisPeaks.add(annotationEnd);	
							}
						}
						
						//For best peak
						if(tmpScores[USE_SCORE]>bestScores[USE_SCORE]){
							for(int i=0;i<bestScores.length;i++){
								bestScores[i] = tmpScores[i];
							}
							bestScoreDistanceFromEnd = is3p? extend: -extend;
							geneToWindowMap.put(duplicateNameMap.get(annotation), annotationEnd);
						}
						if(is3p && tmpScores[USE_SCORE]>downstreamExpr)
							downstreamExpr = tmpScores[USE_SCORE];
						if(!is3p && tmpScores[USE_SCORE]>upstreamExpr)
							upstreamExpr = tmpScores[USE_SCORE];
						extend += (STEP);
					}
					
					/*
					 * Write the result for this annotation to the matrix
					 */
					resultMatrix.set(duplicateNameMap.get(annotation), BEST_EXPR_COL, bestScores[USE_SCORE]);
					resultMatrix.set(duplicateNameMap.get(annotation), BEST_PVAL_COL, bestScores[PVAL_SCORE]);
					resultMatrix.set(duplicateNameMap.get(annotation), ANNOTATED_END_EXPR_COL, annotatedEndExpr);
					if(maxIntoGene > 0)
						resultMatrix.set(duplicateNameMap.get(annotation), UPSTREAM_EXPR_COL, upstreamExpr);
					if(maxExtension > 0)
						resultMatrix.set(duplicateNameMap.get(annotation), DOWNSTREAM_EXPR_COL, downstreamExpr);
					if(maxIntoGene > 0 || maxExtension > 0)
						resultMatrix.set(duplicateNameMap.get(annotation), BEST_DIST_TO_END_COL, bestScoreDistanceFromEnd);
					
					annotation.setBedScore(bestScores[USE_SCORE]); //set RPKM as score
					annotation.setExtraFields(bestScores);
					annotation.addExtraField(duplicateNameMap.get(annotation));
					geneToPeaksMap.put(annotation, thisPeaks);
						
				}
			}
			/**
			 * If there is no data for chromosome in annotation file, 
			 * put zeroes in all fields
			 */
			else{
				for(Gene annotation:annotations.get(chr)){
					
					double [] bestScores = new double[5];
					resultMatrix.set(duplicateNameMap.get(annotation), BEST_EXPR_COL, 0.0);
					resultMatrix.set(duplicateNameMap.get(annotation), BEST_PVAL_COL, 1.0);
					resultMatrix.set(duplicateNameMap.get(annotation), ANNOTATED_END_EXPR_COL, 0.0);
					if(maxIntoGene > 0)
						resultMatrix.set(duplicateNameMap.get(annotation), UPSTREAM_EXPR_COL, 0.0);
					if(maxExtension > 0)
						resultMatrix.set(duplicateNameMap.get(annotation), DOWNSTREAM_EXPR_COL, 0.0);
					if(maxIntoGene > 0 || maxExtension > 0) 
						resultMatrix.set(duplicateNameMap.get(annotation), BEST_DIST_TO_END_COL, 0.0);
					
					annotation.setBedScore(bestScores[USE_SCORE]); //set RPKM as score
					annotation.setExtraFields(bestScores);
					annotation.addExtraField(duplicateNameMap.get(annotation));
				}
			}
		} 

		/*
		 * Write to output File
		 */
		writeOutputFiles(outputFile, resultMatrix, annotations);
		if(collapseIsoforms){
			writeCollapsedOutputForSingle(annotations,outputFile,resultMatrix,geneToWindowMap);
		}
		
		if(debugMode){
			BufferedWriter bedBw = new BufferedWriter(new FileWriter(outputFile+".allPeaks.bed"));
			for(Gene gene: geneToPeaksMap.keySet()){
				for(Gene peak:geneToPeaksMap.get(gene)){
					bedBw.write(peak.toBED());
					//bedBw.append("\t"+duplicateNameMap.get(gene));
					bedBw.newLine();
				}			
			}
			bedBw.close();
		}
	}
	
	private static Annotation reverseStrand(Annotation annotationEnd) {
		if(annotationEnd.getOrientation().equals(Strand.POSITIVE)){
			annotationEnd.setOrientation(Strand.NEGATIVE);
		}
		else{
			if(annotationEnd.getOrientation().equals(Strand.NEGATIVE))
				annotationEnd.setOrientation(Strand.POSITIVE);
		}
		return annotationEnd;
	}

	private static Gene reverseStrand(Gene annotationEnd) {
		if(annotationEnd.getOrientation().equals(Strand.POSITIVE)){
			annotationEnd.setOrientation(Strand.NEGATIVE);
		}
		else{
			if(annotationEnd.getOrientation().equals(Strand.NEGATIVE))
				annotationEnd.setOrientation(Strand.POSITIVE);
		}
		return annotationEnd;
	}

	/**
	 * Calculates the quantification profile for a Multiple 3p sample
	 * @param annotations
	 * @param alignmentFile
	 * @param outputFile
	 * @throws IOException
	 * @throws ParseException
	 */
	private void scoreMultipleSamples(Map<String, Collection<Gene>> annotations,Map<String, String> alignmentFiles, 
			String output,boolean normalizedOutput,Map<String, Collection<String>> conditionMaps,boolean is3p) throws IOException {

		List<Gene> windows = new ArrayList<Gene>();		
		//Map of each gene to its peak
		Map<String,Annotation> geneToWindowMap = new HashMap<String,Annotation>();
			
		/*
		 * The columns of the matrix will be each alignment file name
		 */
		List<String> cols = new ArrayList<String>();
		for(String name: conditionMaps.keySet()){
			for(String sample: conditionMaps.get(name)){
				cols.add(sample);
			}
		}
		
		/* 
		 * Initialize the matrix with row and column names
		 */
		MatrixWithHeaders resultMatrix = new MatrixWithHeaders(rows,cols);
		MatrixWithHeaders pvalueMatrix = new MatrixWithHeaders(rows,cols);
//		MatrixWithHeaders fullGeneResultMatrix = null;
//		if(fullGeneScoreFlag){
			/* 
			 * Initialize the matrix to store score of full gene with row and column names
			 */
/*			fullGeneResultMatrix = new MatrixWithHeaders(rows,cols);
		}
*/		//ELSE DONT INITIALIZE
		logger.info("Scoring using all reads:");

		// Initialize the AlignmentModels
		Map<String,AlignmentModel> libDataModels = new HashMap<String,AlignmentModel>();
		for(String ss: alignmentFiles.keySet()){
			AlignmentModel libmodel = new AlignmentModel(alignmentFiles.get(ss), null, new ArrayList<Predicate<Alignment>>(), pairedFlag,strand,true,maskedRegionFile);
			AlignmentModel model = new AlignmentModel(alignmentFiles.get(ss),new TranscriptInGenomicSpace(libmodel.getRefSequenceLengths()), new ArrayList<Predicate<Alignment>>(), pairedFlag,strand,true,maskedRegionFile);
			if(removePCRDuplicatesFlag)
				model.addFilter(new PCRDuplicateFilter());
			if(filterMultimappers)
				model.addFilter(new UniqueMappedReadsFilter());
			libDataModels.put(ss, model);
		}
		
		//To report all peaks - not just the best one
		Map<Gene,Set<Gene>> geneToPeaksMap = new HashMap<Gene,Set<Gene>>();
						
		// For each chromosome
		for(String chr:annotations.keySet()) {
			
			// If the alignment data has data from that chromosome
			boolean dataForChr = false;
			for(AlignmentModel model:libDataModels.values()){
				//If any alignment file has data for the chromosome, 
				//		set flag to true and exit loop
				if(model.containsReference(chr)){
					dataForChr = true;
					break;
				}
			}
			//If the any alignment file has data from that chromosome
			if(dataForChr){
				IntervalTree<Gene> tree = new IntervalTree<Gene>();
				for(Gene g:annotations.get(chr)){
					tree.put(g.getStart(), g.getEnd(), g);
				}
				logger.info("Processing " + chr);
				//Parse over all annotations on chromosome
				for(Gene annotation:annotations.get(chr)){
					Set<Gene> thisPeaks = new HashSet<Gene>();
					
					Gene annotationEnd = new Gene(annotation);					
					int annotationLength = annotation.getSize();
					//If the length of the annotated transcript > window size being analyzed,
					//get sub annotation for window length
					if(annotationLength>window){
						annotationEnd = is3p? getSubAnnotationFromEnd(annotation,window,0) : getSubAnnotationFromStart(annotation,window,0);

						if(annotationEnd == null){
							logger.warn("Annotation end for " + (annotationLength - window) + "-" + annotationLength + " --> " + annotation.toBED() + " was null.");
						}
					}
					
					Map<String,double[]> bestScores = getScoresForWindow(libDataModels,annotationEnd);
					double bestEnrichment = getEnrichmentForWindow(libDataModels.keySet(),bestScores);
					int bestWindowEnd = 0;
					boolean insideGene = true;
					Gene bestWindow = annotationEnd;

					//USING ENRICHMENT AS MEANS FOR COMPARISON

					//For all peaks
					if(isSignificantInAllModels(libDataModels.keySet(),bestScores)){
						annotationEnd.setScore(bestEnrichment); //set RPKM as score
						annotationEnd.setName(duplicateNameMap.get(annotation));
						thisPeaks.add(annotationEnd);
					}
					
					int intoGene = STEP;
					/*							 
					 * If upstream extension is allowed, use sliding windows with overlaps of STEP
					 * while retreat region length is smaller than (annotation - window), i.e. it lies within the annotation
					 * and it is at distance less than the max region allowed upstream of end of gene
					 */
					while((intoGene<(annotationLength - window)) && (intoGene<maxIntoGene)){
						//get annotation for region of length window, "intoGene" length from end of transcript
						annotationEnd = is3p? getSubAnnotationFromEnd(annotation, window, intoGene): getSubAnnotationFromStart(annotation, window, intoGene);
						/*
						 * calculate score for the window and find max score of windows
						 * Uses scanPRate(Gene, tree) which adds extension factor in AlignmentDataModelStats.
						 * Difference between scoregene and scanPRate from above is that this does not include exons.
						 */ 
						if(annotationEnd !=null){
							Map<String,double[]> tmpScores = getScoresForWindow(libDataModels,annotationEnd);
							double tmpEnrichment = getEnrichmentForWindow(libDataModels.keySet(),tmpScores);
							//For all peaks
							//If significant
							boolean overlaps = false;
							if(isSignificantInAllModels(libDataModels.keySet(), tmpScores)){
								//if overlaps any already added peak
								Set<Gene> tempPeaks = new HashSet<Gene>();
								tempPeaks.addAll(thisPeaks);
								for(Gene peak:thisPeaks){
									if(annotationEnd.overlaps(peak)){
										overlaps = true;
										//if better then replace
										
										if(tmpEnrichment>peak.getBedScore()){
											annotationEnd.setBedScore(tmpEnrichment); //set RPKM as score
											annotationEnd.setName(duplicateNameMap.get(annotation));
											tempPeaks.remove(peak);
											tempPeaks.add(annotationEnd);
										}
									}
								}
								thisPeaks = tempPeaks;
								if(!overlaps){
									annotationEnd.setBedScore(tmpEnrichment); //set RPKM as score
									annotationEnd.setName(duplicateNameMap.get(annotation));
									//does not overlap. Add
									thisPeaks.add(annotationEnd);	
								}
							}
							//For best peak
							if(tmpEnrichment>bestEnrichment){
								bestEnrichment = tmpEnrichment;
								bestWindow = annotationEnd;
								bestWindowEnd = is3p? -intoGene : intoGene;
								insideGene = true;
								bestScores = tmpScores;
								geneToWindowMap.put(duplicateNameMap.get(annotation), annotationEnd);
							}
						}
						intoGene += STEP;
					}
					
					int extend = STEP;
					while(extend < maxExtension) {
						Annotation end = null;
						long chrEnd = 0;
						for(String s:libDataModels.keySet())
							chrEnd = libDataModels.get(s).getRefSequenceLength(chr);
						
						if((is3p && annotation.getOrientation().equals(Strand.NEGATIVE))||(!is3p && annotation.getOrientation().equals(Strand.POSITIVE))){
							if((annotation.getStart() - extend)<0){
								logger.info(annotation.getStart()+" - "+extend +" is less than 0");
								break;
							}
							end = new BasicAnnotation(annotation.getChr(), annotation.getStart() - extend, annotation.getStart() - (extend-window),annotation.getOrientation());
						}else{
							if((annotation.getEnd() + extend)> chrEnd){
								logger.info(annotation.getEnd()+" + "+extend +" is more than "+chrEnd);
								break;
							}
							end = new BasicAnnotation(annotation.getChr(), annotation.getEnd() + (extend - window), annotation.getEnd() + extend,annotation.getOrientation());
						}/*
						 * Get an interval tree for all/any exons that overlap with the extended region
						 */
						Iterator<Gene> endOverlappersIter = tree.overlappingValueIterator(end.getStart(), end.getEnd());
						/*
						 * While there is an overlap with a gene
						 * and gene is same gene
						 */
						boolean overlapperIsSameGene = true;
						while(endOverlappersIter.hasNext() && overlapperIsSameGene){
							
							Gene overlapper = endOverlappersIter.next();
							//If same orientation
							if(overlapper.getOrientation().equals(annotation.getOrientation())){
								//compare the end coordiantes of the gene
								if(is3p){
									if(!(overlapper.getOrientedEnd() == annotation.getOrientedEnd()))
										if(!BEDFileParser.isOverlapCompatible(annotation,overlapper, MIN_OVERLAP))
											overlapperIsSameGene = false;
								}
								else{
									if(!(overlapper.getOrientedStart() == annotation.getOrientedStart()) && !overlapper.overlaps(annotation))
										if(!BEDFileParser.isOverlapCompatible(annotation,overlapper, MIN_OVERLAP))
											overlapperIsSameGene = false;
								}
							}
						}
						if(!overlapperIsSameGene)
							break;
						// Because the extended region cannot overlap another annotation
						
						//No overlap so continue with scoring the region
						annotationEnd = new Gene(end);
						Map<String,double[]> tmpScores = getScoresForWindow(libDataModels,annotationEnd);
						double tmpEnrichment = getEnrichmentForWindow(libDataModels.keySet(),tmpScores);
						
						//For all peaks
						//If significant
						boolean overlaps = false;
						if(isSignificantInAllModels(libDataModels.keySet(), tmpScores)){
							//if overlaps any already added peak
							Set<Gene> tempPeaks = new HashSet<Gene>();
							tempPeaks.addAll(thisPeaks);
							for(Gene peak:thisPeaks){
								if(annotationEnd.overlaps(peak)){
									overlaps = true;
									//if better then replace
									
									if(tmpEnrichment>peak.getBedScore()){
										annotationEnd.setBedScore(tmpEnrichment); //set RPKM as score
										annotationEnd.setName(duplicateNameMap.get(annotation));
										tempPeaks.remove(peak);
										tempPeaks.add(annotationEnd);
									}
								}
							}
							thisPeaks = tempPeaks;
							if(!overlaps){
								annotationEnd.setBedScore(tmpEnrichment); //set RPKM as score
								annotationEnd.setName(duplicateNameMap.get(annotation));
								//does not overlap. Add
								thisPeaks.add(annotationEnd);	
							}
						}
						if(tmpEnrichment>bestEnrichment){
							bestEnrichment = tmpEnrichment;
							bestWindow = annotationEnd;
							bestWindowEnd = extend;
							insideGene = false;
							bestScores = tmpScores;
						}
						extend += (STEP);
					}
					
					
					String rowName = duplicateNameMap.get(annotation);
/*					if(fullGeneScoreFlag){
						if(is3p)
							fullGeneResultMatrix = getScoreFromWindowToGeneStart(libDataModels,bestWindow,annotation,bestWindowEnd,insideGene,fullGeneResultMatrix,rowName);
						else
							fullGeneResultMatrix = getScoreFromWindowToGeneEnd(libDataModels,bestWindow,annotation,bestWindowEnd,insideGene,fullGeneResultMatrix,rowName);
					}
*/					
					for(int i=0;i<cols.size();i++){
						//resultMatrix.set(duplicateNameMap.get(annotation), cols.get(i), scores[i]);
						resultMatrix.set(rowName, cols.get(i), bestScores.get(cols.get(i))[COUNT_SCORE]);
						pvalueMatrix.set(rowName, cols.get(i), bestScores.get(cols.get(i))[PVAL_SCORE]);
					}
					//SET THE NAME FOR THE WINDOW
					bestWindow.setName(duplicateNameMap.get(annotation));
					geneToWindowMap.put(duplicateNameMap.get(annotation), bestWindow);
					logger.debug("Added peak to geneToWindowMap as in the end "+duplicateNameMap.get(annotation)+" "+geneToWindowMap.containsKey(duplicateNameMap.get(annotation)));
					windows.add(bestWindow);
				}
			}
			else{
				logger.info("No data for " + chr);
				for(Gene annotation:annotations.get(chr)){
					for(int i=0;i<cols.size();i++){
						resultMatrix.set(duplicateNameMap.get(annotation), cols.get(i), 0.0);
						pvalueMatrix.set(duplicateNameMap.get(annotation), cols.get(i), 1.0);
/*						if(fullGeneScoreFlag){
							fullGeneResultMatrix.set(duplicateNameMap.get(annotation), cols.get(i), 0.0);
						}
*/					}	
				}
			}
		}
		/*
		 * Write to output File
		 */
		writeOutputFiles(output, resultMatrix,pvalueMatrix, annotations, windows);
/*		if(fullGeneScoreFlag){
			BufferedWriter bw =new BufferedWriter(new FileWriter(output+".fullgene.matrix")); 
			fullGeneResultMatrix.write(bw);
			bw.close();
		}
*/		
		if(collapseIsoforms){
			writeCollapsedOutputForMultiple(annotations,output,resultMatrix,geneToWindowMap);
			writeCollapsedOutputFilesLocalLambda(output,conditionMaps,alignmentFiles,is3p,annotations.keySet(),normalizedOutput);
		}
		
		if(normalizedOutput){
			
			writeNormalizedMatrix(output, resultMatrix);
/*			if(fullGeneScoreFlag){
				writeNormalizedMatrix(output+".fullgene.matrix", fullGeneResultMatrix);
			}
*/		}
/*		if(debugMode){
			BufferedWriter bedBw = new BufferedWriter(new FileWriter(output+".allPeaks.bed"));
			for(Gene gene: geneToPeaksMap.keySet()){
				for(Gene peak:geneToPeaksMap.get(gene)){
					bedBw.write(peak.toBED());
					//bedBw.append("\t"+duplicateNameMap.get(gene));
					bedBw.newLine();
				}			
			}
			bedBw.close();
		}*/
	}
	
	
	
	/**
	 * Calculates the scaling factor for each column of the specified matrix, using the normalization method used by DESeq.
	 * @param mat
	 * @return
	 */
	public static Map<String,Double> calculateNormalizationFactors(MatrixWithHeaders mat){
		
		int M = mat.getNumberColumns();
		int N = mat.getNumberRows();
		
		Map<String,Double> factors = new HashMap<String,Double>();
		double[] means = new double[N];
		//CALCULATE GEOMETRIC MEAN FOR ALL N GENES
		for(int i=0;i<N;i++){
			means[i] = 1.0;
			for(int j=0;j<M;j++){
				means[i] = means[i] * mat.get(i, j);
			}
			means[i] = Math.pow(means[i],((double)1.0/(double)M));
		}
		//For each sample
		for(String c:mat.getColumnNames()){
			
			double[] col = new double[N];
			for(int i=0;i<N;i++){
				col[i] = mat.get(i, c)/means[i]; 
			}
			factors.put(c, (double)1.0/calculateMedian(array2List(col)));
		}
		return factors;
		
	}
	
	/**
	 * This function converts an array of double to List of type Double
	 * @param double[]
	 * @return
	 */
	public static List<Double> array2List(double[] array){
		List<Double> lt=new ArrayList<Double>();

		for(int i=0;i<array.length;i++){
			lt.add(array[i]);
		}
	
		return lt;
	}
	
	/**
	 * Calculates median value of a list of values.
	 * @param values
	 * @return
	 */
	private static double calculateMedian(List<Double> values)
	{
	    List<Double> nonInfValues = new ArrayList<Double>();
	 
	    for(Double v: values){
	    	
	    	if(v.isInfinite()||v.isNaN()){
	    		//Dont add
	    	}
	    	else{
	    		//System.out.println(v);
	    		nonInfValues.add(v);
	    	}
	    }
	    Collections.sort(nonInfValues);
	    if (nonInfValues.size() % 2 == 1){
	    	double median = nonInfValues.get((nonInfValues.size()+1)/2-1);
	    	if(median==0.0)
	    		return (1.0);
	    	else
	    		return median;
	    }
	    else
	    {
	    	double lower;
	    	double upper;
	    	if(nonInfValues.size()==0){
	    		lower= 1.0;
	    		upper = 1.0;
	    	}
	    	else{
	    		lower = nonInfValues.get(nonInfValues.size()/2-1);
	    		upper = nonInfValues.get(nonInfValues.size()/2);
	    	}
			double median = (lower + upper) / 2.0;
			return median;
	    }
	}
	
	/**
	 * Returns the subannotation of length "length" at a distance "distanceFromTranscriptEnd" from the 3' end of the gene
	 * @param annotation
	 * @param length
	 * @param distanceFromTranscriptEnd
	 * @return
	 */
	public static Gene getSubAnnotationFromEnd(Gene annotation, int length, int distanceFromTranscriptEnd) {
		int annotationLength = annotation.getSize();
		Gene subAnnotation = null;
		int relativeStart = 0;
		int relativeEnd   = 0;
		if(annotationLength >= length + distanceFromTranscriptEnd) {
			relativeStart = annotationLength - (distanceFromTranscriptEnd + length);
			relativeEnd   = distanceFromTranscriptEnd;			
//			System.err.println("getSubannotation for " + annotation.getName() + " Annotation_length " + annotationLength + "  Window_size " + length + "("+annotation.getOrientation()+") distFromEnd: " + distanceFromTranscriptEnd + " relstart-relend " + relativeStart+"-"+relativeEnd);
			//TRIM FUNCTION IS STRAND SPECIFIC
			subAnnotation = new Gene(annotation.trim(relativeStart, relativeEnd));
			if(!subAnnotation.getOrientation().equals(annotation.getOrientation())) {
				System.err.println("Annotation and sub annotation orientations don't match");
			}
		}
		else{
			System.err.println("Annotation length "+annotationLength+" is less than length "+length +" distance from end "+ distanceFromTranscriptEnd);
		}
		return subAnnotation;
	}
	
	/**
	 * Returns the subannotation of length "length" at a distance "distanceFromTranscriptStart" from the 5' end of the gene
	 * @param annotation
	 * @param length
	 * @param distanceFromTranscriptEnd
	 * @return
	 */
	public static Gene getSubAnnotationFromStart(Gene annotation, int length, int distanceFromTranscriptStart) {
		int annotationLength = annotation.getSize();
		Gene subAnnotation = null;
		int relativeStart = 0;
		int relativeEnd   = 0;
		if(annotationLength >= length + distanceFromTranscriptStart) {
			relativeStart = distanceFromTranscriptStart;
			relativeEnd   = annotationLength - (distanceFromTranscriptStart + length);
			//System.err.println("getSubannotation for " + annotation.toUCSC() + " Annotation_length " + annotationLength + "  Window_size " + length + "("+annotation.getOrientation()+") distFromEnd: " + distanceFromTranscriptEnd + " relstart-relend " + relativeStart+"-"+relativeEnd);
			subAnnotation = new Gene(annotation.trim(relativeStart, relativeEnd));
			if(!subAnnotation.getOrientation().equals(annotation.getOrientation())) {
				System.err.println("Annotation and sub annotation orientations don't match");
			}
		}
		return subAnnotation;
	}

	/**
	 * This function writes the output to two files
	 * @param outFileName Name of the output file provided by user
	 * @param resultMatrix matrix to be written to the output file
	 * @param allGenes List of all genes in annotation file
	 * @throws IOException
	 */
	private static void writeOutputFiles(String outFileName,MatrixWithHeaders resultMatrix, Map<String,Collection<Gene>> annotations) throws IOException{
		
		// Name of the bed file is outFileName followed by .bed
		String bedFileName = outFileName+".bed";
		
		/*
		 * Write the expression file
		 */
		BufferedWriter outBw = new BufferedWriter(new FileWriter(outFileName));
		resultMatrix.write(outBw);
		outBw.close();
		/*
		 * Write the bed File
		 */
		BufferedWriter bedBw = new BufferedWriter(new FileWriter(bedFileName));
		for(String chr:annotations.keySet()){
			for(Gene gene: annotations.get(chr)){
				bedBw.write(gene.toBED());
				bedBw.append("\t"+duplicateNameMap.get(gene));
				bedBw.newLine();
			}
		}
		bedBw.close();
		
	}
	
	
	/**
	 * This function writes the output to three files
	 * @param outFileName Name of the output file provided by user
	 * @param resultMatrix matrix to be written to the output file
	 * @param allGenes List of all genes in annotation file
	 * @throws IOException
	 */
	private static void writeOutputFiles(String outFileName,MatrixWithHeaders resultMatrix,MatrixWithHeaders pvaluesMatrix, Map<String,Collection<Gene>> annotations,List<Gene> windows) throws IOException{
		
		writeOutputFiles(outFileName,resultMatrix,annotations);
		
		if(debugMode){
			BufferedWriter outBw = new BufferedWriter(new FileWriter(outFileName+".pvalues.matrix"));
			pvaluesMatrix.write(outBw);
			outBw.close();
			// Name of the second bed file is outFileName followed by .windows.bed
			String bedFileName = outFileName+".windows.bed";		
			//Write the bed File
			BufferedWriter bedBw = new BufferedWriter(new FileWriter(bedFileName));
			for(Gene gene: windows){
				bedBw.write(gene.toString());
				//bedBw.append("\t"+duplicateNameMap.get(gene));
				bedBw.newLine();
			}
			bedBw.close();
		}
		
	}

	
	/**
	 * Returns the reconstructed annotation set based on the reconstruction provided
	 * @param completeReconstruction
	 * @param annotations
	 * @param annotationFile
	 * @return
	 * @throws IOException
	 */
	public static BEDFileParser reconstructGeneEnds(BEDFileParser completeReconstruction, BEDFileParser annotations,String annotationFile) throws IOException{
		
		/*
		 *  ARGUMENTS:
		 *  double minSpliceFrequency
		 *  int[] fixedWidth
		 *  String chrToSegment
		 *  boolean trimEnds	FALSE FOR NOW
		 *  boolean filterCanonical
		 *  String sequenceFile
		 *  double alpha
		 *  int chrStart
		 *  int chrEnd
		 *  double trimQuantile
		 */
		//BEDAnnotationCollection completeReconstruction = completeCDAM.segmentChromosome(0.1, fixedWidth, null, true, true, sequenceFile, 0.05, 0, Integer.MAX_VALUE, 0.1);
		completeReconstruction.makeGenes(); //This avoids having to compare to many different isoforms of different sizes usually present in reconstructed transcriptomes.
		completeReconstruction.writeFullBed("reconstructions.togenes.bed");
		annotations.equalizeTranscriptEnds(completeReconstruction);
		annotations.writeFullBed(annotationFile+".3p5pAdjusted.bed");
		
		return annotations;
	}
	
	/**
	 * Writes the collapsed gene output for a single endSeq sample
	 * @param annotationCollection
	 * @param outputFile
	 * @param resultMatrix
	 * @param geneToWindowMap
	 * @throws IOException
	 */
	
	//TODO: NOT WORKING
	private static void writeCollapsedOutputForSingle(Map<String,Collection<Gene>> annotationCollection,String outputFile,MatrixWithHeaders resultMatrix,Map<String, Annotation> geneToWindowMap) throws IOException{

		logger.info("Collapsing isoforms to genes:");
		String collapsedOutput = outputFile+".collapsed.table";
		BufferedWriter bw1 = new BufferedWriter(new FileWriter(outputFile+".collapsedGenes.bed"));
		BufferedWriter bwW = null;
		if(debugMode)
			bwW =  new BufferedWriter(new FileWriter(outputFile+".collapsedGenes.windows.bed"));
		BufferedWriter bw = new BufferedWriter(new FileWriter(collapsedOutput));
		BufferedWriter bwDist = new BufferedWriter(new FileWriter(collapsedOutput+".distance"));
		
		for(String colName:resultMatrix.getColumnNames()){
			if(!colName.equals(BEST_DIST_TO_END_COL))
				bw.write("\t"+colName);
		}
		bw.newLine();
		
		for(String chr:collapsedGenes.keySet()){
			
			IntervalTree<Gene> tree=collapsedGenes.get(chr);
			Iterator<Gene> geneIt=tree.valueIterator();
			while(geneIt.hasNext()) {
				
				//For each collapsed gene				
				Gene gene = geneIt.next();
				
				double scores[] = new double[resultMatrix.columnDimension()];
				for(int i=0;i<scores.length;i++){
					scores[i] = 0.0;//resultMatrix.get(duplicateNameMap.get(gene), i);
				}
					
				List<Annotation> prevPeaks = new ArrayList<Annotation>();
				Collection<Annotation> exons = new ArrayList<Annotation>();
				
				//For all constituent genes
				for(Gene isoform:collapsedGeneMap.get(gene)){
//					logger.info(gene.getName()+" and "+isoform.getName());
					if(geneToWindowMap.containsKey(duplicateNameMap.get(isoform))){
						if(collapsedGeneMap.get(gene).size()==1){
							prevPeaks.add(geneToWindowMap.get(duplicateNameMap.get(isoform)));
							exons.add(geneToWindowMap.get(duplicateNameMap.get(isoform)));
							for(int i=0;i<scores.length;i++){
								scores[i] = resultMatrix.get(duplicateNameMap.get(isoform), i); 
							}
						}
						else{
							boolean overLapsPeak = false;
							for(Annotation peak:prevPeaks){
								if(peak.overlaps(geneToWindowMap.get(duplicateNameMap.get(isoform))))
									overLapsPeak = true;
							}				
							if(!overLapsPeak){
								prevPeaks.add(geneToWindowMap.get(duplicateNameMap.get(isoform)));
								exons.add(geneToWindowMap.get(duplicateNameMap.get(isoform)));
								for(int i=0;i<scores.length;i++)
									scores[i] += resultMatrix.get(duplicateNameMap.get(isoform), i);
							}
						}								
					}
					else{
						logger.info("Does not contain");
					}
				}
				bw.write(gene.getName()+"\t");
				for(int i=0;i<scores.length;i++)
					bw.write(scores[i]+"\t");
				bw.newLine();
				bw1.write(gene.toBED());
		  		bw1.newLine();
		  		Gene windows = (makeGene(exons,gene.getName(),gene.getOrientation()));
		  		int distance = getDistance(gene,windows);
		  		bwDist.write(gene.getName()+"\t"+distance+"\n");
		  		if(debugMode)
		  			bwW.write(windows.toBED()+"\n");
			}
		}	
		bw.close();
		bw1.close();
		if(debugMode)
			bwW.close();

	}

	
	/**
	 * This function calculates all scores for all given alignment models for a given window
	 * @param libDataModels
	 * @param annotationEnd
	 * @return Vector of size libDataModels with enrichment of each 
	 * @throws IOException
	 */
	private static Map<String,double[]> getScoresForWindow(Map<String, AlignmentModel> libDataModels,Gene annotationEnd) throws IOException{
		
		Map<String,double[]> result = new HashMap<String,double[]>();
		if(oppositeStrand)
			annotationEnd = reverseStrand(annotationEnd);
		
		for(String modelName:libDataModels.keySet()){
			
			if(libDataModels.get(modelName).containsReference(annotationEnd.getChr())){	
				
				/*
				 * Returns an array of scores
				 * [0] = count
			 	 * [1] = RPKM
				 * [2] = RPK
			 	 * [3] = region total
			 	 * [4] = total
			 	 * [5] = Annotation length
			 	 * [6] = scan p value
				 * [7] = global lambda
			 	 * [8] = local lambda
			 	 * [9] = region length  
			 	 */		
				double [] scores = new ScanStatisticScore(libDataModels.get(modelName), annotationEnd,false).getScores();
				
				result.put(modelName, scores);
			}
		}
		
		return result;
	}
	
	/**
	 * This function calculates all scores for all given alignment models for a given window
	 * @param libDataModels
	 * @param annotationEnd
	 * @param geneScores 
	 * @param geneLength 
	 * @return Vector of size libDataModels with enrichment of each 
	 * @throws IOException
	 */
	private static Map<String,double[]> getStartEndScoresForWindow(Map<String, AlignmentModel> libDataModels,Gene annotationEnd,boolean is3p, int geneLength, Map<String, Double> geneScores) throws IOException{
		
		Map<String,double[]> result = new HashMap<String,double[]>();
		
		for(String modelName:libDataModels.keySet()){
			
			if(libDataModels.get(modelName).containsReference(annotationEnd.getChr())){	
				
				AbstractWindowScore score = null;
				if((is3p && !oppositeStrand)|| (!is3p && oppositeStrand))
					score = new EndInWindowCountScore(libDataModels.get(modelName), annotationEnd);
				
				else
					score = new StartInWindowCountScore(libDataModels.get(modelName), annotationEnd);
				/*
				 * Returns an array of scores
				 * [0] = count
			 	 * [1] = RPKM
				 * [2] = RPK
			 	 * [3] = region total
			 	 * [4] = total
			 	 * [5] = Annotation length
			 	 * [6] = scan p value
				 * [7] = global lambda
			 	 * [8] = local lambda
			 	 * [9] = region length  
			 	 */		
				double [] scores = new BasicScanStatisticScore(libDataModels.get(modelName), annotationEnd,score.getScore(), false, geneLength, geneScores.get(modelName)).getScores();
				
				result.put(modelName, scores);
			}
		}
		
		return result;
	}
	/**
	 * This function calculates the best enrichment score for a given window
	 * @param libDataModels
	 * @param annotationEnd
	 * @return Vector of size libDataModels with enrichment of each 
	 * @throws IOException
	 */
	private static double getEnrichmentForWindow(Collection<String> allModelNames,Map<String,double[]> modelScores) throws IOException{
		
		double sum = 0.0;
		for(String modelName:allModelNames){
			if(modelScores.containsKey(modelName)){
				double[] scores = modelScores.get(modelName);
				double value = (scores[COUNT_SCORE]/scores[ANNOTATION_LENGTH])/scores[LOCAL_LAMBDA];
				if(!Double.isNaN(value))
					sum += value;
			}
		}
		
		//return max;
		return sum;
	}
	
	
	/**
	 * This function returns true is said window is significant in all models
	 * @param libDataModels
	 * @param annotationEnd
	 * @return Vector of size libDataModels with enrichment of each 
	 * @throws IOException
	 */
	private static boolean isSignificantInAllModels(Collection<String> allModelNames,Map<String,double[]> modelScores) throws IOException{
		
		for(String modelName:allModelNames){
			if(modelScores.containsKey(modelName)){
				if(modelScores.get(modelName)[PVAL_SCORE]>PVAL_THRESHOLD)
					return false;
			}
			else
				return false;
		}
		return true;
	}
	
	
	private static void writeCollapsedOutputForMultiple(Map<String,Collection<Gene>> annotations,String outputFile,MatrixWithHeaders resultMatrix,Map<String, Annotation> geneToWindowMap) throws IOException{

		logger.info("Collapsing isoforms to genes:");
		BufferedWriter bw1 = new BufferedWriter(new FileWriter(outputFile+".collapsedGenes.bed"));		
		
		for(String chr:collapsedGenes.keySet()){
			
			IntervalTree<Gene> tree=collapsedGenes.get(chr);
			Iterator<Gene> geneIt=tree.valueIterator();
			while(geneIt.hasNext()) {
				
				//For each collapsed gene				
				Gene gene = geneIt.next();
					
				bw1.write(gene.toBED());
		  		bw1.newLine();
			}
		}	
		bw1.close();
		
	}
	
	/**
	 * Will calculate isoform and gene level expression for the collapsed genes
	 * @param outputFile
	 * @param alignmentFiles 
	 * @param conditionMaps 
	 * @param is3p 
	 * @param set 
	 * @throws IOException 
	 */
	private static void writeCollapsedOutputFilesLocalLambda(String outputFile, Map<String, Collection<String>> conditionMaps, Map<String, String> alignmentFiles, boolean is3p, Set<String> set,boolean normalizedOutput) throws IOException{
		
		logger.info("Calculating gene level expression");
		// The columns of the matrix will be each alignment file name
		List<String> cols = new ArrayList<String>();
		for(String name: conditionMaps.keySet())
			for(String sample: conditionMaps.get(name))
				cols.add(sample);
				
		List<String> geneNames = new ArrayList<String>();
		for(Gene gene:collapsedGeneMap.keySet())
			geneNames.add(gene.getName());
		
		// Initialize the AlignmentModels
		Map<String,AlignmentModel> libDataModels = new HashMap<String,AlignmentModel>();
		for(String ss: alignmentFiles.keySet()){
			AlignmentModel libmodel = new AlignmentModel(alignmentFiles.get(ss), null, new ArrayList<Predicate<Alignment>>(), pairedFlag,strand,true,maskedRegionFile);
			
			AlignmentModel model = new AlignmentModel(alignmentFiles.get(ss), new TranscriptInGenomicSpace(libmodel.getRefSequenceLengths()), new ArrayList<Predicate<Alignment>>(), pairedFlag,strand,true,maskedRegionFile);
			if(removePCRDuplicatesFlag)
				model.addFilter(new PCRDuplicateFilter());
			if(filterMultimappers)
				model.addFilter(new UniqueMappedReadsFilter());
			libDataModels.put(ss, model);
		}	
		
		//To report all peaks - not just the best one - FOR CHROMOSOME
		BufferedWriter bw = new BufferedWriter(new FileWriter(outputFile+".collapsed.significantEnds"));
		
		BufferedWriter bwBed = null;
		if(debugMode){
			bwBed = new BufferedWriter(new FileWriter(outputFile+".collapsed.significantEnds.bed"));
		}
		BufferedWriter bwFracs = new BufferedWriter(new FileWriter(outputFile+".collapsed.significantEnds.fractions"));
		bw.write("Gene\tPeak\t");
		bwFracs.write("Gene\tPeak\t");
		for(String sample:cols){
			bw.write(sample+"\t");
			bwFracs.write(sample+"\t");
		}
		bw.newLine();
		bwFracs.newLine();

		//Sum of all the significant peaks
		MatrixWithHeaders fullgeneMatrix = new MatrixWithHeaders(geneNames,cols);
		//Sum of all significant peaks divided by the number of peaks
		List<String> peaks = new ArrayList<String>();
		String NUM_PEAKS = "NUM_PEAKS";
		peaks.add(NUM_PEAKS);
		MatrixWithHeaders numPeaksMatrix = new MatrixWithHeaders(geneNames,peaks);
		
		//Map of each gene to its best peak
		Map<String,Annotation> geneToBestWindowMap = new HashMap<String,Annotation>();
		
		//Best peak 
		MatrixWithHeaders genesMatrix = new MatrixWithHeaders(geneNames,cols);
		MatrixWithHeaders genespvalueMatrix = new MatrixWithHeaders(geneNames,cols);

		List<Gene> windows = new ArrayList<Gene>();	
		//For each chromosome
		for(String chr:set){
			
			// If the alignment data has data from that chromosome
			boolean dataForChr = false;
			for(AlignmentModel model:libDataModels.values()){
				//If any alignment file has data for the chromosome, 
				//		set flag to true and exit loop
				if(model.containsReference(chr)){
					dataForChr = true;
					break;
				}
			}
			
			if(dataForChr){
				IntervalTree<Gene> tree=collapsedGenes.get(chr);
				Iterator<Gene> geneIt=tree.valueIterator();
				logger.info("Processing "+chr);
				while(geneIt.hasNext()) {
					//For each collapsed gene				
					Gene gene = geneIt.next();
					Map<Gene,Map<String,Double>> sigPeakToScores = new HashMap<Gene,Map<String,Double>>();
					Gene geneEnd = gene;
					int geneLength = gene.getSize();
					Map<String,Double> geneScores = new HashMap<String,Double>();
					
					for(String modelName:libDataModels.keySet()){
						if((is3p && !oppositeStrand)||(!is3p && oppositeStrand))
							geneScores.put(modelName,new EndInWindowCountScore(libDataModels.get(modelName), gene).getScore());
						else
							geneScores.put(modelName,new StartInWindowCountScore(libDataModels.get(modelName), gene).getScore());
					}
					//If the length of the annotated transcript > window size being analyzed,
					//get sub annotation for window length
					if(geneLength>window){
						geneEnd = is3p? getSubAnnotationFromEnd(gene,window,0) : getSubAnnotationFromStart(gene,window,0);
						if(geneEnd == null){
							logger.warn("Annotation end for " + (geneLength - window) + "-" + geneLength + " --> " + gene.toBED() + " was null.");
						}
					}
					
					//Get the scores
					Map<String,double[]> bestScores = getStartEndScoresForWindow(libDataModels,geneEnd,is3p,geneLength,geneScores);
					//For the bestWindow
					double bestEnrichment = getEnrichmentForWindow(libDataModels.keySet(),bestScores);
					Gene bestWindow = geneEnd;
									
					Gene bestSignificantWindow = null;
					double bestSignificantEnrichment = Double.MIN_VALUE;
					
//					boolean aboveThreshold = (Statistics.max(geneScores.values())>=10) ? true:false;
//					if(aboveThreshold)
//						logger.info("Flag set for: "+gene.getName()+" with max count "+Statistics.max(geneScores.values()));
					//Is the window significant
					if(isSignificantWindow(libDataModels.keySet(),bestScores)){// && aboveThreshold){
						geneEnd.setBedScore(bestEnrichment); //set RPKM as score
						geneEnd.setName(gene.getName());
						sigPeakToScores.put(geneEnd, new HashMap<String,Double>());
						for(String column:cols)
							sigPeakToScores.get(geneEnd).put(column, bestScores.get(column)[USE_SCORE]);
						bestSignificantEnrichment = bestEnrichment;
						bestSignificantWindow = geneEnd;
					}
					
					//Parse the gene
					int intoGene = STEP;
					while((intoGene<(geneLength - window))){
						//get annotation for region of length window, "intoGene" length from end of transcript
						geneEnd = is3p? getSubAnnotationFromEnd(gene, window, intoGene): getSubAnnotationFromStart(gene, window, intoGene);
						
						if(geneEnd !=null){
							Map<String,double[]> tmpScores = getStartEndScoresForWindow(libDataModels,geneEnd,is3p,geneLength,geneScores);
							double tmpEnrichment = getEnrichmentForWindow(libDataModels.keySet(),tmpScores);
							//For all peaks
							//If significant
							boolean overlaps = false;
							//Is the window significant
							//Add to thisPeaks - thats it
							if(isSignificantWindow(libDataModels.keySet(),tmpScores)){//  && aboveThreshold){
									geneEnd.setBedScore(tmpEnrichment); //set RPKM as score
									geneEnd.setName(gene.getName());
									sigPeakToScores.put(geneEnd, new HashMap<String,Double>());
									for(String sample:cols)
										sigPeakToScores.get(geneEnd).put(sample, tmpScores.get(sample)[USE_SCORE]);
									if(tmpEnrichment>bestSignificantEnrichment){
										bestSignificantEnrichment = tmpEnrichment;
										bestSignificantWindow = geneEnd;
									}
							}							
							//For best peak
							if(tmpEnrichment>bestEnrichment){
								bestEnrichment = tmpEnrichment;
								bestWindow = geneEnd;
								bestScores = tmpScores;
							}
						}
						intoGene += STEP;
					}
					
					//EXTENSION
					int extend = STEP;
					while(extend < maxExtension) {
						Annotation end = null;
						long chrEnd = 0;
						for(String s:libDataModels.keySet())
							chrEnd = libDataModels.get(s).getRefSequenceLength(chr);
						
						if((is3p && gene.getOrientation().equals(Strand.NEGATIVE))||(!is3p && gene.getOrientation().equals(Strand.POSITIVE))){
							if((gene.getStart() - extend)<0){
								logger.info(gene.getStart()+" - "+extend +" is less than 0");
								break;
							}
							end = new BasicAnnotation(gene.getChr(), gene.getStart() - extend, gene.getStart() - (extend-window),gene.getOrientation());
						}
						else{
							if((gene.getEnd() + extend)> chrEnd){
								logger.info(gene.getEnd()+" + "+extend +" is more than "+chrEnd);
								break;
							}
							end = new BasicAnnotation(gene.getChr(), gene.getEnd() + (extend - window), gene.getEnd() + extend,gene.getOrientation());

						}
						// Get an interval tree for all/any exons that overlap with the extended region
						Iterator<Gene> endOverlappersIter = tree.overlappingValueIterator(end.getStart(), end.getEnd());
						//While there is an overlap with a gene and gene is same gene
						boolean overlapperIsSameGene = true;
						while(endOverlappersIter.hasNext() && overlapperIsSameGene){
							
							Gene overlapper = endOverlappersIter.next();
							//If same orientation
							if(overlapper.getOrientation().equals(gene.getOrientation())){
								//compare the end coordiantes of the gene
								if(is3p){
									if(!(overlapper.getOrientedEnd() == gene.getOrientedEnd()))
										if(!BEDFileParser.isOverlapCompatible(gene,overlapper, MIN_OVERLAP))
											overlapperIsSameGene = false;
								}
								else{
									if(!(overlapper.getOrientedStart() == gene.getOrientedStart()) && !overlapper.overlaps(gene))
										if(!BEDFileParser.isOverlapCompatible(gene,overlapper, MIN_OVERLAP))
											overlapperIsSameGene = false;
								}
							}
						}
						if(!overlapperIsSameGene)
							break;
						geneEnd = new Gene(end);
						//geneEnd.setOrientation(gene.getOrientation());
						//No overlap so continue with scoring the region
						Map<String,double[]> tmpScores = getStartEndScoresForWindow(libDataModels,geneEnd,is3p,geneLength,geneScores);
						double tmpEnrichment = getEnrichmentForWindow(libDataModels.keySet(),tmpScores);
						
						//For all peaks
						//If significant
						boolean overlaps = false;
						//Is the window significant
						if(isSignificantWindow(libDataModels.keySet(),tmpScores)){// && aboveThreshold){
								geneEnd.setBedScore(tmpEnrichment); //set RPKM as score
								geneEnd.setName(gene.getName());
								sigPeakToScores.put(geneEnd, new HashMap<String,Double>());
								for(String sample:cols)
									sigPeakToScores.get(geneEnd).put(sample, tmpScores.get(sample)[USE_SCORE]);
								if(tmpEnrichment>bestSignificantEnrichment){
									bestSignificantEnrichment = tmpEnrichment;
									bestSignificantWindow = geneEnd;
								}
						}
						if(tmpEnrichment>bestEnrichment){
							bestEnrichment = tmpEnrichment;
							bestWindow = geneEnd;
							bestScores = tmpScores;
						}
						extend += (STEP);
					}
					
					String rowName = gene.getName();
					
					//Best window
					for(int i=0;i<cols.size();i++){
						//resultMatrix.set(duplicateNameMap.get(annotation), cols.get(i), scores[i]);
						genesMatrix.set(rowName, cols.get(i), bestScores.get(cols.get(i))[COUNT_SCORE]);
						genespvalueMatrix.set(rowName, cols.get(i), bestScores.get(cols.get(i))[PVAL_SCORE]);
					}
					//SET THE NAME FOR THE WINDOW
					bestWindow.setName(gene.getName());
					geneToBestWindowMap.put(gene.getName(), bestWindow);
					logger.debug("Added peak to geneToWindowMap as in the end "+gene.getName()+" "+geneToBestWindowMap.containsKey(gene.getName()));
					windows.add(bestWindow);

					//FIND ALL SIGNIFICANT PEAKS
					int numPeaks = 0;
					if(sigPeakToScores.isEmpty() || bestSignificantWindow==null){
						logger.info("No significant for: "+gene.getName()+"\t");
						for(String sample:cols){
							fullgeneMatrix.set(gene.getName(), sample, genesMatrix.get(gene.getName(), sample));
						}
					}else{
						Annotation nextBestWindow = bestSignificantWindow;
						Annotation thisBestWindow = bestSignificantWindow;
						Map<Gene,Map<String,Double>> peakScores = new HashMap<Gene,Map<String,Double>>();
						peakScores.putAll(sigPeakToScores);
						//For fractions
						Map<Annotation,Map<String,Double>> thispeakScores = new HashMap<Annotation,Map<String,Double>>();

						boolean peakFound = true;
						while(!sigPeakToScores.isEmpty() && peakFound){	
							thisBestWindow = nextBestWindow;
							
							//WRITE THE NEXT BEST WINDOW FIRST
							//WRITE TO THE FILE
							logger.info(gene.getName()+"\t");
							logger.info("Peak : "+thisBestWindow.toBED()+"\t"+sigPeakToScores.containsKey(thisBestWindow));
							
							bw.write(gene.getName()+"\t"+thisBestWindow.toUCSC());
							if(debugMode){
								bwBed.write(thisBestWindow.toBED()+"\n");
							}
							thispeakScores.put(thisBestWindow, new HashMap<String,Double>());
							for(String sample:cols){
								bw.write("\t"+sigPeakToScores.get(thisBestWindow).get(sample));
								double value = fullgeneMatrix.get(gene.getName(), sample)+sigPeakToScores.get(thisBestWindow).get(sample);
								fullgeneMatrix.set(gene.getName(), sample, value);	
								thispeakScores.get(thisBestWindow).put(sample, sigPeakToScores.get(thisBestWindow).get(sample));
							}
							bw.write("\n");		
							
							double enrichment = Double.MIN_VALUE;
							peakFound = false;
							for(Gene peak:sigPeakToScores.keySet()){
								if(peak.overlaps(thisBestWindow))
									peakScores.remove(peak);
								else{
									if(peak.getBedScore()>enrichment){
										peakFound=true;
										nextBestWindow = peak;
										enrichment = peak.getBedScore();
									}
								}
							}
							numPeaks++;
							
							sigPeakToScores = new HashMap<Gene,Map<String,Double>>();
							sigPeakToScores.putAll(peakScores);
						}
	
						//Write to fractions files
						for(Annotation p:thispeakScores.keySet()){
							bwFracs.write(gene.getName()+"\t"+p.toUCSC());
							for(String sample:cols){
								double value = fullgeneMatrix.get(gene.getName(), sample);
								bwFracs.write("\t"+(thispeakScores.get(p).get(sample)/value));
							}
							bwFracs.newLine();
						}
					}
					
					//the per peak matrix
					numPeaksMatrix.set(gene.getName(), NUM_PEAKS, numPeaks);
					
				}
			}
			else{
				logger.info("No data for " + chr);
				IntervalTree<Gene> tree=collapsedGenes.get(chr);
				Iterator<Gene> geneIt=tree.valueIterator();
				logger.info("Processing "+chr);
				while(geneIt.hasNext()) {
					Gene gene = geneIt.next();
					for(int i=0;i<cols.size();i++){
						genesMatrix.set(gene.getName(), cols.get(i), 0.0);
						genespvalueMatrix.set(gene.getName(), cols.get(i), 1.0);
						fullgeneMatrix.set(gene.getName(), cols.get(i), 0.0);
						numPeaksMatrix.set(gene.getName(), NUM_PEAKS, 0);
					}
				}
			}
		}
		
		//TODO: write out the peaks in geneToSignificantWindowsMap
		bw.close();
		if(debugMode){
			bwBed.close();
		}
		bwFracs.close();
		
		//Best isoform
		bw =new BufferedWriter(new FileWriter(outputFile+".collapsed.bestIsoform.scores")); 
		genesMatrix.write(bw);
		bw.close();
		if(normalizedOutput)
			writeNormalizedMatrix(outputFile+".collapsed.bestIsoform.scores", genesMatrix);
		if(debugMode){
			bw =new BufferedWriter(new FileWriter(outputFile+".collapsed.bestIsoform.pvalues.matrix")); 
			genespvalueMatrix.write(bw);
			bw.close();	
		////Sum of all significant peaks divided by the number of peaks
			bw =new BufferedWriter(new FileWriter(outputFile+".collapsed.numPeaks")); 
			numPeaksMatrix.write(bw);
			bw.close();		
			
		}
		if(fullGeneScoreFlag){
			//Sum of all the significant peaks
			bw =new BufferedWriter(new FileWriter(outputFile+".fullgene.matrix")); 
			fullgeneMatrix.write(bw);
			bw.close();
			if(normalizedOutput)
				writeNormalizedMatrix(outputFile+".fullgene.matrix", fullgeneMatrix);
		}
		
		if(debugMode){
			BufferedWriter bedBw = new BufferedWriter(new FileWriter(outputFile+".collapsed.bestIsoform.windows.bed"));
			for(Gene gene: windows){
				bedBw.write(gene.toString());
				bedBw.newLine();
			}
			bedBw.close();
		}
	}
	

	/**
	 * Will calculate isoform and gene level expression for the collapsed genes
	 * @param outputFile
	 * @param alignmentFiles 
	 * @param conditionMaps 
	 * @param is3p 
	 * @param set 
	 * @throws IOException 
	 */
	private static void writeCollapsedOutputFilesOld(String outputFile, Map<String, Collection<String>> conditionMaps, Map<String, String> alignmentFiles, boolean is3p, Set<String> set,boolean normalizedOutput) throws IOException{
		
		logger.info("Calculating gene level expression");
		// The columns of the matrix will be each alignment file name
		List<String> cols = new ArrayList<String>();
		for(String name: conditionMaps.keySet())
			for(String sample: conditionMaps.get(name))
				cols.add(sample);
				
		List<String> geneNames = new ArrayList<String>();
		for(Gene gene:collapsedGeneMap.keySet())
			geneNames.add(gene.getName());
		
		// Initialize the AlignmentModels
		Map<String,AlignmentModel> libDataModels = new HashMap<String,AlignmentModel>();
		for(String ss: alignmentFiles.keySet()){
			AlignmentModel libmodel = new AlignmentModel(alignmentFiles.get(ss), null, new ArrayList<Predicate<Alignment>>(), pairedFlag,strand,true,maskedRegionFile);
			
			AlignmentModel model = new AlignmentModel(alignmentFiles.get(ss), new TranscriptInGenomicSpace(libmodel.getRefSequenceLengths()), new ArrayList<Predicate<Alignment>>(), pairedFlag,strand,true,maskedRegionFile);
			if(removePCRDuplicatesFlag)
				model.addFilter(new PCRDuplicateFilter());
			if(filterMultimappers)
				model.addFilter(new UniqueMappedReadsFilter());
			libDataModels.put(ss, model);
		}	
		
		//To report all peaks - not just the best one - FOR CHROMOSOME
		BufferedWriter bw = new BufferedWriter(new FileWriter(outputFile+".collapsed.significantEnds"));
		BufferedWriter bwBed = new BufferedWriter(new FileWriter(outputFile+".collapsed.significantEnds.bed"));
		bw.write("Gene\tPeak\t");
		for(String sample:cols){
			bw.write(sample+"\t");
		}
		bw.newLine();

		//Sum of all the significant peaks
		MatrixWithHeaders fullgeneMatrix = new MatrixWithHeaders(geneNames,cols);
		//Sum of all significant peaks divided by the number of peaks
		List<String> peaks = new ArrayList<String>();
		String NUM_PEAKS = "NUM_PEAKS";
		peaks.add(NUM_PEAKS);
		MatrixWithHeaders numPeaksMatrix = new MatrixWithHeaders(geneNames,peaks);
		
		//Map of each gene to its peak
		Map<String,Annotation> geneToBestWindowMap = new HashMap<String,Annotation>();
		
		//Best peak 
		MatrixWithHeaders genesMatrix = new MatrixWithHeaders(geneNames,cols);
		MatrixWithHeaders genespvalueMatrix = new MatrixWithHeaders(geneNames,cols);

		List<Gene> windows = new ArrayList<Gene>();	
		//For each chromosome
		for(String chr:set){
			
			// If the alignment data has data from that chromosome
			boolean dataForChr = false;
			for(AlignmentModel model:libDataModels.values()){
				//If any alignment file has data for the chromosome, 
				//		set flag to true and exit loop
				if(model.containsReference(chr)){
					dataForChr = true;
					break;
				}
			}
			
			if(dataForChr){
				IntervalTree<Gene> tree=collapsedGenes.get(chr);
				Iterator<Gene> geneIt=tree.valueIterator();
				logger.info("Processing "+chr);
				while(geneIt.hasNext()) {
					//For each collapsed gene				
					Gene gene = geneIt.next();
					Set<Gene> thisPeaks = new HashSet<Gene>();
					Map<Gene,Map<String,Double>> sigPeakToScores = new HashMap<Gene,Map<String,Double>>();
					Gene geneEnd = gene;
					int geneLength = gene.getSize();
					
					//If the length of the annotated transcript > window size being analyzed,
					//get sub annotation for window length
					if(geneLength>window){
						geneEnd = is3p? getSubAnnotationFromEnd(gene,window,0) : getSubAnnotationFromStart(gene,window,0);
						if(geneEnd == null){
							logger.warn("Annotation end for " + (geneLength - window) + "-" + geneLength + " --> " + gene.toBED() + " was null.");
						}
					}
					
					//Get the scores
					Map<String,double[]> bestScores = getScoresForWindow(libDataModels,geneEnd);
					//For the bestWindow
					double bestEnrichment = getEnrichmentForWindow(libDataModels.keySet(),bestScores);
					int bestWindowEnd = 0;
					boolean insideGene = true;
					Gene bestWindow = geneEnd;
					
					Map<String,double[]> primaryScores = null;
					Map<String,double[]> secondaryScores = null;
					double primaryEnrichment=0;
					double secondaryEnrichment=0;
					
					//Is the window significant
					if(isSignificantWindow(libDataModels.keySet(),bestScores)){
						geneEnd.setScore(bestEnrichment); //set RPKM as score
						geneEnd.setName(gene.getName());
						thisPeaks.add(geneEnd);
						
						sigPeakToScores.put(geneEnd, new HashMap<String,Double>());
						for(String sample:cols)
							sigPeakToScores.get(geneEnd).put(sample, bestScores.get(sample)[USE_SCORE]);
						
						primaryEnrichment = bestEnrichment;
						primaryScores = bestScores;
						//if(!geneToSignificantWindowsMap.containsKey(gene))
						//	geneToSignificantWindowsMap.put(gene, new TreeSet<Gene>());
						//geneToSignificantWindowsMap.get(gene).add(geneEnd);
					}
					
					//Parse the gene
					int intoGene = STEP;
					while((intoGene<(geneLength - window))){
						//get annotation for region of length window, "intoGene" length from end of transcript
						geneEnd = is3p? getSubAnnotationFromEnd(gene, window, intoGene): getSubAnnotationFromStart(gene, window, intoGene);
						
						if(geneEnd !=null){
							Map<String,double[]> tmpScores = getScoresForWindow(libDataModels,geneEnd);
							double tmpEnrichment = getEnrichmentForWindow(libDataModels.keySet(),tmpScores);
							//For all peaks
							//If significant
							boolean overlaps = false;
							//Is the window significant
							if(isSignificantWindow(libDataModels.keySet(),tmpScores)){
								//if overlaps any already added peak
								Set<Gene> tempPeaks = new HashSet<Gene>();
								tempPeaks.addAll(thisPeaks);
								for(Gene peak:thisPeaks){
									if(geneEnd.overlaps(peak)){
										overlaps = true;
										//if better then replace
										
										if(tmpEnrichment>peak.getBedScore()){
											geneEnd.setBedScore(tmpEnrichment); //set RPKM as score
											geneEnd.setName(gene.getName());
											tempPeaks.remove(peak);
											sigPeakToScores.remove(peak);
											tempPeaks.add(geneEnd);
											sigPeakToScores.put(geneEnd, new HashMap<String,Double>());
											for(String column:cols)
												sigPeakToScores.get(geneEnd).put(column, tmpScores.get(column)[USE_SCORE]);
										}
									}
								}
								thisPeaks = tempPeaks;
								if(!overlaps){
									geneEnd.setBedScore(tmpEnrichment); //set RPKM as score
									geneEnd.setName(gene.getName());
									//does not overlap. Add
									thisPeaks.add(geneEnd);	
									sigPeakToScores.put(geneEnd, new HashMap<String,Double>());
									for(String sample:cols)
										sigPeakToScores.get(geneEnd).put(sample, tmpScores.get(sample)[USE_SCORE]);
								}
								
								//Better than primary, so this becomes primary and primary becomes secondary
								if(tmpEnrichment>primaryEnrichment){
									secondaryEnrichment = primaryEnrichment;
									secondaryScores = primaryScores;
									primaryEnrichment = tmpEnrichment;
									primaryScores = tmpScores;
								}//if not better than primary but better than secondary, then replace secondary
								else if(tmpEnrichment>secondaryEnrichment){
									secondaryEnrichment = tmpEnrichment;
									secondaryScores = tmpScores;
								}
							}							
							//For best peak
							if(tmpEnrichment>bestEnrichment){
								bestEnrichment = tmpEnrichment;
								bestWindow = geneEnd;
								bestWindowEnd = is3p? -intoGene : intoGene;
								insideGene = true;
								bestScores = tmpScores;
							}
						}
						intoGene += STEP;
					}
					
					//EXTENSION
					int extend = STEP;
					while(extend < maxExtension) {
						Annotation end = null;
						long chrEnd = 0;
						for(String s:libDataModels.keySet())
							chrEnd = libDataModels.get(s).getRefSequenceLength(chr);
						
						if((is3p && gene.getOrientation().equals(Strand.NEGATIVE))||(!is3p && gene.getOrientation().equals(Strand.POSITIVE))){
							if((gene.getStart() - extend)<0){
								logger.info(gene.getStart()+" - "+extend +" is less than 0");
								break;
							}
							end = new BasicAnnotation(gene.getChr(), gene.getStart() - extend, gene.getStart() - (extend-window),gene.getOrientation());
						}
						else{
							if((gene.getEnd() + extend)> chrEnd){
								logger.info(gene.getEnd()+" + "+extend +" is more than "+chrEnd);
								break;
							}
							end = new BasicAnnotation(gene.getChr(), gene.getEnd() + (extend - window), gene.getEnd() + extend,gene.getOrientation());

						}
						// Get an interval tree for all/any exons that overlap with the extended region
						Iterator<Gene> endOverlappersIter = tree.overlappingValueIterator(end.getStart(), end.getEnd());
						//While there is an overlap with a gene and gene is same gene
						boolean overlapperIsSameGene = true;
						while(endOverlappersIter.hasNext() && overlapperIsSameGene){
							
							Gene overlapper = endOverlappersIter.next();
							//If same orientation
							if(overlapper.getOrientation().equals(gene.getOrientation())){
								//compare the end coordiantes of the gene
								if(is3p){
									if(!(overlapper.getOrientedEnd() == gene.getOrientedEnd()))
										if(!BEDFileParser.isOverlapCompatible(gene,overlapper, MIN_OVERLAP))
											overlapperIsSameGene = false;
								}
								else{
									if(!(overlapper.getOrientedStart() == gene.getOrientedStart()) && !overlapper.overlaps(gene))
										if(!BEDFileParser.isOverlapCompatible(gene,overlapper, MIN_OVERLAP))
											overlapperIsSameGene = false;
								}
							}
						}
						if(!overlapperIsSameGene)
							break;
						geneEnd = new Gene(end);
						//geneEnd.setOrientation(gene.getOrientation());
						//No overlap so continue with scoring the region
						Map<String,double[]> tmpScores = getScoresForWindow(libDataModels,geneEnd);
						double tmpEnrichment = getEnrichmentForWindow(libDataModels.keySet(),tmpScores);
						
						//For all peaks
						//If significant
						boolean overlaps = false;
						//Is the window significant
						if(isSignificantWindow(libDataModels.keySet(),tmpScores)){
							//if overlaps any already added peak
							Set<Gene> tempPeaks = new HashSet<Gene>();
							tempPeaks.addAll(thisPeaks);
							for(Gene peak:thisPeaks){
								if(geneEnd.overlaps(peak)){
									overlaps = true;
									//if better then replace
									
									if(tmpEnrichment>peak.getBedScore()){
										geneEnd.setBedScore(tmpEnrichment); //set RPKM as score
										geneEnd.setName(gene.getName());
										tempPeaks.remove(peak);
										sigPeakToScores.remove(peak);
										tempPeaks.add(geneEnd);
										sigPeakToScores.put(geneEnd, new HashMap<String,Double>());
										for(String sample:cols)
											sigPeakToScores.get(geneEnd).put(sample, tmpScores.get(sample)[USE_SCORE]);
									}
								}
							}
							thisPeaks = tempPeaks;
							if(!overlaps){
								geneEnd.setBedScore(tmpEnrichment); //set RPKM as score
								geneEnd.setName(gene.getName());
								//does not overlap. Add
								thisPeaks.add(geneEnd);	
								sigPeakToScores.put(geneEnd, new HashMap<String,Double>());
								for(String sample:cols)
									sigPeakToScores.get(geneEnd).put(sample, tmpScores.get(sample)[USE_SCORE]);
							}
							
							//Better than primary, so this becomes primary and primary becomes secondary
							if(tmpEnrichment>primaryEnrichment){
								secondaryEnrichment = primaryEnrichment;
								secondaryScores = primaryScores;
								primaryEnrichment = tmpEnrichment;
								primaryScores = tmpScores;
							}//if not better than primary but better than secondary, then replace secondary
							else if(tmpEnrichment>secondaryEnrichment){
								secondaryEnrichment = tmpEnrichment;
								secondaryScores = tmpScores;
							}
						}
						if(tmpEnrichment>bestEnrichment){
							bestEnrichment = tmpEnrichment;
							bestWindow = geneEnd;
							bestWindowEnd = extend;
							insideGene = false;
							bestScores = tmpScores;
						}
						extend += (STEP);
					}
					
					String rowName = gene.getName();
					
					for(int i=0;i<cols.size();i++){
						//resultMatrix.set(duplicateNameMap.get(annotation), cols.get(i), scores[i]);
						genesMatrix.set(rowName, cols.get(i), bestScores.get(cols.get(i))[COUNT_SCORE]);
						genespvalueMatrix.set(rowName, cols.get(i), bestScores.get(cols.get(i))[PVAL_SCORE]);
					}
					//SET THE NAME FOR THE WINDOW
					bestWindow.setName(gene.getName());
					geneToBestWindowMap.put(gene.getName(), bestWindow);
					logger.debug("Added peak to geneToWindowMap as in the end "+gene.getName()+" "+geneToBestWindowMap.containsKey(gene.getName()));
					windows.add(bestWindow);
					//full gene matrix
					for(Gene peak:thisPeaks){
						//WRITE TO THE FILE
						bw.write(gene.getName()+"\t"+peak.toUCSC());
						bwBed.write(peak.toBED()+"\n");
						for(String sample:cols){
							bw.write("\t"+sigPeakToScores.get(peak).get(sample));
							double value = fullgeneMatrix.get(gene.getName(), sample)+sigPeakToScores.get(peak).get(sample);
							fullgeneMatrix.set(gene.getName(), sample, value);	
						}
						bw.write("\n");		
					}
					if(thisPeaks.isEmpty()){
						for(String sample:cols){
							fullgeneMatrix.set(gene.getName(), sample, genesMatrix.get(gene.getName(), sample));
						}
					}
					//the per peak matrix
					numPeaksMatrix.set(gene.getName(), NUM_PEAKS, thisPeaks.size());
					
				}
			}
			else{
				logger.info("No data for " + chr);
				IntervalTree<Gene> tree=collapsedGenes.get(chr);
				Iterator<Gene> geneIt=tree.valueIterator();
				logger.info("Processing "+chr);
				while(geneIt.hasNext()) {
					Gene gene = geneIt.next();
					for(int i=0;i<cols.size();i++){
						genesMatrix.set(gene.getName(), cols.get(i), 0.0);
						genespvalueMatrix.set(gene.getName(), cols.get(i), 1.0);
						fullgeneMatrix.set(gene.getName(), cols.get(i), 0.0);
						numPeaksMatrix.set(gene.getName(), NUM_PEAKS, 0);
					}
				}
			}
		}
		
		//TODO: write out the peaks in geneToSignificantWindowsMap
		bw.close();
		bwBed.close();
		
		//Best isoform
		bw =new BufferedWriter(new FileWriter(outputFile+".collapsed.bestIsoform.scores")); 
		genesMatrix.write(bw);
		bw.close();
		if(normalizedOutput)
			writeNormalizedMatrix(outputFile+".collapsed.bestIsoform.scores", genesMatrix);
		if(debugMode){
			bw =new BufferedWriter(new FileWriter(outputFile+".collapsed.bestIsoform.pvalues.matrix")); 
			genespvalueMatrix.write(bw);
			bw.close();	
		////Sum of all significant peaks divided by the number of peaks
			bw =new BufferedWriter(new FileWriter(outputFile+".collapsed.numPeaks")); 
			numPeaksMatrix.write(bw);
			bw.close();		
			
		}
		if(fullGeneScoreFlag){
			//Sum of all the significant peaks
			bw =new BufferedWriter(new FileWriter(outputFile+".fullgene.matrix")); 
			fullgeneMatrix.write(bw);
			bw.close();
			if(normalizedOutput)
				writeNormalizedMatrix(outputFile+".fullgene.matrix", fullgeneMatrix);
		}
		
		if(debugMode){
			BufferedWriter bedBw = new BufferedWriter(new FileWriter(outputFile+".collapsed.bestIsoform.windows.bed"));
			for(Gene gene: windows){
				bedBw.write(gene.toString());
				bedBw.newLine();
			}
			bedBw.close();
		}
	}
	
	/**
	 * Returns true is the window is "Significant" in at least one of the samples.
	 * @param keySet
	 * @param bestScores
	 * @return
	 */
	private static boolean isSignificantWindow(Set<String> allModelNames, Map<String, double[]> modelScores) {
		
		for(String modelName:allModelNames){
			if(modelScores.containsKey(modelName)){
				if(modelScores.get(modelName)[PVAL_SCORE] < PVAL_THRESHOLD)
					return true;
			}
		}
		return false;
	}

	/**
	 * Returns a gene with the specified windows as exons
	 * @param windows
	 * @param name
	 * @param orientation
	 * @return
	 */
	private static Gene makeGene(Collection<Annotation> windows, String name,Strand orientation) {
		
		Collection<Annotation> exons = new ArrayList<Annotation>();
		
		for(Annotation w:windows){
			for(Annotation block:w.getBlocks())
				exons.add(block);
		}
		
		Gene gene = new Gene(exons,name,orientation);		
		return gene;
	}
	
	
	/**
	 * This function calculates the scores for all alignment models from a given window to the start of the gene and adds it to the matrix
	 * @param libDataModels
	 * @param annotationEnd
	 * @param annotation
	 * @param bestWindowEnd : value of this parameter is the end of the best window. If the window is inside the gene of interest, get the annotation upto this distance from the gene end. If the window is outside the gene, get the annotation upto this many bps outside the gene end. 
	 * @return  Vector of size libDataModels with score of each 
	 * @throws IOException
	 */
	private static MatrixWithHeaders getScoreFromWindowToGeneStart(Map<String,AlignmentModel> libDataModels,Gene annotationEnd,Gene annotation,int bestWindowEnd,boolean insideGene,MatrixWithHeaders fullGeneScoreMatrix,String rowName) throws IOException{
		
		Gene subannotation;
		//1.If the annotationEnd is not within the RefSeq Gene
		if(annotationEnd.percentOverlapping(annotation)<1.0){
			//Then, Get subannotation from gene start to gene End, including the alignment from gene end to annotationEnd
			//subannotation = annotation.getExtended3primeIsoform(annotationEnd.absoluteToRelativePosition(annotation.getOrientedEnd())+annotationEnd.getTranscriptLength());
			//EXTEND IS UNSTRANDED
			
			subannotation = annotation.isNegativeStrand() ? annotation.extendAnnotation(bestWindowEnd,0):annotation.extendAnnotation(0,bestWindowEnd);
		}
		//Else get subannotation from gene start to annotationEnd window
		else{
		//	System.out.println("Annotation End is inside the gene:True compared to "+new Boolean(insideGene).toString());
			//subannotation = getSubAnnotationFromEnd(annotation, annotation.getTranscriptLength()-bestWindowEnd,bestWindowEnd-1);
			subannotation=annotation;
		}
		if(oppositeStrand)
			subannotation = reverseStrand(subannotation);
		//Score the subAnnotation
		Map<String,Double> scores = new HashMap<String,Double>();
		for(String st:libDataModels.keySet()){
			//double compareScore = libDataModels[i].scoreGene(annotation)[USE_SCORE];
			if(libDataModels.get(st).containsReference(subannotation.getChr())){	
				scores.put(st,new ScanStatisticScore(libDataModels.get(st), subannotation,false).getScores()[USE_SCORE]);
			}
			else{
				scores.put(st, 0.0);
			}
		}
		
		for(int i=0;i<fullGeneScoreMatrix.columnDimension();i++){
			fullGeneScoreMatrix.set(rowName, fullGeneScoreMatrix.getColoumnName(i), scores.get(fullGeneScoreMatrix.getColoumnName(i)));
		}
		
		return fullGeneScoreMatrix;
	}
	
	private static int getDistance(Gene gene,Gene windows){
			
			if(gene.getOrientedEnd()==windows.getOrientedEnd())
				return 0;
			else{
				if(gene.isNegativeStrand()){
					//Outside gene
					if((gene.getStart()-windows.getStart())>0){
	//					logger.info("Gene is negative and peak outside gene = "+gene.getStart()+"-"+windows.getStart());
						return (gene.getStart()-windows.getStart());
					}
					else{
	//					logger.info("Gene is negative and peak INSIDE gene = "+gene.genomicToTranscriptPosition(windows.getStart())+"-"+gene.genomicToTranscriptPosition(gene.getStart()));
						return (gene.genomicToTranscriptPosition(windows.getStart()) - gene.genomicToTranscriptPosition(gene.getStart()));
					}
				}
				else{
					//Outside gene
					if((windows.getEnd() - gene.getEnd())>0){
	//					logger.info("Gene is positive and peak outside gene = "+windows.getEnd()+"-"+gene.getEnd());
						return (windows.getEnd() - gene.getEnd());
					}
					else{
	//					logger.info("Gene is positive and peak INSIDE gene = "+gene.genomicToTranscriptPosition(windows.getEnd()-1)+"-"+gene.genomicToTranscriptPosition(gene.getEnd()-1));
						return (gene.genomicToTranscriptPosition(windows.getEnd()-1) - gene.genomicToTranscriptPosition(gene.getEnd()-1));
					}
				}
			}
		}
	/**
	 * This function calculates the scores for all alignment models from a given window to the start of the gene and adds it to the matrix
	 * @param libDataModels
	 * @param annotationStart
	 * @param annotation
	 * @param bestWindowEnd : value of this parameter is the end of the best window. If the window is inside the gene of interest, get the annotation upto this distance from the gene end. If the window is outside the gene, get the annotation upto this many bps outside the gene end. 
	 * @return  Vector of size libDataModels with score of each 
	 * @throws IOException
	 */
	private static MatrixWithHeaders getScoreFromWindowToGeneEnd(Map<String,AlignmentModel> libDataModels,Gene annotationStart,Gene annotation,int bestWindowStart,boolean insideGene,MatrixWithHeaders fullGeneScoreMatrix,String rowName) throws IOException{

		Gene subannotation;
		//1.If the annotationStart is not within the RefSeq Gene
		if(annotationStart.percentOverlapping(annotation)<1.0){
			//System.out.println("Annotation End is not inside the gene:False compared to "+new Boolean(insideGene).toString());
			//Then, Get subannotation from gene start to gene End, including the alignment from gene end to annotationStart
			//subannotation = annotation.getExtended3primeIsoform(annotationStart.absoluteToRelativePosition(annotation.getOrientedEnd())+annotationStart.getTranscriptLength());
			subannotation = annotation.isNegativeStrand() ? annotation.extendAnnotation(0,bestWindowStart) : annotation.extendAnnotation(bestWindowStart,0);
		}
		//Else get subannotation from gene start to annotationStart window
		else{
		//	System.out.println("Annotation End is inside the gene:True compared to "+new Boolean(insideGene).toString());
			//subannotation = getSubAnnotationFromStart(annotation, annotation.getTranscriptLength()-bestWindowStart,bestWindowStart-1);
			subannotation=annotation;
		}
	
		//Score the subAnnotation
		Map<String,Double> scores = new HashMap<String,Double>();
		if(oppositeStrand)
			subannotation = reverseStrand(subannotation);
		for(String ss:libDataModels.keySet()){
			
			//double compareScore = libDataModels[i].scoreGene(annotation)[USE_SCORE];
			if(libDataModels.get(ss).containsReference(subannotation.getChr())){	
				scores.put(ss, scores.put(ss,new ScanStatisticScore(libDataModels.get(ss), subannotation,false).getScores()[USE_SCORE]));
				
			}
			else{
				scores.put(ss, 0.0);
			}
		}
		
		for(int i=0;i<fullGeneScoreMatrix.columnDimension();i++){
			fullGeneScoreMatrix.set(rowName, fullGeneScoreMatrix.getColoumnName(i), scores.get(fullGeneScoreMatrix.getColoumnName(i)));
		}
		
		return fullGeneScoreMatrix;
	}

	
	/**
	 * @param args
	 * @throws IOException 
	 * @throws ParseException 
	 */
	public static void main(String[] args) throws ParseException, IOException {
		new ESAT(args);

	}

}
