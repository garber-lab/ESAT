package umms.esat;

import picard.sam.SortSam;

import umms.esat.Window;
import umms.esat.EventCounter;

import java.io.IOException;
import java.text.ParseException;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.lang.annotation.Annotation;
import java.lang.Runtime;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.List;
import java.util.LinkedList;
import java.util.ListIterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.NoSuchElementException;

import org.apache.log4j.Logger;
import org.apache.log4j.LogManager;
import org.apache.log4j.BasicConfigurator;

import umms.core.annotation.Annotation.Strand;
import umms.core.annotation.BEDFileParser;
import umms.core.annotation.BasicAnnotation;
import umms.core.annotation.Gene;
import umms.core.exception.RuntimeIOException;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import net.sf.samtools.*;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import umms.esat.SAMSequenceCountingDict;
import umms.esat.SAMSequenceCountingDictShort;
import umms.esat.SAMSequenceCountingDictFloat;
import umms.core.readers.MappingTableReader;
import umms.core.utils.NexteraPreprocess;
import umms.core.utils.InDropPreprocess;
import umms.core.utils.ExperimentMap;

//import umms.core.utils.ESATUtils;

public class NewESAT {
	
	// This usage message needs to be updated to reflect the actual operation of NewESAT!!!
	static final String usage = "Usage: NewESAT -in <input Bam File> | -alignments <input filelist file>"+
			"\n\t-annotations <reference annotation file [BED file]> | -geneMapping <gene-to-transcript map file>"+
			"\n\t-out <output file basename>"+
			"\n\t**************************************************************"+
			"\n\t\tOPTIONAL arguments"+
			"\n\t**************************************************************"+
			"\n\t-quality <minimum alignment quality [default: no filtering]>"+ 
			"\n\t-task <score3p | score5p> [default: score3p]"+
			"\n\t-unstranded [default: stranded]"+
			"\n\tWindow parameters:"+
			"\n\t\t-wLen <window length [default: 400]>"+
			"\n\t\t-wOlap <window overlap [default: 0]"+
			"\n\t\t-wExt <extension past end of transcript [default: 400]>"+
			"\n\t\t-all [default: disabled]"+
			"\n\tSignificance testing:"+
			"\n\t\t-sigTest <minimum allowable p-value>\n"+
			"\n\tPre-processing alignments from Nextera library reads:"+
			"\n\t\t-nextPrep [default: off]"+
			"\n\tPre-processing alignments from inDrop library reads:"+
			"\n\t\t-inPrep [default: off]"+
			"\n\t\t-umiMin <minimum number of reads per UMI per transcript to be considered valid [default: 10]";
	
	// new comment
	private static HashMap<String,ArrayList<File>> bamFiles;     // key=experiment ID, File[]= list of input files for the experiment
	private static HashMap<String,ArrayList<File>> mmBamFiles;     // key=experiment ID, File[]= list of input files for the experiment (for 'proper' multimap handling)
	private static File outFile;
	private static String annotationFile;
	private static int windowLength;
	private static int windowOverlap;
	private static int windowExtend;
	private static boolean allWindows;  // save all significant windows (default, only single window position with the highest counts 
										// within a set of contiguous overlapping windows.)
	private static String multimap;     // one of "ignore", "normal" or "scale"
	private static boolean qFilter;
	private static int qThresh;         // quality threshold (reads must be GREATER THAN qThresh, if filtering is on 
	private static boolean gMapping;
	private static File gMapFile;       // name of the gene mapping file
	private static String task; 		// 3' or 5' library
	private static float pValThresh;		// minimum allowable p-value for window significance testing
	private static boolean stranded;    // allow for unstranded analysis (defaults to stranded)
	
	/* single-cell parameters */
	private static boolean nextPreprocess;    // Nextera library reads preprocessing flag
											  // NOTE: barcode and UMI are in the read name, separated by "_".
	private static boolean inPreprocess;    // inDrop library reads preprocessing flag
	  										// NOTE: barcode is encoded in filename, UMIs are in the read name, separated by "_".
	private static int umiMin;			// minimum number of reads per UMI that must be mapped to a transcript to be considered a valid UMI 
	private static int bcMin;			// minimum number of reads that must be observed for a barcode to be considered valid (after PCR duplicate removal) 
	
	static final Logger logger = LogManager.getLogger(NewESAT.class.getName());

	private static HashMap<String,HashMap<String,TranscriptCountInfo>> countsMap;
	private static SAMSequenceCountingDict bamDict;
	private static Hashtable<String, Gene> geneTable;
	
	private static InDropPreprocess inDropData;
	private static ExperimentMap expMap;
	
	public NewESAT(String[] args) throws IOException, ParseException, IllegalArgumentException {
	
		/*
		 * @param for ArgumentMap - size, usage, default task
		 * argMap maps the command line arguments to the respective parameters
		 */
	
		// Configure the logger:
		BasicConfigurator.configure();
		
		/* seems like a useful utility that can probably be stripped down */
		ArgumentMap argMap = CLUtil.getParameters(args,usage,"score3p");  /* no default task for now */
		if (!validateArguments(argMap)) {
			logger.error("Please correct input arguments and retry.");
			throw new IOException();
		}

		/* Only create floating-point counters if multimap=="scale" */ 
		if (multimap.equals("scale")) {
			// multimapped read scaling is not implemented yet:
			bamDict = new SAMSequenceCountingDictFloat();	
		} else {
			bamDict = new SAMSequenceCountingDictShort();
		}

		/* START TIMING */
		long startTime = System.nanoTime();
		
		/* Either use the existing gene-to-transcript mapping table, or load in a genomic annotation file */
		Map<String, Collection<Gene>> annotations;
		if (gMapping) {
			/* If collapsing transcripts down to the gene level, load the gene annotation mapping file */
			geneTable = loadGeneTableFromFile(gMapFile); 
			// Create the annotations map, keyed by chromosome:
			annotations = geneMapToAnnotations(geneTable);  
		} else {
			// load the annotations from the annotation (BED) file:
			annotations =  BEDFileParser.loadDataByChr(new File(annotationFile));	
		}
		
		// This should be done regardless of the type of multimap handling, and the rest of the program needs to be
		// re-written to be more efficient, but for now, make an IntervalTree for each chromosome and strand and
		// add all annotations to it. Add extensions to all of the annotations at this point.
		HashMap<String, HashMap<String, IntervalTree<String>>> occupancyTree = new HashMap<String, HashMap<String, IntervalTree<String>>>();
		if (multimap.equals("proper")) {
			fillOccupancyTree(occupancyTree, annotations, windowExtend, task);
		}
		
		/* collect all read start location counts from the input alignments file(s) */
		mmBamFiles = new HashMap<String,ArrayList<File>>();
		bamDict = countReadStartsFromAlignments(bamDict, bamFiles, qFilter, qThresh, multimap, stranded, occupancyTree, mmBamFiles); 
	
		// If handling multimapped reads "properly", call the function again with the multimapped temp files:
		if (multimap.equals("proper")) {
			bamDict = countReadStartsFromAlignments(bamDict, mmBamFiles, qFilter, qThresh, "ignore", stranded, occupancyTree, mmBamFiles);
			// add any files in the mmBamFiles list to the list of bamFiles:
			for (String exp:mmBamFiles.keySet()) {
				Iterator<File> fIter = mmBamFiles.get(exp).iterator();
				while (fIter.hasNext()) {
					File f = fIter.next();
					bamFiles.get(exp).add(f);
				}
			}
			// From here on, ignore multimappers:
			multimap = "ignore";
		}
		
		/*****************************************************************************************************
		 * BEGIN Single-cell data preprocessing 
		 ******************************************************************************************************/
		if (nextPreprocess) {
			NexteraPreprocess nextData = new NexteraPreprocess(bamFiles, annotations, qFilter, qThresh, multimap, windowExtend, stranded, task, umiMin);
			bamFiles = nextData.getPreprocessedFiles();
		} else if (inPreprocess) {
			//InDropPreprocess inDropData = new InDropPreprocess(bamFiles, annotations, qFilter, qThresh, multimap, windowExtend, stranded, task, umiMin);
			/* New version of InDropPreprocess: 
			 * Assumes umiMin=1 and that the barcode and UMI are concatenated with the readID as <readID>:<bc1>:<bc2>:<umi>
			 */
			inDropData = new InDropPreprocess(bamFiles, annotations, qFilter, qThresh, multimap, windowExtend, stranded, task);
			bamFiles = inDropData.getPreprocessedFiles();
			// Fill in barcode counts from preprocessed files, if necessary:
			int rCount = inDropData.fillBarcodeCounts();
			// Remove low-count barcodes, if -bcMin value is given:
			if (bcMin>0) {
				HashMap<String, Integer> bcStats = inDropData.filterLowcountBarcodes(bcMin);
				logger.info((bcStats.get("startCount")-bcStats.get("endCount"))+" low-count barcodes removed. "+
						bcStats.get("endCount")+" remaining.");
			}
		}	
		
		/*****************************************************************************************************
		 * END Single-cell data preprocessing 
		 ******************************************************************************************************/

		/* create the experiment map to be used by makeCountingIntervalTree(), fillExperimentWindowCounter() and writeExperimentCounter(): */
		if (inPreprocess) {
			expMap = new ExperimentMap(bamFiles, inDropData);
		} else {
			expMap = new ExperimentMap(bamFiles);
		}
		
		/* Count all reads beginning within the exons of each of the transcripts in the annotationFile */
		countsMap = bamDict.countWindowedTranscriptReadStarts(annotations, windowLength, windowOverlap, windowExtend, task, pValThresh, allWindows);
		
		/* Make an intervalTree containing only Windows with non-zero counts across ALL experiments */
		//HashMap<String, HashMap<String, IntervalTree<EventCounter>>> windowTree = makeCountingIntervalTree(countsMap, bamFiles.keySet().size());
		HashMap<String, HashMap<String, IntervalTree<EventCounter>>> windowTree = makeCountingIntervalTree(countsMap, expMap.getNexp());
		
		/* re-process the alignments files to count all reads that start within intervals in the windowTree (i.e., within windows in cleanCountsMap) */
		//fillExperimentWindowCounter(windowTree, bamFiles, qFilter, qThresh, multimap, stranded);
		fillExperimentWindowCounter(windowTree, expMap, qFilter, qThresh, multimap, stranded);
		
		/* write the output file */
		//writeExperimentCountsFile(windowTree, bamFiles, outFile);
		writeExperimentCountsFile(windowTree, expMap, outFile);

		/* STOP AND REPORT TIMING */
		long stopTime = System.nanoTime();
		logger.info("Total processing time: "+(stopTime-startTime)/1e9+" sec\n");
	}
									
	public static void main(String[] args) throws ParseException, IOException {
		new NewESAT(args);
	}
	

	private static boolean validateArguments(ArgumentMap argMap) throws IOException {
		/* Validates the input arguments to ensure that all parameters are consistent and
		 * fills in provided and default values for all required parameters.
		 * 
		 * In the case of single/multiple input alignment files, either the -in <inputFile> flag or
		 * -alignments <inputFileListFile> can be used. Multiple files (-alignments) overrides -in. 
		 * The bamFiles argument (a String array) will contain either a single file name (-in), or the
		 * full set of all alignment files (-alignments). NOTE: In the case of multiple alignments files
		 * ONLY the header of the first file is used to determine the alignment segments.ls
		 * 
		 * 
		 * @param	args	an ArgumentMap containing the input parameters from the command line.
		 * @return	success	a boolean flag indicating success/failure of the parameter set validity.
		 */
		
		/* Windowed read count test parameters */
		windowLength = argMap.isPresent("wLen")? argMap.getInteger("wLen") : 50;
		if (windowLength<1) {
			logger.error("Illegal value for wLen: "+windowLength+" (window length must be >= 1.");
			throw new IllegalArgumentException();
		}
		windowOverlap = argMap.isPresent("wOlap")? argMap.getInteger("wOlap") : 0;   // only allow overlap if window significance is being tested
		if (windowOverlap<0) {
			logger.error("Illegal value for wOlap: "+windowOverlap+" (window overlap must be >= 0.");
			throw new IllegalArgumentException();
		}
		windowExtend = argMap.isPresent("wExt")? argMap.getInteger("wExt") : 0;
		if (windowExtend<0) {
			logger.error("Illegal value for wExt: "+windowExtend+" (extension must be >= 0.");
			throw new IllegalArgumentException();
		}
		allWindows = argMap.isPresent("all");
		
		/* Multimapping parameters */
		if (!argMap.isPresent("multimap")) {
			multimap = "normal";     // default to treat multimapped reads as normal reads
		} else {
			multimap = argMap.get("multimap");
			if (!(multimap.equals("ignore") || multimap.equals("normal") || multimap.equals("scale") || multimap.equals("proper"))) {
				logger.error("-multimap flag must be one of ignore, normal, or scale (is set to "+multimap+")");
				throw new IOException();
			}
		}
		
		/* 3' or 5' library (defaults to 3') */
		task = argMap.getTask();  
		
		/* Quality filtering */
		if (!argMap.isPresent("quality")) {
			qFilter = false;
		} else {
			qFilter = true;
			qThresh = argMap.getInteger("quality");   // quality must be GREATER THAN qThresh for read to be processed
			if (qThresh<0) {
				logger.error("Illegal value for quality: "+qThresh+" (quality threshold must be >= 0.");
				throw new IllegalArgumentException();
			}
		}
		
		/* Stranded or unstranded alignments */
		stranded = argMap.isPresent("unstranded")? false : true;
		
		/* single-cell pre-processing? */
		nextPreprocess = argMap.isPresent("nextPrep") ? true : false;
		inPreprocess = argMap.isPresent("inPrep") ? true : false;
		umiMin = argMap.isPresent("umiMin") ? argMap.getInteger("umiMin") : 0;
		bcMin = argMap.isPresent("bcMin") ? argMap.getInteger("bcMin") : 0;
		
		// Allow multiple inputs 
		if (argMap.isPresent("alignments")){
			// The input file file contains a listing of all input files, with an experiment identifier for
			// each file as: <experimentID>\t<input BAM file>
			// fill the File array with the list of input files
			String inputFileFile = argMap.get("alignments");
			bamFiles = loadBamFileList(inputFileFile);    
		}
		if (argMap.hasInputFile()) {
			bamFiles = new HashMap<String, ArrayList<File>>();
			// Set the default experiment name to "Exp1"
			bamFiles.put("Exp1",new ArrayList<File>());
			bamFiles.get("Exp1").add(new File(argMap.getInput()));
		}
		
		// output file
		outFile = new File(argMap.getOutput());

		// Must have either an annotations file or a gene-to-transcript mapping file (gene-mapping takes precedence over annotation file):
		if (argMap.isPresent("geneMapping"))
		{
			//  allow all transcripts for a gene to be read from an input file.
			// The input file is assumed to be a table in the format provided by the UCSC website (genomes.ucsc.edu)
			// with the following features selected:
			// Clade: Mammal
			// genome: human
			// assembly: <appropriate assembly>
			// table: refGene
			// region: genome
			// output format: all fields from selected table
			//
			// The table MUST have the following columns:
			//   name: the transcript RefSeq ID
			//   chrom: chromosome
			//   strand: strand, + or -
			//   txStart, txEnd: transcript start and end location (0-based)
			//   exonStarts, exonEnds: starting and ending location of each exon (paired, 0-based)
			//   name2: gene symbol
			gMapping = true;
			gMapFile = new File(argMap.get("geneMapping"));   // name of the gene mapping file
		} else if (argMap.isPresent("annotations")) {
			annotationFile = argMap.get("annotations");
			if(!annotationFile.endsWith(".bed") && !annotationFile.endsWith(".BED")){
				logger.error("Please supply an annotation file in the BED format",new RuntimeIOException());
				return false;
			}
			gMapping = false;
		} else {
			logger.error("Either an annotation file or gene-to-transcript mapping file must be provided.");
			return false;
		}
		
		// Significance testing:
		pValThresh = argMap.isPresent("sigTest")? argMap.getFloat("sigTest") : 1;   

		return true;   // default return value if all tests pass
	}

//	/*******************************************************************************************************/
//	/****** Methods pulled out of ESATUtils.java ***********************************************************/
//	/*******************************************************************************************************/
	public static void writeOutputBEDFile(HashMap<String,HashMap<String,LinkedList<Window>>> countsMap, File outFile, boolean collapseGenes) throws IOException {
		// Open the output file:
		FileWriter writer = new FileWriter(outFile);

		/* write the data in bedGraph format */
		// Header line:
		writer.write("track type=bedGraph name=\"BedGraph Format\" description=\"BedGraph format\" visibility=full color=200,100,0 altColor=0,100,200 priority=20\n");
		for (String chr:countsMap.keySet()) {
			for (String gene:countsMap.get(chr).keySet()) {
				ListIterator<Window> wIter = countsMap.get(chr).get(gene).listIterator();
				while (wIter.hasNext()) {
					try {
						Window w = wIter.next();
						// bedGraph format: chr\tstart\tend\tvalue:
						String oStr = chr+"\t"+w.getStart()+"\t"+w.getEnd()+"\t"+w.getCount()+"\n";
						writer.write(oStr);
					} catch (NoSuchElementException e) {
						logger.error("NoSuchElementException for "+gene);
					}
				}	
			}
		}
		writer.flush();
		writer.close();
	}

	public void writeExperimentCountsFile(HashMap<String, HashMap<String, IntervalTree<EventCounter>>> windowTree,
			//HashMap<String,ArrayList<File>> bamfiles, File outFile) throws IOException {
			ExperimentMap eMap, File outFile) throws IOException {
				
		String baseName = outFile.getAbsolutePath();
		File wFile = new File(baseName+".window.txt");  // window-level counts file
		File gFile = new File(baseName+".gene.txt");  // gene-level counts file
		
		// Open the output files:
		FileWriter wWriter = new FileWriter(wFile);
		FileWriter gWriter = new FileWriter(gFile);

		// Header line for window file:
		String wStr = "Symbol\tchr\tstart\tend\tstrand";
		// Header line for gene file:
		String gStr = "Symbol\tchr\tstrand";
//		for (String e:bamfiles.keySet()) {      // before single-cell update
		int nCols = eMap.getNexp();             // after single-cell update
		for (int i=0; i<nCols; i++) {				// after single-cell update
			String e = eMap.getName(i);			// after single-cell update
			wStr+="\t"+e;
			gStr+="\t"+e;
		}
		wWriter.write(wStr+"\n");   // write the window file header  
		gWriter.write(gStr+"\n");   // write the gene file header  
		
		//int nExp = bamfiles.keySet().size();   // number of experiments    // before single-cell update
		int nExp = eMap.getNexp();						// after single-cell update
		
		for (String strand:windowTree.keySet()) {
			for (String chr:windowTree.get(strand).keySet()) {
				// Iterate through the interval tree to extract the counts for each window:
				Iterator<EventCounter> eIter = windowTree.get(strand).get(chr).valueIterator();
				while (eIter.hasNext()) {
					EventCounter e = eIter.next();
					// extract the gene name:
					String[] wLoc = e.getName().split("\t");
					String gName = wLoc[0]; 
					
					if (wLoc.length==1) {
						// if the event counter name only has one field, it is a gene-level counter:
						String oStr = gName+"\t"+chr+"\t"+strand;
						float counts=0;
						for (int i=0; i<nExp; i++) {
							oStr += "\t"+e.getCounts(i);
							counts+=e.getCounts(i);
						}
						// don't bother writing genes with no counts <<< NOTE: Now writing ALL genes, even with 0 counts
						//if (counts>0) {
						gWriter.write(oStr+"\n");	// write to gene-level file
						//}
					} else {
						// otherwise, it is a window-level counter:
						String oStr = e.getName();
						oStr += "\t"+strand;
						for (int i=0; i<nExp; i++) {
							oStr+="\t"+e.getCounts(i);
						}
						wWriter.write(oStr+"\n"); 	// write to window-level file
					}
				}
			}
		}

		// flush and close the writers:
		wWriter.flush();
		wWriter.close();
		gWriter.flush();
		gWriter.close();
	}	
	
	public static void writeOutputESATFile(HashMap<String,HashMap<String,LinkedList<Window>>> countsMap, 
			Map<String, Collection<Gene>> annotations, 
			File outFile) throws IOException {

		// Open the output file:
		FileWriter writer = new FileWriter(outFile);

		/* write the data in ESAT format */
		// Header line:
		writer.write("Symbol\ttranscriptIDs\tcounts\n");

		for (String chr:annotations.keySet()) {
			Iterator<Gene> gIter = annotations.get(chr).iterator();
			while (gIter.hasNext()) {
				Gene thisGene = gIter.next();
				// Gene symbol:
				String symbol = thisGene.getName();
				// Isoforms:
				Collection<Gene> isoforms = thisGene.getIsoforms();
				Iterator<Gene> iIter = isoforms.iterator();
				String iStr = null;         // build the list of isoforms
				while (iIter.hasNext()){
					Gene iso = iIter.next();
					if (iStr==null) {
						iStr = iso.getName();
					} else {
						if (!iso.getName().equals(symbol)) {
							iStr += ","+iso.getName();
						}
					}
				}
				// Window counts:
				if (!(countsMap.containsKey(chr) && countsMap.get(chr).containsKey(symbol))) {
					logger.info("No alignments for "+symbol+" ("+chr+")");
					continue;
				}
				Iterator<Window> wIter = countsMap.get(chr).get(symbol).iterator();
				float totalCounts = 0;
				while (wIter.hasNext()) {
					Window w = wIter.next();
					totalCounts+=w.getCount();
				}
				// write the line to the output file:
				writer.write(symbol+"\t"+iStr+"\t"+totalCounts+"\n");
			}
		}
		writer.flush();
		writer.close();
	}
	
	public static HashMap<String, ArrayList<File>> loadBamFileList(String fileListFile) 
		throws IOException {
		HashMap<String, ArrayList<File>> expBamFiles = new HashMap<String, ArrayList<File>>();
		
		BufferedReader br = new BufferedReader(new FileReader(fileListFile));
		String s;
		while((s = br.readLine())!= null){
			//s = s.replaceAll("\\s+", "\t");
			String[] strSplit = s.split("\t");
			//Check for blank lines or comments (start with "#"):
			if (strSplit.length < 2 || strSplit[0].startsWith("#")) {
				continue;
			} 
			//Map of sample name to alignment File
			String exp = strSplit[0];
			File expFile = new File(strSplit[1]);
			if (!expBamFiles.containsKey(exp)) {
				expBamFiles.put(exp, new ArrayList<File>());
			}
			expBamFiles.get(exp).add(expFile);
		}		
		br.close();
		return expBamFiles;
	}

	private static List<Integer> stringToIntList(String[] vals) {
		List<Integer> outList = new ArrayList<Integer>();
		for (int i=0;i<vals.length;i++) {
			outList.add(Integer.parseInt(vals[i]));
		}
		return outList;
	}
	
	public static Map<String, Collection<Gene>> geneMapToAnnotations(Hashtable<String, Gene>gTable) {
		Map<String, Collection<Gene>> annotations = new TreeMap<String, Collection<Gene>>();
		for (String symbol:gTable.keySet()) {
			String chr = gTable.get(symbol).getChr();   // get the chromosome (used as annotation key)
			if (!annotations.containsKey(chr)) {
				annotations.put(chr, new TreeSet<Gene>());
			}
			//annotations.get(chr).add(gTable.get(symbol));
			Gene g = gTable.get(symbol);
			annotations.get(chr).add(g);
		}
		return annotations;
	}
	
	public void fillOccupancyTree(HashMap<String, HashMap<String, IntervalTree<String>>> oTree,
									Map<String, Collection<Gene>> annotations, int wExt, String task) {
		// fill the occupancy tree with gene/transcript intervals:
		for (String chr:annotations.keySet()) {
			// add an entry for this chromosome:
			oTree.put(chr, new HashMap<String, IntervalTree<String>>());
			// make + and - strand IntervalTrees to avoid testing strand of each gene:
			oTree.get(chr).put("+", new IntervalTree<String>());
			oTree.get(chr).put("-", new IntervalTree<String>());
			for (Gene g:annotations.get(chr)) {
				String strand = g.getStrand().toString();   // strand
				String gName = g.getName();    // gene name
				BasicAnnotation[] eSet = g.getExons();
				for (BasicAnnotation e:eSet) {
					oTree.get(chr).get(strand).put(e.getStart(), e.getEnd(), gName);
				}
				// Add one additional "exon" for the extension. It shouldn't matter if it overlaps another gene,
				// since this will be dealt with when the reads are windowed.
				if (wExt>0) {
					int extStart;  // start of extension
					int extEnd;    // end of extension
					
					if ((strand.equals("+") & task.equals("score5p")) || (strand.equals("-") & task.equals("score3p"))) { 
						extEnd = g.getStart();
						extStart = Math.max(extEnd-wExt,0);
					} else {
						extStart = g.getEnd();
						extEnd = extStart+wExt;
					}
					// Add the extension:
					oTree.get(chr).get(strand).put(extStart, extEnd, gName);
				}
			}
		}
		
		return;
	}
								
	public static Hashtable<String, Gene> loadGeneTableFromFile(File gMapFile) throws IOException {
		
		final MappingTableReader mapFile;
		
		long startGmapTime = System.nanoTime();    // load timer
		// open the gene mapping file
		mapFile = new MappingTableReader(gMapFile);
		// make sure it has the mandatory columns:
		String[] reqFields = {"name","chrom","strand","txStart","txEnd","exonStarts","exonEnds","name2"};
		if (!mapFile.hasMandatoryFields(reqFields)) {
			logger.error("Input mapping file "+gMapFile+" is missing required columns.");
		}
		// build gene-to-transcript map
		Hashtable<String, Gene> geneTable = new Hashtable<String, Gene>();
		boolean notDone = true;
		// NOTE: "name" is the transcriptID, "name2" is the gene symbol
		String[] cOrder = {"name2","chrom","txStart","txEnd","name","strand","exonStarts","exonEnds"};  
		while (notDone) {
			String[] mapData = mapFile.readOrderedFieldsFromLine(cOrder);
			if (mapData.length==0) {
				notDone=false;
			} else {
				// Use the following Gene constructor:
				// Gene(String chr, int start, int end, String name, String orientation, List<Integer> exonsStart, List<Integer> exonsEnd)
				List<Integer> eStarts = stringToIntList(mapData[6].split(","));
				List<Integer> eEnds = stringToIntList(mapData[7].split(","));
				Gene newGene = new Gene(mapData[1],Integer.parseInt(mapData[2]),Integer.parseInt(mapData[3]),mapData[4],mapData[5],eStarts,eEnds); 
				String symbol = mapData[0];
				if (!geneTable.containsKey(mapData[0])) {
					// new gene symbol:
					geneTable.put(symbol, newGene);
				} else {
					// add this transcript as an isoform, but only if this isoform is on the same chromosome as the first,
					// as well as on the same strand:
					Strand gStrand = geneTable.get(symbol).getStrand();
					String gChrom = geneTable.get(symbol).getChr();
					Strand iStrand = newGene.getStrand();
					String iChrom = newGene.getChr();
					if (gStrand.equals(iStrand) && gChrom.equals(iChrom)) {
						geneTable.get(symbol).addIsoform(newGene);
					} else {
						logger.warn("New isoform mismatch for "+symbol+" ("+gChrom+gStrand+") with "+mapData[4]+" ("+iChrom+iStrand+")");
					}
				}
			}
		}
		
		long midGmapTime = System.nanoTime();    // end map loading, start isoform collapsing timer
		// Collapse all multi-isoform genes down to a single exon set:
		List<String> keySet = new ArrayList<String>(geneTable.keySet());    // since this loop modifies the Hashtable on the fly, it is necessary
													// to extract the keySet first, rather then (String symbol:geneTable.keySet())
		for (String symbol:keySet) {
			Gene thisGene = geneTable.get(symbol);
			Collection<Gene> isoforms = thisGene.getIsoforms();
			if (isoforms.size()==1) {
				//logger.info(symbol+" has only 1 isoform");
				// Set the top-level gene name to the the gene symbol:
				geneTable.get(symbol).setName(symbol);
			} else {
				//logger.info(symbol+" has "+isoforms.size()+" isoforms");
				boolean firstIsoform = true;
				Gene mergedGene = null;
				// merge all isoforms
				for (Gene iso:isoforms) {
					if (firstIsoform) {
						mergedGene = iso;
						firstIsoform = false;
					} else {
						mergedGene = mergedGene.takeUnion(iso);
					}
				}
				// Create a replacement gene where the top-level gene has the union of all exons of all 
				// isoforms.
				// NOTE:: doing it this way because there is no simple way to reset the exons of a gene.
				//        It might be much faster if that method was available.
				// The top-level gene name is the gene symbol, and all isoforms are named by their RefSeq ID (NM_, NR_, etc.)
				// If a gene has a single isoform, the top-level gene and its isoform will be identical, except for the name.
				// constructor: Gene(String chr, int start, int end, String name, String orientation, List<Integer> exonsStart, List<Integer> exonsEnd)
				/* make the exon start/end lists */
				BasicAnnotation[] exons = mergedGene.getExons();
				List<Integer> exonStarts = new ArrayList<Integer>();
				List<Integer> exonEnds = new ArrayList<Integer>();
				for (BasicAnnotation exon:exons) {
					exonStarts.add(exon.getStart());
					exonEnds.add(exon.getEnd());
				}
				// create the new gene:
				Gene newGene = new Gene(thisGene.getChr(),thisGene.getStart(),thisGene.getEnd(),symbol,thisGene.getStrand().toString(),exonStarts,exonEnds);
				// add all of the isoforms:
				for (Gene iso:isoforms) {
					newGene.addIsoform(iso);
				}
				// replace the old gene with the new one:
				geneTable.remove(symbol);
				geneTable.put(symbol, newGene);
			}
		}
		long endGmapTime = System.nanoTime();    // end isoform collapsing timer
		logger.info("Loading the gene-to-isoform map took "+(midGmapTime-startGmapTime)/1e9+" sec.");
		logger.info("           Collapsing the genes took "+(endGmapTime-midGmapTime)/1e9+" sec.");

		return geneTable;
	}
	
	public SAMSequenceCountingDict countReadStartsFromAlignments (SAMSequenceCountingDict bamDict, HashMap<String,ArrayList<File>> bamFiles,
																	boolean qFilter, int qThresh, String multimap, boolean stranded, 
																	HashMap<String, HashMap<String, IntervalTree<String>>> occupancyTree, 
																	HashMap<String,ArrayList<File>> mmTempFiles) {
		boolean firstFile = true;      // only read the header from the first alignment file
		int goodQualityCount = 0;
		int badQualityCount = 0;
		int bamFileCount = 0;
		int validReadCount = 0;
		int invalidReadCount = 0;		
		int totalValidReadCount = 0;
		int totalInvalidReadCount = 0;
		SAMRecord r;		
		SAMFileWriterFactory sf = new SAMFileWriterFactory();
		SAMFileHeader tempHeader = new SAMFileHeader();
		SAMFileHeader outHeader;
		// to make compilation errors go away:::
		SAMFileWriter bamWriter = null;
		File multimapTempFile = null; 
		HashMap<String, ArrayList<SAMRecord>> mmMap = null;
		File mmFile = null;
		SAMFileWriter mmWriter = null;
		
		// start file loading timer:
		long startTime = System.nanoTime();
		
		/* open the input alignments file */
		for (String exp:bamFiles.keySet()) {
			
			for (int i=0; i<bamFiles.get(exp).size(); i++){

				long loopStartTime = System.nanoTime();    // loop timer
				int mmCount=0;     // count saved multimapped reads
				validReadCount = 0;   // reset read counters
				invalidReadCount = 0;
						
				// open the next bam file in the list:
				File bamFile = (File) bamFiles.get(exp).get(i);
				logger.info("Processing file: "+bamFile+"...");
				SAMFileReader bamReader = new SAMFileReader(bamFile);   // open as a non-eager reader
				SAMFileHeader bamHeader = bamReader.getFileHeader();
				
				//bamReader.setValidationStringency(ValidationStringency.LENIENT);	
				bamReader.setValidationStringency(ValidationStringency.STRICT);	
				SAMRecordIterator bamIterator = bamReader.iterator();
				
				/* Create a HashMap containing all multimapped reads if -multimap == "proper" */
				/* SHOULD: Create and open a temporary BAM file for multimapped reads if -multimap == "proper" */
				if (multimap.equals("proper")) {
					//mmMap = new HashMap<String, ArrayList<SAMRecord>>();
					/* to use less memory, open a temporary BAM file for multimapped reads */
					try {
						mmFile = File.createTempFile("multimapped_valid_", ".bam");   // !!!should probably include the experiment ID
						// delete after exit:
						mmFile.deleteOnExit();
						mmWriter = sf.makeBAMWriter(bamHeader, false, mmFile);
					} catch (Exception e) {
						e.printStackTrace();
					}		
				}

				if (firstFile) {
					// use the header information in the first bam file to create counts storage
					bamDict.setLogger(logger);
					bamDict.copySequences(bamHeader.getSequenceDictionary());    // copy the sequence map from the original dictionary into the counting dict
					firstFile = false;
				}

				// process each read:
				while (bamIterator.hasNext()) {
					try {
						r = bamIterator.next();
					} catch (SAMFormatException e) {
						// skip SAM Format errors but log a warning:
						logger.warn(e.getMessage());
						continue;
					}
					// process the read:
					if (!r.getReadUnmappedFlag()) {
						// if quality filtering is turned on, skip low-quality reads:
						if (qFilter==true) {
							if (r.getMappingQuality()>qThresh){
								goodQualityCount++;
							} else {
								badQualityCount++;
								continue;
							}
						}
						bamDict.updateCount(r, multimap, stranded);
						// proper handling of multimapped reads
						if (multimap.equals("proper") & SAMSequenceCountingDict.getMultimapCount(r)>1) {
							// To reduce the amount of memory required, write the multimapped reads out to a temp file.
							// Then, before processing, sort the temp file by read ID.
							// Write to temp file:     !!!!!!!!!!
							
							// Writing to temp file replaces this:
							// *******
							//String readName = r.getReadName();
							// update the list reads with this read name:
							//if (!mmMap.containsKey(r.getReadName())) {
							//	mmMap.put(readName, new ArrayList<SAMRecord>());
							//}
							//mmMap.get(readName).add(r);
							// *******
							
							// Write the multimapped read out to the temporary BAM file
							mmWriter.addAlignment(r);
							mmCount+=1;   // update multimap count for this file
						}
						// update the read start count
						validReadCount++;
					} else {
						// Skip unmapped reads, but count them 
						invalidReadCount++;
					}
				}

				// close the bam file reader
				bamReader.close();
				if (multimap.equals("proper") & mmCount>0){    // don't bother if there were no multimapped reads
					//System.out.print("Total unique read IDs: "+mmMap.keySet().size()+"\n");
					System.out.print("Total multimapped reads: "+mmCount+"\n");

					// close the temporary BAM file:
					mmWriter.close();
					
					// process the multimapped reads and get the file where the reads were written:
					//multimapTempFile = processMultimapTempFile(mmMap, occupancyTree, bamHeader, sf);  // original
					multimapTempFile = processMultimapTempFile(mmFile, occupancyTree, bamHeader, sf);  // using temp mm file
					// !!! Should give the option of saving the temp file as a command-line argument:
//					if (!saveMMTempFiles) {
//						multimapTempFile.deleteOnExit();
//					}
					// Add this file to the output HashMap:
					if (!mmTempFiles.containsKey(exp)) {
						// initialize, if necessary:
						mmTempFiles.put(exp, new ArrayList<File>());
					}
					mmTempFiles.get(exp).add(multimapTempFile);
					
					// delete the temp multimapped read file:
					mmFile.delete();
				}
				bamFileCount++;
				
				// time the loop:
				long loopEndTime = System.nanoTime();
				logger.info("Experiment "+exp+" BAM file "+bamFile+" processed in "+(loopEndTime-loopStartTime)/1e9+" sec\n");
				logger.info("  "+validReadCount+" valid reads\n");
				logger.info("  "+invalidReadCount+" invalid reads");

				// accumulate counts:
				totalValidReadCount+=validReadCount;
				totalInvalidReadCount+=invalidReadCount;
			}
		}
		
		long stopTime = System.nanoTime();
		logger.info(bamFileCount+" BAM files processed in "+(stopTime-startTime)/1e9+" sec\n");
		logger.info("  "+totalValidReadCount+" total valid reads\n");
		if (qFilter) {
			//logger.info("     "+goodQualityCount+" reads pass the quality threshold\n");
			logger.info("     "+badQualityCount+" reads fail the quality threshold\n");
		}
		logger.info("  "+totalInvalidReadCount+" total invalid reads");
		
		// Return the updated counts dictionary:
		return bamDict;
	}

	//private File processMultimapTempFile(HashMap<String, ArrayList<SAMRecord>> mmMap,
	private File processMultimapTempFile(File mmFile,
											HashMap<String, HashMap<String, IntervalTree<String>>> occupancyTree,											 
											SAMFileHeader bHeader, SAMFileWriterFactory sf) {
		// Process the multimapped reads in this file, save the valid reads in a new temp file,
		// delete this file and return the name of the new file. The new file should have the same header
		// as the original input BAM file.
		File newTempFile = null;
		File mmFileSorted = null;
		SAMFileWriter newWriter = null;
		SAMRecord r;
		int uAligns = 0;
		int mmCount = 0;
		
		SortSam samSorter = new SortSam();
		
		System.out.print("Processing multimappers...");
		// First, sort the input BAM file by read ID:
		try {
			mmFileSorted = File.createTempFile("multimapped_valid_", ".bam");   // !!!should probably include the experiment ID
		} catch (Exception e) {
			e.printStackTrace();
		}		
		// create the command line arg list:
		String[] clOptions = new String[3];
		clOptions[0] = "I="+mmFile.toString();
		clOptions[1] = "SO=queryname";
		clOptions[2] = "O="+mmFileSorted.toString();
		// call the samSorter program:
		int exitCode = samSorter.instanceMain(clOptions);
		// ### should probably test the exit code, in case the sorting fails for some reason....
		
		// Open the (sorted) output file:
		SAMFileReader mmReader = new SAMFileReader(mmFileSorted);   // open as a non-eager reader
		SAMRecordIterator mmIterator = mmReader.iterator();
		
		// create the new output file:
		try {
			newTempFile = File.createTempFile("multimapped_valid_", ".bam");   // !!!should probably include the experiment ID
			// delete after exit:
			newTempFile.deleteOnExit();
			newWriter = sf.makeBAMWriter(bHeader, false, newTempFile);
		} catch (Exception e) {
			e.printStackTrace();
		}		
		
		// START READ PROCESSING:::::::::
		String lastReadID = null;
		ArrayList<SAMRecord> rArray = new ArrayList<SAMRecord>(); 
		while (mmIterator.hasNext()) {
			try {
				r = mmIterator.next();
			} catch (SAMFormatException e) {
				// skip SAM Format errors but log a warning:
				logger.warn(e.getMessage());
				continue;
			}
			
			String thisReadID = r.getReadName();
			mmCount += 1;   // update the multimapped read count
			
			// Process the read array, or add this read to it:
			if (!thisReadID.equals(lastReadID)) {
				if (!rArray.isEmpty()) {
					// Process the ArrayList
					SAMRecord bestRead = findBestReadMapping(rArray, occupancyTree);
					// write the read to the temp file, if it is not null:
					if (bestRead != null) {
						newWriter.addAlignment(bestRead);
						uAligns++;
					}
					// reset rArray, add current read, update lastReadID, continue:
					rArray = new ArrayList<SAMRecord>();
				} 
				// add the current read to the array and continue
				rArray.add(r);
				lastReadID = thisReadID;

			} else {
				// same readID, add to rArray, update lastReadID, continue:
				rArray.add(r);
				lastReadID = thisReadID;
			}
		}

		// Process any reads remaining in rArray:
		if (!rArray.isEmpty()) {
			// Process the ArrayList
			SAMRecord bestRead = findBestReadMapping(rArray, occupancyTree);
			// write the read to the temp file, if it is not null:
			if (bestRead != null) {
				newWriter.addAlignment(bestRead);
				uAligns++;
			}
		}
		// END READ PROCESSING:::::::::::
		
		// close the file:
		newWriter.close();
		//System.out.print("done.\n");
		//System.out.printf(" Number of multimapped reads: %d\n",mmCount);
		System.out.printf(" Number of valid alignments: %d\n",uAligns);
		
		// close the temp multimapped reader:
		mmReader.close();
		// delete the sorted temp multimapped file:
		mmFileSorted.delete();
		
		return newTempFile;
	}

	private SAMRecord findBestReadMapping(ArrayList<SAMRecord> rArray, HashMap<String, HashMap<String, IntervalTree<String>>> occupancyTree) {
		SAMRecord bestRead = null;
		SAMRecord r;
		String chr;
		String strand;
		int readStart; 
		
		// loop over each read ID:
		int bestIdx = -1;    // initialize best index to -1
		int bestCount = 0;   // initialize valid alignment counter
		for (int i=0; i<rArray.size(); i++) {
			r = rArray.get(i);   // get the read
			// check to see if it is mapped to a transcript:
			chr = r.getReferenceName();                     // chromosome
			if (r.getReadNegativeStrandFlag()) {            // strand
				strand = "-";
			} else {
				strand = "+";
			}
			readStart = r.getAlignmentStart();              // start location

			// Test read start against occupancyTree:
			if (occupancyTree.containsKey(chr)) {
				//int oCount = occupancyTree.get(chr).get(strand).numOverlappers(readStart, readStart+1);
				// ignore strandedness of multimapped reads (i.e., test if the read falls within any transcript region, regardless of strand):
				int oCount = occupancyTree.get(chr).get("+").numOverlappers(readStart, readStart+1) +
						occupancyTree.get(chr).get("-").numOverlappers(readStart, readStart+1);
				if (oCount>0) {
					// if this overlaps any intervals in the tree, save this index and increment the count:
					bestIdx = i;
					bestCount++;
				}
			}
		}
		// if there were any reads that map to any gene, but not more than one, write the indexed read to the output file:
		if (bestCount==1) {
			// Change the multimap count to 1:
			r = rArray.get(bestIdx);
			// Finally, make sure that this read maps to this location on the correct strand (this was ignored before):
			chr = r.getReferenceName();
			if (r.getReadNegativeStrandFlag()) {            // strand
				strand = "-";
			} else {
				strand = "+";
			}
			readStart = r.getAlignmentStart();              // start location
			int oCount = occupancyTree.get(chr).get(strand).numOverlappers(readStart, readStart+1);
			if (oCount>0) {
				r.setAttribute("NH", 1);
				bestRead = r;
			}
		}
		
		// Return the best read, if there was one. If there was no best read, bestRead will be null.
		return bestRead;
	}
	
	public HashMap<String, HashMap<String, IntervalTree<EventCounter>>> makeCountingIntervalTree(HashMap<String,HashMap<String,TranscriptCountInfo>> countsMap, int nExp) {
		// Builds a stranded HashMap of IntervalTrees, one per chromosome
		
		HashMap<String, HashMap<String, IntervalTree<EventCounter>>> cleanTree = new HashMap<String, HashMap<String, IntervalTree<EventCounter>>>();
		// create top-level node for the two strands:
		cleanTree.put("+", new HashMap<String, IntervalTree<EventCounter>>());
		cleanTree.put("-", new HashMap<String, IntervalTree<EventCounter>>());
		
		// keep track of how many non-significant windows there are
		int inWindowCount = 0;

		// Iterate over chromosomes:
		for (String chr:countsMap.keySet()) {
			// Iterate over genes/transcripts:
			for (String gene:countsMap.get(chr).keySet()) {
				// First add all intervals of significant windows:
				ListIterator<Window> wIter = countsMap.get(chr).get(gene).getWindows().listIterator();
				String strand = countsMap.get(chr).get(gene).getStrand();
				while (wIter.hasNext()) {
					Window w = wIter.next();
					inWindowCount++;
					// add a new Window to the cleanCountsMap:
					//strand = w.getStrand();
					// tab-delimited node name:
					String nName = gene+"\t"+chr+"\t"+w.getStart()+"\t"+w.getEnd();
					//String nName = gene+"."+listIdx+" "+chr+":"+w.getStart()+"-"+w.getEnd()+" ("+strand+")";
					if (w.getStart()>=w.getEnd()) {
						logger.warn("start>end for "+gene);
					}

					EventCounter e = new EventCounter(nName, nExp, w);
					// *** TEST: initialize the event counter with the total alignments
					e.setSumCounts(w.getCount());
					if (!cleanTree.get(strand).containsKey(chr)) {
						cleanTree.get(strand).put(chr, new IntervalTree<EventCounter>());
					}
					cleanTree.get(strand).get(chr).put(w.getStart(),w.getEnd(), e);   // add the node to the tree
				}
				// next, add an event counter for intervals of the full gene/transcript to allow accumulation of gene-level counts:
				IntervalTree<String> eTree = countsMap.get(chr).get(gene).getITree();
				Window gWindow = new Window(strand, chr, eTree, gene, nExp);
				EventCounter e = new EventCounter(gene, nExp, gWindow);
				if (!cleanTree.get(strand).containsKey(chr)) {
					cleanTree.get(strand).put(chr, new IntervalTree<EventCounter>());
				}
				cleanTree.get(strand).get(chr).put(gWindow.getStart(), gWindow.getEnd(), e);
			}
		}
		
		logger.info("Total window count: "+inWindowCount);
		
		return cleanTree;
	}
	
	public void fillExperimentWindowCounter(HashMap<String, HashMap<String, IntervalTree<EventCounter>>> windowTree, 
											//HashMap<String,ArrayList<File>> bamFiles,
											ExperimentMap eMap,
											boolean qFilter,
											int qThresh,
											String multimap,
											boolean stranded) {
		
		SAMRecord r;		// alignment
		String rStrand;		// alignment strand
		String rName;		// alignment name (chromosome)
		int rStart;			// alignment start location
		
		// Get the list of experiment names:
		//Object[] expList = bamFiles.keySet().toArray();   // before single-cell update
		Object [] expList = eMap.getBamFiles().keySet().toArray();    // after single cell update
		
		// Iterate over each experiment:
		for (int eIdx=0; eIdx<expList.length; eIdx++) {
			
			Object exp = expList[eIdx];
			
			// Iterate over the files in each experiment:
			for (int i=0; i<bamFiles.get(exp).size(); i++){

				//long loopStartTime = System.nanoTime();    // loop timer
				
				// open the next bam file in the list:
				File bamFile = (File) bamFiles.get(exp).get(i);
				SAMFileReader bamReader = new SAMFileReader(bamFile);   // open as a non-eager reader
				bamReader.setValidationStringency(ValidationStringency.STRICT);	
				SAMRecordIterator bamIterator = bamReader.iterator();
				logger.info("Processing file: "+bamFile+"...");

				// process each read:
				while (bamIterator.hasNext()) {
					try {
						r = bamIterator.next();
					} catch (SAMFormatException e) {
						// skip SAM Format errors but log a warning:
						logger.warn(e.getMessage());
						continue;
					}
					// process the read:
					if (!r.getReadUnmappedFlag()) {
						// if quality filtering is turned on, skip low-quality reads:
						if (qFilter==true) {
							if (!(r.getMappingQuality()>qThresh)) {
								// skip bad reads
								continue;
							}
						}
						// Update the counts in cleanCountsMap if the read start location is contained in
						// an interval in the windowTree.
					   	rName = r.getReferenceName();                // chromosome ID
				    	rStart = (int)(r.getAlignmentStart())-1;   // alignments are 1-based, arrays are 0-based
				    	if (stranded & r.getReadNegativeStrandFlag()) {
				    		rStrand = "-";
				    	} else {
				    		rStrand = "+";
				    	}
				    	String cString = r.getCigarString(); 
				    	// Note: if the CigarString is "*", it indicates that the read is unmapped. It would be better 
				    	//       if SAMRecord had a isMapped() method.
				    	if (cString!="*") {
				    		// Deal with multimapped reads:
				    		float fractCount;
			    			int mmCount = SAMSequenceCountingDict.getMultimapCount(r);
				    		if (multimap.equals("normal") || multimap.equals("proper")) {
				    			fractCount=1;
				    		} else if (multimap.equals("ignore")) {
				    			if (mmCount==1) {
				    				fractCount=1;
				    			} else {
				    				fractCount=0;   // hacky way to skip reads... 
				    			}
				    		} else {
				    			// scaled mulitmapped reads:
				    			fractCount=1f/mmCount;
				    		}

				    		// check if this read start is contained in any intervals in the tree:
				    		// NOTE: the reads are checked against the region rStart-1 to rStart so that reads at the right-most position
				    		//  	in a window are placed in the window and not one base past the end of the window. This was causing a problem
				    		//   	where a significant window was found, but the counts in the window were reported as 0, because all of the 
				    		// 		reads were in the last base of the window.
				    		if (windowTree.get(rStrand).containsKey(rName) && windowTree.get(rStrand).get(rName).numOverlappers(rStart-1, rStart)>0) {
			    			Iterator<IntervalTree.Node<EventCounter>> oIter = windowTree.get(rStrand).get(rName).overlappers(rStart-1,rStart);
				    			if (rStart==13857423) {
				    				logger.info("Read at "+rName+":"+rStart+" ("+rStrand+") has "+
				    						windowTree.get(rStrand).get(rName).numOverlappers(rStart-1, rStart)+" overlapping intervals");
				    			}
				    			while (oIter.hasNext()) {
				    				Node<EventCounter> n = oIter.next();
				    				// This node might contain multiple EventCounters. Update them all:
				    				Collection<EventCounter> cvNode = n.getContainedValues();
				    				for (EventCounter e:cvNode) {
				    					// 	update the count for this interval:
				    					if (!eMap.isSingleCell()) {
				    						int cIdx = eMap.getIndex(exp.toString());
				    						//e.addIntervalCount(rStart, rStart+1, eIdx, fractCount);   // add (possibly) fractional counts if read is contained in an interval
				    						e.addIntervalCount(rStart, rStart+1, cIdx, fractCount);   // add (possibly) fractional counts if read is contained in an interval					    								
				    					} else {
				    						String bc = InDropPreprocess.getBarcodeFromRead(r);
				    						String eName = exp+":"+bc;
				    						int cIdx = eMap.getIndex(eName);
				    						if (cIdx>=0) {
				    							// update the count if this is a valid experiment and barcode:
				    							e.addIntervalCount(rStart, rStart+1, cIdx, fractCount);   // add (possibly) fractional counts if read is contained in an interval
				    						}
				    					}
				    				}
				    			}
				    		}
				    	}
					}
				}
				bamReader.close();
			}
		}
	} 	
}	
