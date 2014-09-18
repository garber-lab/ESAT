package umms.core;

import umms.esat.Window;

import java.io.IOException;
import java.text.ParseException;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
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

import umms.core.annotation.Gene;
import umms.core.exception.RuntimeIOException;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.annotation.BEDFileParser;
import net.sf.samtools.*;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import umms.esat.SAMSequenceCountingDict;
import umms.core.annotation.Annotation;

public class BamWindowedReadTest {
	
	static final String usage = "Usage: BamReadTest <input bam file> [-scs [-wellBC <n>] [-UMI <m>]] "+
			"\n\t-in <input bam file>: input bam file (sorting and indexing not required)" + 
			"\n\t-out <output file name>"+
			"\n\t-annotations <reference annotation file [BED file]>"+
			"\n\t**************************************************************"+
			"\n\t\tOPTIONAL arguments"+
			"\n\t**************************************************************"+
			"\n\t-scs indicates that input file is UMI- and well barcode-tagged SCS data (optional)"+
			"\n\t-wellBC <n> well barcode length [default=6]"+
			"\n\t-UMI <m> UMI length [default=10]";
	
	private static HashMap<String,ArrayList> bamFiles;     // key=experiment ID, File[]= list of input files for the experiment
	private static File bamFile;
	private static File outFile;
	private static String annotationFile;
	static Logger logger = Logger.getLogger(BamReadTest.class.getName());
	private static SAMRecord r;
	private static int windowLength;
	private static int windowOverlap;
	private static int windowExtend;
	private static String multimap;     // one of "ignore", "normal" or "scale"
	private static SAMSequenceCountingDict bamDict;
	
	public BamWindowedReadTest(String[] args) throws IOException, ParseException {
		
		/*
		 * @param for ArgumentMap - size, usage, default task
		 * argMap maps the command line arguments to the respective parameters
		 */
		/* seems like a useful utility that can probably be stripped down */
		ArgumentMap argMap = CLUtil.getParameters(args,usage,"dummy");  /* no default task for now */
		if (!validateArguments(argMap)) {
			logger.error("Please correct input arguments and retry.");
			throw new IOException();
		}
		HashMap<String,HashMap<String,LinkedList<Window>>> countsMap;
	
		/* START TIMING */
		long startTime = System.nanoTime();
		int totalReadCount = 0;
		int totalInvalidReadCount = 0;
		int bamFileCount = 0;
		boolean firstFile = true;      // only read the header from the first alignment file
		
		/* open the input alignments file */
		for (String exp:bamFiles.keySet()) {
			
			for (int i=0; i<bamFiles.get(exp).size(); i++){

				long loopStartTime = System.nanoTime();    // loop timer
			
				// open the next bam file in the list:
				bamFile = (File) bamFiles.get(exp).get(i);
				SAMFileReader bamReader = new SAMFileReader(bamFile);   // open as a non-eager reader
				bamReader.setValidationStringency(ValidationStringency.LENIENT);	
				SAMRecordIterator bamIterator = bamReader.iterator();

				if (firstFile) {
					// use the header information in the first bam file to create counts storage
					SAMFileHeader bamHeader = bamReader.getFileHeader();
					bamDict = new SAMSequenceCountingDict();
					bamDict.setLogger(logger);
					bamDict.copySequences(bamHeader.getSequenceDictionary());    // copy the sequence map from the original dictionary into the counting dict
					firstFile = false;
				}

				/* read through the input file and increment the count for the start location of each read */
				// NOTE: Baseline test - simply counting up the number of reads starting at each location for a full
				//                       10M counts file takes ~30sec
				int readCount = 0;
				int invalidReadCount = 0;
				while (bamIterator.hasNext()) {
					r = bamIterator.next();
					if (!r.getReadUnmappedFlag()) {
						bamDict.updateCount(r);
						readCount++;
					} else {
						invalidReadCount++;
					}
				}

				// close the bam file reader
				bamReader.close();
				bamFileCount++;
				
				// time the loop:
				long loopEndTime = System.nanoTime();
				logger.info("Experiment "+exp+" BAM file "+bamFile+" processed in "+(loopEndTime-loopStartTime)/1e9+" sec\n");
				logger.info("  "+readCount+" reads\n");
				logger.info("  "+invalidReadCount+" invalid reads");

				// accumulate counts:
				totalReadCount+=readCount;
				totalInvalidReadCount+=invalidReadCount;
			}
		}
		
		long stopTime = System.nanoTime();
		logger.info(bamFileCount+" BAM files processed in "+(stopTime-startTime)/1e9+" sec\n");
		logger.info("  "+totalReadCount+" reads\n");
		logger.info("  "+totalInvalidReadCount+" invalid reads");
		long startTime2 = System.nanoTime();
	
		/* open the input alignments file */
		Map<String,Collection<Gene>> annotations =  BEDFileParser.loadDataByChr(new File(annotationFile));	

		/* Count all reads beginning within the exons of each of the transcripts in the annotationFile */
		countsMap = bamDict.countWindowedTranscriptReadStarts(annotations, windowLength, windowOverlap, windowExtend);

		stopTime = System.nanoTime();
		logger.info("Mapping reads to  "+bamFile+" took "+(stopTime-startTime2)/1e9+" sec\n");
		
		// for testing:
//		for(String chr:annotations.keySet()){
//			for(Gene gene : annotations.get(chr)) {
//				logger.info(gene.getName()+": "+gene.getScore()+" read starts ("+gene.getNormalizedScore()+" normalized)");
//			}	
//		}

		/* STOP AND REPORT TIMING */
		stopTime = System.nanoTime();
		logger.info("Total processing time: "+(stopTime-startTime)/1e9+" sec\n");
		
		
		// Write the output file:
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

	public static void main(String[] args) throws ParseException, IOException {
		new BamWindowedReadTest(args);
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
		windowLength = argMap.isPresent("wLen")? argMap.getInteger("wLen") : 400;
		windowOverlap = argMap.isPresent("wOlap")? argMap.getInteger("wOlap") : 40;
		windowExtend = argMap.isPresent("wExt")? argMap.getInteger("wExt") : 400;
		
		/* Multimapping parameters */
		if (!argMap.isPresent("multimap")) {
			multimap = "ignore";
		} else {
			multimap = argMap.get("multimap");
			if (!(multimap.equals("ignore") || multimap.equals("normal") || multimap.equals("scale"))) {
				logger.error("-multimap flag must be one of ignore, normal, or scale (is set to "+multimap+")");
				throw new IOException();
			}
		}
		
		// Allow multiple inputs 
		if (argMap.isPresent("alignments")){
			// The input file file contains a listing of all input files, with an experiment identifier for
			// each file as: <experimentID>\t<input BAM file>
			// fill the File array with the list of input files
			String inputFileFile = argMap.get("alignments");
			bamFiles = loadBamFileList(inputFileFile);
		}
		if (argMap.hasInputFile()) {
			bamFiles = new HashMap<String, ArrayList>();
			// Set the default experiment name to "Exp1"
			bamFiles.put("Exp1",new ArrayList());
			bamFiles.get("Exp1").add(new File(argMap.getInput()));
		}
		outFile = new File(argMap.getOutput());
		annotationFile = argMap.getMandatory("annotations");
		if(!annotationFile.endsWith(".bed") && !annotationFile.endsWith(".BED")){
			logger.error("Please supply an annotation file in the BED format",new RuntimeIOException());
			return false;
		}
			
		//  
		return true;   // default return value if all tests pass
	}
	
	private static HashMap<String, ArrayList> loadBamFileList(String fileListFile) 
		throws IOException {
		HashMap<String, ArrayList> expBamFiles = new HashMap<String, ArrayList>();
		
		BufferedReader br = new BufferedReader(new FileReader(fileListFile));
		String s;
		while((s = br.readLine())!= null){
			String[] strSplit = s.split("\t");
			//Check for blank lines or comments (start with "#"):
			if (strSplit.length < 2 || strSplit[0].startsWith("#")) {
				continue;
			} 
			//Map of sample name to alignment File
			String exp = strSplit[0];
			File expFile = new File(strSplit[1]);
			if (!expBamFiles.containsKey(exp)) {
				expBamFiles.put(exp, new ArrayList());
			}
			expBamFiles.get(exp).add(expFile);
		}		
		br.close();
		return expBamFiles;
	}
}