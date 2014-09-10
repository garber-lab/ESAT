package umms.core;

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
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

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

public class BamReadTest {
	
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
	
	private static int wellBcLength;
	private static int umiLength;
	private static boolean scsFlag;
	private static File bamFile;
	private static File outFile;
	private static String annotationFile;
	static Logger logger = Logger.getLogger(BamReadTest.class.getName());
	private static SAMRecord r;
	
	public BamReadTest(String[] args) throws IOException, ParseException {
		
		/*
		 * @param for ArgumentMap - size, usage, default task
		 * argMap maps the command line arguments to the respective parameters
		 */
		/* seems like a useful utility that can probably be stripped down */
		ArgumentMap argMap = CLUtil.getParameters(args,usage,"dummy");  /* no default task for now */

		scsFlag = argMap.isPresent("scs");
		wellBcLength = argMap.isPresent("wellBC")? argMap.getInteger("wellBC") : 6;
		umiLength = argMap.isPresent("UMI")? argMap.getInteger("UMI") : 10;
		
		bamFile = new File(argMap.getInput());
		outFile = new File(argMap.getOutput());
		annotationFile = argMap.getMandatory("annotations");
		if(!annotationFile.endsWith(".bed") && !annotationFile.endsWith(".BED")){
			logger.error("Please supply an annotation file in the BED format",new RuntimeIOException());
		}
		
		/* START TIMING */
		long startTime = System.nanoTime();
				
		/* open the input alignments file */
		SAMFileReader bamReader = new SAMFileReader(bamFile, true);   // open as an eager reader
		bamReader.setValidationStringency(ValidationStringency.LENIENT);	
		SAMRecordIterator bamIterator = bamReader.iterator();

		SAMFileHeader bamHeader = bamReader.getFileHeader();
		SAMSequenceCountingDict bamDict = new SAMSequenceCountingDict();
		bamDict.setLogger(logger);
		bamDict.copySequences(bamHeader.getSequenceDictionary());    // copy the sequence map from the original dictionary into the counting dict
		
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
		long stopTime = System.nanoTime();
		logger.info("BAM file "+bamFile+" processed in "+(stopTime-startTime)/1e9+" sec\n");
		logger.info("  "+readCount+" reads\n");
		logger.info("  "+invalidReadCount+" invalid reads");
		long startTime2 = System.nanoTime();

		/* open the input alignments file */
		Map<String,Collection<Gene>> annotations =  BEDFileParser.loadDataByChr(new File(annotationFile));	
		
		/* Count all reads beginning within the exons of each of the transcripts in the annotationFile */
		bamDict.simpleCountTranscriptReadStarts(annotations);

		stopTime = System.nanoTime();
		logger.info("Mapping reads to  "+bamFile+" took "+(stopTime-startTime2)/1e9+" sec\n");
		
		// for testing:
		for(String chr:annotations.keySet()){
			for(Gene gene : annotations.get(chr)) {
				logger.info(gene.getName()+": "+gene.getScore()+" read starts ("+gene.getNormalizedScore()+" normalized)");
			}	
		}

		/* STOP AND REPORT TIMING */
		stopTime = System.nanoTime();
		logger.info("Total processing time: "+(stopTime-startTime)/1e9+" sec\n");

		/* close the files */
		bamReader.close();
	}

	public static void main(String[] args) throws ParseException, IOException {
		new BamReadTest(args);
	}
}