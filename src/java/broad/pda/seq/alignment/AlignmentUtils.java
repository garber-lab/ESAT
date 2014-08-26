package broad.pda.seq.alignment;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import javax.swing.text.NumberFormatter;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;
import nextgen.core.annotation.Annotation.Strand;
import nextgen.core.annotation.Gene;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;

import cern.colt.Arrays;

import broad.core.annotation.BED;
import broad.core.annotation.BEDReader;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.datastructures.Alignments;
import broad.pda.seq.segmentation.AlignmentDataModelStats;
import broad.pda.seq.segmentation.ContinuousDataAlignmentModel;
import broad.pda.seq.segmentation.GenericAlignmentDataModel;

import org.apache.commons.io.FilenameUtils;

public class AlignmentUtils {
	static final String USAGE = "Tasks: " +
	"\n\tbreakPairs. Takes a paired map BAM or SAM file and writes out one file corresponding to each left and right pairs -in <Alignment> -outdir <Output directory to write pairs>"+
	"\n\ttoWig. Takes a bam or sam file and outputs a wig file where every bases score is the number of reads overlaping the base. " +
	"\n\t\t-alignment - The alignment file in BAM or SAM format" +
	"\n\t\t-window - The somoothing window to use, the score at a base will be the weighted average for the overlaps of all bases in the window (default 1)" +
	"\n\t\t-minMappingQual - Minimum mapping quality for a read to be counted (default 0)"+
	"\n\t\t-ext - Read extension if desired"+
	"\n\t\t-out - Output file, default output goes to standard output "+
	"\n";
	
	static Logger logger = Logger.getLogger(AlignmentUtils.class.getName());
	static final double WEIGHT_DECREASE_RATE = 1.5;
	
	public static void main (String [] args) throws Exception {
		Globals.setHeadless(true);
		ArgumentMap argMap = CLUtil.getParameters(args,USAGE , "breakPairs");
		if(argMap.getTask().equalsIgnoreCase("breakpairs")) {
			
			String in = argMap.getInput();
			File inFile = new File(in);
			String inName = inFile.getName();
			String suffix = inName.endsWith(".sam") ? ".sam" : ".bam";
			
			String outdir = argMap.getOutputDir();
			String firstPairAln = outdir + "/" + inName.replace(suffix, "").concat(".firstPair").concat(suffix);
			String secondPairAln = outdir + "/" + inName.replace(suffix, "").concat(".secondPair").concat(suffix);
			
			final SAMFileReader inputSam = new SAMFileReader(inFile);
			final SAMFileWriter outputFirstSam = new SAMFileWriterFactory().makeSAMOrBAMWriter(inputSam.getFileHeader(),  true,new File(firstPairAln) );
			final SAMFileWriter outputSecondSam = new SAMFileWriterFactory().makeSAMOrBAMWriter(inputSam.getFileHeader(),  true,new File(secondPairAln) );
			int unmapped = 0;
			int first = 0;
			int second = 0;
			int neither = 0;
		    for (final SAMRecord samRecord : inputSam) {
		    	if(samRecord.getReadUnmappedFlag()) {
		    		unmapped++;
		    		continue; //ignore unmapped reads
		    	}
		    	if(samRecord.getFirstOfPairFlag() && !samRecord.getSecondOfPairFlag()) {
		    		outputFirstSam.addAlignment(samRecord);
		    		first++;
		    	} else if(!samRecord.getFirstOfPairFlag() && samRecord.getSecondOfPairFlag()){
		    		outputSecondSam.addAlignment(samRecord);
		    		second++;
		    	} else {
		    		neither++;
		    	}
		    }

		    outputFirstSam.close();
		    outputSecondSam.close();
		    inputSam.close();
		    System.err.println("DONE: Unammped reads " + unmapped + " first ends " + first + " second ends " + second + " neither " + neither);
		} else if ("TOWIG".equalsIgnoreCase(argMap.getTask())) {
			int minMapQual = argMap.containsKey("minMappingQual") ? argMap.getInteger("minMappingQual") : -1;
			int window     = argMap.containsKey("window") ? argMap.getInteger("window") : 1;
			String alnFile = argMap.getMandatory("alignment");
			int ext        = argMap.containsKey("ext") ? argMap.getInteger("ext") : 0;
			String rgbVector = argMap.containsKey("rgb") ? argMap.getMandatory("rgb") : "80,0,40";
			ContinuousDataAlignmentModel cdam = loadAlignmentData(alnFile, false, minMapQual);
			cdam.setExtensionFactor(ext);
			
			Collection<String> chromosomes = cdam.getChromosomeLengths().keySet();
			double [] count = new double[window];
			double [] weights = computeWeights (window);
			logger.info("Using a smoothing window of " + window + " bases, weights: " + Arrays.toString(weights));
			BufferedWriter bw = argMap.getOutputWriter();
			File alnFileFile = new File(alnFile);
			bw.write("track name="+alnFileFile.getName()+".wig type=wiggle_0 color="+rgbVector );
			bw.newLine();
			for(String chr : chromosomes) {
				logger.debug("Processing chromosome " + chr);
				bw.write("variableStep  chrom="+chr+" span=1");
				bw.newLine();
				int position =  (int) Math.ceil(window/2) - 1;
				double maxPositionCount = 0;
				for(int i = 0; i < window; i++) {
					count[i] =  cdam.count(new Alignments(chr, i, i+1));
					maxPositionCount = Math.max(maxPositionCount, count[i]);
				}

				while (position  < cdam.getChromosomeLength(chr)) {
					if(maxPositionCount > 0) {
						double smoothedLastValue = computeSmoothedCount(count, weights);
						double roundedSmoothedValue = Math.round((10*smoothedLastValue))/10.0;
						if(position >= 0) {
							bw.write(String.valueOf(position+1));//Wigs are 1-based
							bw.write("\t");
							bw.write(String.valueOf(roundedSmoothedValue));
							bw.newLine();
						} else {
							logger.warn("Negative position while working with " + chr + " pos: " + position);
						}
					}
					
					Alignments tmp = new Alignments(chr, position, position+1);
					double nextCount = cdam.count(tmp);
					if(maxPositionCount > 0 || nextCount>0) {
						maxPositionCount = nextCount;
						for(int i = 0; i < window - 1 ; i++) {
							count[i] = count[i+1];
							maxPositionCount = Math.max(maxPositionCount, count[i]);
						}
						count[window-1] = nextCount;
					}
					
					position++;
					if(position % 5000000 ==0) {
						logger.debug("Past position " + position);
					}
				}
				cdam.resetTreeCache();
			}
			bw.close();
		} else {
			System.out.println("Invalid call. Ussage:\n" + USAGE);
		}
	}
	
	
	private static double computeSmoothedCount(double[] count, double[] weights) {
		double value = 0;
		for (int i = 0 ; i < count.length; i++) {
			value += count[i] * weights[i];
		}
		return value;
	}


	private static double[] computeWeights(int window) {
		double [] weights = new double[window];
		if(window == 1) {
			weights[0] =1; 
		} else {
			int midPoint = (int) window/2 - 1;
			//Every step from midpoint will have a weight decrease of WEIGHT_DECREASE_RATE. 
			// If WEIGHT_DECREASE_RATE = 2 (every step from midpoint contributes half the weight), and window = 7 
			// Then we first compute p such that 2*p + 2 * WEIGHT_DECREASE_RATE^2 *p + 2*WEIGHT_DECREASE_RATE^4 + WEIGHT_DECREASE_RATE^8 = 1
			// The we can set the weights. p = 1/(sum(WEIGHT_DECREASE_RATE^i) + WEIGHT_DECREASE_RATE^remainingPairs, where i =0..remainingPairs-1.
			double sum = 0;
			for(int i = 0; i <= midPoint;i++) {
				sum += Math.pow(WEIGHT_DECREASE_RATE,i);
			}
			for(int i = midPoint + 1; i < window; i++) {
				sum += Math.pow(WEIGHT_DECREASE_RATE,window - 1 - i);
			}
			double p = 1/sum;


			for(int i = 0; i <= midPoint;i++) {
				weights[i] = Math.pow(WEIGHT_DECREASE_RATE,i)*p;
			}
			for(int i = midPoint + 1; i < window; i++) {
				weights[i] = Math.pow(WEIGHT_DECREASE_RATE,window - 1 - i)*p;
			}
		}
		return weights;
	}


	public  static ContinuousDataAlignmentModel loadAlignmentData(String alignmentFile) throws IOException {
		return loadAlignmentData(alignmentFile, true);
	}
	
	public  static ContinuousDataAlignmentModel loadAlignmentData(String alignmentFile, boolean loadChromosomeStats) throws IOException {
		return loadAlignmentData(alignmentFile, loadChromosomeStats, 0);
	}
	
	public  static ContinuousDataAlignmentModel loadAlignmentData(String alignmentFile, boolean loadChromosomeStats, int minMappingQuality) throws IOException {

		return loadAlignmentData(alignmentFile, loadChromosomeStats, minMappingQuality, false, false); //DEFAULTs: Not to remove dups nor weigh reads
	}
	
	/**
	 * @author: skadri
	 * @param alignmentFile				Name of alignment file
	 * @param loadChromosomeStats		File with chromosome stats
	 * @param minMappingQuality			Mapping quality threshold
	 * @param removeDuplicateFlags		Set to true if PCR duplicates should be ignored
	 * @return
	 * @throws IOException
	 */
	public  static ContinuousDataAlignmentModel loadAlignmentData(String alignmentFile, boolean loadChromosomeStats, int minMappingQuality, boolean removeDuplicateFlags) throws IOException {
		return loadAlignmentData(alignmentFile, loadChromosomeStats, minMappingQuality, removeDuplicateFlags, false); //DEFAULT: Not to  weigh read counts.
	}
	
	/**
	 * @author: skadri
	 * @param alignmentFile				Name of alignment file
	 * @param loadChromosomeStats		File with chromosome stats
	 * @param removeDuplicateFlags		Set to true if PCR duplicates should be ignored
	 * @param weighReadCounts			Set to true if read counts will be weighed by the NH flag
	 * @return
	 * @throws IOException
	 */
	public  static ContinuousDataAlignmentModel loadAlignmentData(String alignmentFile, boolean loadChromosomeStats, boolean removeDuplicateFlags, boolean weighReadCounts) throws IOException {
		return loadAlignmentData(alignmentFile, loadChromosomeStats, 0, removeDuplicateFlags, weighReadCounts); //DEFAULT: Min mapping quality = 0
	}
	
	/**
	 * @author: skadri
	 * @param alignmentFile				Name of alignment file
	 * @param loadChromosomeStats		File with chromosome stats
	 * @param minMappingQuality			Mapping quality threshold = will be set to 0.0 if weighReadCounts is true
	 * @param removeDuplicateFlags		Set to true if PCR duplicates should be ignored
	 * @param weighReadCounts			Set to true if read counts will be weighed by the NH flag
	 * @return
	 * @throws IOException
	 */
	public  static ContinuousDataAlignmentModel loadAlignmentData(String alignmentFile, boolean loadChromosomeStats, int minMappingQuality, boolean removeDuplicateFlags, boolean weighReadCounts) throws IOException {
		return loadAlignmentData(alignmentFile, loadChromosomeStats, minMappingQuality, removeDuplicateFlags, weighReadCounts, null);
	}
	
	public  static ContinuousDataAlignmentModel loadAlignmentData(String alignmentFile, boolean loadChromosomeStats, int minMappingQuality, 
			  boolean removeDuplicateFlags, boolean weighReadCounts, String strand) throws IOException {
		return loadAlignmentData(alignmentFile, loadChromosomeStats, minMappingQuality, removeDuplicateFlags, weighReadCounts, strand, false);  // DEFAULT: treat paired reads as 1 fragment = false
	}

	public  static ContinuousDataAlignmentModel loadAlignmentData(String alignmentFile, boolean loadChromosomeStats, int minMappingQuality, 
																  boolean removeDuplicateFlags, boolean weighReadCounts, String strand, boolean loadPairsAsFragments) throws IOException {
		logger.info("Processing " + alignmentFile);
		AlignmentDataModelStats dataModelStats = null;
		GenericAlignmentDataModel dataModel = new GenericAlignmentDataModel(alignmentFile, null, false, minMappingQuality, removeDuplicateFlags, weighReadCounts, strand, loadPairsAsFragments);
		
		StoredAlignmentStats sas = loadStoredAlignmentStats(alignmentFile);
		
		if(sas.isSafeToLoad()  ) {
			logger.info("Instantiating model stats from previously computed stats");
			dataModelStats = new AlignmentDataModelStats();
			dataModelStats.setAlignmentDataModel(dataModel);
			dataModelStats.setStoredStats(sas);
		} else {
			logger.info("Model stats were unsafe to load becuase " + sas.getReasonForUnsafe());
			dataModelStats  = new AlignmentDataModelStats(dataModel, null, null, 0, loadChromosomeStats);
			StoredAlignmentStats.store(dataModelStats);
		}
		logger.info("\tcomputed or loaded data model stats");
		ContinuousDataAlignmentModel data = new ContinuousDataAlignmentModel(dataModelStats);
		logger.info("\tInstantiated continuous data alignment model");
		return data;
	}
	
	

	public static void scoreAnnotations(BEDReader annotationReader, ContinuousDataAlignmentModel data) throws IOException {
		Iterator<String> chrIt = annotationReader.getChromosomeIterator();
		while(chrIt.hasNext()) {
			String chr = chrIt.next();
			List<BED> chrAnnotations = annotationReader.getChromosomeBEDs(chr);
			for(BED chrAnnotation: chrAnnotations) {
				throw new UnsupportedOperationException("fixme");
				/*Alignments tmp = new Alignments(chrAnnotation);
				double [] scores = data.scoreSegment(tmp);
				chrAnnotation.addExtraScore(scores[1]); //0 pvalue, 1 enrichment, 2 read number, 3 reads/base
				*/
			}
			//TODO: remove this cast - Moran's problem
			//((AlignmentDataModel) data).resetTreeCache();
			data.resetTreeCache();
		}
		
	}
	
	private static StoredAlignmentStats loadStoredAlignmentStats(String alignmentFileLocation) throws IOException {
		File alignmentFile = new File(alignmentFileLocation);
		String idxExtension = alignmentFileLocation.endsWith(".sam") ? ".sai" : ".bai"; //Only supports bam and sam formats 
		File alignmentIdxFile = new File(alignmentFileLocation + idxExtension);
		File alignmentStoredStatFile = new File(alignmentFileLocation + StoredAlignmentStats.SAVED_MODEL_STATS_EXT);
		
		// Try another format
		if (!alignmentIdxFile.exists()) {
			System.out.println(FilenameUtils.removeExtension(alignmentFileLocation) + idxExtension);
			alignmentIdxFile = new File(FilenameUtils.removeExtension(alignmentFileLocation) + idxExtension);
		}
		
		StoredAlignmentStats sas = new StoredAlignmentStats();
		if(!alignmentFile.exists()) {
			sas.setAsUnsafeToLoad();
			sas.setReasonForUnsafe("Alignment file: " + alignmentFile.getAbsolutePath() + " does not exist");
		}
		
		else if (!alignmentIdxFile.exists()) {
			sas.setAsUnsafeToLoad();
			sas.setReasonForUnsafe("Alignment index file: " + alignmentIdxFile.getAbsolutePath() + " does not exist");
		}
		
		else if(!alignmentStoredStatFile.exists()) {
			sas.setAsUnsafeToLoad();
			sas.setReasonForUnsafe("Did not find stored alignment statistics.");
		} else {
			sas.load(alignmentStoredStatFile);
			/** TODO: Fix this - it doesn't seem to be working for me (Jesse Aug 20 2012)
			if(sas.getLastModifiedStoredTimeStamp() < alignmentStoredStatFile.lastModified()) {
				//sas.setAsUnsafeToLoad();
				sas.setReasonForUnsafe("Stored alignment stats file seems to have been modified since its creation");
			} else if (alignmentFile.lastModified() > sas.getLastModifiedAlignmentTimeStampAtCreation()) {
				sas.setAsUnsafeToLoad();
				sas.setReasonForUnsafe("Algnment file "+alignmentFileLocation+ " seems to have been modified since its creation");
			} else if (alignmentIdxFile.lastModified() > sas.getLastModifiedAlignmentIdxTimeStampAtCreation()) {
				sas.setAsUnsafeToLoad();
				sas.setReasonForUnsafe("Algnment index file "+alignmentIdxFile.getAbsolutePath()+ " seems to have been modified since its creation");
			} else if (alignmentIdxFile.length() != sas.getAlignmentIdxSizeAtCreation()) {
				sas.setAsUnsafeToLoad();
				sas.setReasonForUnsafe("Algnment index file "+alignmentIdxFile.getAbsolutePath()+ " has a different size since alignment statistics were last stored");
							**/ 
			if (alignmentFile.length() != sas.getAlignmentSizeAtCreation()) {
				sas.setAsUnsafeToLoad();
				sas.setReasonForUnsafe("Algnment file "+alignmentFileLocation+ " has a different size since alignment statistics were last stored");
			} else {
				logger.debug("Alignmnent stats stored file looks safe to load.");
				sas.setAsSafeToLoad();
				sas.load(alignmentStoredStatFile);
			}

		}		
		
		return sas;
	}
	public static Map<String, Integer> getSizesFromSAM(String samFileName) {
		final SAMFileReader inputSam = new SAMFileReader(new File(samFileName));
		SAMFileHeader header = inputSam.getFileHeader();
		SAMSequenceDictionary dictionary = header.getSequenceDictionary();
		List<SAMSequenceRecord> references = dictionary.getSequences();
		Map<String, Integer> sizes = new LinkedHashMap<String, Integer>(references.size());
		for(SAMSequenceRecord r : references) {
			sizes.put(r.getSequenceName(), r.getSequenceLength());
		}
		inputSam.close();
		
		return sizes;
	}
	
	public static Gene SAMFormatFullBED(SAMRecord sam){
		String chr=sam.getReferenceName();
    	int start=sam.getAlignmentStart()-1; //TODO: Do we need the -1, it is a left over from the method above.
    	Strand strand=sam.getReadNegativeStrandFlag() ? Strand.NEGATIVE : Strand.POSITIVE;
    	String cigar= sam.getCigarString(); //parse cigar and put the breaks into their component parts, in parallel keep track of the splices
    	Collection<Alignments> aligns=parseCigar(cigar, chr, start);
		Gene gene=new Gene( aligns);
		gene.setName(sam.getReadName());
		gene.setOrientation(strand);
		if(sam.getReadString() != null) {gene.setSequence(sam.getReadString()); }
    	return gene;
	}
	
	private static Collection<Alignments> parseCigar(String cigar, String chr, int start){
		ArrayList[] array=parseCigar(cigar); //type , numbers
		
		ArrayList<Character> type=array[0];
		ArrayList<Integer> numbers=array[1];
		
		Collection<Alignments> rtrn=new ArrayList();
		
		int currentStart=start;
		boolean shouldReturn=true;
		
		for(int i=0; i<type.size(); i++){
			char flag=type.get(i);
			Integer num=numbers.get(i);
			//end is currentStart+num-1;
			if(flag==('M')){Alignments align=new Alignments(chr, currentStart, currentStart+num); currentStart=currentStart+num; rtrn.add(align);}
			else if(flag==('N')){currentStart=currentStart+num;} //currentStart=currentStart+num-1;
			else{shouldReturn=false;}
		}
		if(!shouldReturn){return new ArrayList();}
		return rtrn;
	}
	
	private static ArrayList[] parseCigar(String cigar){
		char[] chars=cigar.toCharArray();
		
		ArrayList type=new ArrayList();
		ArrayList numbers=new ArrayList();
		
		ArrayList ordered=new ArrayList();
		String str="";
		for(int i=0; i<chars.length; i++){
			if(chars[i]=='M' || chars[i]=='N' || chars[i]=='I' || chars[i]=='D' || chars[i]=='S' || chars[i]=='H' || chars[i]=='P'){numbers.add(new Integer(str)); type.add(chars[i]); str=str+chars[i]; ordered.add(str);  str="";}
			else{str=str+chars[i];}
		}
		
		//System.err.println(ordered);
		
		ArrayList[] rtrn={type, numbers};
		return rtrn;
	}
	
	

}
