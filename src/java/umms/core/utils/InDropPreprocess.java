package umms.core.utils;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.List;
import java.util.Vector;
import java.util.Iterator;

import broad.core.datastructures.IntervalTree;
import net.sf.samtools.SAMProgramRecord;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMFormatException;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecordIterator;



//import org.apache.commons.math3.util.MultidimensionalCounter.Iterator;
import org.apache.log4j.Logger;

import umms.core.annotation.Gene;
import umms.core.annotation.BasicAnnotation;
import umms.esat.SAMSequenceCountingDict;

public class InDropPreprocess {

	private static final String PROGRAM_VERSION = "0.01";
	static Logger logger = Logger.getLogger(InDropPreprocess.class.getName());
	/* output pre-processed BAM files */
	HashMap<String,ArrayList<File>> bamFiles_prep = new HashMap<String,ArrayList<File>>();
	HashMap<String,Integer> bcCounts = new HashMap<String, Integer>();
	
	public static String getBarcodeFromRead(SAMRecord r) {
		String BC = "";	
		// extract the well barcode and UMI:
		String readName = r.getReadName();
		String[] fields = readName.split(":");
		// New version expects <read name>:<barcode>:<UMI> (3 fields minimum):
		if (fields.length < 3) {
			logger.warn("Improper read name: "+readName);
		} else {
			int fLen = fields.length;
			BC = fields[fLen-2];
		}

		return BC;
	}	
	
	/* this version of InDropPreprocess does not have a umiMin parameter, and saves the first read mapping to 
	 * gene/barcode/umi. The barcode and UMI are assumed to be concatenated with the read ID as <readID>:<bc>:<umi>: 
	 */
	public InDropPreprocess(HashMap<String,ArrayList<File>> bamFiles, 
			Map<String, Collection<Gene>> annotations, 
			boolean qFilter, int qThresh, String multimap, int wExt, boolean stranded, String task,
			boolean filtAT, int filtAtN) throws IOException {
		
		SAMRecord r;
		SAMFileWriterFactory sf = new SAMFileWriterFactory();
		File outFile;
		int readsIn = 0;
		int readsOut = 0;
		
		/* First, build the exon interval map */
		HashMap<String, HashMap<String, IntervalTree<String>>> eMap = buildExonIntervalMap(annotations, wExt, stranded, task);
		
		/* open the input alignments file */
		for (String exp:bamFiles.keySet()) {
			
			/* storage for well barcode/UMI counts (new for each experiment)*/
			/* Keys: strand:chr:gene:barcode:UMI => <count> */
			HashMap<String, HashMap<String, HashMap<String, HashMap<String, HashMap<String, Integer>>>>> umiCount = 
				new HashMap<String, HashMap<String, HashMap<String, HashMap<String, HashMap<String, Integer>>>>>();

			for (int i=0; i<bamFiles.get(exp).size(); i++){

				long loopStartTime = System.nanoTime();    // loop timer

				// open the next BAM file in the list:
				File bamFile = (File) bamFiles.get(exp).get(i);
				logger.info("Processing file: "+bamFile+"...");
				SAMFileReader bamReader = new SAMFileReader(bamFile);   // open as a non-eager reader
				SAMFileHeader bamHeader = bamReader.getFileHeader();    // get the header information
				
				// create the pre-processed output BAM file:
				// The processed file is called <original bam file base name>_inPrep.bam
				String inFile = bamFile.getCanonicalPath();
				int extPos = inFile.length()-4;
				// File outFile = new File(inFile.substring(0,extPos)+"_inPrep"+inFile.substring(extPos));
				
				// Make a temporary file for the PCR de-duplicated reads:
				//logger.info("creating temporary file: "+inFile.substring(0,extPos)+"xxxxxxx.bam");				
				//outFile = File.createTempFile(inFile.substring(0,extPos), ".bam");   
				outFile = File.createTempFile("sc_umiFiltered_", ".bam");   
				logger.info("temporary file created: "+outFile);
				// delete after exit:
				outFile.deleteOnExit();
				
				// add this file to the list of files to be processed by ESAT:
				if (!bamFiles_prep.containsKey(exp)) {
					bamFiles_prep.put(exp, new ArrayList<File>());
				}
				bamFiles_prep.get(exp).add(outFile);

				// copy the header from the input BAM file:
				SAMFileWriter bamWriter = sf.makeBAMWriter(bamHeader, false, outFile);

				//bamReader.setValidationStringency(ValidationStringency.LENIENT);	
				bamReader.setValidationStringency(ValidationStringency.STRICT);	
				SAMRecordIterator bamIterator = bamReader.iterator();

				int readCount = 0;
				int writeCount = 0;
				
				while (bamIterator.hasNext()) {
					try {
						r = bamIterator.next();
						
						// test long A/T stretch reads:
 						if (filtAT && SAMSequenceCountingDict.passATFilt(r, filtAtN)) {
 							continue;
 						}
 						
						readCount+=1;
						int mmCount=SAMSequenceCountingDict.getMultimapCount(r);
						if (mmCount>1 && multimap.equals("ignore")) {
							// skip multimapped reads if "ignore" is selected:
							continue;
						}
						
						// check if read start overlaps any transcript in the annotations:
						Vector<String> oLaps = readStartOverlap(r, eMap); 
						if (!oLaps.isEmpty()) {
							// if so, extract the cell barcode and UMI, and add counts for the overlapping transcript(s)
							if (updateUmiCounts(r,oLaps,umiCount)) {
								/* write this read out as the exemplar read for this cell/transcript/UMI */
								bamWriter.addAlignment(r);
								writeCount+=1;
								/* update the count for this exp:barcode */
								String bc = getBarcodeFromRead(r);
								String expBc = exp+":"+bc;
								if (bcCounts.containsKey(expBc)) {
									bcCounts.put(expBc, bcCounts.get(expBc)+1);
								} else {
									bcCounts.put(expBc, 1);
								}
							}
						}
					} catch (SAMFormatException e) {
						// skip SAM Format errors but log a warning:
						logger.warn(e.getMessage());
						continue;
					}
				}
				
				bamReader.close();
				bamWriter.close();
				
				long loopEndTime = System.nanoTime();    // loop timer
				logger.info("Reads: "+readCount+" writes: "+writeCount);
				logger.info("Preprocessing file: "+bamFile+" took "+(loopEndTime-loopStartTime)/1e9+" sec\n");
				readsIn+=readCount;
				readsOut+=writeCount;
			}
		}
		logger.info("Preprocessing complete: Total reads in: "+readsIn+" Total reads out: "+readsOut);
	}
	
	public HashMap<String,ArrayList<File>> getPreprocessedFiles() {
		return bamFiles_prep;
	}
	
	public int fillBarcodeCounts(boolean filtAT, int filtAtN) {
		SAMRecord r;
		int rCount = 0;
		// Only fill bcCounts if it is empty:
		if (bcCounts.isEmpty()) {
			// If it is empty, fill it with data from PCR-duplicate-cleaned files:
			/* open the input alignments file */
			for (String exp:bamFiles_prep.keySet()) {
				for (int i=0; i<bamFiles_prep.get(exp).size(); i++){
					// open the next BAM file in the list:
					File bamFile = (File) bamFiles_prep.get(exp).get(i);
					logger.info("Processing file: "+bamFile+"...");
					long loopStartTime = System.nanoTime();    // loop timer
					SAMFileReader bamReader = new SAMFileReader(bamFile);   // open as a non-eager reader
					// process all reads:
					bamReader.setValidationStringency(ValidationStringency.STRICT);	
					SAMRecordIterator bamIterator = bamReader.iterator();

					while (bamIterator.hasNext()) {
						try {
							r = bamIterator.next();
							
							// test long A/T stretch reads:
	 						if (filtAT && SAMSequenceCountingDict.passATFilt(r, filtAtN)) {
	 							continue;
	 						}
	 						
							rCount++;
							/* update the count for this exp:barcode */
							String bc = getBarcodeFromRead(r);
							String expBc = exp+":"+bc;
							if (bcCounts.containsKey(expBc)) {
								bcCounts.put(expBc, bcCounts.get(expBc)+1);
							} else {
								bcCounts.put(expBc, 1);
							}
						} catch (SAMFormatException e) {
							// skip SAM Format errors but log a warning:
							logger.warn(e.getMessage());
							continue;
						}
					}
					bamReader.close();
					long loopEndTime = System.nanoTime();    // loop timer
					logger.info("Counting barcodes in file: "+bamFile+" took "+(loopEndTime-loopStartTime)/1e9+" sec\n");
				}
			}
		}
		return rCount;
	}
	
	/* Remove barcodes with observations < bcMin, and return a HashMap with some statistics */
	public HashMap<String, Integer> filterLowcountBarcodes(int bcMin) {
		HashMap<String, Integer> bcStats = new HashMap<String, Integer>();
		bcStats.put("startCount", bcCounts.size());  // number of barcodes before removing low counts
		int rangeMin = Integer.MAX_VALUE;
		int rangeMax = 0;
		List<String> badBc = new ArrayList<String>();
		// Find barcodes with too few reads and add them to a list:
		for (String bc:bcCounts.keySet()) {
			int bcN = bcCounts.get(bc);    // barcode observations
			if (bcN < bcMin) {
				badBc.add(bc);
			} else {
				if (bcN<rangeMin) {
					rangeMin = bcN;           // smallest number of barcode observations > bcMin
				} 	
				if (bcN>rangeMax) {
					rangeMax = bcN;           // largest number of barcode observations
				}
			}
		}
		
		// remove the bad barcodes from the HashMap:
		for (String bc:badBc) {
			bcCounts.remove(bc);
		}
		
		// fill in stats and return:
		bcStats.put("endCount",bcCounts.size());     // number of barcodes after cleaning
		bcStats.put("minCount", rangeMin);
		bcStats.put("maxCount", rangeMax);
		
		return bcStats;
	}

	/* New version without umiMin parameter */
	public boolean updateUmiCounts(SAMRecord r, Vector<String> oLaps, 
			HashMap<String, HashMap<String, HashMap<String, HashMap<String, HashMap<String, Integer>>>>> umiCount) {
		boolean writeExemplar=false;
		
		// HashMap keys:
		String chr = r.getReferenceName();
		String strand = r.getReadNegativeStrandFlag() ? "-" : "+";
		String gName;
		String BC;
		String UMI;
		
		// extract the well barcode and UMI:
		String readName = r.getReadName();
		String[] fields = readName.split(":");
		if (fields.length < 3) {
			logger.warn("Improper read name: "+readName);
			return writeExemplar;
		} else {
			int fLen = fields.length;
			BC = fields[fLen-2];
			UMI = fields[fLen-1];
		}
		
		// update the counts:
		// Strand key:
		if (!umiCount.containsKey(strand)) {
			umiCount.put(strand, new HashMap<String, HashMap<String, HashMap<String, HashMap<String, Integer>>>>());
		}
		// chromosome key:
		if (!umiCount.get(strand).containsKey(chr)) {
			umiCount.get(strand).put(chr, new HashMap<String, HashMap<String, HashMap<String, Integer>>>());
		}
		// gene symbol(s) keys:
		for (String g:oLaps) {
			if (!umiCount.get(strand).get(chr).containsKey(g)) {
				umiCount.get(strand).get(chr).put(g, new HashMap<String, HashMap<String, Integer>>());
			}
			// well barcode keys:
			if (!umiCount.get(strand).get(chr).get(g).containsKey(BC)) {
				umiCount.get(strand).get(chr).get(g).put(BC, new HashMap<String, Integer>());
			}
			if (!umiCount.get(strand).get(chr).get(g).get(BC).containsKey(UMI)) {
				umiCount.get(strand).get(chr).get(g).get(BC).put(UMI, 0);
			}
			// Update the counts for this gene/barcode/UMI
			umiCount.get(strand).get(chr).get(g).get(BC).put(UMI,umiCount.get(strand).get(chr).get(g).get(BC).get(UMI)+1);
			// If any transcript hits the UMI count threshold, set the writeExemplar flag:
			int umiN = umiCount.get(strand).get(chr).get(g).get(BC).get(UMI);   // TEMPORARY
			if (umiCount.get(strand).get(chr).get(g).get(BC).get(UMI)==1) {
				// only write the first read with a BC:UMI mapped to this location:
				writeExemplar=true;
			}
		}
		
		return writeExemplar;
	}	
	
	public static Vector<String> readStartOverlap(SAMRecord r, HashMap<String, HashMap<String, IntervalTree<String>>> eMap) {
		String chr;		// alignment chromosome
		String strand;	// alignment strand   /* TODO: add unstranded */
		int aStart;		// alignment start
		Vector<String> oLaps;	// true=has overlap, false=no overlap
		
		oLaps = new Vector<String>();	// initialize to "no overlap"
		
		chr = r.getReferenceName();
		strand = r.getReadNegativeStrandFlag() ? "-" : "+";   // read strand
		// negative-strand alignments 'start' at the read end:
		if (strand.matches("\\+")) {
			aStart = r.getAlignmentStart();
		} else {
			aStart = r.getAlignmentEnd();
		}
				
		if (eMap.containsKey(chr)) {
			if (eMap.get(chr).containsKey(strand)) {
				IntervalTree<String> gSet = eMap.get(chr).get(strand);
				if (gSet.numOverlappers(aStart, aStart+1)>0) {
					Iterator<IntervalTree.Node<String>> oIter = gSet.overlappers(aStart, aStart+1);
					while (oIter.hasNext()) {
						oLaps.add(oIter.next().getValue());   // add the gene symbol to the list
					}
				}
			}
		}
		return oLaps;
	}
	
	HashMap<String, HashMap<String, IntervalTree<String>>> buildExonIntervalMap(Map<String, Collection<Gene>> annotations,
										int wExt, boolean stranded, String task) {
		// Builds a set of IntervalTrees, one per chromosome, containing all exons, plus all extensions:
		HashMap<String, HashMap<String, IntervalTree<String>>> eTree = new HashMap<String, HashMap<String, IntervalTree<String>>>();
		
		/* first pass: create interval trees with exon coordinates from annotations */
		for (String chr:annotations.keySet()) {
			for (Gene g:annotations.get(chr)) {
				if (!eTree.containsKey(chr)) {
					eTree.put(chr, new HashMap<String, IntervalTree<String>>());
				}
				String gStrand = g.getStrand().toString();
				if (!eTree.get(chr).containsKey(gStrand)) {
					eTree.get(chr).put(gStrand, new IntervalTree<String>());
				}
				
				Iterator<BasicAnnotation> eIter = (Iterator<BasicAnnotation>) g.getBlocks().iterator();  // exon iterator
				String gName = g.getName();  // gene symbol
				while (eIter.hasNext()) {
					BasicAnnotation a = eIter.next();
					int aStart = a.getStart();
					int aEnd = a.getEnd();
					// create a tree for this chromosome, if necessary:
					// Add the interval with the gene symbol as the name:
					eTree.get(chr).get(gStrand).put(aStart,aEnd,gName);
				}
			}
		}
		
		/* second pass: add extensions */
		for (String chr:annotations.keySet()) {
			for (Gene g:annotations.get(chr)) {
				String gStrand = g.getStrand().toString();
				String gName = g.getName();
				int gStart = g.getStart();   // lowest gene coordinate
				int gEnd = g.getEnd();       // highest gene coordinate
				int extStart;
				int extEnd;
				// Compute extended transcript coordinates:
				if ((gStrand.equals("+") & task.equals("score3p")) | (gStrand.equals("-") & task.equals("score5p"))) {
					// extend past the 'right' end as far as possible, up to wExt or collision with the neighboring gene:
					extStart = gEnd;
					extEnd = gEnd+wExt;
					if (eTree.get(chr).get(gStrand).numOverlappers(extStart, extEnd)>0) {
						int nOlap = eTree.get(chr).get(gStrand).numOverlappers(extStart, extEnd);
						IntervalTree.Node<String> minNode = eTree.get(chr).get(gStrand).minOverlapper(extStart, extEnd);
						// Truncate RIGHT end of extension to avoid collision:
						Iterator<IntervalTree.Node<String>> oIter = eTree.get(chr).get(gStrand).overlappers(extStart, extEnd);
						while (oIter.hasNext()) {
							IntervalTree.Node<String> n = oIter.next();
							if (n.getStart()<extEnd) {
								extEnd = Math.max(n.getStart(),extStart);
							}
						}
					}
				} else {
					// extend past the 'left' end as far as possible, up to wExt or collision with the neighboring gene:
					extEnd = gStart;
					extStart = Math.max(0,extEnd-wExt);
					if (eTree.get(chr).get(gStrand).numOverlappers(extStart, extEnd)>0) {
						// Truncate LEFT end of extension to avoid collision:
						Iterator<IntervalTree.Node<String>> oIter = eTree.get(chr).get(gStrand).overlappers(extStart, extEnd);
						while (oIter.hasNext()) {
							IntervalTree.Node<String> n = oIter.next();
							if (n.getEnd()>extStart) {
								extStart = Math.min(n.getEnd(),extEnd);
							}
						}
					}
				}
				/* If the interval is non-zero, add it to the tree */
				if (extEnd>extStart) {
					eTree.get(chr).get(gStrand).put(extStart,  extEnd, gName);
				}
			}
		}
		
		return eTree;
	}
	
	// get the bcCounts map
	public HashMap<String,Integer> getBcCounts() {
		return bcCounts;
	}
	public static void main() {
		
	}
}