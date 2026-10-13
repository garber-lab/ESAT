package umms.core.utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.List;
import java.util.Vector;
import java.util.Iterator;
import java.util.LinkedHashMap;

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

	// State 
	private boolean qFilter;
	private int qThreshold;
	private String multimap;
	private int wExt;
	private boolean stranded;
	private String task;
	private boolean filtAT;
	private int umiMin;
	private int filtAtN;
	
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
	public InDropPreprocess(boolean qFilter, int qThresh, String multimap, int wExt, boolean stranded, String task,
			boolean filtAT, int filtAtN , int umiMin)      // umiMin added 
					throws IOException {
		this.qFilter = qFilter;
		this.qThreshold = qThresh;
		this.multimap = multimap;
		this.wExt = wExt;
		this.stranded = stranded;
		this.task = task;
		this.filtAT = filtAT;
		this.umiMin = umiMin;
		this.filtAtN = filtAtN;
	}

	public void  process(HashMap<String, ArrayList<File>> bamFiles, Map<String, Collection<Gene>> annotations, boolean writeUMISupport) throws IOException {	
		SAMRecord r;
		SAMFileWriterFactory sf = new SAMFileWriterFactory();
		File outFile;
		int readsIn = 0;
		int readsOut = 0;
		int totalAlignedUmis = 0;
		int totalDedupedUmis = 0;
		
		/* First, build the exon interval map */
		HashMap<String, HashMap<String, IntervalTree<String>>> eMap = buildExonIntervalMap(annotations, wExt, stranded, task);
		
		/* open the input alignments file */
		for (String exp:bamFiles.keySet()) {
			logger.debug("exp variable: " + exp);
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
							// increment the totalAlignedUmis counter:
							totalAlignedUmis+=1;
							// if so, extract the cell barcode and UMI, and add counts for the overlapping transcript(s)
							if (updateUmiCounts(r,oLaps,umiCount,umiMin)) {
								// only update the totalDedupedUmis counter when writing an alignment:
								totalDedupedUmis+=1;
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
				if(writeUMISupport) {
					/* Keys: strand:chr:gene:barcode:UMI => <count> */
					File outUmiCntFile = new File( exp+".umi.distributions.txt");
					BufferedWriter  outUmiCntFileBW = new BufferedWriter(new FileWriter(outUmiCntFile));
					logger.info("Writing out UMI support data for experiment: " + exp + " to ");
					long startWriteUMISupport = System.nanoTime();  
					for (String strnd : umiCount.keySet()) {
						 HashMap<String, HashMap<String, HashMap<String, HashMap<String, Integer>>>> strandedUmis = umiCount.get(strnd);
						 for( HashMap<String, HashMap<String, HashMap<String, Integer>>>  chrStrandedUmis : strandedUmis.values()) {
							 for (String g : chrStrandedUmis.keySet()) {
								 HashMap<String, HashMap<String, Integer>> geneUmis = chrStrandedUmis.get(g);
								 for (String bc : geneUmis.keySet()) {
									 HashMap<String, Integer> geneBCUmis = geneUmis.get(bc);
									 String key = bc+"\t"+g;
									 for (String umi : geneBCUmis.keySet()) {
										 outUmiCntFileBW.write(key);
										 outUmiCntFileBW.write("\t"+ umi+"\t"+geneBCUmis.get(umi));
										 outUmiCntFileBW.newLine();
									 }
								 }
							 }
						 }
					}
					outUmiCntFileBW.close();
					long timeTaken = System.nanoTime() - startWriteUMISupport;
					logger.info("Finished Writing  UMI support data for experiment: " + exp + " took: " + (timeTaken)/1e9+" sec\n");
					
					// translated Python code for cleaning low-end UMIs
					
					long startCleanLowEndUmis = System.nanoTime();
					
					File cleanLowEndUmisOut = new File(exp + ".expression.txt");
					File cleanLowEndUmisReadsOut = new File(exp + ".reads.txt");
					int nMin = 2;
					int uMin = 1;
					
					cleanLowEndUmis(outUmiCntFile, cleanLowEndUmisOut, cleanLowEndUmisReadsOut, nMin, uMin);
					
					logger.info("Finished cleaning low end UMIs for experiment " + exp + " in " + (System.nanoTime() - startCleanLowEndUmis)/1e9 + "s");
				}
				
			}
		}
		logger.info("Preprocessing complete: Total reads in: "+readsIn+" Total reads out: "+readsOut);
		logger.info("Total aligned UMIs: "+totalAlignedUmis+"  Total de-duped UMIs: "+totalDedupedUmis);
	}
	
	private void cleanLowEndUmis(File inFile, File outFile, File readsFile, int nMin, int uMin) throws IOException {
		BufferedReader fIn = new BufferedReader(new FileReader(inFile));
		
		LinkedHashMap<String, Integer> umiDict = new LinkedHashMap<String, Integer>();
		LinkedHashMap<String, Map<String, Map<String, Integer>>> umiHist = new LinkedHashMap<String, Map<String, Map<String, Integer>>>();
		LinkedHashMap<String, Integer> gDict = new LinkedHashMap<String, Integer>();
		int nLines = 0;
		
		for (String line = fIn.readLine(); line != null; line = fIn.readLine()) {
			nLines++;
			
			// split into barcode and count
			String[] fields = line.split("\t");
			
			// signal bad input line
			if (fields.length < 4) {
				logger.error("Too few items in line: " + line);
			}
			
			// parse the line
			String bc = fields[0];
			String g = fields[1];
			String umi = fields[2];
			int n = Integer.parseInt(fields[3]);
			
			// update the gene dictionary
			gDict.putIfAbsent(g, 0);
			
			// UMI counts histogram
			umiHist.putIfAbsent(bc, new LinkedHashMap<String, Map<String, Integer>>());
			umiHist.get(bc).putIfAbsent(g, new LinkedHashMap<String, Integer>());
			umiHist.get(bc).get(g).putIfAbsent(umi, 0);
			umiHist.get(bc).get(g).replace(umi, umiHist.get(bc).get(g).get(umi) + n);
			
			// keep stats on all UMIs
			umiDict.putIfAbsent(umi, 0);
			umiDict.replace(umi, umiDict.get(umi) + n);
		}
		
		fIn.close();
		
		// get lists of UMIs with one count and UMIs with >nMin counts:
		int meanMin = 2; // minimum number of UMI counts for use in computing UMI mean
		for (String bc : umiHist.keySet()) {
			for (String g : umiHist.get(bc).keySet()) {
				// for average UMI count calculations
				int nzUmis = 0; // number of UMIs with non-zero counts
				int ntUmis = 0; // number of UMIs with >1 counts
				int umiSum = 0; // sum for UMI mean calculation	
				// make lists of singlets and multis
				ArrayList<String> singlets = new ArrayList<String>();
				ArrayList<String> multis = new ArrayList<String>();
				for (String umi : umiHist.get(bc).get(g).keySet()) {
					int count = umiHist.get(bc).get(g).get(umi);
					if (count == 1) {
						singlets.add(umi);
					} else if (count >= nMin) {
						multis.add(umi);
						nzUmis++;
					} else {
						nzUmis++;
					}
				}
				
				// separate true singlets from ones that have a Hamming distance of 1 from one or more of the multis
				ArrayList<String> rmUmis = new ArrayList<String>();
				HashMap<String, Integer> incUmis = new HashMap<String, Integer>();
				testUmiDist(singlets, multis, rmUmis, incUmis);
				
				// remove bad UMIs from the dictionary
				for (String umi : rmUmis) {
					umiHist.get(bc).get(g).remove(umi); // delete error UMI
				}
				// update the counts for UMIs that were uniquely one-off from one of the singletons
				for (String umi : incUmis.keySet()) {
					umiHist.get(bc).get(g).replace(umi, umiHist.get(bc).get(g).get(umi) + incUmis.get(umi));
				}
			}
		}
		
		// For each gene, the expression value is len(umiHist[g].keys()) (i.e., the number of UMIs remaining with at least one copy)
		BufferedWriter fOut = new BufferedWriter(new FileWriter(outFile));
		
		// open the optional read counts file
		boolean writeReads = (readsFile != null) ? true : false;
		BufferedWriter fReads = new BufferedWriter(new FileWriter(readsFile));
		
		// remove any barcodes/cells with fewer than uMin total UMIs (if uMin>0)
		if (uMin > 0) {
			ArrayList<String> bcRm = new ArrayList<String>();
			for (String bc : umiHist.keySet()) { // loop over barcodes
				int bcSum = 0;
				for (String g : umiHist.get(bc).keySet()) { // loop over genes
					bcSum += umiHist.get(bc).get(g).keySet().size(); // sum total UMIs for this cell
				}
				// remove this barcode if the total count is less than uMin
				if (bcSum < uMin) {
					bcRm.add(bc);
				}
			}
			// remove low-count cells
			for (String bc : bcRm) {
				umiHist.remove(bc);
			}
		}
		
		// file header
		String hStr = "gene"; // header
		ArrayList<String> bcList = new ArrayList<String>();
		bcList.addAll(umiHist.keySet());
		for (String bc : bcList) {
			hStr += "\t" + bc;
		}
		fOut.write(hStr + "\n");
		
		// option reads count file header
		if (writeReads) {
			fReads.write(hStr + "\n");
		}
		
		// per-gene count
		for (String g : gDict.keySet()) { // loop over all observed genes
			String oStr = g;
			String rStr = g;
			for (String bc : bcList) { // UMI counts for each cell
				if (umiHist.get(bc).containsKey(g)) {
					oStr += "\t" + umiHist.get(bc).get(g).keySet().size();
					rStr += "\t" + sumValues(umiHist.get(bc).get(g));
				} else {
					oStr += "\t0";
					rStr += "\t0";
				}
			}
			
			fOut.write(oStr + "\n");
			if (writeReads) {
				fReads.write(rStr + "\n");
			}
		}
		
		fOut.close();
		if (writeReads) {
			fReads.close();
		}
	}
	
	/**
	 * find the Hamming distance between two input strings:
	 * @param x
	 * @param y
	 * @return
	 */
	private int hammingDistance(String x, String y) {
		int hd;
		
		if (x.length() != y.length()) {
			hd = x.length();
		} else {
			hd = 0;
			for (int i = 0; i < x.length(); i++) {
				if (x.charAt(i) != y.charAt(i)) {
					hd += 1; // count mismatches
				}
			}
		}
		return hd;
	}
	
	/**
	 * test each element of singList against all elements of multList to see if any
     * singletons are one base off from any of the multis.
	 * @param singList list of UMIs with only one count
	 * @param multList list of UMIs with count >= nMin
	 * @param badUmis output list of UMIs to be removed
	 * @param repUmis UMIs to be incremented, with the increment count
	 */
	private void testUmiDist(List<String> singList, List<String> multList, List<String> badUmis, Map<String, Integer> repUmis) {
		for (String bc0 : singList) {
			ArrayList<String> bc1off = new ArrayList<String>(); // barcodes with a Hamming distance of 1 from the singelton
			for (String bc1 : multList) {
				if (hammingDistance(bc0, bc1) == 1) {
					bc1off.add(bc1);
				}
			}
			// if a barcode is one off from any multi, remove it:
			if (bc1off.size() > 0) {
				badUmis.add(bc0);
			}
			// if bc0 is one off from only ONE multi, increment the count for that multi:
			if (bc1off.size() == 1) {
				if (repUmis.containsKey(bc1off.get(0))) {
					repUmis.replace(bc1off.get(0), repUmis.get(bc1off.get(0)));
				} else {
					repUmis.put(bc1off.get(0), 1);
				}
			}
		}
	}
	
	/**
	 * sums the values of the values in a key:value list:
	 * @param d
	 * @return
	 */
	private int sumValues(Map<String, Integer> d) {
		int x = 0;
		for (Integer value : d.values()) {
			x += value;
		}
		return x;
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

	/* New version with umiMin parameter */
	public boolean updateUmiCounts(SAMRecord r, Vector<String> oLaps, 
			HashMap<String, HashMap<String, HashMap<String, HashMap<String, HashMap<String, Integer>>>>> umiCount,
			int umiMin) {
		boolean writeExemplar=false;
		
		// HashMap keys:
		String chr = r.getReferenceName();
		String strand = r.getReadNegativeStrandFlag() ? "-" : "+";
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
			if (umiN==umiMin) { 
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
	
	class UMICount {
		/* Keys: strand:chr:gene:barcode:UMI => <count> */
		private HashMap<String, HashMap<String, HashMap<String, HashMap<String, HashMap<String, Integer>>>>> data =
				new HashMap<String, HashMap<String, HashMap<String, HashMap<String, HashMap<String, Integer>>>>>();
		
		
	}

}
