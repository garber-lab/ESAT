package umms.core.utils;

import umms.esat.Window;

import java.io.IOException;
import java.text.ParseException;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.List;
import java.util.LinkedList;
import java.util.ListIterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.NoSuchElementException;

import org.apache.log4j.Logger;

import umms.core.annotation.Annotation.Strand;
import umms.core.annotation.BasicAnnotation;
import umms.core.annotation.Gene;
import net.sf.samtools.*;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import umms.esat.SAMSequenceCountingDict;
import umms.esat.SAMSequenceCountingDictShort;
import umms.esat.SAMSequenceCountingDictFloat;
import umms.core.readers.MappingTableReader;

public class ESATUtils {
	
	public static Logger logger = Logger.getLogger(ESATUtils.class.getName());

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
			annotations.get(chr).add(gTable.get(symbol));
		}
		return annotations;
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
																	boolean qFilter, int qThresh) {
		boolean firstFile = true;      // only read the header from the first alignment file
		int goodQualityCount = 0;
		int badQualityCount = 0;
		int bamFileCount = 0;
		int validReadCount = 0;
		int invalidReadCount = 0;		
		int totalValidReadCount = 0;
		int totalInvalidReadCount = 0;
		SAMRecord r;		
		
		// start file loading timer:
		long startTime = System.nanoTime();

		/* open the input alignments file */
		for (String exp:bamFiles.keySet()) {
			
			for (int i=0; i<bamFiles.get(exp).size(); i++){

				long loopStartTime = System.nanoTime();    // loop timer
			
				// open the next bam file in the list:
				File bamFile = (File) bamFiles.get(exp).get(i);
				SAMFileReader bamReader = new SAMFileReader(bamFile);   // open as a non-eager reader
				//bamReader.setValidationStringency(ValidationStringency.LENIENT);	
				bamReader.setValidationStringency(ValidationStringency.STRICT);	
				SAMRecordIterator bamIterator = bamReader.iterator();
				logger.info("Processing file: "+bamFile+"...");

				if (firstFile) {
					// use the header information in the first bam file to create counts storage
					SAMFileHeader bamHeader = bamReader.getFileHeader();
					//bamDict = new SAMSequenceCountingDict_short();
					
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
						bamDict.updateCount(r);
						// update the read start count
						validReadCount++;
					} else {
						// Skip unmapped reads, but count them 
						invalidReadCount++;
					}
				}

				// close the bam file reader
				bamReader.close();
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
		logger.info("  "+totalValidReadCount+" valid reads\n");
		if (qFilter) {
			//logger.info("     "+goodQualityCount+" reads pass the quality threshold\n");
			logger.info("     "+badQualityCount+" reads fail the quality threshold\n");
		}
		logger.info("  "+totalInvalidReadCount+" invalid reads");
		
		// Return the updated counts dictionary:
		return bamDict;
	}
	
}