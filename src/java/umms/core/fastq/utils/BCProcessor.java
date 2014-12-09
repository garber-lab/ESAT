package umms.core.fastq.utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.regex.Pattern;

import org.apache.log4j.Logger;

import picard.cmdline.CommandLineProgram;
import picard.cmdline.Option;
import picard.cmdline.Usage;
import umms.core.fastq.FastqParser;
import umms.core.fastq.FastqSequence;

/**
 * This class is meant to generally process reads that have barcodes 
 * barcodes can be composed of a Unique Molecular Identifier (UMI) and a sample barcode (BC) 
 * the barcode is assumed to be read always at the same position in each read.
 * UMIs when present will be added to the description and appended to the read name prior to the first space
 * @author mgarber
 *
 */
public class BCProcessor extends CommandLineProgram {
	private static final String PROGRAM_VERSION = "0.01";
	static Logger logger = Logger.getLogger(BCProcessor.class.getName());
	
    //report the region which at least have on fragment score larger than the MINPEAK SCORE
	@Usage 
	public static final String USAGE = "Usage: BCProcessor [options]";
    @Option(doc="Fastq fasta file to process, barcode is assumed to be in this file", shortName="f1", optional=false) 
    public File FASTQ1;
    @Option(doc="If paired reads, this is the fastq file with the other pair", shortName="f2", optional=true) 
    public File FASTQ2;
    @Option(doc="Barcode pattern in the form of SSSSSNNNNOO, where the S indicate the number of nucleotides determining the "+
    "sample barcode and the N indicate the nucleotides defining the UMI, and the O indicate other nucleotides present",
    		shortName="B", optional=false)
    public String BC_FORMAT;
    @Option(doc="Position (zero based) in the read where the barcode is expected to start, default is the first position (0)" , shortName="S", optional=true)
    public int BC_START_POS = 0;
    @Option (doc="2-colum Tab delimited file with barcode to sample mapping" , shortName="M", optional=false)
    public File BC_SAMPLE_MAP;
    @Option (doc="Maximum hamming distance to the sample barcode in order to assign the read to the sample (default 0)", shortName="HD", optional=true)
    public int MAX_HAMMING_DIST = 0;
    @Option (doc="Output directory (default .)", shortName="O", optional=true)
    public File OUTDIR = new File(".");
    
    private int bcLength;
    private int idxSampleBCStart;
	private int idxSampleBCEnd;
	private int idxUMIStart;
	private int idxUMIEnd;
	private int idxOtherStart;
	private int idxOtherEnd;
	private boolean isPaired;
	private Map<String,Integer> sampleBarcodeCounts;
	private Map<String,Integer> randomBarcodeCounts;
	private Map<String,Integer> otherBarcodeCounts;
	
	private Map<String,Integer> sampleNonBarcodeCounts;
	private Map<String,Integer> randomNonBarcodeCounts;
	private Map<String,Integer> otherNonBarcodeCounts;
	
    /** Stock main method. */
    public static void main(final String[] argv) {
        System.exit(new BCProcessor().instanceMain(argv));
    }

	@Override
	protected int doWork() {
		sampleBarcodeCounts = new HashMap<String, Integer>();
		randomBarcodeCounts = new HashMap<String, Integer>();
		otherBarcodeCounts  = new HashMap<String, Integer>();
		
		sampleNonBarcodeCounts = new HashMap<String, Integer>();
		randomNonBarcodeCounts = new HashMap<String, Integer>();
		otherNonBarcodeCounts  = new HashMap<String, Integer>();
		
		HashMap<String, String> bcSampleMap = null;
		try {
			bcSampleMap = readBCMap(BC_SAMPLE_MAP);
		} catch (FileNotFoundException e) {
			System.out.println("Could not access the barcode sample mapping file");
			e.printStackTrace();
			return 1;
		} catch (IOException e) {
			System.out.println("Error accessing barcode sample mapping file");
			e.printStackTrace();
			return 1;
		}
		
		init();
		
		FastqParser parser1 = new FastqParser();
		FastqParser parser2 = new FastqParser();
		
		HashMap<String, BufferedWriter[]> writerMap = null;
		BufferedWriter unassignedP1 = null;
		BufferedWriter unassignedP2 = null;
		
		try {
			try {
				parser1.start(FASTQ1);
			} catch (IOException e) {
				System.out.println("Error reading fastq file " + FASTQ1);
				e.printStackTrace();
				return 1;
			}
			
			if(isPaired) {
				try {
					parser2.start(FASTQ2);
				} catch (IOException e) {
					System.out.println("Error reading fastq file " + FASTQ2);
					e.printStackTrace();
					return 1;
				}
			}			
			try {
				writerMap = buildWriterMap(bcSampleMap) ;
				unassignedP1 = new BufferedWriter(new FileWriter(OUTDIR + "/unassigned.P1.fq"));
				if (isPaired) {
					unassignedP2 = new BufferedWriter(new FileWriter(OUTDIR + "/unassigned.P2.fq"));
				}
			} catch (IOException e) {
				System.out.println("Could not create output files for barcode splitting");
				e.printStackTrace();
				return 1;
			}
				
			
			while (parser1.hasNext()) {
				FastqSequence r1 = parser1.next();
				String seq = r1.getSequence();
				
				String [] bcInfo = parseBarcode(seq); // First position saple, second UMI, third OTHER
				
				BufferedWriter out1 = unassignedP1;
				BufferedWriter out2 = unassignedP2;

				String bc = getClosestBc(bcSampleMap, bcInfo[0]);
				if( bc != null) {
					BufferedWriter [] out = writerMap.get(bc);
					out1 = out[0];
					if(isPaired) { 
						out2 = out[1];
					}
					
					setBarcode(r1, bcInfo);
					r1.excise(this.BC_START_POS,this.bcLength);
					updateCounts(bcInfo);
				} else {
					updateNonCounts(bcInfo);
				}
				
				
				r1.write(out1);
				if( isPaired) {
					FastqSequence r2 = parser2.next();
					if( bc!= null) {
						setBarcode(r2, bcInfo);
					}
					r2.excise(this.BC_START_POS,this.bcLength);
					r2.write(out2);
				}
				
				if(bc != null) {
					updateCounts(bcInfo);
				}
				
			}
			
			writeReports();
			
		} catch (IOException e) {
			System.out.println("Error writing records to file");
			e.printStackTrace();
			return 1;

		} finally {
			try {
				if(parser1 != null )
				parser1.close();
			} catch (IOException e) {
				System.out.println("Could not close fastq1 file");
			}
			if(isPaired) { 
				try {
					parser2.close();
				} catch (IOException e) {
					System.out.println("Could not close fastq2 file");
				}
			}
			
			try {
				if(unassignedP1 != null ) {
					unassignedP1.close();
				}
				if(unassignedP2 != null ) {
					unassignedP2.close();
				}
				if(writerMap != null) {
					closeWriters(writerMap);
				}
			} catch (IOException e) {
				System.err.println("Error closing fastq writers");
				e.printStackTrace();
			}
			
		}
		
		return 0;
	}

	private void writeReports() throws IOException {
		writeReport(OUTDIR + "/sampleBarcodes.rprt", sampleBarcodeCounts);	
		writeReport(OUTDIR + "/randomBarcodes.rprt", randomBarcodeCounts);
		writeReport(OUTDIR + "/extraSeq.rprt", otherBarcodeCounts);
		writeReport(OUTDIR + "/nonSampleBarcodes.rprt", sampleNonBarcodeCounts);
		writeReport(OUTDIR + "/randomBarcodesForNonSamples.rprt", randomNonBarcodeCounts);
		writeReport(OUTDIR + "/extraSeqForNonSamples.rprt", otherNonBarcodeCounts);
	}

	private void writeReport(String outFileName,  Map<String, Integer> countMap) throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(outFileName));
		for (String s : countMap.keySet()) {
			bw.write(s);
			bw.write("\t");
			bw.write(String.valueOf(countMap.get(s)));
			bw.newLine();
		}
		bw.close();
	}

	private void updateCounts(String [] bcInfo) {
		int sampleCounts = sampleBarcodeCounts.containsKey(bcInfo[0]) ? sampleBarcodeCounts.get(bcInfo[0]) + 1 : 1;
		sampleBarcodeCounts.put(bcInfo[0], sampleCounts);
		
		int randomBCCounts = randomBarcodeCounts.containsKey(bcInfo[1]) ? randomBarcodeCounts.get(bcInfo[1]) + 1 : 1;
		randomBarcodeCounts.put(bcInfo[1], randomBCCounts);
		
		int otherCounts = otherBarcodeCounts.containsKey(bcInfo[2]) ? otherBarcodeCounts.get(bcInfo[2]) + 1 : 1;
		otherBarcodeCounts.put(bcInfo[2], otherCounts);
	}
	
	private void updateNonCounts(String [] bcInfo) {
		int sampleCounts = sampleNonBarcodeCounts.containsKey(bcInfo[0]) ? sampleNonBarcodeCounts.get(bcInfo[0]) + 1 : 1;
		sampleNonBarcodeCounts.put(bcInfo[0], sampleCounts);
		
		int randomBCCounts = randomNonBarcodeCounts.containsKey(bcInfo[1]) ? randomNonBarcodeCounts.get(bcInfo[1]) + 1 : 1;
		randomNonBarcodeCounts.put(bcInfo[1], randomBCCounts);
		
		int otherCounts = otherNonBarcodeCounts.containsKey(bcInfo[2]) ? otherNonBarcodeCounts.get(bcInfo[2]) + 1 : 1;
		otherNonBarcodeCounts.put(bcInfo[2], otherCounts);
	}

	private void setBarcode(FastqSequence r, String[] bcInfo) {
		String name = r.getName();
		String compoundBC = bcInfo[0] + "_" + bcInfo[1];
		if (name.contains(" ")) {
			name = name.replace(" ", ":" + compoundBC + " ");
		} else {
			name = name + ":" + compoundBC;
		}
		
		r.setName(name);
		
		String comment = r.getDescription();
		if (comment == null || comment.trim().isEmpty()) {
			comment = compoundBC;
		} else {
			comment = comment + ":" + compoundBC;
		}
		r.setDescription(comment);
		
	}

	protected void init() {
		bcLength = BC_FORMAT.length();
		idxSampleBCStart  = BC_FORMAT.indexOf('S');
		idxSampleBCEnd    = BC_FORMAT.lastIndexOf('S');
		
		idxUMIStart  = BC_FORMAT.indexOf('N');
		idxUMIEnd    = BC_FORMAT.lastIndexOf('N');
		
		idxOtherStart = BC_FORMAT.indexOf('O');
		idxOtherEnd   = BC_FORMAT.lastIndexOf('O');
		isPaired = !(FASTQ2 == null);
	}

	
	private String getClosestBc(HashMap<String, String> bcSampleMap, String observedBC) {
		// TODO Add hemming distance comparison
		return bcSampleMap.containsKey(observedBC) ? observedBC : null;
	}
	
	//	TODO: HANDLE THE CASE WHERE ANY OF THE START/ENDS ARE 0
	private String[] parseBarcode(String seq) {
		String [] info = new String [3];
		
		String bcString = seq.substring(BC_START_POS, BC_START_POS + bcLength +1);

		
		info [0] = bcString.substring(idxSampleBCStart,idxSampleBCEnd + 1);
		info [1] = bcString.substring(idxUMIStart,idxUMIEnd+1);
		// Handle the case where there is no "O" (other) bases:
		if (idxOtherStart>-1) {
			info [2] = bcString.substring(idxOtherStart,idxOtherEnd + 1);
		} else {
			info [2] = "";
		}

		//logger.debug("Read: " + seq + " bc: " + bcString + " sampleBC: " + info[0] + " UMI: " + info[1] + "other "+ info[2]);
		
		return info;
	}

	private void closeWriters(HashMap<String, BufferedWriter[]> writerMap) throws IOException {
		Iterator<BufferedWriter[]> it = writerMap.values().iterator();
		while(it.hasNext()) {
			BufferedWriter[] bw = it.next();
			if(! (bw == null) ) {
				bw[0].flush();
				bw[0].close();
				if(isPaired) {
					bw[1].flush();
					bw[1].close();
				}
			}
		}
		
	}

	private HashMap<String, BufferedWriter[]> buildWriterMap(
			HashMap<String, String> bcSampleMap) throws IOException {

		HashMap<String, BufferedWriter[]> map = new HashMap<String, BufferedWriter[]>();
		Iterator<String> bcIt = bcSampleMap.keySet().iterator();
		while (bcIt.hasNext()) {
			String bc = bcIt.next();
			
			if (! isPaired) {
				BufferedWriter[] writerArray = new BufferedWriter[1];
				String outSampleName = OUTDIR + "/" + bc + "_" + bcSampleMap.get(bc)+ ".fq";
				writerArray[0] = new BufferedWriter(new FileWriter(outSampleName));
				map.put(bc, writerArray);
			} else {
				BufferedWriter[] writerArray = new BufferedWriter[2];
				String outSampleName = OUTDIR + "/" + bc + "_" + bcSampleMap.get(bc)+ ".p1.fq";
				writerArray[0] = new BufferedWriter(new FileWriter(outSampleName));
				String outSampleName2 = OUTDIR + "/" + bc + "_" + bcSampleMap.get(bc)+ ".p2.fq";
				writerArray[1] = new BufferedWriter(new FileWriter(outSampleName2));
				map.put(bc, writerArray);
			}
		}
		
		return map;
	}

	private HashMap<String, String> readBCMap(File file) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(file));
		HashMap<String, String> map = new HashMap<String, String>();
		String line = null;
		while (  (line = br.readLine()) != null) {
			String [] info = line.split("\\s+");
			map.put(info[0], info[1]);
		}
		br.close();
		return map;
	}


}
