package umms.core.fastq.utils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import org.apache.log4j.Logger;

import picard.cmdline.CommandLineProgram;
import picard.cmdline.Option;
import picard.cmdline.Usage;
import umms.core.fastq.FastqParser;
import umms.core.fastq.FastqSequence;
import umms.core.sequence.SequenceUtils;

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
    @Option (doc="Flag to indicate that the barcode is embeded in the read name accordingly to a Casava Specification.", shortName="IRN", optional=true)
    public boolean IN_READ_NAME = false;

    
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
			bcSampleMap = BCEvaluator.readBCMap(BC_SAMPLE_MAP);
		} catch (FileNotFoundException e) {
			System.out.println("Could not access the barcode sample mapping file");
			e.printStackTrace();
			return 1;
		} catch (IOException e) {
			System.out.println("Error accessing barcode sample mapping file");
			e.printStackTrace();
			return 1;
		}

		if(! OUTDIR.exists()) {
			try {
				System.out.println("Output directory does not exists " + OUTDIR+ ", going to create it");
				OUTDIR.mkdir();
			}catch (SecurityException e) {
				System.out.println("Could not create output directory");
				e.printStackTrace();
				return 1;
			}
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
				
				String [] bcInfo = parseBarcode(r1); // First position sample, second UMI, third OTHER
				
				if(bcInfo[0] == null) {
					logger.error("Obtained a null barcode for read " + r1.getName() + " sequence " + r1.getSequence());
				}
				
				BufferedWriter out1 = unassignedP1;
				BufferedWriter out2 = unassignedP2;

				String bc = getClosestBc(bcSampleMap, bcInfo[0]);
				if( bc != null) {
					//Need to update the BCInfo with the closest barcode once it is found
					bcInfo[0] = bc;
					BufferedWriter [] out = writerMap.get(bc);
					out1 = out[0];
					if(isPaired) { 
						out2 = out[1];
					}
					
					setBarcode(r1, bcInfo);
					if( !IN_READ_NAME) {
						r1.excise(this.BC_START_POS,this.bcLength);
					}
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
					//TODO: Handle when/if there is another barcode on R2
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
		
		if(bcInfo[1] != null) {
			int randomBCCounts = randomBarcodeCounts.containsKey(bcInfo[1]) ? randomBarcodeCounts.get(bcInfo[1]) + 1 : 1;
			randomBarcodeCounts.put(bcInfo[1], randomBCCounts);
		}
		
		if(bcInfo[2] != null) {
			int otherCounts = otherBarcodeCounts.containsKey(bcInfo[2]) ? otherBarcodeCounts.get(bcInfo[2]) + 1 : 1;
			otherBarcodeCounts.put(bcInfo[2], otherCounts);
		}
	}
	
	private void updateNonCounts(String [] bcInfo) {
		int sampleCounts = sampleNonBarcodeCounts.containsKey(bcInfo[0]) ? sampleNonBarcodeCounts.get(bcInfo[0]) + 1 : 1;
		sampleNonBarcodeCounts.put(bcInfo[0], sampleCounts);
		
		if(bcInfo[1] != null) {
			int randomBCCounts = randomNonBarcodeCounts.containsKey(bcInfo[1]) ? randomNonBarcodeCounts.get(bcInfo[1]) + 1 : 1;
			randomNonBarcodeCounts.put(bcInfo[1], randomBCCounts);
		}
		
		if(bcInfo[2] != null) {
			int otherCounts = otherNonBarcodeCounts.containsKey(bcInfo[2]) ? otherNonBarcodeCounts.get(bcInfo[2]) + 1 : 1;
			otherNonBarcodeCounts.put(bcInfo[2], otherCounts);
		}
	}

	private void setBarcode(FastqSequence r, String[] bcInfo) {
		String compoundBC = bcInfo[0] + (bcInfo[1] != null && bcInfo[1].length() > 0 ? "_" + bcInfo[1] : "");
		String newName = buildBarcodeFromReadName( r.getName(), compoundBC);
		
		r.setName(newName);
		
		/*String comment = r.getDescription();
		if (comment == null || comment.trim().isEmpty()) {
			comment = compoundBC;
		} else {
			comment = comment + ":" + compoundBC;
		}
		r.setDescription(comment);*/
	}
	
	/**
	 * Cumbersome way to replace the existing barcode with he compound one we built
	 * @param readName
	 * @param compoundBC
	 * @return
	 */
	private String buildBarcodeFromReadName(String readName, String compoundBC) {
		String newName = readName;

		String [] idComponents = readName.split (":");
		
		if (idComponents.length  == 5) { // Casava version 1.4 or earlier
			idComponents[idComponents.length - 1] = idComponents[idComponents.length - 1].replaceFirst("#.+/", "#"+compoundBC+"/");
			newName = join(idComponents);
		} else { //Later Casava version (> 1.8) has barcode at the end.
			idComponents[idComponents.length - 1] = compoundBC;
			newName = join(idComponents);
		}

		/*
		if (name.contains(" ")) {
			name = name.replace(" ", ":" + compoundBC + " ");
		} else {
			name = name + ":" + compoundBC;
		}
		*/

		return newName;
	}

	private String join(String[] idComponents) {
		if(idComponents.length == 0) {
			return "";
		}
		
		StringBuilder sb = new StringBuilder(idComponents[0]);
		for (int i = 1 ; i < idComponents.length; i++) {
			sb.append(":").append(idComponents[i]);
		}
		return sb.toString();
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
		String barcode = bcSampleMap.containsKey(observedBC) ? observedBC : null;
		
		int bestDist = Integer.MAX_VALUE;
		int bcsAtBestDist = 0;
		if (barcode == null) {
			for (String bc : bcSampleMap.keySet()) {
				int dist = SequenceUtils.hamming(bc, observedBC);
				if ( dist <= MAX_HAMMING_DIST && dist <= bestDist) {
					if (bestDist > dist) {
						bestDist = dist;
						bcsAtBestDist = 1;
						barcode = bc;
					} else {
						bcsAtBestDist++;
					}
				}
				
			}
			if( bcsAtBestDist > 1  ) { //TODO: Perhaps be more stringent with other nearby barcodes being farther than a single base from the "one"
				barcode = null;
			}
		} 
		
		return barcode;//TODO: Shall we examine other barcodes that are at the required HD?
	}
	
	//	TODO: HANDLE THE CASE WHERE ANY OF THE START/ENDS ARE 0
	private String[] parseBarcode(FastqSequence r) {
		if(this.IN_READ_NAME ) {
			return parseBarcodeFromName(r.getName());
		} else {
			return parseBarcodeFromSequence( r.getSequence());
		}

	}
	
    /**
     * This now assumes Casava formats
     */
	private String[] parseBarcodeFromName(String name) {
	    String [] idComponents = name.split (":");
	    String bc = idComponents[idComponents.length - 1];
	    
	    if (idComponents.length  == 4) { // Casava version 1.4 or earlier
		bc = bc.replace("#", "");
		bc = bc.replace("/","");
	    } 
	    return splitBarcode(bc);
	}

	private String [] parseBarcodeFromSequence (String seq) {
		
		String bcString = seq.substring(BC_START_POS, BC_START_POS + bcLength +1);		
		return splitBarcode( bcString);	 	
	}

	private String [] splitBarcode( String bcString) {
		String [] info = new String [3];
		
		//This is to handle non-complex barcodes that for example only have the sample id
		if(idxSampleBCStart > -1 ) {
			info [0] = bcString.substring(idxSampleBCStart,idxSampleBCEnd + 1);
			if(idxUMIStart > -1) {
				info [1] = bcString.substring(idxUMIStart,idxUMIEnd+1);
				
			}
			if (idxOtherStart > -1) {
				info [2] = bcString.substring(idxOtherStart,idxOtherEnd + 1);
			}
		} else {
			info[0] = bcString; // assuming that the barcode is the sample id. Not handling the case where the barcode is only the UMI
		}

		return info;
	}

	private void closeWriters(HashMap<String, BufferedWriter[]> writerMap) throws IOException {
		Iterator<BufferedWriter[]> it = writerMap.values().iterator();
		while(it.hasNext()) {
			BufferedWriter[] bw = it.next();
			if(! (bw == null) ) {
				bw[0].close();
				if(isPaired) {
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

}
