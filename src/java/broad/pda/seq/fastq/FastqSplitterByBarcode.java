package broad.pda.seq.fastq;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;

import broad.core.sequence.Sequence;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;

/**
 * Reads through a fastq file and writes records to different files depending on 
 * the barcode found within the start or end of the read. Reads are trimmed to remove
 * the barcode.
 * @author mgarber
 *
 */
public class FastqSplitterByBarcode {
	static Logger logger = Logger.getLogger(FastqSplitterByBarcode.class.getName());
	static private String UNMATCHED = "unmatched";
	static private Pattern nameBarcodePattern = Pattern.compile("#[A,C,G,T,N]+");
	//static private Pattern endOfReadMamePairInfo = Pattern.compile("\\(/[0-9]+$\\)");
	static private Pattern endOfReadMamePairInfo = Pattern.compile("(/[0-9]+$)");

	private HashMap<String, String> sampleToBC;
	private HashMap<String, String> reverseBCs;
	private HashMap<String, String> bcToSample;
	private HashMap<String, Integer> readsPerBarcode;
	private boolean atEnd  = false;
	private boolean paired = false;
	
	private int read2Trim5p = 0;
	private int read2Trim3p = 0;
	private int read1Trim3p = 0;
	private int read1Trim5p = 0;
	private boolean addToReadName;
	
	private FastqParser pair1Parser = new FastqParser();
	private FastqParser pair2Parser = new FastqParser();
	
	private HashMap<String, BufferedWriter> pair1Writers = new HashMap<String, BufferedWriter>();
	private HashMap<String, BufferedWriter> pair2Writers = new HashMap<String, BufferedWriter>();
	
	

	static final String USAGE = "FastqSplitterByBarcode. Reads a fastq file and writes reads to \ndifferent files depending on their barcode: " +
			"\n\t-in <Path to the fastq file (NO standard input is supported must input a file paht)>" +
			"\n\t-outdir <Ouput directory for split files>" +
			"\n\t-bcfile <File of the form sample\tbarcode>" +
			"\n\t-pair2File <If paired end this specifies the path to the second read file. The barcode is assumed to be on the first file>" +
			"\n\t-atEnd <add this flag is the barcode is expected at the end, not the begginging of the read>" +
			"\n\t-read1Trim5p <number of bases to trim the first read (after barcode)>" +
			"\n\t-read1Trim3p <number of bases to trim the first read (after barcode)>" +
			"\n\t-read2Trim5p <number of bases to trim the second read>" +
			"\n\t-read2Trim3p <number of bases to trim the second read>" +
			"\n\t-addToReadName <Add this flag if the trimmed sequence from the second read should added to the read name's end after a _BC_ separator>" +
			"\n";
	
	
	
	
	public FastqSplitterByBarcode(String bcFile) throws IOException {
		sampleToBC = new HashMap<String, String>();
		bcToSample = new HashMap<String, String>();
		readsPerBarcode = new HashMap<String, Integer>();
		reverseBCs  = new HashMap<String, String>();
		
		BufferedReader br = new BufferedReader(new FileReader(bcFile));
		String line = null;
		while ( (line = br.readLine()) != null) {
			String [] info = line.split("\\s+");
			sampleToBC.put(info[0], info[1]);
			bcToSample.put(info[1], info[0]);
			readsPerBarcode.put(info[1], 0);
			reverseBCs.put(info[1], Sequence.reverseSequence(info[1]) );
		}
		br.close();
		readsPerBarcode.put(UNMATCHED, 0);
	}

	public static void main (String [] args) throws Exception {
		Globals.setHeadless(true);
		ArgumentMap argMap = CLUtil.getParameters(args,USAGE , "splitfile");
		File outDir = new File(argMap.getOutputDir());
		if(! outDir.exists()) {
			outDir.mkdir();
		}
		
		FastqSplitterByBarcode fsbb = new FastqSplitterByBarcode(argMap.getMandatory("bcfile"));
		
		// Setting trimming parameters
		fsbb.setTrimmingParams(argMap);
		fsbb.setAddToReadName(argMap.containsKey("addToReadName"));
		fsbb.setAtEnd(argMap.containsKey("atEnd"));

		fsbb.start(argMap.getInput(), argMap.getOutputDir(), argMap.get("pair2File"));
		
		int numReads = 0;
		while(fsbb.hasMore() ) {
			fsbb.processRecord();
			if(numReads++ % 500000 == 0) {
				logger.info(numReads + " processed");
			}
		}
		
		fsbb.stop();
		fsbb.writeReport();
	
	}

	private void writeReport() {
		System.out.println("Barced\tReadNumber");
		for (String bc : bcToSample.keySet()) {
			System.out.println(bc + "\t" + readsPerBarcode.get(bc));
		}
		System.out.println(UNMATCHED + "\t" + readsPerBarcode.get(UNMATCHED));
	}

	public void stop() throws IOException {
		pair1Parser.close();
		
		for (BufferedWriter bw : pair1Writers.values()) {
			bw.close();
		}
		
		if(isPaired()) {
			pair2Parser.close();
			for (BufferedWriter bw : pair2Writers.values()) {
				bw.close();
			}
		}
	}

	private boolean processRecord() throws IOException {
		FastqSequence pair1 = pair1Parser.next();
		String readSeq = pair1.getSequence();
	
		FastqSequence pair2 = null;
		if(isPaired())  {
			pair2 = pair2Parser.next();
			if (! checkPairNames(pair1, pair2)) {
				logger.error("Incompatible read pair names. Read1 name: " + pair1.getName() + " read2 name: " + pair2.getName());
				// SHould break execution
			}
		}
		
		boolean found = false;
		Iterator<String> bcIt = bcToSample.keySet().iterator();
		while (!found & bcIt.hasNext()) {
			String bc = bcIt.next();
			if(atEnd) {
				found = readSeq.endsWith(bc) || readSeq.endsWith(reverseBCs.get(bc));
			} else {
				found = readSeq.startsWith(bc) || readSeq.startsWith(reverseBCs.get(bc));
			}
			if(found) {
				readsPerBarcode.put(bc, readsPerBarcode.get(bc) + 1);
				trim(pair1, pair2, bc);
				updateReadInfoWithBarcode(pair1, bc);
				pair1.write(pair1Writers.get(bc));
				if (isPaired()) {
					updateReadInfoWithBarcode(pair2, bc);
					pair2.write(pair2Writers.get(bc));
				}
			}
		}
		
		if (! found) {
			pair1.write(pair1Writers.get(UNMATCHED));
			if (isPaired()) {
				pair2.write(pair2Writers.get(UNMATCHED));
			}
			readsPerBarcode.put(UNMATCHED, readsPerBarcode.get(UNMATCHED) + 1);
		}
		
		return found;
	}

	private boolean checkPairNames(FastqSequence pair1, FastqSequence pair2) {
	
		boolean rtrn = true;
		
		if (pair2 != null) {
			String [] pair1NameInfo = pair1.getName().split("\\s");
			String [] pair2NameInfo = pair2.getName().split("\\s");
			rtrn = pair1NameInfo[0].equals(pair2NameInfo[0]);
		}
		return rtrn;
	}

	private void trim(FastqSequence pair1, FastqSequence pair2, String bc) {
		if(atEnd) {
			pair1.trimEndBases(bc.length() + getRead1Trim3p());
			pair1.trimStartBases(getRead1Trim5p());
		} else {
			pair1.trimEndBases(getRead1Trim3p());
			pair1.trimStartBases(bc.length() +  getRead1Trim5p());
		}
		if(pair2 != null) {
			pair2.trimStartBases(getRead2Trim5p());
			pair2.trimEndBases(getRead2Trim3p());
		}
	}

	protected void updateReadInfoWithBarcode(FastqSequence read, String bc) {
		if(read.getDescription().equals(read.getName()) || read.description.length() < 2) {
			read.setDescription(bc);
		} else {
			read.setName(read.getName() + "__BC__" + bc);
		}
	}

	public void start(String pair1Source, String outputDir, String pair2Source) throws IOException {
		File pair1SourceFile = checkFile(pair1Source);
		if(pair1SourceFile == null) {
			throw new IOException("Fastq File for first read " + pair1Source + " could not be found or can't be read");
		}
		
		logger.info("pair1 source " + pair1Source + " pair2 source " + pair2Source);
		pair1Parser.start(pair1SourceFile);
		
		if(pair2Source != null & pair2Source.trim().length() > 0) {
			this.paired = true;
			File pair2SourceFile = checkFile(pair2Source);
			pair2Parser.start(pair2SourceFile);
		}	
		
		ArrayList<String> bcs = new ArrayList<String>(bcToSample.keySet());
		bcs.add(UNMATCHED);
		Iterator<String> bcIt = bcs.iterator();
		while(bcIt.hasNext()) {
			String bc = bcIt.next();
			String sample = bcToSample.containsKey(bc) ? bcToSample.get(bc) : bc;
			logger.debug("processing bc " + bc + " sample " + sample);
			String bcPair1Dest = outputDir + "/" + removePathSeparators(removeExt(pair1Source)) + "." + sample.replaceFirst("\\s", "_") + ".fq";
			File bcPair1DestFile = new File(bcPair1Dest);
			if(bcPair1DestFile.exists() ){
				bcPair1DestFile.delete();
			}
 			BufferedWriter bcPair1Writer = new BufferedWriter( new FileWriter(bcPair1Dest));
			this.pair1Writers.put(bc, bcPair1Writer);
			if(isPaired()) {
				String bcPair2Dest = outputDir + "/" + removePathSeparators(removeExt(pair2Source)) +"." + sample.replaceFirst("\\s", "_") + ".fq";
				File bcPair2DestFile = new File(bcPair2Dest);
				if(bcPair2DestFile.exists() ){
					bcPair2DestFile.delete();
				}
				BufferedWriter bcPair2Writer = new BufferedWriter( new FileWriter(bcPair2Dest));
				this.pair2Writers.put(bc, bcPair2Writer);
			}
		}
		
	}
	


	private String removePathSeparators(String fileName) {
		File f = new File(fileName);
		return f.getName();
	}

	private String removeExt(String source) {
		String sourceNoExt = source.replaceFirst("\\.[^\\.]+$", "");
		return sourceNoExt;
	}

	private File checkFile(String source) {
		File sourceFile = null;
		if( source != null) {
			sourceFile = new File(source);
			if(!sourceFile.exists() & !sourceFile.canRead()) {
				sourceFile = null;
			}
		}
		return sourceFile;
	}

	private void setTrimmingParams(ArgumentMap argMap) {

		setRead1Trim5p(argMap.isPresent("read1Trim5p") ? argMap.getInteger("read1Trim5p") : 0);
		setRead1Trim3p(argMap.isPresent("read1Trim3p") ? argMap.getInteger("read1Trim3p") : 0);
		setRead2Trim5p(argMap.isPresent("read2Trim5p") ? argMap.getInteger("read2Trim5p") : 0);
		setRead2Trim3p(argMap.isPresent("read2Trim3p") ? argMap.getInteger("read2Trim3p") : 0);
	}
	
	public boolean hasMore() { return this.pair1Parser.hasNext();}


	public boolean isPaired() {
		return paired;
	}

	public void setPaired(boolean paired) {
		this.paired = paired;
	}

	public boolean isAtEnd() {
		return atEnd;
	}

	public void setAtEnd(boolean atEnd) {
		this.atEnd = atEnd;
	}

	public int getRead2Trim5p() {
		return read2Trim5p;
	}

	public void setRead2Trim5p(int read2Trim5p) {
		this.read2Trim5p = read2Trim5p;
	}

	public int getRead2Trim3p() {
		return read2Trim3p;
	}

	public void setRead2Trim3p(int read2Trim3p) {
		this.read2Trim3p = read2Trim3p;
	}

	public int getRead1Trim3p() {
		return read1Trim3p;
	}

	public void setRead1Trim3p(int read1Trim3p) {
		this.read1Trim3p = read1Trim3p;
	}

	public int getRead1Trim5p() {
		return read1Trim5p;
	}

	public void setRead1Trim5p(int read1Trim5p) {
		this.read1Trim5p = read1Trim5p;
	}

	public boolean isAddToReadName() {
		return addToReadName;
	}

	public void setAddToReadName(boolean addToReadName) {
		this.addToReadName = addToReadName;
	}

	protected FastqParser getPair1Parser() {
		return pair1Parser;
	}

	protected void setPair1Parser(FastqParser pair1Parser) {
		this.pair1Parser = pair1Parser;
	}

	protected FastqParser getPair2Parser() {
		return pair2Parser;
	}

	protected void setPair2Parser(FastqParser pair2Parser) {
		this.pair2Parser = pair2Parser;
	}

	
}
