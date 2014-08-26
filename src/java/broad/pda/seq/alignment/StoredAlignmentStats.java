package broad.pda.seq.alignment;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

import broad.pda.seq.segmentation.AlignmentDataModelStats;


/**
 * Class that manages the storage of an AlignmentDataModelStats to permanent storage
 * Format of the file is simple, the first line stores the timestamp of file creation
 * the next two columns are tab separated lines containing information about the alignment and index file when the stats were stored:
 * Size, last modified date
 * The next line stores the global RPKM constant
 * The next line should be a P if the alignment is paired and anything else, usually U if it is not.
 * Then there is a line per reference sequence containing its statistics: # markers, # reads, masked regions, lambda,  
 * @author mgarber
 *
 */
public class StoredAlignmentStats {
	public static final String SAVED_MODEL_STATS_EXT = ".srptr.aln.stats";
	
	private long lastModifiedStoredTimeStamp;
	private long lastModifiedAlignmentTimeStampAtCreation;
	private long lastModifiedAlignmentIdxTimeStampAtCreation;
	private long alignmentSizeAtCreation;
	private long alignmentIdxSizeAtCreation;
	private boolean isAlignmentPaired;
	private double globalRPKMConstant;
	
	private Map<String, Double> numMarkers;
	private Map<String, Double> numberOfReads;
	private Map<String, Integer> maskedRegions;
	private Map<String, Double> lambda;
	private Map<String, Double> localRpkmConstant;
	
	private boolean isSafeToLoad;
	private String reasonForUnsafe;

	public StoredAlignmentStats() throws IOException {
		super();
		
	}
	
	public void load(File alignmentStoredStatFile) throws NumberFormatException, IOException {
		BufferedReader br = new BufferedReader(new FileReader(alignmentStoredStatFile));
		String line = null;
		
		numMarkers = new HashMap<String, Double>();
		numberOfReads = new HashMap<String, Double>();
		maskedRegions = new HashMap<String, Integer>();
		lambda        = new HashMap<String, Double>();
		localRpkmConstant = new HashMap<String, Double>();
		
		int lineNum = 0;
		while( (line = br.readLine()) != null) {
			String [] info = line.split("\t");
			switch (lineNum){
			case 0:
				lastModifiedStoredTimeStamp = Long.parseLong(line);
				break;
			case 1:
				alignmentSizeAtCreation = Long.parseLong(info[0]);
				lastModifiedAlignmentTimeStampAtCreation =  Long.parseLong(info[1]);
				break;
			case 2:
				alignmentIdxSizeAtCreation = Long.parseLong(info[0]);
				lastModifiedAlignmentIdxTimeStampAtCreation =  Long.parseLong(info[1]);
				break;
			case 3:
				globalRPKMConstant = Double.parseDouble(line);
				break;
			case 4:
				isAlignmentPaired = "P".equals(line);
				break;
			default: 
				String ref = info[0];
				numMarkers.put(ref, Double.parseDouble(info[1]));
				numberOfReads.put(ref, Double.parseDouble(info[2]));
				maskedRegions.put(ref, (int) Double.parseDouble(info[3]));
				lambda.put(ref, Double.parseDouble(info[4]));
				localRpkmConstant.put(ref, globalRPKMConstant);
			}
			
			lineNum++;
		}
		br.close();
	}	
		
	
	public static void store(AlignmentDataModelStats alignmentStats) throws IOException {
		String pairedString = alignmentStats.isPaired() ? "P" : "U";
		String alignmentPath = alignmentStats.getData().getModelFilePath();
		File alignmentFile = new File(alignmentPath);
		File alignmentIdxFile = new File(alignmentPath + (alignmentPath.endsWith(".sam") ? ".sai" : ".bai"));
		
		File tmpFile = File.createTempFile(".temporary.scrptr", ".tmp");
		BufferedWriter bw = new BufferedWriter(new FileWriter(alignmentPath + SAVED_MODEL_STATS_EXT));
		bw.write(String.valueOf(tmpFile.lastModified()));
		bw.newLine();
		tmpFile.delete();
		bw.write(alignmentFile.length() +"\t" + alignmentFile.lastModified());
		bw.newLine();
		bw.write(alignmentIdxFile.length() +"\t" + alignmentIdxFile.lastModified());
		bw.newLine();
		bw.write(String.valueOf(alignmentStats.getGlobalRPKMConstant()));
		bw.newLine();
		bw.write(pairedString);
		bw.newLine();
		
		Collection<String> chrs = alignmentStats.getChromosomeLengths().keySet();
		for(String chr : chrs ) {
			bw.write(chr);
			bw.write("\t");

			bw.write(alignmentStats.getNumberMarkers(chr)+"\t" +alignmentStats.getNumberOfReads(chr) + 
					"\t" + alignmentStats.getMaskedRegionsSize(chr) +"\t" + alignmentStats.getLambda(chr) );
			bw.newLine();
		}
		bw.close();
	}

	public long getLastModifiedStoredTimeStamp() {
		return lastModifiedStoredTimeStamp;
	}

	public long getLastModifiedAlignmentTimeStampAtCreation() {
		return lastModifiedAlignmentTimeStampAtCreation;
	}

	public long getLastModifiedAlignmentIdxTimeStampAtCreation() {
		return lastModifiedAlignmentIdxTimeStampAtCreation;
	}

	public long getAlignmentSizeAtCreation() {
		return alignmentSizeAtCreation;
	}

	public long getAlignmentIdxSizeAtCreation() {
		return alignmentIdxSizeAtCreation;
	}

	public boolean isAlignmentPaired() {
		return isAlignmentPaired;
	}

	public double getGlobalRPKMConstant() {
		return globalRPKMConstant;
	}

	public Map<String, Double> getNumMarkers() {
		return numMarkers;
	}

	public Map<String, Double> getNumberOfReads() {
		return numberOfReads;
	}

	public Map<String, Integer> getMaskedRegions() {
		return maskedRegions;
	}

	public Map<String, Double> getLambda() {
		return lambda;
	}

	public Map<String, Double> getLocalRpkmConstant() {
		return localRpkmConstant;
	}
	

	protected boolean isSafeToLoad() {
		return isSafeToLoad;
	}

	protected void setAsUnsafeToLoad(){
		this.isSafeToLoad = false;
	}
	
	protected void setAsSafeToLoad(){
		this.isSafeToLoad = true;
	}

	protected String getReasonForUnsafe() {
		return reasonForUnsafe;
	}

	protected void setReasonForUnsafe(String reasonForUnsafe) {
		this.reasonForUnsafe = reasonForUnsafe;
	}


}
