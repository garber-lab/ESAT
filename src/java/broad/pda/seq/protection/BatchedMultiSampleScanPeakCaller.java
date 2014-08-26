/**
 * 
 */
package broad.pda.seq.protection;

import java.io.File;
import java.io.IOException;

import org.apache.log4j.Level;

/**
 * @author prussell
 *
 */
public class BatchedMultiSampleScanPeakCaller extends MultiSampleScanPeakCaller {

	private String sampleName;
	private String chr;
	private SampleData sampleData;
	
	private BatchedMultiSampleScanPeakCaller(MultiSampleScanPeakCaller m, String sample, String chrName) throws IOException {
		super(m);
		sampleName = sample;
		chr = chrName;
				
		boolean foundSample = false;
		for(SampleData s : allSamples) {
			if(s.getSampleName().equals(sampleName)) {
				if(foundSample) {
					throw new IllegalArgumentException("Found sample name " + sampleName + " twice.");
				}
				sampleData = s;
				foundSample = true;
			}
		}
		if(foundSample) {
			logger.info("Found sample " + sampleName + ".");
		}
		
	}

	public static String[] extendSuperArgsForSampleAndChr(String[] commandArgs, String sampleName, String chrName) {
		String[] rtrn = new String[commandArgs.length + 2];
		for(int i=0; i < commandArgs.length; i++) {
			rtrn[i] = commandArgs[i];
		}
		rtrn[commandArgs.length] = sampleName;
		rtrn[commandArgs.length + 1] = chrName;
		return rtrn;
	}
	
	private static String[] getSuperCommandArgs(String[] extendedCommandArgs) {
		String[] rtrn = new String[extendedCommandArgs.length - 2];
		for(int i=0; i < rtrn.length; i++) {
			rtrn[i] = extendedCommandArgs[i];
		}
		return rtrn;
	}
	
	private static String getSampleName(String[] extendedCommandArgs) {
		return extendedCommandArgs[extendedCommandArgs.length - 2];
	}
	
	private static String getChrName(String[] extendedCommandArgs) {
		return extendedCommandArgs[extendedCommandArgs.length - 1];
	}
	
	private void writePeaks(String outDir) throws IOException {
		String outFile = getPeakBedFileName(sampleData, outDir, chr);
		writeSingleSampleScanPeaks(sampleData, outFile, chr);
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
		
		String[] superArgs = getSuperCommandArgs(args);
		
		String sampleName = getSampleName(args);
		String chrName = getChrName(args);
		MultiSampleScanPeakCaller m = createFromCommandArgs(superArgs);
		BatchedMultiSampleScanPeakCaller b = new BatchedMultiSampleScanPeakCaller(m, sampleName, chrName);
		String outDir = commandLineOutDir(superArgs);
		File o = new File(outDir);
		@SuppressWarnings("unused")
		boolean madeDir = o.mkdir();
		if(!o.exists()) {
			throw new IOException("Could not create directory " + outDir);
		}
		
		b.initializeFilterRejectWriters(chrName, outDir + "/" + FILTER_REJECT_DIR);
		if(commandLineHasDebugFlag(superArgs)) {
			b.setLoggerLevel(Level.DEBUG);
		}
		b.writePeaks(outDir);
		b.closeFilterRejectWriters();
	}

}
