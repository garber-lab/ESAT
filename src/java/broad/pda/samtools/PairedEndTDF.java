/**
 * 
 */
package broad.pda.samtools;

import java.io.File;
import java.util.Iterator; 


import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.Log;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

import org.apache.commons.lang3.StringUtils;

import nextgen.core.pipeline.LSFJob;

/**
 * @author engreitz
 * September 9, 2012
 */
public class PairedEndTDF  extends CommandLineProgram {

    private static final Log log = Log.getInstance(PairedEndTDF.class);
	
	@Usage
	public String USAGE =
	            "Converts a paired-end SAM or BAM file to a TDF file for IGV,\ncounting coverage by fragments instead of individual reads.\n";

	@Option(doc = "The SAM or BAM file to be converted",
			 optional = false,
			 shortName = StandardOptionDefinitions.INPUT_SHORT_NAME)
	public File INPUT;

	@Option(doc = "TDF file",
	        optional = false, shortName = "O")
	public File OUTPUT;
	
	@Option(doc = "Window size", optional = true, shortName = "W")
	public int WINDOW = 25;
	
	@Option(doc = "Genome build (e.g. mm9)", optional = true)
	public String GENOME = "mm9";
	
	@Option(doc = "Keep intermediate BAM file", optional = true)
	public boolean KEEP_INTERMEDIATE = false;

	//@Option(doc = "Minimum mapping quality", optional = true)
	//public int MIN_MAP_QUALITY = 0;
	
	@Option(doc = "Include duplicates", optional = true)
	public boolean INCLUDE_DUPLICATES = false;
	
	@Option(doc = "Maximum insert size (filter out anything bigger)", optional = true)
	public int MAX_INSERT = 2000;
	
	
	public void extendPairedEndReads(final File inSamOrBam, final File outSamOrBam) throws Exception{

		final SAMPairedEndFileReader inputReader = new SAMPairedEndFileReader(inSamOrBam);
		final SAMFileHeader.SortOrder inputSortOrder = inputReader.getReader().getFileHeader().getSortOrder();
		if (!inputSortOrder.equals(SAMFileHeader.SortOrder.coordinate)) {
			throw new Exception("Input SAM or BAM file must be sorted in coordinate order.");
		}

		final SAMFileWriter outputWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(inputReader.getReader().getFileHeader(), true, outSamOrBam);
		Iterator<SAMRecord> itr = inputReader.getSAMRecordIterator();
		
		int pairCounter = 0;
		while (itr.hasNext()) {
			SAMRecord rec = itr.next();	
			outputWriter.addAlignment(rec);
			pairCounter++;
		}

		log.info(pairCounter + " read pairs found.");
		outputWriter.close();
		inputReader.close();
	}



	@Override
	protected int doWork() {

		File tmpBam = null;

		try {
			IoUtil.assertFileIsReadable(INPUT);
			IoUtil.assertFileIsWritable(OUTPUT);

			tmpBam = new File(INPUT.getAbsolutePath() + ".PairedEndTDF.bam");
			IoUtil.assertFileIsWritable(tmpBam);

			extendPairedEndReads(INPUT, tmpBam);

			String command = "/xchip/igv/tools/igvtools count -w " + WINDOW + " -e 0 ";
			if (INCLUDE_DUPLICATES) command += "--includeDuplicates ";
			command += tmpBam.getAbsolutePath() + " " + OUTPUT.getAbsolutePath() + " " + GENOME;
			System.out.println(command);

			Runtime run=Runtime.getRuntime();
			LSFJob job = new LSFJob(run, command);
			job.submit();
			
			if (!KEEP_INTERMEDIATE) tmpBam.delete();
			
			return 0;

		} catch (Exception e) {
			if (tmpBam.exists() && !KEEP_INTERMEDIATE && !tmpBam.delete()) {
				log.warn("Failed to delete " + tmpBam.getAbsolutePath());
			}

			log.error(e, "Failed to convert " + INPUT.getName());
			return 1;
		}
	}

	@Override
	protected String[] customCommandLineValidation() {
		if (INPUT.equals(OUTPUT)) {
			return new String[]{"INPUT file and OUTPUT file must differ!"};
		}
		return super.customCommandLineValidation();
	}

	/**
	 * Stock main method.
	 *
	 * @param args main arguments
	 */
	public static void main(final String[] args) {
		System.exit(new PairedEndTDF().instanceMain(args));
	}

}
