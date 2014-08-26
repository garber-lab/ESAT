package broad.pda.seq.rap;

import java.io.File;
import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.Log;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class UnmarkDuplicates  extends CommandLineProgram  {
    private static final Log log = Log.getInstance(UnmarkDuplicates.class);
    
    @Usage 
    public final String USAGE = "Remove the duplicate flag from all reads in a SAM or BAM file.";
    
    @Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="Input SAM or BAM files to analyze.")
    public File INPUT;
    
    @Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="The output file to right marked records to")
    public File OUTPUT;
    
	
	/**
	 * Stock main method.
	 *
	 * @param args main arguments
	 */
	public static void main(final String[] args) {
		System.exit(new UnmarkDuplicates().instanceMain(args));
	}
	

	@Override
	protected int doWork() {
		
		try {
			
	        IoUtil.assertFileIsReadable(INPUT);
	        IoUtil.assertFileIsWritable(OUTPUT);
			
	        final SAMFileReader reader = new SAMFileReader(INPUT);
	        final SAMFileHeader header = reader.getFileHeader();
	        final SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, true, OUTPUT);
	        SAMRecordIterator itr = reader.iterator();
	        while (itr.hasNext()) {
	        	SAMRecord record = itr.next();
	        	record.setDuplicateReadFlag(false);
	        	writer.addAlignment(record);
	        }
	        itr.close();
	        writer.close();
	        reader.close();
	        
		} catch (Exception e) {
			log.error(e);
		}
		
		return 0;
	}
}
