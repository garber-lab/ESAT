package broad.pda.samtools;

import java.io.File;
import java.util.Set;
import java.util.HashSet;
import org.apache.commons.io.FileUtils;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.fastq.FastqReader;
import net.sf.picard.fastq.FastqRecord;
import net.sf.picard.fastq.FastqWriter;
import net.sf.picard.fastq.FastqWriterFactory;
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.Log;

public class FilterFastqReads extends CommandLineProgram {
    private static final Log log = Log.getInstance(TrimFastqReads.class);
    
    @Usage
    public String USAGE = "Filters FASTQ files to remove reads matching the provided read name list.";
    
    @Option(doc = "FASTQ file 1", shortName = "I1")
	public File INPUT_1;
    
    @Option(doc = "FASTQ file 2 (assumed paired with file 1)", shortName = "I2", optional=true)
	public File INPUT_2 = null;
    
    @Option(doc = "Output FASTQ 1", shortName = "O1")
    public File OUTPUT_1;
    
    @Option(doc = "Output FASTQ 2", shortName = "O2", optional = true)
    public File OUTPUT_2 = null;
    
    @Option(doc = "Plain text file with one read name per line", shortName = "RLF")
    public File READ_LIST_FILE;  
    
    
	
	/**
	 * Stock main method.
	 *
	 * @param args main arguments
	 */
	public static void main(final String[] args) {
		System.exit(new FilterFastqReads().instanceMain(args));
	}
	
	

    @Override
    protected int doWork() {
    	try {
			IoUtil.assertFileIsReadable(INPUT_1);
			IoUtil.assertFileIsReadable(INPUT_2);
			
			FastqReader reader1 = new FastqReader(INPUT_1);
			FastqReader reader2 = null;
			if (INPUT_2 != null) reader2 = new FastqReader(INPUT_2);
			
			FastqWriter writer1 = new FastqWriterFactory().newWriter(OUTPUT_1);
			FastqWriter writer2 = null;
			if (INPUT_2 != null) writer2 = new FastqWriterFactory().newWriter(OUTPUT_2);

			Set<String> lines = new HashSet<String>();
			lines.addAll(FileUtils.readLines(READ_LIST_FILE, "utf-8"));
			
			int counter = 0;
	
			while (reader1.hasNext()) {
				FastqRecord rec1 = reader1.next();
				FastqRecord rec2 = null;
				if (INPUT_2 != null) rec2 = reader2.next();
				
	            String frec1Name = getReadName(rec1.getReadHeader());
	            String frec2Name = null;
	            if (INPUT_2 != null) frec2Name = getReadName(rec2.getReadHeader());

	            if (INPUT_2 != null && !frec1Name.equals(frec2Name)) {
	            	throw new IllegalStateException("Read names for read 1 and read 2 do not match.");
	            }
	            
				if (!lines.contains(frec1Name)) {
					writer1.write(rec1);
					if (INPUT_2 != null) writer2.write(rec2);
				} else {
					counter++;
				}
			}
			
			log.info("Total reads removed per file: " + counter);
			
			writer1.close();
			if (INPUT_2 != null) writer2.close();
			reader1.close();
			if (INPUT_2 != null) reader2.close();
		
    	} catch (Exception e) {
    		log.error(e);
    		e.printStackTrace();
    		return 1;
    	}
    	return 0;
    }
    
    // COPIED FROM PICARD FastqToSam.jar
    // Read names cannot contain blanks
    public String getReadName(final String fastaqHeader) {
        final int idx = fastaqHeader.indexOf(" ");
        return (idx == -1) ? fastaqHeader : fastaqHeader.substring(0,idx); 
    }
    
}
