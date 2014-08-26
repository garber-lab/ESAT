package broad.pda.samtools;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;

import broad.pda.seq.fastq.FastqParser;
import broad.pda.seq.fastq.FastqSequence;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.Log;

public class TrimFastqReads extends CommandLineProgram {  
    private static final Log log = Log.getInstance(FilterFastqReads.class);
    
    @Usage
    public String USAGE = "Trims FASTQ files to remove bases at each end.";
    
    @Option(doc = "FASTQ file ", shortName = "I")
	public File INPUT;
    
    @Option(doc = "Nucleotides to remove at start", shortName = "S", optional = true)
    public int removeAtStart = 0;
    
    @Option(doc = "Nuclotides to remove at end", shortName = "E", optional = true)
    public int removeAtEnd = 0;
    
    @Option(doc = "Output", shortName = "O")
    public File OUTPUT;  
    
    
	
	/**
	 * Stock main method.
	 *
	 * @param args main arguments
	 */
	public static void main(final String[] args) {
		System.exit(new TrimFastqReads().instanceMain(args));
	}
	
	

    @Override
    protected int doWork() {
    	try {
			IoUtil.assertFileIsReadable(INPUT);
			
			FastqParser reader1 = new FastqParser();

			reader1.start(INPUT);
			
			BufferedWriter writer = new BufferedWriter(new FileWriter(OUTPUT));

			while (reader1.hasNext()) {
				FastqSequence rec = reader1.next();

				if(removeAtEnd > 0)
					rec.trimEndBases(removeAtEnd);
				if(removeAtStart > 0 ) 
					rec.trimStartBases(removeAtStart);
				rec.write(writer);
	
			}
			
			writer.close();
			reader1.close();
    	} catch (Exception e) {
    		log.error(e);
    		e.printStackTrace();
    		return 1;
    	}
    	return 0;
    }
    
 
}
