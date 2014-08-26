package broad.pda.seq.rap;

import java.io.*;
import java.util.Arrays;

import net.sf.picard.cmdline.Usage;
import net.sf.picard.cmdline.Option;
import net.sf.picard.io.IoUtil;
import net.sf.picard.reference.*;
import net.sf.picard.util.Log;
import net.sf.samtools.util.*;
import net.sf.samtools.*;

import nextgen.core.annotation.*;

public class GCContent extends GenomeCommandLineProgram {
    private static final Log log = Log.getInstance(GCContent.class);
    
    
    @Usage public final String USAGE = "Calculates % GC in windows across the genome and outputs in bedgraph format.";
 
    @Option(doc="Reference sequence file (FASTA)", shortName="REF")
    public File REFERENCE_SEQUENCE;
    
    @Option(doc="Output file (bedgraph)", shortName="O")
    public File OUTPUT;
    
    @Option(doc="Window size")
    public Integer WINDOW = 100;
    
 
    private final int CHUNK_SIZE = 10000;  // read this many windows at a time
    
    
	
	/**
	 * Stock main method.
	 *
	 * @param args main arguments
	 */
	public static void main(final String[] args) {
		System.exit(new GCContent().instanceMain(args));
	}
	
    @Override
    protected int doWork() {
        IoUtil.assertFileIsReadable(REFERENCE_SEQUENCE);
        IoUtil.assertFileIsWritable(OUTPUT);

        try {
        	
        	final ReferenceSequenceFile ref = ReferenceSequenceFileFactory.getReferenceSequenceFile(REFERENCE_SEQUENCE);

        	BufferedWriter bw = new BufferedWriter(new FileWriter(OUTPUT));

        	for (SAMSequenceRecord seq : ref.getSequenceDictionary().getSequences()) {
        		int length = seq.getSequenceLength();

        		for (int i = 1; i < length; i += WINDOW * CHUNK_SIZE) {

        			int end = Math.min(i + WINDOW*CHUNK_SIZE - 1, seq.getSequenceLength());
        			ReferenceSequence chunk = ref.getSubsequenceAt(seq.getSequenceName(), i, end);

        			for (int j = 0; j < chunk.length(); j += WINDOW) {
        				byte[] toScore = Arrays.copyOfRange(chunk.getBases(), j, Math.min(j + WINDOW, chunk.length()));
        				double gc = calculateGc(toScore);

        				if (gc >= 0) {
        					Annotation a = new BasicAnnotation(seq.getSequenceName(), i + j - 1, i + j + WINDOW - 1);
        					a.setScore(gc);
        					bw.write(a.toBEDGraph());
        					bw.newLine();
        				}
        			}
        			bw.flush();
        		}
				
        	}
        	
        	bw.close();
        	
        } catch (Exception e) {
        	log.error(e);
        }

        return 0;
    }	
    	
    
    
    /**
     * Returns -1 if no bases are ATCG (e.g. all N)
     * @param bases
     * @return
     */
    public double calculateGc(final byte[] bases) {
    	int gcs = 0, total = 0;
    	for (int i = 0; i < bases.length; i++) {
    		final byte b = bases[i];
    		if (SequenceUtil.isValidBase(b)) {
    			total++;
    			if (b == 'C' || b == 'G' || b == 'c' || b == 'g') ++gcs;
    		}
    	}
    	
    	double result = -1;
    	if (total > 0) {
    		result = gcs / (double) total;
    	}
    	return result;
    }
}
