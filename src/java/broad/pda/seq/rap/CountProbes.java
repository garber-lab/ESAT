package broad.pda.seq.rap;

import java.io.File;
import java.io.BufferedWriter;
import java.io.FileWriter;

import broad.core.datastructures.Pair;
import broad.pda.samtools.PairedEndReadIterator;
import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.Log;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMRecord;


/**
 * Program to find read pairs in which one of the pairs aligns to at least part 
 * of the qPCR tag sequence of the probe.  Input SAM files have been aligned to 
 * the full probe sequences (including the target sequence, e.g. of Xist).
 * @author engreitz
 *
 */

public class CountProbes extends CommandLineProgram {
    private static final Log log = Log.getInstance(CountProbes.class);
    
    @Option(doc="Input SAM file, aligned to Probe FASTA file with Bowtie 2.\n NOTE: Read pairs should occupy consecutive lines.")
    public File INPUT;
    
    @Option(doc="Output SAM file, containing reads that contain probe qPCR tag sequence")
    public File OUTPUT;
    
    @Option(doc="Output tabbed text file containing start and end indices")
    public File STATS;
    
    @Option(doc="Number of extra T7 bases added in the Bowtie index")
    public Integer T7BASES = 26;
    
    @Option(doc="Minimum number of bases overlapping with the qPCR tag to call a hit")
    public Integer OVERLAP = 5;
	
    // Define constants describing how many bases on each end of the probe are qPCR tag sequence
    // (non-target sequence that does not reside in the genome)
    
    public static final int ORIGINAL_ARRAY_LEFT = 30;
    public static final int ORIGINAL_ARRAY_RIGHT = 20;
    public static final int MAIN_ARRAY_LEFT = 30;
    public static final int MAIN_ARRAY_RIGHT = 30;
    public static final int SIMPLE_ARRAY_LEFT = 15;
    public static final int SIMPLE_ARRAY_RIGHT = 15;
    public static final int ARRAY_140105 = 17;
    public static final int MRNA_ARRAY_LEFT = 40;
    public static final int MRNA_ARRAY_RIGHT = 40;
    
    // Keep track of how often the probe sequences map to specific starting / ending locations.
    // Arrays automatically initialize to zero.
    // In practice, this system isn't perfect, because the read pairs might be aligned to different
    //  overlapping probes that contain the same target sequence.
    private int[] startIndex = new int[250];
    private int[] endIndex = new int[250];
    
	/**
	 * Stock main method.
	 *
	 * @param args main arguments
	 */
	public static void main(final String[] args) {
		System.exit(new CountProbes().instanceMain(args));
	}
	

    @Override
    protected int doWork() {
        IoUtil.assertFileIsReadable(INPUT);
        IoUtil.assertFileIsWritable(OUTPUT);
        IoUtil.assertFileIsWritable(STATS);

        try {
        	// Read and filter
        	SAMFileReader reader = new SAMFileReader(INPUT);
        	final SAMFileHeader header = reader.getFileHeader();

        	// In our current pipeline, unsorted BAM files output by BWA alignment have read pairs next to each other.
        	if (header.getSortOrder() == SAMFileHeader.SortOrder.unsorted) {
        		log.warn("Treating unsorted SAM file as queryname-sorted.  Check that this produces correct behavior.");
        	}
        	header.setSortOrder(SAMFileHeader.SortOrder.queryname);
        	
        	PairedEndReadIterator itr = new PairedEndReadIterator(reader);
        	
        	header.setSortOrder(SAMFileHeader.SortOrder.unsorted);
        	SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, false, OUTPUT);
        	
        	while (itr.hasNext()) {
        		Pair<SAMRecord> pair = itr.next();
        		SAMRecord first = pair.getValue1();
        		SAMRecord second = pair.getValue2();
        		
        		if (readMapsToProbe(first, header) || readMapsToProbe(second, header)) {
        			writer.addAlignment(first);
        			writer.addAlignment(second);
        			addStartEndIndices(first, second);
        		}
        	}
        	
        	reader.close();
        	writer.close();
        	
        	BufferedWriter bw = new BufferedWriter(new FileWriter(STATS));
    		bw.write("index\tstart\tend\n");
        	for (int i = 0; i < startIndex.length; i++) {
        		bw.write(i + "\t" + startIndex[i] + "\t"+ endIndex[i] + "\n");
        	}
        	bw.close();
        	
        } catch (Exception e) {
        	log.error(e);
        }

        return 0;
    }	
    
    
    private boolean readMapsToProbe(SAMRecord read, SAMFileHeader header) {
    	if (read.getReadUnmappedFlag()) return false;
    	
    	int left, right;
    	if (read.getReferenceName().contains("OriginalArray")) {
    		left = ORIGINAL_ARRAY_LEFT;
    		right = ORIGINAL_ARRAY_RIGHT;
    	} else if (read.getReferenceName().contains("MainArray")) {
    		left = MAIN_ARRAY_LEFT;
    		right = MAIN_ARRAY_RIGHT;
    	} else if (read.getReferenceName().contains("SimpleArray")) {
    		left = SIMPLE_ARRAY_LEFT;
    		right = SIMPLE_ARRAY_RIGHT;
    	} else if (read.getReferenceName().contains("140105")) {
    		left = ARRAY_140105;
    		right = ARRAY_140105;
    	} else if (read.getReferenceName().contains("Mouse")) {
    		left = MRNA_ARRAY_LEFT;
    		right = MRNA_ARRAY_RIGHT;
    	} else {
    		log.error("This reference name doesn't match: " + read.getReferenceName());
    		return false;
    	}
    	
    	left = T7BASES + left;
    	right = header.getSequence(read.getReferenceIndex()).getSequenceLength() - right - T7BASES;
    	
    	if ((read.getAlignmentStart() <= left - OVERLAP) || read.getAlignmentEnd() >= right + OVERLAP) {
    		return true;
    	}
    	
    	return false;
    }
    	
    
    private void addStartEndIndices(SAMRecord first, SAMRecord second) {
    	if (!first.getReadUnmappedFlag() && !second.getReadUnmappedFlag()) {
    		int start = Math.min(first.getAlignmentStart(), second.getAlignmentStart());
    		int end = Math.max(first.getAlignmentEnd(), second.getAlignmentEnd());
    		startIndex[start]++;
    		endIndex[end]++;
    	}
    }
}
