package broad.pda.samtools;

import java.io.File;
import java.io.FileWriter;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.Log;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import broad.core.datastructures.Pair;
import broad.pda.seq.rap.CountProbes;

/**
 * @author engreitz
 * Remove reads that do not map uniquely to the genome.  If a read's mate is removed this way, 
 * mark that read as unpaired.  
 */
public class FilterPairedReadsByMAQ extends CommandLineProgram {
	private static final Log log = Log.getInstance(CountProbes.class);
	    
	@Option(shortName="I", doc="Input SAM file, aligned to Probe FASTA file with Bowtie 2.\n NOTE: Read pairs should occupy consecutive lines.")
	public File INPUT;

	@Option(shortName="MAQ", doc="Reads that have MAQ scores less than this will be filtered out.")
	public Integer MINIMUM_MAPPING_QUALITY = 30;
	
	@Option(shortName="O", doc="Output SAM file, containing reads that contain probe qPCR tag sequence")
	public File OUTPUT;
	
	@Option(doc="If true, remove both reads if one of the pair does not pass the filter.")
	public Boolean REMOVE_MATE = false;
	
	@Option(doc="Optional metrics file")
	public File METRICS=null;
	
	@Option(shortName=StandardOptionDefinitions.SORT_ORDER_SHORT_NAME, doc="Sort order of output file")
	public SAMFileHeader.SortOrder SORT_ORDER = SAMFileHeader.SortOrder.coordinate;
	
	private int passed = 0;
	private int failed = 0;
	private int passedAsPairs = 0;
	private int total = 0;
	
	/**
	 * Stock main method.
	 *
	 * @param args main arguments
	 */
	public static void main(final String[] args) {
		System.exit(new FilterPairedReadsByMAQ().instanceMain(args));
	}
	
    @Override
    protected int doWork() {
        IoUtil.assertFileIsReadable(INPUT);
        IoUtil.assertFileIsWritable(OUTPUT);

        try {
           	// Read and filter
        	SAMFileReader reader = new SAMFileReader(INPUT);
        	SAMFileHeader header = reader.getFileHeader();
        	
        	// In our current pipeline, unsorted BAM files output by BWA alignment have read pairs next to each other.
        	if (header.getSortOrder() == SAMFileHeader.SortOrder.unsorted) {
        		log.warn("Treating unsorted SAM file as queryname-sorted.  Check that this produces correct behavior.");
        	}
        	header.setSortOrder(SAMFileHeader.SortOrder.queryname);
        	
        	PairedEndReadIterator itr = new PairedEndReadIterator(reader);

        	header.setSortOrder(SORT_ORDER);
        	SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, false, OUTPUT);

        	while (itr.hasNext()) {
        		Pair<SAMRecord> pair = itr.next();
        		SAMRecord first = pair.getValue1();
        		SAMRecord second = pair.getValue2();
        		
        		total += 1;
        		
        		boolean acceptFirst = false;
        		boolean acceptSecond = false;
        		
        		if (REMOVE_MATE) {
        			if (acceptMappingQuality(first) && acceptMappingQuality(second)) {
        				acceptFirst = true;
        				acceptSecond = true;
        			}
        		} else {
        			if (acceptMappingQuality(first)) {
        				acceptFirst = true;
        			}
        			if (acceptMappingQuality(second)) {
        				acceptSecond = true;
        			}
        		}
        		
        		if (acceptFirst) {
        			if (!acceptSecond) {
        				removeMate(first);
        			}
        			writer.addAlignment(first);
        			passed++;
        		} else {
        			failed++;
        		}
        		
        		if (acceptSecond) {
        			if (!acceptFirst) {
        				removeMate(second);
        			}
        			writer.addAlignment(second);
        			passed++;
        		} else {
        			failed++;
        		}
        		
        		if (acceptFirst && acceptSecond) passedAsPairs += 2;
        	}        	
        	
        	log.info("Total single end reads skipped: " + itr.getSingleReadsSkipped());
        	log.info("Total reads accepted: " + passed);
        	log.info("Total reads removed: " + failed);
        	log.info("Total reads accepted as pairs: " + passedAsPairs);
        	itr.close();
        	writer.close();
        	reader.close();
        	
        	if (METRICS != null) {
        		FileWriter fw = new FileWriter(METRICS);
        		fw.write("TOTAL_READS\tTOTAL_PAIRS\tREADS_ACCEPTED\tREADS_REMOVED\tREADS_ACCEPTED_AS_PAIRS\n");
        		fw.write((total*2 + itr.getSingleReadsSkipped()) + "\t" + total + "\t" + passed + "\t" + failed + "\t" + passedAsPairs + "\n");
        		fw.close();
        	}
        	
        } catch (Exception e) {
        	log.error(e);
        }

        return 0;
    }
    
    
    private void removeMate(SAMRecord first) {
		first.setMateUnmappedFlag(true);
		first.setMateAlignmentStart(0);
		first.setMateNegativeStrandFlag(false);
		first.setMateReferenceName("*");
    }
    
    
    private boolean acceptMappingQuality(SAMRecord record) {
    	return record.getMappingQuality() >= MINIMUM_MAPPING_QUALITY;
    }
    
}
