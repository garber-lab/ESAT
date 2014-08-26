package broad.pda.seq.rap.rna;

import java.io.*;

import org.apache.log4j.Logger;

import net.sf.picard.PicardException;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.util.RExecutor;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;
import nextgen.core.alignment.*;

import nextgen.core.annotation.*;
import nextgen.core.annotation.filter.FullyContainedFilter;
import nextgen.core.exception.RuntimeIOException;
import nextgen.core.model.AlignmentModel;

public class CollectRnaSeqMetrics extends net.sf.picard.analysis.CollectRnaSeqMetrics {
	static Logger logger = Logger.getLogger(AlignmentModel.class.getName());
	
    public final String USAGE = getStandardUsagePreamble() +
            "Program to collect metrics about the alignment of RNA to various functional classes of loci in the genome:" +
            " coding, intronic, UTR, intergenic, ribosomal.\n" +
            "Also determines strand-specificity for strand-specific libraries.\n" +
            "For RAP, provide a EXCLUDE_REGION specifying the bounds of the target RNA.  Metrics will be calculated excluding this region.";


    @Option(doc="BED file containing regions to exclude from the calculations.  Reads will be counted in this area and output to a separate file.")
    public File EXCLUDE_BED=null;
	
	@Option(doc="File to write count of reads for excluded region")
	public File EXCLUDE_REGION_FILE=null;
	
	@Option(doc="Buffer around excluded region to also exclude")
	public Integer BUFFER = 1000;

	
	private AnnotationList<Annotation> excludedRegions = null;
	private int excludedCount = 0;
	private int excludedCountAntisense = 0;
	
    /** Required main method implementation. */
    public static void main(final String[] argv) {
        new CollectRnaSeqMetrics().instanceMainWithExit(argv);
    }

  
    @Override
    protected void setup(final SAMFileHeader header, final File samFile) {
    	super.setup(header, samFile);
    	if (EXCLUDE_BED != null) {
    		try {
    			excludedRegions = AnnotationFileReader.load(EXCLUDE_BED, Annotation.class, new BasicAnnotation.Factory());
    		} catch (IOException e) {
    			throw new RuntimeIOException(e);
    		}
    	}
    }
    
    @Override
    protected void acceptRead(final SAMRecord rec, final ReferenceSequence refSeq) {
    	Alignment read = new SingleEndAlignment(rec);
    	Annotation closest = excludedRegions.getClosest(read);
    	if (closest != null && closest.overlaps(read, BUFFER)) {
    		if (closest.getStrand() == read.getStrand() && !rec.getFirstOfPairFlag() ||
    				closest.getStrand() != read.getStrand() && rec.getFirstOfPairFlag()) {
    			++excludedCount;
    		} else {
    			++excludedCountAntisense;
    		}
    	} else {
    		super.acceptRead(rec, refSeq);
    	}
    }
    
    @Override
    protected void finish() {
    	try {
    		super.finish();
    	} catch (IllegalArgumentException e) {
    		// Script will fail to find the R file ... skip
    		logger.warn("Cannot create transcript coverage file ... call the original Picard JAR to do that.");
    	}
    	if (EXCLUDE_REGION_FILE == null && excludedRegions != null) {
    		System.out.println("Found " + (excludedCount + excludedCountAntisense) + " reads mapping to given region; these were excluded from the analysis.");
    	} else if (EXCLUDE_REGION_FILE != null) {
    		try {
    			FileWriter writer = new FileWriter(EXCLUDE_REGION_FILE);
    			writer.write("READS_MAPPING_TO_EXCLUDED_REGION\tREADS_MAPPING_TO_EXCLUDED_REGION_ANTISENSE\n" + excludedCount + "\t" + excludedCountAntisense + "\n");
    			writer.close();
    		} catch (IOException e) {
    			throw new RuntimeIOException(e);
    		}
    		
    	}
    }
}
