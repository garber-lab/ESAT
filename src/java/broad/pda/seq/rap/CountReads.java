package broad.pda.seq.rap;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.List;

import broad.core.annotation.ShortBEDReader;

import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.util.Log;
import nextgen.core.annotation.*;
import nextgen.core.annotation.filter.FullyContainedFilter;
import nextgen.core.model.AlignmentModel;
import nextgen.core.model.score.*;

// test comment
public class CountReads extends GenomeScoringProgram {
    private static final Log log = Log.getInstance(CountReads.class);

    @Usage
    public String USAGE = "Counts reads overlapping the provided BED files.  This is currently written for Genomic Space analyses but should be modified to allow other coordinate systems.";
    
    @Option(doc="Target BED file")
    public File ANNOTATION_FILE;
    
    @Option(doc="Output BED file", shortName="O")
    public File OUTPUT;
    
	/**
	 * Stock main method.
	 *
	 * @param args main arguments
	 */
	public static void main(final String[] args) {
		System.exit(new CountReads().instanceMain(args));
	}
    
	
	@Override
	protected int doWork() {
		
		try {
			if (OUTPUT.exists()) {
				OUTPUT.delete();
			}

			List<Annotation> regions = getRegions();
			AnnotationList<Annotation> targets = AnnotationFileReader.load(ANNOTATION_FILE, Annotation.class, new BasicAnnotation.Factory(), getCoordinateSpace(), new FullyContainedFilter(regions));
			log.info("Loaded " + targets.size() + " annotations.");
			
			BufferedWriter bw = new BufferedWriter(new FileWriter(OUTPUT,true));
			WindowProcessor<? extends WindowScore> processor = getWindowProcessor();
			
			for (Annotation region : regions) {
				log.info("Starting: " + region.toUCSC());
				
				WindowScoreIterator<? extends WindowScore> itr = new WindowScoreIterator(targets.getOverlappingAnnotations(region), processor, region);
				while (itr.hasNext()) {
					bw.write(itr.next().toString());
					bw.newLine();
				}
			}
			
			bw.close();
			
		} catch (Exception e) {
			log.error(e);
		}
		
		return 0;
	}
}
