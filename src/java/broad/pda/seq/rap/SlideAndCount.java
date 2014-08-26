package broad.pda.seq.rap;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.Iterator;
import java.util.List;

import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.util.Log;

import nextgen.core.annotation.*;
import nextgen.core.model.*;
import nextgen.core.model.score.*;


public class SlideAndCount extends GenomeScoringProgram {
    private static final Log log = Log.getInstance(SlideAndCount.class);
	
    @Usage
    public String USAGE = "Slides across the genome and counts reads (or calculates ratios for two models).";
   
	@Option(doc="Window size")
	public Integer WINDOW = null;
	
	@Option(doc="Overlap between windows")
	public Integer OVERLAP = null;

	@Option(doc="Output file", shortName="O")
	public File OUTPUT;

	/**
	 * Stock main method.
	 *
	 * @param args main arguments
	 */
	public static void main(final String[] args) {
		System.exit(new SlideAndCount().instanceMain(args));
	}
	

	@Override
	protected int doWork() {
		
		try {
			if (OUTPUT.exists()) {
				OUTPUT.delete();
			}
			
			List<Annotation> regions = getRegions();
			//log.info("Regions:" + regions);
			
			BufferedWriter bw = new BufferedWriter(new FileWriter(OUTPUT,true));
			WindowProcessor<? extends WindowScore> processor = getWindowProcessor();
					
			for (Annotation region : regions) {
				log.info(region);
				log.info("Starting: " + region.toUCSC());
				
				Iterator<? extends Annotation> windowIterator = getCoordinateSpace().getWindowIterator(region, WINDOW, OVERLAP);
				WindowScoreIterator<? extends WindowScore> itr = new WindowScoreIterator(windowIterator, processor, region);
				
				while (itr.hasNext()) {
					WindowScore curr = itr.next();
					bw.write(curr.toString());
					bw.newLine();
				}
				itr.close();
				
			}
			
			bw.close();
			
		} catch (Exception e) {
			log.error(e);
		}
		
		return 0;
	}
	
}
