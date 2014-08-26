package broad.pda.seq.rap;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.Iterator;
import java.util.List;
import java.io.IOException;

import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.util.Log;
import net.sf.samtools.util.CloseableIterator;

import nextgen.core.annotation.*;
import nextgen.core.general.TabbedReader;
import nextgen.core.model.score.*;


public class SlideAndCalculate extends GenomeScoringProgram {
    private static final Log log = Log.getInstance(SlideAndCalculate.class);
	
    @Usage
    public String USAGE = "Slides across the genome and counts reads (or calculates ratios for two models).  This is currently written for Genomic Space analyses but should be modified to allow other coordinate systems.";
   
	@Option(doc="Window size")
	public int WINDOW;
	
	@Option(doc="Overlap between windows")
	public int OVERLAP;

	@Option(doc="Output file", shortName="O")
	public File OUTPUT;

	@Option(doc="Delta factor to add to RPKM before calculating ratios", shortName="RPKM")
	public Double RPKM_RATIO_OFFSET = RatioScore.RPKM_OFFSET;
	

	/**
	 * Stock main method.
	 *
	 * @param args main arguments
	 */
	public static void main(final String[] args) {
		System.exit(new SlideAndCalculate().instanceMain(args));
	}

	@Override
	protected int doWork() {
		
		try {
			if (OUTPUT.exists()) {
				OUTPUT.delete();
			}

			RatioScore.RPKM_OFFSET = RPKM_RATIO_OFFSET;
			
			BufferedWriter bw = new BufferedWriter(new FileWriter(OUTPUT,true));
			Iterator<? extends WindowScore> itr = getWindowScoreIterator();
			while (itr.hasNext()) {
				WindowScore curr = itr.next();
				
				// Check to see if this window should be masked out using the given mask settings,
				// which might be different from the mask settings in the input SlideAndCount BED file
				if (getCoordinateSpace().isValidWindow(curr.getAnnotation())) {
					bw.write(curr.toString());
					bw.newLine();
				}
			}
				
			bw.close();
			
		} catch (Exception e) {
			log.error(e);
		}
		
		return 0;
	}
	
	public Iterator<? extends WindowScore> getWindowScoreIterator() throws IOException {
		Iterator<? extends WindowScore> itr;
		CloseableIterator<CountScore> target = TabbedReader.read(TARGET, CountScore.class, new CountScore.Factory());
		CloseableIterator<CountScore> control = TabbedReader.read(CONTROL, CountScore.class, new CountScore.Factory());
		if (SCORE.equalsIgnoreCase("ratio")) {
			itr = new RatioScore.RatioScoreIterator(target, control);
		} else {
			throw new IllegalArgumentException("Could not find scoring class for " + SCORE);
		}
		return itr;
	}
	
}
