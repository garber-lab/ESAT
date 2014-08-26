package broad.pda.seq.rap;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Iterator;
import java.util.List;

import broad.core.annotation.ShortBEDReader;

import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.Log;
import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.AnnotationFileReader;
import nextgen.core.annotation.AnnotationList;
import nextgen.core.annotation.BasicAnnotation;
import nextgen.core.annotation.filter.FullyContainedFilter;
import nextgen.core.coordinatesystem.*;
import nextgen.core.model.score.WindowProcessor;
import nextgen.core.model.score.WindowScore;
import nextgen.core.model.score.WindowScoreIterator;

public class PlotAggregateRegions extends GenomeScoringProgram {
    private static final Log log = Log.getInstance(PlotAggregateRegions.class);
    
    @Usage
    public String USAGE = "Plot scores as as a function of position from / in a target.  Scans through introns in the annotation file";
    
    @Option(doc="Output file", shortName="O")
    public File OUTPUT;
    
	@Option(doc="File containing genomic annotations (bed or bedGraph format)")
	public File ANNOTATION_FILE;
	
	@Option(doc="Window size")
	public int WINDOW;

	@Option(doc="Overlap between windows", optional=true)
	public int OVERLAP=0;
	
	@Option(doc="Length of the inner region", optional=true)
	public int INNER_LENGTH = 0;
	
	@Option(doc="Length of the outer region", optional=true)
	public int OUTER_LENGTH = 0;
	
	@Option(doc="Length of the middle region", optional=true)
	public int MIDDLE_LENGTH = 0;
	
	@Option(doc="Set to true to generate one-sided diagrams")
	public boolean SYMMETRIC=false;
	
	@Option(doc="Eliminate windows in the inner and middle subregions that are outside the boundaries of the annotation")
	public boolean CLIP_AT_BOUNDARIES=true;

	
	
	private int regionCounter = 0;
	
	/**
	 * Stock main method.
	 *
	 * @param args main arguments
	 */
	public static void main(final String[] args) {
		System.exit(new PlotAggregateRegions().instanceMain(args));
	}
	
	protected void loadCoordinateSpace() {
		// load coordinate space without masking; apply masking manually later
		coordinateSpace = new GenomicSpace(SIZES);
	}
	
	@Override
	protected int doWork() {
		try {

			List<Annotation> regions = getRegions();
			AnnotationList<Annotation> targets = AnnotationFileReader.load(ANNOTATION_FILE, Annotation.class, new BasicAnnotation.Factory(), getCoordinateSpace(), new FullyContainedFilter(regions));
			log.info("Loaded " + targets.size() + " annotations.");
			
			BufferedWriter bw = new BufferedWriter(new FileWriter(OUTPUT));
			WindowProcessor<? extends WindowScore> processor = getWindowProcessor();
			GenomicSpace maskedSpace = new GenomicSpace(SIZES, MASK_FILE, PCT_MASKED_ALLOWED);
			
			for (Annotation region : regions) {
				log.info("Starting: " + region.toUCSC());
				
				Iterator<Annotation> itr = targets.getOverlappingAnnotations(region);
				while (itr.hasNext()) {
					
					// Use the genomic coordinate space to fill in the target (for Genomic Space)
					Annotation target = itr.next();
					Annotation fragment = maskedSpace.getFragment(target).iterator().next();
					fragment.setOrientation(target.getOrientation());
					fragment.setName(target.getName());
					target = fragment;
					
					PlotRegions subregions = generateSubregions(target);
					regionCounter += 1;

					Integer counter = 0;
					if (SYMMETRIC) {
						counter = scanAndPrint(target, subregions.beginOuter, "outer", processor, maskedSpace, bw, counter, false);
						counter = scanAndPrint(target, subregions.beginInner, "inner", processor, maskedSpace, bw, counter, true);
						counter = scanAndPrint(target, subregions.middle, "middle", processor, maskedSpace, bw, counter, true);
						counter = 0; // reset counter since the end is symmetric to the beginning
						counter = scanAndPrint(target, subregions.endOuter, "outer", processor, maskedSpace, bw, counter, false);
						counter = scanAndPrint(target, subregions.endInner, "inner", processor, maskedSpace, bw, counter, true);
					} else {
						counter = scanAndPrint(target, subregions.beginOuter, "beginOuter", processor, maskedSpace, bw, counter, false);
						counter = scanAndPrint(target, subregions.beginInner, "beginInner", processor, maskedSpace, bw, counter, true);
						counter = scanAndPrint(target, subregions.middle, "middle", processor, maskedSpace, bw, counter, true);
						counter = scanAndPrint(target, subregions.endInner, "endInner", processor, maskedSpace, bw, counter, true);
						counter = scanAndPrint(target, subregions.endOuter, "endOuter", processor, maskedSpace, bw, counter, false);
					}
				}
			}
			
			bw.close();
			
			
		} catch (Exception e) {
			log.error(e);
			System.exit(1);
		}
		
		return 0;
	}
	
	/**
	 * Generate subregions to scan based on the given region.
	 * @param region
	 * @return
	 * @throws IOException
	 */
	public PlotRegions generateSubregions(Annotation region) throws IOException {
		PlotRegions subregions;
		
		Annotation outerLeft = new BasicAnnotation(region.getReferenceName(), region.getStart() - OUTER_LENGTH, region.getStart(), region.getOrientation());
		Annotation outerRight = new BasicAnnotation(region.getReferenceName(), region.getEnd(), region.getEnd() + OUTER_LENGTH, region.getOrientation());
		Annotation innerLeft = new BasicAnnotation(region.getReferenceName(), region.getStart(), region.getStart() + INNER_LENGTH, region.getOrientation());
		Annotation innerRight = new BasicAnnotation(region.getReferenceName(), region.getEnd() - INNER_LENGTH, region.getEnd(), region.getOrientation());
		// TODO:  Adjust so that these regions don't extend beyond the length of a short gene
		
		int midway = (int) Math.floor((region.getEnd() + region.getStart())/2.0);
		int midAdjust = (int) Math.floor(MIDDLE_LENGTH/2.0);
		Annotation middle = new BasicAnnotation(region.getReferenceName(), midway - midAdjust, midway - midAdjust + MIDDLE_LENGTH, region.getOrientation());
		
		if (region.isNegativeStrand()) {
			subregions = new PlotRegions(outerRight, innerRight, middle, innerLeft, outerLeft);
		} else {
			subregions = new PlotRegions(outerLeft, innerLeft, middle, innerRight, outerRight);
		} 
		
		if (SYMMETRIC) {
			reverseOrientation(subregions.endOuter);
			reverseOrientation(subregions.endInner);
		}
		
		return subregions;
	}
	
	
	private void reverseOrientation(Annotation a) {
		if (a.isNegativeStrand()) {
			a.setOrientation(Annotation.Strand.POSITIVE);
		} else {
			a.setOrientation(Annotation.Strand.NEGATIVE);
		}
	}
	
	
	/**
	 * Scan a subregion and print to a file
	 * @param subregion
	 * @param name
	 * @param processor
	 * @param maskedSpace
	 * @param bw
	 * @param counter Keeps track of the relative base from the beginning of the subregions.
	 * @throws IOException
	 */
	private int scanAndPrint(final Annotation target, Annotation subregion, final String name, WindowProcessor<? extends WindowScore> processor, final GenomicSpace maskedSpace, BufferedWriter bw, Integer counter, boolean clip) throws IOException {
		//log.info("starting " + name + " " + subregion.toUCSC());
		Iterator<? extends Annotation> windowIterator = ((GenomicSpace) getCoordinateSpace()).getWindowIterator(subregion, WINDOW, OVERLAP, true);
		WindowScoreIterator<? extends WindowScore> itr = new WindowScoreIterator(windowIterator, processor, subregion);
		while (itr.hasNext()) {
			
			try {
				WindowScore ws = itr.next();

				// Allow window if it is not masked
				if (maskedSpace.isValidWindow(ws.getAnnotation())) {


					// Accept windows only if they do not pass outside boundaries of the target
					if (!clip || !CLIP_AT_BOUNDARIES || target.contains(ws.getAnnotation())) {

						String regionName = target.getName();
						if (regionName.equals("")) regionName = regionCounter + "";
						// Write the subregion name, location, and the result
						bw.write(regionName + "\t" + name + "\t" + counter + "\t" + ws.toString() + "\n");
					}

				}
			} catch (AnnotationOutOfBoundsException e) {
				// this is okay .. just skip this window
			}
			
			// Advance counter even if window is skipped
			counter += (WINDOW - OVERLAP);
		}
		return counter;
	}
	
	
	private class PlotRegions {
		public Annotation beginOuter, beginInner, middle, endInner, endOuter;
		public PlotRegions(Annotation beginOuter, Annotation beginInner, Annotation middle, Annotation endInner, Annotation endOuter) {
			this.beginOuter = beginOuter;
			this.beginInner = beginInner;
			this.middle = middle;
			this.endInner = endInner;
			this.endOuter = endOuter;
		}
	}
}
