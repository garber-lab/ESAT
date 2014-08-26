package broad.pda.seq.rap;

import java.io.*;
import java.util.List;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.util.Log;
import net.sf.samtools.util.CloseableIterator;
import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.AnnotationFileReader;
import nextgen.core.annotation.AnnotationList;
import nextgen.core.annotation.BasicAnnotation;
import nextgen.core.annotation.filter.FullyContainedFilter;
import nextgen.core.model.score.WindowProcessor;
import nextgen.core.model.score.WindowScore;
import nextgen.core.model.score.WindowScoreIterator;

public class DistanceToBED extends CommandLineProgram {
    private static final Log log = Log.getInstance(DistanceToBED.class);
    
	@Option(doc="BED file containing query annotations.")
	public File QUERY;

	@Option(doc="BED file containing reference annotations (anchors).")
	public File REF;
	
	@Option(doc="Exclude reads that overlap the following BED file", optional=true)
	public File EXCLUDE_BED = null;
	
	@Option(doc="Default behavior is to throw out query annotations where the closest reference is not the same strand.  Set to false to allow this")
	public Boolean IGNORE_STRAND = false;
	
	@Option(doc="Default behavior is to match the same strand.  Set to true for matching opposite strand")
	public Boolean MATCH_OPPOSITE_STRAND = false;
	
	@Option(doc="Set to true if you want to calculate distances between single bases - otherwise neighboring bases score as 0 distance", optional=true)
	public Boolean POINT_MATH = false;
	
	@Option(doc="Output BED file")
	public File OUTPUT;

	/**
	 * Stock main method.
	 *
	 * @param args main arguments
	 */
	public static void main(final String[] args) {
		System.exit(new DistanceToBED().instanceMain(args));
	}
    
	
	@Override
	protected int doWork() {
		
		try {
			if (OUTPUT.exists()) {
				OUTPUT.delete();
			}

			if (EXCLUDE_BED != null) {
				AnnotationList<Annotation> exclude = AnnotationFileReader.load(EXCLUDE_BED, Annotation.class, new BasicAnnotation.Factory());
				throw new UnsupportedOperationException("EXCLUDE_BED option not implemented yet");
			}
			
			// Load references into memory
			AnnotationList<Annotation> ref = AnnotationFileReader.load(REF, Annotation.class, new BasicAnnotation.Factory());
			log.info("Loaded " + ref.size() + " reference annotations.");
			
			// Iterate through buffered query annotations to handle large files gracefully
			CloseableIterator<Annotation> itr = AnnotationFileReader.read(QUERY, Annotation.class, new BasicAnnotation.Factory());
			
			BufferedWriter bw = new BufferedWriter(new FileWriter(OUTPUT,true));
			
			while (itr.hasNext()) {
				Annotation curr = itr.next();
				Annotation closestRef = ref.getClosest(curr);
				
				if (closestRef == null) continue;
				if (!IGNORE_STRAND && !MATCH_OPPOSITE_STRAND && curr.getStrand() != closestRef.getStrand()) {
					// skip if the strands do not match
					continue;
				} else if (!IGNORE_STRAND && MATCH_OPPOSITE_STRAND && curr.getStrand() == closestRef.getStrand()) {
					continue;
				}
				
				int dist = 0;
				
				if (POINT_MATH) {
					dist = curr.getStart() - closestRef.getStart();
				} else {
					dist = curr.getDistanceTo(closestRef);
					if (curr.getEnd() < closestRef.getStart()) dist = -1 * dist;
				}
				
				if (curr.getStrand() == Annotation.Strand.NEGATIVE) dist = -1 * dist;
				
				bw.write(curr.toShortBED() + "\t" + dist + "\t" + closestRef.getName() + "\n"); 
			}
			
			itr.close();
			bw.close();
			
		} catch (Exception e) {
			log.error(e);
		}
		
		return 0;
	}
}