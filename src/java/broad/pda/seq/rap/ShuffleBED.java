package broad.pda.seq.rap;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

import broad.core.annotation.ShortBED;
import broad.core.annotation.ShortBEDReader;
import nextgen.core.annotation.AnnotationList;
import nextgen.core.coordinatesystem.GenomicSpace;

import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.util.Log;
import net.sf.samtools.util.CloseableIterator;
import nextgen.core.annotation.Annotation;

public class ShuffleBED extends GenomeCommandLineProgram {
    private static final Log log = Log.getInstance(ShuffleBED.class);

    @Usage
    public String USAGE = "Permute the regions in a BED file.  Use MASK_FILE and PCT_MASKED_ALLOWED to adjust behavior at masked regions.";
   
    @Option(doc="Input BED file containing query regions", shortName="I")
    public File INPUT;

	@Option(doc="Output file", shortName="O")
	public File OUTPUT;
	
	@Option(doc="Number of permutations", shortName="N")
	public Integer PERMUTATIONS = 1;
	
	
	/**
	 * Stock main method.
	 *
	 * @param args main arguments
	 */
	public static void main(final String[] args) {
		System.exit(new ShuffleBED().instanceMain(args));
	}
	
	

	@Override
	protected int doWork() {
		
		try {
			if (OUTPUT.exists()) {
				OUTPUT.delete();
			}
			
			AnnotationList<Annotation> regions = getRegionSet();
			ShortBEDReader peakReader = new ShortBEDReader(INPUT.getAbsolutePath());
			AnnotationList<ShortBED> peaks = new AnnotationList<ShortBED>(getCoordinateSpace(), peakReader.getAnnotationList()).getOverlappingAnnotationList(regions, true);
			
			// Cast to genomic space
			shuffleAndWriteAnnotations(regions, (GenomicSpace) getCoordinateSpace(), peaks, OUTPUT, PERMUTATIONS);
		} catch (Exception e) {
			log.error(e);
		}
		
		return 0;
	}	

	
	public static void shuffleAndWriteAnnotations(AnnotationList<Annotation> regions, GenomicSpace cs, AnnotationList<? extends Annotation> annotations, File output, int permutations) throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(output));
		
		// Permute multiple times and write all to the same file
		for (int i = 1; i <= permutations; i++) {
			// Shuffle the annotations within each region separately, but write all together to the permuted BED file
			for (Annotation region : regions) {
				CloseableIterator<? extends Annotation> overlappers = annotations.getOverlappingAnnotations(region, true);
				while (overlappers.hasNext()) {
					Annotation shuffled = overlappers.next().copy();
					cs.permuteAnnotation(shuffled, region);
					bw.write(shuffled.toShortBED());
					bw.newLine();
				}
			}
		}			
		bw.close();
	}
	
	
	
}
