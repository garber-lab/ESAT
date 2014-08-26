package broad.pda.seq.rap;

import java.io.File;
import java.io.IOException;
import java.util.List;

import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.util.Log;
import net.sf.samtools.util.CloseableIterator;
import nextgen.core.alignment.Alignment;
import nextgen.core.annotation.Annotation;
import nextgen.core.coordinatesystem.GenomicSpace;
import nextgen.core.model.AlignmentModel;
import nextgen.core.writers.PairedEndWriter;

public class ShuffleBAM extends GenomeCommandLineProgram {
    private static final Log log = Log.getInstance(ShuffleBAM.class);

    @Usage
    public String USAGE = "Shuffle reads in a SAM or BAM file, avoiding unmappable regions.";
	
    @Option(doc="Input SAM or BAM file", shortName="I")
	public File INPUT;

	@Option(doc="Output file basename", shortName="O")
	public File OUTPUT;
	
	@Option(doc="Number of permutations", shortName="N")
	public Integer PERMUTATIONS = 1;
	
	@Option(doc="Set the random generator seed for reproducible behavior.", optional=true)
	public Integer RANDOM_SEED = null;
	

	/**
	 * Stock main method.
	 *
	 * @param args main arguments
	 */
	public static void main(final String[] args) {
		System.exit(new ShuffleBAM().instanceMain(args));
	}
	

	@Override
	protected int doWork() {
		
		try {
			for (int i = 0; i < PERMUTATIONS; i++) {
				String output = OUTPUT.getAbsolutePath();
				if (PERMUTATIONS != 1) {
					output = output + "_" + i;
				}
				permuteOnce(output);
			}
		} catch (Exception e) {
			log.error(e);
		}
		
		return 0;
	}
	
	
	private void permuteOnce(String output) throws IOException {
		List<Annotation> regions = getRegions();
		String permutedOutput = output + PairedEndWriter.PAIRED_END_EXTENSION;

		log.info("Setting up writer: " + permutedOutput);
		final PairedEndWriter outputWriter = new PairedEndWriter(INPUT, permutedOutput);

		// Reopen the input file in a model. Fails if not in genomic space.
		((GenomicSpace) getCoordinateSpace()).setPercentMaskedAllowed(0.0);  // don't allow permuting reads to masked regions
		AlignmentModel model = loadAlignmentModel(INPUT);

		if (RANDOM_SEED != null) GenomicSpace.setSeed(RANDOM_SEED);

		// Permute reads in the model and write to the output file
		for (Annotation region : regions) {
			CloseableIterator<Alignment> itr = model.getPermutedAnnotations(region);
			while (itr.hasNext()) {
				Alignment curr = itr.next();
				outputWriter.addAlignment(curr);
			}
			itr.close();
		}
		outputWriter.close(); // this creates a .bai index
	}
}
