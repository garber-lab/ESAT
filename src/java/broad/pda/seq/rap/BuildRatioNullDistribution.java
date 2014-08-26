package broad.pda.seq.rap;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.ggf.drmaa.DrmaaException;

import nextgen.core.pipeline.Job;
import nextgen.core.pipeline.JobUtils;
import nextgen.core.pipeline.LSFJob;

import broad.core.math.EmpiricalDistribution;

import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.util.Log;
import net.sf.samtools.util.CloseableIterator;
import nextgen.core.coordinatesystem.GenomicSpace;
import nextgen.core.model.AlignmentModel;
import nextgen.core.model.score.WindowScoreIterator;
import nextgen.core.alignment.Alignment;
import nextgen.core.model.score.*;
import nextgen.core.annotation.Annotation;
import nextgen.core.writers.PairedEndWriter;
import nextgen.core.readFilters.*;

public class BuildRatioNullDistribution extends GenomeCommandLineProgram {
    private static final Log log = Log.getInstance(BuildRatioNullDistribution.class);
	
    @Usage
    public String USAGE = "Builds a null distribution of window ratios between two SAM or BAM files";
   
	@Option(doc="Window size")
	public int WINDOW;
	
	@Option(doc="Overlap between windows")
	public int OVERLAP;

	@Option(doc="Input SAM or BAM file", shortName="I")
	public File TARGET;

	@Option(doc="Control SAM or BAM file for normalization.")
	public File CONTROL;
	
	@Option(doc="Output file basename")
	public File OUTPUT;
	
	@Option(doc="Number of permutations", shortName="N")
	public Integer PERMUTATIONS = 100;
	
	@Option(doc="Keep intermediate BAM and empirical distribution files for debugging purposes")
	public boolean KEEP_INTERMEDIATES = false;
	
	@Option(doc="Use existing intermediate BAM and empirical distribution files")
	public boolean USE_INTERMEDIATES = false;
	
	@Option(doc="Queue to submit jobs to.  If null, will simply process with one permutation", optional=true)
	public String QUEUE = null;
	
	@Option(doc="Set the random generator seed for reproducible behavior (note: does not work for submitting multiple jobs).", optional=true)
	public Integer RANDOM_SEED = null;
	
	
	private Runtime run = Runtime.getRuntime();
	private String jobID = LSFJob.generateJobID();

	/**
	 * Stock main method.
	 *
	 * @param args main arguments
	 */
	public static void main(final String[] args) {
		System.exit(new BuildRatioNullDistribution().instanceMain(args));
	}
	

	@Override
	protected int doWork() {
		
		try {
			
			if (QUEUE != null) {
				// Load INPUT and CONTROL files to set up pebam and GlobalStats to avoid concurrent writing in the submitted jobs
				AlignmentModel input = loadAlignmentModel(TARGET);
				AlignmentModel control = loadAlignmentModel(CONTROL);
				control.getGlobalCount();  // init global stats
				
				submitJobs();
				collateResults();
				if (!KEEP_INTERMEDIATES) cleanup();
			} else {
				permuteOnce(OUTPUT.getAbsolutePath());
			}
			
		} catch (Exception e) {
			log.error(e);
		}
		
		return 0;
	}

	
	
	/**
	 * Submit one job per permutation
	 * @throws IOException
	 * @throws InterruptedException
	 * @throws DrmaaException 
	 */
	private void submitJobs() throws IOException, InterruptedException, DrmaaException {
		Collection<Job> jobs = new ArrayList<Job>();
		for (int i = 0; i < PERMUTATIONS; i++) {

			String command = "-M 4 -P RAP java -Xmx4g -cp " + 
							 BuildRatioNullDistribution.class.getProtectionDomain().getCodeSource().getLocation().getPath() +
							 " broad.pda.seq.rap.BuildRatioNullDistribution";
			
			String[] argv = getCommandLineParser().getArgv();
			for (String arg : argv) {
				if (arg.indexOf("QUEUE") < 0 && arg.indexOf("OUTPUT") < 0) // don't set QUEUE
					command += " " + arg;
			}
			
			String outfile = OUTPUT.getAbsolutePath() + "_" + i;
			command += " OUTPUT=" + outfile;
			
						//	 " WINDOW=" + WINDOW + " OVERLAP=" + OVERLAP + " INPUT=" + INPUT + " OUTPUT=" + outfile + " REGION=" + REGION + " MASK_FILE=" + MASK_FILE + " MAX_FRAGMENT_LENGTH=" + MAX_FRAGMENT_LENGTH + " SIZES=" + SIZES;
			String bsubOutput = outfile + ".bsub";
			LSFJob job = new LSFJob(run, jobID, command, bsubOutput, QUEUE);
			jobs.add(job);
			job.submit();
		}
		JobUtils.waitForAll(jobs);
	}

	
	/**
	 * Collect the permutation results into a single empirical distribution
	 * @throws IOException
	 */
	private void collateResults() throws IOException {
		EmpiricalDistribution ed = getEmptyEmpiricalDistribution();
		for (int i = 0; i < PERMUTATIONS; i++) {
			String file = OUTPUT.getAbsolutePath() + "_" + i + ".empiricalDistribution.txt";
			EmpiricalDistribution curr = new EmpiricalDistribution(new File(file));
			ed.addDistribution(curr);
		}
		ed.write(OUTPUT.getAbsolutePath() + ".empiricalDistribution.txt");
	}
	
	
	/**
	 * Delete intermediate / temporary files
	 * @throws IOException
	 */
	private void cleanup() throws IOException {
		for (int i = 0; i < PERMUTATIONS; i++) {
			String file = OUTPUT.getAbsolutePath() + "_" + i + ".empiricalDistribution.txt";
			new File(file).delete();
			file = OUTPUT.getAbsolutePath() + "_" + i + ".PairedEnd.bam";
			new File(file).delete();
			file = OUTPUT.getAbsolutePath() + "_" + i + ".PairedEnd.bam.bai";
			new File(file).delete();
			file = OUTPUT.getAbsolutePath() + "_" + i + ".PairedEnd.bam.GenomicSpaceStats";
			new File(file).delete();
			file = OUTPUT.getAbsolutePath() + "_" + i + ".bsub";
			new File(file).delete();
		}
	}
	
	/**
	 * Do the processing for a single permutation
	 * @param output
	 * @throws IOException
	 */
	private void permuteOnce(String output) throws IOException {
		List<Annotation> regions = getRegions();

		String permutedOutput = output + PairedEndWriter.PAIRED_END_EXTENSION;
		
		if (USE_INTERMEDIATES && new File(permutedOutput).exists()) {
			log.info("Using existing intermediate file: " + permutedOutput);
		} else {
			log.info("Setting up writer: " + permutedOutput);
			final PairedEndWriter outputWriter = new PairedEndWriter(TARGET, permutedOutput);

			// Reopen the input file in a model. Casts to genomic space.
			((GenomicSpace) getCoordinateSpace()).setPercentMaskedAllowed(0.0);  // don't allow permuting reads to masked regions
			AlignmentModel model = loadAlignmentModel(TARGET);
			
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
		
		String edOutput = output + ".empiricalDistribution.txt";
		if (USE_INTERMEDIATES && new File(edOutput).exists()) {
			log.info("Using existing empirical distribution file: " + edOutput);
		} else {
			// Now scan over the permuted data and calculate ratios. Casts to genomic space.
			((GenomicSpace) getCoordinateSpace()).setPercentMaskedAllowed(PCT_MASKED_ALLOWED);  // reset pct masked allowed after changing above
			AlignmentModel permuted = loadAlignmentModel(new File(permutedOutput));
			AlignmentModel control = loadAlignmentModel(CONTROL);
			WindowProcessor<RatioScore> processor = new RatioScore.Processor(permuted, control);

			EmpiricalDistribution ed = getEmptyEmpiricalDistribution();

			for (Annotation region : regions) {
				WindowScoreIterator<RatioScore> windowItr = permuted.scan(region, WINDOW, OVERLAP, processor);
				while (windowItr.hasNext()) {
					RatioScore curr = windowItr.next();
					ed.add(curr.getLog2RegionRatio());
					if (curr.getLog2RegionRatio() < -20 || curr.getLog2RegionRatio() > 20) {
						log.warn("Ratio out of bounds: " + curr);
					}
				}
			}

			ed.write(output + ".empiricalDistribution.txt");
		}
	}
	
	
	/**
	 * Generate an empty empirical distribution.  Use this function so that all EDs in the class
	 * have a standard number of bins, etc.
	 * @return
	 */
	public static EmpiricalDistribution getEmptyEmpiricalDistribution() {
		return new EmpiricalDistribution(2000, -20, 20);
	}
	
}
