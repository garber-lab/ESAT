package broad.pda.seq.rap;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.FileOutputStream;
import java.io.FileInputStream;
import java.io.PrintStream;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.SortedMap;
import java.util.TreeMap;

import broad.core.annotation.ShortBEDReader;
import broad.pda.enrichment.EnrichmentUtils;

import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.util.Log;
import net.sf.samtools.util.CloseableIterator;
import nextgen.core.annotation.*;
import nextgen.core.annotation.filter.*;
import nextgen.core.coordinatesystem.CoordinateSpace;
import nextgen.core.coordinatesystem.GenomicSpace;
import nextgen.core.model.score.CountScore;

import broad.core.math.EmpiricalDistribution;

import org.apache.commons.math3.stat.inference.*;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.stat.descriptive.*;
import org.apache.commons.math3.stat.descriptive.moment.*;
import org.apache.commons.math3.stat.descriptive.rank.*;
import org.apache.commons.io.IOUtils;
import org.ggf.drmaa.DrmaaException;

import nextgen.core.pipeline.Job;
import nextgen.core.pipeline.JobUtils;
import nextgen.core.pipeline.LSFJob;

public class CollectAnnotationEnrichments extends GenomeCommandLineProgram {
    private static final Log log = Log.getInstance(RatioPermutationPeakCaller.class);
    
    @Usage
    public String USAGE = "Calculate overlap/count enrichments for regions in the query BED file.  " + 
                          "Significance is assessed by permuting the query BED file.  " + 
                          "Regions in the BED file can be scored against other BED files or against BAM files.";
   
    @Option(doc="Input BED file containing query regions", shortName="I")
    public File INPUT;
    
	@Option(doc="Files or directories containing genomic annotations (BED format) or data files (BAM format)")
	public List<File> ANNOTATIONS;

	@Option(doc="Name of the analysis")
	public String ANALYSIS_NAME;
	
	@Option(doc="Output file")
	public File OUTPUT;
	
	@Option(doc="Number of permutations")
	public Integer PERMUTATIONS = 100;
	
	@Option(doc="Keep intermediate files")
	public boolean KEEP_INTERMEDIATES = false;
	
	@Option(doc="Queue to submit jobs to.  If non-null, will process all submitted annotations in this job.", optional=true)
	public String QUEUE = null;
	
	@Option(doc="Whether to calculate permutations or not.")
	public boolean PERMUTE = true;
	
	@Option(doc="Annotation base name.")
	public String ANNOTATION_BASE = "";
	
	@Option(doc="File extensions for discrete annotations to score.  Defaults to 'bed'.", shortName="DISCRETE_EXT")
	public List<String> DISCRETE_EXTENSION = new ArrayList<String>();
	
	@Option(doc="File extensions for BAM files to score.", shortName="BAM_EXT")
	public List<String> BAM_EXTENSION = new ArrayList<String>();
	
	@Option(doc="File extensions for continuous annotations to score.  Defaults to 'bedgraph'.", shortName="CONTINUOUS_EXT")
	public List<String> CONTINUOUS_EXTENSION = new ArrayList<String>();
	
	@Option(doc="Whether to load BAM files as paired end.")
	public boolean PAIRED_END = true;
	
	
	/**
	 * Stock main method.
	 *
	 * @param args main arguments
	 */
	public static void main(final String[] args) {
		System.exit(new CollectAnnotationEnrichments().instanceMain(args));
	}
	
	@Override
	protected String[] customCommandLineValidation() {
		// Add custom validation here.
		if (DISCRETE_EXTENSION.size() == 0) {
			DISCRETE_EXTENSION.add("bed");
		}
		if (BAM_EXTENSION.size() == 0) {
			BAM_EXTENSION.add("bam");
		}
		if (CONTINUOUS_EXTENSION.size() == 0) {
			CONTINUOUS_EXTENSION.add("bedgraph");
		}
		return super.customCommandLineValidation();
	}
	
	
	@Override
	protected int doWork() {
		
		try {
			if (OUTPUT.exists()) {
				OUTPUT.delete();
			}
			
			AnnotationList<Annotation> regions = getRegionSet();
			//ShortBEDReader peakReader = new ShortBEDReader(INPUT.getAbsolutePath());
			//AnnotationList<Annotation> peaks = new AnnotationList<Annotation>(getCoordinateSpace(), peakReader.getOverlappers(regions.toList(), true));
			AnnotationList<Annotation> peaks = AnnotationFileReader.load(INPUT, Annotation.class, new BasicAnnotation.Factory(), getCoordinateSpace(), new FullyContainedFilter(regions));
			log.info("Loaded " + peaks.size() + " peaks.");
			
			if (PERMUTE) shufflePeaks(regions, peaks, OUTPUT.getAbsolutePath());
			
			SortedMap<String,File> discreteMap = new TreeMap<String,File>();
			SortedMap<String,File> bamMap = new TreeMap<String,File>();
			SortedMap<String,File> continuousMap = new TreeMap<String,File>();
			
			for (File file : ANNOTATIONS) {
				discreteMap.putAll(EnrichmentUtils.findAnnotationFiles(file, DISCRETE_EXTENSION.toArray(new String[DISCRETE_EXTENSION.size()]), ANNOTATION_BASE));
				bamMap.putAll(EnrichmentUtils.findAnnotationFiles(file, BAM_EXTENSION.toArray(new String[BAM_EXTENSION.size()]), ANNOTATION_BASE));
				continuousMap.putAll(EnrichmentUtils.findAnnotationFiles(file, CONTINUOUS_EXTENSION.toArray(new String[CONTINUOUS_EXTENSION.size()]), ANNOTATION_BASE));
			}
			
			SortedMap<String,File> annotationMap = new TreeMap<String,File>();
			annotationMap.putAll(discreteMap);
			annotationMap.putAll(bamMap);
			annotationMap.putAll(continuousMap);
			
			if (QUEUE != null & annotationMap.size() > 1) {
				// This is the head job:  Submit children, collect results, cleanup
				submitJobs(annotationMap);
			} else {
				// This is a processing job - do some calculations
				calculateDiscreteAnnotationEnrichments(regions, peaks, discreteMap);
				calculateBamEnrichments(regions, peaks, bamMap);
				//calculateContinuousEnrichments(regions, peaks, continuousMap);
			}
			
			if (annotationMap.size() > 1) {
				collectResults(annotationMap);
				if (!KEEP_INTERMEDIATES) cleanup(annotationMap);
			}
					
		} catch (Exception e) {
			log.error(e);
		}
		
		return 0;
	}	
	
	
	protected void shufflePeaks(AnnotationList<Annotation> regions, AnnotationList<Annotation> peaks, String output) throws IOException {
		for (int i = 1; i <= PERMUTATIONS; i++) {
			String currOutput = output + "_" + i + ".bed";
			// Casts to genomic space
			ShuffleBED.shuffleAndWriteAnnotations(regions, (GenomicSpace) getCoordinateSpace(), peaks, new File(currOutput), 1);
		}
	}
	
	
	protected void submitJobs(SortedMap<String,File> annotationMap) throws IOException, InterruptedException, DrmaaException {
		Runtime run = Runtime.getRuntime();
		String jobID = LSFJob.generateJobID();
		Collection<Job> jobs = new ArrayList<Job>();
		
		for (String key : annotationMap.keySet()) {
			String command = "-M 4 -P RAP java -Xmx2g -cp " +
							 CollectAnnotationEnrichments.class.getProtectionDomain().getCodeSource().getLocation().getPath() +
							 " broad.pda.seq.rap.CollectAnnotationEnrichments";
			String[] argv = getCommandLineParser().getArgv();
			for (String arg : argv) {
				if (arg.indexOf("QUEUE") < 0 && arg.indexOf("PERMUTE") < 0 && arg.indexOf("ANNOTATIONS") < 0 && arg.indexOf("ANNOTATION_BASE") < 0) 
					command += " " + arg;
			}
			
			command += " PERMUTE=false QUEUE=null ANNOTATIONS=" + annotationMap.get(key).getAbsolutePath();
			command += " ANNOTATION_BASE=" + key.substring(0,key.lastIndexOf('.')+1);
			String bsub = OUTPUT.getAbsolutePath() + "." + key + ".bsub";
			LSFJob job = new LSFJob(run, jobID, command, bsub, QUEUE);
			jobs.add(job);
			job.submit();
		}
		
		JobUtils.waitForAll(jobs);
		if(!JobUtils.allSucceeded(jobs)) {
			log.error("One or more jobs failed. Proceeding anyway.");
		}
		
	}
	
	
	protected void calculateDiscreteAnnotationEnrichments(AnnotationList<Annotation> regions, AnnotationList<Annotation> peaks, SortedMap<String,File> discreteAnnotations) throws IOException {
		for (String key : discreteAnnotations.keySet()) {
			try {
				//ShortBEDReader annotationReader = new ShortBEDReader(discreteAnnotations.get(key).getAbsolutePath());
				//AnnotationList<Annotation> annotations = new AnnotationList<Annotation>(getCoordinateSpace(), annotationReader.getOverlappers(regions.toList(), true));
				AnnotationList<Annotation> annotations = AnnotationFileReader.load(discreteAnnotations.get(key), Annotation.class, new BasicAnnotation.Factory(), getCoordinateSpace(), new OverlapFilter(regions));
				processOneAnnotationSet(key, regions, peaks, annotations);
			} catch (RuntimeException e) {
				log.error("Skipping " + key + " due to error: " + e);
				e.printStackTrace();
			}
		}
	}
	
	
	protected void calculateBamEnrichments(AnnotationList<Annotation> regions, AnnotationList<Annotation> peaks, SortedMap<String,File> bamMap) throws IOException {
		for (String key : bamMap.keySet()) {
			try {
				AnnotationCollection<? extends Annotation> model = loadAlignmentModel(bamMap.get(key));
				processOneAnnotationSet(key, regions, peaks, model);
			} catch (RuntimeException e) {
				log.error("Skipping " + key + " due to error: " + e);
				e.printStackTrace();
			}
		}
	}
	
	
	protected void processOneAnnotationSet(String name, AnnotationList<Annotation> regions, AnnotationList<Annotation> peaks, AnnotationCollection<? extends Annotation> annotations) throws IOException {
		log.info("Processing " + name);
		List<CountScore> nullScores = new ArrayList<CountScore>();
		for (int i = 1; i <= PERMUTATIONS; i++) {
			ShortBEDReader shuffledReader = new ShortBEDReader(OUTPUT.getAbsolutePath() + "_" + i + ".bed");
			AnnotationList<Annotation> shuffled = new AnnotationList<Annotation>(getCoordinateSpace(), shuffledReader.getAnnotationList());
			nullScores.addAll(getPeakScores(regions, shuffled, annotations));
		}
		List<CountScore> peakScores = getPeakScores(regions, peaks, annotations);
		
		FileWriter writer = new FileWriter(new File(OUTPUT.getAbsolutePath() + "." + name + ".txt"));
		writer.write(new OverlapEnrichmentScore(ANALYSIS_NAME, name, peakScores, nullScores).toString());
		writer.write("\n");
		writer.close();
		
		// Write peak and null distributions for debugging or plotting purposes
		BufferedWriter peakWriter = new BufferedWriter(new FileWriter(new File(OUTPUT.getAbsolutePath() + "." + name + ".peaks.bed")));
		for (CountScore score : peakScores) {
			peakWriter.write(score.getAnnotation().toShortBED() + "\t" + (int) score.getCount() + "\n");
		}
		peakWriter.close();
		
		BufferedWriter nullWriter = new BufferedWriter(new FileWriter(new File(OUTPUT.getAbsolutePath() + "." + name + ".null.bed")));
		for (CountScore score : nullScores) {
			nullWriter.write(score.getAnnotation().toShortBED() + "\t" + (int) score.getCount() + "\n");
		}
		nullWriter.close();
	}

	
	protected void collectResults(SortedMap<String,File> discreteAnnotations) throws IOException {
		FileOutputStream output = new FileOutputStream(new File(OUTPUT.getAbsolutePath() + ".txt"));
		PrintStream print = new PrintStream(output);
		print.print(OverlapEnrichmentScore.getColumnNames() + "\n");
		for (String key : discreteAnnotations.keySet()) {
			File file = new File(OUTPUT.getAbsolutePath() + "." + key + ".txt");
			System.out.print("Loading " + file.getPath() + "...");
			if (file.exists()) {
				System.out.println("success.");
				FileInputStream input = new FileInputStream(file);
				IOUtils.copy(input, output);
				input.close();
			} else {
				System.out.println("failure.");
			}
		}
		output.close();
	}
	
	
	protected void cleanup(SortedMap<String,File> discreteAnnotations) {
		for (String key : discreteAnnotations.keySet()) {
			File file = new File(OUTPUT.getAbsolutePath() + "." + key + ".txt");
			if (file.exists()) file.delete();
		}
		for (int i = 1; i <= PERMUTATIONS; i++) {
			File file = new File(OUTPUT.getAbsolutePath() + "_" + i + ".bed");
			if (file.exists()) file.delete();
		}
	}
	
	protected List<CountScore> getPeakScores(AnnotationList<Annotation> regions, AnnotationList<? extends Annotation> peaks, AnnotationCollection<? extends Annotation> annotations) {
		List<CountScore> scores = new ArrayList<CountScore>();

		double total = annotations.getGlobalCount();
		double regionTotal = (double) annotations.getCount(regions);
		CloseableIterator<? extends Annotation> itr = peaks.getOverlappingAnnotationList(regions).iterator();
		while (itr.hasNext()) {
			Annotation peak = itr.next();
			double score = annotations.getCount(peak);
			scores.add(new CountScore(peak, score, regionTotal, total));
		}

		return scores;
		
	}
	
	
	public static class OverlapEnrichmentScore {
		String name, analysisName;
		private int totalOverlap = 0, regionTotal, total;
		private List<Double> peakDensities = new ArrayList<Double>();
		private List<Double> nullDensities = new ArrayList<Double>();
		private double[] peakDensitiesPrimitive, nullDensitiesPrimitive;
		
		public OverlapEnrichmentScore(String analysisName, String name, List<CountScore> peakScores, List<CountScore> nullScores) {
			this.name = name;
			this.analysisName = analysisName;
			regionTotal = (int) peakScores.get(0).getRegionTotal();
			total = (int) peakScores.get(0).getTotal();
			
			for (CountScore peakScore : peakScores) {
				peakDensities.add(peakScore.getCount() / peakScore.getAnnotation().length());
				totalOverlap += peakScore.getCount();
			}
			
			for (CountScore nullScore : nullScores) {
				nullDensities.add(nullScore.getCount() / nullScore.getAnnotation().length());
			}
			
			peakDensitiesPrimitive = ArrayUtils.toPrimitive(peakDensities.toArray(new Double[peakDensities.size()]));
			nullDensitiesPrimitive = ArrayUtils.toPrimitive(nullDensities.toArray(new Double[nullDensities.size()]));
		}
		
		public double getMeanPeakDensity() {
			UnivariateStatistic stat = new Mean();
			return stat.evaluate(peakDensitiesPrimitive);
		}
		
		public double getMeanNullDensity() {
			UnivariateStatistic stat = new Mean();
			return stat.evaluate(nullDensitiesPrimitive);
		}
		
		public double getSdPeakDensity() {
			UnivariateStatistic stat = new StandardDeviation();
			return stat.evaluate(peakDensitiesPrimitive);
		}
		
		public double getSdNullDensity() {
			UnivariateStatistic stat = new StandardDeviation();
			return stat.evaluate(nullDensitiesPrimitive);
		}
		
		public double getMedianPeakDensity() {
			UnivariateStatistic stat = new Median();
			return stat.evaluate(peakDensitiesPrimitive);
		}
		
		public double getMedianNullDensity() {
			UnivariateStatistic stat = new Median();
			return stat.evaluate(nullDensitiesPrimitive);
		}
		
		public double getTTestPvalue() {
			return TestUtils.tTest(peakDensitiesPrimitive, nullDensitiesPrimitive);
		}
		
		public double getMannWhitneyPValue() {
			// NOTE please update your lib directory and add the new commons-math JAR file
			MannWhitneyUTest stat = new MannWhitneyUTest();
			return stat.mannWhitneyUTest(peakDensitiesPrimitive, nullDensitiesPrimitive);
		}
		
		public double getKLDivergence() {
			List<Double> all = new ArrayList<Double>();
			all.addAll(peakDensities);
			all.addAll(nullDensities);
			EmpiricalDistribution ed = new EmpiricalDistribution(all);
			EmpiricalDistribution peakEd = new EmpiricalDistribution(peakDensitiesPrimitive, ed.getBinNumber(), ed.getMin(), ed.getMax());
			EmpiricalDistribution nullEd = new EmpiricalDistribution(nullDensitiesPrimitive, ed.getBinNumber(), ed.getMin(), ed.getMax());
			return peakEd.KLDivergence(nullEd);
		}
		
		
		public String toString() {
			return analysisName + "\t" + name + "\t" + (getMeanPeakDensity()/getMeanNullDensity()) + "\t" + getMeanPeakDensity() + "\t" + getSdPeakDensity() + "\t" + getMeanNullDensity() + "\t" + getSdNullDensity() + "\t" + getMedianPeakDensity() + "\t" + getMedianNullDensity() + "\t" +
			       getTTestPvalue() + "\t" + getMannWhitneyPValue() + "\t" + getKLDivergence() + "\t" + regionTotal + "\t" + total;
		}
		
		public static String getColumnNames() {
			return "analysis\tname\tdensityRatio\tmeanPeakDensity\tsdPeakDensity\tmeanNullDensity\tsdNullDensity\tmedianPeakDensity\tmedianNullDensity\tT-Test\tMannWhitneyTest\tKLDivergence\tregionTotal\ttotal";
		}
	}
	
}
