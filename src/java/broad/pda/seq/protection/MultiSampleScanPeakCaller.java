/**
 * Call protein binding sites on RNA with Protect-seq data from multiple control and multiple signal samples
 */
package broad.pda.seq.protection;

import broad.core.parser.CommandLineParser;
import broad.core.parser.StringParser;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.ggf.drmaa.DrmaaException;

import nextgen.core.pipeline.Job;
import nextgen.core.pipeline.JobUtils;
import nextgen.core.pipeline.LSFJob;

import broad.core.math.MathUtil;
import broad.core.math.Statistics;
import broad.pda.annotation.BEDFileParser;

import nextgen.core.alignment.Alignment;
import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Gene;
import nextgen.core.annotation.Annotation.Strand;
import nextgen.core.coordinatesystem.CoordinateSpace;
import nextgen.core.coordinatesystem.TranscriptomeSpace;
import nextgen.core.model.AlignmentModel;
import nextgen.core.model.TranscriptomeSpaceAlignmentModel;
import nextgen.core.model.score.BinomialEnrichmentScore;
import nextgen.core.model.score.ScanStatisticScore;
import nextgen.core.model.score.WindowProcessor;
import nextgen.core.model.score.WindowScoreIterator;
import nextgen.core.utils.AlignmentUtils;
import nextgen.core.utils.AnnotationUtils;
import nextgen.core.feature.GeneWindow;

/**
 * @author prussell
 *
 */
public class MultiSampleScanPeakCaller {
	
	private TranscriptomeSpace coord;
	protected Map<String, Collection<Gene>> genes;
	private Map<Gene, Map<Annotation, Double>> tStatisticWindowScores;
	private Map<SampleData, Map<Gene, Map<Annotation, Double>>> singleSampleWindowEnrichmentOverGene;
	private Map<SampleData, Map<Gene, Collection<Gene>>> singleSampleScanPeaks;
	protected Map<SampleData, Map<Gene, Map<Annotation, BinomialEnrichmentScore>>> binomialWindowScores;
	private Map<SampleData,WindowProcessor<BinomialEnrichmentScore>> binomialProcessors;
	protected GenomeSpaceSampleData expressionData;
	protected ArrayList<SampleData> controlSamples;
	protected ArrayList<SampleData> signalSamples;
	protected ArrayList<SampleData> allSamples;
	private boolean useBinomialScore;
	private SampleData binomialCtrl;
	protected static Logger logger = Logger.getLogger(MultiSampleScanPeakCaller.class.getName());
	protected int windowSize;
	protected int stepSize;
	private static int DEFAULT_WINDOW_SIZE = 20;
	private static int DEFAULT_STEP_SIZE = 1;
	private static double DEFAULT_PEAK_SCAN_P_VALUE_CUTOFF = 0.001;
	private static double DEFAULT_PEAK_WINDOW_COUNT_CUTOFF = 10;
	private static double DEFAULT_TRIM_PEAK_QUANTILE = 0.6;
	private static int DEFAULT_BATCH_MEM_REQUEST = 8;
	private static double DEFAULT_EXPRESSION_SCAN_P_VALUE_CUTOFF = 0.01;
	private static boolean DEFAULT_FIRST_READ_TRANSCRIPTION_STRAND = false;
	private static double DEFAULT_PEAK_MAX_PCT_DUPLICATES = 0.5;
	private static boolean DEFAULT_FILTER_BY_STRAND = true;
	private static boolean DEFAULT_EXTRA_FIELDS = false;
	private static boolean DEFAULT_USE_BINOMIAL = false;
	private boolean filterByStrand;
	private boolean extraFields;
	protected int numControls;
	protected int numSignals;
	protected int numSamples;
	protected Random random;
	private Collection<SamplePermutation> sampleIdentityPermutations;
	private double peakWindowScanPvalCutoff;
	private double peakWindowCountCutoff;
	private double trimQuantile;
	private boolean permutationScoring;
	private String sampleFile;
	private String bedAnnotationFile;
	protected String sizeFile;
	private boolean firstReadTranscriptionStrand;
	private double peakCutoffMaxReplicatePct;
	private static int RGB_RED_WITH_GENE = 106;
	private static int RGB_GREEN_WITH_GENE = 7;
	private static int RGB_BLUE_WITH_GENE = 205;
	private static int RGB_RED_AGAINST_GENE = 218;
	private static int RGB_GREEN_AGAINST_GENE = 2;
	private static int RGB_BLUE_AGAINST_GENE = 38;
	private static int RGB_RED_UNKNOWN = 62;
	private static int RGB_GREEN_UNKNOWN = 63;
	private static int RGB_BLUE_UNKNOWN = 73;
	private Map<SampleData, FileWriter> expressionFilterRejectWriters;
	private static String REJECT_FILE_PREFIX_EXPRESSION_FILTER = "expression_filter";
	private Map<SampleData, FileWriter> windowCountRejectFileWriters;
	private static String REJECT_FILE_PREFIX_WINDOW_COUNT = "window_count_filter";
	private Map<SampleData, FileWriter> windowScanPvalAllFragmentsRejectFileWriters;
	private static String REJECT_FILE_PREFIX_WINDOW_SCAN_PVAL_ALL_FRAGMENTS = "window_scan_pval_all_fragments_filter";
	private Map<SampleData, FileWriter> windowScanPvalWithFragmentLengthFilterRejectFileWriters;
	private static String REJECT_FILE_PREFIX_WINDOW_SCAN_PVAL_WITH_FRAGMENT_LENGTH = "fragment_length_filter";
	private Map<SampleData, FileWriter> peakScanPvalRejectWriters;
	private static String REJECT_FILE_PREFIX_PEAK_SCAN_PVAL = "peak_scan_pval_filter";
	private Map<SampleData, FileWriter> duplicateRejectWriters;
	private static String REJECT_FILE_PREFIX_DUPLICATE = "duplication_filter";
	private Map<SampleData, FileWriter> strandRejectWriters;
	private static String REJECT_FILE_PREFIX_STRAND = "strand_filter";
	private ArrayList<FileWriter> allRejectFileWriters;
	protected static String FILTER_REJECT_DIR = "filter_rejects";

	
	protected MultiSampleScanPeakCaller(MultiSampleScanPeakCaller other) throws IOException {
		this(other.sampleFile, other.bedAnnotationFile, other.sizeFile, other.windowSize, other.stepSize);
		copyParameters(other);
	}
	
	/**
	 * Instantiate with default parameters
	 * @param sampleListFile File containing sample list
	 * @param bedFile Bed gene annotation
	 * @throws IOException
	 */
	@SuppressWarnings("unused")
	private MultiSampleScanPeakCaller(String sampleListFile, String bedFile, String chrSizeFile) throws IOException {
		this(sampleListFile, bedFile, chrSizeFile, DEFAULT_WINDOW_SIZE, DEFAULT_STEP_SIZE);
	}
	
	
	/**
	 * Instantiate using file of sample information
	 * @param sampleListFile File containing sample list
	 * @param bedFile Bed gene annotation
	 * @param chrSizeFile Chromosome size file
	 * @param window Window size
	 * @param step Step size
	 * @throws IOException
	 */
	public MultiSampleScanPeakCaller(String sampleListFile, String bedFile, String chrSizeFile, int window, int step) throws IOException {
		
		sampleFile = sampleListFile;
		bedAnnotationFile = bedFile;
		sizeFile = chrSizeFile;
		windowSize = window;
		stepSize = step;
		genes = BEDFileParser.loadDataByChr(new File(bedFile));
		coord = new TranscriptomeSpace(genes);
		random = new Random();
		
		if(sampleFile != null) {
			initializeSamplesFromSampleListFile();
			initializeScoreMaps();
			initializeProcessors();
		}
		
		trimQuantile = DEFAULT_TRIM_PEAK_QUANTILE;
		peakWindowCountCutoff = DEFAULT_PEAK_WINDOW_COUNT_CUTOFF;
		peakWindowScanPvalCutoff = DEFAULT_PEAK_SCAN_P_VALUE_CUTOFF;
		firstReadTranscriptionStrand = DEFAULT_FIRST_READ_TRANSCRIPTION_STRAND;
		peakCutoffMaxReplicatePct = DEFAULT_PEAK_MAX_PCT_DUPLICATES;
		filterByStrand = DEFAULT_FILTER_BY_STRAND;
		extraFields = DEFAULT_EXTRA_FIELDS;
		useBinomialScore = DEFAULT_USE_BINOMIAL;
		
		
	}
	
	protected void initializeFilterRejectWriters(String commonSuffix, String outDir) throws IOException {
		File dir = new File(outDir);
		@SuppressWarnings("unused")
		boolean madeDir = dir.mkdir();
		logger.info("");
		logger.info("Writing regions rejected by filters to files in directory " + outDir);
		expressionFilterRejectWriters = new HashMap<SampleData, FileWriter>();
		windowCountRejectFileWriters = new HashMap<SampleData, FileWriter>();
		windowScanPvalAllFragmentsRejectFileWriters = new HashMap<SampleData, FileWriter>();
		windowScanPvalWithFragmentLengthFilterRejectFileWriters = new HashMap<SampleData, FileWriter>();
		peakScanPvalRejectWriters = new HashMap<SampleData, FileWriter>();
		duplicateRejectWriters = new HashMap<SampleData, FileWriter>();
		strandRejectWriters = new HashMap<SampleData, FileWriter>();
		allRejectFileWriters = new ArrayList<FileWriter>();

		for(SampleData sample : allSamples) {
			expressionFilterRejectWriters.put(sample, new FileWriter(outDir + "/" + REJECT_FILE_PREFIX_EXPRESSION_FILTER + "_" + sample.getSampleName() + "_" + commonSuffix + ".bed"));
			windowCountRejectFileWriters.put(sample, new FileWriter(outDir + "/" + REJECT_FILE_PREFIX_WINDOW_COUNT + "_" + sample.getSampleName() + "_" + commonSuffix + ".bed"));
			windowScanPvalAllFragmentsRejectFileWriters.put(sample, new FileWriter(outDir + "/" + REJECT_FILE_PREFIX_WINDOW_SCAN_PVAL_ALL_FRAGMENTS + "_" + sample.getSampleName() + "_" + commonSuffix + ".bed"));
			windowScanPvalWithFragmentLengthFilterRejectFileWriters.put(sample, new FileWriter(outDir + "/" + REJECT_FILE_PREFIX_WINDOW_SCAN_PVAL_WITH_FRAGMENT_LENGTH + "_" + sample.getSampleName() + "_" + commonSuffix + ".bed"));
			peakScanPvalRejectWriters.put(sample, new FileWriter(outDir + "/" + REJECT_FILE_PREFIX_PEAK_SCAN_PVAL + "_" + sample.getSampleName() + "_" + commonSuffix + ".bed"));
			duplicateRejectWriters.put(sample, new FileWriter(outDir + "/" + REJECT_FILE_PREFIX_DUPLICATE + "_" + sample.getSampleName() + "_" + commonSuffix + ".bed"));
			strandRejectWriters.put(sample, new FileWriter(outDir + "/" + REJECT_FILE_PREFIX_STRAND + "_" + sample.getSampleName() + "_" + commonSuffix + ".bed"));
			allRejectFileWriters.add(expressionFilterRejectWriters.get(sample));
			allRejectFileWriters.add(windowCountRejectFileWriters.get(sample));
			allRejectFileWriters.add(windowScanPvalAllFragmentsRejectFileWriters.get(sample));
			allRejectFileWriters.add(windowScanPvalWithFragmentLengthFilterRejectFileWriters.get(sample));
			allRejectFileWriters.add(peakScanPvalRejectWriters.get(sample));
			allRejectFileWriters.add(duplicateRejectWriters.get(sample));
			allRejectFileWriters.add(strandRejectWriters.get(sample));
		}
	}
	
	protected void closeFilterRejectWriters() throws IOException {
		for(FileWriter w : allRejectFileWriters) {
			w.close();
		}
	}
	
	private void copyParameters(MultiSampleScanPeakCaller other) {
		setPeakTrimQuantile(other.trimQuantile);
		setPeakWindowCountCutoff(other.peakWindowCountCutoff);
		setPeakWindowScanPvalCutoff(other.peakWindowScanPvalCutoff);
		setPeakCutoffMostCommonReplicate(other.peakCutoffMaxReplicatePct);
		setFirstReadTranscriptionStrand(other.firstReadTranscriptionStrand);
		setExpressionScanPvalueCutoff(other.expressionData.getExpressionScanPvalueCutoff());
		setFilterByStrand(other.filterByStrand);
		setExtraFields(other.extraFields);
		setBinomialScore(other.useBinomialScore);
	}
	
	/**
	 * Set whether to use strand information in reads to filter peaks
	 * @param useStrandFilter Whether to apply strand filter
	 */
	public void setFilterByStrand(boolean useStrandFilter) {
		filterByStrand = useStrandFilter;
	}
	
	/**
	 * Set whether to print extra fields
	 * @param useExtraFields
	 */
	public void setExtraFields(boolean useExtraFields) {
		extraFields = useExtraFields;
	}
	
	/**
	 * Set whether to print extra fields
	 * @param useExtraFields
	 */
	public void setBinomialScore(boolean binomialScore) {
		useBinomialScore = binomialScore;
	}
	
	/**
	 * Set cutoff for the percentage of fragments overlapping a peak that come from the most common replicate fragment
	 * @param maxPct The max percentage
	 */
	public void setPeakCutoffMostCommonReplicate(double maxPct) {
		peakCutoffMaxReplicatePct = maxPct;
	}
	
	/**
	 * Set quantile for trim max contiguous algorithm for peak calling
	 * @param trimPeakQuantile The quantile of data values
	 */
	public void setPeakTrimQuantile(double trimPeakQuantile) {
		trimQuantile = trimPeakQuantile;
	}
	
	/**
	 * Set minimum number of fragments overlapping a window to be considered for possible peak
	 * @param peakCountCutoff Minimum fragment count for each window
	 */
	public void setPeakWindowCountCutoff(double peakCountCutoff) {
		peakWindowCountCutoff = peakCountCutoff;
	}
	
	/**
	 * Set scan P value cutoff for window within transcript to be considered for possible peak
	 * @param peakScanPvalCutoff Max scan P value
	 */
	public void setPeakWindowScanPvalCutoff(double peakScanPvalCutoff) {
		peakWindowScanPvalCutoff = peakScanPvalCutoff;
	}
	
	/**
	 * Set whether read 1 is transcription strand
	 * @param firstReadIsTranscriptionStrand True if read 1 is transcription strand
	 */
	public void setFirstReadTranscriptionStrand(boolean firstReadIsTranscriptionStrand) {
		firstReadTranscriptionStrand = firstReadIsTranscriptionStrand;
	}
	
	/**
	 * Set genome wide scan P value cutoff for expression of transcript
	 * @param expressionScanPvalCutoff P value cutoff for transcript expression against genomic background
	 */
	public void setExpressionScanPvalueCutoff(double expressionScanPvalCutoff) {
		for(SampleData sample : allSamples) {
			sample.setExpressionScanPvalueCutoff(expressionScanPvalCutoff);
		}
	}
	
	/**
	 * Instantiate with a single sample and an expression sample and default paramters
	 * @param expressionBamFile Bam file for expression sample
	 * @param signalBamFile Bam file for signal sample
	 * @param bedFile Bed gene annotation
	 * @param chrSizeFile Chromosome size file
	 * @throws IOException
	 */
	public MultiSampleScanPeakCaller(String expressionBamFile, String signalBamFile, String bedFile, String chrSizeFile) throws IOException {
		this(expressionBamFile, signalBamFile, bedFile, chrSizeFile, DEFAULT_WINDOW_SIZE, DEFAULT_STEP_SIZE);
	}
	
	/**
	 * Instantiate with a single signal sample and an expression sample
	 * @param expressionBamFile Bam file for expression sample
	 * @param signalBamFile Bam file for signal sample
	 * @param bedFile Bed gene annotation
	 * @param chrSizeFile Chromosome size file
	 * @param window Window size
	 * @param step Step size
	 * @throws IOException
	 */
	public MultiSampleScanPeakCaller(String expressionBamFile, String signalBamFile, String bedFile, String chrSizeFile, int window, int step) throws IOException {
		this(null, bedFile, chrSizeFile, window, step);
		initializeWithSingleSample(expressionBamFile, signalBamFile, chrSizeFile);
		initializeScoreMaps();
		if (useBinomialScore) {
			initializeProcessors();
		}
	}
	
	private void initializeWithSingleSample(String expressionBamFile, String signalBamFile, String chrSizeFile) throws IOException {
		controlSamples = new ArrayList<SampleData>();
		numControls = controlSamples.size();
		signalSamples = new ArrayList<SampleData>();
		SampleData signalSample = new SampleData(signalBamFile, firstReadTranscriptionStrand, genes, windowSize, stepSize, DEFAULT_EXPRESSION_SCAN_P_VALUE_CUTOFF, true, false);
		signalSamples.add(signalSample);
		numSignals = signalSamples.size();
		expressionData = new GenomeSpaceSampleData(expressionBamFile, chrSizeFile, genes, windowSize, stepSize, DEFAULT_EXPRESSION_SCAN_P_VALUE_CUTOFF, true);
		allSamples = new ArrayList<SampleData>();
		allSamples.addAll(controlSamples);
		allSamples.addAll(signalSamples);
		numSamples = allSamples.size();
		binomialCtrl = new SampleData(expressionBamFile, firstReadTranscriptionStrand, genes, windowSize, stepSize, DEFAULT_EXPRESSION_SCAN_P_VALUE_CUTOFF, true, false);
	}
	
	private void initializeSamplesFromSampleListFile() throws IOException {
		SampleFileParser p = new SampleFileParser(sampleFile);
		controlSamples = p.getControlDatasets();
		numControls = controlSamples.size();
		signalSamples = p.getSignalDatasets();
		expressionData = p.getExpressionData();
		numSignals = signalSamples.size();
		allSamples = new ArrayList<SampleData>();
		allSamples.addAll(controlSamples);
		allSamples.addAll(signalSamples);
		numSamples = allSamples.size();
		binomialCtrl = p.getBinomialCtrlData();
	}
	
	private void initializeScoreMaps() {
		tStatisticWindowScores = new TreeMap<Gene, Map<Annotation, Double>>();
		singleSampleWindowEnrichmentOverGene = new HashMap<SampleData, Map<Gene, Map<Annotation, Double>>>();
		singleSampleScanPeaks = new HashMap<SampleData, Map<Gene, Collection<Gene>>>();
		binomialWindowScores = new HashMap<SampleData, Map<Gene, Map<Annotation, BinomialEnrichmentScore>>>();
		for(SampleData sample : allSamples) {
			if(!singleSampleWindowEnrichmentOverGene.containsKey(sample)) {
				Map<Gene, Map<Annotation, Double>> m = new TreeMap<Gene, Map<Annotation, Double>>();
				singleSampleWindowEnrichmentOverGene.put(sample, m);
			}
			if(!singleSampleScanPeaks.containsKey(sample)) {
				Map<Gene, Collection<Gene>> m = new TreeMap<Gene, Collection<Gene>>();
				singleSampleScanPeaks.put(sample, m);
			}
			if (!binomialWindowScores.containsKey(sample)) {
				Map<Gene, Map<Annotation, BinomialEnrichmentScore>> m = new TreeMap<Gene, Map<Annotation, BinomialEnrichmentScore>>();
				binomialWindowScores.put(sample, m);
			}
		}
	}
	
	private void initializeProcessors() throws IOException {
		binomialProcessors = new HashMap<SampleData,WindowProcessor<BinomialEnrichmentScore>>();
		for(SampleData sample : allSamples) {
			if (!binomialProcessors.containsKey(sample)) {
				SampleData ctrl = binomialCtrl;
				try {
					WindowProcessor<BinomialEnrichmentScore> p = new BinomialEnrichmentScore.Processor(sample.getData(),ctrl.getData());
					binomialProcessors.put(sample, p);
				} catch(Exception e) {
					try {
						sample.getData();
					} catch(NullPointerException f) {
						logger.debug("sample data is null");
					}
					try {
						ctrl.getData();
					} catch(NullPointerException g) {
						logger.debug("ctrl data is null");
					}
					
				}
				
			}
		}
	}
	
	/**
	 * Whether the gene is expressed in all control samples at the given significance level
	 * @param gene The gene
	 * @return Whether the gene is expressed by these criteria
	 */
	public boolean isExpressedInAllControlSamples(Gene gene) {
		for(SampleData control : controlSamples) {
			if(!control.isExpressed(gene)) {
				return false;
			}
		}
		return true;
	}
	
	/**
	 * Whether the gene is significantly expressed in the special expression sample
	 * @param gene The gene
	 * @return True iff the gene is expressed in the expression sample
	 */
	public boolean isExpressed(Gene gene) {
		return expressionData.isExpressed(gene);
	}
	
	@SuppressWarnings("unused")
	private void batchWriteSingleSampleScanPeaksAllSamples(String[] commandArgs) throws IOException, InterruptedException, DrmaaException {
		batchWriteSingleSampleScanPeaksAllSamples(commandArgs, null, DEFAULT_BATCH_MEM_REQUEST);
	}
	
	@SuppressWarnings("unused")
	private void batchWriteSingleSampleScanPeaksAllSamples(String[] commandArgs, String chrListFile) throws IOException, InterruptedException, DrmaaException {
		batchWriteSingleSampleScanPeaksAllSamples(commandArgs, chrListFile, DEFAULT_BATCH_MEM_REQUEST);
	}
	
	@SuppressWarnings("unused")
	private void batchWriteSingleSampleScanPeaksAllSamples(String[] commandArgs, int memRequestGb) throws IOException, InterruptedException, DrmaaException {
		batchWriteSingleSampleScanPeaksAllSamples(commandArgs, null, memRequestGb);
	}
	
	private void batchWriteSingleSampleScanPeaksAllSamples(String[] commandArgs, String chrListFile, int memRequestGb) throws IOException, InterruptedException, DrmaaException {
		
		logger.info("");
		logger.info("\nBatching out peak calling by sample and chromosome...\n");
		
		int xmx = (int)Math.floor(0.9 * memRequestGb);
		int xms = (int)Math.floor(0.7 * memRequestGb);
		int xmn = (int)Math.floor(0.5 * memRequestGb);
		
		TreeSet<String> chrs = new TreeSet<String>();
		if(chrListFile == null) {
			chrs.addAll(genes.keySet());
		} else {
			FileReader r = new FileReader(chrListFile);
			BufferedReader b = new BufferedReader(r);
			StringParser s = new StringParser();
			while(b.ready()) {
				String line = b.readLine();
				s.parse(line);
				if(s.getFieldCount() == 0) continue;
				if(!genes.keySet().contains(line)) {
					throw new IllegalArgumentException("Chromosome name " + line + " not recognized.");
				}
				chrs.add(line);
			}
			r.close();
			b.close();
		}
		
		String jar = commandLineBatchJar(commandArgs);
		ArrayList<Job> jobs = new ArrayList<Job>();
		String outDir = commandLineOutDir(commandArgs);
		File o = new File(outDir);
		@SuppressWarnings("unused")
		boolean madeDir = o.mkdir();
		if(!o.exists()) {
			throw new IOException("Could not create directory " + outDir);
		}
		
		Map<String, String> cmmds = new TreeMap<String, String>();
		
		for(SampleData sample : allSamples) {
			for(String chr : chrs) {
				String[] batchedCmmdArgs = BatchedMultiSampleScanPeakCaller.extendSuperArgsForSampleAndChr(commandArgs, sample.getSampleName(), chr);
				String args = "";
				for(int i=0; i < batchedCmmdArgs.length; i++) {
					args += batchedCmmdArgs[i] + " ";
				}
				String cmmd = "java -jar -Xmx" + xmx + "g -Xms" + xms + "g -Xmn" + xmn + "g " + jar + " " + args;
				logger.info("Running command: " + cmmd);
				String jobID = sample.getSampleName() + "_" + chr + "_" + Long.valueOf(System.currentTimeMillis()).toString();
				LSFJob job = new LSFJob(Runtime.getRuntime(), jobID, cmmd, outDir + "/" + jobID + ".bsub", "week", memRequestGb);
				jobs.add(job);
				cmmds.put(jobID, cmmd);
				logger.info("LSF job ID is " + jobID + ".");
				// Submit job
				job.submit();

			}
		}

		
		logger.info("");
		logger.info("Waiting for jobs to finish...");
		JobUtils.waitForAll(jobs);
		
		logger.info("\nAll jobs finished.\n");
		
	}
	
	/**
	 * Get name of bed file to write for peaks
	 * @param sample Sample
	 * @param outDir Output directory name or null if current directory
	 * @param chrName Chromosome name or null if all chromosomes
	 * @return File name
	 */
	protected String getPeakBedFileName(SampleData sample, String outDir, String chrName) {
		String rtrn = "";
		if(outDir != null) {
			rtrn += outDir + "/";
		}
		rtrn += sample.getSampleName() + "_scan_peaks_" + windowSize + "_" + stepSize + "_" + peakWindowScanPvalCutoff + "_" + trimQuantile;
		if(chrName != null) {
			rtrn += "_" + chrName;
		}
		rtrn += ".bed";
		return rtrn;
	}
	
	/**
	 * Write scan peaks for all samples to separate bed files
	 * @throws IOException
	 */
	private void writeSingleSampleScanPeaksAllSamples(String outDir) throws IOException {
		logger.info("");
		logger.info("Writing single sample scan peaks for each sample...");
		File o = new File(outDir);
		@SuppressWarnings("unused")
		boolean madeDir = o.mkdir();
		if(!o.exists()) {
			throw new IOException("Could not create directory " + outDir);
		}
		for(SampleData signal : signalSamples) {
			String outfile = getPeakBedFileName(signal, outDir, null);
			writeSingleSampleScanPeaks(signal, outfile);
		}		
		for(SampleData control : controlSamples) {
			String outfile = getPeakBedFileName(control, outDir, null);
			writeSingleSampleScanPeaks(control, outfile);
		}
		logger.info("Done writing single sample scan peaks for all samples.");
	}
	
	/**
	 * Write all single sample scan peaks for the sample to file
	 * @param sample The sample
	 * @param outFile Output bed file
	 * @param r Red value for bed file color
	 * @param g Green value for bed file color
	 * @param b Blue value for bed file color
	 * @throws IOException 
	 */
	private void writeSingleSampleScanPeaks(SampleData sample, String outFile) throws IOException {
		writeSingleSampleScanPeaks(sample, outFile, null);
	}
	
	/**
	 * Write all single sample scan peaks for the sample and chromosome to file
	 * @param sample The sample
	 * @param outFile Output bed file
	 * @param r Red value for bed file color
	 * @param g Green value for bed file color
	 * @param b Blue value for bed file color
	 * @param chrName Only write peaks for this chromosome
	 * @throws IOException 
	 */
	protected void writeSingleSampleScanPeaks(SampleData sample, String outFile, String chrName) throws IOException {
		logger.info("Writing single sample scan peaks for sample " + sample.getSampleName() + " to file " + outFile + "...");
		FileWriter w = new FileWriter(outFile);
		for(String chr : genes.keySet()) {
			if(chrName != null) {
				if(!chr.equals(chrName)) {
					continue;
				}
			}
			for(Gene gene : genes.get(chr)) {
				Collection<Gene> peaks = getSingleSampleScanPeaks(sample, gene);
				for(Gene window : peaks) {
					GeneWindow geneWindow = new GeneWindow(window); 
					int r = RGB_RED_UNKNOWN;
					int g = RGB_GREEN_UNKNOWN;
					int b = RGB_BLUE_UNKNOWN;
					if(window.getOrientation().equals(Strand.UNKNOWN)) {
						String name = window.getName();
						name += "_STRAND_UNKNOWN";
						geneWindow.setName(name);						
					}
					if(window.getOrientation().equals(gene.getOrientation()) && !window.getOrientation().equals(Strand.UNKNOWN)) {
						r = RGB_RED_WITH_GENE;
						g = RGB_GREEN_WITH_GENE;
						b = RGB_BLUE_WITH_GENE;
					}
					if(!window.getOrientation().equals(gene.getOrientation()) && !window.getOrientation().equals(Strand.UNKNOWN)) {
						String name = window.getName();
						name += "_STRAND_AGAINST_GENE";
						geneWindow.setName(name);
						r = RGB_RED_AGAINST_GENE;
						g = RGB_GREEN_AGAINST_GENE;
						b = RGB_BLUE_AGAINST_GENE;
					}
					geneWindow.setBedScore(window.getScore());
					if (extraFields) {
						w.write(geneWindow.toBED(true, r, g, b) + "\n");
					} else {
						w.write(geneWindow.toBED(r, g, b) + "\n");
					}
				}
			}
		}
		w.close();
		logger.info("Done writing scan peaks for sample " + sample.getSampleName() + ".");
	}
	
	/**
	 * Get scan peaks for the sample and the gene
	 * @param sample The sample
	 * @param gene The gene
	 * @return Significant scan peaks
	 * @throws IOException 
	 */
	public Collection<Gene> getSingleSampleScanPeaks(SampleData sample, Gene gene) throws IOException {
		if(singleSampleScanPeaks.get(sample).containsKey(gene)) {
			return singleSampleScanPeaks.get(sample).get(gene);
		}
		identifySingleSampleScanPeaks(sample, gene, expressionFilterRejectWriters.get(sample));
		return singleSampleScanPeaks.get(sample).get(gene);
	}
	
	/**
	 * Identify significant peaks within the sample and the gene
	 * @param sample The sample
	 * @param gene The gene
	 * @throws IOException
	 */
	private void identifySingleSampleScanPeaks(SampleData sample, Gene gene, FileWriter rejectFileWriterExpression) throws IOException {
		
		TreeSet<Annotation> finalPeaks = new TreeSet<Annotation>();
		TreeSet<Gene> rtrnPeaks = new TreeSet<Gene>();
		TranscriptomeSpaceAlignmentModel data = sample.getData();
		TreeSet<Annotation> scanSignificantWindows = new TreeSet<Annotation>();
		
		// If gene is not expressed, skip
		if(!isExpressed(gene)) {
			logger.info("Gene " + gene.getName() + " (" + gene.getChr() + ":" + gene.getStart() + "-" + gene.getEnd() + ") not expressed in expression dataset.");
			singleSampleScanPeaks.get(sample).put(gene, rtrnPeaks);
			rejectFileWriterExpression.write(gene.toBED() + "\n");
			return;
		}
				
		logger.info("Finding scan peaks for sample " + sample.getSampleName() + " and gene " + gene.getName() + " (" + gene.getChr() + ":" + gene.getStart() + "-" + gene.getEnd() + ")");
		
		// Get fixed size windows with sufficient count and significant scan statistic
		// Note: this is where I need to split for binomial score
		if (!useBinomialScore) {
			scanSignificantWindows = findScanSignificantWindows(sample, gene, windowCountRejectFileWriters.get(sample), windowScanPvalAllFragmentsRejectFileWriters.get(sample), windowScanPvalWithFragmentLengthFilterRejectFileWriters.get(sample));
		} else {
			scanSignificantWindows = findBinomialSignificantWindows(sample, gene, windowCountRejectFileWriters.get(sample), windowScanPvalAllFragmentsRejectFileWriters.get(sample), windowScanPvalWithFragmentLengthFilterRejectFileWriters.get(sample));
		}
		// If no significant windows return
		if(scanSignificantWindows.isEmpty()) {
			singleSampleScanPeaks.get(sample).put(gene, rtrnPeaks);
			return;
		}
		
		// Merge overlapping windows
		Collection<Annotation> mergedWindows = mergePeaks(scanSignificantWindows);
		
		// Trim each window
		TreeSet<Annotation> mergedTree = new TreeSet<Annotation>();
		mergedTree.addAll(mergedWindows);
		TreeSet<Annotation> trimmedMergedWindows = trimWindows(mergedTree, data);
		
		// Filter by scan statistic again
		
		TreeSet<Annotation> scanSigWindows = filterByScanStatistic(trimmedMergedWindows, sample, gene, peakScanPvalRejectWriters.get(sample));
		
		// Filter on strand information
		TreeSet<Annotation> correctStrandWindows = new TreeSet<Annotation>();
		if(filterByStrand) {
			correctStrandWindows.addAll(filterByStrandInformation(scanSigWindows, sample, gene, strandRejectWriters.get(sample)));
		} else {
			correctStrandWindows.addAll(scanSigWindows);
		}
		
		// Filter on percent duplicates
		TreeSet<Annotation> dupsOk = filterByDuplicates(correctStrandWindows, sample, peakCutoffMaxReplicatePct, duplicateRejectWriters.get(sample));
		
		// Final peaks
		finalPeaks.addAll(dupsOk);
		
		double geneAvgCoverage = sample.getGeneAverageCoverage(gene);
		double geneCount = sample.getGeneCount(gene);
		int geneSize = coord.getSize(gene);
		
		// Add finishing touches to peaks
		for(Annotation peak : finalPeaks) {
			
			Gene window = new Gene(peak);
			// Name peaks
			window.setName(gene.getName() + ":" + window.getChr() + ":" + window.getStart() + "-" + window.getEnd());
			
			// Set peak score to enrichment
			double enrichment = sample.getEnrichmentOverGene(gene, window);
			double pval;
			double[] extraFields;
			if (!useBinomialScore) {
				ScanStatisticScore score = sample.scoreWindow(gene, window);
				pval = score.getScanPvalue();
				extraFields = new double[6];
				extraFields[0] = score.getCount();
				extraFields[1] = score.getAverageCoverage(data);
				extraFields[2] = geneCount;
				extraFields[3] = geneAvgCoverage;
				extraFields[4] = (double) geneSize;
				extraFields[5] = pval;
			}
			else {
				BinomialEnrichmentScore score = scoreWindowBinomial(gene,window,sample);
				score.setRegionLength(geneSize);
				score.refreshPvalue();
				pval = score.getPvalue();
				extraFields = new double[6];
				extraFields[0] = score.getSampleCount();
				extraFields[1] = score.getCtrlCount();
				extraFields[2] = score.getSampleRegionCount();
				extraFields[3] = score.getCtrlRegionCount();
				extraFields[4] = score.getRegionLength();
				extraFields[5] = pval;
				enrichment = score.getEnrichmentOverControl();
			}
			
			
			
			logger.debug("FINAL_PEAK\t" + gene.getName());
			logger.debug("FINAL_PEAK\t" + window.toBED());
			//logger.debug("FINAL_PEAK\tname=" + window.getName());
			//logger.debug("FINAL_PEAK\twindow_count=" + windowCount);
			//logger.debug("FINAL_PEAK\twindow_size=" + coord.getSize(window));
			//logger.debug("FINAL_PEAK\twindow_avg_coverage=" + windowAvgCoverage);
			//logger.debug("FINAL_PEAK\tgene_count=" + geneCount);
			//logger.debug("FINAL_PEAK\tgene_size=" + geneSize);
			//logger.debug("FINAL_PEAK\tgene_avg_coverage=" + geneAvgCoverage);
			//logger.debug("FINAL_PEAK\tenrichment_over_transcript=" + enrichment);
			//logger.debug("FINAL_PEAK\tscore=" + window.getScore());
			//logger.debug("FINAL_PEAK\torientation=" + window.getOrientation().toString());
			
			window.setScore(enrichment);
			window.setExtraFields(extraFields);
			
			rtrnPeaks.add(window);

		}
		
		singleSampleScanPeaks.get(sample).put(gene, rtrnPeaks);
	}
	
	private TreeSet<Annotation> findScanSignificantWindows(SampleData sample, Gene gene, FileWriter windowCountRejectFileWriter, FileWriter rejectFileWriterAllFragments, FileWriter rejectFileWriterFragmentLengthFilter) throws IOException {
		TreeSet<Annotation> rtrn = new TreeSet<Annotation>();
		Map<Annotation, ScanStatisticScore> windowScores = sample.getWindowScores(gene);
		for(Annotation window : windowScores.keySet()) {
			ScanStatisticScore score = windowScores.get(window);
			double count = score.getCount();
			if(count < peakWindowCountCutoff) {
				windowCountRejectFileWriter.write(window.toBED() + "\n");
				continue;
			}
			double pval = score.getScanPvalue();
			if(pval < peakWindowScanPvalCutoff) {
				//logger.debug("FIXED_SIZE_WINDOW_IS_SIGNIFICANT_BEFORE_FRAGMENT_LENGTH_FILTER\t" + gene.getName());
				//logger.debug("FIXED_SIZE_WINDOW_IS_SIGNIFICANT_BEFORE_FRAGMENT_LENGTH_FILTER\t" + window.toBED());
				//logger.debug("FIXED_SIZE_WINDOW_IS_SIGNIFICANT_BEFORE_FRAGMENT_LENGTH_FILTER\tglobal_length=" + score.getGlobalLength());
				//logger.debug("FIXED_SIZE_WINDOW_IS_SIGNIFICANT_BEFORE_FRAGMENT_LENGTH_FILTER\tglobal_count=" + score.getTotal());
				//logger.debug("FIXED_SIZE_WINDOW_IS_SIGNIFICANT_BEFORE_FRAGMENT_LENGTH_FILTER\tglobal_lambda=" + score.getGlobalLambda());
				//logger.debug("FIXED_SIZE_WINDOW_IS_SIGNIFICANT_BEFORE_FRAGMENT_LENGTH_FILTER\twindow_size=" + score.getCoordinateSpace().getSize(window));
				//logger.debug("FIXED_SIZE_WINDOW_IS_SIGNIFICANT_BEFORE_FRAGMENT_LENGTH_FILTER\twindow_count=" + score.getCount());
				//logger.debug("FIXED_SIZE_WINDOW_IS_SIGNIFICANT_BEFORE_FRAGMENT_LENGTH_FILTER\tpval=" + score.getScanPvalue());
				ScanStatisticScore fragmentLengthFilterScore = sample.scoreWindowWithFragmentLengthFilter(gene, window);
				double pval2 = fragmentLengthFilterScore.getScanPvalue();
				if(pval2 < peakWindowScanPvalCutoff) {
					//logger.debug("FIXED_SIZE_WINDOW_IS_SIGNIFICANT_AFTER_FRAGMENT_LENGTH_FILTER\tglobal_length=" + fragmentLengthFilterScore.getGlobalLength());
					//logger.debug("FIXED_SIZE_WINDOW_IS_SIGNIFICANT_AFTER_FRAGMENT_LENGTH_FILTER\tglobal_count=" + fragmentLengthFilterScore.getTotal());
					//logger.debug("FIXED_SIZE_WINDOW_IS_SIGNIFICANT_AFTER_FRAGMENT_LENGTH_FILTER\tglobal_lambda=" + fragmentLengthFilterScore.getGlobalLambda());
					//logger.debug("FIXED_SIZE_WINDOW_IS_SIGNIFICANT_AFTER_FRAGMENT_LENGTH_FILTER\twindow_size=" + fragmentLengthFilterScore.getCoordinateSpace().getSize(window));
					//logger.debug("FIXED_SIZE_WINDOW_IS_SIGNIFICANT_AFTER_FRAGMENT_LENGTH_FILTER\twindow_count=" + fragmentLengthFilterScore.getCount());
					//logger.debug("FIXED_SIZE_WINDOW_IS_SIGNIFICANT_AFTER_FRAGMENT_LENGTH_FILTER\tpval=" + fragmentLengthFilterScore.getScanPvalue());				
					rtrn.add(window);
				} else {
					//logger.debug("FIXED_SIZE_WINDOW_NOT_SIGNIFICANT_AFTER_FRAGMENT_LENGTH_FILTER\tglobal_length=" + fragmentLengthFilterScore.getGlobalLength());
					//logger.debug("FIXED_SIZE_WINDOW_NOT_SIGNIFICANT_AFTER_FRAGMENT_LENGTH_FILTER\tglobal_count=" + fragmentLengthFilterScore.getTotal());
					//logger.debug("FIXED_SIZE_WINDOW_NOT_SIGNIFICANT_AFTER_FRAGMENT_LENGTH_FILTER\tglobal_lambda=" + fragmentLengthFilterScore.getGlobalLambda());
					//logger.debug("FIXED_SIZE_WINDOW_NOT_SIGNIFICANT_AFTER_FRAGMENT_LENGTH_FILTER\twindow_size=" + fragmentLengthFilterScore.getCoordinateSpace().getSize(window));
					//logger.debug("FIXED_SIZE_WINDOW_NOT_SIGNIFICANT_AFTER_FRAGMENT_LENGTH_FILTER\twindow_count=" + fragmentLengthFilterScore.getCount());
					//logger.debug("FIXED_SIZE_WINDOW_NOT_SIGNIFICANT_AFTER_FRAGMENT_LENGTH_FILTER\tpval=" + fragmentLengthFilterScore.getScanPvalue());				
					rejectFileWriterFragmentLengthFilter.write(window.toBED() + "\n");
				}
			} else {
				rejectFileWriterAllFragments.write(window.toBED() + "\n");
			}
		}
		return rtrn;
	}
	
	private TreeSet<Annotation> findBinomialSignificantWindows(SampleData sample, Gene gene, FileWriter windowCountRejectFileWriter, FileWriter rejectFileWriterAllFragments, FileWriter rejectFileWriterFragmentLengthFilter) throws IOException {
		TreeSet<Annotation> rtrn = new TreeSet<Annotation>();
		Map<Annotation, BinomialEnrichmentScore> windowScores = getBinomialWindowScores(sample,gene);
		for(Annotation window : windowScores.keySet()) {
			BinomialEnrichmentScore score = windowScores.get(window);
			double count = score.getCount();
			if(count < peakWindowCountCutoff) {
				windowCountRejectFileWriter.write(window.toBED() + "\n");
				continue;
			}
			double pval = score.getPvalue();
			if(pval < peakWindowScanPvalCutoff) {
				logger.debug("FIXED_SIZE_WINDOW_IS_SIGNIFICANT_BEFORE_FRAGMENT_LENGTH_FILTER\t" + gene.getName());
				logger.debug("FIXED_SIZE_WINDOW_IS_SIGNIFICANT_BEFORE_FRAGMENT_LENGTH_FILTER\t" + window.toBED());
				logger.debug("FIXED_SIZE_WINDOW_IS_SIGNIFICANT_BEFORE_FRAGMENT_LENGTH_FILTER\tsampleCounts=" + score.getSampleCount());
				logger.debug("FIXED_SIZE_WINDOW_IS_SIGNIFICANT_BEFORE_FRAGMENT_LENGTH_FILTER\tctrlCounts=" + score.getCtrlCount());
				logger.debug("FIXED_SIZE_WINDOW_IS_SIGNIFICANT_BEFORE_FRAGMENT_LENGTH_FILTER\tsampleRegionCounts=" + score.getSampleRegionCount());
				logger.debug("FIXED_SIZE_WINDOW_IS_SIGNIFICANT_BEFORE_FRAGMENT_LENGTH_FILTER\tctrlRegionCounts=" + score.getCtrlRegionCount());
				logger.debug("FIXED_SIZE_WINDOW_IS_SIGNIFICANT_BEFORE_FRAGMENT_LENGTH_FILTER\tpval=" + score.getPvalue());
				BinomialEnrichmentScore fragmentLengthFilterScore = scoreWindowBinomialWithFragmentLengthFilter(gene, window, sample);
				double pval2 = fragmentLengthFilterScore.getPvalue();
				if(pval2 < peakWindowScanPvalCutoff) {
				rtrn.add(window);
				} else {
				rejectFileWriterFragmentLengthFilter.write(window.toBED() + "\n");
				}
			} else {
				logger.debug("FIXED_SIZE_WINDOW_IS_NOT_SIGNIFICANT\t" + gene.getName());
				logger.debug("FIXED_SIZE_WINDOW_IS_NOT_SIGNIFICANT\t" + window.toBED());
				logger.debug("FIXED_SIZE_WINDOW_IS_NOT_SIGNIFICANT\tsampleCounts=" + score.getSampleCount());
				logger.debug("FIXED_SIZE_WINDOW_IS_NOT_SIGNIFICANT\tctrlCounts=" + score.getCtrlCount());
				logger.debug("FIXED_SIZE_WINDOW_IS_NOT_SIGNIFICANT\tsampleRegionCounts=" + score.getSampleRegionCount());
				logger.debug("FIXED_SIZE_WINDOW_IS_NOT_SIGNIFICANT\tctrlRegionCounts=" + score.getCtrlRegionCount());
				logger.debug("FIXED_SIZE_WINDOW_IS_NOT_SIGNIFICANT\tpval=" + score.getPvalue());

				rejectFileWriterAllFragments.write(window.toBED() + "\n");
			}
		}
		return rtrn;
	}
	
	private Map<Annotation, BinomialEnrichmentScore> getBinomialWindowScores(SampleData sample,Gene gene) {
		if (binomialWindowScores.get(sample).containsKey(gene)) {
			return binomialWindowScores.get(sample).get(gene);
		} else {
			computeBinomialWindowScores(sample,gene);
			return binomialWindowScores.get(sample).get(gene);
		}
	}
	
	private void computeBinomialWindowScores(SampleData sample, Gene gene) {
		Map<Annotation, BinomialEnrichmentScore> scores = new TreeMap<Annotation, BinomialEnrichmentScore>();
		if(gene.getSize() < windowSize) {
			logger.info(gene.getName() + " is smaller than window size. Not computing window binding site scores.");
			binomialWindowScores.get(sample).put(gene, scores);
			return;
		}
		WindowScoreIterator<BinomialEnrichmentScore> iter = null;
		try {
			iter = sample.getData().scan(gene,sample.getWindowSize(),sample.getWindowSize()-sample.getStepSize(),binomialProcessors.get(sample));
		} catch(NullPointerException e) {
			logger.info("Gene: " + gene.getName());
			logger.info("Gene: " + gene.toBED());
		}
		
		double sampleGeneCount = sample.getGeneCount(gene);
		double ctrlGeneCount = binomialCtrl.getGeneCount(gene);
		while (iter.hasNext()) {
			BinomialEnrichmentScore score = iter.next();
			Annotation window = score.getAnnotation();
			score.setSampleRegionCount(sampleGeneCount);
			score.setCtrlRegionCount(ctrlGeneCount);
			score.setRegionLength(gene.size());
			score.refreshPvalue();
			//logger.debug("SCORE_ALL_WINDOWS_IN_GENE\t" + gene.getName());
			//logger.debug("SCORE_ALL_WINDOWS_IN_GENE\t" + window.getChr() + ":" + window.getStart() + "-" + window.getEnd());
			//logger.debug("SCORE_ALL_WINDOWS_IN_GENE\tsample_count=" + score.getSampleCount());
			//logger.debug("SCORE_ALL_WINDOWS_IN_GENE\tctrl_count=" + score.getCtrlCount());
			//logger.debug("SCORE_ALL_WINDOWS_IN_GENE\tsample_gene_count=" + score.getSampleRegionCount());
			//logger.debug("SCORE_ALL_WINDOWS_IN_GENE\tctrl_gene_count=" + score.getCtrlRegionCount());
			//logger.debug("SCORE_ALL_WINDOWS_IN_GENE\tpval=" + score.getPvalue());
			scores.put(window, score);
		}
		binomialWindowScores.get(sample).put(gene, scores);
	}
	
	private BinomialEnrichmentScore scoreWindowBinomial(Gene gene, Annotation window, SampleData sample, double count) {
		BinomialEnrichmentScore score = scoreWindowBinomial(gene,window,sample);
		score.setCount(count);
		score.refreshPvalue();
		return score;
	}
	
	private BinomialEnrichmentScore scoreWindowBinomial(Gene gene, Annotation window, SampleData sample) {
		double sampleGeneCount = sample.getGeneCount(gene);
		double ctrlGeneCount = binomialCtrl.getGeneCount(gene);
		double regionLength = gene.getSize();
		BinomialEnrichmentScore score = new BinomialEnrichmentScore(sample.getData(),binomialCtrl.getData(),window,regionLength);
		score.setSampleRegionCount(sampleGeneCount);
		score.setCtrlRegionCount(ctrlGeneCount);
		score.refreshPvalue();
		return score;
	}
	
	private BinomialEnrichmentScore scoreWindowBinomialWithFragmentLengthFilter(Gene gene, Annotation window, SampleData sample) {
		double sampleGeneCount = sample.getFragmentLengthFilterData().getCount(gene);
		double ctrlGeneCount = binomialCtrl.getFragmentLengthFilterData().getCount(gene);
		double regionLength = gene.getSize();
		BinomialEnrichmentScore score = new BinomialEnrichmentScore(sample.getFragmentLengthFilterData(),binomialCtrl.getFragmentLengthFilterData(),window,regionLength);
		score.setSampleRegionCount(sampleGeneCount);
		score.setCtrlRegionCount(ctrlGeneCount);
		score.refreshPvalue();
		return score;
	}
	
	private TreeSet<Annotation> trimWindows(TreeSet<Annotation> untrimmed, TranscriptomeSpaceAlignmentModel data) throws IOException {
		TreeSet<Annotation> rtrn = new TreeSet<Annotation>();
		for(Annotation window : untrimmed) {
			List<Double> coverageData = data.getPositionCountList(new Gene(window));
			Annotation trimmed = SampleData.trimMaxContiguous(window, coverageData, trimQuantile);
			rtrn.add(trimmed);
			//logger.debug("MERGED_TRIMMED_WINDOW\t" + window.toBED());
		}
		return rtrn;
	}
	
	private TreeSet<Annotation> filterByScanStatistic(TreeSet<Annotation> preFilter, SampleData sample, Gene gene, FileWriter rejectFileWriter) throws IOException {
		TreeSet<Annotation> rtrn = new TreeSet<Annotation>();
		double p;
		for(Annotation window : preFilter) {
			if (useBinomialScore) {
				BinomialEnrichmentScore score = scoreWindowBinomial(gene,window,sample);
				p = score.getPvalue();
			} else {
				ScanStatisticScore score = sample.scoreWindow(gene, window);
				p = score.getScanPvalue();
			}
			if(p < peakWindowScanPvalCutoff) {
				logger.debug("MERGED_TRIMMED_WINDOW_IS_SIGNIFICANT\t" + gene.getName());
				//logger.debug("MERGED_TRIMMED_WINDOW_IS_SIGNIFICANT\t" + window.toBED());
				//logger.debug("MERGED_TRIMMED_WINDOW_IS_SIGNIFICANT\tglobal_length=" + score.getGlobalLength());
				//logger.debug("MERGED_TRIMMED_WINDOW_IS_SIGNIFICANT\tglobal_count=" + score.getTotal());
				//logger.debug("MERGED_TRIMMED_WINDOW_IS_SIGNIFICANT\tglobal_lambda=" + score.getGlobalLambda());
				//logger.debug("MERGED_TRIMMED_WINDOW_IS_SIGNIFICANT\twindow_size=" + score.getCoordinateSpace().getSize(window));
				//logger.debug("MERGED_TRIMMED_WINDOW_IS_SIGNIFICANT\twindow_count=" + score.getCount());
				//logger.debug("MERGED_TRIMMED_WINDOW_IS_SIGNIFICANT\tpval=" + score.getScanPvalue());
				rtrn.add(window);
			} else {
				rejectFileWriter.write(window.toBED() + "\n");
			}
		}
		return rtrn;
	}
	
	private static TreeSet<Annotation> filterByDuplicates(TreeSet<Annotation> preFilter, SampleData sample, double maxPctMostCommonRead, FileWriter rejectFileWriter) throws IOException {
		TreeSet<Annotation> rtrn = new TreeSet<Annotation>();
		AlignmentModel data = sample.getData();
		for(Annotation window : preFilter) {
			Map<Alignment, Integer> replicateCounts = data.getOverlappingReadReplicateCounts(window, false);
			int total = 0;
			int largestReplicate = 0;
			String mostCommon = "";
			for(Alignment read : replicateCounts.keySet()) {
				int count = replicateCounts.get(read).intValue();
				total += count;
				if(count > largestReplicate) {
					largestReplicate = count;
					mostCommon = read.getChr() + ":" + read.getStart() + "-" + read.getEnd();
				}
			}
			double mostCommonPct = (double) largestReplicate / (double) total;
			if(mostCommonPct > maxPctMostCommonRead) {
				//logger.debug("TOO_MANY_DUPLICATE_READS\t" + mostCommon + "\t" + mostCommonPct + " duplicates");
				rejectFileWriter.write(window.toBED() + "\n");
				continue;
			}
			//logger.debug("NUMBER_OF_DUPLICATE_READS_OK\t" + mostCommon + "\t" +  mostCommonPct + " duplicates");
			rtrn.add(window);
		}
		return rtrn;
	}
	
	private static TreeSet<Annotation> filterByStrandInformation(TreeSet<Annotation> preFilter, SampleData sample, Gene gene, FileWriter rejectFileWriter) throws IOException {
		TreeSet<Annotation> rtrn = new TreeSet<Annotation>();
		for(Annotation window : preFilter) {
			Strand orientation = AlignmentUtils.assignOrientationToWindow(sample.getOriginalBamFile(), window, sample.firstReadTranscriptionStrand(), 0.9);
			window.setOrientation(orientation);
			if(orientation.equals(gene.getOrientation())) {
				//logger.debug("STRAND_OK\t" + gene.toBED());
				//logger.debug("STRAND_OK\t" + window.toBED());
				rtrn.add(window);
			} else if(orientation.equals(Strand.UNKNOWN)) {
				//logger.debug("STRAND_UNKNOWN_SKIPPING\t" + gene.toBED());
				//logger.debug("STRAND_UNKNOWN_SKIPPING\t" + window.toBED());				
				rejectFileWriter.write(window.toBED() + "\n");
			} else {
				//logger.debug("WRONG_STRAND_SKIPPING\t" + gene.toBED());
				//logger.debug("WRONG_STRAND_SKIPPING\t" + window.toBED());				
				rejectFileWriter.write(window.toBED() + "\n");
			}
		}
		return rtrn;
	}
	
	private void computeSingleSampleWindowEnrichmentsOverGenes() {
		logger.info("Computing window enrichments for each sample...");
		for(String chr : genes.keySet()) {
			for(Gene gene : genes.get(chr)) {
				if(!isExpressed(gene)) {
					continue;
				}
				computeSingleSampleWindowEnrichmentOverGene(gene);
			}
		}		
		//writeSingleSampleWindowScoresToFileIfNeeded();
		logger.info("Done computing single sample window enrichments.");
	}
	
	/**
	 * Score all genes
	 * @throws IOException 
	 */
	@SuppressWarnings("unused")
	private void scoreGenesTStatisticScore() throws IOException {
		
		if(!permutationScoring) {
			throw new IllegalStateException("Must enable 'do permutation scoring' option");
		}
		
		computeSingleSampleWindowEnrichmentsOverGenes();
		for(String chr : genes.keySet()) {
			logger.info("Scoring genes on chromosome " + chr);
			for(Gene gene : genes.get(chr)) {
				if(!isExpressed(gene)) {
					logger.info(gene.getName() + " is not expressed. Skipping.");
					continue;
				}
				logger.info("Scoring gene " + gene.getName());
				computeTStatisticWindowScores(gene);
			}
		}
	}
	
	
/*	*//**
	 * For samples that didn't read window scores from file, write to files
	 * @throws IOException
	 *//*
	private void writeSingleSampleWindowScoresToFileIfNeeded() throws IOException {
		for(SampleData sample : allSamples) {
			if(!sample.gotWindowScoresFromFile()) {
				sample.writeWindowScoresToFile(this);
			}
		}
		logger.info("Done writing window score files.");
	}
*/	
	/**
	 * For each sample compute the enrichment of each window over the gene average
	 * Cache the window enrichments
	 * @param gene The gene
	 */
	private void computeSingleSampleWindowEnrichmentOverGene(Gene gene) {
		logger.info("Getting enrichment for each window and each sample...");
		for(SampleData sample : allSamples) {
			Map<Annotation,Double> sampleWindowEnrichments = new TreeMap<Annotation,Double>();
			double geneAvgCoverage = sample.getGeneAverageCoverage(gene);
			logger.info("Sample " + sample.getSampleName() + " Gene average coverage = " + geneAvgCoverage);
			Map<Annotation, ScanStatisticScore> scores = sample.getWindowScores(gene);
			
			if(gene.getSize() < windowSize) {
				logger.info(gene.getName() + " is smaller than window size. Not computing single sample window enrichments.");
				singleSampleWindowEnrichmentOverGene.get(sample).put(gene, sampleWindowEnrichments);
				continue;
			}
			for(Annotation window : scores.keySet()) {
				double windowAvgCoverage = scores.get(window).getAverageCoverage(sample.getData());
				double enrichment = windowAvgCoverage / geneAvgCoverage;
				logger.info(sample.getSampleName() + "\t" + gene.getName() + "\t" + window.getChr() + ":" + window.getStart() + "-" + window.getEnd() + "\tavg_coverage=" + windowAvgCoverage + "\twindow_enrichment=" + enrichment);
				sampleWindowEnrichments.put(window, Double.valueOf(enrichment));
			}
			singleSampleWindowEnrichmentOverGene.get(sample).put(gene, sampleWindowEnrichments);
		}
	}
	
	/**
	 * Compute t statistic score for each window of gene and store scores
	 * @param gene The gene
	 */
	private void computeTStatisticWindowScores(Gene gene) {
		
		if(!permutationScoring) {
			throw new IllegalStateException("Must enable 'do permutation scoring' option");
		}
		
		Map<Annotation, Double> tStatisticScoresThisGene = new TreeMap<Annotation, Double>();
		if(gene.getSize() < windowSize) {
			logger.info(gene.getName() + " is smaller than window size. Not computing t statistic window scores.");
			tStatisticWindowScores.put(gene, tStatisticScoresThisGene);
			return;
		}		
		SampleData tmp = allSamples.iterator().next();
		Collection<Annotation> windows = singleSampleWindowEnrichmentOverGene.get(tmp).get(gene).keySet();
		for(Annotation window : windows) {
			double score = tStatisticWindowScore(gene, window, controlSamples, signalSamples);
			tStatisticScoresThisGene.put(window, Double.valueOf(score));
			String e = "";
			for(SampleData control : controlSamples) {
				e += control.getSampleName() + ":" + control.getWindowScores(gene).get(window).getCount() + ":" + singleSampleWindowEnrichmentOverGene.get(control).get(gene).get(window).toString() + "\t";
			}
			for(SampleData signal : controlSamples) {
				e += signal.getSampleName() + ":" + signal.getWindowScores(gene).get(window).getCount() + ":" + singleSampleWindowEnrichmentOverGene.get(signal).get(gene).get(window).toString() + "\t";
			}
			logger.info(gene.getName() + "\t" + window.getChr() + ":" + window.getStart() + "-" + window.getEnd() + "\t" + score + "\t" + e);
		}
		tStatisticWindowScores.put(gene, tStatisticScoresThisGene);
	}
	
	/**
	 * The t statistic score for a single window
	 * T statistic between control and signal samples
	 * Statistic is positive iff mean of signal enrichments is greater than mean of control enrichments
	 * @param gene The parent gene
	 * @param window The window
	 * @param controls Control samples
	 * @param signals Signal samples
	 * @return The score for the window
	 */
	private double tStatisticWindowScore(Gene gene, Annotation window, Collection<SampleData> controls, Collection<SampleData> signals) {
		
		if(!permutationScoring) {
			throw new IllegalStateException("Must enable 'do permutation scoring' option");
		}
		
		List<Double> controlEnrichments = new ArrayList<Double>();
		List<Double> signalEnrichments = new ArrayList<Double>();
		for(SampleData control : controls) {
			controlEnrichments.add(singleSampleWindowEnrichmentOverGene.get(control).get(gene).get(window));
		}
		for(SampleData signal : signals) {
			signalEnrichments.add(singleSampleWindowEnrichmentOverGene.get(signal).get(gene).get(window));
		}
		return Statistics.tstat(signalEnrichments, controlEnrichments);
	}
	
	
	/**
	 * Nominal P value calculated relative to null distribution of t statistic score
	 * Where null distribution is determined by permuting the contol/signal labels of samples
	 * @param gene The gene the window belongs to
	 * @param window The window
	 * @return The nominal P value for the score of the window
	 */
	@SuppressWarnings("unused")
	private double empiricalNominalPvalTStatisticScore(Gene gene, Annotation window) {
		
		if(!permutationScoring) {
			throw new IllegalStateException("Must enable 'do permutation scoring' option");
		}
		
		double windowScore = tStatisticWindowScores.get(gene).get(window).doubleValue();
		double numLess = 0;
		double numMore = 0;
		for(SamplePermutation perm : sampleIdentityPermutations) {
			double score = tStatisticWindowScore(gene, window, perm.getControls(), perm.getSignals());
			if(score > windowScore) numMore++;
			else numLess++;
		}
		return numMore / (numLess + numMore);
	}

	/**
	 * FDR for t statistic score of window
	 * @param gene The gene the window belongs to
	 * @param window The window
	 * @return The corrected P value
	 */
	@SuppressWarnings({ "static-method", "unused" })
	private double tStatisticWindowScoreFDR(Gene gene, Annotation window) {
		throw new UnsupportedOperationException("TODO");
	}

	public Collection<Annotation> mergePeaks(Collection<Annotation> peaks) {
		TreeSet<Annotation> peakTree = new TreeSet<Annotation>();
		peakTree.addAll(peaks);
		Collection<Annotation> mergedWindows = AnnotationUtils.mergeOverlappingBlocks(peakTree);
		for(Annotation window : mergedWindows) {
			//logger.debug("MERGED_WINDOW\t" + window.toBED());
		}
		return mergedWindows;
	}
	
	public void setCoordinateSpace(CoordinateSpace space) {
		coord = (TranscriptomeSpace) space;
	}

	public CoordinateSpace getCoordinateSpace() {
		return coord;
	}
	
	
	/**
	 * Get a systematic iterator over all possible sample identity permutations
	 * @return An iterator over all permutations of sample identity (control or signal)
	 */
	@SuppressWarnings("unused")
	private Iterator<SamplePermutation> getIterAllSampleIdentityPermutations() {
		
		if(!permutationScoring) {
			throw new IllegalStateException("Must enable 'do permutation scoring' option");
		}
		
		return new SystematicSamplePermutationIterator();
	}

	/**
	 * Get an iterator over a fixed number of random sample identity permutations
	 * @param numPermutations The number of permutations to get
	 * @return An iterator over sample identity permutations (control or signal)
	 */
	@SuppressWarnings("unused")
	private Iterator<SamplePermutation> getIterRandomSampleIdentityPermutations(int numPermutations) {
		
		if(!permutationScoring) {
			throw new IllegalStateException("Must enable 'do permutation scoring' option");
		}
		
		return new RandomSamplePermutationIterator(numPermutations);
	}
	
	/**
	 * Get an iterator over a set of sample identity permutations
	 * Either a fixed number of random permutations or all possible permutations, whichever is less
	 * @param maxNumPermutations The max number of permutations to get
	 * @return An iterator over at most the max number of sample identity permutations (control or signal)
	 */
	@SuppressWarnings("unused")
	private Iterator<SamplePermutation> getIterRandomOrAllSampleIdentityPermutations(int maxNumPermutations) {
		
		if(!permutationScoring) {
			throw new IllegalStateException("Must enable 'do permutation scoring' option");
		}
		
		SystematicSamplePermutationIterator systematicIter = new SystematicSamplePermutationIterator();
		if(systematicIter.getTotalNumPermutations() <= maxNumPermutations) {
			return systematicIter;
		}
		return new RandomSamplePermutationIterator(maxNumPermutations);
	}
	
	/**
	 * Get a random permutation of control and signal identities
	 * @return The random permutation
	 */
	protected SamplePermutation getOneRandomSamplePermutation() {
		ArrayList<Integer> controlPositions = new ArrayList<Integer>();
		int numControlsAssigned = 0;
		while(numControlsAssigned < numControls) {
			Integer randPos = Integer.valueOf(random.nextInt(numSamples));
			if(controlPositions.contains(randPos)) {
				continue;
			}
			controlPositions.add(randPos);
			numControlsAssigned++;
		}
		return new SamplePermutation(controlPositions);
	}

	/**
	 * Set logger levels for all samples
	 * @param level Level
	 */
	protected void setLoggerLevel(Level level) {
		logger.setLevel(level);
		for(SampleData sample : signalSamples) {
			sample.getLogger().setLevel(level);
		}
		for(SampleData sample : controlSamples) {
			sample.getLogger().setLevel(level);
		}
		expressionData.getLogger().setLevel(level);
	}
	
	private static CommandLineParser getCommandLineParser(String[] commandArgs) {
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-l", "Sample list file", true);
		p.addStringArg("-b", "Bed file of genes", true);
		p.addIntArg("-w", "Window size", false, DEFAULT_WINDOW_SIZE);
		p.addIntArg("-s", "Step size", false, DEFAULT_STEP_SIZE);
		p.addDoubleArg("-sp", "Scan P value cutoff for peak within gene", false, DEFAULT_PEAK_SCAN_P_VALUE_CUTOFF);
		p.addDoubleArg("-cp", "Window count cutoff for peak", false, DEFAULT_PEAK_WINDOW_COUNT_CUTOFF);
		p.addBooleanArg("-batch", "Batch out peak writing by sample name and chromosome", false, false);
		p.addStringArg("-cl", "Chromosome list file for batched run", false, null);
		p.addIntArg("-m", "Memory request for batched processes", false, DEFAULT_BATCH_MEM_REQUEST);
		p.addStringArg("-bj", "Batched peak caller jar file", false, null);
		p.addStringArg("-o", "Output directory", false, null);
		p.addDoubleArg("-q", "Quantile for peak trimming by trim max contiguous algorithm", false, DEFAULT_TRIM_PEAK_QUANTILE);
		p.addStringArg("-c", "Chromosome size file", true);
		p.addBooleanArg("-d", "Debug logging", false, false);
		p.addDoubleArg("-ep", "Scan P value cutoff for gene expression", false, DEFAULT_EXPRESSION_SCAN_P_VALUE_CUTOFF);
		p.addBooleanArg("-ft", "First read is transcription strand", false, DEFAULT_FIRST_READ_TRANSCRIPTION_STRAND);
		p.addDoubleArg("-r", "Cutoff for percentage of reads in peak coming from the most common replicate fragment", false, DEFAULT_PEAK_MAX_PCT_DUPLICATES);
		p.addBooleanArg("-sf", "Apply strand filter using read strand info", false, DEFAULT_FILTER_BY_STRAND);
		p.addBooleanArg("-ef", "Print additional info in BED file", false, DEFAULT_EXTRA_FIELDS);
		p.addBooleanArg("-binom", "Use binomial score", false,DEFAULT_USE_BINOMIAL);
		p.parse(commandArgs);
		return p;
	}
	
	protected static MultiSampleScanPeakCaller createFromCommandArgs(String[] commandArgs) throws IOException {
		CommandLineParser p = getCommandLineParser(commandArgs);
		String sampleListFile = p.getStringArg("-l");
		String bedFile = p.getStringArg("-b");
		int windowSize = p.getIntArg("-w");
		int stepSize = p.getIntArg("-s");
		double scanPvalCutoff = p.getDoubleArg("-sp");
		double trimQuantile = p.getDoubleArg("-q");
		String chrSizeFile = p.getStringArg("-c");
		double windowCountCutoff = p.getDoubleArg("-cp");
		double expressionScanPvalCutoff = p.getDoubleArg("-ep");
		boolean firstReadIsTranscriptionStrand = p.getBooleanArg("-ft");
		double maxPctMostCommonReplicatePerPeak = p.getDoubleArg("-r");
		boolean useStrandFilter = p.getBooleanArg("-sf");
		boolean extraFields =  p.getBooleanArg("-ef");
		boolean binomialScore = p.getBooleanArg("-binom");
		
		MultiSampleScanPeakCaller m = new MultiSampleScanPeakCaller(sampleListFile, bedFile, chrSizeFile, windowSize, stepSize);
		m.setExpressionScanPvalueCutoff(expressionScanPvalCutoff);
		m.setFirstReadTranscriptionStrand(firstReadIsTranscriptionStrand);
		m.setPeakCutoffMostCommonReplicate(maxPctMostCommonReplicatePerPeak);
		m.setPeakTrimQuantile(trimQuantile);
		m.setPeakWindowCountCutoff(windowCountCutoff);
		m.setPeakWindowScanPvalCutoff(scanPvalCutoff);
		m.setFilterByStrand(useStrandFilter);
		m.setExtraFields(extraFields);
		m.setBinomialScore(binomialScore);
		
		return m;
		 
	}
	
	private static String commandLineBatchChrList(String[] commandArgs) {
		CommandLineParser p = getCommandLineParser(commandArgs);
		return p.getStringArg("-cl");
	}
	
	private static int commandLineBatchMemRequest(String[] commandArgs) {
		CommandLineParser p = getCommandLineParser(commandArgs);
		return p.getIntArg("-m");
	}
	
	protected static boolean commandLineHasDebugFlag(String[] commandArgs) {
		CommandLineParser p = getCommandLineParser(commandArgs);
		return p.getBooleanArg("-d");
	}
	
	private static boolean commandLineHasBatchFlag(String[] commandArgs) {
		CommandLineParser p = getCommandLineParser(commandArgs);
		return p.getBooleanArg("-batch");
	}
	
	protected static String commandLineOutDir(String[] commandArgs) {
		CommandLineParser p = getCommandLineParser(commandArgs);
		return p.getStringArg("-o");		
	}
	
	private static String commandLineBatchJar(String[] commandArgs) {
		CommandLineParser p = getCommandLineParser(commandArgs);
		String jar = p.getStringArg("-bj");
		if(jar == null) {
			throw new IllegalArgumentException("Must provide batch peak caller jar file with option -bj.");
		}
		return jar;
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 * @throws InterruptedException 
	 * @throws DrmaaException 
	 */
	public static void main(String[] args) throws IOException, InterruptedException, DrmaaException {

		MultiSampleScanPeakCaller m = createFromCommandArgs(args);
		
		if(commandLineHasDebugFlag(args)) {
			m.setLoggerLevel(Level.DEBUG);
		}
		
		if(commandLineHasBatchFlag(args)) {
			m.batchWriteSingleSampleScanPeaksAllSamples(args, commandLineBatchChrList(args), commandLineBatchMemRequest(args));
		} else {
			m.initializeFilterRejectWriters("all_chr", FILTER_REJECT_DIR);
			m.writeSingleSampleScanPeaksAllSamples(commandLineOutDir(args));
			m.closeFilterRejectWriters();
		}
		
		
		logger.info("");
		logger.info("All done.");
		
	}

	
	/**
	 * Helper class to parse file containing list of sample types and bam file names
	 * @author prussell
	 *
	 */
	private class SampleFileParser {
		
		private String EXPRESSION_LABEL = "Expression";
		private String EXPRESSION_PVAL_CUTOFF_LABEL = "Expression_pval_cutoff";
		private String EXPRESSION_AVG_COVERAGE_CUTOFF_LABEL = "Expression_avg_coverage_cutoff";
		private String CONTROL_LABEL = "Control";
		private String SIGNAL_LABEL = "Signal";
		private String FIRST_READ_TRANSCRIPTION_STRAND_LABEL = "first_read_transcription_strand";
		private String SECOND_READ_TRANSCRIPTION_STRAND_LABEL = "second_read_transcription_strand";
		private GenomeSpaceSampleData expressionSampleData;
		private ArrayList<SampleData> controlData;
		private ArrayList<SampleData> signalData;

		
		public SampleFileParser(String sampleFileName) throws IOException {
			controlData = new ArrayList<SampleData>();
			signalData = new ArrayList<SampleData>();
			expressionSampleData = null;
			binomialCtrl = null;
			parseFile(sampleFileName);
		}
		
		/**
		 * Get the control datasets
		 * @return Set of control datasets
		 */
		public ArrayList<SampleData> getControlDatasets() {
			return controlData;
		}
		
		/**
		 * Get the signal datasets
		 * @return Set of signal datasets
		 */
		public ArrayList<SampleData> getSignalDatasets() {
			return signalData;
		}
		
		/**
		 * Get the expression dataset
		 * @return Expression dataset
		 */
		public GenomeSpaceSampleData getExpressionData() {
			return expressionSampleData;
		}
		
		/**
		 * Get the binomial control
		 * @return Binomial control dataset
		 */
		public SampleData getBinomialCtrlData() {
			return binomialCtrl;
		}
		
		/**
		 * Parse the sample file and populate data sets
		 * @throws IOException
		 */
		private void parseFile(String sampleFileName) throws IOException {
			boolean foundExpressionData = false;
			FileReader r = new FileReader(sampleFileName);
			BufferedReader b = new BufferedReader(r);
			StringParser s = new StringParser();
			
			if(!b.ready()) {
				crashWithHelpMessage();
			}
			
			String pvalLine = b.readLine();
			s.parse(pvalLine);
			boolean expByScanPval = false;
			if(s.getFieldCount() != 2) {
				crashWithHelpMessage();
			}
			if(s.asString(0).equals(EXPRESSION_AVG_COVERAGE_CUTOFF_LABEL)) {
				expByScanPval = false;
			} else if (s.asString(0).equals(EXPRESSION_PVAL_CUTOFF_LABEL)) {
				expByScanPval = true;
			} else {
				crashWithHelpMessage();
			}
			double cutoff = s.asDouble(1);
			
			if(!b.ready()) {
				crashWithHelpMessage();
			}
			
			while(b.ready()) {
				
				String line = b.readLine();
				s.parse(line);
				
				if(s.getFieldCount() == 0) continue;
				
				String label = s.asString(0);
				String bamFile = s.asString(1);

				if(label.equals(EXPRESSION_LABEL)) {
					if(s.getFieldCount() != 3) crashWithHelpMessage();
					if(foundExpressionData) {
						crashWithHelpMessage();
					}
					@SuppressWarnings("unused")
					boolean firstReadIsTranscriptionStrand = false;
					if(s.asString(2).equals(FIRST_READ_TRANSCRIPTION_STRAND_LABEL)) {
						firstReadIsTranscriptionStrand = true;
					} else {
						if(!s.asString(2).equals(SECOND_READ_TRANSCRIPTION_STRAND_LABEL)) {
							crashWithHelpMessage();
						}
					}
					logger.info("Creating sample data object for gene expression from bam file " + bamFile);
					GenomeSpaceSampleData sample = new GenomeSpaceSampleData(bamFile, sizeFile, genes, windowSize, stepSize, cutoff, true);
					expressionSampleData = sample;
					foundExpressionData = true;
					SampleData sample2 = new SampleData(bamFile, firstReadIsTranscriptionStrand, genes, windowSize, stepSize, cutoff, expByScanPval, false);
					binomialCtrl = sample2;
					continue;
				}
				
				if(label.equals(CONTROL_LABEL)) {
					if(s.getFieldCount() != 3) crashWithHelpMessage();
					logger.info("Creating sample data object for bam file " + bamFile);
					boolean firstReadIsTranscriptionStrand = false;
					if(s.asString(2).equals(FIRST_READ_TRANSCRIPTION_STRAND_LABEL)) {
						firstReadIsTranscriptionStrand = true;
					} else {
						if(!s.asString(2).equals(SECOND_READ_TRANSCRIPTION_STRAND_LABEL)) {
							crashWithHelpMessage();
						}
					}
					SampleData sample = new SampleData(bamFile, firstReadIsTranscriptionStrand, genes, windowSize, stepSize, cutoff, expByScanPval, false);
					controlData.add(sample);
					continue;
				}
				
				if(label.equals(SIGNAL_LABEL)) {
					if(s.getFieldCount() != 3) crashWithHelpMessage();
					logger.info("Creating sample data object for bam file " + bamFile);
					boolean firstReadIsTranscriptionStrand = false;
					if(s.asString(2).equals(FIRST_READ_TRANSCRIPTION_STRAND_LABEL)) {
						firstReadIsTranscriptionStrand = true;
					} else {
						if(!s.asString(2).equals(SECOND_READ_TRANSCRIPTION_STRAND_LABEL)) {
							crashWithHelpMessage();
						}
					}
					SampleData sample = new SampleData(bamFile, firstReadIsTranscriptionStrand, genes, windowSize, stepSize, cutoff, expByScanPval, false);
					signalData.add(sample);
					continue;
				}
				
				crashWithHelpMessage();
				
			}
			
			r.close();
			b.close();
			
			if(!foundExpressionData) {
				crashWithHelpMessage();
			}
			
		}
		
		/**
		 * Crash and print help message if sample file is invalid
		 */
		private void crashWithHelpMessage() {
			logger.error("");
			logger.error("**********");
			logger.error("");
			logger.error("Sample file not valid.");
			logger.error("");
			logger.error("First line must be:");
			logger.error(EXPRESSION_PVAL_CUTOFF_LABEL + "\t<pval_cutoff>");
			logger.error("-OR-");
			logger.error(EXPRESSION_AVG_COVERAGE_CUTOFF_LABEL + "\t<avg_depth_cutoff>");
			logger.error("");
			logger.error("Exactly one line must be of the form:");
			logger.error(EXPRESSION_LABEL + "\t<bam_file_name>\t" + FIRST_READ_TRANSCRIPTION_STRAND_LABEL + " OR " + SECOND_READ_TRANSCRIPTION_STRAND_LABEL);
			logger.error("");
			logger.error("Each additional line must be of the form:");
			logger.error(CONTROL_LABEL + "\t<bam_file_name>\t" + FIRST_READ_TRANSCRIPTION_STRAND_LABEL + " OR " + SECOND_READ_TRANSCRIPTION_STRAND_LABEL);
			logger.error("- or -");
			logger.error(SIGNAL_LABEL + "\t<bam_file_name>" + FIRST_READ_TRANSCRIPTION_STRAND_LABEL + " OR " + SECOND_READ_TRANSCRIPTION_STRAND_LABEL);
			logger.error("");
			logger.error("**********");
			logger.error("");
			throw new IllegalArgumentException("Sample file not valid.");
		}
		
		
	}

	
	/**
	 * Iterator over randomly generated sample identity permutations
	 * @author prussell
	 *
	 */
	private class RandomSamplePermutationIterator implements Iterator<SamplePermutation> {
		
		private int numPermutations;
		private int nextPosition;
		
		/**
		 * Construct with number of permutations to generate
		 * @param numToGenerate
		 */
		public RandomSamplePermutationIterator(int numToGenerate) {
			numPermutations = numToGenerate;
			nextPosition = 0;
		}
		
		/**
		 * Get the number of permutations to generate
		 * @return Number of permutations to generate
		 */
		@SuppressWarnings("unused")
		public int getNumToGenerate() {
			return numPermutations;
		}

		@Override
		public boolean hasNext() {
			return nextPosition < numPermutations;
		}

		@Override
		public SamplePermutation next() {
			nextPosition++;
			return getOneRandomSamplePermutation();
		}

		@Override
		public void remove() {
			throw new UnsupportedOperationException("TODO");
		}
		
	}
	
	/**
	 * Systematic iterator over all possible sample identity permutations
	 * @author prussell
	 *
	 */
	private class SystematicSamplePermutationIterator implements Iterator<SamplePermutation> {
		
		private int n;
		private int r;
		private int[] currentOneBasedPositions;
		private int[] nextOneBasedPositions;
		private int[] lastOneBasedPositions;
		
		/**
		 * Initialize first set of "control" positions to be 0,1,...,numControls
		 * Set last set to n-r,n-r+1,...,n-1
		 */
		public SystematicSamplePermutationIterator() {
			n = numSamples;
			r = numControls;
			nextOneBasedPositions = new int[r];
			for(int i = 0; i < r; i++) {
				nextOneBasedPositions[i] = algorithmI(i);
				lastOneBasedPositions[i] = algorithmI(n - r + i);
			}
		}
		
		/**
		 * Translate zero based position to one based position from combination generation algorithm
		 * @param currentPosition Zero based position
		 * @return Position in combination generation algorithm
		 */
		private int algorithmI(int currentPosition) {
			return currentPosition + 1;
		}
		
		/**
		 * Translate one based position from combination generation algorithm to zero based position
		 * @param algorithmI Position in combination generation algorithm
		 * @return Zero based position
		 */
		private int position(int algorithmI) {
			return algorithmI - 1;
		}
		
		/**
		 * Get the total number of all possible permutations
		 * @return The total number of permutations
		 */
		public long getTotalNumPermutations() {
			return MathUtil.binomialCoefficient(numSamples, numControls);
		}

		@Override
		public boolean hasNext() {
			return currentOneBasedPositions != lastOneBasedPositions;
		}

		@Override
		public SamplePermutation next() {
			
			currentOneBasedPositions = nextOneBasedPositions;
			// Increment next positions according to lexicographic order
			int k = 0;
			for(int i = 0; i < r; i++) {
				if(nextOneBasedPositions[i] < n - r + algorithmI(i)) {
					k = i;
				}
			}
			nextOneBasedPositions[k]++;
			for(int i = k + 1; i < r; i++) {
				nextOneBasedPositions[i] = nextOneBasedPositions[i-1]+1;
			}
			
			ArrayList<Integer> rtrnPositions = new ArrayList<Integer>();
			for(int i = 0; i < r; i++) {
				rtrnPositions.add(new Integer(position(currentOneBasedPositions[i])));
			}
			
			return new SamplePermutation(rtrnPositions);
			
		}

		@Override
		public void remove() {
			throw new UnsupportedOperationException("TODO");
		}
		
	}
	
	/**
	 * A permutation of sample identities (control or signal)
	 * @author prussell
	 *
	 */
	private class SamplePermutation {
		
		private Collection<SampleData> controls;
		private Collection<SampleData> signals;
		private ArrayList<Integer> controlPositionsInRealSampleList;
		private ArrayList<Integer> signalPositionsInRealSampleList;
		
		/**
		 * Construct with positions of "control" and "signal" samples from the real sample list
		 * @param labeledControlsPositionsInRealSampleList Positions in real sample list of samples to label as controls
		 */
		public SamplePermutation(ArrayList<Integer> labeledControlsPositionsInRealSampleList) {
			controlPositionsInRealSampleList = labeledControlsPositionsInRealSampleList;
			signalPositionsInRealSampleList = getSignalPositionsInRealSampleList();
			controls = new ArrayList<SampleData>();
			for(int i=0; i < controlPositionsInRealSampleList.size(); i++) {
				controls.add(allSamples.get(controlPositionsInRealSampleList.get(i).intValue()));
			}
			signals = new ArrayList<SampleData>();
			for(int i=0; i < signalPositionsInRealSampleList.size(); i++) {
				signals.add(allSamples.get(signalPositionsInRealSampleList.get(i).intValue()));
			}
		}
		
		/**
		 * Get positions of "signal" labeled samples in real sample list
		 * Based on the samples labeled "control"
		 * @return Positions of "signal" labeled samples in real sample list
		 */
		private ArrayList<Integer> getSignalPositionsInRealSampleList() {
			ArrayList<Integer> rtrn = new ArrayList<Integer>();
			for(int i=0; i<numSamples; i++) {
				if(!controlPositionsInRealSampleList.contains(Integer.valueOf(i))) {
					rtrn.add(Integer.valueOf(i));
				}
			}
			return rtrn;
		}
		
		/**
		 * Get the samples labeled as controls in this permutation
		 * @return The controls in this permutation
		 */
		public Collection<SampleData> getControls() {
			return controls;
		}
		
		/**
		 * Get the samples labeled as signals in this permutation
		 * @return The signals in this permutation
		 */
		public Collection<SampleData> getSignals() {
			return signals;
		}
		
		
	}


}
