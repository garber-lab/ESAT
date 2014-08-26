package umms.core.model.score;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import umms.core.annotation.Annotation;
import umms.core.annotation.Gene;
import umms.core.coordinatesystem.CoordinateSpace;
import umms.core.model.AlignmentModel;

/**
 * @author shari
 * Stores scores from multiple tests
 * To use additional scores, add score class to scores hash
 * p-values are stored in the other score classes included here
 * "Region" refers to an arbitrary region being scanned, like a gene, or a chromosome, or some subset of the whole coordinate space
 */
public class MultiScore extends CountScore {

	private Map<String,CountScore> scores;
	protected AlignmentModel sampleModel;
	protected AlignmentModel ctrlModel;
	private CoordinateSpace sampleCoordSpace;
	private CoordinateSpace ctrlCoordSpace;
	private double ctrlCount;              
	private double sampleRegionCount;
	private double ctrlRegionCount;
	private double regionLength;
	private Gene gene;
	private static double DEFAULT_REGION_LENGTH = -1;
	private static String BINOMIAL_KEY = "binomial";
	private static String SCAN_STAT_KEY = "scan_statistic";
	private static ArrayList<String> ALL_KEYS = new ArrayList<String>(Arrays.asList(BINOMIAL_KEY,SCAN_STAT_KEY));
	
	/**
	 * Initialize individual scores by name.
	 * To add a new score to MultiScore, add key to ALL_KEYS and add constructor here
	 * @param scoreName
	 * @return
	 */
	private CountScore getScore(String scoreName) {
		CountScore rtrn = null;
		if (scoreName.equals(BINOMIAL_KEY)) {
			rtrn = new BinomialEnrichmentScore(this);
		} else if (scoreName.equals(SCAN_STAT_KEY)) {
			rtrn = new ScanStatisticScore(this);
		}
		return rtrn;
	}
	
	/**
	 * Constructs score without region/gene info
	 * @param sample sample alignment data
	 * @param ctrl control alignment data
	 * @param a annotation window/gene/etc being scored
	 * @param regionLength length of gene
	 * @param fullyContained only count reads fully contained in window
	 */
	public MultiScore(AlignmentModel sample, AlignmentModel ctrl, Annotation a, double regionLength, boolean fullyContained) {
		super(sample, a, fullyContained);
		setSampleCoordSpace(sample.getCoordinateSpace());
		ctrlCoordSpace = ctrl.getCoordinateSpace();
		this.sampleModel = sample;
		this.ctrlModel = ctrl;
		scores = new HashMap<String,CountScore>();
		
		if (!(getSampleCoordSpace().getChromosomeNames().equals(ctrlCoordSpace.getChromosomeNames()) & getSampleCoordSpace().getLength() == getCtrlCoordSpace().getLength())) {
			throw new IllegalArgumentException("Sample coordinate space must match control coordinate space");
		}
		
		gene = new Gene(a);
		setRegionLength(regionLength);
		setCtrlCount(ctrl.getCount(a,fullyContained));
		setSampleRegionCount(getCount());
		setCtrlRegionCount(getCtrlCount());
		
		initializeScores();
		
	}
	
	public MultiScore(AlignmentModel sample, AlignmentModel ctrl, Annotation a, boolean fullyContained) {
		this(sample,ctrl,a,DEFAULT_REGION_LENGTH,fullyContained);
	}
	
	public MultiScore(AlignmentModel sample, AlignmentModel ctrl, Annotation a, double regionLength) {
		this(sample,ctrl,a,regionLength,false);
	}
	
	public MultiScore(AlignmentModel sample, AlignmentModel ctrl, Annotation a) {
		this(sample,ctrl,a,false);
	}
	
	// Constructor with gene information
	public MultiScore(AlignmentModel sample, AlignmentModel ctrl, Annotation a, Gene gene, boolean fullyContained) {
		this(sample,ctrl,a,gene.getSize(),fullyContained);
		this.gene = gene;
		if (!gene.geneSpanContains(new Gene(a))) {
			throw new IllegalArgumentException("Gene must contain annotation window");
		}
		setCtrlRegionCount(ctrl.getCount(gene,fullyContained));
		setSampleRegionCount(sample.getCount(gene,fullyContained));
		setRegionLength(gene.getSize());
		
		updateScores();
	}
	
	public MultiScore(AlignmentModel sample, AlignmentModel ctrl, Annotation a, Gene gene) {
		this(sample,ctrl,a,gene,false);
	}
	
	// If only one alignment model given, only computes scan statistic score
	public MultiScore(AlignmentModel sample, Annotation a, boolean fullyContained) {
		super(sample, a, fullyContained);
		setSampleCoordSpace(sample.getCoordinateSpace());
		this.sampleModel = sample;
		this.gene = new Gene(a);
		setRegionLength(sampleModel.getGlobalLength());
		setSampleRegionCount(getTotal());
		
		initializeScores(SCAN_STAT_KEY);
	}
	
	// Constructors with region/gene counts given (faster than passing gene if counts are precomputed)
	public MultiScore(AlignmentModel sample, AlignmentModel ctrl, Annotation a, double sampleRegCount, double ctrlRegCount) {
		this(sample,ctrl,a,sampleRegCount,ctrlRegCount,DEFAULT_REGION_LENGTH,false);
	}

	public MultiScore(AlignmentModel sample, AlignmentModel ctrl, Annotation a, double sampleRegCount, double ctrlRegCount, double regionLength) {
		this(sample,ctrl,a,sampleRegCount,ctrlRegCount,regionLength,false);
	}
	
	public MultiScore(AlignmentModel sample, AlignmentModel ctrl, Annotation a, double sampleRegCount, double ctrlRegCount, double regionLength, boolean fullyContained) {
		this(sample,ctrl,a,regionLength,fullyContained);
		setSampleRegionCount(sampleRegCount);
		setCtrlRegionCount(ctrlRegCount);
		updateScores();
	}
	
	// Constructors using previous score
	public MultiScore(AlignmentModel sample, AlignmentModel ctrl, Annotation annotation, MultiScore previousScore, double newSampleCount,double newCtrlCount) {
		super(previousScore, annotation, newSampleCount); //Set the new score without computing
		sampleModel = sample;
		ctrlModel = ctrl;
		setSampleCoordSpace(sample.getCoordinateSpace());
		ctrlCount = newCtrlCount;
		this.sampleRegionCount = previousScore.getSampleRegionCount();
		this.ctrlRegionCount = previousScore.getCtrlRegionCount();
		this.regionLength = previousScore.getRegionLength();
		this.gene = previousScore.getGene();
		initializeScores();
	}
	
	private void initializeScores() {
		scores = new HashMap<String,CountScore>();
		for (String scoreName : ALL_KEYS) {
			CountScore score = getScore(scoreName);
			scores.put(scoreName, score);
		}
	}
	
	private void initializeScores(String scoreName) {
		scores = new HashMap<String,CountScore>();
		scores.put(scoreName, getScore(scoreName));
	}
	
	public void updateScores() {
		for (String scoreName : ALL_KEYS) {
			CountScore score = getScore(scoreName);
			scores.put(scoreName, score);
		}
	}
	
	public void updateScores(String scoreName) {
		scores.put(scoreName, getScore(scoreName));
	}
	
	public double getEnrichmentOverControl() {
		double normSampleCounts = getCount()/getSampleRegionCount();
		double normCtrlCounts = getCtrlCount()/getCtrlRegionCount();
		return normSampleCounts/normCtrlCounts;
	}
	
	public double getAverageCoverage(AlignmentModel model) {
		ScanStatisticScore score = (ScanStatisticScore) scores.get(SCAN_STAT_KEY);
		return score.getAverageCoverage(model);
	}
		
	public static class Processor extends WindowProcessor.AbstractProcessor<MultiScore> {
		protected AlignmentModel sample;
		protected AlignmentModel ctrl;
		protected double sampleRegionCount = DEFAULT_REGION_TOTAL;
		protected double ctrlRegionCount = DEFAULT_REGION_TOTAL;
		protected double regionLength = DEFAULT_REGION_LENGTH;
		private boolean fullyContainedReads;
		private Gene gene;
		
		public Processor(AlignmentModel sample, AlignmentModel ctrl, boolean fullyContained) {
			this.sample = sample;
			this.ctrl = ctrl;
			this.fullyContainedReads = fullyContained;
		}
		
		public Processor(AlignmentModel sample, AlignmentModel ctrl) {
			this(sample,ctrl,false);
		}
		
		public Processor(AlignmentModel sample) {
			this.sample = sample;
			this.ctrl = sample;
		}
		
		public void initRegion(Annotation a){
			if (a != null){
				sampleRegionCount = sample.getCount(a);
				ctrlRegionCount = ctrl.getCount(a);
				regionLength = a.size();
			}
		}
		
		public MultiScore processWindow(Annotation a){
			return new MultiScore(sample, ctrl, a, sampleRegionCount, ctrlRegionCount, regionLength, fullyContainedReads);
		}
		
		/**
		 * Compute the count using the previous windowScore
		 * @param nextRegion 
		 * @param previousScore The WindowScore before
		 * @return the count of the current window
		 */
		private double computeSampleCount(Annotation nextRegion, MultiScore previousScore) {
			//else, get the minus region scores and the plus value scores
			//This is not so simple because we'll need to use the fully contained regions
			double subtractVal=sample.getCountExcludingRegion(previousScore.getAnnotation().minus(nextRegion), nextRegion);
			double addVal=sample.getCountExcludingRegion(nextRegion.minus(previousScore.getAnnotation()), previousScore.getAnnotation());
			return (previousScore.getSampleCount()-subtractVal)+addVal;
		}
		
		private double computeCtrlCount(Annotation nextRegion, MultiScore previousScore) {
			//else, get the minus region scores and the plus value scores
			//This is not so simple because we'll need to use the fully contained regions
			double subtractVal=ctrl.getCountExcludingRegion(previousScore.getAnnotation().minus(nextRegion), nextRegion);
			double addVal=ctrl.getCountExcludingRegion(nextRegion.minus(previousScore.getAnnotation()), previousScore.getAnnotation());
			return (previousScore.getCtrlCount()-subtractVal)+addVal;
		}
		
		@Override
		public MultiScore processWindow(Annotation annotation, MultiScore previousScore) {
			//if the previous score is null or they don't overlap
			if(previousScore==null || !annotation.overlaps(previousScore.getAnnotation())){
				//compute the score directly
				return processWindow(annotation);
			}
			double newSampleCount=computeSampleCount(annotation, previousScore);
			double newCtrlCount=computeCtrlCount(annotation, previousScore);
			return new MultiScore(sample, ctrl, annotation, previousScore, newSampleCount, newCtrlCount);
			}
		
	}
	
	public void getDebugInfo(String prefix) {
		logger.debug(prefix + " " + gene.getName()+" "+annotation.toBED()+" ");
		logger.debug(prefix + " sample_count="+Double.toString(getSampleCount()) + " ctrl_count="+Double.toString(getCtrlCount()));
		logger.debug(prefix + " sample_region_count="+Double.toString(getSampleRegionCount())+"  ctrl_region_count="+Double.toString(getCtrlRegionCount()));
		logger.debug(prefix + " window_length="+Double.toString(annotation.getSize())+"  region_length="+Double.toString(getRegionLength()));
		logger.debug(prefix + " binomial_p="+Double.toString(scores.get(BINOMIAL_KEY).getPvalue())+" scan_p="+Double.toString(scores.get(SCAN_STAT_KEY).getPvalue()));
	}
	
	public double getPvalue(String scoreName) { return scores.get(scoreName).getPvalue(); }
	
	public void setSampleCount(double d) { setCount(d); }
	public void setCtrlCount(double d) { ctrlCount = d; }

	public void setSampleRegionCount(double d) { sampleRegionCount = d; }
	public void setCtrlRegionCount(double d) { ctrlRegionCount = d; }
	
	public double getSampleCount() { return getCount(); }
	public double getCtrlCount() { return ctrlCount; }

	public double getCtrlRegionCount() { return ctrlRegionCount; }
	public double getSampleRegionCount() { return sampleRegionCount; }
	
	public void setRegionLength(double d) { regionLength = d; }
	public double getRegionLength() { return regionLength; }
	
	public void setGene(Gene g) { gene = g; }
	public Gene getGene() { return gene; }

	public CoordinateSpace getSampleCoordSpace() { return sampleCoordSpace; }
	public void setSampleCoordSpace(CoordinateSpace sampleCoordSpace) { this.sampleCoordSpace = sampleCoordSpace; }

	public CoordinateSpace getCtrlCoordSpace() { return ctrlCoordSpace; }
	public void setCtrlCoordSpace(CoordinateSpace sampleCoordSpace) { this.ctrlCoordSpace = sampleCoordSpace; }

	public double getGlobalLength() { 
		ScanStatisticScore score = (ScanStatisticScore) scores.get(SCAN_STAT_KEY);
		return score.getGlobalLength();
	}
	
	public double getGlobalLambda() { 
		ScanStatisticScore score = (ScanStatisticScore) scores.get(SCAN_STAT_KEY);
		return score.getGlobalLambda();
	}
	
}
