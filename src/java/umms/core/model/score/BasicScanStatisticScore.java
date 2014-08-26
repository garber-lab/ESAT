package umms.core.model.score;

import broad.core.math.ScanStatistics;
import net.sf.samtools.util.CloseableIterator;
import umms.core.alignment.Alignment;
import umms.core.annotation.Annotation;
import umms.core.annotation.Gene;
import umms.core.coordinatesystem.CoordinateSpace;
import umms.core.model.AlignmentModel;

public class BasicScanStatisticScore extends CountScore {

	/**
	 * Global statistics are stored in the CountScore class
	 *
	 * P-Value statistics are stored here
	 * "Region" refers to an arbitrary region being scanned, like a gene, or a chromosome, or some subset of the whole coordinate space
	 */
	private double scanPvalue;
	private double regionLength; // Length of region for local p-value calculation
	private double fullyContainedNumberOfReads;
	private double globalLength;
	private double windowLength;

	public BasicScanStatisticScore(AlignmentModel model, Annotation annotation,double score,boolean global,double regionLength,double regionScore) {
		super(annotation,score,regionScore,model.getGlobalCount());
		setGlobalLength(model.getGlobalLength());
		setRegionLength(regionLength);
		setWindowLength(model.getCoordinateSpace().getSize(annotation));
		
		try {
			if(global)
				setScanPvalue(calculateGlobalPvalue(model));
			else
				setScanPvalue(calculateLocalPvalue());
		} catch(Exception e) {
			logger.info("Could not set scan P value for annotation " + annotation.getName());
			logger.info(e.toString());
		}
		getAnnotation().setScore(getScanPvalue());
	}
	
	/**
	 * BAD constructor
	 * @param model
	 * @param annotation
	 * @param regionTotal
	 * @param regionLength2
	 * @param fullyContainedReads
	 */
	public BasicScanStatisticScore(AlignmentModel model, Annotation annotation,
			double regionTotal, double regionLength2,
			boolean fullyContainedReads) {
		this(model,annotation,regionTotal,true,regionLength2,regionTotal);
		
	}

	/**
	 * BAD CONSTRUCTOR
	 * @param model
	 * @param annotation
	 * @param previousScore
	 * @param newScore
	 */
	public BasicScanStatisticScore(AlignmentModel model, Annotation annotation,
			BasicScanStatisticScore previousScore, double newScore) {
		this(model, annotation, previousScore.getScore(), annotation.getSize(), false);
	}

	public double calculateGlobalPvalue(AlignmentModel model){
		return ScanStatistics.calculatePVal(new Double(getCount()).intValue(), model.getGlobalLambda(), getWindowLength(), getGlobalLength());
	}
	
	public double calculateLocalPvalue(){
		return ScanStatistics.calculatePVal(new Double(getCount()).intValue(), getLocalLambda(), getWindowLength(), getRegionLength());
	}
		
	public double getAverageCoverage() { 

		double avgCoverage = (double) getCount() / (double) getWindowLength();
		return avgCoverage;
	}
	
	public double getWindowLength(){
		return windowLength;
	}

	public void setWindowLength(double length){
		this.windowLength = length;
	}
	
	public double getGlobalEnrichment() {
		return getAverageCoverage() / getGlobalLambda();
	}

	public double getLocalLambda() {
		return (getRegionTotal()+1) / getRegionLength();
	}

	public void setGlobalLength(double d) {
		globalLength = d;
	}
	
	public double getGlobalLambda() {
		return getTotal() / getGlobalLength();
	}
	
	public double getGlobalLength() {
		return globalLength;
	}
	
	public void setRegionLength(double regionLength) {
		this.regionLength = regionLength;
	}

	public double getRegionLength() {
		return regionLength;
	}

	public void setScanPvalue(double scanPvalue) {
		this.scanPvalue = scanPvalue;
	}
	
	public void setPvalue(double scanPvalue) {
		this.scanPvalue = scanPvalue;
	}

	public double getScanPvalue() {
		return scanPvalue;
	}
	
	@Override
	public double getPvalue() {
		return scanPvalue;
	}

	public void setFullyContainedNumberOfReads(double fullyContainedNumberOfReads) {
		this.fullyContainedNumberOfReads = fullyContainedNumberOfReads;
	}

	public double getFullyContainedNumberOfReads() {
		return fullyContainedNumberOfReads;
	}
	
	/*
	 * Returns an array of scores
	 * [0] = count
	 * [1] = RPKM
	 * [2] = RPK
	 * [3] = region total
	 * [4] = total
	 * [5] = Annotation length
	 * [6] = scan p value
	 * [7] = global lambda
	 * [8] = local lambda
	 * [9] = region length  
	 */
	public double[] getScores(){
		double[] scores = new double[10];
		scores[0] = getCount();
		scores[1] = getRPKM();
		scores[2] = getRPK();
		scores[3] = getRegionTotal();
		scores[4] = getTotal();
		scores[5] = getAnnotation().length();
		scores[6] = getScanPvalue();
		scores[7] = getGlobalLambda();
		scores[8] = getLocalLambda();
		scores[9] = getRegionLength();
		return scores;
	}
	
	public String toString() {
		return super.toString() + "\t" + getScanPvalue() + "\t" 
		+ getGlobalLambda() + "\t" + getLocalLambda() + "\t" + getRegionLength();
	}
	
	
	public static class Processor extends WindowProcessor.AbstractProcessor<BasicScanStatisticScore> {
		protected AlignmentModel model;
		protected double regionTotal = DEFAULT_REGION_TOTAL;
		protected double regionLength = DEFAULT_REGION_TOTAL;
		private boolean fullyContainedReads;
		
		public Processor(AlignmentModel model) {
			this(model, false);
		}
		
		public Processor(AlignmentModel model, boolean fullyContained) {
			this.model = model;
			this.fullyContainedReads = fullyContained;
		}
		
		public BasicScanStatisticScore processWindow(Annotation annotation) {
			logger.info("Processor in BasicScanStatisticScore is called");
			return new BasicScanStatisticScore(model, annotation, regionTotal, regionLength, fullyContainedReads);
		}
		
		public void initRegion(Annotation region) {
			if (region != null) {
				regionTotal = model.getCount(region);
				regionLength = region.length();
			}
		}

		/**
		 * Compute the count using the previous windowScore
		 * @param nextRegion 
		 * @param previousScore The WindowScore before
		 * @return the count of the current window
		 */
		private double computeCount(Annotation nextRegion, CountScore previousScore) {
			//else, get the minus region scores and the plus value scores
			//This is not so simple because we'll need to use the fully contained regions
			double subtractVal=model.getCountExcludingRegion(previousScore.getAnnotation().minus(nextRegion), nextRegion);
			double addVal=model.getCountExcludingRegion(nextRegion.minus(previousScore.getAnnotation()), previousScore.getAnnotation());
			return (previousScore.getCount()-subtractVal)+addVal;
		}
		
		@Override
		public BasicScanStatisticScore processWindow(Annotation annotation, BasicScanStatisticScore previousScore) {
			//if the previous score is null or they don't overlap
			if(previousScore==null || !annotation.overlaps(previousScore.getAnnotation())){
				//compute the score directly
				return processWindow(annotation);
			}
			
			double newScore=computeCount(annotation, previousScore);
			return new BasicScanStatisticScore(model, annotation, previousScore, newScore);
		}
	}

	/**
	 * True iff scan P value and annotation are equal
	 */
	@Override
	public boolean equals(Object o) {
		BasicScanStatisticScore otherScore = (BasicScanStatisticScore) o;
		if(scanPvalue != otherScore.getScanPvalue()) return false;
		if(!getAnnotation().equals(otherScore.getAnnotation())) return false;
		return true;
	}
	

	
}
