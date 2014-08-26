package umms.core.model.score;

import broad.pda.seq.segmentation.AlignmentDataModelStats;
import net.sf.samtools.util.CloseableIterator;
import umms.core.alignment.Alignment;
import umms.core.annotation.Annotation;
import umms.core.annotation.Gene;
import umms.core.coordinatesystem.CoordinateSpace;
import umms.core.model.AlignmentModel;
import jsc.distributions.Binomial;

public class BinomialEnrichmentScore extends CountScore {

	private double Pvalue;
	private CoordinateSpace sampleCoordSpace;
	private CoordinateSpace ctrlCoordSpace;
	private double ctrlCount;  // sampleCount is count            
	private double sampleRegionCount;
	private double ctrlRegionCount; //Total bases in gene/region/chrom in control
	private double regionLength;
	private static double DEFAULT_REGION_LENGTH = -1;
	
	/**
	 * Binomial score for windows
	 * P-Value statistics are stored here
	 * "Region" refers to an arbitrary region being scanned, like a gene, or a chromosome, or some subset of the whole coordinate space
	 */

	public BinomialEnrichmentScore(AlignmentModel sample, AlignmentModel ctrl, Annotation a, boolean fullyContained) {
		super(sample, a, fullyContained);
		sampleCoordSpace = sample.getCoordinateSpace();
		ctrlCoordSpace = ctrl.getCoordinateSpace();
		
		if (!sampleCoordSpace.getChromosomeNames().equals(ctrlCoordSpace.getChromosomeNames())) {
			throw new IllegalArgumentException("Sample coordinate space must match control coordinate space");
		}
		
		setCtrlCount(ctrl.getCount(a,fullyContained));
		setCtrlRegionCount(ctrl.getCount(a,fullyContained));
		setSampleRegionCount(sample.getCount(a,fullyContained));
		setRegionLength(DEFAULT_REGION_LENGTH);
	
		try {
			setPvalue(calculatePVal(getSampleCount(), getCtrlCount(), getSampleRegionCount(), getCtrlRegionCount()));
			} catch(Exception e){
				logger.info("Cound not set P value for annotation " + a.getName());
			}
			getAnnotation().setScore(getPvalue());
	}
	
	public BinomialEnrichmentScore(AlignmentModel sample, AlignmentModel ctrl, Annotation a, double regionLength, boolean fullyContained) {
		this(sample,ctrl,a,fullyContained);
		setRegionLength(regionLength);
		refreshPvalue();
		
	}
	
	public BinomialEnrichmentScore(AlignmentModel sample, AlignmentModel ctrl, Annotation a, double regionLength) {
		this(sample,ctrl,a,false);
		setRegionLength(regionLength);
		refreshPvalue();
		
	}
	
	public BinomialEnrichmentScore(AlignmentModel sample, Annotation a) {
		this(sample, sample, a);
	}
	
	public BinomialEnrichmentScore(AlignmentModel sample, AlignmentModel ctrl, Annotation a) {
		this(sample,ctrl,a,false);
	}
	
	
	public BinomialEnrichmentScore(AlignmentModel sample, AlignmentModel ctrl, Annotation a, double sampleRegCount, double ctrlRegCount) {
		this(sample,ctrl,a,sampleRegCount,ctrlRegCount,false);
	}

	public BinomialEnrichmentScore(AlignmentModel sample, AlignmentModel ctrl, Annotation a, double sampleRegCount, double ctrlRegCount, boolean fullyContained) {
		this(sample,ctrl,a,fullyContained);
		setSampleRegionCount(sampleRegCount);
		setCtrlRegionCount(ctrlRegCount);
		setPvalue(calculatePVal(getSampleCount(), getCtrlCount(), getSampleRegionCount(), getCtrlRegionCount()));
		
		getAnnotation().setScore(getPvalue());
	}
	
	public BinomialEnrichmentScore(AlignmentModel sample, AlignmentModel ctrl, Annotation annotation, BinomialEnrichmentScore previousScore, double newSampleCount,double newCtrlCount) {
		super(previousScore, annotation, newSampleCount); //Set the new score without computing
		sampleCoordSpace = sample.getCoordinateSpace();
		ctrlCount = newCtrlCount;
		this.sampleRegionCount = previousScore.getSampleRegionCount();
		this.ctrlRegionCount = previousScore.getCtrlRegionCount();
		try {
			setPvalue(calculatePVal(new Double(getSampleCount()), new Double(getCtrlCount()), sampleRegionCount, ctrlRegionCount));
		} catch(Exception e) {
			logger.info("Could not set scan P value for annotation " + annotation.getName());
		}
		getAnnotation().setScore(getPvalue());
	}
	
	public BinomialEnrichmentScore(MultiScore other) {
		super(other);
		sampleCoordSpace = other.getSampleCoordSpace();
		ctrlCoordSpace = other.getCtrlCoordSpace();
		setCtrlCount(other.getCtrlCount());
		setSampleRegionCount(other.getSampleRegionCount());
		setCtrlRegionCount(other.getCtrlRegionCount());
		setRegionLength(other.getRegionLength());
		try {
			setPvalue(calculatePVal(new Double(getSampleCount()), new Double(getCtrlCount()), sampleRegionCount, ctrlRegionCount, regionLength, annotation.getSize()));
		} catch(Exception e) {
			logger.info("Could not set scan P value for annotation " + annotation.getName());
			logger.info("Sample count " + getSampleCount() + " control count " + getCtrlCount() + " sample region count " + sampleRegionCount + " control region count " + ctrlRegionCount + " region length " + regionLength + " annotation size " + annotation.getSize());
		}
	}
	
	/**
	 * @param a				Sample reads
	 * @param b	        	Control reads
	 * @param sampleRegionCounts	Sample reads in region
	 * @param ctrlRegionCounts	Control reads in region
	 * @return				Binomial P(X>a)
	 */
	public double calculatePVal(double a, double b, double sampleRegionCounts, double ctrlRegionCounts) {
		if (a+b<=2|ctrlRegionCounts+sampleRegionCounts<2) {return 1;}
		double p = sampleRegionCounts/(sampleRegionCounts+ctrlRegionCounts);
		if (p==0) {return 1;}
		long n = (long) (a + b);
		Binomial C = new Binomial(n,p);
		double pval = 1 - C.cdf(a);
		return pval;

	}
	
	public double calculatePVal(double a, double b, double sampleRegionCounts, double ctrlRegionCounts, double regionLength, double windowSize) {
		double pval1 = calculatePVal(a,b,sampleRegionCounts,ctrlRegionCounts);
		if (regionLength == DEFAULT_REGION_LENGTH | regionLength <= windowSize) {
			return pval1;
		} else {
			double p = windowSize/regionLength;
			if (p==0) {return 1;}
			long n = (long) sampleRegionCounts;
			Binomial C = new Binomial(n,p);
			double pval2 = 1 - C.cdf(a);
			return Math.max(pval1, pval2);
		}
	}
	
	@Override
	public void refreshPvalue() {
		if (regionLength != DEFAULT_REGION_LENGTH) {
			setPvalue(calculatePVal(getCount(),ctrlCount,sampleRegionCount,ctrlRegionCount,regionLength,annotation.getSize()));
		} else {
			setPvalue(calculatePVal(getCount(),ctrlCount,sampleRegionCount,ctrlRegionCount));
		}
	}
	
	public double getEnrichmentOverControl() {
		double normSampleCounts = getCount()/getSampleRegionCount();
		double normCtrlCounts = getCtrlCount()/getCtrlRegionCount();
		return normSampleCounts/normCtrlCounts;
	}
	
	public static class Processor extends WindowProcessor.AbstractProcessor<BinomialEnrichmentScore> {
		protected AlignmentModel sample;
		protected AlignmentModel ctrl;
		protected double sampleRegionCount = DEFAULT_REGION_TOTAL;
		protected double ctrlRegionCount = DEFAULT_REGION_TOTAL;
		private boolean fullyContainedReads;
		
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
			}
		}
		
		public BinomialEnrichmentScore processWindow(Annotation a){
			return new BinomialEnrichmentScore(sample, ctrl, a, sampleRegionCount, ctrlRegionCount, fullyContainedReads);
		}
		
		/**
		 * Compute the count using the previous windowScore
		 * @param nextRegion 
		 * @param previousScore The WindowScore before
		 * @return the count of the current window
		 */
		private double computeSampleCount(Annotation nextRegion, BinomialEnrichmentScore previousScore) {
			//else, get the minus region scores and the plus value scores
			//This is not so simple because we'll need to use the fully contained regions
			double subtractVal=sample.getCountExcludingRegion(previousScore.getAnnotation().minus(nextRegion), nextRegion);
			double addVal=sample.getCountExcludingRegion(nextRegion.minus(previousScore.getAnnotation()), previousScore.getAnnotation());
			return (previousScore.getSampleCount()-subtractVal)+addVal;
		}
		
		private double computeCtrlCount(Annotation nextRegion, BinomialEnrichmentScore previousScore) {
			//else, get the minus region scores and the plus value scores
			//This is not so simple because we'll need to use the fully contained regions
			double subtractVal=ctrl.getCountExcludingRegion(previousScore.getAnnotation().minus(nextRegion), nextRegion);
			double addVal=ctrl.getCountExcludingRegion(nextRegion.minus(previousScore.getAnnotation()), previousScore.getAnnotation());
			return (previousScore.getCtrlCount()-subtractVal)+addVal;
		}
		
		@Override
		public BinomialEnrichmentScore processWindow(Annotation annotation, BinomialEnrichmentScore previousScore) {
			//if the previous score is null or they don't overlap
			if(previousScore==null || !annotation.overlaps(previousScore.getAnnotation())){
				//compute the score directly
				return processWindow(annotation);
			}
				
			double newSampleCount=computeSampleCount(annotation, previousScore);
			double newCtrlCount=computeCtrlCount(annotation, previousScore);
			return new BinomialEnrichmentScore(sample, ctrl, annotation, previousScore, newSampleCount, newCtrlCount);
			}
		
	}
	
	public void setPvalue(double scanPvalue) { this.Pvalue = scanPvalue; }
	@Override
	public double getPvalue() { return Pvalue; }

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

}
