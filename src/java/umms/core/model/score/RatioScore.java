package umms.core.model.score;

import java.util.Iterator;

import umms.core.annotation.Annotation;
import umms.core.annotation.AnnotationCollection;

import org.apache.commons.math3.distribution.BinomialDistribution;

public class RatioScore extends CountScore {
	
	public static double RPKM_OFFSET = 0.001; // reads per Kilobase per Million
	
	// the following variables are for the denominator
	private double dCount;               // use a double to allow for weighted read counting
	private double dTotal;                // whole "genome"
	private double dRegionTotal = DEFAULT_REGION_TOTAL;          // whole chromosome / gene / region under consideration / etc.

	public RatioScore(Annotation a) {
		super(a);
	}
	
	public RatioScore(Annotation numerator, Annotation denominator, double nTotal, double dTotal) {
		this(numerator, denominator, nTotal, dTotal, DEFAULT_REGION_TOTAL, DEFAULT_REGION_TOTAL);
	}
	
	public RatioScore(Annotation numerator, Annotation denominator, double nTotal, double dTotal, double nRegionTotal, double dRegionTotal) {
		super(numerator, numerator.getScore(), nTotal, nRegionTotal);
		if (!numerator.equals(denominator)) {
			throw new IllegalArgumentException("Cannot create RatioScore from numerator and denominator with different coordinates.");
		}
		dCount = denominator.getScore();
		this.dTotal = dTotal;
		this.dRegionTotal = dRegionTotal;
	}
	
	public RatioScore(CountScore numerator, CountScore denominator) {
		super(numerator);
		if (!numerator.getAnnotation().equals(denominator.getAnnotation())) {
			throw new IllegalArgumentException("Cannot create RatioScore from numerator and denominator with different coordinates.");
		}
		dCount = denominator.getCount();
		dTotal = denominator.getTotal();
		dRegionTotal = denominator.getRegionTotal();
	}
	
	public RatioScore(AnnotationCollection<? extends Annotation> numerator, AnnotationCollection<? extends Annotation> denominator, Annotation annotation) {
		super(numerator, annotation);
		setDenominatorCount(denominator.getCount(annotation));
		setDenominatorTotal(denominator.getGlobalCount());
		getAnnotation().setScore(getRatio());
	}
	
	public RatioScore(AnnotationCollection<? extends Annotation> numerator, AnnotationCollection<? extends Annotation> denominator, Annotation annotation, double numeratorRegionTotal, double denominatorRegionTotal) {
		this(numerator, denominator, annotation);
		setNumeratorRegionTotal(numeratorRegionTotal);
		setDenominatorRegionTotal(denominatorRegionTotal);
	}
	
	public void setNumeratorCount(double d) { setCount(d); }
	public void setDenominatorCount(double d) { dCount = d; }
	public void setNumeratorTotal(double d) { setTotal(d); }
	public void setDenominatorTotal(double d) { dTotal = d; }
	public void setNumeratorRegionTotal(double d) { setRegionTotal(d); }
	public void setDenominatorRegionTotal(double d) { dRegionTotal = d; }
	
	public double getNumeratorCount() { return getCount(); }
	public double getDenominatorCount() { return dCount; }
	public double getNumeratorTotal() { return getTotal(); }
	public double getDenominatorTotal() { return dTotal; }
	public double getNumeratorRegionTotal() { return getRegionTotal(); }
	public double getDenominatorRegionTotal() { return dRegionTotal; }
	
	public double getNumeratorRPKM() { return CountScore.asRPKM(getNumeratorCount(), getNumeratorTotal(), getAnnotation().size()); }
	public double getDenominatorRPKM() { return CountScore.asRPKM(getDenominatorCount(), getDenominatorTotal(), getAnnotation().size()); }
	
	public double getRatio() {
		return (getNumeratorRPKM() + RPKM_OFFSET) / (getDenominatorRPKM() + RPKM_OFFSET);
	}
	
	@Override
	public double getScore() { 
		return getRatio();
	}
	
	/**
	 * Binomial test: combine numerator and denominator counts, then calculate cumulative probability of seeing >= this proportion
	 * @return -log10 p-value for the numerator having more reads than expected.
	 */
	public double getEnrichmentBinomialScore() {
		double p=getEnrichmentBinomialPValue();
		return -1.0 * Math.log10(p);
	}
	
	/**
	 * Binomial test: combine numerator and denominator counts, then calculate cumulative probability of seeing >= this proportion
	 * @return p-value for the numerator having more reads than expected
	 */
	public double getEnrichmentBinomialPValue(){
		int total = (int) Math.round(getNumeratorCount() + getDenominatorCount());
		double p = getNumeratorTotal() / (getNumeratorTotal() + getDenominatorTotal());
		BinomialDistribution b = new BinomialDistribution(total, p);
		double pGreater = 1 - b.cumulativeProbability((int) Math.round(getNumeratorCount()));
		double pEqual = b.probability((int) Math.round(getNumeratorCount()));
		return pGreater + pEqual;
	}
	
	/**
	 * Binomial test: combine numerator and denominator counts, then calculate cumulative probability of seeing less than this proportion
	 * @return -log10 p-value for the numerator having less reads than expected.
	 */
	public double getDepletionBinomialScore() {
		double p=getDepletionBinomialPValue();
		return -1.0 * Math.log10(p);
	}
	
	/**
	 * Binomial test: combine numerator and denominator counts, then calculate cumulative probability of seeing less than this proportion
	 * @return p-value for the numerator having less reads than expected.
	 */
	public double getDepletionBinomialPValue(){
		int total = (int) Math.round(getNumeratorCount() + getDenominatorCount());
		double p = getNumeratorTotal() / (getNumeratorTotal() + getDenominatorTotal());
		BinomialDistribution b = new BinomialDistribution(total, p);
		return b.cumulativeProbability((int) Math.round(getNumeratorCount()));
	}
	
	public double getLog2Ratio() {
		return Math.log(getRatio()) / Math.log(2);
	}
	
	public double getRegionRatio() {
		return (getRegionNumeratorRPKM() + RPKM_OFFSET) / (getRegionDenominatorRPKM() + RPKM_OFFSET);
	}
	
	public double getLog2RegionRatio() {
		return Math.log(getRegionRatio()) / Math.log(2);
	}
	
	public double getRegionNumeratorRPKM() { return CountScore.asRPKM(getNumeratorCount(), getNumeratorRegionTotal(), getAnnotation().size()); }
	public double getRegionDenominatorRPKM() { return CountScore.asRPKM(getDenominatorCount(), getDenominatorRegionTotal(), getAnnotation().size()); }

	
	public String toString() {
		annotation.setScore(getScore());
		return annotation.toBED() + "\t" + getRatio() + "\t" + getLog2Ratio() + "\t" + getNumeratorCount() + "\t" + getNumeratorRPKM() + 
				"\t" + getNumeratorRegionTotal() + "\t" + getNumeratorTotal() + "\t" + getDenominatorCount() + "\t" + getDenominatorRPKM() + 
				"\t" + getDenominatorRegionTotal() + "\t" + getDenominatorTotal() + "\t" + getEnrichmentBinomialScore() + "\t" + getDepletionBinomialScore();
	}
	
	public static class Processor extends WindowProcessor.AbstractProcessor<RatioScore> {
		protected AnnotationCollection<? extends Annotation> numerator, denominator;
		protected double numeratorRegionTotal = DEFAULT_REGION_TOTAL;
		protected double denominatorRegionTotal = DEFAULT_REGION_TOTAL;
		
		public Processor(AnnotationCollection<? extends Annotation> numerator, AnnotationCollection<? extends Annotation> denominator) {
			this.numerator = numerator;
			this.denominator = denominator;
		}
		
		public RatioScore processWindow(Annotation annotation) {
			return new RatioScore(numerator, denominator, annotation, numeratorRegionTotal, denominatorRegionTotal);
		}
		
		public void initRegion(Annotation region) {
			if (region != null) {
				numeratorRegionTotal = numerator.getCount(region);
				denominatorRegionTotal = denominator.getCount(region);
			}
		}

		public double getNumeratorRegionTotal() { return numeratorRegionTotal; }
		public double getDenominatorRegionTotal() { return denominatorRegionTotal; }


		//TODO MG: I dont really understand what is happening here so I will just compute directly
		//Using the cache may speed this up, we should discuss
		public RatioScore processWindow(Annotation annotation, RatioScore previousScore) {
			return processWindow(annotation);
		}
	}
	
	
	
	/**
	 * @author engreitz
	 * Create a RatioScore from two CountScores
	 */
	public static class RatioScoreIterator implements Iterator<RatioScore> {
		protected Iterator<CountScore> numerator, denominator;
		
		public RatioScoreIterator(Iterator<CountScore> numeratorItr, Iterator<CountScore> denominatorItr) {
			numerator = numeratorItr;
			denominator = denominatorItr;
		}
		
		@Override
		public boolean hasNext() {
			return numerator.hasNext() && denominator.hasNext();
		}
		
		@Override
		public RatioScore next() {
			return new RatioScore(numerator.next(), denominator.next());
		}

		@Override
		public void remove() {
			throw new UnsupportedOperationException("not supported");
		}
	}
	
}
