package umms.core.model.score;

import java.io.IOException;

import umms.core.annotation.Annotation;
import umms.core.coordinatesystem.GenomicSpace;
import umms.core.model.AlignmentModel;

import org.apache.commons.math3.distribution.BinomialDistribution;

import broad.core.math.Distribution;
import broad.core.math.ScanStatistics;
import broad.pda.seq.segmentation.AlignmentDataModelStats;

public class WindowAllScore extends CountScore {
	
	public final static double RPKM_OFFSET = 0.001; // reads per Kilobase per Million
	
	// the following variables are for the denominator
	private double dCount;               // use a double to allow for weighted read counting
	private double dTotal;                // whole "genome"
	private double dRegionTotal = DEFAULT_REGION_TOTAL;          // whole chromosome / gene / region under consideration / etc.
	private double numeratorLambda;
	private double denominatorLambda;
	private double numeratorScanPVal;
	private double denominatorScanPVal;

	public WindowAllScore(Annotation a) {
		super(a);
	}
	
	public WindowAllScore(AlignmentModel numerator, AlignmentModel denominator, Annotation annotation) {
		super(numerator, annotation);
		setDenominatorCount(denominator.getCount(annotation));
		setDenominatorTotal(denominator.size());
		getAnnotation().setScore(getRatio());

	}
	
	public WindowAllScore(AlignmentModel numerator, AlignmentModel denominator, Annotation annotation, double numeratorRegionTotal, double denominatorRegionTotal,double rapLambda, double backgroundLambda) {
		this(numerator, denominator, annotation);
		setNumeratorRegionTotal(numeratorRegionTotal);
		setDenominatorRegionTotal(denominatorRegionTotal);
		setNumeratorLambda(rapLambda);
		setDenominatorLambda(backgroundLambda);
		setNumeratorScanPVal(ScanStatistics.calculatePVal(new Double(getCount()).intValue(), numerator.getGlobalLambda(), numerator.getCoordinateSpace().getSize(annotation), numerator.getGlobalLength()));
		setDenominatorScanPVal(ScanStatistics.calculatePVal(new Double(getCount()).intValue(), denominator.getGlobalLambda(), denominator.getCoordinateSpace().getSize(annotation), denominator.getGlobalLength()));

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
	public double getDenominatorLambda(){ return denominatorLambda;}
	public double getNumeratorLambda(){ return numeratorLambda;}
	
	private void setDenominatorLambda(double backgroundLambda) {
		this.denominatorLambda = backgroundLambda;
	}

	private void setNumeratorLambda(double rapLambda) {
		this.numeratorLambda = rapLambda;
	}

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
	 * Based on prussell's function in PairedSampleCoverage
	 * Compute Skellam P-value of read counts in a region given the two parameters
	 * denominatorLambda Poisson lambda for background sample
	 *  numeratorLambda Poisson lambda for signal sample
	 *  denominatorCount The count for background sample
	 *  numeratorCount The count for signal sample
	 * @return The probability under the null hypothesis of observing a greater difference
	 * @throws IOException
	 */

	public double getSkellamPValue() {
		double pval;
		if(getNumeratorLambda() <=0.0 || denominatorLambda<=0.0){
			pval = 1.0;
		}
		else if(getNumeratorCount() < numeratorLambda){
			pval = 1.0;
		}
		else{
			pval = Distribution.skellamRightTail(new Double(getNumeratorCount() - getDenominatorCount()).intValue(), numeratorLambda, denominatorLambda);
		}
		return -1.0 * Math.log10(pval);
	}
	
	public void setNumeratorScanPVal(double pval) {
		numeratorScanPVal = pval;
	}
	public void setDenominatorScanPVal(double pval) {
		denominatorScanPVal = pval;
	}
	
	public double getNumeratorScanPVal() {
		return numeratorScanPVal;
	}
	public double getDenominatorScanPVal() {
		return denominatorScanPVal;
	}
	
	/**
	 * Binomial test: combine numerator and denominator counts, then calculate cumulative probability of seeing >= this proportion
	 * @return -log10 p-value for the numerator having more reads than expected.
	 */
	public double getEnrichmentBinomialPValue() {
		int total = (int) Math.round(getNumeratorCount() + getDenominatorCount());
		double p = getNumeratorTotal() / (getNumeratorTotal() + getDenominatorTotal());
		BinomialDistribution b = new BinomialDistribution(total, p);
		double pGreater = 1 - b.cumulativeProbability((int) Math.round(getNumeratorCount()));
		double pEqual = b.probability((int) Math.round(getNumeratorCount()));
		return -1.0 * Math.log10(pGreater + pEqual);
	}
	
	/**
	 * Binomial test: combine numerator and denominator counts, then calculate cumulative probability of seeing less than this proportion
	 * @return -log10 p-value for the numerator having less reads than expected.
	 */
	public double getDepletionBinomialPValue() {
		int total = (int) Math.round(getNumeratorCount() + getDenominatorCount());
		double p = getNumeratorTotal() / (getNumeratorTotal() + getDenominatorTotal());
		BinomialDistribution b = new BinomialDistribution(total, p);
		return -1.0 * Math.log10(b.cumulativeProbability((int) Math.round(getNumeratorCount())));
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
		return annotation.toBED() + 
				"\t" + getRatio() + 			//12
				"\t" + getLog2Ratio() + 		//13
				"\t" + getNumeratorCount() + 	//14
				"\t" + getNumeratorRPKM() + 	//15
				"\t" + getNumeratorRegionTotal() + //16
				"\t" + getNumeratorTotal() + 	//17
				"\t" + getDenominatorCount() + 	//18
				"\t" + getDenominatorRPKM() + 	//19
				"\t" + getDenominatorRegionTotal() + //20
				"\t" + getDenominatorTotal() + //21
				"\t" + getEnrichmentBinomialPValue()+	//22
				"\t" + getDepletionBinomialPValue()+	//23 
				"\t" + getSkellamPValue() +			//24
				"\t" + getNumeratorScanPVal() +		//25
				"\t" + getDenominatorScanPVal();	//26
	}
	
	public static class Processor extends WindowProcessor.AbstractProcessor<WindowAllScore> {
		protected AlignmentModel numerator, denominator;
		protected double numeratorRegionTotal = DEFAULT_REGION_TOTAL;
		protected double denominatorRegionTotal = DEFAULT_REGION_TOTAL;
		protected double numeratorLambda = 0.0;
		protected double denominatorLambda = 0.0;
		
		public Processor(AlignmentModel numerator, AlignmentModel denominator) {
			this.numerator = numerator;
			this.denominator = denominator;
		}
		
		public WindowAllScore processWindow(Annotation annotation) {
			return new WindowAllScore(numerator, denominator, annotation, numeratorRegionTotal, denominatorRegionTotal,numeratorLambda,denominatorLambda);
		}
		
		public void initRegion(Annotation region) {
			if (region != null) {
				numeratorRegionTotal = numerator.getCount(region);
				denominatorRegionTotal = denominator.getCount(region);
				//Calculate lambdas
				String chr = region.getChr();
				denominatorLambda = denominator.getRefSequenceLambda(chr);
				denominatorLambda = (denominatorLambda*denominator.getRefSequenceLength(chr))/
						((GenomicSpace)denominator.getCoordinateSpace()).getUnmaskedLength(chr);
				numeratorLambda = numerator.getRefSequenceLambda(chr);
				numeratorLambda = (numeratorLambda*numerator.getRefSequenceLength(chr))/
						((GenomicSpace)numerator.getCoordinateSpace()).getUnmaskedLength(chr);
				System.out.println("Numerator lambda: "+numeratorLambda+" denominatorLambda "+ denominatorLambda);
			}
		}

		public double getNumeratorRegionTotal() { return numeratorRegionTotal; }
		public double getDenominatorRegionTotal() { return denominatorRegionTotal; }


		//TODO MG: I dont really understand what is happening here so I will just compute directly
		//Using the cache may speed this up, we should discuss
		public WindowAllScore processWindow(Annotation annotation, WindowAllScore previousScore) {
			return processWindow(annotation);
		}
	}
	
}
