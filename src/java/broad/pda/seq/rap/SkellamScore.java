package broad.pda.seq.rap;

import java.io.IOException;

import broad.core.math.Distribution;
import nextgen.core.annotation.Annotation;
import nextgen.core.model.AlignmentModel;
import nextgen.core.model.score.ScanStatisticScore;
import nextgen.core.model.score.WindowProcessor;
import nextgen.core.model.score.WindowScore;

public class SkellamScore extends WindowScore.AbstractWindowScore {

	//numerator is rap and denominator is control
	private int numeratorCount;
	private int denominatorCount;
	private double numeratorLambda;
	private double denominatorLambda;
	private AlignmentModel numerator;
	private AlignmentModel denominator;
	
	public SkellamScore(Annotation a) {
		super(a);
	}
	
	public SkellamScore(AlignmentModel numerator, AlignmentModel denominator, Annotation annotation, double rapLambda, double backgroundLambda) {
		super(annotation);
		this.numerator = numerator;
		this.denominator = denominator;
		setDenominatorCount((int) Math.round(denominator.getCount(annotation)));
		setNumeratorCount((int) Math.round(numerator.getCount(annotation)));
		setNumeratorLambda(rapLambda);
		setDenominatorLambda(backgroundLambda);
		getAnnotation().setScore(getSkellamPValue());
	}
	
	private void setNumeratorCount(int round) {
		this.numeratorCount = round;
	}

	private void setDenominatorCount(int round) {
		this.denominatorCount = round;
	}

	private void setDenominatorLambda(double backgroundLambda) {
		this.denominatorLambda = backgroundLambda;
	}

	private void setNumeratorLambda(double rapLambda) {
		this.numeratorLambda = rapLambda;
	}
	
	private int getNumeratorCount(){
		return numeratorCount;
	}

	private int getDenominatorCount(){
		return denominatorCount;
	}

	@Override
	public double getScore() {
		
		return getSkellamPValue();
	}
	
	public double getNumeratorScanPVal() {
		return new ScanStatisticScore(numerator, getAnnotation()).getScanPvalue();
	}
	public double getDenominatorScanPVal() {
		return new ScanStatisticScore(denominator, getAnnotation()).getScanPvalue();
	}
	
	public String toString() {
		annotation.setScore(getScore());
		return annotation.toBED() + "\t" + getScore() + "\t" + getNumeratorCount() + "\t" + getDenominatorCount() + "\t" + getNumeratorScanPVal()+ "\t" + this.getDenominatorScanPVal() +  "\n";
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
		if(numeratorLambda <=0.0 || denominatorLambda<=0.0){
			return 1.0;
		}
		else if(numeratorCount < numeratorLambda){
			return 1.0;
		}
		else{
			return Distribution.skellamRightTail(numeratorCount - denominatorCount, numeratorLambda, denominatorLambda);
		}
	}

	/**
	 * Returns true if the signal count is less than the chromosome lambda
	 * @return
	 */
	public boolean countLessThanLambda(){
		return(numeratorCount < numeratorLambda);
	}
	
	
	public static class Processor extends WindowProcessor.AbstractProcessor<SkellamScore> {
		protected AlignmentModel numerator, denominator;
		protected double numeratorLambda;
		protected double denominatorLambda;
		
		public Processor(AlignmentModel numerator, AlignmentModel denominator,double numeratorLambda,double denominatorLambda) {
			this.numerator = numerator;
			this.denominator = denominator;
			this.numeratorLambda = numeratorLambda;
			this.denominatorLambda = denominatorLambda;
		}
		
		
		@Override
		public SkellamScore processWindow(Annotation annotation) {
			return new SkellamScore(numerator, denominator, annotation, numeratorLambda,denominatorLambda);
		}

		@Override
		//TODO: Check
		public SkellamScore processWindow(Annotation annotation,
				SkellamScore previousScore) {
			return processWindow(annotation);
		}
		
		public void initRegion(Annotation region) {
			if (region != null) {
				//TODO:CHECK
			}
		}
	}

}
