package umms.core.model.score;

import java.util.Iterator;

import net.sf.samtools.util.CloseableIterator;
import umms.core.alignment.Alignment;
import umms.core.annotation.Annotation;
import umms.core.annotation.AnnotationCollection;
import umms.core.coordinatesystem.CoordinateSpace;
import umms.core.model.AlignmentModel;
import jsc.distributions.Poisson;

/**
 * @author shari
 * This class represents Poisson enrichment score for ChIP/CLIP from Mikkelsen et al. 2010.
 * Can take either sample and control models or just sample.
 */
public class PoissonEnrichmentScore extends CountScore {

	public static double DEFAULT_REGION_TOTAL = -1.0;
	
	private double Pvalue;
	private double regionLength;
	private double globalLength;
	private CoordinateSpace sampleCoordSpace;
	private CoordinateSpace ctrlCoordSpace;
	private double ctrlCount;               
	private double ctrlTotal;       //Total bases in control coord space         
	private double ctrlRegionTotal; 
	
	public PoissonEnrichmentScore(AlignmentModel sample, AlignmentModel ctrl, Annotation a) {
		super(sample, a);
		sampleCoordSpace = sample.getCoordinateSpace();
		ctrlCoordSpace = ctrl.getCoordinateSpace();
		
		if (!sampleCoordSpace.getChromosomeNames().equals(ctrlCoordSpace.getChromosomeNames())) {
			throw new IllegalArgumentException("Sample coordinate space must match control coordinate space");
		}
		
		setGlobalLength(sample.getGlobalLength());
		setCtrlCount(ctrl.getCount(a));
		setCtrlTotal(ctrl.size());
		
		try {
		setPvalue(calculatePVal(getSampleCount(), getCtrlCount(), getSampleLambda(), getCtrlLambda(), a.getSize()));
		} catch(Exception e){
			logger.info("Cound not set P value for annotation " + a.getName());
		}
		getAnnotation().setScore(getPvalue());
	}
	
	public PoissonEnrichmentScore(AlignmentModel sample, Annotation a) {
		this(sample, sample, a);
	}

	/**
	 * @param a				Sample reads
	 * @param b	        	Control reads
	 * @param lambdaSample	Sample reads total / Genome size (ie expected reads per base)
	 * @param lambdaCtrl	Control reads total / Genome size (ie expected reads per base)
	 * @param w				window size
	 * @return				Poisson pvalue
	 */
	public double calculatePVal(double a, double b, double lambdaSample, double lambdaCtrl, double w) {
		double lambdaWsample = lambdaSample*w;
		double lambdaWctrl = lambdaCtrl*w;
		double ctrlE = b/lambdaWctrl;
		double x = Math.max(1, ctrlE)*lambdaWsample;
		Poisson C = new Poisson(x);
		double p = 1 - C.cdf(a);
		return p;
	}
	
	public static class Processor extends WindowProcessor.AbstractProcessor<PoissonEnrichmentScore> {
		protected AlignmentModel sample;
		protected AlignmentModel ctrl;
		protected double sampleRegionTotal = DEFAULT_REGION_TOTAL;
		protected double ctrlRegionTotal = DEFAULT_REGION_TOTAL;
		protected double regionLength = DEFAULT_REGION_TOTAL;
		private boolean fullyContainedReads;
		
		public Processor(AlignmentModel sample, AlignmentModel ctrl) {
			this.sample = sample;
			this.ctrl = ctrl;
		}
		
		public Processor(AlignmentModel sample) {
			this.sample = sample;
			this.ctrl = sample;
		}
		
		public void initRegion(Annotation a){
			if (a != null){
				sampleRegionTotal = sample.getCount(a);
				ctrlRegionTotal = ctrl.getCount(a);
				regionLength = a.length();
			}
		}
		
		public PoissonEnrichmentScore processWindow(Annotation a){
			return new PoissonEnrichmentScore(sample, ctrl, a);
		}
		
	// Todo: implement calculation from previous score
	public PoissonEnrichmentScore processWindow(Annotation a, PoissonEnrichmentScore previousScore){
			return processWindow(a);
		}
	
	}
	
	public double getAverageCoverage(AlignmentModel data, CoordinateSpace coordSpace) { 
		int regionSize = coordSpace.getSize(annotation);
		CloseableIterator<Alignment> readsIter = data.getOverlappingReads(getAnnotation(), false);
		int basesInReads = 0;
		while(readsIter.hasNext()) {
			Alignment read = readsIter.next();
			basesInReads += read.getOverlap(annotation);
		}
		double avgCoverage = (double) basesInReads / (double)regionSize;
		return avgCoverage;
	}
	
	public void setPvalue(double scanPvalue) { this.Pvalue = scanPvalue; }
	public double getPvalue() { return Pvalue; }
	
	public void setGlobalLength(double d) { globalLength = d; }
	public void setRegionLength(double regionLength) { this.regionLength = regionLength; }
	
	public double getGlobalLength() { return globalLength; }
	public double getRegionLength() { return regionLength; }
	
	public double getSampleLambda() { return getSampleTotal() / getGlobalLength(); }
	public double getCtrlLambda() { return getCtrlTotal() / getGlobalLength(); }
	
	public double getLocalSampleLambda() { return getSampleRegionTotal() / getRegionLength(); }
	public double getCtrlSampleLambda() { return getCtrlRegionTotal() / getRegionLength(); }
	
	public void setSampleCount(double d) { setCount(d); }
	public void setCtrlCount(double d) { ctrlCount = d; }
	public void setSampleTotal(double d) { setTotal(d); }
	public void setCtrlTotal(double d) { ctrlTotal = d; }
	public void setSampleRegionTotal(double d) { setRegionTotal(d); }
	public void setCtrlRegionTotal(double d) { ctrlRegionTotal = d; }
	
	public double getSampleCount() { return getCount(); }
	public double getCtrlCount() { return ctrlCount; }
	public double getSampleTotal() { return getTotal(); }
	public double getCtrlTotal() { return ctrlTotal; }
	public double getSampleRegionTotal() { return getRegionTotal(); }
	public double getCtrlRegionTotal() { return ctrlRegionTotal; }
	
}
