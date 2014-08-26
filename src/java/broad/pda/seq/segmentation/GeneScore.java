package broad.pda.seq.segmentation;

import nextgen.core.annotation.Gene;
import nextgen.core.feature.Window;

public class GeneScore extends Gene{

	private double scanPvalue;
	private double enrichment;
	private double numberOfReads;
	private double avergeNumberOfReads;
	private double RPKM;
	private double localLambda;
	private double geneLength;
	private double nominalPValue;
	private double fullyContainedNumberOfReads;
	
	
	//From ContinuousDataAlignmentModel.score()
	/**
	 * A class to represent the standard scores assigned to a window or Gene in the ContinuousDataAlignmentModel
	 */
	public GeneScore(Gene gene, double[] scores){
		super(gene);
		this.setScanPvalue(scores[0]);
		this.setEnrichment(scores[1]);
		this.setNumberOfReads(scores[2]);
		this.setAvergeNumberOfReads(scores[3]);
		this.setRPKM(scores[4]);
		this.setLocalLambda(scores[5]);
		this.setGeneLength(scores[6]);
		this.setNominalPValue(scores[7]);
		this.setFullyContainedNumberOfReads(scores[8]);
	}
	
	public GeneScore(Gene gene){
		super(gene);
	}
		
	public void setEnrichment(double enrichment) {
		this.enrichment = enrichment;
	}

	
	
	public double getEnrichment() {
		return enrichment;
	}

	public void setNumberOfReads(double numberOfReads) {
		this.numberOfReads = numberOfReads;
	}

	public double getNumberOfReads() {
		return numberOfReads;
	}

	public void setAvergeNumberOfReads(double avergeNumberOfReads) {
		this.avergeNumberOfReads = avergeNumberOfReads;
	}

	public double getAvergeNumberOfReads() {
		return avergeNumberOfReads;
	}

	public void setRPKM(double rPKM) {
		RPKM = rPKM;
	}

	public double getRPKM() {
		return RPKM;
	}

	public void setLocalLambda(double localLambda) {
		this.localLambda = localLambda;
	}

	public double getLocalLambda() {
		return localLambda;
	}

	public void setGeneLength(double geneLength) {
		this.geneLength = geneLength;
	}

	public double getGeneLength() {
		return geneLength;
	}

	public void setNominalPValue(double nominalPValue) {
		this.nominalPValue = nominalPValue;
	}

	public double getNominalPValue() {
		return nominalPValue;
	}


	public void setScanPvalue(double scanPvalue) {
		this.scanPvalue = scanPvalue;
	}


	public double getScanPvalue() {
		return scanPvalue;
	}

	public void setFullyContainedNumberOfReads(double fullyContainedNumberOfReads) {
		this.fullyContainedNumberOfReads = fullyContainedNumberOfReads;
	}

	public double getFullyContainedNumberOfReads() {
		return fullyContainedNumberOfReads;
	}

	
	
}
