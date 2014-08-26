package broad.pda.seq.segmentation;

public class GeneCounts {
	private double exonicCounts;
	private double doubleronicCounts;
	private double utr3Counts;
	private double utr5Counts;
	
	public GeneCounts(double exonicCounts, double doubleronicCounts, double utr3Counts, double utr5Counts) {
		this.exonicCounts = exonicCounts;
		this.doubleronicCounts = doubleronicCounts;
		this.utr3Counts = utr3Counts;
		this.utr5Counts = utr5Counts;
	}
	
	
	public double getExonicCounts() {
		return exonicCounts;
	}
	public double getdoubleronicCounts() {
		return doubleronicCounts;
	}
	public double getUtr3Counts() {
		return utr3Counts;
	}
	public double getUtr5Counts() {
		return utr5Counts;
	}
	

}
