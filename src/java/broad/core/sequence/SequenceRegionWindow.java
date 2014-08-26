package broad.core.sequence;

public abstract class SequenceRegionWindow  {
	private SequenceRegion region;
	private double score;
	
	public SequenceRegionWindow(SequenceRegion region) {
		this.region = region;
	}
	
	
	public abstract void computeScore();
	
	public double getScore() {return score; }
	protected void  setScore(double score) { this.score = score; }
	
	public SequenceRegion getRegion() {return region;}
}
