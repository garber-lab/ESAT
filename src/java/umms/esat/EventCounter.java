package umms.esat;

public class EventCounter {
	private String name;
	private float[] counts;  // float to handle partial counts from multimapped reads
	private float sumCounts;  // test feature to make sure that the total number of reads is the sum of counts[]
	
	public EventCounter(String nodeName, int n) {
		name = nodeName;
		counts = new float[n];
	}
	
	public void setSumCounts(float x) {
		sumCounts = x;
	}
	
	public void incrementCount(int n) {
		counts[n]++;
	}
	
	public void addCount(int n, float val) {
		counts[n]+=val;
	}
	
	public String getName() {
		return name;
	}
	
	public float[] getAllCounts() {
		return counts;
	}
	
	public float getCounts(int n) {
		return counts[n];
	}
}