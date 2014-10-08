package umms.esat;

import umms.esat.Window;

public class EventCounter {
	private String name;
	private float[] counts;  // float to handle partial counts from multimapped reads
	private float sumCounts;  // test feature to make sure that the total number of reads is the sum of counts[]
	private Window w;

	public EventCounter(String nodeName, int n, Window inW) {
		name = nodeName;
		counts = new float[n];
		w = inW;
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
	
	public boolean hasIntervalTree() {
		return w.hasITree();
	}
	
	public void incrementIntervalCount(int iStart, int iEnd, int n) {
		// If the node's window contains an interval tree, need to do a second-level check to see if
		// the read actually falls into one of the exon ranges. If there is no interval tree, it means that 
		// this genomic range does not span more than one exon, so just update the count without any further checking.

		if (!w.hasITree() || w.hasOverlap(iStart, iEnd)) { 
			counts[n]++;
		} 
	}

	public void addIntervalCount(int iStart, int iEnd, int n, float fractCount) {
		// If the node's window contains an interval tree, need to do a second-level check to see if
		// the read actually falls into one of the exon ranges. If there is no interval tree, it means that 
		// this genomic range does not span more than one exon, so just update the count without any further checking.

		if (!w.hasITree() || w.hasOverlap(iStart, iEnd)) { 
			counts[n]+=fractCount;
		} 
	}
}