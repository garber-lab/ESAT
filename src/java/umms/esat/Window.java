package umms.esat;

public class Window {
	private String strand;	// strand (+ or -)
	private String chr;		// chromosome (or segment name)
	private int start;		// genomic start coordinate  
	private int end;        // genomic end coordinate  
	private String name;    // transcript name
	private float count;    // total reads starting in this window
	private float[] counts;  // allows for collection of counts from multiple conditions
	
	// empty constructor:
	public Window() {
	}

	// main constructor:
	public Window(String wStrand, String wChr, int wStart, int wEnd) {
		strand = wStrand;
		chr = wChr;
		start = wStart;
		end = wEnd;
	}
	
	// constructor with name:
	public Window(String wStrand, String wChr, int wStart, int wEnd, String tName) {
		strand = wStrand;
		chr = wChr;
		start = wStart;
		end = wEnd;
		name = tName;
	}
	
	// constructor with name and multiple experiment/condition counters:
	public Window(String wStrand, String wChr, int wStart, int wEnd, String tName, int nExp) {
		strand = wStrand;
		chr = wChr;
		start = wStart;
		end = wEnd;
		name = tName;
		counts = new float[nExp];
	}
	
	// METHODS //

	// getters and setters:
	public void setStrand(String wStr) {
		if (wStr.equals("+") || wStr.equals("-")) {
			strand = wStr;
		}  // should throw an exception here
	}
	
	public void setChr(String wChr) {
		chr = wChr;
	}

	public void setStart(int wStart) {
		start = wStart;
	}

	public void setEnd(int wEnd) {
		end = wEnd;
	}

	public void setCount(float wCount) {
		count = wCount;
	}

	public void incrementCounts(int eIdx) {
		// update the counts for experiment eIdx:
		counts[eIdx]++;
	}

	public String getStrand() {
		return strand;
	}
	
	public String getChr() {
		return chr;
	}

	public int getStart() {
		return start;
	}

	public int getEnd() {
		return end;
	}

	public float getCount() {
		return count;
	}
	
	// get counts for all experiments:
	public float[] getCounts() {
		return counts;
	}
	
	// get counts for experiment eIdx:
	public float getCounts(int eIdx) {
		return counts[eIdx];
	}
	
	public String getName() {
		return name;
	}

	// OTHER METHODS //
	public boolean isPositive() {
		if (strand.equals("+")) {
			return true;
		} else {
			return false;
		}
	}
	
	public float addCounts(float newCounts) {
		count += newCounts;
		return count;
	}
}