package umms.esat;

public class Window {
	private String strand;	// strand (+ or -)
	private String chr;		// chromosome (or segment name)
	private int start;		// genomic start coordinate  
	private int end;        // genomic end coordinate  
	private String name;    // transcript name
	private float count;    // total reads starting in this window  
	
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