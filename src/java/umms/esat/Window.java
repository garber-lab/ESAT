package umms.esat;

import java.util.Iterator;

import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;

public class Window {
	private String strand;	// strand (+ or -)
	private String chr;		// chromosome (or segment name)
	private int start;		// genomic start coordinate  
	private int end;        // genomic end coordinate  
	private String name;    // transcript name
	private float count;    // total reads starting in this window
	private float[] counts;  // allows for collection of counts from multiple conditions
	private IntervalTree<String> iTree;    // Interval tree for windows that span multiple exons
	
	// empty constructor:
	public Window() {
	}

	// main constructor:
	public Window(String wStrand, String wChr, int wStart, int wEnd) {
		strand = wStrand;
		chr = wChr;
		start = wStart;
		end = wEnd;
		iTree = null; 
	}
	
	// constructor with name:
	public Window(String wStrand, String wChr, int wStart, int wEnd, String tName) {
		strand = wStrand;
		chr = wChr;
		start = wStart;
		end = wEnd;
		name = tName;
		iTree = null;
	}
	
	// constructor with name and multiple experiment/condition counters:
	public Window(String wStrand, String wChr, int wStart, int wEnd, String tName, int nExp) {
		strand = wStrand;
		chr = wChr;
		start = wStart;
		end = wEnd;
		name = tName;
		counts = new float[nExp];
		iTree = null;
	}
	
	// constructor for pre-existing IntervalTrees:
	public Window(String wStrand, String wChr, IntervalTree<String> eTree, String tName, int nExp) {
		strand = wStrand;
		chr = wChr;
		name = tName;
		counts = new float[nExp];
		iTree = eTree;
		start = eTree.min().getStart();   	// start of earliest interval
		end = eTree.max().getEnd();			// end of latest interval
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
	
	// Interval tree methods:
	// Check to see if this window spans multiple exons (i.e., has an interval tree):
	public boolean hasITree() {
		if (iTree==null) {
			return false;
		} else {
			return true;
		}
	}
	
	// Create an interval tree:
	public void addITree() {
		iTree = new IntervalTree<String>();
	}
	
	// add intervals to he interval tree for this window, if necessary:
	public void addIntervals(int iStart, int iEnd, IntervalTree<String>refTree) {
		int nOver = refTree.numOverlappers(iStart, iEnd);
		// if the window range spans only one exon, no need to create a tree for this window
		if (refTree.numOverlappers(iStart, iEnd)==1) {
			return;     // don't create a tree
		} else if (refTree.numOverlappers(iStart, iEnd)==2) {
			// if the two exons are one of the edge exons and the extension, the genomic range
			// should still be contiguous, so don't create a tree. If the name of either of the 
			// overlappers ends with .ext, this should be the case:
			Iterator<String> iIter = refTree.overlappingValueIterator(iStart, iEnd);
			while (iIter.hasNext()) {
				if (iIter.next().endsWith("ext")) {
					return;
				}
			}
		}
		// Otherwise, there are multiple exons and the window spans them:
		iTree = new IntervalTree<String>();   // create the tree
		Iterator<Node<String>> iIter = refTree.overlappers(iStart, iEnd);
		while (iIter.hasNext()) {
			Node<String> n = iIter.next();
			int eStart = n.getStart();  // exon start
			int eEnd = n.getEnd();      // exon end
			if (iStart>eStart) {
				// window begins after the beginning of the exon
				eStart = iStart;
			} 
			if (iEnd<eEnd) {
				// window ends before the end of the exon
				eEnd = iEnd;
			}
			iTree.put(eStart, eEnd, n.getValue());
		}
		return;		
	}
	
	// return true if any interval in the tree overlaps the input interval:
	public boolean hasOverlap(int iStart, int iEnd) {
		if (iTree.numOverlappers(iStart, iEnd)>0) {
			return true;
		} else {
			return false;
		}
	}
}



