package umms.esat;

import java.util.LinkedList;
import broad.core.datastructures.IntervalTree;

public class TranscriptCountInfo {
	
	/* This object stores the significant windows spanning all transcripts of a gene, as well
	 * as an interval tree covering ALL exons of all transcripts of the gene and any additional
	 * extension past the gene boundary. (This is required for determining total counts for each 
	 * gene for each experimental condition.)
	 */
	String name;
	String strand;
	LinkedList<Window> sigWindows;
	IntervalTree<String> exonTree;
	
	public TranscriptCountInfo(String gName, String gStrand, LinkedList<Window> wList, IntervalTree<String> iTree) {
		name = gName;
		strand = gStrand;
		sigWindows = wList;
		exonTree = iTree;
	}
	
	public LinkedList<Window> getWindows() {
		return sigWindows;
	}

	public IntervalTree<String> getITree() {
		return exonTree;
	}
	
	public String getName() {
		return name;
	}
	
	public String getStrand() {
		return strand;
	}
}
