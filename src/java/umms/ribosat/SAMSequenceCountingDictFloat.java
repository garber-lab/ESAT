package umms.ribosat;

import java.util.HashMap;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceRecord;

public class SAMSequenceCountingDictFloat extends SAMSequenceCountingDict {
	
	protected HashMap<String, HashMap<String,float[]>> startCounts = new HashMap<String, HashMap<String, float[]>>(); 	

	// **** FLOATING-POINT version simply increments the count by x
    public void incrementStartCounts(String refName, String strand, int alignStart, float fractCount) {
    	startCounts.get(refName).get(strand)[alignStart]+=fractCount;   // increment the counter
    }

    public void copyToLocalCounts(String chr, String strand, int eStart, int cStart, int eLen, float[] floatCounts) {
    	System.arraycopy(startCounts.get(chr).get(strand), eStart, floatCounts, cStart, eLen);
    }
    
    public void updateCount(final SAMRecord r, String multimap, boolean stranded) {
    	/** 
    	 * increments the counter for how many reads had alignments beginning at this position.
    	 * The total counts are stored as short ints used as unsigned short ints. If the count is
    	 * negative, it indicates that the total count is >32767 and <65536, and should be converted
    	 * to the correct int value by adding 65536.
    	 * 
    	 * @param	r	a SAMRecord, a single alignment record
    	 * @see		SAMRecord
    	 */
    	String refName;
    	int alignStart;
    	float fractCount = 1; 
    	String strand;
    	
    	if (stranded & r.getReadNegativeStrandFlag()) {
    		strand = "-";
    	} else  {
    		strand = "+";
    	}
     	
    	// Check to see if storage has already been created for this reference sequence:
    	refName = r.getReferenceName();
    	alignStart = (int)(r.getAlignmentStart())-1;   // alignments are 1-based, arrays are 0-based
    	String cString = r.getCigarString(); 
    	// Note: if the CigarString is "*", it indicates that the read is unmapped. It would be better 
    	//       if SAMRecord had a isMapped() method.
    	if (cString!="*" && !startCounts.containsKey(refName)) {
    		// Find the maximum coordinate of the refName in the dictionary
    		SAMSequenceRecord seq = this.getSequence(refName);
    		// Allocate a float array for storage of the number of reads starting at each location on each strand
    		startCounts.put(refName, new HashMap<String, float[]>());
    		startCounts.get(refName).put("+", new float[seq.getSequenceLength()]);   // forward strand
    		startCounts.get(refName).put("-", new float[seq.getSequenceLength()]);   // forward strand
    		
    	}
    	// Skip unaligned reads:
    	// (multimap should only ever be "scale" for this method)
    	if (cString!="*") {
    		int n = getMultimapCount(r);
    		fractCount = 1f/n;
    		incrementStartCounts(refName, strand, alignStart, fractCount);
    	}
    }

    public boolean startCountsHasKey(String chr) {
    	return startCounts.containsKey(chr);
    }
    
    public float getStartCounts(String chr, String strand, int i) {
		return startCounts.get(chr).get(strand)[i];
    }
}