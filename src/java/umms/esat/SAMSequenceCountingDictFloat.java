package umms.esat;

import java.util.HashMap;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceRecord;

public class SAMSequenceCountingDictFloat extends SAMSequenceCountingDict {
	
	protected HashMap<String, float[]> startCounts = new HashMap<String, float[]>(); 	

	// **** FLOATING-POINT version simply increments the count by x
    public void incrementStartCounts(String refName, int alignStart, float fractCount) {
    	startCounts.get(refName)[alignStart]+=fractCount;   // increment the counter
    }

    public void copyToLocalCounts(String chr, int eStart, int cStart, int eLen, float[] floatCounts) {
    	System.arraycopy(startCounts.get(chr),  eStart, floatCounts, cStart, eLen);
    }
    public void updateCount(final SAMRecord r) {
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
    	
    	// Check to see if storage has already been created for this reference sequence:
    	refName = r.getReferenceName();
    	alignStart = (int)(r.getAlignmentStart())-1;   // alignments are 1-based, arrays are 0-based
    	String cString = r.getCigarString(); 
    	// Note: if the CigarString is "*", it indicates that the read is unmapped. It would be better 
    	//       if SAMRecord had a isMapped() method.
    	if (cString!="*" && !startCounts.containsKey(refName)) {
    		// Find the maximum coordinate of the refName in the dictionary
    		SAMSequenceRecord seq = this.getSequence(refName);
    		// Allocate a short int array for storage of the number of reads starting at each location
    		startCounts.put(refName, new float[seq.getSequenceLength()]);
    	}
    	// Skip unaligned reads:
    	if (cString!="*") {
    		incrementStartCounts(refName, alignStart, fractCount);
    	}
    }

    public boolean startCountsHasKey(String chr) {
    	return startCounts.containsKey(chr);
    }
    
    public float getStartCounts(String chr, int i) {
		return startCounts.get(chr)[i];
    }
    
    public int getChrLength(String chr) {
		return startCounts.get(chr).length;
    }
    
}