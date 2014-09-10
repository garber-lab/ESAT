package umms.esat;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.Iterator;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.SAMSequenceDictionary;




import net.sf.samtools.SAMFileReader.ValidationStringency;

//import org.apache.commons.math3.util.MultidimensionalCounter.Iterator;
import org.apache.log4j.Logger;

import umms.core.annotation.Annotation;
import umms.core.annotation.Gene;

public class SAMSequenceCountingDict extends SAMSequenceDictionary {
/**
 *     Extends the SAMSequenceDictionary class to add startCounts, a set of simple counter arrays that 
 *     keep track of the number of reads beginning at each genomic location. The startCounts arrays are
 *     short int[] to save space, since there is one entry for every location in the genome. The arrays 
 *     are not constructed until the "reference" (generally a chromosome ID) is not found as a key in the
 *     startCounts HashMap. This is done to avoid allocating unnecessary storage. Although the startCounts
 *     arrays are short ints, they are treated as unsigned short ints, so counts >32767 and <65536 appear
 *     as negative numbers. Negative short ints can be converted to the correct int value by adding 65536.
 *     The maximum allowable count is -1, which converts to 65535. If the count exceeds 65535 for any 
 *     location, a warning message is generated and the count is left at 65535.
 *     
 *     @param	none, but requires sequences from an existing SAMSequenceDictionary via the method copySequences()
 *     @see		SAMSequenceDictionary
 */
	/* start counts holds the count of the number of reads that start at each base within each of the 
	 * dictionary header entries */
	private HashMap<String, short[]> startCounts = new HashMap();
//	private long allocatedMem = 0;
//	private int allocatedArrays = 0;
	public Logger logger;
    
    public SAMSequenceCountingDict () {
    	super();
    }
    
    public void setLogger(Logger logger) {
    	this.logger = logger;	
    }
    public void copySequences(final SAMSequenceDictionary dict) {
    	/**
    	 * Copies the sequence information from an existing SAMSequenceDictionary, which contains all
    	 * possible alignment references for this SAM/BAM file and their maximum lengths.
    	 *  
    	 * Suggested code fragment:
    	 *  // open the file reader with eager decoding (true) so that each line is parsed on reading:	
    	 *  	SAMFileReader bamReader = new SAMFileReader(bamFile, true);
		 *	// set the ValidationStringency flag to LENIENT to generate warnings, rather than errors, 
		 *  // if MAPQ is not 0 for unmapped reads:
		 *  	bamReader.setValidationStringency(ValidationStringency.LENIENT);	
		 *	// create the new counting dictionary:
		 * 		SAMSequenceCountingDict bamDict = new SAMSequenceCountingDict();
    	 * 	// access the header from the open reader:	
    	 * 		SAMFileHeader bamHeader = bamReader.getFileHeader();
    	 *  // copy the "sequence" information from the header (gets the list of available possible
    	 *  // alignment references and their lengths):
		 *		bamDict.copySequences(bamHeader.getSequenceDictionary());
		 *  // give the dictionary access to the logger, to allow info, warnings, and errors to be reported: 
		 *		bamDict.setLogger(logger);
		 *
		 *	@see	SAMSequenceDictionary
    	 */
    	this.setSequences(dict.getSequences());
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
    	// Check to see if storage has already been created for this reference sequence:
    	refName = r.getReferenceName();
    	alignStart = (int)(r.getAlignmentStart());
    	String cString = r.getCigarString(); 
    	// Note: if the CigarString is "*", it indicates that the read is unmapped. It would be better 
    	//       if SAMRecord had a isMapped() method.
    	if (cString!="*" && !this.startCounts.containsKey(refName)) {
    		// Find the maximum coordinate of the refName in the dictionary
    		SAMSequenceRecord seq = this.getSequence(refName);
    		// Allocate a short int array for storage of the number of reads starting at each location
    		this.startCounts.put(refName, new short[seq.getSequenceLength()]);
    	}
    	// Skip unaligned reads:
    	if (cString!="*") {
    		// NOTE: This effectively treats the counts as short, unsigned ints. Negative numbers indicate
    		//       that the total count is >32767 and less than 65535. The test for counts!=1 limits the
    		//       total counts to 65535 after converting back to an int.
    		if (this.startCounts.get(refName)[alignStart]!=-1) {
    			this.startCounts.get(refName)[alignStart]++;   // increment the counter
    		} else {
    			logger.warn("location "+alignStart+" in "+refName+" has >65535 counts.");
    			// NOTE: since the counts arrays are short ints, counts between 32767 and 65535 will be negative.
    			//       To adjust the negative values, <correct int value> = 65536+<negative count value>)
    		}
    	}
    }
    
    public int countExonReadStarts(final Set<Annotation> eSet) {
    	/**
    	 * sums the count of all reads starting within all of the exons in this Set.
    	 * 
    	 * @param	eSet	a Set of Annotations defining the boundaries of a set of exons
    	 * @return			the (int) sum of reads
    	 */
    	int countSum = 0;
    	int eStart = 0;
    	int eEnd = 0;
    	
    	Iterator<Annotation> eIter = eSet.iterator();    // iterator over the exons
    	while (eIter.hasNext()) {
    		Annotation exon = eIter.next();    // get the exon
    		String chr = exon.getChr();        // exon chromosome
    		if (this.startCounts.containsKey(chr)) {
    			eStart = exon.getStart();          // exon start
    			eEnd = exon.getEnd();              // exon end
    			/* compute the sum read starts over the exon interval (convert from 0-based to 1-based limits) */
    			for (int i=eStart+1; i<eEnd; i++) {
    				countSum+=this.startCounts.get(chr)[i];
    			}
    		} else {
    			logger.warn("Counts array does not contain key: "+chr+"\n");
    			break;    // assume that all exons are from the same chromosome (which wasn't found)
    		}
    	}
    	return countSum;
    }
    
    public void simpleCountTranscriptReadStarts (final Map<String,Collection<Gene>> annotations) {
    	/**
    	 * counts the number of reads starting within the boundaries of the exons of a transcript/gene
    	 * where the exon boundaries are specified in the Set of Annotations returned by the getExonSet()
    	 * method on the Gene. The sum of reads count is stored in the Gene "score" field, via the
    	 * setScore() method on the Gene. (NOTE: it might be better to modify the Gene object to contain
    	 * a new field for this count, in case the "score" field is used for something else.
    	 * 
    	 * @param	annotation	a Map of all available Gene Annotations
    	 */
    	int eCount = 0;  
    	
    	// Iterate over all "chromosomes":
    	for(String chr:annotations.keySet()){
			// Iterate over all transcripts:
    		for(Gene gene : annotations.get(chr)) {
				Set eSet = gene.getExonSet();
				// sum all counts starting within this exon set:
				eCount = this.countExonReadStarts(eSet);
				// set the gene "score" to the number of reads starting within these exons:
				gene.setScore(eCount);
			}	
		}
    	
    }
}
