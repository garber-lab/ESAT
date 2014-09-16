package umms.esat;

import umms.esat.Window;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.LinkedList;
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

import broad.core.datastructures.IntervalTree;
import umms.core.annotation.Annotation;
import umms.core.annotation.Annotation.Strand;
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
    	 *  	X SAMFileReader bamReader = new SAMFileReader(bamFile, true);
    	 *    	> SAMFileReader bamReader = new SAMFileReader(bamFile);  non-eager reader is faster
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
    	alignStart = (int)(r.getAlignmentStart())-1;   // alignments are 1-based, arrays are 0-based
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
    			eStart = exon.getStart()-1;          // exon start (1-based annotation, 0-based arrays)
    			eEnd = exon.getEnd()-1;              // exon end (1-based annotation, 0-based arrays)
    			/* compute the sum read starts over the exon interval */
    			for (int i=eStart; i<eEnd; i++) {
    				// NOTE: the value of getEnd() is both 1-based and "right-open", so it points to the first 
    				// base AFTER the end of the exon (which makes the arithmetic work better, apparently). For this reason, 
    				// the loop uses "i<eEnd", rather than "i<=eEnd".
    				countSum+=this.startCounts.get(chr)[i];
    			}
    		} else {
    			logger.warn("Counts array does not contain key: "+chr+"\n");
    			break;    // assume that all exons are from the same chromosome (which wasn't found)
    		}
    	}
    	return countSum;
    }
    
    public LinkedList<Window> countWindowedReadStarts(final Set<Annotation> eSet, 
    									final IntervalTree<String> iTree_fwd,
    									final IntervalTree<String> iTree_rev,
    									final int window, 
    									final int overlap, 
    									final int extend,
    									final Gene gene) {
    	/**
    	 * sums the count of all reads starting within sliding windows across all of the 
    	 * exons in this Set. An array is created by concatenating the exons and extending by
    	 * "extend" bases, then a window is stepped across the array overlapping "overlap" bases 
    	 * with the previous window, counting the number of reads in each window.  
    	 * 
    	 * @param	eSet	a Set of Annotations defining the boundaries of a set of exons
    	 * @return			a list of windows and counts. 
    	 */
    	int aLen = 0;
    	short[] counts; 
    	float[] floatCounts; 
    	int[] gCoords;    // genomic coordinates of each location
    	int eStart;
    	int eEnd;
    	String gStrand;   // strand of the gene
    	int localExtend = extend;    // default extension
    	
    	LinkedList<Window> wList = new LinkedList<Window>();    // Is this the right type?...
    	
    	Iterator<Annotation> eIter = eSet.iterator();    // iterator over the exons
    	
    	/* get the BED file segment to which this segment is aligned: */
    	String chr = eSet.iterator().next().getChr();
    	
    	// reset the iterator:
    	eIter = eSet.iterator();

    	/* Adjust the extension for this gene:
    	 * If negative strand and extension will end before base 0, or
    	 * If positive strand and extension will end after the last base, 
    	 * shorten the extension.
    	 */
    	
    	if (gene.isNegativeStrand()) {
    		int gStart = gene.getStart();
    		if ((gStart-extend)<0) {
    			localExtend = gStart;     // limit extension to the beginning of the chromosome/segment
    			logger.warn(extend+"-base extension for gene "+gene.getName()+" reduced to "+localExtend+" bases");
    		}
    	} else {
    		int gEnd = gene.getEnd();
    		int segEnd = this.startCounts.get(chr).length;
    		if (((gEnd+extend)-1)>segEnd) {
    			localExtend = segEnd-gEnd;   // limit extension to the end of the chromosome/segment
    			logger.warn(extend+"-base extension for gene "+gene.getName()+" reduced to "+localExtend+" bases");
    		}
    	}
    	
    	/* build the array containing all counts for this exon set: */
    	/* first, scan through the exon set to determine how long the composite array needs to be */
    	while (eIter.hasNext()) {
    		Annotation exon = eIter.next();
    		// Arithmetic works better with 0-based, exclusive (native) coordinates
    		eStart = exon.getStart();  
    		eEnd = exon.getEnd();      
    		aLen += (exon.getEnd()-exon.getStart());  // sum the exon lengths
    	}
    	/* add the extension */
    	aLen += localExtend;
    	
    	/* allocate the storage */
    	counts = new short[aLen];
    	floatCounts = new float[aLen];
    	gCoords = new int[aLen];
    	
    	/* reset the iterator */
    	eIter = eSet.iterator();
    	/* copy the counts from startCounts */
    	
    	int cStart;
    	if (gene.isNegativeStrand()) {
    		cStart = localExtend;      // leave room to copy counts from lower genomic coordinates for negative-strand genes
    	} else {
    		cStart = 0;	      // otherwise, start at the beginning of the startCounts array
    	}
    	eStart = 0;    // save the start position of the last exon
    	eEnd = 0;      // save the end position of the last exon
    	int eLen = eEnd-eStart; 
    	Strand strand = gene.getStrand();
    	while (eIter.hasNext()) {
    		Annotation exon = eIter.next();    // get the exon
    		chr = exon.getChr();        // exon chromosome
    		if (!this.startCounts.containsKey(chr)) {
    			logger.warn("Counts array does not contain key: "+chr);
    			return new LinkedList<Window>();
    		}
    		// use the native coordinate system (0-based, exclusive):
    		eStart = exon.getStart();          // exon start
			eEnd = exon.getEnd();              // exon end
			eLen = eEnd-eStart;        // exon length
    		System.arraycopy(this.startCounts.get(chr),eStart, counts, cStart, eLen);
    		// fill the gCoords array with genomic coordinates:
    		for (int i=0; i<eLen; i++) {
    			gCoords[cStart+i]=eStart+i;
    		}
    		cStart+=eLen;    // point to the start position of the next copied segment
    	}
    	/* extend past the end of the last exon by "localExtend" bases. For now, just report a warning 
    	 * if this collides with another exon */
    	if (gene.isNegativeStrand()) {
    		//eEnd = gene.getFirstExon().getStart()+1;
    		eEnd = gene.getFirstExon().getStart();
    		eStart = eEnd-localExtend;
    		cStart = 0;               // start copying into the beginning of the array
    	} else {
    		eStart = eEnd;            // Since the coordinate system is 0-based, exclusive, eEnd points to the first base AFTER the end of the last exon
    		eEnd = eStart+localExtend;     // exclusive index of the end of the extended segment
    	}
    	eLen = localExtend;            // number of bases to extend past the end of the last exon
    	
    	/* Copy the extension counts */
    	System.arraycopy(this.startCounts.get(chr),eStart, counts, cStart, eLen);
		// fill the gCoords array with genomic coordinates:
		for (int i=0; i<eLen; i++) {
			gCoords[cStart+i]=eStart+i;
		}

    	// check if extension overlaps other genes:
    	if (strand.toString().equals("+")) {
    		if (iTree_fwd.numOverlappers(eStart, eEnd)>0) {
    			logger.warn(localExtend+" base extension of "+gene.getName()+" has "+iTree_fwd.numOverlappers(eStart, eEnd)+" overlaps (fwd)");
    			Iterator<String> oLapIter = iTree_fwd.overlappingValueIterator(eStart, eEnd);
    			while (oLapIter.hasNext()) {
    				logger.warn("   "+oLapIter.next());
    			}
    		}
    	} else {
    		if (iTree_rev.numOverlappers(eStart, eEnd)>0) {
    			logger.warn(localExtend+" base extension of "+gene.getName()+" has "+iTree_rev.numOverlappers(eStart, eEnd)+" overlaps (rev)");
    			Iterator<String> oLapIter = iTree_rev.overlappingValueIterator(eStart, eEnd);
    			while (oLapIter.hasNext()) {
    				logger.warn("   "+oLapIter.next());
    			}
    		}
    	}
    	
    	/* Copy the short (unsigned) int counts array to floating-point */
    	for (int i=0;i<aLen; i++) {
    		/* need to figure out how to do this same thing for floats... */
    		if (counts[i]<0) {
    			floatCounts[i] = (float) (counts[i]+65536);    // adjust for counts >32767 and <65536
    		} else {
    			floatCounts[i] = (float) counts[i];
    		}
    	}

    	/* Iterate a sliding window across the counts array and sum the counts over the window. Create a new Window object
    	 * for each step and add them to the output list.
    	 */
    	int sumStart = 0;					// start of summing window
    	int sumEnd = sumStart+window;       // end of summing window (exclusive)
    	if (gene.isNegativeStrand()) {
    		gStrand = "-";
    	} else {
    		gStrand = "+";
    	}
    	
    	while (sumEnd < aLen) {
    		Window thisWindow = new Window(gStrand, chr, gCoords[sumStart], gCoords[sumEnd], gene.getName());
    		float countSum = 0;
    		for (int i=sumStart; i<sumEnd; i++) {
    			countSum += floatCounts[i];
    		}
    		thisWindow.setCount(countSum);    // update count for this window
    		wList.add(thisWindow);            // add this window to the output list
    		/* update start and end for the next window */
    		sumStart += (window-overlap);
    		sumEnd = sumStart+window;
    	}
    	
    	return wList;
    }
        
    
    public void simpleCountTranscriptReadStarts (final Map<String,Collection<Gene>> annotations) {
    	/**
    	 * counts the number of reads starting within the boundaries of the exons of a transcript/gene
    	 * where the exon boundaries are specified in the Set of Annotations returned by the getExonSet()
    	 * method on the Gene. The sum of reads count is stored in the Gene "score" field, via the
    	 * setScore() method on the Gene. (NOTE: it might be better to modify the Gene object to contain
    	 * a new field for this count, in case the "score" field is used for something else.
    	 * 
    	 * @param	annotations	a Map of all available Gene Annotations
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
    
    public HashMap<String,HashMap<String, LinkedList<Window>>> countWindowedTranscriptReadStarts (final Map<String,Collection<Gene>> annotations, 
    												final int window, 
    												final int overlap,
    												final int extend) {
    	/**
    	 * counts the number of reads starting within each sliding window of "window" bases with and
    	 * overlap of "overlap" bases. The transcript is extended past the last exon "extend" bases past
    	 * the end of the last exon, unless it collides with the next gene in the Annotation Set.
    	 * 
    	 * @param	annotations	a Map of all available Gene Annotations
    	 * @param	window	the length of the sliding window, in bases
    	 * @param	overlap	as the window slides, the overlap of each window with the previous one
    	 * @param	extend	number of bases past the last exon to extend the counting
    	 */
    	HashMap<String,HashMap<String,LinkedList<Window>>> countsMap = new HashMap();   // countsMap[chr][transID][Window]
    	LinkedList<Window> eCount;  
    	
    	// Iterate over all "chromosomes":
    	for(String chr:annotations.keySet()){
			// Iterate over all transcripts:
    		Collection<Gene> aSet = annotations.get(chr);
    		
    		/* Build strand-specific IntervalTrees containing all exons for this segment */
    		IntervalTree<String> iTree_fwd = new IntervalTree();   // forward strand tree
    		IntervalTree<String> iTree_rev = new IntervalTree();   // reverse strand tree
    		for(Gene gene : annotations.get(chr)) {
				Set eSet = gene.getExonSet();
				String gName = gene.getName();   // gene/transcript name
				Iterator<Annotation> eIter = eSet.iterator();
				int eID = 0;    // exon ID
				while (eIter.hasNext()) {
					Annotation exon = eIter.next();
					int eStart = exon.getStart();
					int eEnd = exon.getEnd();
					Strand strand = exon.getStrand();
					String intervalName = gName+"_"+eID;
					if (strand.toString().equals("+")) {
						// forward strand interval tree
						iTree_fwd.put(eStart, eEnd, intervalName);
					} else {
						iTree_rev.put(eStart, eEnd, intervalName);
					}						
					eID++;	// increment the exon ID
				}
    		}

    		for(Gene gene : annotations.get(chr)) {
    			String gName = gene.getName();
				Set eSet = gene.getExonSet();
				// sum all counts starting within this exon set:
				eCount = this.countWindowedReadStarts(eSet, iTree_fwd, iTree_rev, window, overlap, extend, gene);
				
				// add this window set to the HashMap:
				if (countsMap.containsKey(chr)) {
					if (countsMap.get(chr).containsKey(gName)) {
						logger.warn("Duplicate entry for transcript "+gName+". Skipping...");
					} else {
						countsMap.get(chr).put(gName, eCount);   // add the counts windows to the hashmap
					}
				} else {
					countsMap.put(chr, new HashMap());           // add an entry in the main map for the chromosome/segment
					countsMap.get(chr).put(gName, eCount);       // add the counts windows to the hashmap
				}
			}	
		}
    	
    return countsMap;
    	
    }    
}
