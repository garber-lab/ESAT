package umms.esat;

import umms.esat.Window;

import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Map;
import java.util.Set;
import java.util.Iterator;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceDictionary;

//import org.apache.commons.math3.util.MultidimensionalCounter.Iterator;
import org.apache.log4j.Logger;

import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import umms.core.annotation.Annotation;
import umms.core.annotation.Annotation.Strand;
import umms.core.annotation.Gene;

abstract public class SAMSequenceCountingDict extends SAMSequenceDictionary {
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
	//protected HashMap<String, short[]> startCounts = new HashMap<String, short[]>();      // **** short or float 
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
    	setSequences(dict.getSequences());
    }
    
 
   public float countExonReadStarts(final Set<Annotation> eSet) {
    	/**
    	 * sums the count of all reads starting within all of the exons in this Set.
    	 * 
    	 * @param	eSet	a Set of Annotations defining the boundaries of a set of exons
    	 * @return			the (int) sum of reads
    	 */
    	float countSum = 0;
    	int eStart = 0;
    	int eEnd = 0;
    	
    	Iterator<Annotation> eIter = eSet.iterator();    // iterator over the exons
    	while (eIter.hasNext()) {
    		Annotation exon = eIter.next();    // get the exon
    		String chr = exon.getChr();        // exon chromosome
    		String strand = exon.getStrand().toString();
    		if (startCountsHasKey(chr)) {
    			eStart = exon.getStart()-1;          // exon start (1-based annotation, 0-based arrays)
    			eEnd = exon.getEnd()-1;              // exon end (1-based annotation, 0-based arrays)
    			/* compute the sum read starts over the exon interval */
    			for (int i=eStart; i<eEnd; i++) {
    				// NOTE: the value of getEnd() is both 1-based and "right-open", so it points to the first 
    				// base AFTER the end of the exon (which makes the arithmetic work better, apparently). For this reason, 
    				// the loop uses "i<eEnd", rather than "i<=eEnd".
    				countSum+=getStartCounts(chr, strand, i);
    			}
    		} else {
    			logger.warn("Counts array does not contain key: "+chr+"\n");
    			break;    // assume that all exons are from the same chromosome (which wasn't found)
    		}
    	}
    	return countSum;
    }
    
    public LinkedList<Window> countWindowedReadStarts(final Set<Annotation> eSet, 
    		final HashMap<String, HashMap<String,IntervalTree<String>>> iTree,
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
    	float[] floatCounts; 
    	int[] gCoords;    // genomic coordinates of each location
    	int eStart = 0;
    	int eEnd = 0;
    	int eLen = 0; 
    	String gStrand;   // strand of the gene
    	int localExtend = extend;    // default extension
    	IntervalTree<String> exonTree = new IntervalTree<String>();
		String oStr = null;

    	LinkedList<Window> wList = new LinkedList<Window>();    // Is this the right type?...

    	// get the gene strand
    	if (gene.isNegativeStrand()) {
    		gStrand = "-";
    	} else {
    		gStrand = "+";
    	}

    	/* get the BED file segment to which this segment is aligned: */
    	Iterator<Annotation> eIter = eSet.iterator();    // iterator over the exons
    	Annotation e = eIter.next(); 
    	String chr = e.getChr();
    	String strand = e.getStrand().toString();
    	int nExons = eSet.size();      // number of exons for this gene
    	
    	// reset the iterator:
    	eIter = eSet.iterator();

    	/* Step 1: determine the start and end coordinate range for this gene, and the total number
    	 * of bases over all exons.
    	 */
    	int gMin = Integer.MAX_VALUE;
    	int gMax = 0;
    	int tLen = 0;   // transcript length, in bases
    	int eID = 0;    // exon ID

    	while (eIter.hasNext()) {
    		e = eIter.next();
    		eStart = e.getStart();
    		eEnd = e.getEnd();

    		// add an interval to the interval tree:
    		exonTree.put(eStart,  eEnd, gene.getName()+".e"+eID);
    		eID++;
    		
    		/* find earliest start and latest end */
    		if (eStart<gMin) {
    			gMin = eStart;
    		}
    		if (eEnd>gMax) {
    			gMax = eEnd;
    		}
    		
    		// update the transcript length
    		tLen += eEnd-eStart;
    	}
    	
    	/* Step 2: Determine the maximum length, up to "extend" bases, that the transcript can be extended without
    	 * overlapping a nearby gene. This will be the new value of "localExtend"
    	 */
    	// if ((gene.isNegativeStrand() && <task is 3prime>) || (!gene.isNegativeStrand() && <task is 5prime>))  // eventual test 
    	if (gStrand.equals("-")) {
    		/* "-" strand (with 3prime libraries, or "+" strand with 5prime libraries) */
    		eStart = exonTree.min().getStart();           // "leftmost" transcript genomic coordinate
    		int minExt = Math.max(0, eStart-extend);      // initial minimum of extension
    		int maxExt = eStart;
    		if (iTree.get(strand).get(chr).numOverlappers(minExt, maxExt) > 0) {
    			// trim overlap
    			Iterator<Node<String>> iIter = iTree.get(strand).get(chr).overlappers(minExt, maxExt);
    			oStr = null;
    			while (iIter.hasNext()) {
    				Node<String> n = iIter.next();
    				// construct string of names of overlappers
    				String iName = n.getValue();
    				if (oStr==null) {
    					oStr = iName;
    				} else {
    					oStr += ","+iName;
    				}
    				int iEnd = n.getEnd();
    				if (iEnd>minExt) {
    					minExt = Math.min(iEnd,maxExt);   
    				}
    			}
    		}
			localExtend = maxExt-minExt;   // final extension value
			if (localExtend<extend) {
				logger.warn("Gene "+gene.getName()+" ("+strand+") extension overlaps "+oStr+". "+extend+"-base extension shortened to "+localExtend);
			}
			// if the there is an extension past the transcript, andd an interval to the tree:
			if (localExtend>0) {
				exonTree.put(minExt, maxExt, gene.getName()+".ext");
			}
    	} else {
    		/* "+" strand (with 3prime libraries, or "-" strand with 5prime libraries) */
    		eEnd = exonTree.max().getEnd();             // "rightmost" transcript genomic coordinate
    		int minExt = eEnd;
    		int maxExt = Math.min(eEnd+extend, getChrLength(chr));  // initial maximum of extension 
    		if (iTree.get(strand).get(chr).numOverlappers(minExt, maxExt) > 0) {
    			// trim overlap
    			Iterator<Node<String>> iIter = iTree.get(strand).get(chr).overlappers(minExt, maxExt);
    			oStr = null;
    			while (iIter.hasNext()) {
    				Node<String> n = iIter.next();
    				// construct string of names of overlappers
    				String iName = n.getValue();
    				if (oStr==null) {
    					oStr = iName;
    				} else {
    					oStr += ","+iName;
    				}
    				int iStart = n.getStart();
    				if (iStart<maxExt) {
    					maxExt = Math.max(iStart,minExt);   
    				}
    			}
    		}
			localExtend = maxExt-minExt;   // final extension value
			if (localExtend<extend) {
				logger.warn("Gene "+gene.getName()+" ("+strand+") extension overlaps "+oStr+". "+extend+"-base extension shortened to "+localExtend);
			}
			// if the there is an extension past the transcript, andd an interval to the tree:
			if (localExtend>0) {
				exonTree.put(minExt, maxExt, gene.getName()+".ext");
			}
    	}
    	
    	/* Step 3: Copy the counts and genomic coordinates for all locations within
    	 * this gene's transcript, plus the extension.
    	 */
    	/* allocate the storage */
    	floatCounts = new float[tLen+localExtend];
    	gCoords = new int[tLen+localExtend];
    	int cStart = 0;     // index of where to start copying

    	/* If 3' library and negative strand, or 5' library and positive strand copy counts and genomic coordinates 
    	 * for the extension to the "left" of the transcript start: 
    	 */
    	// if ((gene.isNegativeStrand() && <task is 3prime>) || (!gene.isNegativeStrand() && <task is 5prime>))  // eventual test 
    	if (gStrand.equals("-")) {
    		/* "-" strand (with 3prime libraries, or "+" strand with 5prime libraries) */
    		eStart = gMin-localExtend;
    		eLen = localExtend;
			copyToLocalCounts(chr, gStrand, eStart, cStart, eLen, floatCounts);
			for (int i=0; i<localExtend; i++) {
				gCoords[cStart+i] = eStart+i;      // corresponding genomic coordinates of the counts
    		}
			cStart+=localExtend;
    	}
    	
    	/* For all exons, copy the counts and genomic coordinates */
    	eIter = eSet.iterator();
    	while (eIter.hasNext()) {
    		e = eIter.next();
    		eStart = e.getStart();
    		eLen = e.getEnd()-e.getStart();
			copyToLocalCounts(chr, gStrand, eStart, cStart, eLen, floatCounts);
    		for (int i=0; i<eLen; i++) {
    			gCoords[cStart+i] = eStart+i;
    		}
    		cStart+=eLen;
    	}
    	
    	/* If 3' library and positive strand, or 5' library and negative strand copy counts and genomic coordinates 
    	 * for the extension to the "right" of the transcript end: 
    	 */
    	// if ((!gene.isNegativeStrand() && <task is 5prime>) || (gene.isNegativeStrand() && <task is 3prime>))  // eventual test 
    	if (!gene.isNegativeStrand()) {
    		/* "+" strand (with 3prime libraries, or "-" strand with 5prime libraries) */
    		eStart = gMax;   // +1?
    		eLen = localExtend;
			copyToLocalCounts(chr, gStrand, eStart, cStart, eLen, floatCounts);
			for (int i=0; i<localExtend; i++) {
				gCoords[cStart+i] = eStart+i;      // corresponding genomic coordinates of the counts
    		}
			cStart+=localExtend;
    	}
    	
    	/* Iterate a sliding window across the counts array and sum the counts over the window. Create a new Window object
    	 * for each step and add them to the output list.
    	 */
    	int sumStart = 0;					// start of summing window
    	int sumEnd = sumStart+window;       // end of summing window (exclusive)

    	aLen = tLen+localExtend;    // this is the total number of count locations being considered
    	sumStart = 0;
    	sumEnd = sumStart+window;

    	if (aLen<=window) {
    		// deal with special case of the total gene length (plus extension) being shorter than a window length:
    		sumEnd = aLen-1;     // number of valid reads, minus one to compensate for 0-based index;
			float countSum = 0;
			for (int i=sumStart; i<sumEnd; i++) {
				countSum += floatCounts[i];
			}
			// don't bother saving windows with zero counts:
			if (countSum>0.0) {
	    		Window thisWindow = new Window(gStrand, chr, gCoords[sumStart], gCoords[sumEnd], gene.getName());
				thisWindow.setCount(countSum);    // update count for this window
				// if the gene has more than one exon, check to see if the window spans more than one
				if (nExons > 1) {
					thisWindow.addIntervals(gCoords[sumStart], gCoords[sumEnd], exonTree);
				}
				wList.add(thisWindow);            // add this window to the output list
			}
    	} else {
    		while (sumEnd < aLen) {
    			float countSum = 0;
    			for (int i=sumStart; i<sumEnd; i++) {
    				countSum += floatCounts[i];
    			}
    			// don't bother saving windows with zero counts:
    			if (countSum>0.0) {
        			Window thisWindow = new Window(gStrand, chr, gCoords[sumStart], gCoords[sumEnd], gene.getName());
        			thisWindow.setCount(countSum);    // update count for this window
    				// if the gene has more than one exon, check to see if the window spans more than one
    				if (nExons > 1) {
    					thisWindow.addIntervals(gCoords[sumStart], gCoords[sumEnd], exonTree);
    				}
    				wList.add(thisWindow);            // add this window to the output list
    			}
    			/* update start and end for the next window */
    			sumStart += (window-overlap);
    			sumEnd = sumStart+window;
    			/* NB: Any partial windows will be skipped! */
    		}
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
    	float eCount = 0;  
    	
    	// Iterate over all "chromosomes":
    	for(String chr:annotations.keySet()){
			// Iterate over all transcripts:
    		for(Gene gene : annotations.get(chr)) {
				Set eSet = gene.getExonSet();
				// sum all counts starting within this exon set:
				eCount = countExonReadStarts(eSet);
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
    	// countsMap[chr][transID][Window]
    	HashMap<String,HashMap<String,LinkedList<Window>>> countsMap = 
    			new HashMap<String,HashMap<String,LinkedList<Window>>>();   
    	LinkedList<Window> eCount;  
    	// iTree[strand][chromosome]<tree>
    	HashMap<String, HashMap<String, IntervalTree<String>>> iTree = new HashMap<String, HashMap<String, IntervalTree<String>>>();
    	iTree.put("+", new HashMap<String, IntervalTree<String>>());
    	iTree.put("-", new HashMap<String, IntervalTree<String>>());

		// Iterate over all "chromosomes":
    	for(String chr:annotations.keySet()){
    		/* Build strand-specific IntervalTrees containing all exons for this segment */
    		iTree.get("+").put(chr, new IntervalTree<String>());
    		iTree.get("-").put(chr, new IntervalTree<String>());
			// Iterate over all transcripts:
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
						iTree.get("+").get(chr).put(eStart, eEnd, intervalName);
					} else {
						iTree.get("-").get(chr).put(eStart, eEnd, intervalName);
					}						
					eID++;	// increment the exon ID
				}
    		}

    		// Count windowed read starts over all genes:
    		for(Gene gene : annotations.get(chr)) {
    			String gName = gene.getName();
				Set eSet = gene.getExonSet();
				
				// skip any chromosomes/segments that were not observed in the alignments files:
				if (!startCountsHasKey(chr)) {
					continue;
				}
				
				// test to see if this is a duplicate transcript:
				if (countsMap.containsKey(chr)) {	
					if (countsMap.get(chr).containsKey(gName)) {
						logger.warn("Duplicate entry for transcript "+gName+". Skipping...");	
						continue;
					}
				}

				// sum all counts starting within this exon set:
				eCount = countWindowedReadStarts(eSet, iTree, window, overlap, extend, gene);
				
				// add this window set to the HashMap:
				if (!countsMap.containsKey(chr)) {
					countsMap.put(chr, new HashMap<String,LinkedList<Window>>());           // add an entry in the main map for the chromosome/segment
				}
				countsMap.get(chr).put(gName, eCount);       // add the counts windows to the hashmap
			}	
		}
    	
    return countsMap;
    	
    }    
    
    public int getChrLength(String chr) {
    	//return startCounts.get(chr).length;
    	return getSequence(chr).getSequenceLength();
    }
    

    abstract void incrementStartCounts(String refName, String strand, int alignStart, float fractCount);
	abstract void copyToLocalCounts(String chr, String strand, int eStart, int cStart, int eLen, float[] floatCounts);
    abstract public void updateCount(SAMRecord r);
    abstract boolean startCountsHasKey(String chr);
    abstract float getStartCounts(String chr, String strand, int i);
}
