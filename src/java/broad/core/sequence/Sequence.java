package broad.core.sequence;

import broad.core.parser.StringParser;

import java.security.InvalidParameterException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.BasicAnnotation;
import nextgen.core.annotation.Annotation.Strand;
import nextgen.core.annotation.Gene;

import org.apache.log4j.Logger;

import Jama.Matrix;
import broad.core.annotation.BasicGenomicAnnotation;
import broad.core.annotation.GenomicAnnotation;
import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.pda.datastructures.Alignments;
import broad.pda.seq.graph.ChromosomeWithBubblesJGraphT.WindowIterator;

public class Sequence {
	static Logger logger = Logger.getLogger(Sequence.class.getName());
	private static final Pattern SOFT_MASKED_PAT = Pattern.compile("[acgt]+");
	private static final Pattern SEQUENCE_GAP_PAT = Pattern.compile("[Nn]+");
	private static final Pattern UNGAPPED_SEQUENCE_PAT = Pattern.compile("[^Nn]+");
	
	private String id;
	private StringBuilder sequenceBases;
	private short[] encodedSequence;
	private Matrix vectorEncodedSequence;
	private boolean forwardStrand = true;
	private boolean encodeIgnoreCase = false;
	
	public static final char[] SHORT_READ_FLOW = {'T','A','C','G'};
	public static final char[] LONG_READ_FLOW  = {'G','A','T','C'};
	
	public static final short SHORT_ENCODED_A = 0;
	public static final short SHORT_ENCODED_C = 1;
	public static final short SHORT_ENCODED_G = 2;
	public static final short SHORT_ENCODED_T = 3;
	public static final short SHORT_ENCODED_GAP = 4;
	public static final short SHORT_ENCODED_a = 5;
	public static final short SHORT_ENCODED_c = 6;
	public static final short SHORT_ENCODED_g = 7;
	public static final short SHORT_ENCODED_t = 8;
	public static final short SHORT_ENCODED_N = 9;

	public Sequence(String id) {
		super();
		this.id = id;
		sequenceBases = new StringBuilder();
	}
	
	public Sequence(String id, boolean isLarge) {
		this(id);
		if(isLarge) {
			sequenceBases = new StringBuilder(850000000);
		}
	}
	
	public Sequence(String id, int expectedSize) {
		this(id);
		sequenceBases = new StringBuilder(expectedSize);

	}
	


	public String getId() {
		return id;
	}
	
	public void setForwardStrand(boolean isForwardStrand) { this.forwardStrand = isForwardStrand;}

	public void setId(String id) {
		this.id = id;
	}
	
	public void unloadSequence() {
		sequenceBases.delete(0, sequenceBases.length());
		sequenceBases.trimToSize();
	}
	
	public void unloadEncodedSequence() {
		encodedSequence = null;
	}
	
	public void unloadAllSequences () {
		unloadSequence();
		unloadEncodedSequence();
	}
	
	public int getLength() {
		return sequenceBases.length();
	}
	
	public void setCharAt (int position, char newCharacter) throws IllegalAccessException {
		if (sequenceBases == null) {
			throw new IllegalAccessException("This methods is only implemented for un encoded sequences. sequenceBases is null, it must be non null.");
		}
		
		sequenceBases.setCharAt(position, newCharacter);
	}

	public String getSequenceBases() {
	    String bases = "";
		if(sequenceBases.length() > 0) {
			bases =  sequenceBases.toString();
		} else if (encodedSequence != null && encodedSequence.length > 0) {
			StringBuffer buf = new StringBuffer(encodedSequence.length);
			if(!encodeIgnoreCase) {
				for(int i = 0 ; i < encodedSequence.length; i++) {
					switch (encodedSequence[i]) {
					case SHORT_ENCODED_c : buf.append("c"); break;
					case SHORT_ENCODED_C : buf.append("C"); break;
					case SHORT_ENCODED_G : buf.append("G"); break;
					case SHORT_ENCODED_g : buf.append("g"); break;
					case SHORT_ENCODED_a : buf.append("a"); break;
					case SHORT_ENCODED_A : buf.append("A"); break;
					case SHORT_ENCODED_T : buf.append("T"); break;
					case SHORT_ENCODED_t : buf.append("t"); break;
					case SHORT_ENCODED_GAP : buf.append("-"); break;
					default : buf.append("N"); break;
					}
				}

			}else {
				for(int i = 0 ; i < encodedSequence.length; i++) {
					switch (encodedSequence[i]) {
					case SHORT_ENCODED_A : buf.append("A"); break;
					case SHORT_ENCODED_C : buf.append("C"); break;
					case SHORT_ENCODED_G : buf.append("G"); break;
					case SHORT_ENCODED_T : buf.append("T"); break;
					case SHORT_ENCODED_GAP : buf.append("-"); break;
					default : buf.append("N"); break;
					}
				}
			}
			bases = buf.toString();
		} else if (vectorEncodedSequence != null) {
			Random random = new Random();
			int seqLength = vectorEncodedSequence.getColumnDimension();
			StringBuffer buf = new StringBuffer(seqLength);
			for(int j = 0; j < seqLength; j++) {
				double draw = random.nextDouble();
				int drawedBase = -1;
				double cummulativeProbability = 0;
				for(int i = 0; i < 4; i++) {
					cummulativeProbability += vectorEncodedSequence.get(i,j);
					if(draw <= cummulativeProbability ) {
						drawedBase = i;
						break;
					}
				}
				
				switch (drawedBase) {
					case SHORT_ENCODED_A : buf.append("A"); break;
					case SHORT_ENCODED_C : buf.append("C"); break;
					case SHORT_ENCODED_G : buf.append("G"); break;
					case SHORT_ENCODED_T : buf.append("T"); break;
					default : buf.append("-"); break;				
				}
			}
			bases = buf.toString();
		}
		
		return bases;
	}

	public void setSequenceBases(String sequence) {
		this.sequenceBases = new StringBuilder(sequence);
	}
	
	public void append(String partialSequence) {
		sequenceBases.append(partialSequence);
	}
	
	public Sequence shuffle() {
		Sequence shuffle = new Sequence(getId()+"_shuffle");
		List<Integer> indeces = new ArrayList<Integer>(sequenceBases.length());
		for(int i = 0; i < sequenceBases.length(); i++) {
			indeces.add(i);
		}
		Random r = new Random();
		while(indeces.size() > 0) {
			int idx = indeces.remove(r.nextInt(indeces.size()));
			shuffle.append(sequenceBases.charAt(idx));
		}
		return shuffle;
	}
	
	public void append(char c) {
		sequenceBases.append(c);
	}
	
	public boolean isGap(int position) {
		boolean isGap = false;
		if(encodedSequence != null && encodedSequence.length >= position) {
			isGap =  encodeIgnoreCase && (encodedSequence[position] == 4);
		} else if (vectorEncodedSequence != null && vectorEncodedSequence.getColumnDimension() >= position) {
			double maxProb = 0;
			
			for(int i = 0; i < 4; i++) {
				double baseProbability = vectorEncodedSequence.get(i,position);
				if(maxProb < baseProbability ) {
					maxProb = baseProbability;
				}
			}
			
			isGap = maxProb == 0;
		} else if (sequenceBases != null & sequenceBases.length() >= position) {
			isGap = '-' == sequenceBases.charAt(position);
		}
		
		return isGap;
	}
	
	public List<GenomicAnnotation> getSoftmaskedRegions() {
		List<GenomicAnnotation> softMaskedRegions = new ArrayList<GenomicAnnotation>();
		if(sequenceBases == null) {
			return softMaskedRegions;
		}
		
		Matcher m = SOFT_MASKED_PAT.matcher(sequenceBases);
		int i = 1;
		while(m.find()) {
			BasicGenomicAnnotation softMaskedReg = new BasicGenomicAnnotation(getId() + "_SoftmaskedReg_" + i++);
			softMaskedReg.setStart(m.start() + 1);
			softMaskedReg.setEnd(m.end() + 1);
			softMaskedReg.setChromosome(getId());
			softMaskedRegions.add(softMaskedReg);
		}
		return softMaskedRegions;
	}
	
	public int countSoftMaskedBases() {
		int count = 0;
		Matcher m = SOFT_MASKED_PAT.matcher(sequenceBases);
		while(m.find()) {
			count += m.end() - m.start() + 1;
		}
		return count;
	}
	
	public void maskSoftmaskedRegions() {
		if(sequenceBases != null) {
			Matcher m = SOFT_MASKED_PAT.matcher(sequenceBases);
			while(m.find()) {
				StringBuilder ns = new StringBuilder(m.end() - m.start());
				for(int j = 0; j < m.end() - m.start(); j++) { ns.append("N");} 
				sequenceBases.replace(m.start(), m.end(), ns.toString());
			}
		}
	}
	
	public short[] getEncodedSequence () {
		return encodedSequence;
	}
	
	public Matrix getVectorEncodedSequence() {
		return vectorEncodedSequence;
	}
	
	public Matrix encodeSequenceAsVector() {
		
		vectorEncodedSequence = new Matrix(4, sequenceBases.length());
		
		for(int j = 0; j < sequenceBases.length(); j++) {
			char c = sequenceBases.charAt(j);
			if('a' == c || 'A' == c) {
				vectorEncodedSequence.set(SHORT_ENCODED_A, j, 1);
			} else if ('C' == c || 'c' == c) {
				vectorEncodedSequence.set(SHORT_ENCODED_C, j, 1);
			}else if ('G' == c || 'g' == c) {
				vectorEncodedSequence.set(SHORT_ENCODED_G, j, 1);
			}else if ('T' == c || 't' == c) {
				vectorEncodedSequence.set(SHORT_ENCODED_T, j, 1);
			} 
		}	

		return vectorEncodedSequence;
	}
	
	public short[] encodeSequence() {
		if(encodedSequence != null) {
			return encodedSequence;
		}

		this.encodedSequence = encodeSequence(sequenceBases);
		encodeIgnoreCase = false;
		return this.encodedSequence;
	}

	public static short[] encodeSequence(StringBuilder bases) {
		short [] encodedSeq = new short[bases.length()];
		
		for(int i = 0; i < bases.length(); i++) {
			char c = bases.charAt(i);
			if('c' == c) {
				encodedSeq[i] = SHORT_ENCODED_c;
			}else if ('C' == c) {
				encodedSeq[i] = SHORT_ENCODED_C;
			}else if ('G' == c) {
				encodedSeq[i] = SHORT_ENCODED_G;
			}else if ('g' == c) {
				encodedSeq[i] = SHORT_ENCODED_g;
			}else if ('a' == c) {
				encodedSeq[i] = SHORT_ENCODED_a;
			}else if ('A' == c) {
				encodedSeq[i] = SHORT_ENCODED_A;
			}else if('t' == c) {
				encodedSeq[i] = SHORT_ENCODED_t;
			}else if('T' == c) {
				encodedSeq[i] = SHORT_ENCODED_T;
			}else {
				encodedSeq[i] = SHORT_ENCODED_N;
			}
		}			
		return encodedSeq;
	}
	
	public short[] encodeSequenceIgnoreCase( ) {
		return encodeSequenceIgnoreCase(false);
	}
	
	public short[] encodeSequenceIgnoreCase(boolean distinguishMissingSequence ) {
		if(encodedSequence != null) {
			return encodedSequence;
		}
		this.encodedSequence = encodeSequenceIgnoreCase(sequenceBases, distinguishMissingSequence);
		encodeIgnoreCase = true;
		return this.encodedSequence;
		
	}
	
	public static short [] encodeSequenceIgnoreCase(StringBuilder sequence) {
		return encodeSequenceIgnoreCase(sequence, false);
	}
	
	public static short [] encodeSequenceIgnoreCase(StringBuilder sequence, boolean distinguishMissingSequence) {
		short [] encodedSeq = new short[sequence.length()];
		
		for(int i = 0; i < sequence.length(); i++) {
			char c = sequence.charAt(i);
			if('a' == c || 'A' == c) {
				encodedSeq[i] = SHORT_ENCODED_A;
			}else if ('C' == c || 'c' == c) {
				encodedSeq[i] = SHORT_ENCODED_C;
			}else if ('G' == c || 'g' == c) {
				encodedSeq[i] = SHORT_ENCODED_G;
			}else if ('T' == c || 't' == c) {
				encodedSeq[i] = SHORT_ENCODED_T;
			} else if ('-' ==  c) {
				encodedSeq[i] = SHORT_ENCODED_GAP;
			} else if (distinguishMissingSequence && ('N' == c || 'n' == c) ){
				encodedSeq[i] = SHORT_ENCODED_N;
			} else {
				encodedSeq[i] = SHORT_ENCODED_GAP;  //All other letter codes are cosidered gaps.
				//System.err.println("Unexpected letter " + Character.toString(c) );
			}
		}			
		return encodedSeq;
	}
	
	
	public List<Short> compute454Flow(char[] flowOrder) {
		return compute454Flow(flowOrder, 1, sequenceBases.length());
	}

	public List<Short> compute454Flow(char[] flowOrder, int start, int end) {
		char [] subSeqChrs = new char[end - start];
		sequenceBases.getChars(start - 1, end - 1, subSeqChrs, 0);
		ArrayList<Short> flow = new ArrayList<Short>();
		int seqIdx = 0;
		int flowIdx = 0;
		char flowLetter = 0;
		short flowRead = 0;
		while(seqIdx < subSeqChrs.length) {
			flowLetter = flowOrder[flowIdx % 4];
			char seqChar = Character.toUpperCase(subSeqChrs[seqIdx]);
			if(seqChar != 'A' && seqChar != 'T' && seqChar != 'G' && seqChar != 'C') {
				System.err.println("Sequence Character " + seqChar + " is not a base, ignoring it");
				seqIdx++;
			}
			if(flowLetter == seqChar) {
				flowRead++;
				seqIdx++;
			} else {
				flow.add(flowRead);
				flowIdx++;
				flowRead = 0;
			}

		}
		// TODO Auto-generated method stub
		return flow;
	}
	
	public float gcContent() {
		return computeGCContent(sequenceBases.toString());
	}
	
		
	public boolean contains(String motif) {
		Sequence revMotifSeq = new Sequence("rev_motif");
		revMotifSeq.setSequenceBases(motif);
		revMotifSeq.reverse();
		String revMotif = revMotifSeq.getSequenceBases();
		
		return getSequenceBases().contains(motif) || getSequenceBases().contains(revMotif);
	}
	
	public static float computeGCContent(String dnaString) {
		int gcs = 0;
		for(int i = 0; i < dnaString.length(); i++) {
			char c = Character.toUpperCase(dnaString.charAt(i));
			if('C' == c || 'G' == c) {
				gcs++;
			}
		}
		
		return ((float)gcs)/((float)dnaString.length());
	}
	
	
	public static float [] computeSequenceComposition(String dnaString) {
		int gs = 0;
		int cs = 0;
		int as = 0;
		int ts = 0;
		int ns = 0;
		for(int i = 0; i < dnaString.length(); i++) {
			char c = Character.toUpperCase(dnaString.charAt(i));
			switch (c) {
			case 'A' : as++;
				break;
			case 'C' : cs++;
				break;
			case 'G':  gs++;
				break;
			case 'T' : cs++;
				break;
			default :  ns++;
			}
		}
		int total = gs + cs + as + ts ;
		float [] composition = {as/(float) total, cs/(float)total, gs/(float)total, ts/(float) total };
		return composition;
	}
	
	public void reverse() {
		String seqBases = getSequenceBases();
		StringBuilder reversedSeq = new StringBuilder(seqBases.length());
		for(int j = seqBases.length() - 1; j >= 0 ; j--) {
			char c = seqBases.charAt(j);
			if('c' == c) {
				reversedSeq.append('g');
			}else if ('C' == c) {
				reversedSeq.append('G');
			}else if ('G' == c) {
				reversedSeq.append('C');
			}else if ('g' == c) {
				reversedSeq.append('c');
			}else if ('a' == c) {
				reversedSeq.append('t');
			}else if ('A' == c) {
				reversedSeq.append('T');
			}else if('t' == c) {
				reversedSeq.append('a');
			}else if('T' == c) {
				reversedSeq.append('A');
			}else if('N'==c){
				reversedSeq.append('N');
			}else if('n'==c){
				reversedSeq.append('n');
			}else {
				reversedSeq.append(c);
			}
		}
		
		sequenceBases = reversedSeq;
		if(encodedSequence != null) {
			encodedSequence = null;
			encodeSequenceIgnoreCase();
		} else if(vectorEncodedSequence != null) {
			vectorEncodedSequence = null;
			encodeSequenceAsVector();
		}
		forwardStrand = false;
	}
	
	public static Sequence reverseSequence(Sequence seq) {
		String r = reverseSequence(seq.getSequenceBases());
		Sequence rtrn = new Sequence(seq.getId() + "_reverse");
		rtrn.setSequenceBases(r);
		return rtrn;
	}
	
	public static String reverseSequence(String seq) {
		Sequence tmpSeq = new Sequence("tmp");
		tmpSeq.setSequenceBases(seq);
		tmpSeq.reverse();
		return tmpSeq.getSequenceBases();
	}
	
	public void uppercase() {
		StringBuilder uppercasedSeq = new StringBuilder(sequenceBases.toString().toUpperCase());
		sequenceBases = uppercasedSeq;
	}
	
	/**
	 * Mask part of sequence
	 * @param region Region to mask. Masks blocks only.
	 * @param softmask If true, mask to lower case. If false, mask to Ns.
	 */
	public void mask(Annotation region, boolean softmask) {
		if(!region.getChr().equals(id)) {
			throw new IllegalArgumentException("Region chromosome (" + region.getChr() + " does not match sequence name " + id);
		}
		for(Annotation block : region.getBlocks()) {
			int start = block.getStart();
			int end = block.getEnd();
			if(softmask) {
				String origStr = getSequenceBases().substring(start, end);
				String lower = origStr.toLowerCase();
				sequenceBases.replace(start, end, lower);
				//logger.info("Changed " + origStr + " to " + lower + ". Now: " + getSubSequence("", start, end).getSequenceBases());
			} else {
				int len = end - start;
				String ns = "";
				for(int i = 0; i < len; i++) {
					ns += "N";
				}
				sequenceBases.replace(start, end, ns);
			}
			
		}
	}
	
	/**
	 * 
	 * @param start - start of the region to extract the starting base will be included starting at 0
	 * @param end - the end of the region to extract, the base at position end will not be included.
	 * @return
	 */
	public SequenceRegion getRegion(int start, int end) {
		SequenceRegion region = new SequenceRegion(getId());
		region.setRegionStart(start);
		region.setRegionEnd(end);
		String bases = getSequenceBases();
		if(bases != null && bases.length() > 0) {
			region.setSequenceBases(bases.substring(start, end));
		}
		
		return region;
	}
	
	public void getRegion(SequenceRegion region) {
		region.setSequenceBases(sequenceBases.substring(region.getStart(), region.getEnd()));
		if(region.inReversedOrientation()) {
			region.reverse();
			//System.err.println("Sequence reversed");
		}
	}
	
	public void getRegion(SequenceRegion region, boolean softmask) {
		//if(region!=null){region.setSequenceBases(""); System.err.println("NULL");}
		region.setSequenceBases(sequenceBases.substring(region.getStart(), region.getEnd()));
		if(!softmask){return;}
		String seq=region.getSequenceBases();
		char[] chars=seq.toCharArray();
		String temp="";
		for(int i=0; i<chars.length; i++){
			if(chars[i]=='a' || chars[i]=='c' || chars[i]=='g' || chars[i]=='t'){temp+="N";}
			else{temp+=chars[i];}
		}
		region.setSequenceBases(temp);
	}
	
	public void getRegion(SequenceRegion region, boolean softmask, Map<String, IntervalTree<Alignments>> okRepeats) {
		//if(region!=null){region.setSequenceBases(""); System.err.println("NULL");}
		region.setSequenceBases(sequenceBases.substring(region.getStart(), region.getEnd()));
		
		if(!softmask){return;}
		
		//System.err.println(region.getChromosome()+"\t"+region.getStart()+"\t"+region.getEnd());
		
		String seq=region.getSequenceBases();
		char[] chars=seq.toCharArray();
		String temp="";
		for(int i=0; i<chars.length; i++){
			if(chars[i]=='a' || chars[i]=='c' || chars[i]=='g' || chars[i]=='t'){
				boolean mask=checkIfMask(i+region.getStart(), okRepeats.get(region.getChromosome())); if(mask){temp+="N";} else{temp+=chars[i];}
			}
			else{temp+=chars[i];}
		}
		region.setSequenceBases(temp);
	}
	
	private boolean checkIfMask(int index, IntervalTree<Alignments> okRepeats){
		boolean mask=true;
		if(okRepeats==null || okRepeats.isEmpty()){return true;}
		
		Iterator<Node<Alignments>> overlappers=okRepeats.overlappers(index, index+1);
		mask=!overlappers.hasNext();
		
		return mask;
	}
	
	public void getRegions(List<? extends SequenceRegion> regions) {
		Iterator<? extends SequenceRegion> it = regions.iterator();
		while(it.hasNext()) {
			SequenceRegion reg = it.next();
			getRegion(reg);
		}
		
	}
	
	public SequenceRegion extractRegionBasedOnGC(float targetGC, int size, int buffer) {
		SequenceRegion theRegion = null;
		float closestGC = 0;

		System.out.println("extracting region based on GC for sequence " + getId() + " of size " + getSequenceBases().length() + " buffer " + buffer);
		
		if(getSequenceBases() == null || getSequenceBases().length() < buffer) {
			return theRegion;
		}

		for(int i = buffer; i < (getSequenceBases().length() - buffer - size); i += buffer) {
			String currentSequence = getSequenceBases().substring(i - buffer, i + size + buffer - 1);
			float currentGC = computeGCContent(currentSequence);
			System.out.println("position " + i + " currentGC " + currentGC + ", targetGC " + targetGC + ",  closestGC distance to target " + Math.abs(closestGC - targetGC) + " currentGC distance to target " + Math.abs(targetGC - currentGC) );
			if(Math.abs(closestGC - targetGC) > Math.abs(targetGC - currentGC) ){
				closestGC = currentGC;
				theRegion = new SequenceRegion(getId());
				theRegion.setSequenceBases(currentSequence);
				theRegion.setRegionStart(i - buffer);
				theRegion.setEnd(i + size + buffer - 1);
			}
		}
		System.out.println("Returning closest GC region: " + theRegion + " GC% " + closestGC);
		return theRegion;
	}

	public boolean isForwardStrand() { return forwardStrand;}

	public StringBuilder getSequenceBuilder() {
		return sequenceBases == null || sequenceBases.length() == 0 
			? new StringBuilder(getSequenceBases()) 
			: sequenceBases;
	}
	
	public int getGapsSize() {
		int totalGaps = 0;
		if(sequenceBases != null && sequenceBases.length() > 0) {
			char [] baseArray = sequenceBases.toString().toUpperCase().toCharArray();
			for(int i = 0; i < baseArray.length; i++) {
				if(baseArray[i] == 'N') {
					totalGaps++;
				}
			}
		}
		return totalGaps;
	}
	
	public WindowSlider getSlider(int windowSize, int overlap) {
		return WindowSlider.getSlider(this, windowSize, overlap);
	}
	
	public List<SequenceRegion> chunk(int chunkSize, int chunkOverlap) {
		int numOfChunks = (int) Math.floor(getLength()/((float)(chunkSize - chunkOverlap)));
		int chunkStart = 0;
		int chunkEnd   = 0;
		System.out.println("\tChunking " + getId() + " number of chunks " + numOfChunks);
		List<SequenceRegion> chunks = new ArrayList<SequenceRegion>(numOfChunks);
		
		for(int i = 0; i< numOfChunks; i++) {
			chunkStart = i * (chunkSize - chunkOverlap);
			chunkEnd   = chunkStart + chunkSize;
			SequenceRegion chunk = getRegion(chunkStart, chunkEnd);
			chunks.add(chunk);
		}
		if(chunkEnd < getLength()) {
			chunkStart = chunkEnd - chunkOverlap;
			SequenceRegion lastChunk = getRegion(chunkStart, getLength());
			chunks.add(lastChunk);
		}
		
		System.out.println("\tObtained " + chunks.size() + " chunks"); 
		return chunks;
	}

	public int getEnd() {
		return getLength();
	}
	
	/**
	 * Finds the list of ungapped regions in the sequence.
	 * @return A list of two sized integer arrays with the start and end of the ungapped regions. The convention here is semi closed
	 * 		   intervals, each list items is of the form [start, end).
	 */
	public List<int []> findUngappedSequenceChunks() {
		Matcher m = UNGAPPED_SEQUENCE_PAT.matcher(getSequenceBases());
		ArrayList<int []> gapUngapp = new ArrayList<int []>();
		//System.out.println("Sequence for " + getId() + ": " + getSequenceBases());
		while(m.find()) {
			int [] ungappedReg = {m.start(), m.end()};
			//System.out.println("\treg["+m.start()+"-"+m.end()+"]: " + getSequenceBases().substring(m.start(), m.end()));
			gapUngapp.add(ungappedReg);
		}
		//System.out.println("");

		return gapUngapp;			
	}
	
	public static String decode (short [] shortEncodedSequence) {
		StringBuilder buf = new StringBuilder();
		for(int i = 0 ; i < shortEncodedSequence.length; i++) {
			switch (shortEncodedSequence[i]) {
			case SHORT_ENCODED_c : buf.append("c"); break;
			case SHORT_ENCODED_C : buf.append("C"); break;
			case SHORT_ENCODED_G : buf.append("G"); break;
			case SHORT_ENCODED_g : buf.append("g"); break;
			case SHORT_ENCODED_a : buf.append("a"); break;
			case SHORT_ENCODED_A : buf.append("A"); break;
			case SHORT_ENCODED_T : buf.append("T"); break;
			case SHORT_ENCODED_t : buf.append("t"); break;
			case SHORT_ENCODED_GAP : buf.append("-"); break;
			default : buf.append("N"); break;
			}
		}
		return buf.toString();
	}
	
	/**
	 * Get all matches for a string as annotations with respect to this sequence
	 * Search a window
	 * Case is ignored
	 * @param s The string
	 * @param windowStart Start position of window to search
	 * @param windowEnd End of window to search
	 * @param includeRC Also search for reverse complement of string
	 * @return Collection of annotations representing string matches
	 */
	public Collection<Annotation> getMatches(String s, int windowStart, int windowEnd, boolean includeRC) {
		Sequence window = getSubSequence("", windowStart, windowEnd);
		String seqString = window.getSequenceBases().toUpperCase();
		String upper = s.toUpperCase();
		Collection<Annotation> rtrn = new ArrayList<Annotation>();
		int i = 0;
		while(i < seqString.length()) {
			int pos = seqString.indexOf(upper, i);
			if(pos == -1) {
				break;
			}
			if(pos + s.length() > seqString.length()) {
				// Won't be fully contained
				break;
			}
			// Create annotation
			Annotation annot = new BasicAnnotation(getId(), windowStart + pos, windowStart + pos + s.length(), Strand.POSITIVE);
			rtrn.add(annot);
			i = pos + 1;
		}
		if(includeRC) {
			String rc = Sequence.reverseSequence(upper);
			int j = 0;
			while(j < seqString.length()) {
				int pos = seqString.indexOf(rc, j);
				if(pos == -1) {
					break;
				}
				if(pos + s.length() > seqString.length()) {
					// Won't be fully contained
					break;
				}
				// Create annotation
				Annotation annot = new BasicAnnotation(getId(), windowStart + pos, windowStart + pos + s.length(), Strand.POSITIVE);
				rtrn.add(annot);
				j = pos + 1;
			}
		}
		return rtrn;
	}
	
	public Sequence getSubSequence(SequenceRegion region, int extension) {
		String subSeq=sequenceBases.substring(Math.max(region.getStart()-extension, 0), Math.min(region.getEnd()+extension, sequenceBases.length()));
		Sequence seq=new Sequence(region.getName());
		seq.setSequenceBases(subSeq);
		return seq;
	}

	/**
	 * Get subsequence with extension
	 * @param name Name of new sequence to return
	 * @param start Start position of subsequence
	 * @param end Position after last position to include
	 * @param extension Extension
	 * @return The subsequence
	 */
	public Sequence getSubSequence(String name, int start, int end, int extension) {
		String subSeq=sequenceBases.substring(Math.max(start-extension, 0), Math.min(end+extension, sequenceBases.length()));
		Sequence seq=new Sequence(name);
		seq.setSequenceBases(subSeq);
		return seq;
	}
	
	/**
	 * Get subsequence
	 * @param name Name of new sequence to return
	 * @param start Start position of subsequence
	 * @param end Position after last position to include
	 * @return The subsequence
	 */
	public Sequence getSubSequence(String name, int start, int end) {
		return getSubSequence(name, start, end,0);
	}

	/**
	 * Get the spliced transcribed sequence of an annotation
	 * Bases are reported in 5' to 3' direction
	 * @param annot The annotation
	 * @return Sequence with same name as annotation containing the transcribed sequence
	 */
	public Sequence getSubsequence(Gene annot) {
		List<? extends Annotation> blocks = annot.getBlocks();
		Sequence seq = new Sequence(annot.getName());
		for(Annotation block : blocks) {
			Sequence blockSequence = getSubSequence("", block.getStart(), block.getEnd());
			seq.append(blockSequence.getSequenceBases());
		}
		if(annot.getOrientation().equals(Strand.NEGATIVE)) {
			seq.reverse();
		}
		return seq;
	}
	
	/**
	 * Get the spliced transcribed sequence of an annotation
	 * Bases are reported in 5' to 3' direction
	 * @param annot The annotation
	 * @return Sequence with same name as annotation containing the transcribed sequence
	 */
	public Sequence getSubsequence(Annotation annot) {
		if(annot.getOrientation().equals(Strand.UNKNOWN)) {
			throw new IllegalArgumentException("Strand must be known");
		}
		List<? extends Annotation> blocks = annot.getBlocks();
		Sequence seq = new Sequence(annot.getName());
		for(Annotation block : blocks) {
			Sequence blockSequence = getSubSequence("", block.getStart(), block.getEnd());
			seq.append(blockSequence.getSequenceBases());
		}
		if(annot.getOrientation().equals(Strand.NEGATIVE)) {
			seq.reverse();
		}
		return seq;
	}
	
	public void setCapacity(int size) {
		sequenceBases = new StringBuilder(size);
	}

	public Sequence scramble() {
		ArrayList<Character> list=new ArrayList();
		
		char[] bases=this.getSequenceBases().toCharArray();
		for(int i=0; i<bases.length; i++){list.add(bases[i]);}
		
		char[] scrambled=new char[bases.length];
		for(int i=0; i<scrambled.length; i++){
			int index=new Double(Math.random()*list.size()).intValue();
			scrambled[i]=list.remove(index);
		}
		
		String scrSeq=new String(scrambled);
		Sequence rtrn=new Sequence(this.getId()+"_scrambled");
		rtrn.setSequenceBases(scrSeq);
		
		return rtrn;
	}

	public static String generateRandomSequence(int oligoSize) {
		char[] bases=new char[oligoSize];
		for(int i=0; i<bases.length; i++){
			//pick a number from 1-4
			double rand=Math.random();
			if(rand<.25){bases[i]='A';}
			else if(rand<.5){bases[i]='C';}
			else if(rand<.75){bases[i]='G';}
			else {bases[i]='T';}
		}
		return new String(bases);
	}

	public Sequence getAntisense() {
		String anti=reverseSequence(getSequenceBases());
		Sequence antisense=new Sequence(this.getId()+"_antisense");
		antisense.setSequenceBases(anti);
		return antisense;
	}

	public static String get3Prime(String seq, int numBases) {
		char[] bases=seq.toCharArray();
		
		
		int start=Math.max(0, bases.length-numBases);
		char[] rtrn=new char[bases.length-start];
		
		for(int i=start; i<bases.length; i++){
			rtrn[i-start]=bases[i];
		}
		return new String(rtrn);
	}

	public static Collection<String> generateAll5Mers() {
		Collection<String> kmers=new ArrayList<String>();
		
		char[] bases={'A','C','G','T'};
		
		for(int i1=0; i1<bases.length; i1++){
			for(int i2=0; i2<bases.length; i2++){
				for(int i3=0; i3<bases.length; i3++){
					for(int i4=0; i4<bases.length; i4++){
						for(int i5=0; i5<bases.length; i5++){
							char[] random={bases[i1], bases[i2], bases[i3], bases[i4], bases[i5]};
							String seq=new String(random);
							kmers.add(seq);
						}
					}
				}
			}
		}
		return kmers;
	}
	
	/**
	 * Get all kmers composed of the bases A, C, G, T
	 * @param k K
	 * @return All 4^k kmers
	 */
	public static Collection<String> generateAllKmers(int k) {
		if(k < 1) {
			throw new IllegalArgumentException("K must be >= 1");
		}
		if(k == 1) {
			Collection<String> rtrn = new ArrayList<String>();
			rtrn.add("A");
			rtrn.add("C");
			rtrn.add("G");
			rtrn.add("T");
			return rtrn;
		}
		Collection<String> prefixes = generateAllKmers(k-1);
		Collection<String> rtrn = new TreeSet<String>();
		for(String prefix : prefixes) {
			rtrn.add(prefix + "A");
			rtrn.add(prefix + "C");
			rtrn.add(prefix + "G");
			rtrn.add(prefix + "T");
		}
		return rtrn;
	}

	public static String complement(String seq) {
		char[] bases=seq.toCharArray();
		char[] rtrn=new char[bases.length];
		
		for(int i=0; i<bases.length; i++){
			if(bases[i]=='A'){rtrn[bases.length - 1 - i]='T';}
			else if(bases[i]=='C'){rtrn[bases.length - 1 - i]='G';}
			else if(bases[i]=='G'){rtrn[bases.length - 1 - i]='C';}
			else if(bases[i]=='T'){rtrn[bases.length - 1 - i]='A';}
			else if(bases[i]=='a'){rtrn[bases.length - 1 - i]='t';}
			else if(bases[i]=='c'){rtrn[bases.length - 1 - i]='g';}
			else if(bases[i]=='g'){rtrn[bases.length - 1 - i]='c';}
			else if(bases[i]=='t'){rtrn[bases.length - 1 - i]='a';}
			else{rtrn[i]='N';}
		}
		
		return new String(rtrn);
	}
	
	/**
	 * Get the Refseq gene name which appears between the last pair of parentheses in fasta header
	 * @return gene name
	 */
	public String getRefseqGeneName() {
		// Check that sequence ID is in RefSeq format
		if(this.id.indexOf(")") - this.id.indexOf("(") <= 0) {
			throw new InvalidParameterException("Sequence ID " + this.id + " is not a valid RefSeq name");
		};
		// Parse sequence ID and extract name between parentheses
		// Last pair of parentheses contains the gene name
		StringParser p = new StringParser();
		p.parse(this.id,"\\(");
		// Parse the last field again
		int fields = p.getFieldCount();
		p.parse(p.asString(fields - 1),"\\)");
		return p.asString(0);
	}
	
	/**
	 * Get the species name that appears right after the transcript ID or the PREDICTED tag
	 * @return
	 */
	public String getRefseqSpeciesName() {
		// Check that sequence ID is in RefSeq format
		if(this.id.indexOf("|") < 0) {
			throw new InvalidParameterException("Sequence ID " + this.id + " is not a valid RefSeq name");
		};		
		// Parse sequence ID and extract two words after last "|"
		StringParser p = new StringParser();
		p.parse(this.id);
		if(p.asString(1).contentEquals("PREDICTED:")) return p.asString(2) + " " + p.asString(3);
		return p.asString(1) + " " + p.asString(2);
		
	}

	public String getRefseqRnaClassName() {
		if(this.id.contains(", mRNA")) return "mRNA";
		if(this.id.contains(", non-coding RNA")) return "lincRNA";
		if(this.id.contains(", small nuclear RNA")) return "snRNA";
		if(this.id.contains(", small nucleolar RNA")) return "snoRNA";
		if(this.id.contains(", microRNA")) return "microRNA";
		if(this.id.contains(", guide RNA")) return "guideRNA";
		if(this.id.contains(", anstisense RNA")) return "antisenseRNA";
		if(this.id.contains(", RNase P RNA")) return "RNasePRNA";
		if(this.id.contains(", telomerase")) return "telomeraseRNA";
		if(this.id.contains(", small cytoplasmic RNA")) return "smallCytoplasmicRNA";
		if(this.id.contains(", RNase MRP RNA")) return "RNaseMRPRNA";
		if(this.id.contains(", ribosomal RNA")) return "rRNA";
		if(this.id.contains(", ribosomal precursor")) return "rRNAprecursor"; // not a real RefSeq label
		throw new IllegalArgumentException("RefSeq RNA class is not recognized for sequence ID " + this.id);
	}

	/**
	 * Gvien a kmer size, it builds a map of all kmers found
	 * @param kmer
	 */
	public Map<String, KmerInfo> buildKmerMap(int kmer) {
		Map<String, KmerInfo> kmerData = new HashMap<String, KmerInfo>();
		String bases = getSequenceBases();
		logger.debug("Got sequence bases");
		float [] composition = computeSequenceComposition(bases);
		int pos = 0;
		while (pos < getLength() - kmer + 1) {
			String kmerStr = bases.substring(pos,pos + kmer).toUpperCase();
			KmerInfo ki = kmerData.get(kmerStr);
			if(ki == null) {
				//logger.debug("New kmer " + kmerStr + " at " + pos);
				ki = new KmerInfo(kmerStr);
				ki.setSequenceComposition(composition);
				kmerData.put(kmerStr, ki);
			} else {
				//logger.debug("Already seen " + kmerStr + " has " + ki.instances() + " positions " + ki.positions + " adding " + pos);
			}
			ki.addPosition(pos);
			pos++;
		}
		//logger.debug("kmerData size " + kmerData.size());
		return kmerData;
		
	}
	
	public static class KmerInfo {
		String kmerString;
		float []  parentSequenceComposition;
		List<Integer>  positions = new ArrayList<Integer>();
		List<Integer>  genomicPositions = new ArrayList<Integer>();
		public KmerInfo(String kmerStr) {
			this.kmerString = kmerStr;
		}
		public void setSequenceComposition(float[] composition) {
			this.parentSequenceComposition = composition;
		}
		public void addPosition(int pos) {
			positions.add(pos);
		}
		public int instances() {
			return positions.size();
		}
		
		public String toString() {
			StringBuilder sb = new StringBuilder(kmerString);
			sb.append("\t").append(String.valueOf(instances())).append("\t");
			Iterator<Integer> posIt = positions.iterator();
			while(posIt.hasNext()) {
				sb.append(String.valueOf(posIt.next()));
				if(posIt.hasNext()) {
					sb.append(",");
				}
			}
			
			
			sb.append("\t");
			posIt = genomicPositions.iterator();
			while(posIt.hasNext()) {
				sb.append(String.valueOf(posIt.next()));
				if(posIt.hasNext()) {
					sb.append(",");
				}
			}
			return sb.toString();
		}
		public void mapPositionsToGenome(Gene g) {
			genomicPositions = new ArrayList<Integer>();
			for(int p : positions)  {
				genomicPositions.add(g.transcriptToGenomicPosition(p));
			}
			
		}
	}
	
}


