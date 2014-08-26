package broad.core.multiplealignment;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Stack;

import nextgen.core.annotation.Gene;

import org.apache.log4j.Logger;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.parsers.nhx.NHXParser;

import Jama.Matrix;
import broad.core.annotation.AnnotationReader;
import broad.core.annotation.AnnotationReaderFactory;
import broad.core.annotation.BasicGenomicAnnotation;
import broad.core.annotation.GenomicAnnotation;
import broad.core.annotation.LightweightGenomicAnnotation;
import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.core.error.ParseException;
import broad.core.sequence.SequenceRegion;
import broad.core.siphy.ConservationUtils;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.GeneWithIsoforms;

public class MAFAlignment extends MultipleAlignment {
	static Logger logger = Logger.getLogger(MultipleAlignment.class.getName());
	private IntervalTree<MAFMultipleAlignmentBlock> alignmentBlockTree;
	private List<String> sequenceIds;

	static final String HG17_ALIGN_DIR="/seq/genome/ucsc/multiz17way";
	static final String HG17_ALIGN_REF_PREFIX="hg17.";
	static final int indexPositionJump = 100000;
	
	public static String USAGE = "Usage: MAFAlignment TASK=<task_num> <task_args>\n" +
	"\tTasks:" +
	"\n\t\t1. Extract alignments  -in <input file> -out <output file, default writes to standard out>  -start <chromosome start> -end <chromosome end> [-compress <if this flag is present, the alignment output will contain no reference alignment gaps> -outformat <Specify FASTA, PHYLIP or PHYLIPSEQ if a format other than MAF is decired for the output> -noPadding (include this flag if no padding betweein non contiguous MAF blocks is desired, not that any reference to alignment positions will be lost) ]\n" +
	"n\t\t2. Extract regions -in <alignment file in MAF format> -out <output file> -annotations <Anntation file (BED format by default)> [-format <[BED], SIMPLE,GFF>] -seqsToLoad <List of sequences to load or all if non specified> -compress <if the result should be devoid of reference gaps> -outformat <FASTA, PHYLIP, SEQPHYLIP, if none is specifed the default is MAF>]\n" +
	"n\t\t3. Create index file -in <alignment file in MAF format> -out <Index file name or standard output>\n"+
	"n\t\tExtractBed. Extract segements from a bed from different chrs : similar to extract region but will extract seq of each BED line and will generate an outfile for each  exon specified in the BED   -in <Directory of chr alignment file in MAF format ; dir/chrX.maf> -out <output prefix (will generate multiple out files, one for each record in the BED)> -annotations <Anntation file (BED format by default)> [-format <[BED], SIMPLE,GFF>] -seqsToLoad <List of sequences to load or all if non specified> -compress <if the result should be devoid of reference gaps> -outformat <FASTA, PHYLIP, SEQPHYLIP, if none is specifed the default is MAF>]  -fullBed <flag if you want each exon to be reported>\n" ;
	
	
	MAFHeader header;
	//Map<Integer, Long> index;
	IntervalTree<Long> index;
	private String referenceChromosome;
	
	public MAFAlignment() {
		super();
		header = new MAFHeader();
		alignmentBlockTree = new IntervalTree<MAFMultipleAlignmentBlock>();
		sequenceIds = new ArrayList<String>();
	}
	
	public MAFAlignment(String indexFile) throws IOException, ParseException {
		this();
		loadIndex(indexFile);
		//System.out.println("Loaded index file " + index.keySet().size());
	}
	
	public MAFHeader getHeader() { return header;}
	
	// MultipleAlignment overloading.
	public void encode() {
		//Iterator<MAFMultipleAlignmentBlock> alnBlockIt = alignmentBlocks.iterator();
		Iterator<Node<MAFMultipleAlignmentBlock>> alnBlockIt = alignmentBlockTree.iterator();
		while(alnBlockIt.hasNext()) {
			//MultipleAlignment block = alnBlockIt.next();
			Node<MAFMultipleAlignmentBlock> blockNode = alnBlockIt.next();
			blockNode.getValue().encode();
		}
		
	}
	
	public void encodeAsMatrix() {
		//Iterator<MAFMultipleAlignmentBlock> alnBlockIt = alignmentBlocks.iterator();
		Iterator<Node<MAFMultipleAlignmentBlock>> alnBlockIt = alignmentBlockTree.iterator();
		while(alnBlockIt.hasNext()) {
			//MultipleAlignment block = alnBlockIt.next();
			Node<MAFMultipleAlignmentBlock> blockNode = alnBlockIt.next();
			blockNode.getValue().encodeAsMatrix();
		}		
	}
	
	public void reverse() {
		//Iterator<MAFMultipleAlignmentBlock> alnBlockIt = alignmentBlocks.iterator();
		Iterator<Node<MAFMultipleAlignmentBlock>> alnBlockIt = alignmentBlockTree.iterator();
		while(alnBlockIt.hasNext()) {
			//MultipleAlignment block = alnBlockIt.next();
			Node<MAFMultipleAlignmentBlock> blockNode = alnBlockIt.next();
			blockNode.getValue().reverse();
		}
		/*
		Stack<MAFMultipleAlignmentBlock> reversedBlocks = new Stack<MAFMultipleAlignmentBlock>();
		while(alignmentBlocks.size() > 0) {
			MAFMultipleAlignmentBlock block = alignmentBlocks.pop();
			block.reverse();
			reversedBlocks.push(block);
		}			
		alignmentBlocks = reversedBlocks;
		*/
	}
	
	public int length() {
		int length = 0;
		//if(!alignmentBlocks.isEmpty()) {
		if(alignmentBlockTree.size() > 0) {
			MAFMultipleAlignmentBlock first = alignmentBlockTree.min().getValue();//alignmentBlocks.get(0);
			MAFMultipleAlignmentBlock last  = alignmentBlockTree.max().getValue();//alignmentBlocks.get(alignmentBlocks.size() - 1);
			
			length = last.getReferenceEnd()  - first.getReferenceStart();
		}
		
		return length;
	}
	
	/**
	 * Assumes that encode has been called and each aligned sequence is encoded
	 * @param int - alignment column
	 */
	public Map<String, Short> getColumn(int i) {
		Map<String, Short> col = getGapColumn();
		Iterator<Node<MAFMultipleAlignmentBlock>> overlappingNodeIt = alignmentBlockTree.overlappers(i, i+1);

		if(overlappingNodeIt.hasNext()) {
			MAFMultipleAlignmentBlock containingBlock = overlappingNodeIt.next().getValue();
			Map<String, Short> closestCol = containingBlock.getColumn(i);
			Iterator<String> closestAlignSeqIdIt = closestCol.keySet().iterator();
			while(closestAlignSeqIdIt.hasNext()) {
				String seqId = closestAlignSeqIdIt.next();
				col.put(seqId, closestCol.get(seqId));
			}
		}

		return col;
	}
	
	public Map<String, Matrix> getColumnsAsVector(int start, int number) {
		LinkedHashMap<String, Matrix> cols = new LinkedHashMap<String, Matrix>(getAlignedSequenceIds().size());
		Iterator<String> seqIdIt = getAlignedSequenceIds().iterator();

		while(seqIdIt.hasNext()) {
			cols.put(seqIdIt.next(), new Matrix(UNGAPPED_ALPHABET_SIZE, number));
		}
		
		for(int i = start; i < start + number; i++) {
			Map<String, Short> col = getColumn(i);
			seqIdIt = col.keySet().iterator();
			while(seqIdIt.hasNext()) {
				String seqId = seqIdIt.next();
				short base = col.get(seqId);
				Matrix seqMatrix = cols.get(seqId);
				if(base != GAP_CODE) {
					seqMatrix.set(base, i - start,1);
				}
			}
		}
		return cols;
	}

	public void addShortEncodedColumn(Map<String, Short> col) {
		throw new RuntimeException("Not yet implemented, to edit an alignment use a base MultipleAlignment MAF are ReadOnly");
	}
	
	public void addShortEncodedRegion(Map<String, short[]> region) {
		throw new RuntimeException("Not yet implemented, to edit an alignment use a base MultipleAlignment MAF are ReadOnly");
	}
	
	public void addSequence(AlignedSequence seq) {
		throw new RuntimeException("Not yet implemented, to edit an alignment use a base MultipleAlignment MAF are ReadOnly");
	}

	public void addSequences(List<AlignedSequence> sequences) {
		throw new RuntimeException("Not yet implemented, to edit an alignment use a base MultipleAlignment MAF are ReadOnly");		
	}
	
	public boolean isEmpty() { return alignmentBlockTree.isEmpty();}
	
	public void write(BufferedWriter bw) throws IOException {
		if(getIOHelper() != null) {
			super.write(bw);
		} else {
			header.addVariableValuePair("ref", getReferenceId());
			header.write(bw);
			Iterator<Node<MAFMultipleAlignmentBlock>> blockIt = alignmentBlockTree.iterator();
			while(blockIt.hasNext()) {
				blockIt.next().getValue().write(bw);
			}
		}
		
	}
	
	public MAFAlignment getSubAlignment(int refStart, int refEnd, boolean reverse) {
		MAFAlignment subAln = new MAFAlignment();
		subAln.setReferenceId(getReferenceId());
		subAln.sequenceIds = this.sequenceIds;
		Iterator<MAFMultipleAlignmentBlock> overlapperIt = 
			new IntervalTree.ValuesIterator<MAFMultipleAlignmentBlock>(alignmentBlockTree.overlappers(refStart, refEnd));
		LightweightGenomicAnnotation target = new BasicGenomicAnnotation(getReference());
		target.setStart(refStart);
		target.setEnd(refEnd);
		while(overlapperIt.hasNext()) {
			MAFMultipleAlignmentBlock overlappingBlock = overlapperIt.next();
			subAln.addBlock(overlappingBlock.trim(target));
		}
		return subAln;
	}
	
	public AlignedSequence getAlignedSequence(String sequenceId, boolean fillInGapsBetweenBlocks) {
		//System.out.println("Is alignmentBlockTree empty? " + alignmentBlockTree.isEmpty());
		if(!sequenceIds.contains(sequenceId) || alignmentBlockTree.isEmpty() ) {
			//System.out.println("seqids " + sequenceIds + " does not contain <" + sequenceId +"> is contained in sequenceIds " + sequenceIds.contains(sequenceId));
			return null;
		}
		
		AlignedSequence seq = new AlignedSequence(sequenceId);
		seq.setId(sequenceId);
		if(getReferenceId().equals(sequenceId)) {
			//System.out.println("Dealing with ref seq");
			seq.setStart(alignmentBlockTree.min().getValue().getReferenceStart());
			seq.setChromosome(alignmentBlockTree.min().getValue().getAlignedSequence(sequenceId).getChromosome());
			seq.setEnd(alignmentBlockTree.max().getValue().getReferenceEnd());
		}
		Iterator<MAFMultipleAlignmentBlock> it = new IntervalTree.ValuesIterator<MAFMultipleAlignmentBlock>(alignmentBlockTree.iterator());
		MAFMultipleAlignmentBlock lastBlock = null;
		while(it.hasNext()) {
			MAFMultipleAlignmentBlock block = it.next();
			AlignedSequence segment = block.getAlignedSequence(sequenceId);
			if(lastBlock!= null && lastBlock.getReferenceEnd() < block.getReferenceStart() && fillInGapsBetweenBlocks) {
				//System.out.println("last block " + (lastBlock != null ? lastBlock.getReferenceStart()+"-"+lastBlock.getReferenceEnd():" not yet ") + " new block " + block.getReferenceStart()+"-"+block.getReferenceEnd() + ", appending " + (block.getReferenceStart() - lastBlock.getReferenceEnd()) + " Ns");
				for(int i = 0; i < block.getReferenceStart() - lastBlock.getReferenceEnd(); i++) {
					//System.out.println("\t\tappending gaps N");
					seq.append('N');
				}
			}
			seq.append(segment.getSequenceBases());
			/*	
			if(segment != null) {
				seq.appendToSequence(segment.getSequenceBases());
			} else {
				for(int i = 0; i < blockNode.getValue().length(); i++) {
					seq.appendToSequence('-');
				}
			}
			*/
			lastBlock = block;
		}
		return seq;
	}
	
	public AlignedSequence getAlignedSequence(String sequenceId) { 
		return getAlignedSequence(sequenceId, true);
	}
	
	public List<AlignedSequence> getAlignedSequences(boolean fillGapsBetweenBlocks) {
		ArrayList<AlignedSequence> sequences = new ArrayList<AlignedSequence>(sequenceIds.size());
		Iterator<String> idIt = sequenceIds.iterator();
		while(idIt.hasNext()) {
			String id = idIt.next();
			AlignedSequence alignedSeq = getAlignedSequence(id, fillGapsBetweenBlocks);
			sequences.add(alignedSeq);
			
		}
		return sequences;		
	}
	public List<AlignedSequence> getAlignedSequences() {
		return getAlignedSequences(true);
	}
	
	public boolean overlaps(LightweightGenomicAnnotation annotation) {return index.overlappers(annotation.getStart(), annotation.getEnd()).hasNext();}
	
	public boolean contains(LightweightGenomicAnnotation annotation) {
		return index.min().getStart() < annotation.getStart() && index.max().getEnd() > annotation.getEnd();
	}
	
	public List<String> getAlignedSequenceIds() { return sequenceIds;}
	
	public int getReferenceStart() {
		Node<MAFMultipleAlignmentBlock> first = alignmentBlockTree.min();
		return first == null ? -1 : first.getValue().getReferenceStart();
	}
	
	public int getReferenceEnd() {
		Node<MAFMultipleAlignmentBlock> last = alignmentBlockTree.max();
		return last == null ? -1 : last.getValue().getReferenceEnd();
	}
	
	public void compress() {
		Iterator<Node<MAFMultipleAlignmentBlock>> nodeIt =  alignmentBlockTree.iterator();
		while(nodeIt.hasNext()) {
			nodeIt.next().getValue().compress();
		}
	}
	
	public MultipleAlignment toMultipleAlignment(boolean fillGapsBetweenBlocks) {
		MultipleAlignment ma = new MultipleAlignment();
		ma.setReferenceId(getReferenceId());
		ma.setRefGapped(true);
		ma.addSequences(getAlignedSequences(fillGapsBetweenBlocks));
		ma.getReference().setStart(getReferenceStart());
		ma.getReference().setEnd(getReferenceEnd());
		//System.out.println("to multiplealn #aln seqs: " + ma.length() + " ref? " + ma.getReferenceId());
		return ma;		
	}
	
	public MultipleAlignment toMultipleAlignment() {
		return toMultipleAlignment(true);
	}
	
	public MultipleAlignment sampleColumns(int colNum, int numConsecutiveCols) {
		Random r = new Random();
		int treeSize = alignmentBlockTree.size();
		MultipleAlignment sampledAlignment = new MultipleAlignment();
		sampledAlignment.setReferenceId(getReferenceId());
		for(int i = 0; i < colNum - numConsecutiveCols + 1; i = numConsecutiveCols + i) {
			int blockId = r.nextInt(treeSize);
			System.out.println("randomly selected block " + blockId + " out of " + treeSize + " blocks" );
			MAFMultipleAlignmentBlock block = alignmentBlockTree.findByIndex(blockId).getValue();
			MultipleAlignment sampledCols = block.sampleColumns(1,numConsecutiveCols);
			Iterator<String> seqIdIt = getAlignedSequenceIds().iterator();
			while(seqIdIt.hasNext()) {
				String seqId = seqIdIt.next();
				if(sampledCols.getAlignedSequence(seqId) == null) {
					AlignedSequence missingSeq = new AlignedSequence(seqId);
					for(int j = 0; j < numConsecutiveCols; j++) {
						missingSeq.append('-');
					}
					sampledCols.addSequence(missingSeq);
				}
			}
			sampledAlignment.append(sampledCols); //changed without test 5/23/08
		}

		return sampledAlignment;
	}
	
	public void load (String fileName) throws IOException, ParseException {
		load (fileName, 0, SequenceRegion.INF);
	}
	
	public void load(RandomAccessFile handle, List<String> sequencesToLoad) throws IOException, ParseException {
		load(handle, 1, SequenceRegion.INF, sequencesToLoad);
	}
	
	public void load(RandomAccessFile handle, int referenceStart, int referenceEnd) throws IOException, ParseException {
		load(handle, referenceStart, referenceEnd, new ArrayList<String>());
	}
	
	public void load(RandomAccessFile handle, int referenceStart, int referenceEnd, List<String> sequencesToLoad) throws IOException, ParseException {
		long offset = getClosestOffset(referenceStart);
		logger.debug("Starting to read file from " + offset + " from refstart " + referenceStart + " to refend " + referenceEnd);
		alignmentBlockTree = new IntervalTree<MAFMultipleAlignmentBlock>();
		handle.seek(offset);
		boolean okToAdd = true;
		BasicGenomicAnnotation reference = new BasicGenomicAnnotation("reference");
		reference.setStart(referenceStart);
		reference.setEnd(referenceEnd);
		sequenceIds = new ArrayList<String>();
		
		if(sequencesToLoad != null) {
			Iterator<String> seqIt = sequencesToLoad.iterator();
			while(seqIt.hasNext()) {
				sequenceIds.add(seqIt.next());
			}
		}
		
		String line = null;
		String [] lineInfo = null;
		Stack<MAFMultipleAlignmentBlock> alignmentBlockStack = new Stack<MAFMultipleAlignmentBlock>();
		//ArrayList<String> lastAlignmentLines = new ArrayList<String>(); //FOR easier debugging
		while((line = handle.readLine()) != null) {
			//Ignore all comment lines
			if(line.startsWith("#") || line.trim().length() == 0){
				continue;
			}
			
			//System.out.println(line);
			if(line.startsWith("a ")) {
				//System.err.println("New alignment: " + line);
				//First check last alignment to see if it should be kept.
				MAFMultipleAlignmentBlock lastMA = alignmentBlockStack.isEmpty()  ? null : alignmentBlockStack.pop();
				if(lastMA != null) {
					alignmentBlockTree.put(lastMA.getReferenceStart(), lastMA.getReferenceEnd(), lastMA.trim(reference));
				}
				/*
				if(lastMA != null && referenceEnd > 0 && !lastMA.overlaps(reference)) {
					alignments.pop(); // multiple alignment just built does overlap desired region.
					AlignedSequence refAln = lastMA.getAlignedSequence(reference.getContainingSequenceId());
					System.out.println("Reference seq<"+reference.getContainingSequenceId() +"> lastMA aligned seqs : " + lastMA.getAlignedSequenceIds());
					if (refAln.getStart() > referenceEnd) {
						break;
					}
				}
				*/
				MAFMultipleAlignmentBlock ma = new MAFMultipleAlignmentBlock();
				alignmentBlockStack.push(ma);
				lineInfo = line.substring(2).split("\\s");
				ma.setAlignmentInfoFromRawData(lineInfo);
				okToAdd = true;
				//lastAlignmentLines.clear();
			} else if(line.startsWith("s "))  {
				if(alignmentBlockStack.isEmpty() || !okToAdd) {
					continue;
				}
				//System.err.println("\tAlignment aligned seq " + line);
				MAFMultipleAlignmentBlock ma = alignmentBlockStack.peek();
				lineInfo = line.substring(2).split("\\s+");
				AlignedSequence seq = ma.createSequenceFromRawData(lineInfo);
				
				if(getReferenceId() == null) {
					setReferenceId(seq.getId());
					setReferenceChromosome(seq.getChromosome());
				} 
				
				if(ma.getReferenceId() == null) {
					ma.setReferenceId(seq.getId());
				}
				
				if(sequencesToLoad == null || sequencesToLoad.size() == 0 || sequencesToLoad.contains(seq.getId())) {
					//System.err.println("YES: sequence " + seq.getId() + " chr " + seq.getChromosome() + " sequenceToLoad is null? " + (sequencesToLoad == null) + " is empty? " + (sequencesToLoad.size() == 0 ) + " or contains sequence " + sequencesToLoad.contains(seq.getId()));
					ma.addSequence(seq);
				} else {
					//System.err.println("NO: sequence " + seq.getId() + " chr " + seq.getChromosome() + " sequenceToLoad is null? " + (sequencesToLoad == null) + " is empty? " + (sequencesToLoad.size() == 0 ) + " or contains sequence " + sequencesToLoad.contains(seq.getId()));
					continue;
				}
				
				//A bit wasteful, and should fix, but finally check if reference location is OK to add
				if(getReferenceId().equals(seq.getId()) && seq.getEnd() <= referenceStart) {
					alignmentBlockStack.pop();
					okToAdd = false;
				} else if(getReferenceId().equals(seq.getId()) && seq.getStart() >= referenceEnd) {
					alignmentBlockStack.pop();
					break;
				}else {
					//alignmentBlockStack.push(ma);
					if(!getAlignedSequenceIds().contains(seq.getId())) {
						sequenceIds.add(seq.getId());
					}
				}
				//lastAlignmentLines.add(line);
			} else if (line.startsWith("i ")) {
				//We do not handle information lines yet.
				continue;
			}else if (line.startsWith("q ")) {
				//We do not handle quality lines yet.
				continue;
			}else if (line.startsWith("e ") ) {
				//We do not support e lines yet.
				continue;
			} else {
				throw new ParseException("Invalid alignment line <"+ line +">");
			}
			
		}
		//Handle last alignment block
		MAFMultipleAlignmentBlock lastMA = alignmentBlockStack.isEmpty()  ? null : alignmentBlockStack.pop();
		if(lastMA != null && lastMA.getReferenceStart() < referenceEnd && lastMA.getReferenceEnd() > referenceStart) {
			alignmentBlockTree.put(lastMA.getReferenceStart(), lastMA.getReferenceEnd(), lastMA.trim(reference));
			//System.err.println("put last alignment: ");
		}
	}

	public void load(String fileName, int referenceStart, int referenceEnd) 
		throws IOException, ParseException {
		RandomAccessFile raf = new RandomAccessFile(fileName, "r");
		load(raf, referenceStart, referenceEnd);
		raf.close();
	}
	
	public IntervalTree<Long> getIndex() { return index;}
	
	public void writeIndex(String indexFileName) throws IOException {
		Iterator<Node<Long>> idxEntryIt = index.iterator();
		BufferedWriter bw = new BufferedWriter(new FileWriter(indexFileName));
		while(idxEntryIt.hasNext()) {
			Node<Long> entry = idxEntryIt.next();
			bw.write(String.valueOf(entry.getStart()));
			bw.write("\t" + String.valueOf(entry.getEnd() - entry.getStart()));
			bw.write("\t" + String.valueOf(entry.getValue()));
			bw.newLine();
		}
		bw.close();
	}
	
	public void createIndex(String alignmentFile) throws IOException {
		RandomAccessFile raf = new RandomAccessFile(alignmentFile,"r");
		index = new IntervalTree<Long>();//LinkedHashMap<Integer, Long>();
		String line;
		String [] lineInfo= null;
		long lastOffset = 0;
		try {
			boolean readNext = false;
			while((line = raf.readLine()) != null) {
				//Ignore all comment lines
				if(line.startsWith("#") || line.trim().length() == 0){
					continue;
				}
				
				if(line.startsWith("a ")) {
					readNext = true;
					lastOffset = raf.getFilePointer() - line.getBytes().length - 1;
				} else if(line.startsWith("s ")) {
					if(readNext) {
						lineInfo = line.split("\\s+");
						int start = Integer.parseInt(lineInfo[2]);
						int end   = Integer.parseInt(lineInfo[3]) + start;
						index.put(start, end, lastOffset);
					}
					readNext = false;
				}else if (line.startsWith("i ")) {
					//We do not handle information lines yet.
					continue;
				}else if (line.startsWith("q ")) {
					//We do not handle quality lines yet.
					continue;
				}else {
					readNext = false;
				}
				
			}
		} finally {
			try {
				raf.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}		
		
	}
	
	public void setBlocks(List <MAFMultipleAlignmentBlock> blocks) {
		alignmentBlockTree = new IntervalTree<MAFMultipleAlignmentBlock>();
		Iterator<MAFMultipleAlignmentBlock> blockIt = blocks.iterator();
		while(blockIt.hasNext()) {
			addBlock(blockIt.next());
		}

	}
	
	public void addBlock(MAFMultipleAlignmentBlock block) throws IllegalArgumentException{
		Iterator<Node<MAFMultipleAlignmentBlock>> overlappers = 
			alignmentBlockTree.overlappers(block.getReferenceStart(), block.getReferenceEnd());
		if(overlappers != null && overlappers.hasNext() ) {
			throw new  IllegalArgumentException("A block in the alignment to append already overlaps exisiting alignment block. Can't append");
		}
		Iterator<String> alnSeqIdIt = block.getAlignedSequenceIds().iterator();
		while(alnSeqIdIt.hasNext()) {
			String seqId = alnSeqIdIt.next();
			if(!sequenceIds.contains(seqId)) {
				sequenceIds.add(seqId);
			}
		}
		
		alignmentBlockTree.put(block.getReferenceStart(), block.getReferenceEnd(), block);
	}
	
	public void addBlocks(List<MAFMultipleAlignmentBlock> blocks) throws IllegalArgumentException {
		addBlocks(blocks.iterator());
	}
	
	public void addBlocks(Iterator<MAFMultipleAlignmentBlock> blockIt) throws IllegalArgumentException{
		while(blockIt.hasNext()) {
			addBlock(blockIt.next());
		}
	}	
	
	public void append(MAFAlignment mafAln) throws IllegalArgumentException{
		addBlocks(mafAln.getBlockIterator());
	}
	
	public void clear () {
		alignmentBlockTree = new IntervalTree<MAFMultipleAlignmentBlock>();
	}

	long getClosestOffset(int position) {
		Node<Long> node = index.max(position, position + 1);
		return node == null ? 0 : node.getValue();
	}
	
	public Iterator<MAFMultipleAlignmentBlock> getBlockIterator() {
		final Iterator<Node<MAFMultipleAlignmentBlock>> nodeIt = alignmentBlockTree.iterator();
		return new Iterator<MAFMultipleAlignmentBlock>() {

			public boolean hasNext() {
				return nodeIt.hasNext();
			}

			public MAFMultipleAlignmentBlock next() {
				return nodeIt.next().getValue();
			}

			public void remove() {
				nodeIt.remove();
			}
			
		};
		

	}
	
	private void setReferenceChromosome(String chromosome) {
		this.referenceChromosome = chromosome;
		
	}
	
	private String getReferenceChromosome() {return this.referenceChromosome;	}

	private void padN(int start, int end) {
		if(getReference() == null) {
			return;
		}
		
		int startPadLength = getReferenceStart() - start;
		int endPadLength   = end - getReferenceEnd();
		if(startPadLength > 0) {
			MAFMultipleAlignmentBlock pad = new MAFMultipleAlignmentBlock();
			pad.setReferenceId(getReferenceId());
			AlignedSequence ref = new AlignedSequence(getReference().getContainingSequenceId());
			ref.setStart(start);
			ref.setEnd(getReferenceStart());
			ref.setId(getReferenceId());
			pad.addAlignedSequence(ref);
			

			ref.setSequenceBases(buildPad(startPadLength));
			addBlock(pad);
		}
		
		if(endPadLength > 0) {
			MAFMultipleAlignmentBlock pad = new MAFMultipleAlignmentBlock();
			pad.setReferenceId(getReferenceId());
			AlignedSequence ref = new AlignedSequence(getReference().getContainingSequenceId());
			ref.setStart(getReferenceEnd() + 1);
			ref.setEnd(end);
			ref.setId(getReferenceId());
			pad.addAlignedSequence(ref);
			

			ref.setSequenceBases(buildPad(endPadLength));
			addBlock(pad);
		}

	}
	
	private String buildPad(int size) {
		StringBuilder padSeq = new StringBuilder(size );
		for(int i = 0; i <  size; i++) {
			padSeq.append("N");
		}
		
		return padSeq.toString();
	}
	public static class MAFMultipleAlignmentBlock extends MultipleAlignment {
		int pass;
		
		public MAFMultipleAlignmentBlock() {
			super();
			setRefGapped(true);
		}
		
		public MAFMultipleAlignmentBlock(MultipleAlignment ma) {
			this();
			setReferenceId(ma.getReferenceId());
			Iterator<AlignedSequence> it = ma.getAlignedSequences().iterator();
			while(it.hasNext()) {
				addAlignedSequence(it.next());
			}
			
		}

		public int getPass() { return pass;}

		
		public AlignedSequence createSequenceFromRawData(String[] data) {		
			String [] seqNameInfo = data[0].split("\\.");
			AlignedSequence aln = new AlignedSequence(seqNameInfo[0].intern());
			if(seqNameInfo.length > 1) {
				aln.setChromosome(seqNameInfo[1]);
			}
			aln.setId(seqNameInfo[0]);
			aln.setName(seqNameInfo[0]);
			aln.setRegionStart(Integer.parseInt(data[1]));
			aln.setRegionEnd(aln.getRegionStart() + Integer.parseInt(data[2]));
			aln.setStrand(data[3]);
			aln.setTotalLength(Integer.parseInt(data[4]));
			aln.setSequenceBases(data[5]);
			return aln;			
		}
		
		public AlignedSequence addSequenceFromRawData(String[] data) {
			AlignedSequence aln = createSequenceFromRawData(data);
			addSequence(aln);
			return aln;
		}
		
		public void addAlignedSequence(AlignedSequence aln) {
			addSequence(aln);
		}
		
		public AlignedSequence getAlignedSequence(String sequenceId) {
			AlignedSequence seq = super.getAlignedSequence(sequenceId);
			if(seq == null) {
				seq = new AlignedSequence(sequenceId);
				int length = length();
				for(int i = 0; i < length; i++) {
					seq.append('-');
				}
			}
			
			return seq;
		}
		
		public void setAlignmentInfoFromRawData(String[] data) {
			for(int i = 0; i < data.length; i++) {
				String[] nameValPair = data[i].split("=");
				if(nameValPair[0].equalsIgnoreCase("pass")) {
					setPass(Integer.parseInt(nameValPair[1]));
				} else if(nameValPair[0].equalsIgnoreCase("score")) {
					setScore(Float.parseFloat(nameValPair[1]));
				} else {
					System.err.println("Unsuported alignment attribute: " + nameValPair[0] + " .... ignoring");
				}
			}
		}
		
		public void setPass(int pass) { this.pass = pass;}
		
		public MAFMultipleAlignmentBlock trim(LightweightGenomicAnnotation target) {
			//System.out.println("target: " + target.getLocationString() + " reference: " + getReference().getLocationString());
			if(target.getStart() <= getReferenceStart() && target.getEnd() >= getReferenceEnd()) {
				//System.out.println("No trim is nencesary");
				return this;
			}
			
			return new MAFMultipleAlignmentBlock(getSubAlignment(target.getStart(), target.getEnd(), false));
			
		}
		
		public void write(BufferedWriter bw) throws IOException {
			if(getIOHelper() != null) {
				super.write(bw);
			} else {
				int maxSequenceLength = 0;
				int maxSequenceIdLength = 0;
				int maxAlignmentStartLength = 0;
	
				Iterator<AlignedSequence> seqIt = getAlignedSequences().iterator();
				while(seqIt.hasNext()) {
					AlignedSequence seq = seqIt.next();
					if(maxSequenceLength < seq.getLength()) {
						maxSequenceLength = seq.getLength();
					}
					String seqChr = seq.getChromosome() != null ? "."+seq.getChromosomeString() : "";
					if(seq.getContainingSequenceId().length() + seqChr.length()> maxSequenceIdLength) {
						maxSequenceIdLength = seq.getContainingSequenceId().length() + seqChr.length();
					}
					
					if(seq.getStart() > maxAlignmentStartLength) {
						maxAlignmentStartLength = seq.getStart();
					}
				}
				int seqLengthFieldSize = (int) Math.floor(Math.log10(maxSequenceLength)) + 1;
				int seqStartFieldSize = (int) Math.floor(Math.log10(maxAlignmentStartLength)) + 1;
				
				bw.write("a");
				bw.write(" score=");
				bw.write(String.valueOf(getScore()));
				if (getPass() > 0) {
					bw.write(" pass=");
					bw.write(String.valueOf(getPass()));
				}
				bw.newLine();
				
				seqIt = getAlignedSequences().iterator();
				AlignedSequence seq = null;
				while(seqIt.hasNext()) {
					seq = seqIt.next();
					String seqChr = seq.getChromosome() != null ? "."+seq.getChromosomeString() : "";
					bw.write("s ");
					CLUtil.writeLeftJustifiedField(bw,  seq.getContainingSequenceId() + seqChr, maxSequenceIdLength);
					bw.write(" ");
					CLUtil.writeRightJustifiedField(bw, String.valueOf(seq.getRegionStart()), seqStartFieldSize);
					bw.write(" ");
					CLUtil.writeRightJustifiedField(bw, String.valueOf(seq.getUngappedLength()), seqStartFieldSize);
					bw.write(" ");
					bw.write(seq.inReversedOrientation() ? "- " : "+ ");
					CLUtil.writeRightJustifiedField(bw, String.valueOf(seq.getTotalLength()), seqLengthFieldSize);
					bw.write(" ");
					bw.write(seq.getSequenceBases());
					bw.newLine();
				}
				bw.newLine();
			}
		}
	}
	
	public static class MAFHeader {
		MAFAlignment alignment;
		String version;
		String scoring;
		String program;
		
		String runParameters;
		
		Hashtable<String, String> otherVariables = new Hashtable<String, String>();
		
		private static final long serialVersionUID = 239451013421586L;

		protected void setVariablesFromRawData(String [] data) {
			for(int i = 0; i < data.length; i++) {
				String [] variableValuePair = data[i].split("=");
				if(variableValuePair[0].equalsIgnoreCase("version")) {
					version = variableValuePair[1];
				} else if (variableValuePair[0].equalsIgnoreCase("scoring")) {
					scoring = variableValuePair[1];
				} else if (variableValuePair[0].equalsIgnoreCase("program")) {
					program = variableValuePair[1];
				} else {
					otherVariables.put(variableValuePair[0].toLowerCase(), variableValuePair[1]);
				}
				
			}
		}
		
		public String getRunParameters() {
			return runParameters;
		}

		public void addVariableValuePair(String variable, String value) {
			otherVariables.put(variable, value);
		}
		public void write(BufferedWriter bw) throws IOException {
			bw.write("##maf");
			if(version != null && version.length() > 0) {
				bw.write(" version=");
			    bw.write(version);
			}
			if(scoring != null && scoring.length() > 0) {
				bw.write(" scoring=");
				bw.write(scoring);
			}
			if(program != null && program.length() > 0) {
				bw.write(" program=");
				bw.write(program);
			}
			
			Iterator<String> varIt = otherVariables.keySet().iterator();
			while(varIt.hasNext()) {
				String var = varIt.next();
				bw.write(" ");
				bw.write(var);
				bw.write("=");
				bw.write(otherVariables.get(var));
			}
			bw.newLine();
			
			if(runParameters != null) {
				bw.write("# ");
				bw.write(runParameters);
			}

			bw.newLine();
			bw.newLine();
		}

		public void setRunParameters(String runParameters) {
			this.runParameters = runParameters;
		}



		public void setProgram(String program) {
			this.program = program;
		}

		public void setScoring(String scoring) {
			this.scoring = scoring;
		}

		public void setVersion(String version) {
			this.version = version;
		}

		public String getProgram() {
			return program;
		}

		public String getScoring() {
			return scoring;
		}

		public String getVersion() {
			return version;
		}
		
	}
	
	
	
	
	
	protected void setAlignedSequences(List<String> seqIds) {
		sequenceIds = seqIds;
	}

	public void loadIndex(String idxFile) throws IOException {
		index = new IntervalTree<Long>();//new LinkedHashMap<Integer, Long>();;
		BufferedReader br = new BufferedReader(new FileReader(idxFile));
		//System.err.println("Loading index " + idxFile);
		String line = null;
		int l = 0;
		while((line = br.readLine()) != null) {
			String [] info = line.split("\t");
			int start = Integer.parseInt(info[0]);
			int end   = Integer.parseInt(info[1]) + start;
			long offset = Long.parseLong(info[2]);
			//System.out.println("Read line "+ l++);
			index.put(start, end, offset);

		}
		
	}
	
	private Map<String, Short> getGapColumn() {
		HashMap<String, Short> gapCol = new HashMap<String, Short>(sequenceIds.size());
		Iterator<String> sequenceIdIt = sequenceIds.iterator();
		while(sequenceIdIt.hasNext()) {
			gapCol.put(sequenceIdIt.next(), GAP_CODE);
		}
		
		return gapCol;
	}
	
	//Given a bed file and the directory of MAF chr alignments, will extract the specific 
	//aligned region from the MSA and output it to a separate file (one file per record).
	//will also output a files listing all output files
	//spFilter - only report alignment with that species
	//updated to extract by exons .. 
	private static void extractBedRegions(String inDir, String outDir,
			String annotations, boolean noPadding, String outformat,
			boolean compress , List<String> seqsToLoad,String spFilter, boolean fullBed) throws IOException, ParseException {
		
		
		BufferedWriter bw = new BufferedWriter(new FileWriter(outDir+"/MafOutFiles.txt"));
		
	
		BEDFileParser bed = new BEDFileParser (annotations);
		for (Gene g:bed.GetGenes()){
			
			String pos = g.getChr()+"_"+g.getStart()+"_"+g.getEnd();
			String outf = outDir+"/"+pos;
			String inf = inDir+"/"+g.getChr()+".maf";
			boolean wrote=false;
			if (fullBed)
			{
				Alignments [] exons = g.getExons();
				for (int ind = 0; ind < exons.length; ind ++) {
					boolean wrote1 =getRegion (exons[ind].getStart(),exons[ind].getEnd(),inf,noPadding, outformat,compress,outf+"_"+ind ,spFilter);
					wrote = wrote | wrote1;
				}
				 
			}
			else
				wrote =getRegion (g.getStart(),g.getEnd(),inf,noPadding, outformat,compress,outf,spFilter);
			if(wrote)
				bw.write(outf+"\n");
		}
		bw.close();
	}

	private static boolean getRegion (int start, int end, String chrMAF,boolean  noPadding,
			String outformat,boolean compress,String outfile,String spFilter) throws IOException, ParseException{
		
		
		MultipleAlignmentIO maio = null; 
		maio = MultipleAlignmentIOFactory.create(outformat);
		MAFIO mafio = new MAFIO();
		MAFAlignment maf = mafio.load(chrMAF, new ArrayList<String>(), start, end) ;
		MultipleAlignment outAln = maf;
		if(compress) {
			System.err.println("compressing alignment");
			maf.compress();
		}
		//System.err.println("Start " + start + " end " + end + " loaded start " + maf.getReferenceStart() + " loaded end " + maf.getReferenceEnd());
		if(!noPadding) {
			maf.padN(start, end);
		} else {
			outAln = maf.toMultipleAlignment(false);
		}
		List<String>  seqNames= maf.getAlignedSequenceIds();
		if(spFilter=="" || seqNames.contains( spFilter)){
				
			outAln.setIOHelper(maio);
			//System.err.println(maf.getReference().getSequenceBases());
			if (outAln.isEmpty())
				return false;
			else{
				BufferedWriter bw = new BufferedWriter(new FileWriter(outfile));
				outAln.write(bw);
				bw.close();
			}
			return true;
		}
		return false;
	}
	
	
	
	
	public static void main(String [] args) throws IOException, ParseException {
		ArgumentMap argMap = CLUtil.getParameters(args, USAGE);
		if ("1".equals(argMap.getTask())) {
			int start = argMap.getInteger("start");
			int end   = argMap.getInteger("end");
			boolean  noPadding = argMap.containsKey("noPadding");
			MultipleAlignmentIO maio = null; 
			if(argMap.containsKey("outformat")) {
				maio = MultipleAlignmentIOFactory.create(argMap.get("outformat"));
			}
			String in  = argMap.getInput();
			MAFIO mafio = new MAFIO();
			
			MAFAlignment maf = mafio.load(in, new ArrayList<String>(), start, end) ;
			MultipleAlignment outAln = maf;
			if(argMap.isPresent("compress")) {
				System.err.println("compressing alignment");
				maf.compress();
			}
			//System.err.println("Start " + start + " end " + end + " loaded start " + maf.getReferenceStart() + " loaded end " + maf.getReferenceEnd());
			if(!noPadding) {
				maf.padN(start, end);
			} else {
				outAln = maf.toMultipleAlignment(true);
			}
			
			outAln.setIOHelper(maio);
			//System.err.println(maf.getReference().getSequenceBases());
			BufferedWriter bw = argMap.getOutputWriter();
			outAln.write(bw);
			bw.close();
		} else if("2".equals(argMap.getTask())) {
			String format = argMap.containsKey("format") ? argMap.get("format") : "BED";
			AnnotationReader<? extends GenomicAnnotation> annotationReader = 
				AnnotationReaderFactory.create(argMap.getMandatory("annotations"), format);
			annotationReader.merge();
			List<String> seqsToLoad = argMap.getAll("seqsToLoad");
			if(argMap.isPresent("guideTree")) {
				File guideTreeFile = new File(argMap.getMandatory("guideTree"));
				NHXParser nhxParser = new NHXParser();
				nhxParser.setSource(guideTreeFile);
				Phylogeny [] phylogenyArr = nhxParser.parse();
				if(phylogenyArr == null || phylogenyArr.length ==0) {
					System.err.println("The guide tree file provided contained no tree");
				} else if (phylogenyArr.length > 1) {
					System.err.println("Guide tree file contains more than one tree. Using first");
				}
				Phylogeny guideTree = phylogenyArr[0];
				String [] seqs = guideTree.getAllExternalSeqNames();
				seqsToLoad = new ArrayList<String>();
				for(int i = 0; i < seqs.length; i++) {
					seqsToLoad.add(seqs[i]);
				}
			}
			MAFIO mafio = new MAFIO(argMap.getInput(), true);
			MAFAlignment combined = new MAFAlignment();
			
			Iterator<? extends GenomicAnnotation> it = annotationReader.getAnnotationList().iterator();
			while(it.hasNext()) {
				GenomicAnnotation annot = it.next();
				MAFAlignment annotationAlignment = mafio.load(seqsToLoad, annot.getStart(), annot.getEnd());
				if(annotationAlignment == null) {
					System.err.println("AnnotationAlignment for " + annot.getLocationString() + " is empty");
				} else {
					//System.err.println("annotation " + annot.getLocationString() + " annotationAlignmentReference " + ( annotationAlignment.getReference() == null ? "NO REFERENCE SEQUENCE :-(": annotationAlignment.getReference().getLocationString()));
				}
				if( annotationAlignment != null && !annotationAlignment.isEmpty() && annotationAlignment.getReferenceChromosome().equals(annot.getChromosome()) ) {
					if(combined.getReferenceId() == null) {
						combined.setReferenceId(annotationAlignment.getReferenceId());
					}
					combined.append(annotationAlignment);
				}
				
			}
			mafio.destroyFileHandle();
			MultipleAlignment outAln = combined;
			if(argMap.containsKey("outformat") ) {
				outAln = combined.toMultipleAlignment(false);
				outAln.setIOHelper(MultipleAlignmentIOFactory.create(argMap.get("outformat")));
			} 
			if(argMap.containsKey("compress")) {
				outAln.compress();
			}
			BufferedWriter bw = argMap.getOutputWriter();
			//System.err.println("Before writing");
			
			outAln.write(bw, seqsToLoad);
			//System.err.println("Wrote alignment");
			bw.close();
			
		} else if("3".equals(argMap.getTask())) {
			String in = argMap.getInput();
			String out = argMap.getOutput();

			MAFAlignment maf = new MAFAlignment();
			maf.createIndex(in);
			maf.writeIndex(out);
			
		}
		   
		else if ("ExtractBed".equals(argMap.getTask())){
			String inDir= argMap.getInput();
			String outDir = argMap.getOutput();
			String annotations =argMap.getMandatory("annotations");
			boolean  noPadding = argMap.containsKey("noPadding");
			String outformat = argMap.containsKey("outformat")? argMap.get("outformat"): "MAF";
			boolean compress = argMap.containsKey("compress");
			List<String> seqsToLoad = argMap.getAll("seqsToLoad");
			String spFilter = argMap.containsKey("spFilter")? argMap.get("spFilter"): "";
			boolean  fullBed = argMap.containsKey("fullBed");
			extractBedRegions (inDir,outDir,annotations,noPadding,outformat,compress,seqsToLoad,spFilter,fullBed);	
			
		} else {
			System.err.println(USAGE);
		}
	}

	

	
	
	
}
