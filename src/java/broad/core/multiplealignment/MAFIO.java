package broad.core.multiplealignment;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.List;
import java.util.Stack;


import broad.core.error.ParseException;
import broad.core.multiplealignment.MAFAlignment.MAFHeader;
import broad.core.multiplealignment.MAFAlignment.MAFMultipleAlignmentBlock;
import broad.core.multiplealignment.MultipleAlignment.AlignedSequence;

public class MAFIO implements MultipleAlignmentIO {
	RandomAccessFile fileHandle;
	MAFAlignment alignment;
	String alignmentFile;

	public MAFIO(String alignmentFile, boolean keepFileHandle) throws IOException, ParseException {
		alignment = createUnloadedAlignment(alignmentFile);
		this.alignmentFile = alignmentFile;
		if(keepFileHandle) {
			fileHandle = new RandomAccessFile(alignmentFile, "r");
		}
	}
	
	public MAFIO() {
		super();
	}
	public String getPreferredFileExtension() {
		return ".maf";
	}
	
	public void destroyFileHandle() throws IOException {
        if(fileHandle != null) {
            fileHandle.close();
        }
    }

	public MAFAlignment load(List<String> sequencesToLoad, int start, int end) throws IOException, ParseException {
		//RandomAccessFile handle = fileHandle;
		boolean closeFile = false;
		if(fileHandle == null ) { // It would be nice to test whether the handle is closed if not null....
            closeFile = true;
			fileHandle = new RandomAccessFile(alignmentFile, "r");
		}
		try {
			alignment.load(fileHandle, start, end, sequencesToLoad);
		} finally {
			if(closeFile) {
                fileHandle.close();
                fileHandle = null;
            }
		}
		return alignment;
	}
	
	public MAFAlignment createUnloadedAlignment(String fileName) throws IOException, ParseException {
		MAFAlignment aln = new MAFAlignment();
		String idxFileName = fileName + ".index";
		File idxFile = new File(idxFileName);
		if(!idxFile.exists()) {
			aln.createIndex(fileName);
			aln.writeIndex(idxFileName);
		} else {
			aln = new MAFAlignment(idxFileName);
		}
		
		return aln;
	}

	public MAFAlignment load(String fileName) throws IOException, ParseException {
		return load(fileName, new ArrayList<String>());	
	}
	
	public MAFAlignment load(String fileName, List<String> sequencesToLoad) throws IOException, ParseException {
		MAFAlignment aln = new MAFAlignment();
		String idxFileName = fileName + ".index";
		File idxFile = new File(idxFileName);
		if(!idxFile.exists()) {
			aln.createIndex(fileName);
			aln.writeIndex(idxFileName);
		} else {
			aln = new MAFAlignment(idxFileName);
		}
		
		RandomAccessFile raf = new RandomAccessFile(fileName,"r");
		aln.load(raf, sequencesToLoad);
		raf.close();
		return aln;
	}
	
	public MultipleAlignment load(InputStream in) throws IOException, ParseException {
		return load(in, new ArrayList<String>());
	}
	
	public MAFAlignment load(String fileName, List<String> sequencesToLoad, int start, int end) throws IOException, ParseException {
		MAFAlignment aln = new MAFAlignment();
		String idxFileName = fileName + ".index";
		File idxFile = new File(idxFileName);
		if(!idxFile.exists()) {
			System.out.print("Index file not exists, creating and writing it: " + idxFileName);
			aln.createIndex(fileName);
			aln.writeIndex(idxFileName);
			System.out.println("   Done writing index");
		} else {
			aln = new MAFAlignment(idxFileName);
		}
		
		RandomAccessFile raf = new RandomAccessFile(fileName,"r");
		aln.load(raf,start, end,sequencesToLoad);
		raf.close();
		return aln;
	}
	
	public MultipleAlignment load(InputStream in, List<String> sequencesToLoad) throws IOException, ParseException {
		MAFAlignment aln = new MAFAlignment();
		MAFHeader header = aln.getHeader();
		Stack<MAFMultipleAlignmentBlock> blocks = new Stack<MAFMultipleAlignmentBlock>();
		ArrayList<String> seqIds = new ArrayList<String>();
		
		BufferedReader br = new BufferedReader(new InputStreamReader(in));
		String [] lineInfo = null;
		String line = br.readLine();
		if(line == null || !line.startsWith("##maf ")) {
			throw new ParseException("First line in the file " + line + " was either null or did not start with \"##maf \" as specified");
		}			
		
		header.setVariablesFromRawData(line.substring(6).split("\\s"));
		
		line = br.readLine();
		if(line == null || !line.startsWith("# ")) {
			throw new ParseException("Second line in the file " + line + " was either null or did not start with \"# \" as specified");
		}
		header.setRunParameters(line.substring(2));
		
		while((line = br.readLine()) != null) {
			//Ignore all comment lines
			if(line.startsWith("#") || line.trim().length() == 0){
				continue;
			}
			
			if(line.startsWith("a ")) {		
				MAFMultipleAlignmentBlock ma = new MAFMultipleAlignmentBlock();
				blocks.push(ma);
				lineInfo = line.substring(2).split("\\s");
				ma.setAlignmentInfoFromRawData(lineInfo);
			} else if(line.startsWith("s ")) {
				//System.out.println("\tAlignment aligned seq " + line);
				MAFMultipleAlignmentBlock ma = blocks.peek();
				lineInfo = line.substring(2).split("\\s+");
				AlignedSequence alnSeq = ma.addSequenceFromRawData(lineInfo);
				if(sequencesToLoad !=null && sequencesToLoad.size() > 0 && !sequencesToLoad.contains(alnSeq.getId()) ) {
					//System.out.println("sequence " + alnSeq.getId() + " is not in sequencesToLoad, ignoring"); 
					continue;
				}
				
				//System.out.println("ma  getReferenceId? " + ma.getReferenceId());
				
				if(aln.getReferenceId() == null) {
					System.err.println("Ref is " + ma.getReferenceId());
					aln.setReferenceId(alnSeq.getId());
				} 
				
				if(ma.getReferenceId() != null) {
					ma.setReferenceId(aln.getReferenceId());
				}

				if(!seqIds.contains(alnSeq.getId())) {
					seqIds.add(alnSeq.getId());
				}
			} else if (line.startsWith("i ") || line.startsWith("e ")) {
				System.err.println("Information blocks (i/e) are not yet supported");
			} else {
				throw new ParseException("Invalid alignment line <"+ line +">");
			}
			
		}
		
		aln.setBlocks(blocks);
		aln.setAlignedSequences(seqIds);
		return aln.toMultipleAlignment();
	}

	public void write(BufferedWriter bw, MultipleAlignment ma)
			throws IOException {
		// TODO Auto-generated method stub

	}
	
	public void write(BufferedWriter bw, MultipleAlignment ma, List<String> sequenceOrder)
		throws IOException {
		// TODO Auto-generated method stub
	}


}
