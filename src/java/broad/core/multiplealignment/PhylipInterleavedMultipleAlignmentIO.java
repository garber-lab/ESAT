package broad.core.multiplealignment;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;


import broad.core.error.ParseException;
import broad.core.multiplealignment.MultipleAlignment.AlignedSequence;
import broad.core.util.CLUtil;

public class PhylipInterleavedMultipleAlignmentIO implements MultipleAlignmentIO {
	static final int DEFAULT_SPACING    = 9;
	static final int MAX_CHARS_PER_LINE = 90;
	static final int SEQ_NAME_LENGTH    = 10;
	
	private int spacing = DEFAULT_SPACING;
	
	public MultipleAlignment load(String source) throws IOException, ParseException {
		FileInputStream fis = new FileInputStream(new File(source));
		MultipleAlignment m = load(fis);
		fis.close();
		return m;
	}
	
	public MultipleAlignment load(String fileName, List<String> sequencesToLoad) throws IOException, ParseException {
		FileInputStream fis = new FileInputStream(new File(fileName));
		MultipleAlignment m = load(fis, sequencesToLoad);
		fis.close();
		return m;
	}
	
	public MultipleAlignment load(InputStream in) throws IOException, ParseException {
		return load(in, new ArrayList<String>());
	}
	
	public MultipleAlignment load(InputStream in,  List<String> sequencesToLoad) throws IOException {
		BufferedReader br = new BufferedReader(new InputStreamReader(in));
		String line = br.readLine();
		String [] headerInfo = line.split("\\s+");
		int seqNum = Integer.parseInt(headerInfo[1]);
		boolean isFirstBatch = true;
		int seqCounter = 0;
		MultipleAlignment ma = new MultipleAlignment();
		String [] alignedSeqIds = new String[seqNum];
		
		while( (line = br.readLine()) != null) {
			if(line.trim().length() == 0) {
				//System.out.println("empty line, skipping");
				continue;
			}

			if(isFirstBatch) {
				String seqId = line.substring(0,10).trim();
				line = line.substring(10);
				alignedSeqIds[seqCounter] = seqId;
				AlignedSequence alignedSequence = new AlignedSequence(seqId);
				alignedSequence.setId(seqId);
				ma.addSequence(alignedSequence);
				isFirstBatch = seqCounter < seqNum - 1;
				//System.out.println("Seq " + seqCounter +": " + seqId);
			}
			String [] seqAlignedBlocks = line.split("\\s+");			
			AlignedSequence alignedSequence = ma.getAlignedSequence(alignedSeqIds[seqCounter]);
			for(int i = 0 ; i < seqAlignedBlocks.length; i++) {
				alignedSequence.append(seqAlignedBlocks[i]);
			}
			
			//System.out.print("Seqcounter " + seqCounter );
			seqCounter = (seqCounter + 1) % seqNum;    
			//System.out.println(" next: " +seqCounter);
				
		}
		
		return ma;
	}
	
	public void setSpacing(int spacing) {
		this.spacing = spacing;
	}

	public void write(BufferedWriter bw, MultipleAlignment ma) throws IOException {
		List<AlignedSequence> seqs = ma.getAlignedSequences();
		
		int seqNum = seqs.size();
		if(seqs.size() == 0) {
			return;
		}
		int nucleotideNum = seqs.get(0).getSequenceBases().length();
		
		bw.write("\t"+seqNum+"\t"+nucleotideNum);
		bw.newLine();
		
		// Print first line
		for(int i = 0; i < seqs.size(); i++) {
			AlignedSequence seq = seqs.get(i);
			CLUtil.writeLeftJustifiedField(bw, seq.getName(), SEQ_NAME_LENGTH);
			writeLine(bw, seq.getSequenceBases(), 0);
			bw.newLine();
		}
		bw.newLine();
		
		int line = 1;
		int lineNum = (int) Math.ceil(nucleotideNum/(float)MAX_CHARS_PER_LINE);
		while(line < lineNum) {
			for(int i = 0; i < seqs.size(); i++) {
				AlignedSequence seq = seqs.get(i);
				//System.out.print("\t"+seq.getContainingSequenceId());
				writeLine(bw, seq.getSequenceBases(), line);
				bw.newLine();
			}	
			bw.newLine();
			line++;
		}
		
	}
	
	public void write(BufferedWriter bw, MultipleAlignment ma, List<String> order) throws IOException {
		if(order == null || order.size() == 0) {
			write(bw, ma);
			return;
		}
		List<AlignedSequence> seqs = ma.getAlignedSequences();
		
		int seqNum = order.size();
		int nucleotideNum = seqs.get(0).getSequenceBases().length();
		bw.write("\t"+seqNum+"\t"+nucleotideNum);
		bw.newLine();
		
		// Print first line
		for(int i = 0; i < seqNum; i++) {
			AlignedSequence seq = ma.getAlignedSequence(order.get(i));
			CLUtil.writeLeftJustifiedField(bw, seq.getName(), SEQ_NAME_LENGTH);
			writeLine(bw, seq.getSequenceBases(), 0);
			bw.newLine();
		}
		bw.newLine();
		
		int line = 1;
		int lineNum = (int) Math.ceil(nucleotideNum/(float)MAX_CHARS_PER_LINE);
		while(line < lineNum) {
			for(int i = 0; i < seqNum; i++) {
				AlignedSequence seq =  ma.getAlignedSequence(order.get(i));
				writeLine(bw, seq.getSequenceBases(), line);
				bw.newLine();
			}	
			bw.newLine();
			line++;
		}
		
	}
	
	private int  writeLine(BufferedWriter bw, String sequence, int lineNum) throws IOException {
		int start = lineNum * MAX_CHARS_PER_LINE;
		int end   = (int) Math.min(start + MAX_CHARS_PER_LINE, sequence.length());
		//System.out.print("writing extracting sequence "+ sequence + " "+start+ "-"+end);
		String seqToWrite = sequence.substring(start, end);
		int chunkNum = (int) Math.ceil(seqToWrite.length()/(float)spacing);
		int pos = 0;
		for(int i = 0; i < chunkNum; i++) {		
			if(i < chunkNum - 1) {
				bw.write(seqToWrite.substring(pos, pos + spacing));
				bw.write(" ");
			} else {
				bw.write(seqToWrite.substring(pos));
			}
			pos = pos + spacing;
		}
		
		return pos;
	}

	public String getPreferredFileExtension() {
		return "iphy";
	}



}
