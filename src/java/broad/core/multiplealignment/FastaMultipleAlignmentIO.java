package broad.core.multiplealignment;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.Iterator;
import java.util.List;

import broad.core.multiplealignment.MultipleAlignment.AlignedSequence;

import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;

public class FastaMultipleAlignmentIO implements MultipleAlignmentIO {

	public MultipleAlignment load(String source) throws IOException {
		MultipleAlignment ma = new MultipleAlignment();
		FastaSequenceIO fsio = new FastaSequenceIO(new File(source));
		Iterator<Sequence> it = fsio.loadAll().iterator();
		while(it.hasNext()) {
			AlignedSequence seq = new AlignedSequence(it.next());
			ma.addSequence(seq);
		}
		
		return ma;
	}
	
	public MultipleAlignment load(String source, List<String> sequencesToLoad) throws IOException {
		MultipleAlignment ma = new MultipleAlignment();
		FastaSequenceIO fsio = new FastaSequenceIO(new File(source));
		Iterator<Sequence> it = fsio.extractRecords(sequencesToLoad).iterator();
		while(it.hasNext()) {
			AlignedSequence seq = new AlignedSequence(it.next());
			ma.addSequence(seq);
		}
		
		return ma;
	}
	
	public void write(BufferedWriter bw, MultipleAlignment ma) throws IOException {
		FastaSequenceIO fsio = new FastaSequenceIO();
		fsio.write(ma.getAlignedSequences(), bw);
	}
	
	public void write(BufferedWriter bw, MultipleAlignment ma, List<String> orderOfSequences) throws IOException {
		FastaSequenceIO fsio = new FastaSequenceIO();
		if(orderOfSequences != null && orderOfSequences.size() > 0) {
			Iterator<String> orderIt = orderOfSequences.iterator();
			while(orderIt.hasNext()) {
				String seqId = orderIt.next();
				fsio.write(ma.getAlignedSequence(seqId),bw);
			}
		}else {
			fsio.write(ma.getAlignedSequences(), bw);
		}
	}	

	public MultipleAlignment load(InputStream in) throws IOException {
		MultipleAlignment ma = new MultipleAlignment();
		FastaSequenceIO fsio = new FastaSequenceIO();
		Iterator<Sequence> it = fsio.loadAll(in).iterator();
		while(it.hasNext()) {
			AlignedSequence seq = new AlignedSequence(it.next());
			ma.addSequence(seq);
		}
		
		return ma;
	}
	
	public MultipleAlignment load(InputStream in, List<String> sequencesToLoad) throws IOException {
		MultipleAlignment ma = new MultipleAlignment();
		FastaSequenceIO fsio = new FastaSequenceIO();
		Iterator<Sequence> it = fsio.extractRecords(sequencesToLoad, in).iterator();
		while(it.hasNext()) {
			AlignedSequence seq = new AlignedSequence(it.next());
			ma.addSequence(seq);
		}
		
		return ma;
	}

	public String getPreferredFileExtension() {
		return "fa";
	}

}
