package broad.core.multiplealignment;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import broad.core.multiplealignment.MultipleAlignment.AlignedSequence;

import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;

public class PhylipSequencialMultipleAlignmentIO implements MultipleAlignmentIO {

	public MultipleAlignment load(String fileName) throws IOException {
		// TODO Auto-generated method stub
		return null;
	}

	public MultipleAlignment load(InputStream in) throws IOException {
		// TODO Auto-generated method stub
		return null;
	}
	
	public MultipleAlignment load(String fileName, List<String> sequencesToLoad) throws IOException {
		// TODO Auto-generated method stub
		return null;
	}

	public MultipleAlignment load(InputStream in, List<String> sequencesToLoad) throws IOException {
		// TODO Auto-generated method stub
		return null;
	}

	public void write(BufferedWriter bw, MultipleAlignment ma) throws IOException {
		int seqNum = ma.getAlignedSequenceIds().size();
		int alnLength = ma.getReferenceId() != null ? 
				ma.getAlignedSequence(ma.getReferenceId()).getLength() :
				ma.getAlignedSequences().get(0).getLength();
		bw.write("\t" + seqNum + "\t" + alnLength);
		bw.newLine();
		PhylipSequenceIO psio = new PhylipSequenceIO();
		psio.write(ma.getAlignedSequences(), bw);
		
	}
	
	public void write(BufferedWriter bw, MultipleAlignment ma, List<String> order) throws IOException {
		if(order == null || order.size() ==0) {
			write(bw, ma);
			return;
		}
		int seqNum = ma.getAlignedSequenceIds().size();
		int alnLength = ma.getReferenceId() != null ? 
				ma.getAlignedSequence(ma.getReferenceId()).getLength() :
				ma.getAlignedSequences().get(0).getLength();
		bw.write("\t" + seqNum + "\t" + alnLength);
		bw.newLine();
		PhylipSequenceIO psio = new PhylipSequenceIO();
		List<AlignedSequence> orderedSeqs = new ArrayList<AlignedSequence>(order.size());
		Iterator<String> orderIt = order.iterator();
		while(orderIt.hasNext()) {
			orderedSeqs.add(ma.getAlignedSequence(orderIt.next()));
		}
		psio.write(orderedSeqs, bw);
		
	}
	
	static class PhylipSequenceIO extends FastaSequenceIO {
		protected void writeSequenceId(Sequence seq, BufferedWriter bw) throws IOException {
			bw.write(seq.getId());
		}
	}

	public String getPreferredFileExtension() {
		return "phy";
	}

}
