package broad.core.siphy;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.apache.log4j.Logger;
import org.forester.phylogeny.Phylogeny;


import Jama.Matrix;
import broad.core.error.ParseException;
import broad.core.multiplealignment.MAFAlignment;
import broad.core.multiplealignment.MAFIO;
import broad.core.multiplealignment.MultipleAlignment;
import broad.core.multiplealignment.MultipleAlignmentFactory;
import broad.core.multiplealignment.MultipleAlignment.AlignedSequence;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.feature.genome.Chromosome;

public class ConservationUtils {
	static Logger logger = Logger.getLogger(ConservationUtils.class.getName());
	
	public static void setUninformativeNodes(Map<String, Matrix> column, List<String> sequences) {
		setUninformativeNodes(column, sequences, 0);
		
	}
	
	public static void setUninformativeNodes(Map<String, Matrix> column, List<String> sequences, int colIdx ) {
		Iterator<String> seqIt = sequences.iterator();
		
		while(seqIt.hasNext()) {
			String seq = seqIt.next();
			Matrix leaf = column.get(seq);
			
			for(int i = 0; i < leaf.getRowDimension(); i++) {
				leaf.set(i,colIdx,1);
			}
		}
		
	}
	
	public static Phylogeny pruneTree(List<String> toPrune, Phylogeny tree) {
		Phylogeny prunned = tree.copy();
		
		//System.out.println("Original tree external nodes: " + tree.getAllExternalSeqNames() + "\nprunning...\n" + toPrune);
		Iterator<String> nodeNameIt = toPrune.iterator();
		while(nodeNameIt.hasNext()) {
			String nodeName = nodeNameIt.next();
			//System.out.println("ToPrune: " + nodeName);
			Collection<String> nodesForName = prunned.getNodes(nodeName);
			if(nodesForName != null && !nodesForName.isEmpty()) {
				prunned.removeExtNode(prunned.getNode(nodeName));
			} else {
				logger.debug("Node " + nodeName + " was not in tree");
			}
		}
		//System.out.println("prunned external sequences " + prunned.toNewHampshire(false));
		return prunned;
	}
	
	public static Phylogeny pruneTree(String toPrune, Phylogeny tree) {
		List<String> oneSeqList = new ArrayList<String>(1);
		oneSeqList.add(toPrune);
		
		return pruneTree(oneSeqList, tree);
	}
	
	public static List<String> getGappedSeqsInWindowMatrix(int window, Map<String, Matrix> alignmentMatrix, int start) {
		List<String> gappedSeqs = new ArrayList<String>();
		// If total tree branch length is too short already, we can skip to full alignment it aint goint to get any longer.
		
			Iterator<String> seqIdIt = alignmentMatrix.keySet().iterator();
			while(seqIdIt.hasNext()) {
				String seqId = seqIdIt.next();
				Matrix seqBases = alignmentMatrix.get(seqId);
				for(int j = start; j < window + start; j++) {
					boolean isGap = true;
					for( int l = 0; l < seqBases.getRowDimension(); l++ ) {
						if(seqBases.get(l,j) > 0) {
							isGap = false;
						}
					}
					if(isGap && !gappedSeqs.contains(seqId)) {
						gappedSeqs.add(seqId);
						break;
					}
				}
			}
		
		
		return gappedSeqs;
	}
	
	public static MultipleAlignment setUpAlignment(ArgumentMap argMap, String alnFile, 
			String alnFileFormat, List<String> ignoreList, EvolutionaryModel model) 
	throws IOException, ParseException,	FileNotFoundException {
		MultipleAlignment alignment;
		if(!"MAF".equalsIgnoreCase(alnFileFormat)){
			alignment = MultipleAlignmentFactory.create(alnFile, alnFileFormat);
			//scaler.alignment.encodeAsMatrix();
			alignment.setReferenceId(argMap.getMandatory("ref"));
			alignment.compress();
			Iterator<AlignedSequence> it = alignment.getAlignedSequences().iterator();
			while(it.hasNext()) {
				AlignedSequence seq = it.next();
				seq.setEnd(seq.getLength());
			}
		} else {
			int start = argMap.getInteger("start");
			int end   = argMap.getInteger("end");
			logger.debug("Getting alignment from " + start + " to " + end);
			alignment = setUpMAF(alnFile, ignoreList, model,start,end);
		}
		alignment.remove(ignoreList);
		
		return alignment;
	}
	
	

	public static MultipleAlignment setUpMAF(String alnFile,
			List<String> ignoreList, EvolutionaryModel model, 
			int start, int end) throws IOException, ParseException,
			FileNotFoundException {
		MultipleAlignment alignment;
		String mafIndex = alnFile + ".index";
		MAFAlignment mafAln = new MAFAlignment( mafIndex);
		String [] seqs = model.getTree().getAllExternalSeqNames();
		List<String> seqsToLoad = new ArrayList<String>();
		for(int i = 0; i < seqs.length; i++) {
			//if(!ignoreList.contains(seqs[i])) { //TODO: Reinstate this condition only commented out to test
				seqsToLoad.add(seqs[i]);
			//}
		}
		RandomAccessFile alnRaf = new RandomAccessFile(alnFile , "r");		
		mafAln.load(alnRaf, start, end, seqsToLoad);
		mafAln.compress();
		//mafAln.getReference().setEnd(end);
		//BufferedWriter tbw = new BufferedWriter(new OutputStreamWriter(System.err));
		//mafAln.setIOHelper(MultipleAlignmentIOFactory.create("PHYLIP"));
		//mafAln.write(tbw);
		//tbw.flush();
		alnRaf.close();
		
		if(mafAln.isEmpty()) {
			return mafAln;
		}
		alignment = mafAln.toMultipleAlignment();
		return alignment;
	}
	
	public static MultipleAlignment setUpMAF(MAFIO mafio,
			List<String> ignoreList, EvolutionaryModel model, 
			int start, int end) throws IOException, ParseException,
			FileNotFoundException {
		return setUpMAF(mafio, ignoreList, model.getTree(), start, end);
	}
	
	public static MultipleAlignment setUpMAF(MAFIO mafio,
			List<String> ignoreList, Phylogeny guideTree, 
			int start, int end) throws IOException, ParseException,
			FileNotFoundException {
		MultipleAlignment alignment;
		String [] seqs = guideTree.getAllExternalSeqNames();
		List<String> seqsToLoad = new ArrayList<String>();
		for(int i = 0; i < seqs.length; i++) {
			//if(!ignoreList.contains(seqs[i])) {
				seqsToLoad.add(seqs[i]);
			//}
		}
		MAFAlignment mafAln = mafio.load(seqsToLoad, start, end);
		mafAln.compress();
		
		//System.err.println("Is aln empty? " + mafAln.isEmpty());
		if(mafAln.isEmpty()) {
			return mafAln;
		}
		alignment = mafAln.toMultipleAlignment();
		return alignment;
	}
	
	public static List<String> commaSeparatedStringToList(String csvListStr) {
		List<String> list = new ArrayList<String>();
		if(csvListStr != null && csvListStr.trim().length() > 0) {
			String [] ignoreListArr = csvListStr.split(",");
			for(int i = 0; i < ignoreListArr.length; i++ ) {
				list.add(ignoreListArr[i]);
			}
		}
		return list;
	}

}
