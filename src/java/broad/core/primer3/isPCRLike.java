package broad.core.primer3;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import nextgen.core.alignment.SmatchLike;
import nextgen.core.annotation.BasicAnnotation;
import nextgen.core.annotation.Gene;
import nextgen.core.annotation.Annotation.Strand;


import broad.core.datastructures.Pair;
import broad.core.motif.SearchException;
import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.core.sequence.SequenceRegion;
import broad.pda.datastructures.Alignments;

public class isPCRLike {

	Pair<String> primers;
	List<Sequence> sequences;
	int minPerfectMatch=10;
	
	public isPCRLike(Pair<String> primers, List<Sequence> sequences){
		this.primers=primers;
		this.sequences=sequences;
	}
	
	/**
	 * Get all possible amplicons as subsequences of the original sequences
	 * Name of subsequence is equal to the name of the parent sequence
	 * @return Collection of possible amplicons
	 * @throws SearchException
	 */
	public Collection<Sequence> getAllPossibleAmplicons() throws SearchException {
		Collection<Alignments> a = getAllPossibleAmpliconsAsAlignments();
		Map<String, Sequence> seqsByName = new TreeMap<String, Sequence>();
		for(Sequence seq : sequences) {
			seqsByName.put(seq.getId(), seq);
		}
		Collection<Sequence> rtrn = new ArrayList<Sequence>();
		for(Alignments al : a) {
			BasicAnnotation b = new BasicAnnotation(al.getChr(), al.getStart(), al.getEnd(), Strand.POSITIVE);
			Sequence seq = seqsByName.get(b.getChr());
			Sequence subseq = seq.getSubsequence(b);
			subseq.setId(seq.getId());
			rtrn.add(subseq);
		}
		return rtrn;
	}
	
	private Collection<Alignments> getAllPossibleAmpliconsAsAlignments() throws SearchException{
		Collection<Alignments> rtrn=new TreeSet<Alignments>();
		
		String kmer=Sequence.get3Prime(primers.getValue1(), this.minPerfectMatch);
		SmatchLike smatch=new SmatchLike(kmer, sequences, kmer.length());
		Map<String, List<SequenceRegion>> matches=smatch.getForwardTargetRegions(0);
		
		
		kmer=Sequence.reverseSequence(Sequence.get3Prime(primers.getValue2(), this.minPerfectMatch));
		smatch=new SmatchLike(kmer, sequences, kmer.length());
		Map<String, List<SequenceRegion>> reverseMatches=smatch.getForwardTargetRegions(0);
		
		for(String forward: matches.keySet()){
			if(reverseMatches.containsKey(forward)){
				List<SequenceRegion> forwardRegions=matches.get(forward);
				List<SequenceRegion> reverseRegions=reverseMatches.get(forward);
				for(SequenceRegion forwardRegion: forwardRegions){
					for(SequenceRegion reverseRegion: reverseRegions){
						Alignments amplicon=getAmplicon(forwardRegion, reverseRegion);
						if(amplicon!=null){rtrn.add(amplicon);}
					}
				}
			}
		}
		return rtrn;
	}

	private Alignments getAmplicon(SequenceRegion forwardRegion, SequenceRegion reverseRegion) {
		//make sure that its a valid amplicon
		//if forward is before reverse--> valid
		
		Alignments forwardAlign=getAlignments(forwardRegion.getContainingSequenceId());
		Alignments reverseAlign=getAlignments(reverseRegion.getContainingSequenceId());
		
		if(forwardAlign.getStart()>reverseAlign.getStart() || forwardAlign.getEnd()>reverseAlign.getEnd()){return null;}
				
		Alignments rtrn=new Alignments(forwardAlign.getChr(), forwardAlign.getStart(), reverseAlign.getEnd());
		
		return rtrn;
	}
	
	//public Collection<Alignments> getAllLikelyAmplicons(){}
	
	private Alignments getAlignments(String containingSequenceId) {
		String[] tokens=containingSequenceId.split(":");
		String name=tokens[0];
		String other=tokens[1];
		tokens=other.split("-");
		int start=new Integer(tokens[0]);
		int end=new Integer(tokens[1]);
		return new Alignments(name, start, end);
	}

	private static List<Sequence> initializeGeneSequence(String geneSequence) throws IOException {
		FastaSequenceIO fsio = new FastaSequenceIO(geneSequence);
		return fsio.loadAll();
	}
	
	public static void main(String[] args) throws IOException, SearchException{
		Pair<String> primers=new Pair<String>("GAGTTCTGCGGAGGGATGGCA", Sequence.reverseSequence("ATCCTTTCCTCTGCCCCCAGGTCC"));
		List<Sequence> geneSequences=initializeGeneSequence(args[0]);
		isPCRLike isPCR=new isPCRLike(primers, geneSequences);
		Collection<Alignments> amplicons=isPCR.getAllPossibleAmpliconsAsAlignments();
		
		for(Alignments amplicon: amplicons){
			System.err.println(amplicon);
		}
	}
	
}
