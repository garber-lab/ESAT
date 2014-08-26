package umms.core.alignment;


import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.core.math.CombinationGenerator;
import broad.core.motif.SearchException;
import broad.core.motif.SequenceMotif;
import broad.core.primer3.ComputeOriginalScore;
import broad.core.primer3.HairpinKmer;
import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.core.sequence.SequenceRegion;
import broad.core.datastructures.Pair;


//Goal: Use Motif scan to simulate (less efficiently) SMATCH
public class SmatchLike {

	List<Sequence> geneSequences; 
	Collection<String> possibleTargets;
	int seedNumber;
	String kmer;
	//int numMismatches;
	
	public SmatchLike(String kmer, List<Sequence> geneSequences, int seedNum){
		this.kmer=kmer;
		this.geneSequences=geneSequences;
		this.seedNumber=seedNum;
	}
	
	public SmatchLike(String kmer, String primer, int seedNum){
		this.kmer=kmer;
		
		List<Sequence> seqs=new ArrayList<Sequence>();
		Sequence seq=new Sequence("temp");
		seq.setSequenceBases(primer);
		seqs.add(seq);
		
		this.geneSequences=seqs;
		this.seedNumber=seedNum;
	}
	
	public SmatchLike(String kmer2, Collection<Pair<String>> tailPrimers, Pair<String> subprimer1) {
		this.kmer=kmer2;
		this.seedNumber=kmer.length();
		List<Sequence> seqs=new ArrayList<Sequence>();
		
		for(Pair<String> primer: tailPrimers){
			if(!primer.equals(subprimer1)){
				Sequence seq=new Sequence(primer.getValue1());
				seq.setSequenceBases(primer.getValue1());
				seqs.add(seq);
			}
		}
		
		this.geneSequences=seqs;
		
	}
	
//	public SmatchLike(String kmer2, Collection<Pair<String>> tailPrimers, PrimerPair majorPrimer) {
//		this.kmer=kmer2;
//		this.seedNumber=kmer.length();
//		List<Sequence> seqs=new ArrayList<Sequence>();
//		
//		for(Pair<String> primer: tailPrimers){
//			if(!primer.getValue2().equals(majorPrimer.getRightPrimer())){
//				Sequence seq=new Sequence(primer.getValue1());
//				seq.setSequenceBases(primer.getValue1());
//				seqs.add(seq);
//			}
//		}
//		
//		this.geneSequences=seqs;
//		
//	}
	
	
	/*private Collection<String> smatch(String kmer, List<Sequence> geneSequences, int seedNum) throws SearchException{
		Collection<String> possibleTargets=new TreeSet();
		//Seed by enumerating kmers and list regions that have these seeds
		//go through each gene and make a seed by enumerating kmers
		Map<String, List<SequenceRegion>> seededRegion=seed(geneSequences, kmer, seedNum);
		
		//TODO Collapse overlapping SequenceRegions to avoid searching twice
		
		int extension=kmer.toCharArray().length-seedNum;
		
		//for each stored region score full 2 mismatched motif
		possibleTargets=findPossibleTargets(kmer, geneSequences, seededRegion, extension);
		return possibleTargets;
	}*/
	
	

	private Collection<String> smatch(String kmer, List<Sequence> geneSequences, int seedNum, int numMismatch) throws SearchException{
		Collection<String> possibleTargets=new TreeSet();
		Map<String, List<SequenceRegion>> seededRegion=seed(geneSequences, kmer, seedNum);
		//System.err.println("Done seeding...");
		
		int extension=kmer.toCharArray().length-seedNum;
		
		//for each stored region score full 2 mismatched motif
		possibleTargets=findPossibleTargets(kmer, geneSequences, seededRegion, extension, numMismatch);
		
		//System.err.println("done extending ...");
		return possibleTargets;
	}
	
	private Map<String, List<SequenceRegion>> smatchRegions(String kmer, List<Sequence> geneSequences, int seedNum, int numMismatch) throws SearchException{
		Map<String, List<SequenceRegion>> possibleTargets=new TreeMap<String, List<SequenceRegion>>();
		Map<String, List<SequenceRegion>> seededRegion=seed(geneSequences, kmer, seedNum);
		
		int extension=kmer.toCharArray().length-seedNum;
		
		//for each stored region score full 2 mismatched motif
		possibleTargets=findPossibleTargetRegions(kmer, geneSequences, seededRegion, extension, numMismatch);
		
		return possibleTargets;
	}
	
	public Collection<String> getForwardTargets(int numMismatches) throws SearchException{
		Collection<String> forwardMatches=smatch(kmer, geneSequences, seedNumber, numMismatches);
		return forwardMatches;
	}
	
	public Map<String, List<SequenceRegion>> getForwardTargetRegions(int numMismatches) throws SearchException{
		Map<String, List<SequenceRegion>> forwardMatches=smatchRegions(kmer, geneSequences, seedNumber, numMismatches);
		return forwardMatches;
	}
	
	public Collection<String> getReverseTargets(int numMismatches) throws SearchException{
		Collection<String> reverseMatches=smatch(Sequence.reverseSequence(kmer), geneSequences, seedNumber, numMismatches);
		return reverseMatches;
	}
	
	public Collection<String> getAllPossibleTargets(int numMismatches) throws SearchException{
		if(this.possibleTargets!=null){
			return this.possibleTargets;
		}
		else{
			Collection<String> forwardMatches=smatch(kmer, geneSequences, seedNumber, numMismatches);
			Collection<String> reverseMatches=smatch(Sequence.reverseSequence(kmer), geneSequences, seedNumber, numMismatches);
			//System.err.println(kmer+" Forward Scan "+forwardMatches.size()+" Reverse Scan "+reverseMatches.size());
			//if(forwardMatches.size()>0){System.err.println("Forward matches "+forwardMatches);}
			//if(reverseMatches.size()>0){System.err.println("Reverse matches "+reverseMatches);}
			
			this.possibleTargets=(forwardMatches);
			this.possibleTargets.addAll(reverseMatches);
			return this.possibleTargets;
		}
	}
	
	public int numPossibleTargets(){return this.possibleTargets.size();}
	
	private Collection<String> findPossibleTargets(String kmer, List<Sequence> geneSequences, Map<String, List<SequenceRegion>> seededRegion, int extension) throws SearchException {
		Collection<String> rtrn=new TreeSet();
		
		Collection<SequenceMotif> mismatchedMotif=makeSequenceMotifWithMismatches(kmer, 2);
		
		System.err.print("Extending.... ");
		
		//As a first pass test only genes with any hit
		double i=0;
		for(SequenceMotif motif: mismatchedMotif){
			System.err.print(" "+(i/mismatchedMotif.size())*100+"...");
			for(Sequence geneSeq: geneSequences){
				if(seededRegion.containsKey(geneSeq.getId())){
					//for all seeded regions extend by kmer size and store
					List<SequenceRegion> matches = motif.match(geneSeq, seededRegion.get(geneSeq.getId()), extension);
					if(matches.size()>0){rtrn.add(geneSeq.getId());}
				}
			}
			i++;
		}
		
		System.err.println();
		
		return rtrn;
	}
	
	private Collection<String> findPossibleTargets(String kmer, List<Sequence> geneSequences, Map<String, List<SequenceRegion>> seededRegion, int extension, int numMismatch) throws SearchException {
		Collection<String> rtrn=new TreeSet();
		
		Collection<SequenceMotif> mismatchedMotif=makeSequenceMotifWithMismatches(kmer, numMismatch);
		
		//System.err.print("Extending.... ");
		
		//As a first pass test only genes with any hit
		double i=0;
		for(SequenceMotif motif: mismatchedMotif){
			//System.err.print(" "+(i/mismatchedMotif.size())*100+"...");
			for(Sequence geneSeq: geneSequences){
				if(seededRegion.containsKey(geneSeq.getId())){
					//for all seeded regions extend by kmer size and store
					List<SequenceRegion> matches = motif.match(geneSeq, seededRegion.get(geneSeq.getId()), extension);
					if(matches.size()>0){rtrn.add(geneSeq.getId());}
				}
			}
			i++;
		}
		
		//System.err.println();
		
		return rtrn;
	}
	
	private Map<String, List<SequenceRegion>> findPossibleTargetRegions(String kmer, List<Sequence> geneSequences, Map<String, List<SequenceRegion>> seededRegion, int extension, int numMismatch) throws SearchException {
		Map<String, List<SequenceRegion>> rtrn=new TreeMap<String, List<SequenceRegion>>();
		
		Collection<SequenceMotif> mismatchedMotif=makeSequenceMotifWithMismatches(kmer, numMismatch);
		
		//System.err.print("Extending.... ");
		
		//As a first pass test only genes with any hit
		double i=0;
		for(SequenceMotif motif: mismatchedMotif){
			//System.err.print(" "+(i/mismatchedMotif.size())*100+"...");
			for(Sequence geneSeq: geneSequences){
				if(seededRegion.containsKey(geneSeq.getId())){
					//for all seeded regions extend by kmer size and store
					List<SequenceRegion> matches = motif.match(geneSeq, seededRegion.get(geneSeq.getId()), extension);
					if(matches.size()>0){rtrn.put(geneSeq.getId(), matches);}
				}
			}
			i++;
		}
		
		//System.err.println();
		
		return rtrn;
	}
	
	//TODO: Replace with SMATCH which might be much faster
	private static Collection<SequenceMotif> makeSequenceMotifWithMismatches(String kmerSequence, int num) throws SearchException {
		Collection<SequenceMotif> rtrn=new ArrayList();
		
		char[] chars=kmerSequence.toCharArray();
		CombinationGenerator comb=new CombinationGenerator(chars.length, num);
		
		
		while(comb.hasMore()){
			int[] combination=comb.getNext();
			String degenerate=string(chars, combination);
			rtrn.add(new SequenceMotif(degenerate, 1));
		}
				
		return rtrn;
	}

	private static String string(char[] chars, int[] combination) {
		String rtrn="";
		
		for(int i=0; i<chars.length; i++){
			boolean done=false;
			for(int j=0; j<combination.length; j++){
				if(i==combination[j]){rtrn+="N"; done=true;}
			}
			if(!done){rtrn+=chars[i];}
		}
		
		return rtrn;
	}


	private static Collection<SequenceMotif> makeSequenceMotifWithMismatches(String kmerSequence) throws SearchException {
		Collection<SequenceMotif> rtrn=new ArrayList();
		
		char[] chars=kmerSequence.toCharArray();
		for(int i=0; i<chars.length; i++){
			for(int j=i; j<chars.length; j++){
				if(i!=j){
				String degenerate=string(chars, i, j);
				rtrn.add(new SequenceMotif(degenerate, 1));
				}
			}
		}
		
		return rtrn;
	}
	
	
	private static String string(char[] chars, int p1, int p2) {
		String rtrn="";
		
		for(int i=0; i<chars.length; i++){
			if(i==p1 || i==p2){rtrn+="N";}
			else{rtrn+=chars[i];}
		}
		
		return rtrn;
	}

	private void write(String save,	Map<String, List<SequenceRegion>> seededRegion) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String gene: seededRegion.keySet()){
			writer.write(gene+"\t"+seededRegion.get(gene).size()+"\n");
		}
		
		writer.close();
	}

	private Map<String, List<SequenceRegion>> seed(List<Sequence> geneSequences,	String kmer, int seedNum) throws SearchException {
		List<HairpinKmer> seeds=ComputeOriginalScore.enumerateAllKMers(kmer, seedNum); 
		
		Map<String, List<SequenceRegion>> rtrn=new HashMap();
		
		
		
		double i=0;
		for(HairpinKmer seed: seeds){
			//System.err.print(" "+(i/seeds.size())*100+"...");
			SequenceMotif perfectMotif=new SequenceMotif(seed.getKmerSequence(), 1);
			Map<String, List<SequenceRegion>> match=match(perfectMotif, geneSequences);
			rtrn=addMap(rtrn, match);
			i++;
		}
		
		
		
		return rtrn;
	}
	
	private Map<String, List<SequenceRegion>> addMap(Map<String, List<SequenceRegion>> rtrn, Map<String, List<SequenceRegion>> match) {
		for(String gene: match.keySet()){
			List<SequenceRegion> list=new ArrayList();
			if(rtrn.containsKey(gene)){
				list=rtrn.get(gene);
			}
			list.addAll(match.get(gene));
			rtrn.put(gene, list);
		}
		
		return rtrn;
	}

	private Map<String,List<SequenceRegion>> match(SequenceMotif motif, List<Sequence> geneSequence){
		Map<String, List<SequenceRegion>> rtrn=new TreeMap();
		
		for(Sequence seq: geneSequence) {
			List<SequenceRegion> matches = motif.match(seq);
			if(matches.size()>0){rtrn.put(seq.getId(), matches);}
			
		}
		return rtrn;
	}
	
	
	private static List<Sequence> initializeGeneSequence(String geneSequence) throws IOException {
		FastaSequenceIO fsio = new FastaSequenceIO(geneSequence);
		return fsio.loadAll();
	}
	
	public static void main(String[] args)throws IOException, SearchException{
		if(args.length>4){
			String kmer=args[0];
			List<Sequence> geneSequences=initializeGeneSequence(args[1]);
			int seedNum=new Integer(args[2]);
			String save=args[3];
			int numMismatches=new Integer(args[4]);
			SmatchLike smatch=new SmatchLike(kmer, geneSequences, seedNum);
			smatch.getAllPossibleTargets(numMismatches);
		}
		else{System.err.println(usage);}
	}
	
	
	
	private static void write(String save, Collection<String> possibleTargets,	Collection<String> rcPossibleTargets) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String kmer: possibleTargets){
			writer.write(kmer+"\n");
		}
		
		for(String kmer: rcPossibleTargets){
			writer.write(kmer+"\tRC+\n");
		}
		
		writer.close();
	}

	static String usage=" args[0]=kmer \n args[1]=gene sequence \n args[2]=seed size \n args[3]=save \n args[4]=num mismatches";
	
}
