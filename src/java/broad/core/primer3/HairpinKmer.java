package broad.core.primer3;


import java.io.IOException;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import broad.core.motif.SearchException;
import broad.core.motif.SequenceMotif;
import broad.core.sequence.Sequence;
import broad.core.sequence.SequenceRegion;
import broad.core.sequence.WindowSlider;
import broad.pda.rnai.RNAiGeneAnnotation;

public class HairpinKmer implements Comparable<HairpinKmer>{

	String kmer;
	int relativeStart;
	int relativeEnd;
	double mirScore;
	double originalScore;
	double rs8Score;
	boolean blastHit;
	double percentVariantsCovered;
	Collection<RNAiGeneAnnotation> transcripts;
	
	public HairpinKmer(String kmer, int relativeStart, int relativeEnd){
		this.kmer=kmer.toUpperCase();
		this.relativeStart=relativeStart;
		this.relativeEnd=relativeEnd;
		this.rs8Score=-999;
		this.originalScore=-999;
		this.mirScore=-999;
		this.blastHit=false;
		this.percentVariantsCovered=0;
		this.transcripts=new HashSet();
	}
	
	public HairpinKmer(String kmer, int relativeStart, int relativeEnd, Map<String, Double> mirLookup, Set<String> lookup2, Set<String> rcLookup2){
		this.kmer=kmer.toUpperCase();
		this.relativeStart=relativeStart;
		this.relativeEnd=relativeEnd;
		this.computeScore(mirLookup, lookup2, rcLookup2);
		this.blastHit=false;
		this.percentVariantsCovered=0;
		this.transcripts=new HashSet();
	}
	
	public void addTranscript(RNAiGeneAnnotation transcript){this.transcripts.add(transcript);}
	public Collection<RNAiGeneAnnotation> getTranscripts(){return this.transcripts;}
	
	public void setPercentVariantsCovered(double percent){this.percentVariantsCovered=percent;}
	public double getPercentVariantsCovered(){return this.percentVariantsCovered;}
	
	private void computeScore(Map<String, Double> mirLookup, Set<String> lookup2, Set<String> rcLookup2) {
		this.originalScore=ComputeOriginalScore.computeOriginalScorePerKMer(kmer)[0];
		this.mirScore=ComputeMIRScore.computeMIRScore(kmer, mirLookup, lookup2, rcLookup2);
		this.rs8Score=originalScore*mirScore;
		
	}

	private void computeScore() throws IOException{
		this.originalScore=ComputeOriginalScore.computeOriginalScorePerKMer(kmer)[0];
		this.mirScore=ComputeMIRScore.computeMIRScore(kmer);
		this.rs8Score=originalScore*mirScore;
	}
		
	public String getKmerSequence(){return this.kmer;}
	public int getKmerStartPosition(){return this.relativeStart;}
	public int getKmerEndPosition(){return this.relativeEnd;}



	
	public double getRS8Score(){
		if(this.rs8Score==-999){try {
			this.computeScore();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}}
		return this.rs8Score;
	}
	public double getOriginalScore(){
		if(this.originalScore==-999){try {
			this.computeScore();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}}
		return this.originalScore;
	}
		
	public double getMirScore(){
		if(this.mirScore==-999){try {
			this.computeScore();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}}
		return this.mirScore;
	}

	public int compareTo(HairpinKmer kmer) {
		
		//Compare percent of variants covered
		if(kmer.getPercentVariantsCovered()!=getPercentVariantsCovered()){return new Double(kmer.getPercentVariantsCovered()).compareTo(getPercentVariantsCovered());}
		
		//Compare scores
		if(kmer.getRS8Score()!=getRS8Score()){return new Double(kmer.getRS8Score()).compareTo(getRS8Score());}
		
		//first compare starts
		if(kmer.getKmerStartPosition()!=getKmerStartPosition()){return getKmerStartPosition()-kmer.getKmerStartPosition();}
		//if they are the same compare ends
		else if(kmer.getKmerEndPosition()!=getKmerEndPosition()){return getKmerEndPosition()-kmer.getKmerEndPosition();}
		//else compare sequence
		else{return getKmerSequence().compareTo(kmer.getKmerSequence());}
	}

	public boolean setRS8Score(double score){
		if(this.rs8Score==-999){this.rs8Score=score; return true;}
		return false;
	}
	
	public boolean setOriginalScore(double score){
		if(this.originalScore==-999){this.originalScore=score; return true;}
		return false;
	}
	
	public boolean setMirScore(double score){
		if(this.mirScore==-999){this.mirScore=score; return true;}
		return false;
	}
	
	public void setBlastHit(boolean b) {
		this.blastHit=b;
	}
	
	public boolean getBlastHit(){return this.blastHit;}

	public List<SequenceRegion> getRelativePositions(RNAiGeneAnnotation transcript) throws SearchException {
		SequenceMotif motif=new SequenceMotif(this.kmer,1);
		return motif.match(transcript.getSequence());
	}

	//return the overlap in bps between 2 kmers
	public int getOverlap(HairpinKmer hp) {
		int startOverlap=Math.max(hp.getKmerStartPosition(), getKmerStartPosition());
		int endOverlap=Math.min(hp.getKmerEndPosition(), getKmerEndPosition());
		return endOverlap-startOverlap;
	}

	public double getPercentGC() {
		return ComputeOriginalScore.computeGCPercent(kmer);
	}

	public int maxSelfComplementarity() {
		int starting=this.kmer.length()/2;
		for(int i=starting; i>=1; i--){
			Sequence kmerseq=new Sequence("kmer");
			kmerseq.setSequenceBases(kmer);
			WindowSlider slider=WindowSlider.getSlider(kmerseq, i, i-1);
			while(slider.hasNext()){
				SequenceRegion region1=slider.next();
				String seq=region1.getSequenceBases();
				WindowSlider slider2=WindowSlider.getSlider(kmerseq, i, i-1);
				while(slider2.hasNext()){
					SequenceRegion region2=slider2.next();
					String seq2=region2.getSequenceBases();
					if(region1.getStart()!=region2.getStart() && region2.getEnd()!=region1.getEnd()){
					if(seq.equalsIgnoreCase(Sequence.reverseSequence(seq2)) || seq.equalsIgnoreCase(Sequence.complement(seq2))){
						System.err.println(i+" "+seq+" "+seq2);
						return i;
					}
					}
				}
			}
		}
		return 0;
	}
	
}
