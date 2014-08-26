package broad.core.primer3;


import java.io.*;
import java.util.*;

//TODO: Make sure the original score is case insensitive
public class ComputeOriginalScore {
	
	//add option to not target UTR
	//static int k=21;
	
	public static Collection<HairpinKmer> computeOriginalScore(String seq, int k){
		//main
		Collection<HairpinKmer> kmers=enumerateAllKMers(seq, k);
		//Collection<HairpinKmer> kmerScores=computeOriginalScoreForAllKMers(kmers);
		return kmers;
	}
	
	/*private static Collection<HairpinKmer> computeOriginalScoreForAllKMers(ArrayList<HairpinKmer> kmers)throws IOException{
		Map<HairpinKmer, double[]> rtrn=new TreeMap<HairpinKmer, double[]>();
		for(HairpinKmer kmer: kmers){
			//System.err.println(kmer+" "+kmer.toCharArray().length);
			double[] scores=computeOriginalScorePerKMer(kmer.getKmerSequence());
			double mirScore=ComputeMIRScore.ComputeMIRScore(kmer.getKmerSequence());
			//System.err.println(kmer+" "+score);
			double newScore=scores[0]*mirScore;
			double[] array={newScore, scores[0], mirScore};
			kmer.setScores(array);
			rtrn.add(kmer);
		}
		return rtrn;
	}*/
	
	
	
	//enumerate all kmers
	public static List<HairpinKmer> enumerateAllKMers(String sequence, int k){
		List<HairpinKmer> rtrn=new ArrayList<HairpinKmer>();
		char[] seqChar=sequence.toCharArray();
		for(int i=0; i<=seqChar.length-k; i++){
			String kmer="";
			for(int j=0; j<k; j++){
				kmer+=seqChar[i+j];
			}
			HairpinKmer hp=new HairpinKmer(kmer, i, i+k);
			rtrn.add(hp);
		}
		return rtrn;
	}
	
	public static double[] computeOriginalScorePerKMer(String kmer){
		double clampScore=computeClampScore(kmer);
		double GCScore=computeGCScore(kmer);
		double fourInARow=penalty4Mer(kmer);
		double sevenGCsInARow=penalty7GCs(kmer);
		double AAStart=computeAAStart(kmer);
		double maskScore=penalizeMaskedRegions(kmer);
		
		double score=clampScore*GCScore*fourInARow*sevenGCsInARow*AAStart*maskScore;
		
		//original score, clamp score, GCScore, fourInARow, 7GCs, AAStart, maskScore
		double[] originalScores={score, clampScore, GCScore, fourInARow, sevenGCsInARow, AAStart, maskScore};
		return originalScores;
	}
	
	
	
	private static double penalizeMaskedRegions(String kmer){
		char[] seqChar=kmer.toCharArray();
		int count=0;
		for(int i=0; i<seqChar.length; i++){
			if(seqChar[i]!='A' && seqChar[i]!='C' && seqChar[i]!='G' && seqChar[i]!='T'){count++;}
		}
		if(count>0){return (0.00000000000000001);}
		else{return 1;}
	}
	
	private static double computeAAStart(String kmer){
		/*pos0/1=AA -> 0.000000000000001*/
		
		char[] seqChar=kmer.toCharArray();
		if(seqChar[0]=='A' && seqChar[1]=='A'){return 0.000000000000001;}
		return 1;
	}
	
	private static double penalty7GCs(String kmer){
		/*Look for 7 GCs in a row*/
		List<HairpinKmer> subKmers=enumerateAllKMers(kmer, 7);
		int GC7MersCount=0;
		for(HairpinKmer submer: subKmers){
			double GCPercent=computeGCPercent(submer.getKmerSequence());
			if(GCPercent==1){GC7MersCount++;}
		}
		//if(GC7MersCount>0){return .01;}
		if(GC7MersCount>0){return .01;}
		else{return 1;}
	}
	
	private static double penalty4Mer(String kmer){
		/*Look for 4 repeated nt in a row*/
		List<HairpinKmer> subKmers=enumerateAllKMers(kmer, 4);
		int penaltyKmerCount=0;
		for(HairpinKmer hp: subKmers){
			String submer=hp.getKmerSequence();
			if(submer.equalsIgnoreCase("AAAA")){penaltyKmerCount++;}
			else if(submer.equalsIgnoreCase("CCCC")){penaltyKmerCount++;}
			else if(submer.equalsIgnoreCase("GGGG")){penaltyKmerCount++;}
			else if(submer.equalsIgnoreCase("TTTT")){penaltyKmerCount++;}
		}
		
		//if(penaltyKmerCount>0){return .01;}
		if(penaltyKmerCount>0){return 0.01;}
		else{return 1;}
	}
	
	private static double computeGCScore(String kmer){
		/*GC Score
		Rules:
		%GC<=25% -> 0.01
		25%<%GC<=55% -> 3
		55%<%GC<=60% -> 1
		60%<=%GC -> 0.01
		*/
		
		double percentGC=computeGCPercent(kmer);
		if(percentGC<=.25){return .01;}
		else if(percentGC<=.55){return 3;}
		else if(percentGC<=.6){return 1;}
		else{return .01;}
		
	}
	
	public static double computeGCPercent(String kmer){
		char[] seqChar=kmer.toCharArray();
		double GCCount=0;
		for(int i=0; i<seqChar.length; i++){
			if(seqChar[i]=='G' || seqChar[i]=='C' || seqChar[i]=='g' || seqChar[i]=='c'){GCCount++;}
		}
		double percentGC=GCCount/seqChar.length;
		return percentGC;
	}
	
	private static double computeClampScore(String kmer){
		/*3' Clamp Score
		Rules:
		pos17/18/19=3A/T -> 4
		pos17/18/19=2A/T -> 1.5
		pos17/18/19=1A/T -> 0.8
		pos17/18/19=0A/T -> 0.2 */
		
		char[] seqChar=kmer.toCharArray();
		int clampCount=0;
		if(seqChar[17]=='A' || seqChar[17]=='T'){clampCount++;}
		if(seqChar[18]=='A' || seqChar[18]=='T'){clampCount++;}
		if(seqChar[19]=='A' || seqChar[19]=='T'){clampCount++;}
		
		double clampScore=0;
		if(clampCount==3){clampScore=4;}
		else if(clampCount==2){clampScore=1.5;}
		else if(clampCount==1){clampScore=.8;}
		else{clampScore=.2;}
		
		/*
		Multiplied by:
		pos15/16/20=3A/T -> 1.25
		pos15/16/20=2A/T -> 1.1
		pos15/16/20=1A/T -> 0.9
		pos15/16/20=0A/T -> 0.5*/

		int clampCountMult=0;
		double clampScoreMult=0;
		
		if(seqChar[15]=='A' || seqChar[15]=='T'){clampCountMult++;}
		if(seqChar[16]=='A' || seqChar[16]=='T'){clampCountMult++;}
		if(seqChar[20]=='A' || seqChar[20]=='T'){clampCountMult++;}
		
		if(clampCountMult==3){clampScoreMult=1.25;}
		else if(clampCountMult==2){clampScoreMult=1.1;}
		else if(clampCountMult==1){clampScoreMult=.9;}
		else{clampScoreMult=.5;}
		
		clampScore=clampScore*clampScoreMult;
		return clampScore;
	}

	public static List<HairpinKmer> enumerateAllKMers(String sequence, int k, Map<String, Double> mirLookup, Set<String> lookup2, Set<String> rcLookup2) {
		List<HairpinKmer> rtrn=new ArrayList<HairpinKmer>();
		char[] seqChar=sequence.toCharArray();
		for(int i=0; i<=seqChar.length-k; i++){
			String kmer="";
			for(int j=0; j<k; j++){
				kmer+=seqChar[i+j];
			}
			HairpinKmer hp=new HairpinKmer(kmer, i, i+k, mirLookup, lookup2, rcLookup2);
			rtrn.add(hp);
		}
		return rtrn;
	}
		
}
