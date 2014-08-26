package broad.core.primer3;



import java.util.ArrayList;
import java.util.*;
import java.util.Set;
import java.util.TreeSet;
import java.io.*;

public class ComputeMIRScore {
	
	static File mirSeedLookup=new File("/seq/mguttman/RNAi/shRNAs/Designer/MIRSeedLookup.txt");
	static File miRNALookup=new File("/seq/mguttman/RNAi/shRNAs/Designer/HumanMIRNASeeds.txt");
	
	public static double computeMIRScore(String sequence, File mirSeedLookup, File miRNALookup)throws IOException{
		Map<String, Double> mirLookup=parseMirScores(mirSeedLookup);
		Set<String> lookup2=parse(miRNALookup, false);
		Set<String> rcLookup2=parse(miRNALookup, true);
		double mirScore=computeMIRSeedScorePerKMer(sequence, mirLookup, lookup2, rcLookup2);
		return mirScore;
	}
	
	public static double computeMIRScore(String sequence, Map<String, Double> mirLookup, Set<String> lookup2, Set<String> rcLookup2){
		double mirScore=computeMIRSeedScorePerKMer(sequence, mirLookup, lookup2, rcLookup2);
		return mirScore;
	}
	
	public static double computeMIRScore(String sequence)throws IOException{
		return computeMIRScore(sequence, mirSeedLookup, miRNALookup);
	}

	static Map parseMirScores(File file)throws IOException{
		Map rtrn=new TreeMap();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
        String nextLine;
        while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
        	String[] tokens=nextLine.split("\t");
			try{rtrn.put(tokens[0], new Double(tokens[2]));}
			catch(NumberFormatException ex){}	
        }
        reader.close();
        return rtrn;
		
	}
	
	static Map parseMirScores()throws IOException{
		return parseMirScores(mirSeedLookup);
	}
	
	
	
	public static Set parse(File file, boolean rc)throws IOException{
		Set rtrn=new TreeSet();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
        String nextLine;
        while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
        	String[] tokens=nextLine.split("\t");
			if(!rc){rtrn.add(reverseComplement(tokens[3]));}
			else{rtrn.add(tokens[3]);}
        }
        reader.close();
        return rtrn;
		
	}
	
	
	static Set parse(boolean rc)throws IOException{
		return parse(miRNALookup, rc);
	}
	
	private static double computeMIRSeedScorePerKMer(String seq, Map<String, Double> mirLookup, Set lookup2, Set rcLookup2){
		char[] seqChar=seq.toCharArray();
		
		//F11-17
		String f1=seqChar[10]+""+seqChar[11]+""+seqChar[12]+""+seqChar[13]+""+seqChar[14]+""+seqChar[15]+""+seqChar[16];
		
		//F12-18
		String f2=seqChar[11]+""+seqChar[12]+""+seqChar[13]+""+seqChar[14]+""+seqChar[15]+""+seqChar[16]+""+seqChar[17];
				
		//F13-19
		String f3=seqChar[12]+""+seqChar[13]+""+seqChar[14]+""+seqChar[15]+""+seqChar[16]+""+seqChar[17]+""+seqChar[18];
		
		double maxSenseScore=-Double.MAX_VALUE;
		try{	
			maxSenseScore=mirLookup.get(f1);
			maxSenseScore=Math.max(maxSenseScore, mirLookup.get(f2));
			maxSenseScore=Math.max(maxSenseScore, mirLookup.get(f3));
		}catch(NullPointerException ex){maxSenseScore=100;}
		
		if(lookup2.contains(f1)){maxSenseScore=100;}
		if(lookup2.contains(f2)){maxSenseScore=100;}
		if(lookup2.contains(f3)){maxSenseScore=100;}
		
		String reverseSeq=reverseComplement(seq);
		seqChar=reverseSeq.toCharArray();
		
		
		//R14-20
		String r1=seqChar[13]+""+seqChar[14]+""+seqChar[15]+""+seqChar[16]+""+seqChar[17]+""+seqChar[18]+""+seqChar[19];
		
		//R15-21
		String r2=seqChar[14]+""+seqChar[15]+""+seqChar[16]+""+seqChar[17]+""+seqChar[18]+""+seqChar[19]+""+seqChar[20];
		
		//R16-21 + 3' C 
		String r3=seqChar[15]+""+seqChar[16]+""+seqChar[17]+""+seqChar[18]+""+seqChar[19]+""+seqChar[20]+"C";
			
		double maxASenseScore=-Double.MAX_VALUE;
		try{
			maxASenseScore=mirLookup.get(r1);
			maxASenseScore=Math.max(maxASenseScore, mirLookup.get(r2));
			maxASenseScore=Math.max(maxASenseScore, mirLookup.get(r3));
		}catch(NullPointerException ex){maxASenseScore=100;}
		
		if(lookup2.contains(r1)){maxASenseScore=100;}
		if(lookup2.contains(r2)){maxASenseScore=100;}
		if(lookup2.contains(r3)){maxASenseScore=100;}
		
		//System.err.println(seq+" "+maxSenseScore);
		//System.err.println(reverseSeq+" "+maxASenseScore);
		
		double rtrn=forwardScore(maxSenseScore)* reverseScore(maxASenseScore);
		return rtrn;
	}
	
	private static double forwardScore(double score){
		if(score<=.35){return 1.25;}
		if(score<=1.2){return 1.2;}
		if(score<=2.0){return .85;}
		if(score<=3){return .7;}
		return .2;
	}
	
	private static double reverseScore(double score){
		if(score<=.35){return 1.1;}
		if(score<=1.2){return 1;}
		if(score<=2.0){return .9;}
		if(score<=3){return .8;}
		return .7;
	}
	
	
	public static String reverseComplement(String seq){
		char[] seqChar=seq.toCharArray();
		char[] reverse=new char[seqChar.length];
		
		int j=0;
		for(int i=seqChar.length-1; i>=0; i--){
			reverse[j]=seqChar[i];
			j++;
		}
		
		String rtrn="";
		for(int i=0; i<reverse.length; i++){
			if(reverse[i]=='A' || reverse[i]=='a'){rtrn+='T';}
			if(reverse[i]=='T' || reverse[i]=='t'){rtrn+='A';}
			if(reverse[i]=='C' || reverse[i]=='c'){rtrn+='G';}
			if(reverse[i]=='G' || reverse[i]=='g'){rtrn+='C';}
			if(reverse[i]=='N' || reverse[i]=='n'){rtrn+='N';}
		}
		
		return rtrn;
	}
	
	public static String reverse(String seq){
		char[] seqChar=seq.toCharArray();
		char[] reverse=new char[seqChar.length];
		
		int j=0;
		for(int i=seqChar.length-1; i>=0; i--){
			reverse[j]=seqChar[i];
			j++;
		}
		
		String rtrn="";
		for(int i=0; i<reverse.length; i++){
			rtrn+=reverse[i];
		}
		
		return rtrn;
	}
	
	public static String complement(String seq){
		char[] seqChar=seq.toCharArray();
				
		String rtrn="";
		for(int i=0; i<seqChar.length; i++){
			if(seqChar[i]=='A' || seqChar[i]=='a'){rtrn+='T';}
			if(seqChar[i]=='T' || seqChar[i]=='t'){rtrn+='A';}
			if(seqChar[i]=='C' || seqChar[i]=='c'){rtrn+='G';}
			if(seqChar[i]=='G' || seqChar[i]=='g'){rtrn+='C';}
			if(seqChar[i]=='N' || seqChar[i]=='n'){rtrn+='N';}
		}
		
		return rtrn;
	}
	
	
	public static void main(String[] args)throws IOException{
		String seq=args[0];
		double mirScore=ComputeMIRScore.computeMIRScore(seq);
		//System.err.println(mirScore);
	}

	
}
