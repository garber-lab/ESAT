package broad.core.primer3;

import broad.core.annotation.MaximumContiguousSubsequence;
import jaligner.Alignment;

public class PrimerUtils {
	
	public static double computeMaxTM(Alignment align){
		return Math.max(computeTM(align), computeTM2(align));
	}

	public static double computeTM(String seq){
		//Tm = 64.9C + 41C x (number of Gs and Cs in the primer  16.4)/N
		char[] bases=seq.toCharArray();
		double TM=64.9;
		int numGC=0;
		for(int i=0; i<bases.length; i++){
			if(bases[i]=='G'|| bases[i]=='C' || bases[i]=='c' || bases[i]=='g'){numGC++;}
		}
		TM=TM+41*((numGC-16.4)/(double)bases.length);
		
		return TM;
	}
	
	public static void main(String[] args){
		String seq="TAATACGACTCACTATAGGG";
		System.err.println(seq+" "+computeTM(seq));
	}

	public static double percentGC(String seq) {
		char[] bases=seq.toCharArray();
		int numGC=0;
		for(int i=0; i<bases.length; i++){
			if(bases[i]=='G'|| bases[i]=='C' || bases[i]=='c' || bases[i]=='g'){numGC++;}
		}
				
		return 100.0*(numGC/(double)bases.length);
	}

	public static double computeTM(Alignment align) {
		//COmpute the TM only across the matches
		String concat="";
		//Tm = 64.9C + 41C x (number of Gs and Cs in the primer  16.4)/N
		double TM=64.9;
		int numGC=0;
		int numMismatch=0;
		double length=0.0;
		for(int i=0; i<align.getSequence1().length; i++){
			char base1=align.getSequence1()[i];
			char base2=align.getSequence2()[i];
			if(base1==base2 && base1!='-'){
				if(base1=='G'|| base1=='C' || base1=='c' || base1=='g'){numGC++;}
				concat=concat+base1;
				length++;
			}
			else if(base1!=base2 && base1!='-' && base2!='-'){
				numMismatch++;
			}
		}
		
		TM=TM+41*((numGC-16.4)/length);
		double percentMismatch=((numMismatch/(length+numMismatch)))*100;
		double TM1=Math.max(0, TM-percentMismatch);
		//System.err.println(TM+" "+TM1+" "+numMismatch+" "+length+" "+percentMismatch);
		return TM1;
	}
	
	public static double computeTM2(Alignment align) {
		//COmpute the TM only across the matches
		String concat="";
		//Tm = 64.9C + 41C x (number of Gs and Cs in the primer  16.4)/N
		double TM=64.9;
		int numGC=0;
		int numMismatch=0;
		double length=0.0;
		int[] vals=new int[align.getSequence1().length];
		for(int i=0; i<align.getSequence1().length; i++){
			char base1=align.getSequence1()[i];
			char base2=align.getSequence2()[i];
			if(base1==base2 && base1!='-'){
				if(base1=='G'|| base1=='C' || base1=='c' || base1=='g'){numGC++;}
				concat=concat+base1;
				length++;
				vals[i]=1;
			}
			else if(base1!=base2 && base1!='-' && base2!='-'){
				numMismatch++;
				vals[i]=-1000;
			}
			else{vals[i]=-2000;}
		}
		

		int[] maxSub=MaximumContiguousSubsequence.maxSubSum3(vals);
		
		
		if(maxSub[0]>0){
			char[] sub=new char[maxSub[2]-maxSub[1]];
			int counter=0;
			for(int i=maxSub[1]; i<maxSub[2]; i++){
				sub[counter++]=align.getSequence1()[i];
			}
			if(sub.length>10){
				double subTM=computeTM(new String(sub));
				return subTM;
			}
			else{return 0;}
		}
		
		return 0;
	}


}
