package broad.pda.seq.fastq;

import java.io.BufferedWriter;
import java.io.IOException;

public class FastqSequence {

	private String sequence;
	private String quality;
	String name;
	String description;
	
	public FastqSequence(String name, String sequence, String description, String quality){
		this.sequence=sequence;
		this.quality=quality;
		this.name=name;
		this.description=description;
	}
	
	public String getName(){return this.name;}
	public String getSequence(){return this.sequence;}
	public String getDescription() {return this.description;}
	
	public void setName(String name) {this.name = name;}
	public void setDescription(String descr) {this.description = descr;}
	
	public String toString(){
		String rtrn=name+"\n"+sequence+"\n"+description+"\n"+quality;
		return rtrn;
	}
	
	public void write(BufferedWriter bw) throws IOException{
		bw.write( (name.startsWith("@") ) ? name : "@" + name );  // TODO: Fix, store name without @, add @ only for writing
		bw.newLine();
		bw.write(sequence);
		bw.newLine();
		bw.write(description.startsWith("+") ? description : "+" + description); // TODO: Fix, store description without +, add + only for writing
		bw.newLine();
		bw.write(quality);
		bw.newLine();
	}
	
	public void removeAtSymbolFromName() {
		if(name.charAt(0) == '@') {
			name = name.substring(1);
		}
	}
	
	public FastqSequence trimEnds(char letter){
		//get the last occurrence of the specified letter starting from the end and remove it from the sequence
		char[] chars=sequence.toCharArray();
		int endIndex=chars.length-1;
		for(int i=chars.length-1; i>=0; i--){
			if(chars[i]!=letter){endIndex=i; break;}
		}
		String newSequence=sequence.substring(0, endIndex+1);
		String newQuality=quality.substring(0,endIndex+1);
		FastqSequence newSeq=new FastqSequence(name, newSequence, description, newQuality);
		return newSeq;
	}
	
	/**
	 * Trims the last numBasesToTrim bases of the fastq record and returns the nucleotides trimmed.
	 * @param numOfBasesToTrim
	 * @return
	 */
	public FastqSequence trimEndBases(int numOfBasesToTrim){
		String trimmedSequence = sequence.substring(sequence.length() - numOfBasesToTrim);
		String trimmedQual     = quality.substring(quality.length() - numOfBasesToTrim);
		sequence=sequence.substring(0, sequence.length() - numOfBasesToTrim);
		quality=quality.substring(0, quality.length() - numOfBasesToTrim);
		return new FastqSequence(name, trimmedSequence, description, trimmedQual);
	}
	
	/**
	 * Trims the first numBasesToTrim bases of the fastq record and returns the nucleotides trimmed.
	 * @param numOfBasesToTrim
	 * @return
	 */
	public FastqSequence trimStartBases(int numOfBasesToTrim){
		String trimmedSequence = sequence.substring(0, numOfBasesToTrim);
		String trimmedQual     = quality.substring(0, numOfBasesToTrim);
		sequence=sequence.substring(numOfBasesToTrim);
		quality=quality.substring(numOfBasesToTrim);
		return new FastqSequence(name, trimmedSequence, description, trimmedQual);
	}
	
	
	public FastqSequence trimBeginning(char letter){
		//get the last occurrence of the specified letter starting from the end and remove it from the sequence
		char[] chars=sequence.toCharArray();
		int startIndex=0;
		for(int i=0; i<chars.length; i++){
			if(chars[i]!=letter){startIndex=i; break;}
		}
		String newSequence=sequence.substring(startIndex, chars.length);
		String newQuality=quality.substring(startIndex,chars.length);
		FastqSequence newSeq=new FastqSequence(name, newSequence, description, newQuality);
		return newSeq;
	}
	
	public FastqSequence trimFirstNBPs(int n){
		char[] chars=sequence.toCharArray();
		int startIndex=n;
		
		String newSequence=sequence.substring(startIndex, chars.length);
		String newQuality=quality.substring(startIndex,chars.length);
		FastqSequence newSeq=new FastqSequence(name, newSequence, description, newQuality);
		return newSeq;
	}
	
	public String getFirstNBPs(int n){
		return sequence.substring(0, n);
	}
	
	public int getLength(){return sequence.toCharArray().length;}

	public String toFasta() {
		String rtrn="";
		rtrn+=">"+name+"\n"+sequence+"\n";
		return rtrn;
	}

	public String toFasta(int readNum) {
		String rtrn="";
		String name="seq."+readNum+"a";
		rtrn+=">"+name+"\n"+sequence+"\n";
		return rtrn;
	}

	public String toFastq() {
		return this.toString();
	}

	private boolean isPolyA(String seq) {
		//if sequence is a homopolymer of A or T
		char[] bases=seq.toCharArray();
		int aCount=0;
		int tCount=0;
		for(int i=0; i<bases.length; i++){
			if(bases[i]=='A' || bases[i]=='a'){aCount++;}
			if(bases[i]=='T' || bases[i]=='t'){tCount++;}
		}
		if(aCount==bases.length || tCount==bases.length){return true;}
		return false;
	}

	public boolean isPolyA(){return isPolyA(this.sequence);}
	
	public boolean isPartialPolyA(int polyN) {
		String lastNBps=getLastBps(sequence, polyN);
		String firstNBps=getFirstBps(sequence, polyN);
		
		if(isPolyA(lastNBps) || isPolyA(firstNBps)){return true;}
		return false;
	}
	
	private String getLastBps(String sequence, int num){
		return sequence.substring(sequence.toCharArray().length-num);
	}
	
	private String getFirstBps(String sequence, int num){
		return sequence.substring(0, num);
	}

	public FastqSequence trimPolyA() {
		FastqSequence s1=this.trimBeginning('T');
		FastqSequence s2=this.trimBeginning('A');
		FastqSequence s3=this.trimEnds('A');
		FastqSequence s4=this.trimEnds('T');
		
		int minSize=s1.getLength();
		minSize=Math.min(minSize, s2.getLength());
		minSize=Math.min(minSize, s3.getLength());
		minSize=Math.min(minSize, s4.getLength());
		
		if(s1.getLength()==minSize){return s1;}
		if(s2.getLength()==minSize){return s2;}
		if(s3.getLength()==minSize){return s3;}
		if(s4.getLength()==minSize){return s4;}
		return s1;
	}

	public FastqSequence trimLinker(String linkerSeq) {
		for(int i=0; i<linkerSeq.length(); i++){
			String linker=linkerSeq.substring(i, linkerSeq.length());
			String fivePrime=this.sequence.substring(0, linkerSeq.length()-i);
			if(fivePrime.equalsIgnoreCase(linker)){
				String seq=this.sequence.substring((linkerSeq.length()-i), this.sequence.length());
				String qual=this.quality.substring((linkerSeq.length()-i), this.sequence.length());
				//System.err.println(linkerSeq.length()-i +" match "+fivePrime+" "+linker+" "+this.sequence+" TRIMMED "+seq);
				FastqSequence fastq=new FastqSequence(name, seq, description, qual);
				return fastq;
			}
		}
		return this;
	}

	
}
