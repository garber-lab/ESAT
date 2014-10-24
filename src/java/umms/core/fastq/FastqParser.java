package umms.core.fastq;

import java.io.*;
import java.util.*;

import org.apache.log4j.Logger;


public class FastqParser implements Iterator<FastqSequence>{
	static Logger logger = Logger.getLogger(FastqParser.class.getName());
	Collection<FastqSequence> sequences;
	File fastqFile;
	BufferedReader reader;
	private int numberOfSeq;
	String nextLine = null;
		
	public FastqParser() {
		super();
	}
	
	public void start(File fastqFile) throws IOException {
		this.fastqFile = fastqFile;
		logger.debug(fastqFile.exists() + fastqFile.getAbsolutePath() );
		reader=new BufferedReader(new FileReader(fastqFile));
		nextLine = reader.readLine();
	}
	
	public void start (BufferedReader br) throws IOException {
		reader=br;
		nextLine = reader.readLine();
	}
	
	public void convertToFasta(String save)throws IOException{
		FileWriter writer=new FileWriter(save);
	
		String nextLine;
        while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
        	if(nextLine.startsWith("@")){
        		String firstLine=nextLine;
        		String secondLine=reader.readLine();
        		String thirdLine=reader.readLine();
        		String fourthLine=reader.readLine();
        		
        		FastqSequence seq=new FastqSequence(firstLine, secondLine, thirdLine, fourthLine);
    			writer.write(seq.toFasta());
        	}
        }
        writer.close();
	}
	
	public void convertToNumberedFasta(String save)throws IOException{
		FileWriter writer=new FileWriter(save);
	
		String nextLine;
    	int i=1;
        while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
        	if(nextLine.startsWith("@")){
        		String firstLine=nextLine;
        		String secondLine=reader.readLine();
        		String thirdLine=reader.readLine();
        		String fourthLine=reader.readLine();
        		
        		FastqSequence seq=new FastqSequence(firstLine, secondLine, thirdLine, fourthLine);
    			writer.write(seq.toFasta(i));
    			i++;
        	}
        }
        writer.close();
	}
	
	public Collection<FastqSequence> parse(File file)throws IOException{
		Collection<FastqSequence> rtrn=new ArrayList<FastqSequence>();
		String nextLine;
        while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
        	if(nextLine.startsWith("@")){
        		String firstLine=nextLine;
        		String secondLine=reader.readLine();
        		String thirdLine=reader.readLine();
        		String fourthLine=reader.readLine();
        		
        		FastqSequence seq=new FastqSequence(firstLine, secondLine,thirdLine, fourthLine);
  			
        		rtrn.add(seq);
        	}
        	
        	
        }
        return rtrn;
	}
	
		
	public Collection<FastqSequence> getSequences() throws IOException{
		if(this.sequences!=null){
		return this.sequences;
		}
		else{this.sequences=this.parse(fastqFile); return this.sequences;}
	}
	
	public int getNumberOfSequences(){
		if(this.numberOfSeq>0){return this.numberOfSeq;}
		else{
			try{
			BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(this.fastqFile)));
			String nextLine;
			int i=0;
	        while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
	        	if(nextLine.startsWith("@")){i++;}
	        }
	        this.numberOfSeq=i;
	        reader.close();
			}catch(IOException ex){}
			return this.numberOfSeq;
		}
	}

	public boolean hasNext() {
		return nextLine != null;
	}

	public FastqSequence next() {
		FastqSequence seq = null;
		try{

			
	        if (nextLine  != null) {
        		String firstLine=nextLine;
        		String secondLine=reader.readLine();
				String thirdLine=reader.readLine();
        		String fourthLine=reader.readLine();
        		seq=new FastqSequence(firstLine, secondLine, thirdLine, fourthLine);
	        }
	        nextLine = reader.readLine() ;

		}catch(Exception ex){ 
			logger.error("Exception thrown while reading fastq file",ex);
		}

		return seq;
	}

	public void close() throws IOException{
		reader.close();
	}

	public void remove() {}
	
}
