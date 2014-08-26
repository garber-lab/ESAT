package broad.pda.seq.fastq;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;

import broad.core.datastructures.Pair;
import broad.core.sequence.Sequence;
import broad.pda.annotation.BEDFileParser;
import broad.pda.seq.fastq.SplitFastqByIndex;

//For use with internal barcodes

public class SplitAndTrimIndex {

	public SplitAndTrimIndex(File leftFastqFile, File rightFastqFile, File barcodeInfoFile, String saveDir) throws IOException{
		Map<String, String> barcodeInfo=parseBarcodeInfo(barcodeInfoFile);
		
		Map<String, Integer> barcodeCount = trimReads(leftFastqFile, rightFastqFile, barcodeInfo, saveDir);
		SplitFastqByIndex.writeBarcodeCount(saveDir+"/barcodeCount.txt", barcodeCount, true);	
	}

	
	
	/**
	 * Create map associating reverse of barcode with sample ID
	 * @param barcodeInfoFile flat file with barcode sequence in first column, sample name in second column, headers in first line
	 * @return map with key set the RC'ed barcodes, value the sample ID
	 * @throws IOException
	 */
	private Map<String, String> parseBarcodeInfo(File barcodeInfoFile) throws IOException {
		Map<String, String> rtrn=new TreeMap<String, String>();
		Collection<String> lines=BEDFileParser.loadList(barcodeInfoFile.getAbsolutePath(), true);
		
		for(String line: lines){
			String[] tokens=line.split("\t");
			String seq=Sequence.reverseSequence(tokens[0]);
			rtrn.put(seq, tokens[1]);
		}
		
		return rtrn;
	}



	/**
	 * 
	 * @param file1 read 1 fastq file
	 * @param file2 read 2 fastq file
	 * @param barcodeInfo barcde info map returned by parseBarcodeInfo
	 * @param saveDir directory to write to
	 * @throws IOException
	 */
	private Map<String, Integer> trimReads(File file1, File file2, Map<String, String> barcodeInfo, String saveDir) throws IOException {
		Map<String, Pair<FileWriter>> writers=new TreeMap<String, Pair<FileWriter>>();
		Map<String, Integer> barcodeCounter=new HashMap<String, Integer>();
		
		for(String barcode: barcodeInfo.keySet()){
			String name=barcodeInfo.get(barcode);
			FileWriter writer1=new FileWriter(saveDir+"/"+name+"_1.fq");
			FileWriter writer2=new FileWriter(saveDir+"/"+name+"_2.fq");
			Pair<FileWriter> pair=new Pair<FileWriter>(writer1, writer2);
			writers.put(barcode, pair);
		}
		writers.put("remainder", new Pair<FileWriter>(new FileWriter(saveDir+"/remainder_1.fq"), new FileWriter(saveDir+"/remainder_2.fq")));
		
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file1)));
		BufferedReader readerIndex=new BufferedReader(new InputStreamReader(new FileInputStream(file2)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
	    	if(nextLine.startsWith("@")){
	    		FastqSequence seq=new FastqSequence(nextLine, reader.readLine(), reader.readLine(), reader.readLine());
	        	String barcode=seq.getFirstNBPs(8);
	        	FastqSequence trimmed=seq.trimFirstNBPs(9);
	    		Pair<FileWriter> writer=writers.get("remainder");
	    		if(writers.containsKey(barcode)){
	    			writer=writers.get(barcode);
	    		}
	    		
	    		FastqSequence index=new FastqSequence(readerIndex.readLine(), readerIndex.readLine(), readerIndex.readLine(), readerIndex.readLine());
	        	if(!seq.getName().split(" ")[0].equalsIgnoreCase(index.getName().split(" ")[0])){
	        		throw new IllegalArgumentException("Read1 and Read2 dont match");
	        	}
	        	
	    		writer.getValue1().write(trimmed+"\n");
	    		writer.getValue2().write(index+"\n");
	    		
	    		int count = 0;
	            if (barcodeCounter.containsKey(barcode)) {
	            	count = barcodeCounter.get(barcode);
	            }
	            count++;
	            barcodeCounter.put(barcode, count);
	    	}
	    }
		
		for(String barcode: writers.keySet()){
			Pair<FileWriter> writer=writers.get(barcode);
			writer.getValue1().close();
			writer.getValue2().close();
		}
		
		return barcodeCounter;
	}
	
	/**
	 * 
	 * @param args read 1 file, read 2 file, barcode info file, directory to write to
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException{
		if(args.length>3){
			File read1=new File(args[0]);
			File read2=new File(args[1]);
			File barcodeInfo=new File(args[2]);
			String saveDir=args[3];
			new SplitAndTrimIndex(read1, read2, barcodeInfo, saveDir);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=Read 1 \n args[1]=Read 2 \n args[2]=barcode info with headers in first line \n args[3]=save directory";
	
}
