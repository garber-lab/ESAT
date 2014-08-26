package broad.pda.seq.fastq;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.annotation.BEDFileParser;
//import broad.pda.rap.CountProbeKmers;

import net.sf.picard.fastq.*;


public class SplitFastqByIndex {

	public SplitFastqByIndex(String fastq, String indexInfoFname, String index, String saveDir, boolean allowMismatch) throws IOException {
		Map<String, FastqWriter> indexInfo=parseIndexInfo(indexInfoFname, saveDir, allowMismatch);
		
		FastqReader fq = new FastqReader(new File(fastq));
		
		Map<String, Integer> barcodeCount;
		if (index == null) {
			barcodeCount = readAndWriteByBarcode(fq, indexInfo, allowMismatch);
		} else {
			FastqReader indexReader = new FastqReader(new File(index));
			barcodeCount = readAndWriteByBarcodeSeparateIndexFile(fq, indexReader, indexInfo, allowMismatch);
		}

		//writeBarcodeCount(saveDir+"/barcodeCount.txt", barcodeCount, false);
		writeBarcodeCount(saveDir+"/barcodeCount.txt", barcodeCount, true);		
	}

	
	private Map<String, FastqWriter> parseIndexInfo(String indexInfoFile, String saveDir, boolean allowMismatch) throws IOException {
		Map<String, FastqWriter> rtrn=new HashMap<String, FastqWriter>();
		
		Collection<String> lines=BEDFileParser.loadList(indexInfoFile, true);
		
        final FastqWriterFactory factory = new FastqWriterFactory();
		rtrn.put("remainder", factory.newWriter(new File(saveDir+"/remainder.fq")));
		
		for(String line: lines){
			String[] tokens=line.split("\t");
			String seq=tokens[0];
			FastqWriter writer=factory.newWriter(new File(saveDir+"/"+tokens[1]+".fq"));
			rtrn.put(seq, writer);
			
			// If we're allowing mismatches in the barcode, add degenerate sequences (e.g. for CAT add XAT, CXT, CAX)
			if (allowMismatch) {
				throw new UnsupportedOperationException("Need to copy back CountProbeKmers package from broad.pda.rap");
				/*
				Set<String> mismatchStrings = CountProbeKmers.getMismatchStrings(seq);
				for (String str : mismatchStrings) {
					rtrn.put(str, writer);
				}
				*/
			}
		}

		return rtrn;
	}


	private Map<String, Integer> readAndWriteByBarcode(FastqReader reader, Map<String, FastqWriter> indexInfo, boolean allowMismatch) throws IOException {
		Map<String, Integer> barcodeCounter=new HashMap<String, Integer>();
		
		while (reader.hasNext()) {
            final FastqRecord seq = reader.next();
            String barcode = extractBarcodeFromFastq(seq);
            
            writeFastqRecord(seq, barcode, indexInfo, barcodeCounter, allowMismatch);
	    }
		closeAll(indexInfo);
		return barcodeCounter;
	}
	
	
	public static String extractBarcodeFromFastq(FastqRecord record) {
		String id = record.getReadHeader().split("#")[1];
		return id.split("/")[0];
	}
	
	
	private void writeFastqRecord(final FastqRecord seq, final String barcode, final Map<String, FastqWriter> indexInfo, Map<String, Integer> barcodeCounter, boolean allowMismatch) {
    	FastqWriter writer;

    	if (indexInfo.containsKey(barcode)) {
    		writer = indexInfo.get(barcode);
    	} else if (allowMismatch) {
    		throw new UnsupportedOperationException("Need to copy back CountProbeKmers package from broad.pda.rap");
    		/*
    		Set<String> mismatchStrings = CountProbeKmers.getMismatchStrings(barcode);
    		Set<String> intersection = new HashSet<String>(mismatchStrings);
    		intersection.retainAll(indexInfo.keySet());
    		if (intersection.isEmpty()) {
    			writer = indexInfo.get("remainder");
    		} else {
    			writer = indexInfo.get(intersection.toArray()[0]);
    		}
    		*/
    	} else {
    		writer = indexInfo.get("remainder");
    	}
    	
        writer.write(seq);
        
        int count = 0;
        if (barcodeCounter.containsKey(barcode)) {
        	count = barcodeCounter.get(barcode);
        }
        count++;
        barcodeCounter.put(barcode, count);
	}
	
	
	private Map<String, Integer> readAndWriteByBarcodeSeparateIndexFile(FastqReader reader, FastqReader readerIndex, Map<String, FastqWriter> indexInfo, boolean allowMismatch) throws IOException {
		Map<String, Integer> barcodeCounter=new HashMap<String, Integer>();
		
		while (reader.hasNext() && readerIndex.hasNext()) {
            final FastqRecord seq = reader.next();
            final FastqRecord index = readerIndex.next();
            final String barcode = index.getReadString();
            
        	if(!seq.getReadHeader().split(" ")[0].equalsIgnoreCase(index.getReadHeader().split(" ")[0])){
        		throw new IllegalArgumentException("Indexes and Reads");
        	}
        	
            writeFastqRecord(seq, barcode, indexInfo, barcodeCounter, allowMismatch);
	    }
		 closeAll(indexInfo);
		 return barcodeCounter;
	}
	
	
	private void closeAll(Map<String, FastqWriter> indexInfo) throws IOException {
		for(String barcode: indexInfo.keySet()){
			FastqWriter writer=indexInfo.get(barcode);
			writer.close();
		}
	}
	

	/*
	private Map<String, Integer> countBarcodes(File indexFastq) throws IOException{
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(indexFastq)));
		String nextLine;
		Map<String, Integer> barcodeCounter=new TreeMap<String, Integer>();
	    while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
	    	if(nextLine.startsWith("@")){
	    		String firstLine=nextLine;
	        	String secondLine=reader.readLine();
	        	String thirdLine=reader.readLine();
	        	String fourthLine=reader.readLine();
	        	FastqSequence seq=new FastqSequence(firstLine, secondLine, thirdLine, fourthLine);
	        	String barcode=seq.getSequence();
	        	int count=0;
	        	if(barcodeCounter.containsKey(barcode)){count=barcodeCounter.get(barcode);}
	        	count++;
	        	barcodeCounter.put(barcode, count);
	    	}
	    }
	   return barcodeCounter;
	}
	*/
	
	public static void writeBarcodeCount(String save, Map<String, Integer> barcodeCount, boolean sortedByCount) throws IOException {
		BufferedWriter writer = new BufferedWriter(new FileWriter(save));
		
		if (!sortedByCount) {
			for(String barcode: barcodeCount.keySet()){
				writer.write(barcode+"\t"+barcodeCount.get(barcode)+"\n");
			}
		} else {
			List<Map.Entry<String, Integer>> list = new LinkedList(barcodeCount.entrySet());
			Collections.sort(list, new Comparator() {
				public int compare(Object o1, Object o2) {
					return ((Comparable) ((Map.Entry) (o1)).getValue()).compareTo(((Map.Entry) (o2)).getValue());
				}
			});
			Collections.reverse(list);

			for (Map.Entry<String, Integer> entry : list) {
				writer.write(entry.getKey() + "\t" + entry.getValue() + "\n");
			}
		}
		writer.close();
	}

	
	
	private static String USAGE="\n\nSplitFastqByIndex supports two file formats: one where indices are contained in a separate fastq file (old MiSeq format), and one where the indices are contained inline in the FASTQ file (i.e. HiSeq from Koch Institute).\n\n"+
								"\tjava -jar SplitFastqByIndex -in <fastq> -outdir <output directory> -indexInfo <tab-delimited file: barcodes in col1 and sample in col2, with a header row>\n"+
								"\t\t[-index <index fastq file>] [-allowMismatch]\n\n";
		
	public static void main(String[] args)throws IOException{
		ArgumentMap argmap = CLUtil.getParameters(args, USAGE, "full");
		String fastq = argmap.getInput();
		String out = argmap.getOutputDir();
		String indexInfo = argmap.getMandatory("indexInfo");
		String index = argmap.containsKey("index") ? argmap.getMandatory("index") : null;
		boolean allowMismatch = argmap.containsKey("allowMismatch");
		new SplitFastqByIndex(fastq, indexInfo, index, out, allowMismatch);
	}
}
