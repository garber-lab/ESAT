package umms.core.fastq.utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import picard.cmdline.CommandLineProgram;
import picard.cmdline.Option;
import umms.core.sequence.SequenceUtils;

public class BCEvaluator extends CommandLineProgram {
	int minDistance = Integer.MAX_VALUE;
	int maxDistance = 0;
	int [] histogram = null;
	
    @Option (doc="2-colum Tab delimited file with barcode to sample mapping" , shortName="M", optional=false)
    public File BC_SAMPLE_MAP;
    
    @Option (doc="Output File or standard output if non is specified", shortName="O", optional=true)
    public File OUT = null;
	
	static HashMap<String, String> readBCMap(File file) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(file));
		HashMap<String, String> map = new HashMap<String, String>();
		String line = null;
		while (  (line = br.readLine()) != null) {
			String [] info = line.split("\\s+");
			map.put(info[0], info[1]);
		}
		br.close();
		return map;
	}
	
	private void evaluateErrorDetectionAndCoorrection(HashMap<String, String> bcSampleMap) {
		List<String> barcodes = new ArrayList<String>(bcSampleMap.keySet());
		
		histogram = new int [barcodes.get(0).length()+1];
		for (int i = 0; i < barcodes.size() - 1; i++) {
			for (int j = i + 1; j < barcodes.size(); j++) {
				int dist = SequenceUtils.hamming(barcodes.get(i), barcodes.get(j));
				histogram[dist]++;
				minDistance = dist < minDistance ? dist : minDistance;
				maxDistance = dist > maxDistance ? dist : maxDistance;
			}
 		}
		
	}
	
    private void writeOut(File o) throws IOException {
    	BufferedWriter bw = o == null ? 
    				new BufferedWriter(new PrintWriter(System.out)) :
    				new BufferedWriter(new FileWriter(OUT));
    	
    				
    	bw.write("minimum distance: " + minDistance);
    	bw.newLine();
    	bw.write("maximum distance: " + maxDistance);
    	bw.newLine();
    	bw.newLine();
    	bw.write("Histogram:");
    	bw.newLine();
    	for (int i = 0; i <= maxDistance; i++) {
    		bw.write(i + spaces(i, maxDistance) + stars(histogram[i]));
    		bw.newLine();
    	}
    	bw.close();
	}
	
	private String stars(int n) {
		StringBuffer stars = new StringBuffer();
		for (int i = 0; i < n; i++) {
			stars.append("*");
		}
		return stars.toString();
	}

	private String spaces(int number, int maxnumber) {
		int spaces = String.valueOf(maxnumber).length() - String.valueOf(number).length();
		StringBuffer spaceString = new StringBuffer("  ");
		for (int i = 0; i < spaces; i++) {
			spaceString.append(" ");
		}
		return spaceString.toString();
	}

	protected int doWork() {
		
		try {
			HashMap<String, String> bcSampleMap = readBCMap(BC_SAMPLE_MAP);
			evaluateErrorDetectionAndCoorrection(bcSampleMap);
			writeOut(OUT);
		} catch (FileNotFoundException e) {
			System.out.println("Could not access the barcode sample mapping file or write output file");
			e.printStackTrace();
			return 1;
		} catch (IOException e) {
			System.out.println("Error accessing barcode sample mapping file");
			e.printStackTrace();
			return 1;
		}
		
		return 0;
	}



	/** Stock main method. */
    public static void main(final String[] argv) {
        System.exit(new BCEvaluator().instanceMain(argv));
    }

    
}
