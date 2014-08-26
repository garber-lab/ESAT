package broad.core.motif.meme;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class AMEResultReader {
	
	private List<AMEMatch> results;
	private static final Pattern dataPattern = Pattern.compile("^[0-9]+\\.");
	public AMEResultReader(BufferedReader br) throws IOException {
		results = new ArrayList<AMEMatch>(); 
		String line = null;
	//	System.err.println("line " + line);
		while( (line = br.readLine())!= null) {
			line = line.trim();
			Matcher m = dataPattern.matcher(line);
			if(m.find()) {
				//System.err.println("Data line: " + line);
				String [] info = line.split("\\s+");
				String motif = info[5];
				double pval = Double.parseDouble(info[info.length - 4]);
				double correctedPval = Double.parseDouble(info[info.length - 1].replace(")",""));
				results.add(new AMEMatch(motif, pval, correctedPval));
			}
		}
		
	}
	
	public AMEResultReader(InputStream is) throws IOException {
		this(new BufferedReader(new InputStreamReader(is)));
	}
	
	public AMEResultReader(String inputFile) throws IOException {
		this(new BufferedReader(new FileReader(new File(inputFile))));
	}
	
	public AMEResultReader(File inputFile) throws IOException {
		this(new BufferedReader(new FileReader(inputFile)));
	}
	
	public static class AMEMatch {
		String motif;
		double pval;
		double correctedPval;
		

		public AMEMatch(String motif, double pval, double correctedPval) {
			this.motif = motif;
			this.pval = pval;
			this.correctedPval = correctedPval;
		}
		
		public String getMotif() {
			return motif;
		}

		public double getPval() {
			return pval;
		}

		public double getCorrectedPval() {
			return correctedPval;
		}
	}
}
