package broad.core.motif.meme;

import jaligner.Sequence;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.text.NumberFormat;
import java.util.List;
import java.util.Stack;

public class MemeResultReader {
	Stack<MemeMotif> motifs;
	public MemeResultReader(BufferedReader br) throws IOException {
		boolean inMatches = false;
		boolean inMatrix = false;
		int distToData = 0;
		motifs = new Stack<MemeMotif>();
		String line = null;
		int motifNum = 0;
		
		while( (line = br.readLine())!= null) {
			line = line.trim();
			
			if(line.isEmpty() || (line.startsWith("---") && (distToData<=0))) {
				inMatches = false;
				distToData  = 3;
				inMatrix = false;
			}else if(line.startsWith("MOTIF")) { //At general motif information
				String [] info = line.split("\t");
				MemeMotif newMotif = new MemeMotif(info[0]);
				motifs.push(newMotif);
				String [] data = info[1].split("\\s+");
				newMotif.setWidth(Integer.parseInt(data[2]));
				newMotif.setSites(Integer.parseInt(data[5]));
				newMotif.setLlr(Float.parseFloat(data[8]));
				newMotif.setEval(Double.parseDouble(data[11]));
				
				motifNum++;
			}else if(line.contains("Motif "+motifNum+" sites sorted by position p-value")) {
				inMatches = true;
			} else if (inMatches && distToData-- <=0) {
				String [] info = line.split("\\s+");
				//System.err.println("info[0] "+ info[0]+" Line: " +line );
				Sequence match = new Sequence(info[5]);
				match.setId(info[0]);
				motifs.peek().addMatch(match, Integer.parseInt(info[2]), Double.parseDouble(info[3]));
				//System.err.println("Sequence added: " + match.getSequence() + " named " + match.getId());
				distToData--;
			} else if(line.contains("Motif "+motifNum+" position-specific probability matrix")) { 
				inMatrix = true;
			} else if (inMatrix && distToData-- <=0) {
				String [] info = line.split("\\s+");
				double [] column = new double[4];
				for(int i = 0; i < 4; i++ ) {
					column[i] =Double.parseDouble(info[i]);
				}
				motifs.peek().addPWMColumn(column);
			} else {
				//System.err.println(inMatches + ","+distToData+","+inMatrix+" Ignoring line " + line);
			}
		}
		
	}

	
	public MemeResultReader(InputStream is) throws IOException {
		this(new BufferedReader(new InputStreamReader(is)));
	}
	
	public MemeResultReader(String inputFile) throws IOException {
		this(new BufferedReader(new FileReader(new File(inputFile))));
	}
	
	public MemeResultReader(File inputFile) throws IOException {
		this(new BufferedReader(new FileReader(inputFile)));
	}
	
	public List<MemeMotif> getMotifs() {return motifs;}

}
