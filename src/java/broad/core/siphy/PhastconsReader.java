package broad.core.siphy;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import broad.core.annotation.BasicGenomicAnnotation;

public class PhastconsReader {
	ArrayList<PhastConsElement> fastconList;
	File source;

	public PhastconsReader(String fileName) throws FileNotFoundException {
		super();
		fastconList = new ArrayList<PhastConsElement>();
		source = new File(fileName);
		BufferedReader br = new BufferedReader(new FileReader(source));
		String line;
		try {
			while((line = br.readLine()) != null) {
				if(line.startsWith("#")){
					continue;
				}
				PhastConsElement summary = new PhastConsElement(line.split("\t"));
				fastconList.add(summary);
			}
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			try {
				System.out.print("Closing "+fileName);
				br.close();
				System.out.print(" ..... Closed\n");
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}		
	}
	
	public List<PhastConsElement> getElements() {
		return fastconList;
	}
	
	public static class PhastConsElement extends BasicGenomicAnnotation {
		private int score;
		private String chrom;

		public PhastConsElement(String [] rawData) {
			chrom = rawData[0].substring(3);
			setStart(Integer.parseInt(rawData[1]));
			setEnd(Integer.parseInt(rawData[2]));
			setName(rawData[3]);
			score = Integer.parseInt(rawData[4]);
		}
		
		public String toString() {
			StringBuffer b = new StringBuffer("0\t");
			b.append(chrom).append("\t")
				.append(getStart()).append("\t")
				.append(getEnd()).append("\t")
				.append(getName()).append("\t")
				.append(score);
			
			return b.toString();
		}
	}

}
