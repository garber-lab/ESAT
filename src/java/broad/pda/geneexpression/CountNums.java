package broad.pda.geneexpression;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Collection;
import java.util.TreeSet;

import broad.core.datastructures.MatrixWithHeaders;
import broad.core.error.ParseException;

public class CountNums {

	public CountNums(MatrixWithHeaders fdr, Collection<String> directTargets, String save, double cutoff)throws IOException{
		FileWriter writer=new FileWriter(save);
		writer.write("Experiment\tNumber of Direct Hits\tNumber of total hits\tTotal number of direct targets\tPercent direct identified\t% of identified that are direct\n");
		
		for(String experiment: fdr.getColumnNames()){
			int numDirectHits=0;
			int numTotalHits=0;
			for(String gene: fdr.getRowNames()){
				double val=fdr.get(gene, experiment);
				if(val<cutoff){
					numTotalHits++;
					if(directTargets.contains(gene)){numDirectHits++;}
				}
			}
			double percentDirectId=(double)numDirectHits/directTargets.size();
			double percentTotalDirect=(double)numDirectHits/numTotalHits;
			writer.write(experiment+"\t"+numDirectHits+"\t"+numTotalHits+"\t"+directTargets.size()+"\t"+percentDirectId+"\t"+percentTotalDirect+"\n");
		}
		writer.close();
	}
	
	private static Collection<String> parseNames(String file) throws IOException {
		Collection<String> rtrn=new TreeSet<String>();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String nextLine; 
		while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {rtrn.add(nextLine.split("\t")[0].toUpperCase().trim());}
		reader.close();
		return rtrn;
	}
	
	public static void main(String[] args)throws IOException, ParseException{
		if(args.length>3){
			MatrixWithHeaders fdr=new MatrixWithHeaders(args[0]);
			Collection<String> directTargets=parseNames(args[1]);
			String save=args[2];
			double cutoff=new Double(args[3]);
			new CountNums(fdr, directTargets, save, cutoff);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=fdr matrix \n args[1]=direct targets \n args[2]=save \n args[3]=cutoff";
}
