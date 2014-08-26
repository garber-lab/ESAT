package broad.pda.geneexpression.agilent;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.core.datastructures.MatrixWithHeaders;
import broad.core.error.ParseException;
import broad.core.math.Statistics;
import broad.pda.geneexpression.ExpressionExperimentInfo;

public class CollapseProbesByMinPValue {
	
	public CollapseProbesByMinPValue(File file, String save) throws IOException{
		Map<String, Collection<String>> expressionInfo=parseInfo(file);
		Map<String, Double> expressionPValue=parsePValue(file);
		
		writeMin(save, expressionPValue, expressionInfo);
		
	}
	
	private Map<String, Double> parsePValue(File file) throws NumberFormatException, IOException {
		Map<String, Double> rtrn=new TreeMap<String, Double>();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
    	String nextLine;
    	int counter=0;
    	while ((nextLine = reader.readLine()) != null ) {
    		if(counter>0){
    		String[] tokens=nextLine.split("\t");
    		rtrn.put(nextLine, new Double(tokens[14]));
    		}
    		counter++;
    	}
    	
		return rtrn;
	}

	private Map<String, Collection<String>> parseInfo(File file) throws IOException {
		Map<String, Collection<String>> rtrn=new TreeMap<String, Collection<String>>();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
    	String nextLine;
    	int counter=0;
    	while ((nextLine = reader.readLine()) != null ) {
    		if(counter>0){
    		String[] tokens=nextLine.split("\t");
    		Collection<String> set=new TreeSet();
    		if(rtrn.containsKey(tokens[0])){set=rtrn.get(tokens[0]);}
    		set.add(nextLine);
    		rtrn.put(tokens[0], set);
    		}
    		counter++;
    	}
    	
		return rtrn;
	}

	private void writeMin(String save, Map<String, Double> expressionPValue,Map<String, Collection<String>> expressionInfo) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String linc: expressionInfo.keySet()){
			Collection<String> lines=expressionInfo.get(linc);
			double minP=1;
			String minLinc=null;
			for(String line: lines){
				double p=expressionPValue.get(line);
				if(p<=minP){minP=p; minLinc=line;}
			}
			writer.write(minLinc+"\n");
		}
		
		writer.close();
	}

	public static void main(String[] args)throws IOException{
		File file=new File(args[0]);
		String save=args[1];
		new CollapseProbesByMinPValue(file, save);
	}
	
}