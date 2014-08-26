package broad.pda.differentialExpression;

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

import broad.core.datastructures.MatrixWithHeaders;

public class MergeGSEAOutputIntoMatrix {

	public MergeGSEAOutputIntoMatrix(File[] files, String save) throws IOException{
		
		ArrayList<String> names=new ArrayList();
		Map<String, Double>[] maps=new Map[files.length];
		
		for(int i=0; i<files.length; i++){
			String name=files[i].getName();
			maps[i]=parse(files[i]);
			names.add(name);
		}
		
		write(save, maps, names);
	}

	private Map<String, Double> parse(File file) throws IOException {
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
    	
    	Map<String, Double> rtrn=new TreeMap<String, Double>();
        String nextLine;
        
        int i=0;
        while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
        	if(i>2){
        		String[] tokens=nextLine.split("\t");
        		rtrn.put(tokens[0], new Double(tokens[2]));
        	}
        	i++;
        }
        
        reader.close();
        return rtrn;
	}

	private void write(String save, Map<String, Double>[] maps, ArrayList<String> names) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		writer.write("GeneSet Name");
		for(String name: names){writer.write("\t"+name);}
		writer.write("\n");
		
		for(String gene: maps[0].keySet()){
			writer.write(gene);
			for(int i=0; i<maps.length; i++){
				double val=maps[i].get(gene);
				writer.write("\t"+val);
			}
			writer.write("\n");
		}
		
		
		writer.close();
	}
	
	public static void main(String [] args)throws IOException{
		if(args.length>1){
			File[] files=new File(args[0]).listFiles();
			String save=args[1];
			new MergeGSEAOutputIntoMatrix(files, save);
		}
		else{System.err.println(usage);}
	}

	static String usage=" args[0]=files \n args[1]=save";
	
}
