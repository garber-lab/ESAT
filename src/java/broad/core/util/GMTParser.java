package broad.core.util;

import java.io.*;
import java.util.*;

public class GMTParser {

	public static Map<String, Collection<String>> ParseGMTFile(File gmtFile, int minNum)throws IOException{
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(gmtFile)));
    	
    	Map<String, Collection<String>> data=new TreeMap<String, Collection<String>>();
        String nextLine;
        
        while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
          String[] tokens=nextLine.split("\t");
          String key=tokens[0];
          Collection<String> set=new TreeSet();
          for(int i=2; i<tokens.length; i++){set.add(tokens[i].toUpperCase());}
          data.put(key, set);
        }
        
        
        Map<String, Collection<String>> rtrn=new TreeMap();
        for(String name: data.keySet()){
        	Collection<String> set=data.get(name);
        	if(set.size()>minNum){rtrn.put(name, set);}
        }
        
        reader.close();
        return rtrn;
		
	}
	
	
	public static Map<String, Collection<String>> parseGMXFile(File gmxFile, int minNum)throws IOException{
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(gmxFile)));
    	
    	Map<String, Collection<String>> data=new TreeMap();
        String nextLine;
        
        int i=0;
        
        ArrayList<String> names=new ArrayList();
        while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
          String[] tokens=nextLine.split("\t");
          if(i==0){
        	  for(int j=0; j<tokens.length; j++){data.put(tokens[j], new TreeSet()); names.add(tokens[j]);}
          }
          if(i>1){
        	  for(int j=0; j<tokens.length; j++){
        		String name=names.get(j);
        		Set set=(Set)data.get(name);
        		if(tokens[j].length()>0){set.add(tokens[j]);}
        		data.put(name, set);
        		  
        	 }
          }
          i++;
        }
        
        Map<String, Collection<String>> rtrn=new TreeMap();
        for(String name: data.keySet()){
        	Collection<String> set=data.get(name);
        	if(set.size()>minNum){rtrn.put(name, set);}
        }
        
        reader.close();
        return rtrn;
	}
	
	public static Map<String, String> parseCHIPFile(File file)throws IOException{
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
    	Map list=new TreeMap();
		String nextLine;
        while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
          String[] tokens=nextLine.split("\t");
          String key=tokens[0];
          String val=tokens[1];
          list.put(key, val.toUpperCase());
        }
        reader.close();
        return list;	
	}
	
	public static void writeGMT(String save, Map<String, Collection<String>> geneSets)throws IOException{
		FileWriter writer=new FileWriter(save);
		
		for(String geneSet: geneSets.keySet()){
			Collection<String> genes=geneSets.get(geneSet);
			writer.write(geneSet+"\t"+geneSet);
			for(String gene: genes){writer.write("\t"+gene);}
			writer.write("\n");
		}
		
		writer.close();
	}


	public static Map<String, Collection<String>> ParseGeneSetFile(File gmtFile, int minNum) throws IOException {
		if(gmtFile.getName().endsWith("gmt") || gmtFile.getName().endsWith("GMT")){return ParseGMTFile(gmtFile, minNum);}
		else{return parseGMXFile(gmtFile, minNum);}
	}


	public static Map<String, Collection<String>> parseCHIPFileByName(File chipFile) throws IOException {
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(chipFile)));
    	Map<String, Collection<String>> rtrn=new TreeMap<String, Collection<String>>();
		String nextLine;
        while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
          String[] tokens=nextLine.split("\t");
          String key=tokens[0];
          String val=tokens[1].toUpperCase();
          Collection<String> set=new TreeSet<String>();
          if(rtrn.containsKey(val)){set=rtrn.get(val);}
          set.add(key);
          rtrn.put(val, set);
        }
        reader.close();
        return rtrn;	
	}
	
}
