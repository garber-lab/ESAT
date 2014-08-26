package broad.pda.geneexpression.agilent;


import java.io.*;
import java.util.*;

import broad.core.datastructures.MatrixWithHeaders;
import broad.core.math.Statistics;
import broad.core.util.ParseGCTFile;

public class AgilentToGCT {
	static String logHeaderName="gProcessedSignal";
	static String uidHeaderName="ProbeName";
	static String geneHeaderName="SystematicName";
	//static String accessionName="SystematicName";
	static String quantileMetric="gPercentileIntensityProcessedSignal";
	
	double floor=0;
	Map<String, String> chipInfo;
	double[] quantileMetrics;
	

	public AgilentToGCT(File[] files, String save, File experimentInfo)throws IOException{
		//Map<String, Double>[] maps=new Map[files.length];
		MatrixWithHeaders data=null;
		
		
		Map<String, String> infoMap=parse(experimentInfo);
		this.chipInfo=new TreeMap();
		
		this.quantileMetrics=new double[files.length];
		for(int i=0; i<files.length; i++){
			Map<String, Double> map=parseAgilent(files[i], i);
			addToMatrix(data, map, files, files[i]);
			System.err.println(files[i]);
			System.gc();
		}
		
		
		//maps=medianCenter(maps);
		//Map<String, String> chipMap=ParseGCTFile.parseChipFile(chipFile);
		
		data.setPIDToName(chipInfo);
		data.writeGCT(save);
		
		//write(save+".norm.gct", maps, chipInfo, files, infoMap, quantileMetrics);
		//write(save+".Raw.gct", maps, chipInfo, files, infoMap, null);
	}
	
	private void addToMatrix(MatrixWithHeaders data, Map<String, Double> map, File[] files, File file) {
		if(data==null){
			data=new MatrixWithHeaders(new ArrayList(map.keySet()), getNames(files));
		}
		for(String geneName: map.keySet()){
			data.set(geneName, file.getName(), map.get(geneName));
		}
	}

	private List<String> getNames(File[] files) {
		ArrayList<String> rtrn=new ArrayList<String>();
		
		for(int i=0; i<files.length; i++){
			rtrn.add(files[i].getName());
		}
		
		return rtrn;
	}

	private Map<String, Double>[] medianCenter(Map<String, Double>[] maps) {
		Map<String, Double>[] rtrn=new Map[maps.length];
		for(int i=0; i<maps.length; i++){
			rtrn[i]=new TreeMap();
			double median=Statistics.medianCollection(maps[i].values());
			double avg=Statistics.average(maps[i].values());
			System.err.println(i+" "+median+" "+avg);
			for(String gene: maps[i].keySet()){
				double val=maps[i].get(gene);
				//rtrn[i].put(gene, val/avg);
				rtrn[i].put(gene, val);
			}
		}
		return rtrn;
	}

	public Map<String, String> parse(File file)throws IOException{
		if(file==null){return new TreeMap();}
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
    	
    	Map data=new TreeMap();
        String nextLine;
        int count=0;
        while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
           String[] tokens=nextLine.split("\t");
           String pid=tokens[0].trim();
           String name=tokens[1].trim();
           data.put(pid, name);
        }
        
        
        reader.close();
        return data;
	}
	
	
	private Map<String, Double> parseAgilent(File file, int index)throws IOException{
		Map<String, Double> rtrn=new TreeMap();
		
		int nameIndex=6;
		int valIndex=12;
		int descIndex=13;
		int accessionIndex=-1;
		boolean start=false;
		
		//TODO Grab the Norm factor and divide
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
    	String nextLine;
        while ((nextLine = reader.readLine()) != null ) {
        	
        	if(nextLine.startsWith("DATA") && start){
        		String[] tokens=nextLine.split("\t");
				//if(rtrn.containsKey(tokens[6])){System.err.println("Duplicate: "+tokens[6]);}
				String name=tokens[nameIndex];
				double val=new Double(tokens[valIndex]);
				this.chipInfo.put(name, tokens[descIndex]);
        		//String accession=tokens[accessionIndex];
				
        		//if(!accession.isEmpty()){
				if(val<this.floor){val=floor;}
				rtrn.put(name.toUpperCase(), val);
        		//}
        	}
        	else if(nextLine.split("\t")[0].equalsIgnoreCase("Features")){
        		start=true;
        		System.err.println("Feature: "+nextLine);
        		String[] tokens=nextLine.split("\t");
        		for(int i=0; i<tokens.length; i++){
        			if(tokens[i].equalsIgnoreCase(logHeaderName)){valIndex=i;}
        			else if(tokens[i].equalsIgnoreCase(uidHeaderName)){nameIndex=i;}
        			else if(tokens[i].equalsIgnoreCase(geneHeaderName)){descIndex=i;}
        			//else if(tokens[i].equalsIgnoreCase(accessionName)){accessionIndex=i;}
        		}
        	//	System.err.println(valIndex+" "+nameIndex);
        		
        	}
        	
        	else if(nextLine.split("\t")[0].equalsIgnoreCase("Stats")){
        		String[] tokens=nextLine.split("\t");
        		int normIndex=0;
        		for(int i=0; i<tokens.length; i++){
        			if(tokens[i].equalsIgnoreCase(quantileMetric)){normIndex=i;}
        		}
        		nextLine=reader.readLine();
        		tokens=nextLine.split("\t");
        		this.quantileMetrics[index]=new Double(tokens[normIndex]);
        		System.err.println(index+" "+this.quantileMetrics[index]);
        	}
        	
        	
        //	else{System.err.println(nextLine.split("\t")[0]);}
        }
        
        
        reader.close();
        return rtrn;
	}
	
	private Map parseFile(File file)throws IOException{
		Map rtrn=new TreeMap();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
    	String nextLine;
        while ((nextLine = reader.readLine()) != null ) {
        	
        	String[] tokens=nextLine.split("\t");
			rtrn.put(tokens[0], new Double(tokens[1]));
        	
        	
        }
        
        
        reader.close();
        return rtrn;
	}
	
	private void write(String save, Map<String, Double>[] maps, Map<String, String> info, File[] files, Map<String, String> experimentMap, double[] quantiles)throws IOException{
		FileWriter writer=new FileWriter(save);
		writer.write("#1.2\n");
		writer.write(count(maps[0])+"\t"+maps.length+"\n");
		writer.write("PID\tName");
		for(int i=0; i<files.length; i++){
			String name=experimentMap.get(files[i].getName().split("\\.")[0]);
			if(name==null || name.isEmpty()){name=(files[i].getName().split("\\.")[0]);}
			writer.write("\t"+name);	
		}
		writer.write("\n");
		
		for(String pid: maps[0].keySet()){
			if(info.containsKey(pid)){
			String str="";
			boolean write=true;
			str+=(pid+"\t"+info.get(pid));
			for(int i=0; i<maps.length; i++){
				if(maps[i].get(pid)!=null){
					double val=maps[i].get(pid);
					if(quantiles!=null){val=val/quantiles[i];}
					str+=("\t"+(Math.log(val)/Math.log(2)));
					//str+=("\t"+val);
				}
				else{write=false;}
			}
			if(write){writer.write(str+"\n");}
			}
		}
		writer.close();
	}
	
	
	private int count(Map<String, Double> map) {
		int counter=0;
		
		for(String pid: map.keySet()){
			Double val=map.get(pid);
			if(val!=null){counter++;}
		}
		return counter;
	}

	private void write(String save, Map<String, Double>[] maps, Map<String, String> info, File[] files)throws IOException{
		FileWriter writer=new FileWriter(save);
		writer.write("#1.2\n");
		writer.write(maps[0].size()+"\t"+maps.length+"\n");
		writer.write("PID\tName");
		for(int i=0; i<files.length; i++){writer.write("\t"+files[i].getName());}
		writer.write("\n");
		
		
		for(String pid: maps[0].keySet()){
			if(info.containsKey(pid)){
			String str="";
			boolean write=true;
			str+=(pid+"\t"+info.get(pid));
			for(int i=0; i<maps.length; i++){
				if(maps[i].get(pid)!=null){
				str+=("\t"+Math.log(maps[i].get(pid)));
				}
				else{write=false;}
			}
			if(write){writer.write(str+"\n");}
			}
		}
		writer.close();
	}
	
	
	public static void main(String[] args)throws IOException{
		if(args.length>1){
		File[] files=new File(args[0]).listFiles();
		String save=args[1];
		if(args.length>2){File experimentInfo=new File(args[2]);
		new AgilentToGCT(files, save, experimentInfo);
		}
		else{
			new AgilentToGCT(files, save, null);
		}
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=files \n args[1]=save \n args[2]=Chip File";
}