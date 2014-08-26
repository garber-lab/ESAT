package broad.pda.geneexpression.agilent;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.core.datastructures.MatrixWithHeaders;
import broad.pda.geneexpression.ExpressionExperimentInfo;


public class AgilentUtils {

	public static Map<String, AgilentProbe> parseAgilentExperiment(File file) throws IOException{
		Map<String, AgilentProbe> rtrn=new TreeMap<String, AgilentProbe>();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
    	String nextLine;
    	boolean start=false;
    	String header="";
        while ((nextLine = reader.readLine()) != null ) {
        	if(nextLine.startsWith("FEATURES")){header=nextLine; start=true;}
        	if(nextLine.startsWith("DATA") && start){
        		AgilentProbe probe=new AgilentProbe(nextLine, header);
        		rtrn.put(probe.getProbename(), probe);
        	}
        }
		
		
		return rtrn;
	}

	public static AgilentArrayStats parseAgilentExperimentGlobalStats(File file) throws IOException {
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
    	String nextLine;
    	boolean start=false;
    	String header="";
        while ((nextLine = reader.readLine()) != null ) {
        	if(nextLine.startsWith("STATS")){header=nextLine; start=true;}
        	if(nextLine.startsWith("DATA") && start){
        		AgilentArrayStats array=new AgilentArrayStats(nextLine, header);
        		return array;
        	}
        }
		
		
		return null;
	}

	public static Map<String, ExpressionExperimentInfo> parseExperimentInfoFile(File experimentInfoFile) throws IOException {
		Map<String, ExpressionExperimentInfo> rtrn=new TreeMap<String, ExpressionExperimentInfo>();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(experimentInfoFile)));
    	String nextLine;
    	while ((nextLine = reader.readLine()) != null ) {
    		if(!nextLine.startsWith("Barcode")){
    			ExpressionExperimentInfo info=new ExpressionExperimentInfo(nextLine);
    			rtrn.put(info.getExperimentBarcode(), info);
    		}
    	}
    	
		return rtrn;
	}
	
	public static Map<String, String> parseExperimentInfoFileToName(File experimentInfoFile) throws IOException {
		return parseExperimentInfoFileToName(experimentInfoFile, false);
	}
	
	public static Map<String, String> parseExperimentInfoFileToName(File experimentInfoFile, boolean passingFilter) throws IOException {
		Map<String, String> rtrn=new TreeMap<String, String>();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(experimentInfoFile)));
    	String nextLine;
    	while ((nextLine = reader.readLine()) != null ) {
    		if(!nextLine.startsWith("Barcode")){
    			ExpressionExperimentInfo info=new ExpressionExperimentInfo(nextLine);
    			if(!passingFilter || info.passedQC()){rtrn.put(info.getExperimentBarcode(), info.getExperimentName());}
    		}
    	}
    	
		return rtrn;
	}
	
	public static Map<String, Collection<ExpressionExperimentInfo>> parseExperimentInfoFileToGroup(File experimentInfoFile) throws IOException {
		Map<String, Collection<ExpressionExperimentInfo>> rtrn=new TreeMap<String, Collection<ExpressionExperimentInfo>>();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(experimentInfoFile)));
    	String nextLine;
    	while ((nextLine = reader.readLine()) != null ) {
    		if(!nextLine.startsWith("Barcode")){
    			ExpressionExperimentInfo info=new ExpressionExperimentInfo(nextLine);
    			if(info.passedQC()){
    				String group="";
    				if(!info.isControl()){group=info.getExperimentName();}
    				else{group="Control";}
    				Collection<ExpressionExperimentInfo> set=new ArrayList<ExpressionExperimentInfo>();
    				if(rtrn.containsKey(group)){set=rtrn.get(group);}
    				set.add(info);
    				rtrn.put(group, set);
    			}
    		}
    	}
    	
		return rtrn;
	}
	
	public static Map<String, Collection<String>> parseExperimentInfoFileToGroups(File experimentInfoFile) throws IOException {
		Map<String, Collection<String>> rtrn=new TreeMap<String, Collection<String>>();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(experimentInfoFile)));
    	String nextLine;
    	while ((nextLine = reader.readLine()) != null ) {
    		if(!nextLine.startsWith("Barcode")){
    			ExpressionExperimentInfo info=new ExpressionExperimentInfo(nextLine);
    			if(info.passedQC()){
    				String group="";
    				if(!info.isControl()){group=info.getExperimentName();}
    				else{group="Control";}
    				Collection<String> set=new ArrayList<String>();
    				if(rtrn.containsKey(group)){set=rtrn.get(group);}
    				set.add(info.getExperimentBarcode());
    				rtrn.put(group, set);
    			}
    		}
    	}
    	
		return rtrn;
	}
	
	public static Map<String, Double> makeDataMatrix(File agilentFile, Collection<String> nonFlaggedGenes) throws IOException{
		Map<String, Double> rtrn=new TreeMap<String, Double>();
		
		
		String name=agilentFile.getName().split("\\.")[0];
		AgilentArrayStats stats =parseAgilentExperimentGlobalStats(agilentFile);
		Map<String, AgilentProbe> data=AgilentUtils.parseAgilentExperiment(agilentFile);
		
		
		//columns
		Map<String, AgilentProbe> map=data;
		for(String pid: map.keySet()){
			AgilentProbe probe=map.get(pid);
			if(probe.isProbeWellAboveBackground()){nonFlaggedGenes.add(pid);}
			rtrn.put(pid, probe.getNormalizedSignal(stats.spikeInDetectionLimit, stats.saturationValue));//MG I'm removing the normalization to 75%
		}
				
		return rtrn;
	}
	
	public static MatrixWithHeaders makeDataMatrix(File[] agilentFiles) throws IOException{
		Collection<String> nonFlaggedGenes=new TreeSet<String>();
		
		System.err.println(agilentFiles.length);
		
		Map<String, Double>[] maps=new Map[agilentFiles.length];
		List<String> names=new ArrayList<String>();
		
		for(int i=0; i<agilentFiles.length; i++){
			String name=agilentFiles[i].getName().split("\\.")[0];
			System.err.println("Sample "+i+" "+agilentFiles[i]);
			maps[i]=makeDataMatrix(agilentFiles[i], nonFlaggedGenes);
			names.add(name);
		}
		
		MatrixWithHeaders rtrn=new MatrixWithHeaders(new ArrayList(nonFlaggedGenes), names);
		for(int i=0; i<agilentFiles.length; i++){
			String name=agilentFiles[i].getName().split("\\.")[0];
			for(String gene: nonFlaggedGenes){
				rtrn.set(gene, name, maps[i].get(gene));
			}
		}
				
		return rtrn;
	}
		
	
	

	public static MatrixWithHeaders makeDataMatrix(File[] agilentFiles, Map<String, ExpressionExperimentInfo> experimentInfo) throws IOException{
		ArrayList<File> files=new ArrayList<File>();
		
		for(int i=0; i<agilentFiles.length; i++){
			String name=agilentFiles[i].getName().split("\\.")[0];
			ExpressionExperimentInfo info=experimentInfo.get(name);
			System.err.println("Sample #"+i+" "+name+" "+info);
			if(info!=null && info.passedQC()){files.add(agilentFiles[i]);}
		}
		
		File[] filteredFiles=files.toArray(new File[1]);
				
		return makeDataMatrix(filteredFiles);
		
	}

	public static Map<String, Collection<String>> parseExperimentInfoGroups(String experimentInfoFile) throws IOException {
		Map<String, Collection<String>> rtrn=new TreeMap<String, Collection<String>>();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(experimentInfoFile)));
    	String nextLine;
    	while ((nextLine = reader.readLine()) != null ) {
    		if(!nextLine.startsWith("Barcode")){
    			ExpressionExperimentInfo info=new ExpressionExperimentInfo(nextLine);
    			if(info.passedQC()){
    				String group=info.getExperimentName();
    				Collection<String> set=new TreeSet<String>();
    				if(rtrn.containsKey(group)){
    					set=rtrn.get(group);
    				}
    				set.add(info.getExperimentBarcode());
    				rtrn.put(group, set);
    			}
    		}
    	}
    	
		return rtrn;
	}

	public static Map<String, Collection<String>> parseExperimentInfoNonControls(String experimentInfoFile) throws IOException {
		Map<String, Collection<String>> rtrn=new TreeMap<String, Collection<String>>();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(experimentInfoFile)));
    	String nextLine;
    	while ((nextLine = reader.readLine()) != null ) {
    		if(!nextLine.startsWith("Barcode")){
    			ExpressionExperimentInfo info=new ExpressionExperimentInfo(nextLine);
    			if(info.passedQC() && !info.isControl()){
    				String group=info.getExperimentName();
    				Collection<String> set=new TreeSet<String>();
    				if(rtrn.containsKey(group)){
    					set=rtrn.get(group);
    				}
    				set.add(info.getExperimentBarcode());
    				rtrn.put(group, set);
    			}
    		}
    	}
    	
		return rtrn;
	}
	
	public static Collection<String> parseControlSamples(String experimentInfoFile) throws IOException {
		Collection<String> rtrn=new TreeSet<String>();
		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(experimentInfoFile)));
    	String nextLine;
    	while ((nextLine = reader.readLine()) != null ) {
    		if(!nextLine.startsWith("Barcode")){
    			ExpressionExperimentInfo info=new ExpressionExperimentInfo(nextLine);
    			if(info.passedQC() && info.isControl()){
    				rtrn.add(info.getExperimentBarcode());
    			}
    		}
    	}
    	
		return rtrn;
	}
	
}
