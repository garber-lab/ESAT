package broad.core.util;

import java.io.*;
import java.util.*;

import broad.core.datastructures.MatrixWithHeaders;

import Jama.Matrix;



public class ParseGCTFile {

	public static Map loadData(File file)throws IOException{
			BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
	    	
	    	Map data=new TreeMap();
	        String nextLine;
	        int count=0;
	        while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
	           
	        	if(count>2){
				String[] tokens=nextLine.split("\t");
				String key=tokens[0];
				ArrayList list=new ArrayList();
				for(int i=2; i<tokens.length; i++){
					list.add(new Double(tokens[i]));
				}
				data.put(key, list);
	        	
	        	}
	        	else if(count==2){
	        		String[] tokens=nextLine.split("\t");
	    			String key="header";
	    			ArrayList list=new ArrayList();
	    			for(int i=2; i<tokens.length; i++){
	    				list.add(tokens[i]);
	    			}
	    			data.put(key, list);
	        	}
	        	count++;
	        }
	        
	        
	        reader.close();
	        return data;
			
		}
		
		
		public static Map<String, String> loadPIDName(File file)throws IOException{
			BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
	    	
	    	Map data=new TreeMap();
	        String nextLine;
	        int count=0;
	        while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
	           
	        	if(count>2){
				String[] tokens=nextLine.split("\t");
				
				data.put(tokens[0], tokens[1]);
	        	
	        	}
	        	
	        	count++;
	        }
	        
	        
	        reader.close();
	        return data;
			
		}
		
		

		
		
		private static double[] convert(ArrayList<Double> list){
			double[] rtrn=new double[list.size()];
			
			int i=0;
			for(Double val: list){
				rtrn[i++]=val;
			}
			
			return rtrn;
		}
		
		private static ArrayList convert(double[] vals){
			ArrayList list=new ArrayList();
			
			
			for(int i=0; i<vals.length; i++){
				list.add(vals[i]);
			}
			
			return list;
		}
		
		public static Map<String, ArrayList> parseCDT(File file)throws IOException{
			BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
	    	
	    	Map data=new TreeMap();
	        String nextLine;
	        int count=0;
	        while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
	           
	        	if(count>2){
				String[] tokens=nextLine.split("\t");
				String key=tokens[1];
				ArrayList list=new ArrayList();
				for(int i=4; i<tokens.length; i++){
					list.add(new Double(tokens[i]));
				}
				data.put(key, list);
	        	
	        	}
	        	else if(count==0){
	        		String[] tokens=nextLine.split("\t");
	    			String key="header";
	    			ArrayList list=new ArrayList();
	    			for(int i=4; i<tokens.length; i++){
	    				list.add(tokens[i]);
	    			}
	    			data.put(key, list);
	        	}
	        	count++;
	        }
	        
	        
	        reader.close();
	        return data;
		}
		
		public static Map<String, String> parseChipFile(File file)throws IOException{
			BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
	    	
	    	Map data=new TreeMap();
	        String nextLine;
	        int count=0;
	        while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
	           String[] tokens=nextLine.split("\t");
	           String pid=tokens[0].toUpperCase().trim();
	           String name=tokens[1].toUpperCase().trim();
	           data.put(pid, name);
	        }
	        
	        
	        reader.close();
	        return data;
		}
		
		public static ArrayList<String> getOrderedListCDT(File file)throws IOException{
			BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
	    	
	    	ArrayList list=new ArrayList();
	        String nextLine;
	        int count=0;
	        while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
	           
	        	if(count>2){
				String[] tokens=nextLine.split("\t");
				String key=tokens[1];
				//ArrayList list=new ArrayList();
				list.add(key);
	        	}
	        	
	        	count++;
	        }
	        
	        
	        reader.close();
	        return list;
		}
		
		
		public static void writeGCT(String save, Map<String, ArrayList> gct)throws IOException{
			FileWriter writer=new FileWriter(save);
			
			writer.write("#1.2\n");
			writer.write(gct.size()-1+"\t"+gct.get("header").size()+"\n");
			writer.write("PID\tName");
			
			ArrayList<String> header=gct.get("header");
			for(String str: header){writer.write("\t"+str);}
			writer.write("\n");
			
			for(String gene: gct.keySet()){
				if(!gene.equalsIgnoreCase("header")){
				writer.write(gene+"\t"+gene);
				ArrayList list=gct.get(gene);
				for(Object val: list){writer.write("\t"+val);}
				writer.write("\n");
				}
			}
			
			writer.close();
		}
		
		public static ArrayList parseCLS(File file)throws IOException{
			BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
	    	
	    	ArrayList list=new ArrayList();
	        String nextLine;
	        int count=0;
	        while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
	           
	        	if(count==2){
	        		String[] tokens=nextLine.split("\t");
	        		if(tokens.length==1){tokens=nextLine.split(" ");}
					for(int i=0; i<tokens.length; i++){list.add(tokens[i]);}
	        	}
	        	
	        	count++;
	        }
	        
	        
	        reader.close();
	        return list;
		}
		
		
		// Adding an optional line in CLS files that indicates the control sample
		public static String getControlClass(File file)throws IOException{
			String control=null;
			BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
	    	
	    	ArrayList list=new ArrayList();
	        String nextLine;
	        int count=0;
	        while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
	           
	        	if(count==3){
	        		control=nextLine;
	        	}
	        	
	        	count++;
	        }
	        
	        
	        reader.close();
	        return control;
		}
		
		public static Set parseCLSGroups(File file)throws IOException{
			BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
	    	
	    	Set list=new TreeSet();
	        String nextLine;
	        int count=0;
	        while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
	           
	        	if(count==2){
	        		String[] tokens=nextLine.split("\t");
	        		if(tokens.length==1){tokens=nextLine.split(" ");}
					for(int i=0; i<tokens.length; i++){list.add(tokens[i]);}
	        	}
	        	
	        	count++;
	        }
	        
	        
	        reader.close();
	        return list;
		}
		
		public static Map getGroupIndexes(ArrayList<String> groups){
			Map<String, ArrayList> rtrn=new TreeMap();
		
			int index=0;
			for(String group: groups){
				ArrayList temp=new ArrayList();
				if(rtrn.containsKey(group)){temp=rtrn.get(group);}
				temp.add(index);
				rtrn.put(group, temp);
				index++;
			}
			
			return rtrn;
		}
		
		public static Map<String, ArrayList> getGroupIndexes(File clsFile)throws IOException{
			Map<String, ArrayList> rtrn=new TreeMap();
			ArrayList<String> groups=parseCLS(clsFile);
			
			int index=0;
			for(String group: groups){
				ArrayList temp=new ArrayList();
				if(rtrn.containsKey(group)){temp=rtrn.get(group);}
				temp.add(index);
				rtrn.put(group, temp);
				index++;
			}
			
			return rtrn;
		}
		
		
		
		public static Map<Integer, String> getIndexGroup(ArrayList<String> groups){
			Map<Integer, String> rtrn=new TreeMap();
		
			int index=0;
			for(String group: groups){
				//ArrayList temp=new ArrayList();
				//if(rtrn.containsKey(group)){temp=rtrn.get(group);}
				//temp.add(index);
				rtrn.put(index, group);
				index++;
			}
			
			return rtrn;
		}


		public static void writeGMTFile(String save, Map<String, Collection<String>> gmtFile) throws IOException {
			FileWriter writer=new FileWriter(save);
			
			for(String name: gmtFile.keySet()){
				Collection<String> genes=gmtFile.get(name);
				writer.write(name+"\t"+genes.size());
				for(String gene: genes){writer.write("\t"+gene);}
				writer.write("\n");
			}
			
			writer.close();
		}


		public static MatrixWithHeaders parseNumericCls(File file, List<String> header) throws NumberFormatException, IOException {
            BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
	      	
	      	ArrayList<String> row= new ArrayList<String> ();
	      	//Map<String, ArrayList> data = new TreeMap();
	      	List<List<Double>> rawData = new ArrayList<List<Double>>();
	        String nextLine;
	        int count=1;
	        String key="";
	        int lineNum=1;
	        nextLine = reader.readLine();
	        while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
	            if (count%2==1){
	            	key=nextLine.substring(1);
	            	row.add(key);}
	        	if(count%2==0){
	        		String[] tokens=nextLine.split("\t");
	        		List<Double> lineData = new ArrayList<Double>(tokens.length);
	    			rawData.add(lineData);
	    			for(int i = 0 ; i < tokens.length; i++) {
	    				lineData.add(Double.parseDouble(tokens[i]));
	    			}
	    			
	    			lineNum++;
				}
	        	count++;
	        }
	        
	        reader.close();	        
	        MatrixWithHeaders mat=new MatrixWithHeaders(rawData,row,header);
	        
	        return mat;
		}


		public static Collection<String> parseGeneList(File file) throws IOException {
			BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
	      	
	      	Collection<String> rtrn= new TreeSet<String> ();
	      	String nextLine;
	        while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
	            rtrn.add(nextLine.trim().toUpperCase());
	        }
	        
	        reader.close();
			return rtrn;	    
		}
		
		public static Map<String, ArrayList<Double>> logTransform (Map<String, ArrayList<Double>> data, double addFactor){
			
			Map<String, ArrayList<Double>> newData=new HashMap<String, ArrayList<Double>>();
			for (String S:data.keySet())
			{
				ArrayList<Double> oldarr=data.get(S);
				ArrayList<Double> a=new ArrayList<Double>();
				for (int i=0 ; i<oldarr.size(); i++){
					double d=(double)oldarr.get(i)+addFactor;
					a.add(i,Math.log(d));
				}
				newData.put(S,a);
			}
			return newData;
		}
		
}