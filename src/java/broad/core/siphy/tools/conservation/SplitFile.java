package broad.core.siphy.tools.conservation;

import java.io.*;
import java.util.*;

public class SplitFile {

	public SplitFile(File file, String saveDir, int num)throws IOException{
		Vector<String> lines=readFile(file);
		writeSplit(saveDir, lines, num);
	}
	
	private Vector readFile(File file) throws IOException{
		FileInputStream fileInput;
		BufferedReader buf = null;
		Vector temp = new Vector(2000, 500);
		String aLine = new String("");
	    try{
	      fileInput = new FileInputStream(file);
	      buf = new BufferedReader(new InputStreamReader (fileInput));
	    } catch (IOException ex){
	       return temp;
	    }

	    while(true) {
	      try {
	        aLine = buf.readLine();
	        if(aLine == null) break;
	        temp.add(aLine);
	      } catch (IOException e) {
	        temp.removeAllElements();
	        return temp;
	      }
	    }
	    
	    return temp;
	  }
	
	private void writeSplit(String saveDir, Vector<String> lines, int num)throws IOException{
		int counter=0;
		
		FileWriter writer=null;
		for(String line: lines){
			if(counter % num ==0){int i=counter/num; if(writer!=null){writer.close();} writer=new FileWriter(saveDir+"/"+i+".txt");}
			writer.write(line+"\n");
			counter++;
		}
		
	}
	
	public static void main(String[] args)throws IOException{
		File file=new File(args[0]);
		String saveDir=args[1];
		int num=new Integer(args[2]);
		new SplitFile(file, saveDir, num);
	}
}
