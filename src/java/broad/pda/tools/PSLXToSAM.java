package broad.pda.tools;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Writer;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import nextgen.core.annotation.Gene;

import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;

public class PSLXToSAM {
	Pattern hitPattern = Pattern.compile("^[0-9]+.*"); //I am really not sure why .* is required to match desired lines, but it is 
	public void parseAndWrite(InputStream is, Writer w)throws IOException{
		BufferedReader reader=new BufferedReader(new InputStreamReader(is));
		String nextLine;
        while ((nextLine = reader.readLine()) != null) {
        	Matcher m = hitPattern.matcher(nextLine);
        	if(m.matches()) {	
        		Gene gene=new Gene(nextLine);    	 
        		w.write(gene.toSAM()+"\n");
        	}
        }
 	}
	
	public void parseAndWrite(File file, String save)throws IOException{
		InputStream fis = new FileInputStream(file);
		FileWriter fwriter=new FileWriter(save);
		try {
			parseAndWrite(fis, fwriter);
		} finally {    
            fis.close();
            fwriter.close();
		}
	}
	
	
	/*private void parseAndWrite(File file, File fastaFile, String save)throws IOException{
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		FastaSequenceIO faReader=new FastaSequenceIO(fastaFile);
		Set<Utils.RefSeqGene> alignments=new TreeSet<Utils.RefSeqGene>();
		Set<String> names=new TreeSet();
		
		String nextLine;
    	FileWriter writer=new FileWriter(save);
        while ((nextLine = reader.readLine()) != null) {
          try{
        	  Utils.RefSeqGene gene=new Utils.RefSeqGene(nextLine);
        	  //gene.addSequence(faReader.extractRecord(gene.getName()));
        	 
        	  alignments.add(gene);
        	  names.add(gene.getName());
        	  
          }catch(Exception ex){System.err.println(nextLine);}
        	
         }
        
        Map<String, Sequence> sequences=faReader.getMap(names);
		for(Utils.RefSeqGene gene: alignments){
			//System.err.println(gene);
			gene.addSequence(sequences.get(gene.getName()));
			writer.write(gene.toSAM()+"\n");
			
		}
        
            writer.close();
            reader.close();
	}*/
	
	public static void main(String[] args)throws IOException{
		ArgumentMap argMap = CLUtil.getParameters(args, usage);
		
		InputStream is = argMap.getInputStream();
		Writer ow = argMap.getOutputWriter();
		PSLXToSAM parser = new PSLXToSAM();
		try {
			parser.parseAndWrite(is, ow);
		} finally {
			is.close();
			ow.close();
		}
		
	}
	
	static String usage=" args[0] -in <psl file (or just standard input)> -out <save file (or standard output if none is specified)>\n";
}
