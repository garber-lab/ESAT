package broad.pda.differentialExpression;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

import broad.core.datastructures.Pair;
import broad.core.util.GMTParser;
import broad.pda.annotation.BEDFileParser;

public class GrabProbes {

	public GrabProbes(File chipFile, File neighborFile, String save) throws IOException{
		Map<String, Pair<String>> neighbors=getNeighbors(neighborFile);
		Map<String, Collection<String>> chipInfo=GMTParser.parseCHIPFileByName(chipFile);
		
		FileWriter writer=new FileWriter(save);
		for(String linc: neighbors.keySet()){
			Pair<String> genes=neighbors.get(linc);
			Collection<String> left=chipInfo.get(genes.getValue1().toUpperCase());
			if(left==null){left=chipInfo.get(genes.getValue1());}
			Collection<String> right=chipInfo.get(genes.getValue2().toUpperCase());
			if(right==null){right=chipInfo.get(genes.getValue2());}
			if(left!=null){
				for(String leftGene: left){
					writer.write(linc+"\t"+leftGene+"\t"+genes.getValue1()+"\tLeft\n");
				}
			}
			else{System.err.println("Skipped "+genes.getValue1());}
			if(right!=null){
				for(String rightGene: right){
					writer.write(linc+"\t"+rightGene+"\t"+genes.getValue2()+"\tRight\n");
				}
			}
			else{System.err.println("Skipped "+genes.getValue2());}
		}
		writer.close();
	}
	
	private Map<String, Pair<String>> getNeighbors(File neighborFile) throws IOException {
		Collection<String> lines=BEDFileParser.loadList(neighborFile.getAbsolutePath(), true);
		Map<String, Pair<String>> rtrn=new TreeMap<String, Pair<String>>();
		
		for(String line: lines){
			String[] tokens=line.split("\t");
			String name=tokens[1];
			Pair<String> genes=new Pair<String>(tokens[2], tokens[3]);
			rtrn.put(name, genes);
		}
		
		return rtrn;
	}
	
	
	public static void main(String[] args) throws IOException{
		if(args.length>2){
		File chipfile=new File(args[0]);
		File neighborFile=new File(args[1]);
		String save=args[2];
		new GrabProbes(chipfile, neighborFile, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=chipFile \n args[1]=neighbor file \n args[2]=save";
	
}
