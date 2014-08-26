package broad.pda.geneexpression.agilent;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.core.util.GMTParser;
import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;

public class ClassifyAgilentLincProbes {

	public ClassifyAgilentLincProbes(File lincProbesFile, File lincsByNameFile, String save) throws IOException{
		Map<String, String> lincProbes=GMTParser.parseCHIPFile(lincProbesFile);
		Map<String, IntervalTree<String>> lincsByName=parse(lincsByNameFile);
		
		write(save, lincProbes, lincsByName);
	}
	
	private Map<String, IntervalTree<String>> parse(File lincsByNameFile) throws IOException {
		Map<String, IntervalTree<String>> rtrn=new TreeMap<String, IntervalTree<String>>();
		Collection<String> lines=BEDFileParser.loadList(lincsByNameFile.getAbsolutePath(), true);
		
		for(String line: lines){
			String[] tokens=line.split("\t");
			Alignments align=new Alignments(tokens[3], tokens[4], tokens[5]);
			String name=tokens[0];
			IntervalTree<String> tree=new IntervalTree<String>();
			if(rtrn.containsKey(align.getChr())){
				tree=rtrn.get(align.getChr());
			}
			tree.put(align.getStart(), align.getEnd(), name);
			rtrn.put(align.getChr(), tree);
		}
		
		return rtrn;
	}

	private void write(String save, Map<String, String> lincProbes, Map<String, IntervalTree<String>> lincsByName) throws IOException{
		FileWriter writer=new FileWriter(save);
		
		for(String probe: lincProbes.keySet()){
			Alignments position=parsePosition(lincProbes.get(probe));
			if(position!=null){
				if(lincsByName.containsKey(position.getChr())){
				Iterator<Node<String>> overlappers=lincsByName.get(position.getChr()).overlappers(position.getStart(), position.getEnd());
				if(overlappers.hasNext()){
				writer.write(probe);
				Collection<String> set=new TreeSet<String>();
				while(overlappers.hasNext()){
					set.add(overlappers.next().getValue());
				}
				for(String val: set){writer.write("\t"+val);}
				writer.write("\n");
				}
			}else{System.err.println(position);}
			}
		}
		
		writer.close();
	}

	private Alignments parsePosition(String string) {
		try{
		String chr="chr"+string.split(":")[0].split("R")[1];
		String start=string.split(":")[1].split("-")[0];
		String end=string.split(":")[1].split("-")[1].split("_")[0];
		return new Alignments(chr, new Integer(start), new Integer(end));
		}catch(Exception ex){return null;}
	}
	
	
	public static void main(String[] args) throws IOException{
		if(args.length>3){
			File probe=new File(args[0]);
			File nanostring=new File(args[1]);
			String save=args[2];
			new ClassifyAgilentLincProbes(probe, nanostring, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=linc probes file (agilent chip file) \n args[1]=lincs by name (nanostring report file) \n args[2]=save";

}
