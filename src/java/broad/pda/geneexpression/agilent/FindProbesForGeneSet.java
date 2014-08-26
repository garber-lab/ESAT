package broad.pda.geneexpression.agilent;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeSet;

import broad.core.util.ParseGCTFile;

public class FindProbesForGeneSet {

	
	public FindProbesForGeneSet(Map<String, String> geneSet, Map<String, String> probes, String save) throws IOException{
		Collection<String> specificProbes=new TreeSet<String>();
		
		for(String gene: geneSet.keySet()){
			Collection<String> probesForGene=findProbes(probes, gene);
			specificProbes.addAll(probesForGene);
		}
		
		write(save, specificProbes, probes);
	}

	private Collection<String> findProbes(Map<String, String> probes, String gene) {
		Collection<String> rtrn=new TreeSet<String>();
		
		for(String probe: probes.keySet()){
			String name=probes.get(probe);
			if(name.equalsIgnoreCase(gene)){rtrn.add(probe);}
		}
		
		return rtrn;
	}

	private void write(String save, Collection<String> specificProbes, Map<String, String> annotations) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String probe: specificProbes){writer.write(probe+"\t"+annotations.get(probe)+"\n");}
		
		writer.close();
	}
	
	public static void main(String[] args)throws IOException{
		if(args.length>2){
			Map<String, String> geneSet=ParseGCTFile.parseChipFile(new File(args[0]));
			Map<String, String> probes=ParseGCTFile.parseChipFile(new File(args[1]));
			String save=args[2];
			new FindProbesForGeneSet(geneSet, probes, save);
		}
		else{System.err.println(usage);}
	}
	
	private static String usage=" args[0]=gene set \n args[1]=probes \n args[2]=save";

	
}
