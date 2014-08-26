package broad.pda.differentialExpression;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import nextgen.core.annotation.Gene;

import broad.core.util.GMTParser;
import broad.pda.annotation.BEDFileParser;

public class MakePositionalGeneSets {

	public MakePositionalGeneSets(File bedFile, String save, int size) throws IOException{
		Map<String, Collection<Gene>> genesByChr=BEDFileParser.loadDataByChr(bedFile);
		
		Map<String, Collection<String>> geneSets=makeGeneSets(genesByChr, size);
		
		GMTParser.writeGMT(save, geneSets);
	}

	private Map<String, Collection<String>> makeGeneSets(Map<String, Collection<Gene>> genesByChr, int size) {
		Map<String, Collection<String>> rtrn=new TreeMap<String, Collection<String>>();
		
		//i will determine the initial position in the list to start from
		for(String chr: genesByChr.keySet()){
			Collection<Gene> genes=genesByChr.get(chr);
			ArrayList<Gene> geneList=new ArrayList<Gene>(genes);
			for(int i=0; i<geneList.size()-size; i++){
				Collection<Gene> set=new ArrayList<Gene>();
				for(int j=0; j<size; j++){
					set.add(geneList.get(i+j));
				}
				String name=name(set);
				Collection<String> nameSet=set(set);
				rtrn.put(name, nameSet);
			}
		}
		
		
		/*for(int i=0; i<size; i++){
			for(String chr: genesByChr.keySet()){
				Collection<RefSeqGene> genes=genesByChr.get(chr);
				ArrayList<RefSeqGene> geneList=new ArrayList<RefSeqGene>(genes);
				for(int index=i; index<geneList.size(); index+=size){
					Collection<RefSeqGene> set=new ArrayList<RefSeqGene>();
					for(int j=0; j<size; j++){
						try{set.add(geneList.get(index+j));}catch(IndexOutOfBoundsException ex){ex.printStackTrace();}
					}
					String name=name(set);
					Collection<String> nameSet=set(set);
					rtrn.put(name, nameSet);
				}
			}
		}*/
		
		return rtrn;
	}

	private Collection<String> set(Collection<Gene> set) {
		Collection<String> rtrn=new TreeSet<String>();
		
		for(Gene gene: set){
			rtrn.add(gene.getName().toUpperCase());
		}
		
		return rtrn;
	}

	private String name(Collection<Gene> set) {
		int min=set.iterator().next().getStart();
		int max=((Gene)set.toArray()[set.size()-1]).getEnd();
		String chr=set.iterator().next().getChr();
		return chr+":"+min+"-"+max;
	}
	
	public static void main(String[] args)throws IOException{
		if(args.length>2){
			File bedFile=new File(args[0]);
			String save=args[1];
			int size=new Integer(args[2]);
			new MakePositionalGeneSets(bedFile, save, size);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=BED file \n args[1]=save \n args[2]=size";
}
