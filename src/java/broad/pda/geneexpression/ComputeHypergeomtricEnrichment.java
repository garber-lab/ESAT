package broad.pda.geneexpression;

import java.io.*;
import java.util.*;

import broad.core.datastructures.MatrixWithHeaders;
import broad.core.error.ParseException;
import broad.core.math.FDRDistribution;
import broad.core.math.Statistics;
import broad.core.util.GMTParser;
import broad.core.util.ParseGCTFile;

public class ComputeHypergeomtricEnrichment {

	double alpha=.05;
	private int numRandom=100;
	//TODO Compute enrichment using two metrics
	//First, differentially expressed genes
	//Second, top N genes (Up and Down)
	
	public ComputeHypergeomtricEnrichment(MatrixWithHeaders diffExpressed, Map<String, Collection<String>> geneSets, Map<String, String> chipInfo, String save, double alpha, int cutoff) throws IOException{
		this.alpha=alpha;
		MatrixWithHeaders geneSetMatrix=new MatrixWithHeaders(new ArrayList(geneSets.keySet()), diffExpressed.getColumnNames());
		List<String> geneNames=new ArrayList<String>();
		for(String pid: chipInfo.keySet()){geneNames.add(chipInfo.get(pid));}
		
		//geneSets=filterGeneSets(geneSets, diffExpressed, chipInfo);
		//System.err.println("Filtered gene sets "+geneSets.size());
		
		for(String experiment: diffExpressed.getColumnNames()){
			System.err.println(experiment);
			Collection<String> up=new TreeSet();
			Collection<String> down=new TreeSet();
			for(String gene: diffExpressed.getRowNames()){
				double val=diffExpressed.get(gene, experiment);
				String geneName=chipInfo.get(gene);
				if(val>cutoff){up.add(geneName);}
				if(val<-cutoff){down.add(geneName);}
			}
			
			enrichment(up, down, geneSets, chipInfo, geneSetMatrix, experiment, geneNames);
			//write(save+"."+experiment+".gmt", enrichedGeneSets);
		}
	
		geneSetMatrix.writeGCT(save);
		
		
	}
	
	private Map<String, Collection<String>> filterGeneSets(Map<String, Collection<String>> geneSets, MatrixWithHeaders diffExpressed, Map<String, String> chipInfo) {
		Map<String, Collection<String>> rtrn=new TreeMap<String, Collection<String>>();
		
		for(String geneSet: geneSets.keySet()){
			int counter=count(geneSets.get(geneSet), diffExpressed, chipInfo);
			if(counter>0){rtrn.put(geneSet, geneSets.get(geneSet));}
		}
		
		return rtrn;
	}

	private int count(Collection<String> genes, MatrixWithHeaders diffExpressed, Map<String, String> chipInfo) {
		int counter=0;
		for(String gene: genes){
			if(diffExpressed.containsRow(gene)){
				int notZero=countNonZero(diffExpressed.getRow(gene));
				if(notZero>0){counter++;}
			}
		}
		return counter;
	}

	private int countNonZero(double[] row) {
		int counter=0;
		for(int i=0; i<row.length; i++){
			if(row[i]!=0){counter++;}
		}
		return counter;
	}

	private void write(String save, Map<String, Double> enrichedGeneSets) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String geneSet: enrichedGeneSets.keySet()){
			writer.write(geneSet+"\t"+enrichedGeneSets.get(geneSet)+"\n");
		}
		
		writer.close();
	}

	private void enrichment(Collection<String> up, Collection<String> down, Map<String, Collection<String>> geneSets, Map<String, String> chipInfo, MatrixWithHeaders geneSetMatrix, String experiment, List<String> genes) {
		int universe=chipInfo.size();
		//Map<String, Double> rtrn=new TreeMap();
		
		//ArrayList<String> perms=new ArrayList<String>();
		//for(int i=0; i<numRandom; i++){perms.add("P"+i);}
		
		/*MatrixWithHeaders random=new MatrixWithHeaders(new ArrayList(geneSets.keySet()), perms);
		for(int i=0; i<numRandom; i++){
			System.err.print(" "+i);
			Map<String, Collection<String>> randomGeneSets=generateRandom(geneSets, genes);
			Map<String, Double> randomPValues=computeEnrichment(randomGeneSets, universe, up, down);
			for(String geneSet: random.getRowNames()){
				double val=1;
				if(randomPValues.containsKey(geneSet)){val=randomPValues.get(geneSet);}
				random.set(geneSet, ("P"+i), val);
			}
		}
		System.err.println();*/
		
				
		Map<String, Double> vals=computeEnrichment(geneSets, universe, up, down);
		double[] observed=new double[geneSets.size()];
		int counter=0;
		for(String geneSet: geneSets.keySet()){
			double p=1;
			if(vals.containsKey(geneSet)){p=vals.get(geneSet);}
			/*if(p>0){p=Math.min(1.0, geneSets.size()*p);}
			else{
				p=Math.max(-1.0, geneSets.size()*p);
			}*/
			observed[counter++]=p;
			double log=0;
			if(p>0){log=-Math.log10(p);}
			if(p<0){log=Math.log10(-p);}
			geneSetMatrix.set(geneSet, experiment, log);
		}
		
		//TODO Compute FDR normalization
		/*FDRDistribution dist=new FDRDistribution(observed, random.getData(), alpha);
		
		//assign FDR
		for(String geneSet: geneSets.keySet()){
			double fdr=1;
			double val=1;
			if(vals.containsKey(geneSet)){
				val=vals.get(geneSet);
				fdr=dist.getFDR(val);
				if(val<0){fdr=-fdr;}
			}
			
		}*/
		
		//return rtrn;
	}
	
	private Map<String, Collection<String>> generateRandom(Map<String, Collection<String>> geneSets, List<String> genes) {
		Map<String, Collection<String>> rtrn=new TreeMap<String, Collection<String>>();
		
		for(String geneSet: geneSets.keySet()){
			Collection<String> permuted=permute(geneSets.get(geneSet), genes);
			rtrn.put(geneSet, permuted);
		}
		
		return rtrn;
	}
	
	private Collection<String> permute(Collection<String> geneSet, List<String> genes1){
		Collection<String> rtrn=new HashSet<String>();
						
		List<String> genes=new ArrayList<String>(genes1);
		for(int i=0; i<geneSet.size(); i++){
			int index=new Double(Math.random()*genes.size()).intValue();
			String randomGene=(String)genes.remove(index);
			rtrn.add(randomGene);
		}
		
		return rtrn;
	}

	private Map<String, Double> computeEnrichment(Map<String, Collection<String>> geneSets, int universe, Collection<String> up, Collection<String> down) {
		Map<String, Double> rtrn=new TreeMap<String, Double>();
		
		for(String geneSet: geneSets.keySet()){
			int set1=up.size();
			int set2=geneSets.get(geneSet).size();
			int overlap=overlap(up, geneSets.get(geneSet));
			double p=Statistics.hypergeometric(universe, set1, set2, overlap);
			//rtrn.put(geneSet, p);
			
			if(p<this.alpha){
				//if(p==0){p=.00000001;}
				rtrn.put(geneSet, p);
			}
			//else{geneSetMatrix.set(geneSet, experiment, 1);}
			//geneSetMatrix.set(geneSet, experiment, p);
			
			set1=down.size();
			overlap=overlap(down, geneSets.get(geneSet));
			p=Statistics.hypergeometric(universe, set1, set2, overlap);
			//rtrn.put(geneSet, -p);
			if(p<this.alpha){
				//if(p==0){p=.00000001;}
				rtrn.put(geneSet, -p);
			}
			
		}
		
		return rtrn;
	}

	public ComputeHypergeomtricEnrichment(Collection<String> set1, Collection<String> set2, Collection<String> universe) {
		int numberOfGenesTotal=universe.size();
		int diffExpressedGenes=set1.size();
		int genesInSet=set2.size();
		int diffExpressedGenesInSet=overlap(set1, set2);
		System.err.println(numberOfGenesTotal+" "+diffExpressedGenes+" "+genesInSet+" "+diffExpressedGenesInSet);
		//System.err.println(Statistics.hypergeometric(numberOfGenesTotal, diffExpressedGenes, genesInSet, diffExpressedGenesInSet));
	}

	private int overlap(Collection<String> set1, Collection<String> set2) {
		int overlap=0;
		
		for(String g1: set1){
			if(set2.contains(g1) || set2.contains(g1.toUpperCase())){overlap++;}
		}
		
		return overlap;
	}
	
	public static void main(String[] args)throws IOException, ParseException{
		if(args.length>5){
			MatrixWithHeaders diff=new MatrixWithHeaders(args[0]);
			Map<String, Collection<String>> geneSets=GMTParser.ParseGeneSetFile(new File(args[1]), 5);
			Map<String, String> chipInfo=GMTParser.parseCHIPFile(new File(args[2]));
			String save=args[3];
			double alpha=new Double(args[4]);
			int cutoff=new Integer(args[5]);
			new ComputeHypergeomtricEnrichment(diff, geneSets, chipInfo, save, alpha, cutoff);
		}
		else{
			System.err.println(usage);
		}
	}
	
	static String usage=" args[0]=GCT file (discretized FDR matrix) \n args[1]=gene sets \n args[2]=chip info \n args[3]=save \n args[4]=alpha \n args[5]=cutoff";
	
}
