package broad.pda.geneexpression.clustering;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.core.datastructures.MatrixWithHeaders;
import broad.core.error.ParseException;
import broad.core.math.Statistics;

public class RankOrder {

	public RankOrder(MatrixWithHeaders data, String save) throws IOException{
		//sort experiments by the number of hits
		Map<String, Integer> orderedExperiments=orderExperiments(data);
		
		//sort genes by the average rank of experiment hits
		//break ties by taking minimum
		Map<String, Integer> orderedGenes=orderGenes(orderedExperiments, data);
		
		//Reorder based on the new ranks
		data=reorderMatrix(data, orderedExperiments, orderedGenes);
		
		data.writeGCT(save);
	}
	
	private Map<String, Integer> orderGenes(Map<String, Integer> orderedExperiments, MatrixWithHeaders data) {
		//sort genes by the average rank of experiment hits
		//break ties by taking minimum
		
		
		
		Map<String, Integer> rtrn=new TreeMap<String, Integer>();
		
		Map<String, Double> tmp=new TreeMap<String, Double>();
		TreeSet<Double> set=new TreeSet<Double>();
		
		for(String gene: data.getRowNames()){
			//first rank genes by number of hits and then break ties with rank
			double[] vals=data.getRow(gene);
			int count=countZero(vals);
			double topRank=1000.0*count;
			
			double averageRank=averageRank(data, orderedExperiments, gene);
			double signRank=signRank(data, orderedExperiments, gene);
			//double minRank=minRank(data, orderedExperiments, gene);
			//double rank=topRank+averageRank+(signRank/10.0)+(minRank/1000.0);
			double rank=averageRank+(signRank);
			System.err.println(averageRank+" "+signRank+" "+rank);
			tmp.put(gene, rank);
			set.add(rank);
		}
		
		Map<Double, Integer> mapping=new TreeMap<Double, Integer>();
		
		int rank=0;
		for(Double num: set){
			mapping.put(num, rank);
			rank++;
		}
		
		for(String gene: tmp.keySet()){
			Double num=tmp.get(gene);
			int rankOrder=mapping.get(num);
			rtrn.put(gene, rankOrder);
		}
		
		
		return rtrn;
	}

	private double signRank(MatrixWithHeaders data,	Map<String, Integer> orderedExperiments, String gene) {
		for(String experiment: data.getColumnNames()){
			double val=data.get(gene, experiment);
			if(val==-4){return .4;}
			if(val==-3){return .5;}
			if(val==-2){return .6;}
			if(val==-1){return .7;}
			if(val==1){return .3;}
			if(val==2){return .2;}
			if(val==3){return .1;}
			if(val==4){return 0;}
		}
		return .9;
	}

	private int countZero(double[] vals) {
		int counter=0; 
		
		for(int i=0; i<vals.length; i++){
			if(vals[i]==0){counter++;}
		}
		
		return counter;
	}

	private double averageRank(MatrixWithHeaders data, Map<String, Integer> orderedExperiments, String gene) {
		Collection<Integer> ranks=new ArrayList<Integer>();
		for(String experiment: data.getColumnNames()){
			double val=data.get(gene, experiment);
			if(val!=0){ranks.add(orderedExperiments.get(experiment));}
		}
		return Statistics.average(ranks);
	}

	//TODO Need tie breakers, each experiment must have unique rank
	private Map<String, Integer> orderExperiments(MatrixWithHeaders data) {
		Map<String, Integer> rtrn=new TreeMap<String, Integer>();
		
		TreeSet<Integer> nums=new TreeSet<Integer>();
		
		for(String experiment: data.getColumnNames()){
			double[] vals=data.getColumn(experiment);
			int count=countNonZero(vals);
			nums.add(-count);
			rtrn.put(experiment, -count);
		}
		
		Map<Integer, Integer> mapping=new TreeMap<Integer, Integer>();
		
		int counter=0;
		for(Integer num: nums){
			mapping.put(num, counter);
			counter++;
		}
		
		Map<Integer, Collection<String>> test=new TreeMap<Integer, Collection<String>>();
		
		for(String experiment: rtrn.keySet()){
			Integer pos=rtrn.get(experiment);
			Collection<String> set=new TreeSet<String>();
			if(test.containsKey(pos)){set=test.get(pos);}
			set.add(experiment);
			test.put(pos, set);
		}
		
		Map<String, Integer> tmp=new TreeMap<String, Integer>();
		
		int count=0;
		for(Integer pos: test.keySet()){
			Collection<String> set=test.get(pos);
			for(String exp: set){
				tmp.put(exp, count);
				count++;
			}			
		}
		
		rtrn=tmp;
		return rtrn;
	}

	private int countNonZero(double[] vals) {
		int counter=0; 
		
		for(int i=0; i<vals.length; i++){
			if(vals[i]!=0){counter++;}
		}
		
		return counter;
	}

	private MatrixWithHeaders reorderMatrix(MatrixWithHeaders data, Map<String, Integer> orderedExperiments, Map<String, Integer> orderedGenes) {
		Map<Integer, Collection<String>> experiments=new TreeMap<Integer, Collection<String>>();
		Map<Integer, Collection<String>> genes=new TreeMap<Integer, Collection<String>>();
	
		for(String experiment: orderedExperiments.keySet()){
			Integer pos=orderedExperiments.get(experiment);
			Collection<String> set=new TreeSet<String>();
			if(experiments.containsKey(pos)){set=experiments.get(pos);}
			set.add(experiment);
			experiments.put(pos, set);
		}
		
		for(String gene: orderedGenes.keySet()){
			Integer pos=orderedGenes.get(gene);
			Collection<String> set=new TreeSet<String>();
			if(genes.containsKey(pos)){set=genes.get(pos);}
			set.add(gene);
			genes.put(pos, set);
		}
		
		ArrayList<String> geneList=new ArrayList<String>();
		for(Integer genePos: genes.keySet()){
			for(String gene: genes.get(genePos)){
				geneList.add(gene);
			}
		}
		
		ArrayList<String> experimentList=new ArrayList<String>();
		for(Integer experimentPos: experiments.keySet()){
			for(String experiment: experiments.get(experimentPos)){
				experimentList.add(experiment);
			}
		}
		
		MatrixWithHeaders rtrn=new MatrixWithHeaders(geneList, experimentList);
		for(String gene: data.getRowNames()){
			for(String experiment: data.getColumnNames()){
				rtrn.set(gene, experiment, data.get(gene, experiment));
			}
		}
		
		return rtrn;
	}

	public static void main(String[] args)throws IOException, ParseException{
		if(args.length>1){
			MatrixWithHeaders data=new MatrixWithHeaders(args[0]);
			String save=args[1];
			new RankOrder(data, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=data \n args[1]=save";
	
}
