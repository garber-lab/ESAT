package broad.core.math;

import java.util.*;

public class MannWhitneyPermutations {

	double p;

	
	public MannWhitneyPermutations(List<Double> x, List<Double> y, int numPerm){
		double score=computeScore(x,y);
		double[] randomScores=new double[numPerm+1];
		for(int i=0; i<numPerm; i++){
			randomScores[i]=shuffleProbesBetweenGroups(x,y);
		}
		randomScores[numPerm]=score;
		
		double[] minMax=minMax(randomScores);
		EmpiricalDistribution dist=new EmpiricalDistribution(200, minMax[0], minMax[1]);
		dist.addAll(randomScores);
		
		this.p=dist.getCumulativeProbability(score);
	}

	private double shuffleProbesBetweenGroups(List<Double> x, List<Double> y){
		List<Double> xShuffle=new ArrayList<Double>();
		List<Double> yShuffle=new ArrayList<Double>();
		
		for(int i=0; i<x.size(); i++){
			double rand=Math.random();
			if(rand<.5){
				yShuffle.add(x.get(i)); 
				xShuffle.add(y.get(i));//flip
			}else {
				yShuffle.add(y.get(i)); 
				xShuffle.add(x.get(i));
			}//retain
		}
		
		double score=computeScore(xShuffle, yShuffle);
		return score;
	}
	
	
	private double[] minMax(double[] vals){
		double [] rtrn =  {Statistics.min(vals), Statistics.max(vals)};
		return rtrn;
	}
	
	public double getP(){return this.p;}
	
	
	private double computeScore(List<Double> x, List<Double> y){
		//diff sum
		double sumX=Statistics.sum(x);
		double sumY=Statistics.sum(y);
		return sumX-sumY;
	}

	
}
