package broad.core.math;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.*;


public class ComputeFDR {
	
	static int numBins=1000;
	
	public static double FDR(EmpiricalDistribution observed, EmpiricalDistribution[] random, double k){
		double Rk=observed.getNumberOfValuesLargerThan(k);
		double[] Vk=new double[random.length];
		for(int i=0; i<random.length; i++){Vk[i]=random[i].getNumberOfValuesLargerThan(k);}
		double EVk=expected(Vk);
		double FDRk=EVk/Rk;
		if(Rk==0){FDRk=EVk; }
		if(FDRk>1){FDRk=1;}
		//System.err.println("EVk: " + EVk + " Rk: " + Rk + " k: " + k + " fdr " + FDRk);
		return FDRk;
	}
	
	public static double FDR(Map<String, Double> observed, EmpiricalDistribution[] random, double k){
		double Rk=countNumGreaterThanThresh(observed, k);
		double[] Vk=new double[random.length];
		for(int i=0; i<random.length; i++){Vk[i]=random[i].getNumberOfValuesLargerThan(k);}
		double EVk=expected(Vk);
		double FDRk=EVk/Rk;
		if(Rk==0){FDRk=0; }
		if(FDRk>1){FDRk=1;}
		//System.err.println("EVk: " + EVk + " Rk: " + Rk + " k: " + k + " fdr " + FDRk);
		return FDRk;
	}
	
	
	public static double FDRExact(Map observed, Map[] expected, double k){
		double Rk=countNumGreaterThanThresh(observed, k);
		double[] Vk=countNumGreaterThanThresh(expected, k);
		double EVk=expected(Vk);
		double FDRk=EVk/Rk;
		return FDRk;		
	}
	
	public static double FDRExactPrint(Map observed, Map[] expected, double k){
		EmpiricalDistribution[] randomTDist=makeDistributions(observed, expected);
		double RkEstimate=randomTDist[0].getNumberOfValuesLargerThan(Math.abs(k));
		double[] VkVector=new double[randomTDist.length];
		for(int i=0; i<randomTDist.length; i++){VkVector[i]=randomTDist[i].getNumberOfValuesLargerThan(Math.abs(k));}
		double VkEstimate=expected(VkVector);
		
		double Rk=countNumGreaterThanThresh(observed, k);
		double[] Vk=countNumGreaterThanThresh(expected, k);
		double EVk=expected(Vk);
		
		double FDRkEstimate=VkEstimate/RkEstimate;
		double FDRk=EVk/Rk;
		
		System.err.println(k+" "+FDRkEstimate+" "+FDRk+" "+RkEstimate+" "+VkEstimate+" "+Rk+" "+EVk);
		
		//double FDRk=EVk/Rk;
		return FDRk;		
	}
	
	private static EmpiricalDistribution[] makeDistributions(Map<String, Double> tstats, Map[] randomTStats){
		EmpiricalDistribution tDist=new EmpiricalDistribution(tstats.values(), numBins, true);
		EmpiricalDistribution[] randomTDist=new EmpiricalDistribution[randomTStats.length+1];
		randomTDist[0]=tDist;
		for(int i=0; i<randomTStats.length; i++){
			randomTDist[i+1]=new EmpiricalDistribution(randomTStats[i].values(), numBins, true);
		}
		return randomTDist;
	}
	
	private static double[] countNumGreaterThanThresh(Map<String, Double>[] tstats, double k){
		double[] counts=new double[tstats.length];
		k=Math.abs(k);
		for(int i=0; i<counts.length; i++){
			for(String gene: tstats[i].keySet()){
				double t=Math.abs(tstats[i].get(gene));
				if(t>k){counts[i]=counts[i]+1;}
			}
		}
		return counts;
	} 
	
	private static double countNumGreaterThanThresh(Map<String, Double> tstats, double k){
		double count=0.0;
		k=Math.abs(k);
		for(String gene: tstats.keySet()){
			double t=Math.abs(tstats.get(gene));
			if(t>=k){count=count+1;}
		}
		return count;
	} 
	
	private static double expected(double[] values){
		double avg=0;
		for(int i=0; i<values.length; i++){avg+=values[i];}
		return avg/values.length;
	}

}
