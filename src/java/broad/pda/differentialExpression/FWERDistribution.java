package broad.pda.differentialExpression;

import java.util.ArrayList;
import java.util.Collection;

import Jama.Matrix;
import broad.core.math.EmpiricalDistribution;
import broad.core.math.Statistics;

public class FWERDistribution {

	public EmpiricalDistribution minDist;
	public EmpiricalDistribution maxDist;
	
	public FWERDistribution(double[] observed, Matrix permutations){
		
		//TODO: Add observed to distribution
		
		Collection<Double> neg=new ArrayList<Double>();
		Collection<Double> pos=new ArrayList<Double>();
		
		for(int i=0; i<permutations.getColumnDimension(); i++){
			double[] vals=permutations.getColumn(i);
			double min=Statistics.min(vals);
			double max=Statistics.max(vals);
			neg.add(min);
			pos.add(max);
		}
		
		minDist=new EmpiricalDistribution(neg,200);
		maxDist=new EmpiricalDistribution(pos, 200);	
	}

	public double getFWER(double d) {
		double fwer = 0;
		if(d<0){
			fwer =  minDist.getCumulativeProbability(d);
		}
		else{
			fwer =  1-maxDist.getCumulativeProbability(d);
		}
		
		return fwer;
	}
	
}
