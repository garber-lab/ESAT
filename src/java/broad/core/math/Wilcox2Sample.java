package broad.core.math;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

public class Wilcox2Sample {

	
	double Z;
	double p;
	
	public Wilcox2Sample(double[] x, double[] y){
		double[] d=new double[x.length];
		double[] dAbs=new double[x.length];
		TreeMap<Double,List<String>> map=new TreeMap<Double,List<String>>();
		for(int i=0; i<d.length; i++){
			d[i]=x[i]-y[i];
			dAbs[i]=Math.abs(x[i]-y[i]);
			String group= d[i] <= 0 ? "x" : "y";
			List<String>  list = map.get(Math.abs(d[i]));
			if(list == null) {
				list=new ArrayList<String>();
			}
			list.add(group);
			map.put(Math.abs(d[i]), list);
		}
		
		List<Double>[] xyRanked=MannWhitney.rankOrder(map);
		this.Z=calculateZ(xyRanked, d.length);
		this.p=calculateP(Z);
	}
	
	public Wilcox2Sample(List<Double> xList, List<Double> yList){
		
		
		List<Double> xList2=new ArrayList<Double>();
		List<Double> yList2=new ArrayList<Double>();
		
		for(int i=0; i<xList.size(); i++){
			if(!Double.isInfinite(xList.get(i)) &&  !Double.isInfinite(yList.get(i))){
				xList2.add(xList.get(i)); yList2.add(yList.get(i));
			}
		}
		
		
		double[] x=new double[xList2.size()];
		double[] y=new double[yList2.size()];
		
		for(int i=0; i<x.length; i++){x[i] = xList2.get(i); y[i] = yList2.get(i);}
		
		this.p=new Wilcox2Sample(x,y).getP();
	}
	
	private double calculateP(double Z){
		cern.jet.random.Normal norm=new cern.jet.random.Normal(0,1, new cern.jet.random.engine.DRand());
		double cdf=norm.cdf(Z);
		return Math.min(1, Math.min((1-cdf), cdf)*2);
	}
	
	public double getZ(){return this.Z;}
	public double getP(){return this.p;}
	
	private double calculateZ(List<Double>[] xyRanked, int n){
		double U=Statistics.sum(xyRanked[0]);
		double mu=(n*(n+1))/4.0;
		double var=(n*(n+1)*((2*n)+1))/24.0;
		double sigma=Math.sqrt(var);
		double Z=(U-mu)/sigma;
		return Z;
	}
	
}
