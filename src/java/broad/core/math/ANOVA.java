package broad.core.math;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

public class ANOVA {
	double F;
	double BGMS;
	double WGMS;
	double BGSS;
	double WGSS;
	int N;
	int a; //number of groups
	
	/*public ANOVA(Collection<List<Double>> data){
		this.a=data.size();
		this.N=computeN(data);
		this.BGSS=computeBGSS(data, N);
		this.WGSS=computeWGSS(data);
		
		this.BGMS=BGSS/(a-1);
		this.WGMS=WGSS/(N-a);
		this.F=BGMS/WGMS;
		
		
		
		for(List<Double> list: data){
			if(list.size()<2){this.F=0;}
		}
		
	}*/
	
	public ANOVA(Collection<double[]> data){
		this.a=data.size();
		this.N=computeN(data);
		this.BGSS=computeBGSS(data, N);
		this.WGSS=computeWGSS(data);
		
		this.BGMS=BGSS/(a-1);
		this.WGMS=WGSS/(N-a);
		this.F=BGMS/WGMS;
	}
	
	
	public double F(){return this.F;}
	
	
	/*private double computeWGSS(Collection<List<Double>> data){
		double rtrn=0;
		double subtract=0;
		for(List<Double> vals: data){
			rtrn+=sumSquares(vals);
			subtract+=(vals.size()*Math.pow(Statistics.average(vals),2));
		}
		rtrn-=subtract;
		return rtrn;
	}*/
	
	private double computeWGSS(Collection<double[]> data){
		double rtrn=0;
		double subtract=0;
		for(double[] vals: data){
			rtrn+=sumSquares(vals);
			subtract+=(vals.length*Math.pow(Statistics.average(vals),2));
		}
		rtrn-=subtract;
		return rtrn;
	}
	
	
	private double sumSquares(List<Double> vals){
		double rtrn=0;
		for(int i=0; i<vals.size(); i++){
			rtrn+= Math.pow(vals.get(i),2);
		}
		return rtrn;
	}
	
	private double sumSquares(double[] vals){
		double rtrn=0;
		for(int i=0; i<vals.length; i++){
			rtrn+= Math.pow(vals[i],2);
		}
		return rtrn;
	}
	
	/*private double computeBGSS(Collection<List<Double>> data, int N){
		double rtrn=0;
		for(List<Double> vals: data){
			rtrn+=(Math.pow(Statistics.average(vals),2)*vals.size());
		}
		rtrn-=(N*Math.pow(Statistics.metaAvg(data),2));
		return rtrn;
	}*/
	
	private double computeBGSS(Collection<double[]> data, int N){
		double rtrn=0;
		for(double[] vals: data){
			rtrn+=(Math.pow(Statistics.average(vals),2)*vals.length);
		}
		//System.err.println(rtrn);
		rtrn-=(N*Math.pow(Statistics.metaAvg(data),2));
		//System.err.println(Statistics.metaAvg(data)+ " "+N+" "+(N*Statistics.metaAvg(data)));
		return rtrn;
	}
	
	
	/*private int computeN(Collection<List<Double>> data){
		int rtrn=0;
		for(List<Double> vals: data){
			rtrn+=vals.size();
		}
		return rtrn;
	}*/
	
	private int computeN(Collection<double[]> data){
		int rtrn=0;
		for(double[] vals: data){
			rtrn+=vals.length;
		}
		return rtrn;
	}
	
	
	public static void main(String[] args) {
		double[] list1={5957.4,	4054.4,	3914.1,	3691.4,	4282.8,	4340.1,	5490.3,	6093.6};
		double[] list2={5417.8,	4399.3,	3380,	2816.9,	4280.6,	3389.6,	3454.3,	2468.5};
		double[] list3={5661.9,	4887.8,	4242.7,	5638.7,	3552.2,	4079.2,	2867.5};
		double[] list4={7375,	4357.1,	4041.7,	3818.2,	4013.3,	2624.4,	1663.6,	1989.8};
		double[] list5={6700,	3964.4,	3379.4,	2217.5,	2682.2,	2153.6,	3361.4,	3183.2};
		double[] list6={1462.7,	909,	655.6,	496.1,	536,	632.5,	1384.7,	1889.8};
		double[] list7={5759,	3772.4,	1332.2,	621.5,	495.8,	453.8,	232.9};
		
		List<double[]> list=new ArrayList<double []>();
		list.add(list1);
		list.add(list2);
		list.add(list3);
		list.add(list4);
		list.add(list5);
		list.add(list6);
		list.add(list7);
		ANOVA anova=new ANOVA(list);
		System.err.println("F: "+anova.F);
		System.err.println("BGSS: "+anova.BGSS+" BGMS: "+anova.BGMS);
		System.err.println("WGSS: "+anova.WGSS+" WGMS: "+anova.WGMS);
		System.err.println("N: "+anova.N+" a: "+anova.a);
		System.err.println("df: "+(anova.N-anova.a)+","+(anova.a-1));
		System.err.println("TSS: "+(anova.WGSS+anova.BGSS));
	}
}
