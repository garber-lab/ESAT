package broad.core.math;

public class ScanStatistics {
	/**
	 * Params written by Jesse Aug 20, 2012 ... are these definitions correct?
	 * @param k			Observed count
	 * @param lambda	# reads on chromosome / # non-masked bases on chromosome
	 * @param w			window size
	 * @param T			# non-masked bases on chromosome
	 * @return
	 */
	public static double calculatePVal(int k, double lambda, double w, double T){
		if(k<=2){return 1;}
		double lambdaW=lambda*w;   // parameter for Poisson distribution
		double a=((k-lambdaW)/k)*(lambda*(T-w)*poisson(k-1, lambdaW));     // poisson function = Poisson PDF
		double result=Fp(k-1, lambdaW)*Math.exp(-a);					   // Fp = Poisson CDF
		double p=1-result;
		p=Math.abs(p);
		p=Math.min(1, p);
		//p=Math.max(0, p);
		return p;
	}
	
	public static double calculateApproximatePVal(int k, double lambda, double w, double T, double alpha){
		if(k<=2){return 1;}
		double lambdaW=lambda*w;
		double a=((k-lambdaW)/k)*(lambda*(T-w)*poisson(k-1, lambdaW));
		double p=FpWithBreaking(k-1, lambdaW, a, alpha);
		return p;
	}

	public static double poisson(int k, double lambda){
		cern.jet.random.Poisson poiss=new cern.jet.random.Poisson(lambda, new cern.jet.random.engine.DRand());
		return poiss.pdf(k);
	}
	
	public static double Fp(int k,double lambdaW){
		double sum=0;
		for(int i=0; i<=k; i++){
			sum+=poisson(i, lambdaW);
		}
		return sum;
	}
	
	public static double FpWithBreaking(int k,double lambdaW, double a, double stopPoint){
		double sum=0;
		for(int i=0; i<=k; i++){
			sum+=poisson(i, lambdaW);
			double result=sum*Math.exp(-a);
			double p=1-result;
			p=Math.abs(p);
			p=Math.min(1, p);
			if(p<stopPoint){return p;}
		}
		double result=sum*Math.exp(-a);
		double p=1-result;
		p=Math.abs(p);
		p=Math.min(1, p);
		return p;
	}
	
	
}
