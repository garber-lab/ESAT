package broad.pda.geneexpression.clustering;

import broad.core.math.Statistics;

public class PearsonDistance implements ClusterDistanceFunction {

	public double measure(double[] a, double[] b) {
		return 1-Statistics.pearsonDistance(a, b);
	}
	
	public String getName(){return "PearsonDistance";}

}
