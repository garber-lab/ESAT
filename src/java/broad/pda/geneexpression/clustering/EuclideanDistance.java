package broad.pda.geneexpression.clustering;

import broad.core.math.Statistics;

public class EuclideanDistance implements ClusterDistanceFunction {

	@Override
	public double measure(double[] a, double[] b) {
		return Statistics.euclideanDistance(a, b);
	}
	
	public String getName(){return "EuclideanDistance";}

}
