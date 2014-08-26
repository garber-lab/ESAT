package broad.pda.geneexpression.clustering;

import broad.core.math.Statistics;

public class AbsolutePearsonFunction implements ClusterDistanceFunction{

	public double measure(double[] a, double[] b) {
		return 1 - Math.abs(Statistics.pearsonDistance(a, b));
	}

	public String getName(){return "AbsolutePearsonDistance";}
}
