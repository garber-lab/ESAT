package broad.pda.geneexpression.clustering;

public interface ClusterDistanceFunction {
	String getName();
	double measure(double [] a, double [] b);
}
