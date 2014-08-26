package broad.core.siphy;

/**
 * Simple interface that represents a model fit
 * @author mgarber
 *
 */
public interface Fit {
	double getLogLikelihoodRatio();
	int getPosition();

}
