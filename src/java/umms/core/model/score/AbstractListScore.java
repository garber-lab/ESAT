package umms.core.model.score;

import java.util.List;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.stat.descriptive.rank.*;
import org.apache.commons.math3.stat.descriptive.moment.*;
import org.apache.commons.math3.stat.descriptive.summary.*;

import umms.core.annotation.*;

public abstract class AbstractListScore extends WindowScore.AbstractWindowScore {

	public static double DEFAULT_VALUE = -1.0;
	
	double[] scores;

	public AbstractListScore(Annotation a) {
		super(a);
	}
	
	public void setScores(List<Double> scores) {
		this.scores = ArrayUtils.toPrimitive(scores.toArray(new Double[scores.size()]));
	}

	@Override
	public double getScore() {
		return getMean();
	}
	
	public double getCount() {
		return scores.length;
	}
	
	public double getSum() {
		return new Sum().evaluate(scores);
	}
	
	public double getMedian() {
		if (scores.length == 0) return DEFAULT_VALUE;
		return new Median().evaluate(scores);
	}
	
	public double getMean() {
		if (scores.length == 0) return DEFAULT_VALUE;
		return new Mean().evaluate(scores);
	}
	
	public double getStandardDeviation() {
		if (scores.length == 0) return DEFAULT_VALUE;
		return new StandardDeviation().evaluate(scores);
	}
	
	public double[] getScores() {
		return scores;
	}
	
	public String toString() {
		getAnnotation().setScore(getScore());
		return getAnnotation().toBED() + "\t" + getScore() + "\t" + getMean() + "\t" + getStandardDeviation() + "\t" + getMedian() + "\t" + getCount();
	}
}
