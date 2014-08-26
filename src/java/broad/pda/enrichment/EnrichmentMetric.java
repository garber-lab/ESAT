package broad.pda.enrichment;

import java.util.HashSet;
import java.util.Set;
//import net.sf.picard.metrics.MetricBase;


/**
 * All enrichment metrics should inherit this class.
 * Could use the PICARD base for potentially easy output
 * @author engreitz
 *
 */
public class EnrichmentMetric {

	//public enum DATA_TYPE { DISCRETE, CONTINUOUS };
	protected String name, dataType1, dataType2;
	protected String statisticName, method;
	protected double statistic, pvalue;
	protected int totalCount, regionCount;

	public EnrichmentMetric(String name, String dataType1, String dataType2, String statisticName, String method, double statistic, double pvalue, int totalCount, int regionCount) {
		this.name = name;
		this.dataType1 = dataType1;
		this.dataType2 = dataType2;
		this.statisticName = statisticName;
		this.method = method;
		this.statistic = statistic;
		this.pvalue = pvalue;
		this.totalCount = totalCount;
		this.regionCount = regionCount;
	}
	
	public Set<String> getDataTypes() { 
		java.util.Set<String> both = new HashSet<String>();
		both.add(dataType1);
		both.add(dataType2);
		return both;
	}
	
	public String getDataType1() { return dataType1; }
	public String getDataType2() { return dataType2; }
	
	public String getMethod() { return statisticName; }
	public String getSubmethod() { return method; }

	public double getStatistic() { return statistic; }
	public double getPValue() { return pvalue; }
	
	public String toString() {
		return name + "\t" + getId() + "\t" + dataType1 + "\t" + dataType2 + "\t" + statisticName + "\t" + method + "\t" + statistic + "\t" + pvalue + "\t" + totalCount + "\t" + regionCount;
	}
	
	public String getId() {
		return name + "." + statisticName + "." + method;
	}
}
