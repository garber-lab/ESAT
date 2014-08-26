package broad.pda.differentialExpression;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;

import broad.core.datastructures.MatrixWithHeaders;
import broad.core.math.Statistics;

public class DifferentialScoring {

	public static MatrixWithHeaders computeAbsFoldMatrix(MatrixWithHeaders data, Collection<String> group1, Collection<String> group2){
		
		return computeFoldMatrix(data, group1, group2, true);
	}
	
	public static MatrixWithHeaders computeFoldMatrix(MatrixWithHeaders data, Collection<String> group1, Collection<String> group2){
		return computeFoldMatrix(data, group1, group2, false);
	}
	
	private static MatrixWithHeaders computeFoldMatrix(MatrixWithHeaders data, Collection<String> group1, Collection<String> group2, boolean useAbsoluteFold){
		List<String> foldColumn=new ArrayList<String>();
		foldColumn.add("Fold");
		
		MatrixWithHeaders foldMatrix=new MatrixWithHeaders(data.getRowNames(), foldColumn);
		MatrixWithHeaders subset1=data.submatrixByColumnNames(group1);
		MatrixWithHeaders subset2=data.submatrixByColumnNames(group2);
		
		for(String gene: data.getRowNames()){
			double[] gr1=subset1.getRow(gene);
			double[] gr2=subset2.getRow(gene);
			double s= useAbsoluteFold ? DifferentialScoring.computeAbsFold(gr1, gr2) : DifferentialScoring.computeFold(gr1, gr2) ;
			foldMatrix.set(gene, "Fold", s);
		}
		
		return foldMatrix;
	}
	
	public static MatrixWithHeaders computeTestStatistics(MatrixWithHeaders data, Collection<String> group1, Collection<String> group2, double[] fudgeFactors, boolean useFold, boolean paired) {
		List<String> columns=new ArrayList<String>();
		for(int i=0; i<fudgeFactors.length; i++){columns.add(""+fudgeFactors[i]);}
		if(useFold){columns.add("Fold");}
				
		MatrixWithHeaders statisticMatrix=new MatrixWithHeaders(data.getRowNames(), columns);
		
		
		MatrixWithHeaders subset1=data.submatrixByColumnNames(group1);
		MatrixWithHeaders subset2=data.submatrixByColumnNames(group2);
		
		for(int i=0; i<fudgeFactors.length; i++){
			for(String gene: data.getRowNames()){
				double[] gr1=subset1.getRow(gene);
				double[] gr2=subset2.getRow(gene);				
				double s=DifferentialScoring.computeTestStatistic(gr1, gr2, fudgeFactors[i], paired);
				statisticMatrix.set(gene, i, s);
			}
		}
		
		if(useFold){
			for(String gene: data.getRowNames()){
				double[] gr1=subset1.getRow(gene);
				double[] gr2=subset2.getRow(gene);
				double s=DifferentialScoring.computeFold(gr1, gr2);
				statisticMatrix.set(gene, "Fold", s);
			}
		}
		
		
		return statisticMatrix;
	}
	
	public static MatrixWithHeaders computeTestStatistics(MatrixWithHeaders data, Map<String, Collection<String>> groups, double[] fudgeFactors, boolean useFold, boolean paired) {
		List<String> columns=new ArrayList<String>();
		for(int i=0; i<fudgeFactors.length; i++){columns.add(""+fudgeFactors[i]);}
		if(useFold){columns.add("Fold");}
				
		MatrixWithHeaders statisticMatrix=new MatrixWithHeaders(data.getRowNames(), columns);
		
		MatrixWithHeaders[] subset=new MatrixWithHeaders[groups.keySet().size()];
		int counter=0;
		for(String name: groups.keySet()){
			subset[counter++]=data.submatrixByColumnNames(groups.get(name));
		}
				
		for(int i=0; i<fudgeFactors.length; i++){
			for(String gene: data.getRowNames()){
				Collection<double[]> vals=new ArrayList<double[]>();
				for(int k=0; k<subset.length; k++){
					vals.add(subset[k].getRow(gene));
				}
				double s=DifferentialScoring.computeTestStatistic(vals, fudgeFactors[i], paired);
				statisticMatrix.set(gene, i, s);
				if(useFold){statisticMatrix.set(gene, "Fold", DifferentialScoring.computeFold(vals));}
			}
		}
		
		/*if(useFold){
			for(String gene: data.getRowNames()){
				double[] gr1=subset1.getRow(gene);
				double[] gr2=subset2.getRow(gene);
				double s=DifferentialScoring.computeFold(gr1, gr2);
				statisticMatrix.set(gene, "Fold", s);
			}
		}*/
		
		
		return statisticMatrix;
	}
	
	
	//Max fold between any 2 groups
	private static double computeFold(Collection<double[]> vals) {
		double maxFold=-1;
		for(double[] vals1: vals){
			for(double[] vals2: vals){
				double val=Statistics.absFold(vals1, vals2, true);
				if(val>maxFold){maxFold=val;}
			}
		}
		return maxFold;
	}

	private static double computeAbsFold(double[] gr1, double[] gr2) {
		double fold=Statistics.absFold(gr1, gr2, true);
		return fold;
	}
	
	private static double computeFold(double[] gr1, double[] gr2) {
		double fold=Statistics.fold(gr1, gr2, true);
		return fold;
	}

	private static double computeTestStatistic(double[] gr1, double[] gr2, double fudgeFactor, boolean paired) {
		Collection<double[]> groups=new ArrayList<double[]>();
		groups.add(gr1);
		groups.add(gr2);
		return computeTestStatistic(groups, fudgeFactor, paired);
	}
	
	private static double computeTestStatistic(Collection<double[]> groups, double fudgeFactor, boolean paired) {
		if(groups.size()>2){
			//do an ANOVA
			//TODO: add fudge factor to the computation
			double f=Statistics.anovaFStat(groups);
			return f;
		}
		else if(groups.size()==2){
			double[] gr1=(double[])groups.toArray()[0];
			double[] gr2=(double[]) groups.toArray()[1];
			if(gr1.length==gr2.length && paired){
				//TODO Compute paired t-stat
				double pairedT=Statistics.pairedTStat(gr1, gr2);
				return pairedT;
			}
			if(gr1.length==1 && gr2.length==1){
				//return fold change
				double fold=gr1[0]-gr2[0];
				return fold;
			}
			else if(gr1.length==1 || gr2.length==1){
				//do a z-test
				//System.err.println("Using Z-statistic");
				double z;
				if(gr1.length==1){z=Statistics.zScore(gr1[0], gr2, fudgeFactor);}
				else{z=Statistics.zScore(gr2[0], gr1, fudgeFactor);}
				return z;
			}
			else{
				//do a t-test
				double t=Statistics.tstat(gr1, gr2, fudgeFactor);
				return t;
			}
		}
		else{throw new IllegalArgumentException("Must be at least 2 groups to compute differences");}
		
	}
	
}
