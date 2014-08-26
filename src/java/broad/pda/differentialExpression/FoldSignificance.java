package broad.pda.differentialExpression;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.core.datastructures.MatrixWithHeaders;
import broad.core.error.ParseException;
import broad.core.math.Statistics;
import broad.pda.geneexpression.agilent.AgilentUtils;

public class FoldSignificance {
	double alpha=.05;

	public FoldSignificance(MatrixWithHeaders data, MatrixWithHeaders fdrMatrix, Map<String, Collection<String>> groups, String save, double alpha) throws IOException{
		this.alpha=alpha;
		Map<String, String> experimentToGroup=experimentToGroup(groups);
		
		MatrixWithHeaders fold=convertToFold(data, groups);
		MatrixWithHeaders zscore=convertToZScore(data, groups);
		
		for(String gene: data.getRowNames()){
			//For each gene against the controls
			//compute ZScore
			//Compute fold
			
			for(String experiment: data.getColumnNames()){
				String group=experimentToGroup.get(experiment);
				double fdr=1;
				try{fdr=Math.abs(fdrMatrix.get(gene, group));}catch(NullPointerException ex){}
				if(fdr>alpha){
					fold.set(gene, experiment, 0); zscore.set(gene, experiment, 0);
				}
			}
			
		}
		
		fold=filter(fold);
		zscore=filter(zscore);
		
		fold.writeGCT(save+".fold.gct");
		zscore.writeGCT(save+".zscore.gct");
	}

	private MatrixWithHeaders filter(MatrixWithHeaders fold) {
		Collection<String> rtrn=new TreeSet<String>();
		
		for(String gene: fold.getRowNames()){
			int count=countNonZero(fold.getRow(gene));
			if(count>0){rtrn.add(gene);}
		}
		
		System.err.println(fold.getRowNames().size()+" "+rtrn.size());
		return fold.submatrixByRowNames(rtrn);
	}

	private int countNonZero(double[] row) {
		int counter=0;
		
		for(int i=0; i<row.length; i++){
			if(row[i]!=0){counter++;}
		}
		
		return counter;
	}

	private Map<String, String> experimentToGroup(Map<String, Collection<String>> groups) {
		Map<String, String> rtrn=new TreeMap<String, String>();
		
		for(String group: groups.keySet()){
			for(String sample: groups.get(group)){
				rtrn.put(sample, group);
			}
		}
		
		
		return rtrn;
	}

	private MatrixWithHeaders convertToFold(MatrixWithHeaders data, Map<String, Collection<String>> groups) {
		Collection<String> controls=groups.get("Control");
		
		MatrixWithHeaders rtrn=new MatrixWithHeaders(data.getRowNames(), data.getColumnNames());
		
		MatrixWithHeaders controlData=data.submatrixByColumnNames(controls);
		
		for(String group: groups.keySet()){
			MatrixWithHeaders sampleData=data.submatrixByColumnNames(groups.get(group));
			for(String gene: data.getRowNames()){
				for(String experiment: groups.get(group)){
					double fold=fold(sampleData.get(gene, experiment), controlData.getRow(gene));
					rtrn.set(gene, experiment, fold);
				}
			}
			
		}
		return rtrn;
	}
	
	private double fold(double d, double[] row) {
		return d-Statistics.median(row);
	}
	
	private double zscore(double d, double[] row) {
		return Statistics.zScore(d, row, .1);
	}

	
	private MatrixWithHeaders convertToZScore(MatrixWithHeaders data, Map<String, Collection<String>> groups) {
		Collection<String> controls=groups.get("Control");
		
		MatrixWithHeaders rtrn=new MatrixWithHeaders(data.getRowNames(), data.getColumnNames());
		
		MatrixWithHeaders controlData=data.submatrixByColumnNames(controls);
		
		for(String group: groups.keySet()){
			MatrixWithHeaders sampleData=data.submatrixByColumnNames(groups.get(group));
			for(String gene: data.getRowNames()){
				for(String experiment: groups.get(group)){
					double fold=zscore(sampleData.get(gene, experiment), controlData.getRow(gene));
					rtrn.set(gene, experiment, fold);
				}
			}
			
		}
		return rtrn;
	}
	
	public static void main(String[] args)throws IOException, ParseException{
		if(args.length>4){
			MatrixWithHeaders data=new MatrixWithHeaders(args[0]);
			MatrixWithHeaders fdr=new MatrixWithHeaders(args[1]);
			Map<String, Collection<String>> groups=AgilentUtils.parseExperimentInfoFileToGroups(new File(args[2]));
			String save=args[3];
			double alpha=new Double(args[4]);
			new FoldSignificance(data, fdr, groups, save, alpha);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=data \n args[1]=fdr \n args[2]=groups \n args[3]=save \n args[4]=alpha";
}
