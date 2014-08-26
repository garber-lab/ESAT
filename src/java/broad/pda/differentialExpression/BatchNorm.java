package broad.pda.differentialExpression;

import java.io.*;
import java.util.*;

import broad.core.datastructures.MatrixWithHeaders;
import broad.core.error.ParseException;
import broad.pda.geneexpression.ExpressionExperimentInfo;
import broad.pda.geneexpression.agilent.AgilentUtils;

public class BatchNorm {

	public BatchNorm(File gctFile, File batchDescription, String save) throws IOException, ParseException{
		MatrixWithHeaders data=new MatrixWithHeaders(gctFile.getAbsolutePath());
		Map<String, ExpressionExperimentInfo> info=AgilentUtils.parseExperimentInfoFile(batchDescription);
		Map<String, Collection<String>> infoByGroups=AgilentUtils.parseExperimentInfoFileToGroups(batchDescription);
		
		computeFoldByBatch(data, info, infoByGroups);
		data.writeGCT(save);
	}
	
	
	//TODO If there are no controls in the batch then simply skip that batch
	private static void computeFoldByBatch(MatrixWithHeaders data, Map<String, ExpressionExperimentInfo> fullInfo, Map<String, Collection<String>> experimentInfo) throws IOException {
		Map<String, Collection<String>> experimentsByBatch=new TreeMap<String, Collection<String>>();
		for(String group: experimentInfo.keySet()){
			Collection<String> experiments=experimentInfo.get(group);
			for(String experiment: experiments){
				String batch=fullInfo.get(experiment).getBatch();
				Collection<String> all=new TreeSet<String>();
				if(experimentsByBatch.containsKey(batch)){all=experimentsByBatch.get(batch);}
				all.add(experiment);
				experimentsByBatch.put(batch, all);
			}
		}
		
		for(String batch: experimentsByBatch.keySet()){
			MatrixWithHeaders batchMatrix=data.submatrixByColumnNames(experimentsByBatch.get(batch));
			
			Collection<String> controls=new TreeSet<String>();
			for(String experiment: batchMatrix.getColumnNames()){
				if(fullInfo.get(experiment).isControl()){controls.add(experiment);}
			}
			
			if(!controls.isEmpty()){
				for(String experiment: batchMatrix.getColumnNames()){
					for(String gene: batchMatrix.getRowNames()){
						double val=batchMatrix.get(gene, experiment);
						double zScore=fold(val, data.getValues(gene, controls));
						data.set(gene, experiment, zScore);
					}
				}
			}
			
			batchMatrix.writeGCT("batch"+batch+".gct");
		}
	}

	private static double fold(double val, double[] values) {
		return val-median(values);
	}
	
	public static double median(double[] list){
		Arrays.sort(list);
		if(list.length%2==0){
			double val1=list[list.length/2];
			double val2=list[(list.length/2)-1];
			double median=(val1+val2)/2.0;
			return median;
		}
		else{return list[list.length/2];}
	}
	
	public static void main(String[] args) throws IOException, ParseException{
		if(args.length>2){
			File file=new File(args[0]);
			File infoFile=new File(args[1]);
			String save=args[2];
			new BatchNorm(file, infoFile, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=data \n args[1]=description \n args[2]=save";
	
}
