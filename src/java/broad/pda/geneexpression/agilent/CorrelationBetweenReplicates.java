package broad.pda.geneexpression.agilent;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;

import broad.core.datastructures.MatrixWithHeaders;
import broad.core.error.ParseException;
import broad.core.math.Statistics;
import broad.pda.geneexpression.ExpressionExperimentInfo;

public class CorrelationBetweenReplicates {

	
	//Compute the correlation between replicates
	//compared to the correlation between random samples
	//compared to correlation between random samples from the same hyb day
	public CorrelationBetweenReplicates(MatrixWithHeaders data, File experimentInfoFile, String save) throws IOException{
		Map<String, Collection<ExpressionExperimentInfo>> experimentInfo=AgilentUtils.parseExperimentInfoFileToGroup(experimentInfoFile);
		
		FileWriter writer=new FileWriter(save);
		
		for(String group: experimentInfo.keySet()){
			System.err.println(group);
			//writer.write(group);
			Collection<ExpressionExperimentInfo> samples=experimentInfo.get(group);
			Object[] array=samples.toArray();
			for(int i=0; i<array.length; i++){
				for(int j=i; j<array.length; j++){
					if(i!=j){
						double corr=correlation(data, array[i], array[j]);
						writer.write(group+"\t"+((ExpressionExperimentInfo)array[i]).getExperimentBarcode()+"\t"+((ExpressionExperimentInfo)array[j]).getExperimentBarcode()+"\t"+corr+"\n");		
						System.err.println(group+"\t"+((ExpressionExperimentInfo)array[i]).getExperimentBarcode()+"\t"+((ExpressionExperimentInfo)array[j]).getExperimentBarcode()+"\t"+corr);
					}
				}
			}
		}
		writer.close();
	}

	private double correlation(MatrixWithHeaders data, Object o1, Object o2) {
		ExpressionExperimentInfo s1=(ExpressionExperimentInfo) o1;
		ExpressionExperimentInfo s2=(ExpressionExperimentInfo) o2;
		
		double[] val1=data.getColumn(s1.getExperimentBarcode());
		double[] val2=data.getColumn(s2.getExperimentBarcode());
	
		return Statistics.pearsonDistance(val1, val2);
	}
	
	public static void main(String[] args)throws IOException, ParseException{
		if(args.length>2){
			MatrixWithHeaders data=new MatrixWithHeaders(args[0]);
			File experimentInfo=new File(args[1]);
			String save=args[2];
			new CorrelationBetweenReplicates(data, experimentInfo, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=data \n args[1]=experiment info \n args[2]=save";
	
}
