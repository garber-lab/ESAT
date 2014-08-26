package broad.pda.geneexpression;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map;

import broad.pda.geneexpression.agilent.AgilentArrayStats;
import broad.pda.geneexpression.agilent.AgilentUtils;
import broad.pda.geneexpression.agilent.QCAgilentResults;

public class GlobalStats {

	public GlobalStats(File[] files, String save, File experimentInfoFile) throws IOException{
		Map<String, ExpressionExperimentInfo> experimentInfo=AgilentUtils.parseExperimentInfoFile(experimentInfoFile);
		FileWriter writer=new FileWriter(save);
		for(int i=0; i<files.length; i++){
			AgilentArrayStats stats =AgilentUtils.parseAgilentExperimentGlobalStats(files[i]);
			
			double floor=Math.pow(10, stats.getSpikeInDetectionLimit());
			double rtrn=Math.log(floor/stats.getNormalizationFactor())/Math.log(2);
			writer.write(files[i].getName()+"\t"+experimentInfo.get(files[i].getName().split("\\.")[0]).getExperimentHybDate()+"\t"+stats.getNonControl99PercentileSignal()+"\t"+stats.getNonControl50PercentileSignal()+"\t"+rtrn+"\t"+stats.getSpikeInDetectionLimit()+"\t"+stats.getSaturationSignalValue()+"\t"+QCAgilentResults.validate(stats)+"\t"+experimentInfo.get(files[i].getName().split("\\.")[0]).passedQC+"\n");
		}
		writer.close();
	}
	
	public static void main(String[] args)throws IOException{
		if(args.length>2){
			File[] files=new File(args[0]).listFiles();
			String save=args[1];
			File experimentInfo=new File(args[2]);
			new GlobalStats(files, save, experimentInfo);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=files \n args[1]=save \n args[2]=experiment Info";
	
}
