package broad.pda.geneexpression.agilent;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import java.util.TreeSet;

import broad.core.datastructures.MatrixWithHeaders;
import broad.core.error.ParseException;
import broad.core.math.Statistics;
import broad.core.util.GMTParser;
import broad.pda.geneexpression.ExpressionExperimentInfo;

public class CorrelationAgilentNanostring {

	public CorrelationAgilentNanostring(MatrixWithHeaders nanostringData, MatrixWithHeaders agilentData, File nanostringInfo, File agilentInfo, String save, File probeMappingFile) throws IOException{
		Map<String, String> probeMapping=GMTParser.parseCHIPFile(probeMappingFile);
		
		MatrixWithHeaders normAgilent=norm(agilentData, agilentInfo);
		MatrixWithHeaders normNanostring=norm(nanostringData, nanostringInfo);
		
		
		makeJointMatrix(save, normAgilent, normNanostring, probeMapping);
	}

	private MatrixWithHeaders norm(MatrixWithHeaders agilentData, File agilentInfo) throws IOException {
		
		Map<String, Collection<String>> experimentInfo=AgilentUtils.parseExperimentInfoFileToGroups(agilentInfo);
		
		ArrayList<String> list=new ArrayList<String>();
		for(String group: experimentInfo.keySet()){
			if(!group.equalsIgnoreCase("Control")){list.add(group);}
		}
		
		
		MatrixWithHeaders rtrn=new MatrixWithHeaders(agilentData.getRowNames(), list); //all rows and columns are just non-controls
		
		Collection<String> controls=experimentInfo.get("Control");
		MatrixWithHeaders controlData=agilentData.submatrixByColumnNames(controls);
		for(String group: list){
			MatrixWithHeaders groupData=agilentData.submatrixByColumnNames(experimentInfo.get(group));
			for(String gene: agilentData.getRowNames()){
				double t=computeTestStat(controlData, groupData, gene);
				rtrn.set(gene, group, t);
			}
		}
		return rtrn;
	}

	private double computeTestStat(MatrixWithHeaders controlData, MatrixWithHeaders groupData, String gene) {
		return Statistics.tstat(controlData.getRow(gene), groupData.getRow(gene), 0.1);
	}

	private void makeJointMatrix(String save, MatrixWithHeaders normAgilent, MatrixWithHeaders normNanostring, Map<String, String> probeMapping) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		//TODO add header
		
		for(String probe: probeMapping.keySet()){
			String name=probeMapping.get(probe);
			double[] agilentVals=normAgilent.getRow(probe);
			double[] nanostringVals=normNanostring.getRow(name);
			writer.write(probe);
			for(int i=0; i<agilentVals.length; i++){writer.write("\t"+agilentVals[i]);}
			writer.write(name);
			for(int i=0; i<nanostringVals.length; i++){writer.write("\t"+nanostringVals[i]);}
			writer.write("\n");
		}
		
		
		writer.close();
	}
	
	public static void main(String[] args) throws IOException, ParseException{
		if(args.length>5){
			MatrixWithHeaders nanostringData=new MatrixWithHeaders(args[0]);
			MatrixWithHeaders agilentData=new MatrixWithHeaders(args[1]);
			File nanostringInfo=new File(args[2]);
			File agilentInfo=new File(args[3]);
			String save=args[4];
			File probeMapping=new File(args[5]);
			new CorrelationAgilentNanostring(nanostringData, agilentData, nanostringInfo, agilentInfo, save, probeMapping);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=nanostring data \t args[1]=agilent data \t args[2]=nanostring Info \n args[3]=agilentInfo \n args[4]=save \n args[5]=probeMapping";
	
}
