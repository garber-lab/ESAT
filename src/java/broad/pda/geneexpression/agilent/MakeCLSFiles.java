package broad.pda.geneexpression.agilent;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.core.datastructures.MatrixWithHeaders;
import broad.core.error.ParseException;
import broad.pda.geneexpression.ExpressionExperimentInfo;

public class MakeCLSFiles {

	public MakeCLSFiles(String gctFile, File experimentInfoFile, String save) throws IOException, ParseException{
		Map<String, ExpressionExperimentInfo> experimentInfo=AgilentUtils.parseExperimentInfoFile(experimentInfoFile);
		MatrixWithHeaders data=new MatrixWithHeaders(gctFile);
		writeCLSFileBySample(save, data, experimentInfo);
	}
	
	private void writeCLSFileBySample(String save, MatrixWithHeaders data, Map<String, ExpressionExperimentInfo> experimentInfo) throws IOException {
		Collection<String> columns=new ArrayList<String>();
		Collection<String> groups=new TreeSet<String>();
		Collection<String> controls=new TreeSet();
		for(String columnName : data.getColumnNames()) {
			ExpressionExperimentInfo info=experimentInfo.get(columnName);
			if(info.passedQC() && !info.isControl()){columns.add(columnName); groups.add(info.getExperimentName());}
			else if(info.passedQC() && info.isControl()){controls.add(info.getExperimentName());}
		}
		
		FileWriter writer=new FileWriter(save);
		
		Map<String, Integer> clsMap=new TreeMap<String, Integer>();
		
		writer.write(columns.size()+"\t"+(groups.size()+1)+"\t"+1+"\n");
		writer.write("#\tControl");
		int i=0;
		for(String group: groups){
			if(!clsMap.containsKey(group)){clsMap.put(group, i); i++;}
			writer.write("\t"+group);
		}
		writer.write("\n");
		
		for(String column: columns){
			ExpressionExperimentInfo info=experimentInfo.get(column);
			if(info.isControl()){writer.write("0\t");}
			else{writer.write((clsMap.get(info.getExperimentName())+1)+"\t");}
		}
		writer.write("\n");
		
		writer.close();
		
	}
	
	public static void main(String[] args)throws IOException, ParseException{
		if(args.length>2){
		String gctFile=args[0];
		File experimentFile=new File(args[1]);
		String save=args[2];
		new MakeCLSFiles(gctFile, experimentFile, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=gct file \n args[1]=experiment Info \n args[2]=save";
	
}
