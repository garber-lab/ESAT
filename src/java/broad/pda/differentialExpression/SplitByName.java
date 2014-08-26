package broad.pda.differentialExpression;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.core.datastructures.MatrixWithHeaders;
import broad.core.error.ParseException;
import broad.pda.geneexpression.ExpressionExperimentInfo;
import broad.pda.geneexpression.agilent.AgilentUtils;

public class SplitByName {

	public SplitByName(File gctFile, File experimentInfoFile, File chipFile, String save, boolean splitByUnderScore, boolean useControls) throws IOException, ParseException{
		MatrixWithHeaders data=new MatrixWithHeaders(gctFile.getAbsolutePath(), chipFile.getAbsolutePath());
		Map<String, ExpressionExperimentInfo> experimentInfo=AgilentUtils.parseExperimentInfoFile(experimentInfoFile);
		splitByName(data, gctFile, experimentInfo, save, splitByUnderScore, useControls);
	}

	private void splitByName(MatrixWithHeaders data, File gctFile, Map<String, ExpressionExperimentInfo> experimentInfo, String save, boolean splitByUnderScore, boolean useControls) throws IOException {
		Map<String, Collection<String>> samplesByName=new TreeMap<String, Collection<String>>();
		
		Collection<String> controls=new TreeSet();
		
		for(String barcode: data.getColumnNames()){
			if(experimentInfo.containsKey(barcode)){
				ExpressionExperimentInfo info=experimentInfo.get(barcode);
				String group=info.getSampleType();
				String name=info.getExperimentName();
				if(splitByUnderScore){name=info.getExperimentName().split("_")[0];}
				Collection<String> nameSamples=new ArrayList();
				if(samplesByName.containsKey(name)){nameSamples=samplesByName.get(name);}
				nameSamples.add(barcode);
				samplesByName.put(name, nameSamples);
				if(group.equalsIgnoreCase("Control")){controls.add(barcode);}
			}
		}
		
		for(String name: samplesByName.keySet()){
			System.err.println(name);
			Collection<String> samples=samplesByName.get(name);
			if(useControls){samples.addAll(controls);}
			MatrixWithHeaders sub=data.submatrixByColumnNames(samples);
			sub.setPIDToName(data.getPIDToName());
			sub.writeGCT(save+"."+name+".gct");
		}
		
	}
	
	public static void main(String[] args) throws IOException, ParseException{
		if(args.length>5){
			File gctFile=new File(args[0]);
			File experimentInfoFile=new File(args[1]);
			File chipFile=new File(args[2]);
			String save=args[3];
			boolean splitByUnderScore=new Boolean(args[4]);
			boolean useControls=new Boolean(args[5]);
			new SplitByName(gctFile, experimentInfoFile, chipFile, save, splitByUnderScore, useControls);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=gct file \n args[1]=experiment info \n args[2]=chip file \n args[3]=save \n args[4]=split by underscore \n args[5]=use controls";
}
