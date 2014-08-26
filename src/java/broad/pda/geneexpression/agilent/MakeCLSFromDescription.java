package broad.pda.geneexpression.agilent;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.core.datastructures.MatrixWithHeaders;
import broad.core.error.ParseException;
import broad.pda.geneexpression.ExpressionExperimentInfo;

public class MakeCLSFromDescription {

	
	public MakeCLSFromDescription(File gctFile, File descriptionFile, String save) throws IOException, ParseException {
		Map<String, ExpressionExperimentInfo> experimentInfo=AgilentUtils.parseExperimentInfoFile(descriptionFile);
		MatrixWithHeaders data=new MatrixWithHeaders(gctFile.getAbsolutePath());
		data.writeGCTAndCLS(save, experimentInfo);
	}
	
	
	public void writeCLS(String save, Map<String, ExpressionExperimentInfo> experimentInfo, File gctFile) throws IOException, ParseException{
		MatrixWithHeaders data=new MatrixWithHeaders(gctFile.getAbsolutePath());
		FileWriter writer=new FileWriter(save);
		
		Collection<String> passedQC=new TreeSet<String>();
		Map<String, Collection<String>> groups=new TreeMap<String, Collection<String>>();
		
		for(String barcode: experimentInfo.keySet()){
			ExpressionExperimentInfo info=experimentInfo.get(barcode);
			if(info.passedQC()){
				passedQC.add(barcode);
				if(info.isControl()){
					Collection<String> samples=new TreeSet<String>();
					if(groups.containsKey("Controls")){samples=groups.get("Controls");}
					samples.add(barcode);
					groups.put("Controls", samples);
				}
				else{
					Collection<String> samples=new TreeSet<String>();
					if(groups.containsKey(info.getExperimentName())){samples=groups.get(info.getExperimentName());}
					samples.add(barcode);
					groups.put(info.getExperimentName(), samples);
				}
			}	
		}
		
		
		writer.write(passedQC.size()+"\t"+groups.keySet().size()+"\t1\n");
		writer.write("#\tControls");
		for(String group: groups.keySet()){
			writer.write("\t"+group);
		}
		writer.write("\n");
		
		for(String barcode: data.getColumnNames()){
			ExpressionExperimentInfo info=experimentInfo.get(barcode);
			if(info.passedQC()){
				if(info.isControl()){writer.write("Controls");}
				else{writer.write(info.getExperimentName());}
				writer.write("\t");
			}
		}
		writer.write("\n");
		
		writer.close();
	}

	public static void main(String[] args) throws IOException, ParseException{
		if(args.length>2){
		File gctFile=new File(args[0]);
		File descriptionFile=new File(args[1]);
		String save=args[2];
		new MakeCLSFromDescription(gctFile, descriptionFile, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=gct file \n args[1]=description file \n args[2]=save";
}
