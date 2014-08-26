package broad.pda.geneexpression.agilent;

import java.io.File;
import java.io.IOException;
import java.util.Map;

import broad.core.datastructures.MatrixWithHeaders;
import broad.core.error.ParseException;
import broad.core.util.ParseGCTFile;
import broad.pda.geneexpression.ExpressionExperimentInfo;

public class ExtractSamplesAndLabel {

	public ExtractSamplesAndLabel(String gctFile, File info, String save, String chipFile) throws IOException, ParseException{
		MatrixWithHeaders data=new MatrixWithHeaders(gctFile);
		Map<String, String> experimentInfo=AgilentUtils.parseExperimentInfoFileToName(info, true);
		MatrixWithHeaders subset=data.submatrixByColumnNames(experimentInfo.keySet());
		
		Map<String, String> chip=ParseGCTFile.parseChipFile(new File(chipFile));
		subset.setPIDToName(chip);
		subset.writeGCT(save);
	}
	
	public static void main(String[] args)throws IOException, ParseException{
		if(args.length>3){
			String gctFile=args[0];
			File info=new File(args[1]);
			String save=args[2];
			String chipFile=args[3];
			new ExtractSamplesAndLabel(gctFile, info, save, chipFile);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=gct file \n args[1]=experiment info \n args[2]=save \n args[3]=chip file";
	
}
