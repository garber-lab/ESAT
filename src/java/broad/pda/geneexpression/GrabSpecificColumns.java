package broad.pda.geneexpression;

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

public class GrabSpecificColumns {

	public static void main(String[] args) throws IOException, ParseException{
		if(args.length>2){
			MatrixWithHeaders data=new MatrixWithHeaders(args[0]);
			String infoFile=args[1];
			Map<String, ExpressionExperimentInfo> info=AgilentUtils.parseExperimentInfoFile(new File(infoFile));
			//Map<String, Collection<String>> experimentInfo=AgilentUtils.parseExperimentInfoFileToGroups(new File(infoFile));
			//computeFoldByBatch(data, info, experimentInfo);
			//data.writeGCT(args[2]);
			
			data.submatrixByColumnNames(info.keySet()).writeGCT(args[2]);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=data \n args[1]=description \n args[2]=save";
	
	
	
}
