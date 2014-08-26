package broad.projection.gsa;

import java.io.File;
import java.io.IOException;
import java.util.Map;

import broad.core.datastructures.MatrixWithHeaders;
import broad.core.error.ParseException;
import broad.core.util.ParseGCTFile;

public class ExtractGenes {

	
	public static void main(String[] args)throws IOException, ParseException{
		if(args.length>2){
		MatrixWithHeaders gct=new MatrixWithHeaders((args[0]));
		Map<String, String> geneSet=ParseGCTFile.parseChipFile(new File(args[1]));
		String save=args[2];
		
		gct.submatrixByRowNames(geneSet.keySet()).writeGCT(save);
		}
		else{System.err.println(usage);}
		
	}
	static String usage="ExtractGenes args[0]=expression file \n args[1]=chip probe-name key file \n args[2]=save file \n";
	
}
