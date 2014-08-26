package broad.pda.geneexpression.clustering;

import java.io.File;
import java.io.IOException;

import broad.core.datastructures.MatrixWithHeaders;
import broad.core.error.ParseException;

public class ComputeCorrelationOfReplicates {

	public ComputeCorrelationOfReplicates(File file, String save) throws IOException, ParseException{
		MatrixWithHeaders data=new MatrixWithHeaders(file.getAbsolutePath());
		
		for(String experiment: data.getColumnNames()){
			System.err.println(experiment);
		}
	}
	
	public static void main(String[] args) throws IOException, ParseException{
		File file=new File(args[0]);
		String save=args[1];
		new ComputeCorrelationOfReplicates(file, save);
	}
	
}
