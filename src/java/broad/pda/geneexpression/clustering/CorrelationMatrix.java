package broad.pda.geneexpression.clustering;

import broad.core.datastructures.MatrixWithHeaders;
import broad.core.math.Statistics;

public class CorrelationMatrix {

	public CorrelationMatrix(MatrixWithHeaders data, String save){
		
		for(int i=0; i<data.columnDimension(); i++){
			for(int j=0; j<data.columnDimension(); j++){
				double[] val1=data.getColumn(i);
				double[] val2=data.getColumn(j);
				double cor=Statistics.pearsonDistance(val1, val2);
				System.err.println();
			}
		}
		
	}
	
}
