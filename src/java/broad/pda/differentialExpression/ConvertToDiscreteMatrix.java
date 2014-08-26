package broad.pda.differentialExpression;

import java.io.IOException;
import java.util.Collection;
import java.util.TreeSet;

import broad.core.datastructures.MatrixWithHeaders;
import broad.core.error.ParseException;

public class ConvertToDiscreteMatrix {
	double[] alphas={.5,.25,.05,.01};

	public ConvertToDiscreteMatrix(MatrixWithHeaders fdrMatrix, String save, double alpha) throws IOException{
		MatrixWithHeaders discreteMatrix=new MatrixWithHeaders(fdrMatrix.getRowNames(), fdrMatrix.getColumnNames(), fdrMatrix.getPIDToName());
		
		for(String row: fdrMatrix.getRowNames()){
			for(String column: fdrMatrix.getColumnNames()){
				double fdr=fdrMatrix.get(row, column);
				double val=0;
				if(Math.abs(fdr)<alphas[3]){
					if(fdr<0){val=-4;}
					else{val=4;}
				}
				else if(Math.abs(fdr)<alphas[2]){
					if(fdr<0){val=-3;}
					else{val=3;}
				}
				else if(Math.abs(fdr)<alphas[1]){
					if(fdr<0){val=-2;}
					else{val=2;}
				}
				else if(Math.abs(fdr)<alphas[0]){
					if(fdr<0){val=-1;}
					else{val=1;}
				}
				discreteMatrix.set(row, column, val);
			}
		}
		
		//TODO Filter genes that no perturbation affects
		discreteMatrix=filterZeros(discreteMatrix);
		discreteMatrix.writeGCT(save);
	}
	
	
	private MatrixWithHeaders filterZeros(MatrixWithHeaders discreteMatrix) {
		Collection<String> genes=new TreeSet<String>();
		for(String gene: discreteMatrix.getRowNames()){
			double[] vals=discreteMatrix.getRow(gene);
			int count=countNonZero(vals);
			if(count>0){genes.add(gene);}
		}
		return discreteMatrix.submatrixByRowNames(genes);
	}


	private int countNonZero(double[] vals) {
		int counter=0;
		
		for(int i=0; i<vals.length; i++){
			if(vals[i]!=0){counter++;}
		}
		
		return counter;
	}


	public static void main(String[] args) throws IOException, ParseException{
		if(args.length>2){
			MatrixWithHeaders fdrMatrix=new MatrixWithHeaders(args[0]);
			String save=args[1];
			double alpha=new Double(args[2]);
			new ConvertToDiscreteMatrix(fdrMatrix, save, alpha);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=fdrMatrix \n args[1]=save \n args[2]=alpha";
	
}
