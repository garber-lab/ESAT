package broad.pda.geneexpression;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.TreeSet;

import broad.core.datastructures.MatrixWithHeaders;
import broad.core.error.ParseException;

public class FilterGCTFile {

	double floorVal=7;
	
	public FilterGCTFile(String gctFile, String save, double fold) throws IOException, ParseException{
		MatrixWithHeaders expression=new MatrixWithHeaders(gctFile);
		//MatrixWithHeaders filtered=expression.filterInvariantGenes(expression, fold);
		
		MatrixWithHeaders filtered=floorAndFilter(expression, floorVal);
		
		filtered.writeGCT(save);
	}
	
	private MatrixWithHeaders floorAndFilter(MatrixWithHeaders filtered,double floorVal2) {
		Collection<String> rtrn=new TreeSet<String>();
		
		for(String gene: filtered.getRowNames()){
			int count=countLessThan(filtered.getRow(gene), floorVal2);
			if(count>0){rtrn.add(gene);}
		}
		
		System.err.println("Filtered from "+filtered.getRowNames().size()+" to "+rtrn.size());
		return filtered.submatrixByRowNames(rtrn);
	}

	private int countLessThan(double[] row, double floorVal2) {
		int counter=0;
		for(int i=0; i<row.length; i++){
			if(row[i]<floorVal2){counter++;}
		}
		return counter;
	}

	public static void main(String[] args)throws IOException, ParseException{
		if(args.length>2){
			String gct=args[0];
			String save=args[1];
			double fold=new Double(args[2]);
			new FilterGCTFile(gct, save, fold);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=gct file \n args[1]=save \n args[2]=fold";
	
}
