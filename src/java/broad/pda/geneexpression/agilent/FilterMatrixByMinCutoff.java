package broad.pda.geneexpression.agilent;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;

import broad.core.datastructures.MatrixWithHeaders;
import broad.core.error.ParseException;
import broad.core.math.Statistics;
import broad.pda.geneexpression.ExpressionExperimentInfo;

public class FilterMatrixByMinCutoff {

	public FilterMatrixByMinCutoff(MatrixWithHeaders data, String experimentInfoFile, String save, double cutoff) throws IOException{
		Map<String, Collection<String>> experimentInfo=AgilentUtils.parseExperimentInfoNonControls(experimentInfoFile);
		Collection<String> genes=new TreeSet<String>();
		
		Collection<String> controls=AgilentUtils.parseControlSamples(experimentInfoFile);
		for(String group: experimentInfo.keySet()){
			System.err.println(group);
			Collection<String> samples=experimentInfo.get(group);
			Collection<String> passingRows=this.getPassingRows(data, controls, samples, cutoff, true);
			genes.addAll(passingRows);
		}
		
		MatrixWithHeaders submatrix=data.submatrixByRowNames(genes);
		System.err.println("Original matrix has: "+data.getRowNames().size()+" New one: "+submatrix.getRowNames().size());
		
		submatrix.writeGCT(save);
	}
	
	private Collection<String> getPassingRows(MatrixWithHeaders submatrix, Collection<String> group1, Collection<String> group2, double cutoff, boolean preprocessByGroups) {
		Collection<String> rtrn=new TreeSet<String>();
		
		System.err.println(group1.size()+" "+group2.size());
		
		for(String row: submatrix.getRowNames()){
			//TODO Filter by fold change in groups
			if(preprocessByGroups){
				double[] gr1=new double[group1.size()];
				double[] gr2=new double[group2.size()];
				int i=0;
				for(String name: group1){
					gr1[i]=submatrix.get(row, name);
					i++;
				}
				i=0;
				for(String name: group2){
					gr2[i]=submatrix.get(row, name);
					i++;
				}
				double fold=Statistics.absFold(gr1, gr2, true);
				if(fold>cutoff){rtrn.add(row);}
			}
			else{
				double[] vals=submatrix.getRow(row);
				double fold=Statistics.fold(vals, true);
				if(fold>cutoff){rtrn.add(row);}
			}
		}
		
		return rtrn;
	}
	
	public static void main(String[] args)throws IOException, ParseException{
		if(args.length>3){
			MatrixWithHeaders data=new MatrixWithHeaders(args[0]);
			String experimentInfo=args[1];
			String save=args[2];
			double cutoff=new Double(args[3]);
			new FilterMatrixByMinCutoff(data, experimentInfo, save, cutoff);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=gct file \n args[1]=experiment info file \n args[2]=save \n args[3]=cutoff";
	
}
