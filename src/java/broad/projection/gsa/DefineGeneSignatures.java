package broad.projection.gsa;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import broad.core.datastructures.MatrixWithHeaders;
import broad.core.error.ParseException;
import broad.core.math.Statistics;
import broad.core.util.ParseGCTFile;
import broad.pda.geneexpression.ScoreDiffGenes;

public class DefineGeneSignatures {
	
	boolean isLog=true;

	//Compute Differentially expressed genes
	//Also use top n genes
	public DefineGeneSignatures(File gctFile, File clsFile, String save, double fold) throws IOException, ParseException{
		MatrixWithHeaders expression=new MatrixWithHeaders(gctFile.getAbsolutePath());
		Map<String, ArrayList> groupIndexes=ParseGCTFile.getGroupIndexes(clsFile);
		MatrixWithHeaders diffGenes=computeDiffBtwnGroups(expression, groupIndexes);
		diffGenes.writeGCT(save+".diff.gct");
		diffGenes.writeGMT(save+".gmt", fold);
	}
	
	private MatrixWithHeaders computeDiffBtwnGroups(MatrixWithHeaders expression, Map<String, ArrayList> groupIndexes){
		MatrixWithHeaders diff=null;
		
		for(int i=0; i<groupIndexes.size(); i++){
			for(int j=(i+1); j<groupIndexes.size(); j++){
				String group=groupIndexes.keySet().toArray()[i].toString();
				String group2=groupIndexes.keySet().toArray()[j].toString();
				System.err.println(group+" "+group2);
				String name=group+"_vs_"+group2;
								
				Map<String, double[]> g1=expression.submatrixByColumnIndex(groupIndexes.get(group)).toMap();
				Map<String, double[]> g2=expression.submatrixByColumnIndex(groupIndexes.get(group2)).toMap();
				int n1=groupIndexes.get(group).size();
				int n2=groupIndexes.get(group2).size();
				
				
				Map<String, Double> fold=fold(g1, g2, isLog);
				MatrixWithHeaders foldsMatrix=convertToMatrix(fold, name);
				if(diff==null){diff=foldsMatrix;}
				else{diff.appendColumns(foldsMatrix);}
			
			}
		}
		return diff;
	}
	
	private MatrixWithHeaders convertToMatrix(Map<String, Double> map, String name){
		List<String> columns=new ArrayList();
		columns.add(name);
		MatrixWithHeaders rtrn=new MatrixWithHeaders(new ArrayList(map.keySet()), columns);
		
		for(String gene: map.keySet()){
			rtrn.set(gene, 0, map.get(gene));
		}
		
		return rtrn;
	}
	
	private Map<String, Double> fold(Map<String, double[]> g1, Map<String, double[]> g2, boolean isLog){
		Map<String, Double> rtrn=new TreeMap();
		for(String gene: g1.keySet()){
			if(!gene.equalsIgnoreCase("header")){
				double[] gr1=g1.get(gene);
				double[] gr2=g2.get(gene);
				double fold=Statistics.average(gr1)-Statistics.average(gr2);
				if(!isLog){fold=Statistics.average(gr1)/Statistics.average(gr2);}
				//System.err.println(gene+" "+Statistics.average(gr1)+" "+Statistics.average(gr2)+" "+fold+" "+gr1.length+" "+gr2.length);
				rtrn.put(gene, fold);
			}
		}
		return rtrn;
	}
	
	public static void main(String[] args)throws Exception{
		if(args.length>3){
		File gctFile=new File(args[0]);
		File clsFile=new File(args[1]);
		String save=args[2];
		double fold=new Double(args[3]);
		new DefineGeneSignatures(gctFile, clsFile, save, fold);
		}
		else{System.err.println(usage);}
	}
	static String usage=" args[0]=gct file \n args[1]=cls file \n args[2]=save \n args[3]=fold change";
}
