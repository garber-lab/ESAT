package broad.pda.geneexpression;

import java.util.ArrayList;
import java.util.Map;

import broad.core.datastructures.MatrixWithHeaders;

public class ScoreDiffGenes {

	public static MatrixWithHeaders scoreListByTStatistic(MatrixWithHeaders expression, Map<String, ArrayList> groupIndexes, String controlClass, double fudgeFactor){
		MatrixWithHeaders rankedList=null;
		
		boolean useControlClass=false;
		if(controlClass!=null && groupIndexes.containsKey(controlClass)){
			useControlClass=true;
		}
		
		for(String group1: groupIndexes.keySet()){
			if(!useControlClass || group1.equalsIgnoreCase(controlClass)){
			for(String group2: groupIndexes.keySet()){
				String name=group1+"_vs_"+group2;
				System.err.println(name);
				MatrixWithHeaders vals=null;
				if(!group1.equalsIgnoreCase(group2)){
					vals=expression.tScoresByRow(groupIndexes.get(group1), groupIndexes.get(group2), name, fudgeFactor);
				}
					else{
						name=group1+"_vs_REST";
						System.err.println(name);
						vals=expression.tScoresByRow(groupIndexes.get(group1), getREST(group1, groupIndexes), name, fudgeFactor);
					}
				if(rankedList==null){rankedList=vals;}
				else{rankedList.appendColumns(vals);}
			}
			}
			}
		
		
		return rankedList;
	}
	
	

	
	
	public static ArrayList getREST(String group1, Map<String, ArrayList> groupIndexes){
		ArrayList rtrn=new ArrayList();
		
		for(String name: groupIndexes.keySet()){
			if(!name.equalsIgnoreCase(group1)){
				rtrn.addAll(groupIndexes.get(name));
			}
		}
		
		return rtrn;
	}
	
}
