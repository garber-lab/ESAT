package broad.core.math;

import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.Set;

import broad.core.datastructures.MatrixWithHeaders;


public class KSTest {

	
	
	public static double[] KSScores(MatrixWithHeaders rankedList, Collection<String> geneSet, int columnNumber) {

		/**Declaring variables**/
		int rankAtMaxEs=0;
		double scoresAtMaxEs=0;
		//HitIndices hitIndices=new HitIndices(geneSet.size()); //size of geneSet
		double ess_maxdev=0.0;
		double runningScores=0.0;
		int hitCnt=0;
		double ess_pos_list=0.0;
		int rankAtMaxEs_pos_list=0;
		double scoresAtMaxEs_pos_list=0;
		double ess_pos_list_maxdev=0.0;
		int rankAtMaxEs_pos_list_maxdev=0;
		double scoresAtMaxEs_pos_list_maxdev=0;
		double ess_neg_list=0.0;
		int rankAtMaxEs_neg_list=0;
		double scoresAtMaxEs_neg_list=0;
		double ess_neg_list_maxdev=0.0;
		int rankAtMaxEs_neg_list_maxdev=0;
		double scoresAtMaxEs_neg_list_maxdev=0;
		double totalWeight=computeTotalWeight(geneSet, rankedList, columnNumber);
		double missScore=computeMissScore(geneSet, rankedList);
		List<String> geneNames=rankedList.getRowNames();
		/************************/
		
		
		//iterate through each gene
		for (int r = 0; r < rankedList.getNumberRows(); r++) {
		
		
		
		//get gene name
		String rowName =geneNames.get(r); 
		//get score for gene in original ranked list
		double corr = rankedList.get(r, columnNumber);
		//check if its positive
		boolean posList = (corr>=0);
				
		
		/********Need to return and implement scoring of hits********/	
		//if gene (row name) is in the gene set
		double sr=missScore;
		if(geneSet.contains(rowName)){sr=hitScore(totalWeight, rankedList, columnNumber, rowName);}
		runningScores += sr;
		/***********************************************************/
		
		//if ess_maxdev is less than the running score (which means we just updated in this iteration)
		if (Math.abs(ess_maxdev) < Math.abs(runningScores)) { // @note abs here
			//set ess_maxdev to runningScore
			ess_maxdev =  runningScores; // @note no abs here!
			//store the rank at max ES to the current gene
			rankAtMaxEs = r;
			//score at max ES is the score of the ranked list for the gene
			scoresAtMaxEs = rankedList.get(r, columnNumber);
		}
	
		
		//if positive
		if (posList) {
			//if ess_pos less than running score
			if (ess_pos_list < runningScores) { // @note NO abs
				//then update ess_pos to running score
				ess_pos_list = runningScores;
				//set rank at max to current gene
				rankAtMaxEs_pos_list = r;
				//set score at max to current score
				scoresAtMaxEs_pos_list = rankedList.get(r, columnNumber);
			}
		
			//if abs value of ess_pos_maxdev less than absolute running score
			if (Math.abs(ess_pos_list_maxdev) < Math.abs(runningScores)) { // @note YES abs
				//then update the ess_pos_maxdev to running score
				ess_pos_list_maxdev =  runningScores;
				//set rank to current gene
				rankAtMaxEs_pos_list_maxdev = r;
				//set score to current score
				scoresAtMaxEs_pos_list_maxdev = rankedList.get(r, columnNumber);
			}
		} 
		
		//else if negative
		else {
			//if ess_neg greater than running score
			if (ess_neg_list > runningScores) { // Note NO abs
				//then set ess_neg to running score
				ess_neg_list =  runningScores;
				//set rank at maxES_neg to current gene
				rankAtMaxEs_neg_list = r;
				//set score to current gene
				scoresAtMaxEs_neg_list = rankedList.get(r, columnNumber);
			}
		
			//if abs(ess_neg_maxdev)<runningScore
			if (Math.abs(ess_neg_list_maxdev) < Math.abs(runningScores)) { // @note abs
				//then update
				ess_neg_list_maxdev = runningScores;
				rankAtMaxEs_neg_list_maxdev = r;
				scoresAtMaxEs_neg_list_maxdev = rankedList.get(r, columnNumber);;
			}
		}
		
		
		}
		
		
		//System.err.println(ess_maxdev+" "+ess_pos_list_maxdev+" "+ess_neg_list_maxdev);
		// ES, rank at ES, value at ES
		double[] rtrn={ess_maxdev, rankAtMaxEs, scoresAtMaxEs};
		return rtrn;
	
	}
	
	
	
	private static double computeTotalWeight(Collection<String> geneSet, MatrixWithHeaders rankedList, int columnIndex){
		double totalWeight=0;
		for (String gene: geneSet) {
			double score = rankedList.get(gene, columnIndex);
            totalWeight += _abs(score);
        }
		return totalWeight;
	}
	
	private static double computeMissScore(Collection<String> geneSet, MatrixWithHeaders rankedList){
		double nhExpected = geneSet.size();
    	

        final float nTotal = rankedList.getNumberRows();
        double miss_score = 1.0f / (nTotal - nhExpected);
        return -miss_score;
	}

		private static double hitScore(double totalWeight, MatrixWithHeaders rankedList, int columnIndex, String geneName) {
	        	
	        	   	double score = rankedList.get(geneName, columnIndex);
		            score = _abs(score);
		            return score / totalWeight;
	            
	    }

		 // Needed as cdna give some nans for the class metric
	    private static double _abs(double score) {
	        if (Double.isNaN(score) || Double.isInfinite(score)) {
	            return 0.01;
	        } else {
	            return Math.abs(score);
	        }
	    }

	    
}
