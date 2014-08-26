package broad.pda.differentialExpression;

import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.log4j.Logger;

import Jama.Matrix;
import broad.core.datastructures.MatrixWithHeaders;
import broad.core.math.EmpiricalDistribution;
import broad.core.math.FDRDistribution;
import broad.core.math.Statistics;

public class DifferentialExpression {
	static Logger logger = Logger.getLogger(DifferentialExpression.class.getName());
	static double[] fudgeFactorDefaults={0.01};
	double[] fudgeFactors;
	double alpha=.05;
	int numberPermutations=1000;
	
	MatrixWithHeaders testStatistics;
	MatrixWithHeaders fdrMatrix;
	MatrixWithHeaders absFdrMatrix;
	MatrixWithHeaders fwerMatrix;
	MatrixWithHeaders pvalueMatrix;
	MatrixWithHeaders[] permutationMatrix;
	boolean paired=false;
	Collection<String> group1;
	Collection<String> group2;
	MatrixWithHeaders data;
	Map<String, FDRDistribution> fdrDists;
	
	//TODO Cache the permutations for the same number of replicates
	
	
	public DifferentialExpression(MatrixWithHeaders data, Collection<String> group1, Collection<String> group2, boolean useFold, boolean paired, int perm, double[] fudgeFactors, MatrixWithHeaders[] permutationMatrix, Map<String, FDRDistribution> fdrDists){
		this.fudgeFactors=fudgeFactors;
		this.numberPermutations=perm;
		this.group1=group1;
		this.group2=group2;
		this.data=data;
		this.paired=paired;
		this.permutationMatrix=permutationMatrix;
		this.pvalueMatrix=null;
		this.fwerMatrix=null;
		this.fdrDists=fdrDists;
		
		//Step 1: Compute test statistics
		logger.debug("Computing test stats for group1 " + group1 +" and group2 " + group2);
		this.testStatistics=DifferentialScoring.computeTestStatistics(data, group1, group2, fudgeFactors, useFold, paired);
				
		//Step 2: Generate permutations
		if(permutationMatrix==null){
			logger.trace("Permutation matrix is null");
			this.permutationMatrix=Permutations.permuteData(data, group1, group2, numberPermutations, fudgeFactors, useFold, paired);
		}
	
		//Step 3: Assess significance of each gene using the FDR`
		if(this.fdrDists==null){
			MatrixWithHeaders[] matrix=assignFDRValues(testStatistics, this.permutationMatrix);
			this.fdrMatrix=matrix[0];
			this.absFdrMatrix=matrix[1];
		}
		else{
			MatrixWithHeaders[] matrix=assignFDRValues(testStatistics, this.fdrDists);
			this.fdrMatrix=matrix[0];
			this.absFdrMatrix=matrix[1];
		}
		
		//System.err.println("Completed computing FDRs for all genes...");
		
	}
	
	public DifferentialExpression(MatrixWithHeaders data, Map<String, Collection<String>> groups, boolean useFold, boolean paired, int perm, double[] fudgeFactors){
		this.fudgeFactors=fudgeFactors;
		this.numberPermutations=perm;
		this.data=data;
		this.paired=paired;
		this.pvalueMatrix=null;
		this.fwerMatrix=null;
		
		//Step 1: Compute test statistics
		long start=System.currentTimeMillis();
		this.testStatistics=DifferentialScoring.computeTestStatistics(data, groups, fudgeFactors, useFold, paired);
		long end=System.currentTimeMillis();
		logger.debug("Test stat took "+(end-start));
		
		
		//Step 2: Generate permutations
		if(permutationMatrix==null){
			this.permutationMatrix=Permutations.permuteData(data, groups, numberPermutations, fudgeFactors, useFold, paired);
		}
	
		//Step 3: Assess significance of each gene using the FDR
		if(this.fdrDists==null){
			MatrixWithHeaders[] matrix=assignFDRValues(testStatistics, this.permutationMatrix);
			this.fdrMatrix=matrix[0];
			this.absFdrMatrix=matrix[1];
		}
		else{
			MatrixWithHeaders[] matrix=assignFDRValues(testStatistics, this.fdrDists);
			this.fdrMatrix=matrix[0];
			this.absFdrMatrix=matrix[1];
		}
		
		//System.err.println("Completed computing FDRs for all genes...");
		
	}
	
	public DifferentialExpression(MatrixWithHeaders data, Collection<String> group1, Collection<String> group2, boolean useFold, boolean paired, int perm, double[] fudgeFactors, MatrixWithHeaders[] permutationMatrix){
		this(data, group1, group2, useFold, paired, perm, fudgeFactors, permutationMatrix, null);
	}
	
	

	private Map<String, Matrix> makeMatrixMap(MatrixWithHeaders data, MatrixWithHeaders[] permutations) {
		Map<String, Matrix> rtrn=new TreeMap<String, Matrix>();
		
		for(String column: data.getColumnNames()){
			Matrix matrix=new Matrix(data.rowDimension(), permutations.length);
			for(int i=0; i<permutations.length; i++){matrix.setColumn(i, permutations[i].getColumn(column));}
			rtrn.put(column, matrix);
		}
		
		return rtrn;
	}

	public DifferentialExpression(MatrixWithHeaders data, Collection<String> group1, Collection<String> group2, boolean useFold, boolean paired, int perm){
		this(data, group1, group2, useFold, paired, perm, fudgeFactorDefaults, null);
	}
	
	private MatrixWithHeaders assignNominalPValues(MatrixWithHeaders data, MatrixWithHeaders[] permutations) {
		MatrixWithHeaders rtrn=new MatrixWithHeaders(data.getRowNames(), data.getColumnNames()); 
		
		//Loop through each test-stat
		for(String column: data.getColumnNames()){
			Matrix matrix=new Matrix(data.rowDimension(), permutations.length);
			//For each permutation
			for(int i=0; i<permutations.length; i++){matrix.setColumn(i, permutations[i].getColumn(column));}
			EmpiricalDistribution dist=new EmpiricalDistribution(matrix, 200);
			
			int counter=0;
			for(String row: data.getRowNames()){
				double val=data.get(row, column);
				double pvalue=dist.getCumulativeProbability(val); 
				if(val>0){pvalue=1-dist.getCumulativeProbability(val);}
				rtrn.set(row, column, pvalue);
				counter++;
			}
		}
			
			return rtrn;
	}
	
	

	public DifferentialExpression(MatrixWithHeaders data, Collection<String> group1, Collection<String> group2, boolean useFold, boolean paired){
		this(data, group1, group2, useFold, paired, 200);
	}
	
	public DifferentialExpression(MatrixWithHeaders data, Collection<String> group1, Collection<String> group2){ 
		this(data, group1, group2, false, false, 200);
	}
	
	public DifferentialExpression(MatrixWithHeaders data, Collection<String> group1, Collection<String> group2, int perm){
		this(data, group1, group2, false, false, perm);
	}
		

	

	
	
	/*private MatrixWithHeaders[] permuteData(MatrixWithHeaders data, Collection<String> group1, Collection<String> group2, int numberPermutations) {
		MatrixWithHeaders[] rtrn;
		CombinationGenerator comb=new CombinationGenerator((group1.size()+group2.size()), group1.size());
		ArrayList<String> samples=new ArrayList<String>();
		samples.addAll(group1);
		samples.addAll(group2);
		
		int totalNumberPossible=comb.getTotal();
		
		if(totalNumberPossible<numberPermutations){
			//do all
			rtrn=new MatrixWithHeaders[totalNumberPossible];
			int i=0;
			while(comb.hasMore()){
				int[] next=comb.getNext();
				rtrn[i]=generatePermutation(data, next, samples);
				i++;
			}
		}
		
		else{
			//do numPermutations randomly
			rtrn=new MatrixWithHeaders[numberPermutations];
			for(int i=0; i<numberPermutations; i++){
				int[] next=comb.getNextRandom();
				rtrn[i]=generatePermutation(data, next, samples);
				i++;
			}
		}
		
		return rtrn;
	}

	private MatrixWithHeaders generatePermutation(MatrixWithHeaders data, int[] permutedIndeces, ArrayList<String> samples) {
		Collection<String> group1=new TreeSet<String>();
		Collection<String> group2=new TreeSet<String>();
		for(int i=0; i<permutedIndeces.length; i++){group1.add(samples.get(i));}
		for(String sample: samples){if(!group1.contains(sample)){group2.add(sample);}}
		
		return DifferentialScoring.computeTestStatistics(data, group1, group2, fudgeFactors);
	}*/

	public Map<String, FDRDistribution> getFDRDistribution(){
		return this.fdrDists;
	}
	
	private MatrixWithHeaders[] assignFDRValues(MatrixWithHeaders data, MatrixWithHeaders[] permutations) {
		MatrixWithHeaders fdr=new MatrixWithHeaders(data.getRowNames(), data.getColumnNames()); 
		MatrixWithHeaders absFDR=new MatrixWithHeaders(data.getRowNames(), data.getColumnNames()); 
		
		this.fdrDists=new TreeMap<String, FDRDistribution>();
		
		for(String column: data.getColumnNames()){
			logger.debug("processing colname " + column);
			Matrix matrix=new Matrix(data.rowDimension(), permutations.length);
			for(int i=0; i<permutations.length; i++){matrix.setColumn(i, permutations[i].getColumn(column));}
			double[] observed=data.getColumn(column);
			FDRDistribution dist=new FDRDistribution(observed, matrix, alpha);
			fdrDists.put(column, dist);
			int counter=0;
			for(String row: data.getRowNames()){
				double val=data.get(row, column);
				double fdrVal=dist.getFDR(val); //TODO: Consider whether to do this or not
				double absFDRVal=dist.getAbsFDR(val);
				absFDR.set(row, column, absFDRVal);
				fdr.set(row, column, fdrVal);
				counter++;
			}
		}
		
		MatrixWithHeaders[] rtrn={fdr, absFDR};
		
		return rtrn;
	}
	
	private MatrixWithHeaders[] assignFDRValues(MatrixWithHeaders data, Map<String, FDRDistribution> fdrDists) {
		MatrixWithHeaders fdr=new MatrixWithHeaders(data.getRowNames(), data.getColumnNames()); 
		MatrixWithHeaders absFDR=new MatrixWithHeaders(data.getRowNames(), data.getColumnNames()); 
		
		for(String column: data.getColumnNames()){
			double[] observed=data.getColumn(column);
			FDRDistribution dist=new FDRDistribution(observed, fdrDists.get(column));
			int counter=0;
			for(String row: data.getRowNames()){
				double val=data.get(row, column);
				double fdrVal=dist.getFDR(val); //TODO: Consider whether to do this or not
				double absFDRVal=dist.getAbsFDR(val);
				absFDR.set(row, column, absFDRVal);
				fdr.set(row, column, fdrVal);
				counter++;
			}
		}
		
		MatrixWithHeaders[] rtrn={fdr, absFDR};
		
		return rtrn;
	}
	
	private MatrixWithHeaders assignAbsFDRValues(MatrixWithHeaders data, MatrixWithHeaders[] permutations) {
		MatrixWithHeaders fdr=new MatrixWithHeaders(data.getRowNames(), data.getColumnNames()); 
		
		for(String column: data.getColumnNames()){
			Matrix matrix=new Matrix(data.rowDimension(), permutations.length);
			for(int i=0; i<permutations.length; i++){matrix.setColumn(i, permutations[i].getColumn(column));}
			double[] observed=data.getColumn(column);
			FDRDistribution dist=new FDRDistribution(observed, matrix, alpha);
			
			int counter=0;
			for(String row: data.getRowNames()){
				double val=data.get(row, column);
				double fdrVal=dist.getAbsFDR(val); //TODO: Consider whether to do this or not
				//System.err.println(val+" "+dist.getFDR(val)+" abs fdr "+dist.getAbsFDR(val));
				fdr.set(row, column, fdrVal);
				counter++;
				//if(counter % 500 ==0){System.err.println(counter+" "+row+" "+fdrVal);}
			}
		}
		
		return fdr;
	}
	
	private MatrixWithHeaders assignFWERValues(MatrixWithHeaders data, MatrixWithHeaders[] permutations) {
		MatrixWithHeaders fwer=new MatrixWithHeaders(data.getRowNames(), data.getColumnNames()); 
		
		//get experiment random distribution
		for(String column: data.getColumnNames()){
			Matrix randomMatrix=new Matrix(data.rowDimension(), permutations.length);
			for(int i=0; i<permutations.length; i++){randomMatrix.setColumn(i, permutations[i].getColumn(column));}
			double[] observed=data.getColumn(column);
			FWERDistribution dist=new FWERDistribution(observed, randomMatrix);
			
			for(String row: data.getRowNames()){
				double fwerVal=dist.getFWER(data.get(row, column));
				//System.err.println("val " +data.get(row, column) +" fwer: " + fwerVal );
				fwer.set(row, column, fwerVal);
			}
		}
		
		return fwer;
	}

	public MatrixWithHeaders getTestStatisticMatrix(){return this.testStatistics;}
	public MatrixWithHeaders getFDRMatrix(){return this.fdrMatrix;}
	public MatrixWithHeaders getAbsFDRMatrix(){return this.absFdrMatrix;}
	public MatrixWithHeaders[] getPermutationMatrix(){return this.permutationMatrix;}

	public MatrixWithHeaders getFoldMatrix() {
		return DifferentialScoring.computeAbsFoldMatrix(data, group1, group2);
	}

	public MatrixWithHeaders getFWERMatrix() {
		if(this.fwerMatrix==null){this.fwerMatrix=assignFWERValues(testStatistics, this.permutationMatrix);}
		return this.fwerMatrix;
	}

	public MatrixWithHeaders getNominalPValueMatrix() {
		if(this.pvalueMatrix==null){this.pvalueMatrix=assignNominalPValues(testStatistics, this.permutationMatrix);}
		return this.pvalueMatrix;
	}

	public MatrixWithHeaders getGenesPassingFDR(double alpha) {
		Collection<String> list=new TreeSet<String>();
		for(String gene: this.fdrMatrix.getRowNames()){
			double[] fdrs=this.fdrMatrix.getRow(gene);
			double fdr=Statistics.min(fdrs);
			if(fdr<alpha){list.add(gene);}
		}
		
		logger.trace("Before "+this.fdrMatrix.rowDimension()+" after "+list.size());
		
		return this.data.submatrixByRowNames(list);
	}

	/*private MatrixWithHeaders computeTestStatistics(MatrixWithHeaders data, Collection<String> group1, Collection<String> group2, double[] fudgeFactors) {
		List<String> columns=new ArrayList<String>();
		for(int i=0; i<fudgeFactors.length; i++){columns.add(""+fudgeFactors[i]);}
		MatrixWithHeaders statisticMatrix=new MatrixWithHeaders(data.getRowNames(), columns);
		
		MatrixWithHeaders subset1=data.submatrixByColumnNames(group1);
		MatrixWithHeaders subset2=data.submatrixByColumnNames(group2);
		
		for(int i=0; i<fudgeFactors.length; i++){
			for(String gene: data.getRowNames()){
				double[] gr1=subset1.getRow(gene);
				double[] gr2=subset2.getRow(gene);
				double s=DifferentialScoring.computeTestStatistic(gr1, gr2, fudgeFactors[i]);
				statisticMatrix.set(gene, i, s);
			}
		}
		
		return statisticMatrix;
	}*/

	/*private double computeTestStatistic(double[] gr1, double[] gr2, double fudgeFactor) {
		Collection<double[]> groups=new ArrayList<double[]>();
		groups.add(gr1);
		groups.add(gr2);
		return computeTestStatistic(groups, fudgeFactor);
	}

	private double computeTestStatistic(Collection<double[]> groups, double fudgeFactor) {
		if(groups.size()>2){
			//do an ANOVA
			//TODO: add fudge factor to the computation
			double f=Statistics.anovaFStat(groups);
			return f;
		}
		else if(groups.size()==2){
			double[] gr1=(double[])groups.toArray()[0];
			double[] gr2=(double[]) groups.toArray()[1];
			if(gr1.length==1 && gr2.length==1){
				//return fold change
				double fold=gr1[0]-gr2[0];
				return fold;
			}
			else if(gr1.length==1 || gr2.length==1){
				//do a z-test
				double z;
				if(gr1.length==1){z=Statistics.zScore(gr1[0], gr2, fudgeFactor);}
				else{z=Statistics.zScore(gr2[0], gr1, fudgeFactor);}
				return z;
			}
			else{
				//do a t-test
				double t=Statistics.tstat(gr1, gr2, fudgeFactor);
				return t;
			}
		}
		else{throw new IllegalArgumentException("Must be at least 2 groups to compute differences");}
		
	}*/
	
}
