package broad.pda.geneexpression.agilent.analysis;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.core.datastructures.MatrixWithHeaders;
import broad.core.error.ParseException;
import broad.core.math.Statistics;
import broad.pda.differentialExpression.DifferentialExpression;
import broad.pda.differentialExpression.RunDifferentialExpression;
import broad.pda.geneexpression.ExpressionExperimentInfo;
import broad.pda.geneexpression.agilent.AgilentUtils;
import broad.pda.geneexpression.clustering.HierarchicalClustering;

public class ProcessDataset {
	
	private double threshold=2.0;
	private double alpha=0.05;
	private int minSamplesPerGroup=2; //TODO Need to be able to set this value
	private int numPerm=500;
	double[] fudgeFactors={0, 0.01, 0.1, 1};
	boolean clusterRows=false;
	boolean clusterColumns=true;
	
	
	public ProcessDataset(MatrixWithHeaders data, Map<String, Collection<String>> groups, Map<String, String> sampleMapping, String save) throws IOException, ParseException{
		//Filter the number of samples per group
		groups=filterGroups(groups, minSamplesPerGroup);
		
		//Filter by fold differences from controls
		MatrixWithHeaders filteredData=filter(data, groups, threshold);
		
		//TODO Filter variable genes
		//Step 2: Filter all genes with highly variable negatives
		
		//Step 3: Compute differences to median and filter unchanging genes
		
		//Compute differentially expressed genes and filter by this
		RunDifferentialExpression diffExpression=new RunDifferentialExpression(filteredData, groups,100,0.05);
		filteredData=diffExpression.getGenesPassingFDR(alpha);
		
		//TODO speed up ANOVA
		//DifferentialExpression anova=new DifferentialExpression(filteredData, groups, true, false, numPerm, fudgeFactors);
		//filteredData=anova.getGenesPassingFDR(alpha);
		
		/* TODO: uncomment, HierarchicalClustering(filteredData, groups, null)
		 //was not available
		 */
		
		//Step 4: Cluster by keeping groups together
		HierarchicalClustering cluster=new HierarchicalClustering(filteredData, groups, null);
		cluster.cluster(clusterColumns, clusterRows);
		
		filteredData=cluster.getMatrix();
		
		filteredData.writeGCTWithHeaders(save, sampleMapping);
		
		//Compute diff expression and filter genes and cluster by these genes
		
	}

	private Map<String, Collection<String>> filterGroups(Map<String, Collection<String>> groups, int minSamplesPerGroup) {
		Map<String, Collection<String>> rtrn=new TreeMap<String, Collection<String>>();
		for(String group: groups.keySet()){
			Collection<String> vals=groups.get(group);
			if(vals.size()>=minSamplesPerGroup){rtrn.put(group, vals);}
		}
		return rtrn;
	}

	private MatrixWithHeaders filter(MatrixWithHeaders data, Map<String, Collection<String>> groups, double threshold) {
		Collection<String> filteredGenes=new TreeSet<String>();
		Collection<String> controls=groups.get("Control");
		MatrixWithHeaders controlVals=data.submatrixByColumnNames(controls);
		
		for(String group: groups.keySet()){
			if(!group.equalsIgnoreCase("Control")){
				Collection<String> group1=groups.get(group);
				MatrixWithHeaders groupVals=data.submatrixByColumnNames(group1);
				
				//Loop through each gene
				for(String gene: data.getRowNames()){
					if(!filteredGenes.contains(gene)){
						double[] c1=controlVals.getRow(gene);
						double[] g1=groupVals.getRow(gene);
						double fold=Statistics.absFold(c1, g1, true);
						if(g1.length==2){
							double[] v1={g1[0]};
							double[] v2={g1[1]};
							double f1=Statistics.absFold(c1, v1, true);
							double f2=Statistics.absFold(c1, v2, true);
							fold=Math.min(f1, f2);
						}
						if(fold>threshold){filteredGenes.add(gene);}
					}
				}
				
			}
		}
		
		System.err.println("Started with "+data.getRowNames().size()+" Ended with "+filteredGenes.size());
		return data.submatrixByRowNames(filteredGenes);
	}
	
	public static void main(String[] args) throws IOException, ParseException{
		if(args.length>2){
			MatrixWithHeaders data=new MatrixWithHeaders(args[0]);
			Map<String, Collection<String>> groups=AgilentUtils.parseExperimentInfoFileToGroups(new File(args[1]));
			Map<String, String> mapping=AgilentUtils.parseExperimentInfoFileToName(new File(args[1]));
			String save=args[2];
			new ProcessDataset(data, groups, mapping, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=gct \n args[1]=groups \n args[2]=save";
}
