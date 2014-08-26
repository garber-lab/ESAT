package broad.pda.differentialExpression;

import java.io.*;
import java.util.*;

import org.apache.log4j.Logger;

import broad.core.datastructures.MatrixWithHeaders;
import broad.core.error.ParseException;
import broad.core.math.FDRDistribution;
import broad.core.math.Statistics;
import broad.core.util.GMTParser;
import broad.core.util.ParseGCTFile;
import broad.pda.geneexpression.ExpressionExperimentInfo;
import broad.pda.geneexpression.agilent.AgilentUtils;
import broad.projection.gsa.GeneSetEnrichment;
 
//TODO: Make sure absolute FDR is correct

public class RunDifferentialExpression {
	
	double alpha=.05;
	private double cutoff=2.0;
	private boolean preprocessByGroups=true;
	private boolean preprocess=false;
	int minNumPerGroup=2;
	private int numPerm;
	boolean useFold=false;
	double[] fudgeFactors={0, 0.01, 0.1, 1};
	boolean useAbsoluteFDR=true;
	private MatrixWithHeaders fdrMatrix;
	private MatrixWithHeaders data;
	static Logger logger = Logger.getLogger(RunDifferentialExpression.class.getName());
	
	
	public MatrixWithHeaders getResults() { return fdrMatrix;}
	public void writeResults(String outFile) throws IOException {
		fdrMatrix.write(outFile);
	}
	
	//TODO Global permutations
	public RunDifferentialExpression(File gctFile, File experimentInfoFile, File chipFile, String save, boolean preprocess, double alpha, int perm) throws IOException, ParseException{
		this.numPerm=perm;
		this.alpha=alpha;
		data=new MatrixWithHeaders(gctFile.getAbsolutePath());
		Map<String, ExpressionExperimentInfo> experimentInfo=AgilentUtils.parseExperimentInfoFile(experimentInfoFile);
		Map<String, String> annotations=ParseGCTFile.parseChipFile(chipFile);
		data.setPIDToName(annotations);
		
		//Filter data matrix by annotation file
		data=data.submatrixByRowNames(annotations.keySet());
		
		//Filter by columns in info file
		data=data.submatrixByColumnNames(experimentInfo.keySet());
		
		//print the submatrix
		//data.writeGCT(gctFile.getName()+"."+experimentInfoFile.getName()+".gct");
		
		System.err.println("Starting...");
		
		Collection<String> controls=getControls(data, experimentInfo);
		
		Map<String, Collection<String>> groups=getGroups(data, experimentInfo, minNumPerGroup);
		
		//Run individual comparisons
		MatrixWithHeaders fdrMatrix=runComparisons(data, groups, controls, preprocess);
		
		//write output
		fdrMatrix.setPIDToName(annotations);
		fdrMatrix.writeGCT(save);
		
		
		//MatrixWithHeaders geneSetsAll=runComparisonsAgainstAll(data, experimentInfo, annotations, save, preprocess);
		
		//Run all groups against the gene set
		//runGeneSetComparisons(data, experimentInfo, geneSets, save);
	}
	
	/**
	 * 
	 * @param data Data set in which to find differentially expressed genes
	 * @param group1 The list of samples (columns in the data matrix) of group 1
	 * @param group2 The list of samples (columns in the data matrix) of group 2
 	 * @param numPerm Number of permutations to perform. Specify 0 to perform the maximal number of permutations
	 * @param alpha Significance level
	 * @throws IOException
	 */
	public RunDifferentialExpression(MatrixWithHeaders data, Collection<String> group1, Collection<String> group2, int numPerm, double alpha)  {
		this.numPerm=numPerm;//TODO These should be params
		this.alpha=alpha; //TODO These should be params
		this.preprocess=false; //TODO Consider deleting or setting
		this.data=data;
			
		System.err.println("Computing diff expression... group1 " + group1 + " group2 " + group2 );
		
		Collection<String> controls=group1;
		
		Map<String, Collection<String>> mapOfGroup2 = new HashMap<String, Collection<String>>();
		mapOfGroup2.put("group2", group2);		//Run individual comparisons
		fdrMatrix=runComparisons(data, mapOfGroup2, controls, preprocess);
	}
	
	public RunDifferentialExpression(MatrixWithHeaders data, Map<String, Collection<String>> groups, int numPerm, double alpha) throws IOException {
		this.numPerm=numPerm;//TODO These should be params
		this.alpha=alpha; //TODO These should be params
		this.preprocess=false; //TODO Consider deleting or setting
		this.data=data;
			
		System.err.println("Computing diff expression...");
		
		Collection<String> controls=groups.get("Control");
		
		//Run individual comparisons
		fdrMatrix=runComparisons(data, groups, controls, preprocess);
		
		
		fdrMatrix.writeGCT("matrix.fdr");
		//MatrixWithHeaders geneSetsAll=runComparisonsAgainstAll(data, experimentInfo, annotations, save, preprocess);
		
		//Run all groups against the gene set
		//runGeneSetComparisons(data, experimentInfo, geneSets, save);
	}
	
	public MatrixWithHeaders getGenesPassingFDR(double alpha){
		Collection<String> rtrn=new TreeSet<String>();
		
		for(String gene: fdrMatrix.getRowNames()){
			for(String column: fdrMatrix.getColumnNames()){
				double fdr=fdrMatrix.get(gene, column);
				if(Math.abs(fdr)<alpha){rtrn.add(gene);}
			}
		}
		
		System.err.println("After diff: started with "+this.fdrMatrix.getRowNames().size()+" ended with "+rtrn.size());
		
		return data.submatrixByRowNames(rtrn);
	}
	
	public MatrixWithHeaders getFDR(){return this.fdrMatrix;}

	private Collection<String> getControls(MatrixWithHeaders data, Map<String, ExpressionExperimentInfo> experimentInfo) {
		Collection<String> rtrn=new TreeSet<String>();
		
		for(String barcode: data.getColumnNames()){
			if(experimentInfo.containsKey(barcode)){
				ExpressionExperimentInfo info=experimentInfo.get(barcode);
				String group=info.getSampleType();
				if(group.equalsIgnoreCase("Control")){rtrn.add(barcode);}
			}
		}
		
		return rtrn;
	}

	private Map<String, Collection<String>> getGroups(MatrixWithHeaders data, Map<String, ExpressionExperimentInfo> experimentInfo, int minSampleNumbers) {
		Map<String, Collection<String>> samplesByName=new TreeMap<String, Collection<String>>();
		
		Map<String, Collection<String>> diffGenes=new TreeMap<String, Collection<String>>();
		Collection<String> controls=new ArrayList();
		
		for(String barcode: data.getColumnNames()){
			if(experimentInfo.containsKey(barcode)){
				ExpressionExperimentInfo info=experimentInfo.get(barcode);
				String group=info.getSampleType();
				String name=info.getExperimentName();
				Collection<String> nameSamples=new ArrayList();
				if(samplesByName.containsKey(name)){nameSamples=samplesByName.get(name);}
				nameSamples.add(barcode);
				if(!group.equalsIgnoreCase("Control")){samplesByName.put(name, nameSamples);}
				else{controls.add(barcode);}
			}
		}
		
		Collection<String> remove=new TreeSet();
		for(String name: samplesByName.keySet()){
			int numSamples=samplesByName.get(name).size();	
			if(numSamples<minSampleNumbers){
				remove.add(name);
				System.err.println("Skipping "+name+" which only had "+numSamples+" where a min of "+minSampleNumbers+" was required");
			}
		}
		
		for(String name: remove){samplesByName.remove(name);}
				
		samplesByName.put("Control", controls);
		//System.err.println(samplesByName.keySet());
		return samplesByName;
	}

	private void runGeneSetComparisons(MatrixWithHeaders data, Map<String, ExpressionExperimentInfo> experimentInfo, Map<String, Collection<String>> geneSets, String save) throws IOException {
		Collection<String> negatives=new TreeSet<String>();
		Map<String, Collection<String>> samplesByName=new TreeMap<String, Collection<String>>();
				
		for(String barcode: data.getColumnNames()){
			if(experimentInfo.containsKey(barcode)){
				ExpressionExperimentInfo info=experimentInfo.get(barcode);
				String group=info.getSampleType();
				String name=info.getExperimentName();
				Collection<String> nameSamples=new ArrayList();
				if(samplesByName.containsKey(name)){nameSamples=samplesByName.get(name);}
				nameSamples.add(barcode);
				if(group.equalsIgnoreCase("Control")){negatives.add(barcode);}
				else{samplesByName.put(name, nameSamples);}
			}
		}
		
		Collection<String> genes=new TreeSet<String>();
		
		for(String name: samplesByName.keySet()){
			System.err.println(name);
			GeneSetEnrichment diff=new GeneSetEnrichment(data, negatives, samplesByName.get(name), geneSets);
			MatrixWithHeaders ksFDR=diff.getKSFDR();
			ksFDR.writeGCT(save+"."+name+".ksFDR.gct");
			MatrixWithHeaders maxMeanFDR=diff.getMaxMeanFDR();
			maxMeanFDR.writeGCT(save+"."+name+".maxMeanFDR.gct");
			diff.getNormalizedKSEnrichments().writeGCT(save+"."+name+".normalizedKS.gct");
			diff.getNormalizedMaxMeanEnrichments().writeGCT(save+"."+name+".normalizedMaxMean.gct");
		}
	}

	//TODO Need to get better permutations if we are going to reuse them
	private MatrixWithHeaders runComparisons(MatrixWithHeaders data, Map<String, Collection<String>> groups, Collection<String> controls, boolean preprocess)  {
		MatrixWithHeaders rtrn=initializeMatrix(data.getRowNames(), groups.keySet());
							
		Map<Integer, MatrixWithHeaders[]> permsByGroupSize=new TreeMap<Integer, MatrixWithHeaders[]>();
		Map<Integer, MatrixWithHeaders[]> permsByGroupSizeAll=new TreeMap<Integer, MatrixWithHeaders[]>();
		Map<Integer, Map<String, FDRDistribution>> fdrDistsByGroup=new TreeMap<Integer, Map<String, FDRDistribution>>();
		Map<Integer, Map<String, FDRDistribution>> fdrDistsByGroupAll=new TreeMap<Integer, Map<String, FDRDistribution>>();
		
		for(String name: groups.keySet()){
			if(!name.equalsIgnoreCase("Control")){
				Collection<String> samples=groups.get(name);
				logger.debug("Computing for "+name + " with samples: " + samples);
				MatrixWithHeaders dataMatrix=data.submatrixByColumnNames(samples);
				dataMatrix.appendColumns(data.submatrixByColumnNames(controls));
					
				//Step 0: Preprocess the results
				MatrixWithHeaders preprocessedDataMatrix=preprocess(data, controls, samples, preprocess);
				
				//Step 1: Compute differential expression
				DifferentialExpression diff=new DifferentialExpression(preprocessedDataMatrix, controls, samples, this.useFold, false, this.numPerm, this.fudgeFactors, null,null);
				
				Collection<String> allControls=getAllNonSelf(groups, name, controls);
				DifferentialExpression diffAll=new DifferentialExpression(preprocessedDataMatrix, allControls, samples, this.useFold, false, this.numPerm, this.fudgeFactors, null, null);
				permsByGroupSizeAll.put(samples.size(), diffAll.getPermutationMatrix());
				fdrDistsByGroupAll.put(samples.size(), diffAll.getFDRDistribution());
				
				//Step 3: Assign FDR to all genes, defined as the min of each column
				assignFDRs(rtrn, diff, diffAll, name);
			}
		}
		
		return rtrn;
	}
	
	
	private MatrixWithHeaders initializeMatrix(List<String> rowNames, Set<String> keySet) {
		MatrixWithHeaders rtrn=new MatrixWithHeaders(rowNames, new ArrayList(keySet));
	
		for(String row: rtrn.getRowNames()){
			for(String column: rtrn.getColumnNames()){
				rtrn.set(row, column, 1.0);
			}
		}
		
		return rtrn;
	}

	private void assignFDRs(MatrixWithHeaders rtrn,	DifferentialExpression diff, DifferentialExpression diffAll, String name) {
		MatrixWithHeaders fdr=diff.getFDRMatrix();
		MatrixWithHeaders test=diff.getTestStatisticMatrix();
		
		for(String gene: fdr.getRowNames()){
			double fdrVal=Statistics.min(fdr.getRow(gene));
			double testVal=test.getRow(gene)[0];
			double fdrAll=Statistics.min(diffAll.getFDRMatrix().getRow(gene));				
			double minFDR=Math.min(fdrVal, fdrAll);
			if(this.useAbsoluteFDR){
				minFDR=Math.min(minFDR, Statistics.min(diffAll.getAbsFDRMatrix().getRow(gene)));
				minFDR=Math.min(minFDR, Statistics.min(diff.getAbsFDRMatrix().getRow(gene)));
				//writer.write(gene+"\t"+testVal+"\t"+fdr.get(gene, 0)+"\t"+diff.getAbsFDRMatrix().get(gene, 0)+"\n");
			}
			
			if(testVal<0){minFDR=-minFDR;}
			
			rtrn.set(gene, name, minFDR);
		}
				
	}

	private Collection<String> getAllNonSelf(Map<String, Collection<String>> groups, String name,Collection<String> controls) {
		Collection<String> rtrn=new TreeSet();
		
		for(String groupName: groups.keySet()){
			if(!groupName.equalsIgnoreCase(name)){rtrn.addAll(groups.get(groupName));}
		}
		
		rtrn.addAll(controls);
		
		return rtrn;
	}

	/*private MatrixWithHeaders runComparisons(MatrixWithHeaders data, Map<String, ExpressionExperimentInfo> experimentInfo, Map<String, String> annotations, String save, boolean preprocess) throws IOException {
		Collection<String> negatives=new TreeSet<String>();
		Map<String, Collection<String>> samplesByName=new TreeMap<String, Collection<String>>();
		
		Map<String, Collection<String>> diffGenes=new TreeMap<String, Collection<String>>();
		
		for(String barcode: data.getColumnNames()){
			if(experimentInfo.containsKey(barcode)){
				ExpressionExperimentInfo info=experimentInfo.get(barcode);
				String group=info.getSampleType();
				String name=info.getExperimentName();
				Collection<String> nameSamples=new ArrayList();
				if(samplesByName.containsKey(name)){nameSamples=samplesByName.get(name);}
				nameSamples.add(barcode);
				if(group.equalsIgnoreCase("Control")){negatives.add(barcode);}
				else{samplesByName.put(name, nameSamples);}
			}
		}
		
		Collection<String> groupsToUse=new TreeSet();
		for(String name: samplesByName.keySet()){
			int numSamples=samplesByName.get(name).size();	
			if(numSamples>=this.minNumPerGroup){groupsToUse.add(name);}
			else{System.err.println("Skipping "+name+" which only had "+numSamples+" where a min of "+this.minNumPerGroup+" was required");}
		}
		
		MatrixWithHeaders rtrn=new MatrixWithHeaders(data.getRowNames(), new ArrayList(groupsToUse));
		MatrixWithHeaders fwerRtrn=new MatrixWithHeaders(data.getRowNames(), new ArrayList(groupsToUse));
		MatrixWithHeaders testStat=new MatrixWithHeaders(data.getRowNames(), new ArrayList(groupsToUse));
		
		//Collection<String> genes=new TreeSet<String>();
		
		//Map<String, Collection<String>> geneSets=new TreeMap();
		
		for(String name: groupsToUse){
			MatrixWithHeaders dataMatrix=data.submatrixByColumnNames(samplesByName.get(name));
			dataMatrix.appendColumns(data.submatrixByColumnNames(negatives));
				
			//Step 0: Preprocess the results
			MatrixWithHeaders preprocessedDataMatrix=preprocess(data, negatives, samplesByName.get(name), preprocess);
				
			DifferentialExpression diff=new DifferentialExpression(preprocessedDataMatrix, negatives, samplesByName.get(name), this.useFold, false, this.numPerm);
			MatrixWithHeaders fdr=diff.getFDRMatrix();
			MatrixWithHeaders fwer=diff.getFWERMatrix();
			MatrixWithHeaders pvalues=diff.getNominalPValueMatrix();
				
			MatrixWithHeaders test=diff.getTestStatisticMatrix();
						
			//Setting significant genes
			setSignificantGenes(rtrn, getSignificantGenes(fdr), name, test);
			setSignificantGenes(fwerRtrn, getSignificantGenes(fwer), name, test);
			set(testStat, test, name);
				
			//Write FDRs
			fdr.setPIDToName(annotations);
			fdr.appendColumns(fwer);
			fdr.appendColumns(test);
			fdr.writeGCT(save+"."+name+".controls.fdr");
			pvalues.writeGCT(save+"."+name+".controls.pvalues");
		}
		
		rtrn.setPIDToName(annotations);
		fwerRtrn.setPIDToName(annotations);
		testStat.setPIDToName(annotations);
		//GMTParser.writeGMT(save+".gmt", geneSets);
		
		rtrn.writeGCT(save+".controls.fdr");
		testStat.writeGCT(save+".testStat.gct");
		//fwerRtrn.writeGCT(save+".fwer.gct");
		
		return rtrn;
	}*/
	
	private void set(MatrixWithHeaders testStat, MatrixWithHeaders test, String name) {
		for(String geneName: test.getRowNames()){
			double val=test.get(geneName, 0);
			testStat.set(geneName, name, val);
		}
	}

	private MatrixWithHeaders runComparisonsAgainstAll(MatrixWithHeaders data, Map<String, ExpressionExperimentInfo> experimentInfo, Map<String, String> annotations, String save, boolean preprocess) throws IOException {
		Collection<String> negatives=new TreeSet<String>();
		Map<String, Collection<String>> samplesByName=new TreeMap<String, Collection<String>>();
		
		Map<String, Collection<String>> diffGenes=new TreeMap<String, Collection<String>>();
		
		for(String barcode: data.getColumnNames()){
			if(experimentInfo.containsKey(barcode)){
				ExpressionExperimentInfo info=experimentInfo.get(barcode);
				String group=info.getSampleType();
				String name=info.getExperimentName();
				Collection<String> nameSamples=new ArrayList();
				if(samplesByName.containsKey(name)){nameSamples=samplesByName.get(name);}
				nameSamples.add(barcode);
				if(group.equalsIgnoreCase("Control")){negatives.add(barcode);}
				else{samplesByName.put(name, nameSamples);}
			}
		}
		
		Collection<String> groupsToUse=new TreeSet();
		for(String name: samplesByName.keySet()){
			int numSamples=samplesByName.get(name).size();	
			if(numSamples>=this.minNumPerGroup){groupsToUse.add(name);}
			else{System.err.println("Skipping "+name+" which only had "+numSamples+" where a min of "+this.minNumPerGroup+" was required");}
		}
		
		MatrixWithHeaders rtrn=new MatrixWithHeaders(data.getRowNames(), new ArrayList(groupsToUse));
		MatrixWithHeaders fwerRtrn=new MatrixWithHeaders(data.getRowNames(), new ArrayList(groupsToUse));
		
		for(String name: groupsToUse){
			
			//int numSamples=samplesByName.get(name).size();
			
			//if(numSamples>=this.minNumPerGroup){
			//System.err.println(name);
			//MatrixWithHeaders dataMatrix=data.submatrixByColumnNames(samplesByName.get(name));
			//dataMatrix.appendColumns(data.submatrixByColumnNames(negatives));
			
			Collection<String> allNegatives=new TreeSet<String>();
			Collection<String> samples=new TreeSet<String>();
			
			for(String barcode: data.getColumnNames()){
				ExpressionExperimentInfo experiment=experimentInfo.get(barcode);
				if(experiment.passedQC()){
					if(!experiment.getExperimentName().equalsIgnoreCase(name)){allNegatives.add(barcode);}
					else{samples.add(barcode);}
				}
			}
			
			//Step 0: Preprocess the results
			MatrixWithHeaders preprocessedDataMatrix=preprocess(data, allNegatives, samples, preprocess);
			
			
			System.err.println("Negatives: "+allNegatives.size()+" Samples: "+samples.size());
			
			DifferentialExpression diff=new DifferentialExpression(preprocessedDataMatrix, allNegatives, samples, this.numPerm);
			MatrixWithHeaders fdr=diff.getFDRMatrix();
			MatrixWithHeaders fwer=diff.getFWERMatrix();
			
			MatrixWithHeaders test=diff.getTestStatisticMatrix();
			
			//write test stats
			//test.setPIDToName(annotations);
			//test.writeGCT(save+"."+name+".testStat.gct");
			
			
			
			//write permutation test stats
			/*FileWriter writer=new FileWriter(save+"."+name+".perms.gct");
			MatrixWithHeaders[] perms=diff.getPermutationMatrix();
			for(String row: perms[0].getRowNames()){
				writer.write(row+"\t"+row);
				for(int i=0; i<perms.length; i++){
					writer.write("\t"+perms[i].getRow(row)[0]);
				}
				writer.write("\n");
			}
			writer.close();*/
			
			//write FDR test stats
			/*writer=new FileWriter(save+"."+name+".FDRTest.gct");
			for(String row: fdr.getRowNames()){
				writer.write(row+"\t"+fdr.getRow(row)[0]+"\t"+test.getRow(row)[0]+"\n");
			}	
			writer.close();*/
			
			//Collection<String> significantGenes=getSignificantGenes(fdr);
			
			//Setting significant genes
			setSignificantGenes(rtrn, getSignificantGenes(fdr), name, test);
			setSignificantGenes(fwerRtrn, getSignificantGenes(fwer), name, test);
			
			//for(String row: test.getRowNames()){
			//	testStat.set(row, name, test.get(row, name));
			//}
						
			//MatrixWithHeaders filteredData=data.submatrixByRowNames(significantGenes);
			//MatrixWithHeaders testStat=diff.getTestStatisticMatrix().submatrixByRowNames(significantGenes);
			//testStat.setPIDToName(annotations);
			//filteredData.setPIDToName(annotations);
			
			//Map<String, Collection<String>> geneSet=getGeneSets(testStat, name);
			//geneSets.putAll(geneSet);
			
			//filteredData.writeGCT(save+"."+name+".significantGenes.gct");
			//testStat.writeGCT(save+"."+name+".testStatistic.gct");
			//genes.addAll(filteredData.getRowNames());
			//Collection<String> diffSet=new TreeSet<String>();
			//diffSet.addAll(filteredData.getRowNames());
			//diffGenes.put(name, diffSet);
			
			//Write FDRs
			fdr.setPIDToName(annotations);
			fdr.appendColumns(fwer);
			fdr.appendColumns(test);
			fdr.writeGCT(save+"."+name+".all.fdr");
		//}
		}
		
		rtrn.setPIDToName(annotations);
		fwerRtrn.setPIDToName(annotations);
		//testStat.setPIDToName(annotations);
		//GMTParser.writeGMT(save+".gmt", geneSets);
		
		rtrn.writeGCT(save+".all.fdr");
		//testStat.writeGCT(save+".testStat.gct");
		//fwerRtrn.writeGCT(save+".fwer.gct");
		
		return rtrn;
	}
	
	private void setSignificantGenes(MatrixWithHeaders rtrn, Collection<String> significantGenes, String name, MatrixWithHeaders test) {
		
		for(String geneName: significantGenes){
			double val=test.get(geneName, 0);
			double flag=0;
			if(val>0){flag=1;}
			else{flag=-1;}
			rtrn.set(geneName, name,flag);
		}
		
	}

	/*private Map<String, Collection<String>> runComparisons(MatrixWithHeaders data,	Map<String, ExpressionExperimentInfo> experimentInfo, Map<String, String> annotations, String save, boolean preprocess) throws IOException {
		Collection<String> negatives=new TreeSet<String>();
		Map<String, Collection<String>> samplesByName=new TreeMap<String, Collection<String>>();
		
		Map<String, Collection<String>> diffGenes=new TreeMap<String, Collection<String>>();
		
		for(String barcode: data.getColumnNames()){
			if(experimentInfo.containsKey(barcode)){
				ExpressionExperimentInfo info=experimentInfo.get(barcode);
				String group=info.getSampleType();
				String name=info.getExperimentName();
				Collection<String> nameSamples=new ArrayList();
				if(samplesByName.containsKey(name)){nameSamples=samplesByName.get(name);}
				nameSamples.add(barcode);
				if(group.equalsIgnoreCase("Control")){negatives.add(barcode);}
				else{samplesByName.put(name, nameSamples);}
			}
		}
		
		Collection<String> genes=new TreeSet<String>();
		
		Map<String, Collection<String>> geneSets=new TreeMap();
		
		for(String name: samplesByName.keySet()){
			System.err.println(name);
			MatrixWithHeaders dataMatrix=data.submatrixByColumnNames(samplesByName.get(name));
			dataMatrix.appendColumns(data.submatrixByColumnNames(negatives));
			dataMatrix.writeGCT(save+"."+name+".data.gct");
			
			//Step 0: Preprocess the results
			MatrixWithHeaders preprocessedDataMatrix=preprocess(data, negatives, samplesByName.get(name), preprocess);
			
			
			DifferentialExpression diff=new DifferentialExpression(preprocessedDataMatrix, negatives, samplesByName.get(name));
			MatrixWithHeaders fdr=diff.getFDRMatrix();
			
			MatrixWithHeaders test=diff.getTestStatisticMatrix();
			
			//write test stats
			test.setPIDToName(annotations);
			test.writeGCT(save+"."+name+".testStat.gct");
			
			//Write FDRs
			fdr.setPIDToName(annotations);
			fdr.writeGCT(save+"."+name+".fdr");
			
			//write permutation test stats
			FileWriter writer=new FileWriter(save+"."+name+".perms.gct");
			MatrixWithHeaders[] perms=diff.getPermutationMatrix();
			for(String row: perms[0].getRowNames()){
				writer.write(row+"\t"+row);
				for(int i=0; i<perms.length; i++){
					writer.write("\t"+perms[i].getRow(row)[0]);
				}
				writer.write("\n");
			}
			writer.close();
			
			//write FDR test stats
			writer=new FileWriter(save+"."+name+".FDRTest.gct");
			for(String row: fdr.getRowNames()){
				writer.write(row+"\t"+fdr.getRow(row)[0]+"\t"+test.getRow(row)[0]+"\n");
			}	
			writer.close();
			
			Collection<String> significantGenes=getSignificantGenes(fdr);
			
			MatrixWithHeaders filteredData=data.submatrixByRowNames(significantGenes);
			MatrixWithHeaders testStat=diff.getTestStatisticMatrix().submatrixByRowNames(significantGenes);
			testStat.setPIDToName(annotations);
			filteredData.setPIDToName(annotations);
			
			Map<String, Collection<String>> geneSet=getGeneSets(testStat, name);
			geneSets.putAll(geneSet);
			
			//filteredData.writeGCT(save+"."+name+".significantGenes.gct");
			//testStat.writeGCT(save+"."+name+".testStatistic.gct");
			genes.addAll(filteredData.getRowNames());
			Collection<String> diffSet=new TreeSet<String>();
			diffSet.addAll(filteredData.getRowNames());
			diffGenes.put(name, diffSet);
		}
		
		GMTParser.writeGMT(save+".gmt", geneSets);
	
		return diffGenes;
	}*/
	
	//TODO Split by up and down
	//TODO Use gene names not PIDs
	private Map<String, Collection<String>> getGeneSets(MatrixWithHeaders testStat, String name) {
		Map<String, Collection<String>> rtrn=new TreeMap<String, Collection<String>>();
		
		Collection<String> up=new TreeSet<String>();
		Collection<String> down=new TreeSet<String>();
		
		for(String PID: testStat.getRowNames()){
			double val=testStat.getRow(PID)[0];
			if(val<0){down.add(testStat.getPIDToName().get(PID));}
			else if(val>0){up.add(testStat.getPIDToName().get(PID));}
		}
		
		String upName=name+"_UP";
		String downName=name+"_DOWN";
		
		rtrn.put(upName, up);
		rtrn.put(downName, down);
		
		return rtrn;
	}

	private Collection<String> getSignificantGenes(MatrixWithHeaders fdr) {
		Collection<String> rtrn=new TreeSet<String>();
		
		for(String row: fdr.getRowNames()){
			double[] fdrs=fdr.getRow(row);
			double min=Statistics.min(fdrs);
			if(min<alpha){rtrn.add(row);}
		}
		
		return rtrn;
	}
	
	private MatrixWithHeaders preprocess(MatrixWithHeaders data, Collection<String> group1, Collection<String> group2, boolean preprocess) {
		if(!preprocess){return data;}
		Collection<String> columnNames=new TreeSet<String>();
		columnNames.addAll(group1);
		columnNames.addAll(group2);
		
		//Get submatrix by group names
		MatrixWithHeaders submatrix=data.submatrixByColumnNames(columnNames);
		
		//compute fold change
		Collection<String> passingRows=getPassingRows(submatrix, group1, group2, cutoff, preprocessByGroups);
		
		//get all passing rows
		MatrixWithHeaders rtrn=data.submatrixByRowNames(passingRows);
		
		System.err.println("Preprocessed from "+data.rowDimension()+" To "+rtrn.rowDimension());
		
		return rtrn;
	}

	private Collection<String> getPassingRows(MatrixWithHeaders submatrix, Collection<String> group1, Collection<String> group2, double cutoff, boolean preprocessByGroups) {
		Collection<String> rtrn=new TreeSet<String>();
		
		for(String row: submatrix.getRowNames()){
			
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
		if(args.length>6){
		File gctFile=new File(args[0]);
		File experimentInfo=new File(args[1]);
		File chipFile=new File(args[2]);
		String save=args[3];
		boolean preprocess=new Boolean(args[4]);
		double alpha=new Double(args[5]);
		int perm=new Integer(args[6]);
		new RunDifferentialExpression(gctFile, experimentInfo, chipFile, save, preprocess, alpha, perm);
		}
		else{System.err.println(usage);}
	}

	static String usage=" args[0]=gct file \n args[1]=experimentInfo \n args[2]=chip file \n args[3]=save \n args[4]=preprocess \n args[5]=alpha \n args[6]=numPerm";
}
