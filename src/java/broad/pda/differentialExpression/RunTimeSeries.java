package broad.pda.differentialExpression;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.core.datastructures.MatrixWithHeaders;
import broad.core.error.ParseException;
import broad.core.math.Statistics;
import broad.core.util.ParseGCTFile;
import broad.pda.geneexpression.ExpressionExperimentInfo;
import broad.pda.geneexpression.agilent.AgilentUtils;

public class RunTimeSeries {

	//This class takes an expression data set that is a time series
	//and tries to identify genes that are differentially expressed across timed conditions
	boolean preprocess=true;
	private boolean preprocessByGroups=false;
	private double cutoff=2.0;
	Map<String, String> annotations;
	private double alpha=.05;;
	
	public RunTimeSeries(File gctFile, File experimentInfoFile, String save, boolean preprocess, File chipFile, double cutoff) throws IOException, ParseException{
		this.preprocess=preprocess;
		MatrixWithHeaders data=new MatrixWithHeaders(gctFile.getAbsolutePath(), chipFile.getAbsolutePath());
		Map<String, ExpressionExperimentInfo> experimentInfo=AgilentUtils.parseExperimentInfoFile(experimentInfoFile);
		this.annotations=ParseGCTFile.parseChipFile(chipFile);
		
		//Idea 1: Paired t-stat overall points
		//Idea 2: Paired t-stat over subset
		//Idea 3: Group t-stat for individual time points against control
		//Idea 4: subsets against all controls
		
		//TODO: Run the subset comparisons as paired t-tests
		runComparisons(data, experimentInfo, save, preprocess, cutoff);
		//runPairedComparisons(data, experimentInfo, save, preprocess, cutoff);
	}

	private void runPairedComparisons(MatrixWithHeaders data,Map<String, ExpressionExperimentInfo> experimentInfo, String save, boolean preprocess, double cutoff) throws IOException {
		Map<String, Collection<String>> samplesByName=this.getSamplesByName(data, experimentInfo);
		Map<String, Collection<String>> gmtFile=new TreeMap();
		
		Collection<String> negatives=samplesByName.remove("Control");
		Collection<String> allPassingGenes=new TreeSet<String>();
		
		for(String name: samplesByName.keySet()){
			//TODO Get paired samples
			Collection<String> samples=samplesByName.get(name);
			//samples.addAll(negatives);
			//Step 0: Preprocess the results
			MatrixWithHeaders preprocessedDataMatrix=preprocess(data, negatives, samplesByName.get(name), preprocess);
			//preprocessedDataMatrix=preprocessedDataMatrix.submatrixByColumnNames(samples);
			System.err.println(name+" "+samplesByName.get(name));
			Map<String, Collection<String>[]> timeSubsets=this.getPairedLinearSubsets(samples, negatives, experimentInfo);
			
			MatrixWithHeaders fdrMatrix=new MatrixWithHeaders(preprocessedDataMatrix.getRowNames(), new ArrayList(timeSubsets.keySet()));
			MatrixWithHeaders statMatrix=new MatrixWithHeaders(preprocessedDataMatrix.getRowNames(), new ArrayList(timeSubsets.keySet()));
			MatrixWithHeaders foldMatrix=new MatrixWithHeaders(preprocessedDataMatrix.getRowNames(), new ArrayList(timeSubsets.keySet()));
			
			for(String timeSubset: timeSubsets.keySet()){
				System.err.println(name+" "+timeSubset);
				Collection<String>[] subsets=timeSubsets.get(timeSubset);
				
				MatrixWithHeaders dataMatrix=data.submatrixByColumnNames(samplesByName.get(name));
				dataMatrix.appendColumns(data.submatrixByColumnNames(negatives));
				DifferentialExpression diff=new DifferentialExpression(preprocessedDataMatrix, subsets[0], subsets[1], false, true);
				MatrixWithHeaders fdr=diff.getFDRMatrix();
				//fdr.setPIDToName(annotations);
				//fdr.writeGCT(save+"."+name+"."+timeSubset+".paired.fdr");
				
				MatrixWithHeaders stats=diff.getTestStatisticMatrix();
				fdrMatrix.setColumn(fdr.getColumn(0), timeSubset);
				statMatrix.setColumn(stats.getColumn(0), timeSubset);
				foldMatrix.setColumn(diff.getFoldMatrix().getColumn(0), timeSubset);
			}
			fdrMatrix.setPIDToName(annotations);
			fdrMatrix.writeGCT(save+"."+name+".paired.fdr");
			statMatrix.setPIDToName(annotations);
			statMatrix.writeGCT(save+"."+name+".paired.stat");
			Collection<String> passingGenes=this.getSignificantGenes(fdrMatrix, foldMatrix, cutoff);
			allPassingGenes.addAll(passingGenes);
			//MatrixWithHeaders passingSubset=preprocessedDataMatrix.submatrixByRowNames(passingGenes);
			//passingSubset.setPIDToName(annotations);
			//passingSubset.writeGCT(save+"."+name+".gct");
			Map<String, Collection<String>> geneSets=getGeneSets(fdrMatrix, statMatrix, foldMatrix, cutoff, name);
			gmtFile.putAll(geneSets);
			writeGMT(save+"."+name+".paired.gmt", fdrMatrix, statMatrix, foldMatrix, cutoff);
		}
		
		MatrixWithHeaders subset=data.submatrixByRowNames(allPassingGenes);
		subset.setPIDToName(annotations);
		subset.writeGCT(save+".diff.gct");
		ParseGCTFile.writeGMTFile(save+".paired.gmt", gmtFile);
		
		//TODO Write gene set file for all time points
		
		//return diffGenes;
		
		
	}

	private Map<String, Collection<String>> getGeneSets(MatrixWithHeaders fdrMatrix,	MatrixWithHeaders statMatrix, MatrixWithHeaders foldMatrix, double cutoff2, String name) {
		Map<String, Collection<String>> rtrn=new TreeMap();
		
		for(String experiment: fdrMatrix.getColumnNames()){
			Collection<String> up=getUp(fdrMatrix, statMatrix, foldMatrix, experiment, cutoff);
			Collection<String> down=getDown(fdrMatrix, statMatrix, foldMatrix, experiment, cutoff);
			String upName=name+"_"+experiment+"_UP";
			String downName=name+"_"+experiment+"_DOWN";
			rtrn.put(downName, down);
			rtrn.put(upName, up);
		}
		
		return rtrn;
	}

	private Map<String, Collection<String>> getLinearSubsets(Collection<String> names, Map<String, ExpressionExperimentInfo> experimentInfo) {
		Map<Integer, Collection<String>> indexByTime=indexByTime(names, experimentInfo);
		int[] times=getTimes(indexByTime);
		int[] subsetSizes=getPossibleSubsetSizes(times);
		
		Collection<List<Integer>> subsetOfTimes=new HashSet<List<Integer>>();
		
		for(int i=0; i<subsetSizes.length; i++){
			int size=subsetSizes[i];
			Collection<List<Integer>> subsets=getAllSubsets(times, size);
			subsetOfTimes.addAll(subsets);
		}
		
		Map<String, Collection<String>> groups=getGroups(subsetOfTimes, indexByTime);
		return groups;
	}
	
	private Map<String, Collection<String>[]> getPairedLinearSubsets(Collection<String> samples, Collection<String> controls, Map<String, ExpressionExperimentInfo> experimentInfo) {
		Map<Integer, Collection<String>> indexByTime=indexByTime(samples, experimentInfo);
		Map<Integer, Collection<String>> controlsByTime=indexByTime(controls, experimentInfo);
		
		int[] times=getTimes(indexByTime);
		int[] subsetSizes=getPossibleSubsetSizes(times);
		
		Collection<List<Integer>> subsetOfTimes=new HashSet<List<Integer>>();
		
		for(int i=0; i<subsetSizes.length; i++){
			int size=subsetSizes[i];
			Collection<List<Integer>> subsets=getAllSubsets(times, size);
			subsetOfTimes.addAll(subsets);
		}
		
		Map<String, Collection<String>[]> groups=getGroups(subsetOfTimes, indexByTime, controlsByTime);
		return groups;
	}
	
	private Map<String, Collection<String>[]> getGroups(Collection<List<Integer>> subsetOfTimes,Map<Integer, Collection<String>> indexByTime,Map<Integer, Collection<String>> controlsByTime) {
		Map<String, Collection<String>[]> rtrn=new TreeMap<String, Collection<String>[]>();
		
		for(List<Integer> subset: subsetOfTimes){
			String name=makeName(subset);
			Collection<String> group=getGroupAvg(subset, indexByTime);
			Collection<String> controls=getGroupAvg(subset, controlsByTime);
			Collection<String>[] array=new Collection[2];
			array[0]=controls;
			array[1]=group;
			rtrn.put(name, array);
		}
		
		return rtrn;
	}

	private Collection<String> getGroupAvg(List<Integer> subset, Map<Integer, Collection<String>> indexByTime) {
		Collection<String> rtrn=new ArrayList<String>();
		
		for(Integer index: subset){
			Collection<String> vals=indexByTime.get(index);
			if(vals.size()>1){System.err.println("WARN: Had more than one replicate at time "+index+" so took the first");}
			rtrn.add(vals.iterator().next());
		}
		
		return rtrn;
	}

	//slide a window of fixed size, forward and reverse across times
	private static Collection<List<Integer>> getAllSubsets(int[] times, int size) {
		Collection<List<Integer>> rtrn=new ArrayList();
		
		//Forward
		for(int i=0; i<times.length-size; i++){
			int start=i;
			int end=i+size;
			List<Integer> subset=getSubset(start, end, times);
			rtrn.add(subset);
		}
		
		//Reverse
		for(int i=times.length; i>=size; i--){
			int start=i-size;
			int end=i;
			List<Integer> subset=getSubset(start, end, times);
			rtrn.add(subset);
		}
		
		
		return rtrn;
	}

	private static List<Integer> getSubset(int start, int end, int[] times) {
		List<Integer> rtrn=new ArrayList();
		for(int i=start; i<end; i++){
			rtrn.add(times[i]);
		}
		return rtrn;
	}

	private Map<String, Collection<String>> getGroups(Collection<List<Integer>> subsetOfTimes, Map<Integer, Collection<String>> indexByTime) {
		Map<String, Collection<String>> rtrn=new TreeMap<String, Collection<String>>();
		
		for(List<Integer> subset: subsetOfTimes){
			String name=makeName(subset);
			Collection<String> group=getGroup(subset, indexByTime);
			rtrn.put(name, group);
		}
		
		return rtrn;
	}

	private Collection<String> getGroup(List<Integer> subset, Map<Integer, Collection<String>> indexByTime) {
		Collection<String> rtrn=new ArrayList<String>();
		
		for(Integer index: subset){
			Collection<String> vals=indexByTime.get(index);
			rtrn.addAll(vals);
		}
		
		return rtrn;
	}

	private String makeName(List<Integer> subset) {
		String rtrn="t";
	
		for(Integer time: subset){rtrn+=time+"_";}
		
		return rtrn;
	}

	private void runComparisons(MatrixWithHeaders data,	Map<String, ExpressionExperimentInfo> experimentInfo, String save, boolean preprocess, double cutoff) throws IOException {
		Map<String, Collection<String>> samplesByName=this.getSamplesByName(data, experimentInfo);
		Map<String, Collection<String>> gmtFile=new TreeMap();
		Collection<String> negatives=samplesByName.remove("Control");
		Collection<String> allPassingGenes=new TreeSet<String>();
		
		for(String name: samplesByName.keySet()){
			//TODO Get paired samples
			Collection<String> samples=samplesByName.get(name);
			//samples.addAll(negatives);
			//Step 0: Preprocess the results
			MatrixWithHeaders preprocessedDataMatrix=preprocess(data, negatives, samplesByName.get(name), preprocess);
			//preprocessedDataMatrix=preprocessedDataMatrix.submatrixByColumnNames(samples);
			System.err.println(name+" "+samplesByName.get(name));
			Map<String, Collection<String>> timeSubsets=this.getLinearSubsets(samplesByName.get(name), experimentInfo);
			
			MatrixWithHeaders fdrMatrix=new MatrixWithHeaders(preprocessedDataMatrix.getRowNames(), new ArrayList(timeSubsets.keySet()));
			MatrixWithHeaders statMatrix=new MatrixWithHeaders(preprocessedDataMatrix.getRowNames(), new ArrayList(timeSubsets.keySet()));
			MatrixWithHeaders foldMatrix=new MatrixWithHeaders(preprocessedDataMatrix.getRowNames(), new ArrayList(timeSubsets.keySet()));
			
			for(String timeSubset: timeSubsets.keySet()){
				System.err.println(name+" "+timeSubset);
				Collection<String> subset=timeSubsets.get(timeSubset);
				MatrixWithHeaders dataMatrix=data.submatrixByColumnNames(samplesByName.get(name));
				dataMatrix.appendColumns(data.submatrixByColumnNames(negatives));
				DifferentialExpression diff=new DifferentialExpression(preprocessedDataMatrix, negatives, subset);
				MatrixWithHeaders fdr=diff.getFDRMatrix();
				//fdr.setPIDToName(annotations);
				//fdr.writeGCT(save+"."+name+"."+timeSubset+".fdr");
				
				MatrixWithHeaders stats=diff.getTestStatisticMatrix();
				fdrMatrix.setColumn(fdr.getColumn(0), timeSubset);
				statMatrix.setColumn(stats.getColumn(0), timeSubset);
				foldMatrix.setColumn(diff.getFoldMatrix().getColumn(0), timeSubset);
			}
			fdrMatrix.setPIDToName(annotations);
			fdrMatrix.writeGCT(save+"."+name+".fdr");
			statMatrix.setPIDToName(annotations);
			statMatrix.writeGCT(save+"."+name+".stat");
			Collection<String> passingGenes=this.getSignificantGenes(fdrMatrix, foldMatrix, cutoff);
			allPassingGenes.addAll(passingGenes);
			MatrixWithHeaders passingSubset=preprocessedDataMatrix.submatrixByRowNames(passingGenes);
			passingSubset.setPIDToName(annotations);
			passingSubset.writeGCT(save+"."+name+".gct");
			Map<String, Collection<String>> geneSets=getGeneSets(fdrMatrix, statMatrix, foldMatrix, cutoff, name);
			gmtFile.putAll(geneSets);
			writeGMT(save+"."+name+".gmt", fdrMatrix, statMatrix, foldMatrix, cutoff);
		}
		
		MatrixWithHeaders subset=data.submatrixByRowNames(allPassingGenes);
		subset.setPIDToName(annotations);
		subset.writeGCT(save+".diff.gct");
		ParseGCTFile.writeGMTFile(save+".gmt", gmtFile);
		//TODO Write gene set file for all time points
		
		//return diffGenes;
	}
	
	private void writeGMT(String save, MatrixWithHeaders fdrMatrix, MatrixWithHeaders statMatrix, MatrixWithHeaders foldMatrix, double cutoff) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String experiment: fdrMatrix.getColumnNames()){
			Collection<String> up=getUp(fdrMatrix, statMatrix, foldMatrix, experiment, cutoff);
			Collection<String> down=getDown(fdrMatrix, statMatrix, foldMatrix, experiment, cutoff);
			writer.write(experiment+"_UP\t"+up.size());
			for(String gene: up){writer.write("\t"+gene);}
			writer.write("\n");
			writer.write(experiment+"_DOWN\t"+down.size());
			for(String gene: down){writer.write("\t"+gene);}
			writer.write("\n");
		}
		
		writer.close();
	}

	private Collection<String> getUp(MatrixWithHeaders fdrMatrix, MatrixWithHeaders statMatrix, MatrixWithHeaders foldMatrix, String experiment, double cutoff) {
		Collection<String> rtrn=new TreeSet<String>();
	
		for(String name: fdrMatrix.getRowNames()){
			double fdr=fdrMatrix.get(name, experiment);
			double stat=statMatrix.get(name, experiment);
			double fold=(foldMatrix.get(name, experiment));
			if(fdr<alpha && stat>0 && fold>cutoff){
				rtrn.add(annotations.get(name));
			}
		}
		
		return rtrn;
	}

	private Collection<String> getDown(MatrixWithHeaders fdrMatrix, MatrixWithHeaders statMatrix, MatrixWithHeaders foldMatrix, String experiment, double cutoff) {
		Collection<String> rtrn=new TreeSet<String>();
	
		for(String name: fdrMatrix.getRowNames()){
			double fdr=fdrMatrix.get(name, experiment);
			double stat=statMatrix.get(name, experiment);
			double fold=(foldMatrix.get(name, experiment));
			if(fdr<alpha && stat<0 && fold>cutoff){rtrn.add(annotations.get(name));}
		}
		
		return rtrn;
	}

	private Collection<String> getSignificantGenes(MatrixWithHeaders fdr, MatrixWithHeaders foldMatrix, double cutoff2) {
		Collection<String> rtrn=new TreeSet<String>();
		
		for(String row: fdr.getRowNames()){
			double[] fdrs=fdr.getRow(row);
			double[] folds=foldMatrix.getRow(row);
			double min=Statistics.min(fdrs);
			double max=Statistics.max(folds);
			if(min<alpha && max >cutoff2){rtrn.add(row);}
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
	
	private Map<String, Collection<String>> getSamplesByName(MatrixWithHeaders data, Map<String, ExpressionExperimentInfo> experimentInfo){
		Map<String, Collection<String>> samplesByName=new TreeMap<String, Collection<String>>();
		Collection<String> negatives=new TreeSet<String>();	
		
		for(String barcode: data.getColumnNames()){
			if(experimentInfo.containsKey(barcode)){
				ExpressionExperimentInfo info=experimentInfo.get(barcode);
				String group=info.getSampleType();
				if(!group.equalsIgnoreCase("Control")){
					//System.err.println(barcode+" "+group);
					String name=info.getExperimentName();
					Collection<String> nameSamples=new ArrayList();
					if(samplesByName.containsKey(name)){nameSamples=samplesByName.get(name);}
					nameSamples.add(barcode);
					samplesByName.put(name, nameSamples);
				}
				else{negatives.add(barcode);}
			}
		}
		
		samplesByName.put("Control", negatives);
		return samplesByName;
	}

	private static int[] getPossibleSubsetSizes(int[] times) {
		int[] rtrn=new int[times.length-1];
		
		for(int i=0; i<rtrn.length; i++){
			rtrn[i]=(i+2);
			System.err.println(rtrn[i]);
		}
		
		return rtrn;
	}

	private int[] getTimes(Map<Integer, Collection<String>> indexByTime) {
		int[] rtrn=new int[indexByTime.size()];
		
		int counter=0;
		for(Integer val: indexByTime.keySet()){rtrn[counter++]=val;}
		
		return rtrn;
	}

	private Map<Integer, Collection<String>> indexByTime(Collection<String> names, Map<String, ExpressionExperimentInfo> experimentInfo) {
		Map<Integer, Collection<String>> index=new TreeMap<Integer, Collection<String>>();
		
		for(String name: names){
			int time=experimentInfo.get(name).getTimePoint();
			Collection<String> list=new ArrayList<String>();
			if(index.containsKey(time)){list=index.get(time);}
			list.add(name);
			index.put(time, list);
		}
		return index;
	}
	
	public static void main(String[] args) throws IOException, ParseException{
		if(args.length>5){
			File gctFile=new File(args[0]);
			File experimentInfo=new File(args[1]);
			File chipFile=new File(args[2]);
			String save=args[3];
			boolean preprocess=new Boolean(args[4]);
			double cutoff=new Double(args[5]);
			new RunTimeSeries(gctFile, experimentInfo, save, preprocess, chipFile, cutoff);
		}
		else{System.err.println(usage);}
	}
	static String usage=" args[0]=gct file \n args[1]=experimentInfo \n args[2]=chip file \n args[3]=save \n args[4]=preprocess \n args[5]=fold cutoff";
}
