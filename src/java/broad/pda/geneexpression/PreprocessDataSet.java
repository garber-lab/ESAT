package broad.pda.geneexpression;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.core.datastructures.MatrixWithHeaders;
import broad.core.error.ParseException;
import broad.core.math.Statistics;
import broad.core.util.GMTParser;
import broad.pda.geneexpression.agilent.AgilentUtils;

public class PreprocessDataSet {

	double foldChange=2.0;
	
	public PreprocessDataSet(File gctFile, File experimentInfoFile, String save) throws IOException, ParseException{
		MatrixWithHeaders data=new MatrixWithHeaders(gctFile.getAbsolutePath());
		Map<String, ExpressionExperimentInfo> experimentInfo=AgilentUtils.parseExperimentInfoFile(experimentInfoFile);
		
		writePassingSamples(save, data, experimentInfo);
		
		//TODO: write cls file from the experiment info file
		writeCLSFileByGroup(save+".group.cls", data, experimentInfo);
		writeCLSFileBySample(save+".sample.cls", data, experimentInfo);
		
			
		//TODO: Preprocess the results by certain criteria
	}
	
	private void splitFiles(String save, MatrixWithHeaders data, Map<String, ExpressionExperimentInfo> experimentInfo) throws IOException {
		Map<String, String> negativeColumns=new TreeMap<String, String>();
		for(String name: data.getColumnNames()){
			ExpressionExperimentInfo info=experimentInfo.get(name);
			if(info.getSampleType().equalsIgnoreCase("control") && info.passedQC){negativeColumns.put(name, info.getExperimentName());}
		}
		data.writeGCT(save+".negatives.gct", negativeColumns);
	}

	public PreprocessDataSet(File[] agilentFiles, File experimentInfoFile, String save, double foldChange, File chipFile) throws IOException, ParseException{
		Map<String, String> mapping=null;
		if(chipFile!=null){
			mapping=GMTParser.parseCHIPFile(chipFile);
		}
		
		this.foldChange=foldChange;
		Map<String, ExpressionExperimentInfo> experimentInfo=AgilentUtils.parseExperimentInfoFile(experimentInfoFile);
		
		MatrixWithHeaders data=AgilentUtils.makeDataMatrix(agilentFiles, experimentInfo);
		
		//System.err.println("Replacing -999 flags...");
		globalFloor(data);
		
		
		//Filter all genes not differing by at least n fold
		data=filterInvariantGenes(data);
		
		/*norm=filterInvariantGenes(norm);
		quantile=filterInvariantGenes(quantile);*/
		
		MatrixWithHeaders norm=data.scaleNorm();	
		MatrixWithHeaders quantile=data.quantileNormalize();
		
		
		
		
		
		//Pull out negative controls and replicate experiments
		//writeDataSubsets(save+".raw", data, experimentInfo);
		//writeDataSubsets(save+".scaled", norm, experimentInfo);
		//writeDataSubsets(save+".quantile", quantile, experimentInfo);
		
		
		//write GCT files and CLS files
		//writeCLSFileByGroup(save+".group.cls", data, experimentInfo);
		//writeCLSFileBySample(save+".group.cls", data, experimentInfo);
				
		//data.writeBox(save+".raw.box");
		//quantile.writeBox(save+".quantile.box");
		
		if(mapping!=null){
			data=data.submatrixByRowNames(mapping.keySet());
			quantile=quantile.submatrixByRowNames(mapping.keySet());
			norm=norm.submatrixByRowNames(mapping.keySet());
			data.setPIDToName(mapping);
			quantile.setPIDToName(mapping);
			norm.setPIDToName(mapping);
		}
		
		data.writeGCT(save+".raw.gct");
		norm.writeGCT(save+".scaled.gct");
		quantile.writeGCT(save+".quantile.gct");
		
	}

	

	private void globalFloor(MatrixWithHeaders data) {
		
		//Go through all values and get the maximum of the column minimums
		double maxMin=-Double.MAX_VALUE;
		for(String column: data.getColumnNames()){
			double[] vals=data.getColumn(column);
			double min=this.minNonFlagged(vals);
			System.err.println(column+" "+min);
			maxMin=Math.max(maxMin, min);
		}
		
		//replace all -999 with this number
		//replace all values less than this with this number
		for(String column: data.getColumnNames()){
			for(String row: data.getRowNames()){
				double val=data.get(row, column);
				if(val<maxMin){data.set(row, column, maxMin);}
			}
		}
		
		
	}

	private void floor(MatrixWithHeaders quantile) {
		//Find min non -999 value
		double[] vals1=quantile.getColumn(0);
		double min=minNonFlagged(vals1);
		
		for(String column: quantile.getColumnNames()){
			System.err.println(column+" "+minNonFlagged(quantile.getColumn(column)));
		}
		
		//Replace all -999 with min non-flagged
		/*for(String column: quantile.getColumnNames()){
			double[] vals=quantile.getColumn(column);
			for(int i=0; i<vals.length; i++){
				if(vals[i]==-999){quantile.set(i, column, min);}
			}
		}*/
		
	}

	private double minNonFlagged(double[] vals) {
		double min=Double.MAX_VALUE;
		
		for(int i=0; i<vals.length; i++){
			if(vals[i]!=-999){
				min=Math.min(min, vals[i]);
			}
		}
		//System.err.println("Minimum non-zero value is "+min);
		return min;
	}

	private void writeDataSubsets(String save, MatrixWithHeaders data, Map<String, ExpressionExperimentInfo> experimentInfo) throws IOException {
		Map<String, Collection<String>> samplesByGroup=new TreeMap();
		Map<String, Collection<String>> samplesByName=new TreeMap();
				
		for(String barcode: data.getColumnNames()){
			ExpressionExperimentInfo info=experimentInfo.get(barcode);
			String group=info.getSampleType();
			String name=info.getExperimentName();
			Collection<String> groupSamples=new ArrayList();
			Collection<String> nameSamples=new ArrayList();
			if(samplesByGroup.containsKey(group)){groupSamples=samplesByGroup.get(group);}
			if(samplesByName.containsKey(name)){nameSamples=samplesByName.get(name);}
			groupSamples.add(barcode);
			nameSamples.add(barcode);
			samplesByGroup.put(group, groupSamples);
			samplesByName.put(name, nameSamples);
		}
		
		writeSubsets(save, data, samplesByGroup, samplesByName);
	}

	private void writeSubsets(String save, MatrixWithHeaders data, Map<String, Collection<String>> samplesByGroup, Map<String, Collection<String>> samplesByName) throws IOException {
		Collection<String> negatives=samplesByGroup.get("Control");
		
		for(String groupName: samplesByGroup.keySet()){
			if(!groupName.equals("Control")){
				Collection<String> subset=samplesByGroup.get(groupName);
				data.writeCLS(save+"."+groupName+".cls", subset, negatives);
				subset.addAll(negatives);
				data.writeGCT(save+"."+groupName+".gct", subset);
				
			}
		}
		
		for(String sampleName: samplesByName.keySet()){
			Collection<String> subset=samplesByName.get(sampleName);
			data.writeCLS(save+"."+sampleName+".cls", subset, negatives);
			subset.addAll(negatives);
			data.writeGCT(save+"."+sampleName+".gct", subset);		
		}
	
	}

	private MatrixWithHeaders filterInvariantGenes(MatrixWithHeaders data) {
		List<String> rows=new ArrayList<String>();
		
		for(String row: data.getRowNames()){
			double[] vals=data.getRow(row);
			double foldChange=Statistics.fold(vals, true);
			if(foldChange>this.foldChange){rows.add(row);}
		}
		
		return data.submatrixByRowNames(rows);
	}

	private void writeCLSFileBySample(String save, MatrixWithHeaders data, Map<String, ExpressionExperimentInfo> experimentInfo) throws IOException {
		Collection<String> columns=new ArrayList<String>();
		Collection<String> groups=new TreeSet<String>();
		for(String columnName : data.getColumnNames()) {
			ExpressionExperimentInfo info=experimentInfo.get(columnName);
			if(info.passedQC()){columns.add(columnName); groups.add(info.getExperimentName());}
		}
		
		FileWriter writer=new FileWriter(save);
		
		Map<String, Integer> clsMap=new TreeMap<String, Integer>();
		
		writer.write(columns.size()+"\t"+groups.size()+"\t"+1+"\n");
		writer.write("#");
		int i=0;
		for(String group: groups){
			if(!clsMap.containsKey(group)){clsMap.put(group, i); i++;}
			writer.write("\t"+group);
		}
		writer.write("\n");
		
		for(String column: columns){
			ExpressionExperimentInfo info=experimentInfo.get(column);
			writer.write(clsMap.get(info.getExperimentName())+"\t");
		}
		writer.write("\n");
		
		writer.close();
		
	}

	private void writeCLSFileByGroup(String save, MatrixWithHeaders data, Map<String, ExpressionExperimentInfo> experimentInfo) throws IOException {
		Collection<String> columns=new ArrayList<String>();
		Collection<String> groups=new TreeSet<String>();
		for(String columnName : data.getColumnNames()) {
			ExpressionExperimentInfo info=experimentInfo.get(columnName);
			if(info.passedQC()){columns.add(columnName); groups.add(info.getSampleType());}
		}
		
		FileWriter writer=new FileWriter(save);
		
		Map<String, Integer> clsMap=new TreeMap<String, Integer>();
		
		writer.write(columns.size()+"\t"+groups.size()+"\t"+1+"\n");
		writer.write("#");
		int i=0;
		for(String group: groups){
			if(!clsMap.containsKey(group)){clsMap.put(group, i); i++;}
			writer.write("\t"+group);
		}
		writer.write("\n");
		
		for(String column: columns){
			ExpressionExperimentInfo info=experimentInfo.get(column);
			writer.write(clsMap.get(info.getSampleType())+"\t");
		}
		writer.write("\n");
		
		writer.close();
	}

	private void writePassingSamples(String save, MatrixWithHeaders data, Map<String, ExpressionExperimentInfo> experimentInfo) throws IOException {
		
		
		Collection<String> columns=new ArrayList<String>();
		for(String columnName : data.getColumnNames()) {
			ExpressionExperimentInfo info=experimentInfo.get(columnName);
			if(info.passedQC()){columns.add(columnName);}
		}
		
		data.writeGCT(save, columns);
	}
	
	
	public static void main(String[] args)throws IOException, ParseException{
		if(args.length>3){
			File[] gctFiles=new File(args[0]).listFiles();
			File experimentInfo=new File(args[1]);
			String save=args[2];
			double foldChange=new Double(args[3]);
			
			File chipFile=null;
			if(args.length>4){
				chipFile=new File(args[4]);
			}
			
			System.err.println("Agilent files: "+gctFiles.length);
			new PreprocessDataSet(gctFiles, experimentInfo, save, foldChange, chipFile);
			
			//TODO: Add flags to support, 1)making cls files, 2) extracting experimentNames, 3) gene names or probe names, 4) filter genes
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=agilent raw files \n args[1]=experimentinfo file \n args[2]=save \n args[3]=fold change \n args[4]=chip file";
}
