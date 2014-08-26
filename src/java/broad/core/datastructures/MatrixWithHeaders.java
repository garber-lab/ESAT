package broad.core.datastructures;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;



import broad.core.error.ParseException;
import broad.core.math.MathUtil;
import broad.core.math.Statistics;
import broad.core.util.ParseGCTFile;
import broad.pda.geneexpression.ExpressionExperimentInfo;


import Jama.Matrix;

public class MatrixWithHeaders {
	private Matrix data;
	private LinkedHashMap<String, Integer> columnIndexMap;
	private LinkedHashMap<String, Integer> rowIndexMap;
	private LinkedHashMap<String, List<Integer>> rowDescrIndexMap;
	private List<String> columns;
	private List<String> rows;
	private Map<String, String> pidToName;
	private TreeMap<String, String> probeClass;
	private TreeMap<String, Collection<Integer>> replicateMap;
	
	protected MatrixWithHeaders() {
		super();
	}
	
	
	public MatrixWithHeaders(String inputFile)  throws IOException, ParseException {
		BufferedReader br = new BufferedReader(new FileReader(inputFile));
		initFromReader(br);
		br.close();
	}
	
	
	public MatrixWithHeaders(BufferedReader br) throws IOException, ParseException {
		initFromReader(br);
	}
	
	public MatrixWithHeaders(Matrix data, List<String> rowNames, List<String> colNames) {
		this.data = data;
		initNameIndexMaps(rowNames, colNames);
	}
	
	public MatrixWithHeaders(Matrix data, List<String> rowNames, List<String> colNames, List<String> rowDescriptions) {
		this.data = data;
		initNameIndexMaps(rowNames, colNames, rowDescriptions);
	}
	
	public MatrixWithHeaders(List<String> rows, List<String> columns, Map<String, String> pidToName) {
		this.data = new Matrix(rows.size(), columns.size());
		this.pidToName=pidToName;
		initNameIndexMaps(rows, columns);
	}
	
	/**
	 * Create a new matrix that has a subset of columns and rows of the original matrix
	 * @param matrix
	 * @param colsInRows
	 * @param colsInRows2
	 */
	public MatrixWithHeaders(MatrixWithHeaders matrix, List<String> rows, List<String> columns) {
		this(rows, columns);
		for(int i = 0; i < rows.size(); i++) {
			for(int j = 0; j < columns.size(); j++) {
				data.set(i,j, matrix.get(rows.get(i), columns.get(j)));
			}
		}
	}
	
	public MatrixWithHeaders(List<String> rows, List<String> columns) {
		this.data = new Matrix(rows.size(), columns.size());
		initNameIndexMaps(rows, columns);
	}
	
	public MatrixWithHeaders(List<String> rows, List<String> columns, List<String> rowDescriptions) {
		this.data = new Matrix(rows.size(), columns.size());
		initNameIndexMaps(rows, columns, rowDescriptions);
	}
	
	
	public MatrixWithHeaders(String gctFile, String chipFile) throws IOException, ParseException {
		BufferedReader br = new BufferedReader(new FileReader(gctFile));
		initFromReader(br);
		br.close();
		
		Map<String, String> annotations=ParseGCTFile.parseChipFile(new File(chipFile));
		this.setPIDToName(annotations);
	}

	
	public MatrixWithHeaders(List<List<Double>> rawData, ArrayList<String> rowNameList,
			List<String> columnNameList) {
		init(rawData, rowNameList, columnNameList);
	}


	public MatrixWithHeaders(ArrayList<String> rows,
			ArrayList<String> columns, int DefMatVal) {
		
		this(rows, columns);
		for(int i = 0; i < rows.size(); i++) {
			for(int j = 0; j < columns.size(); j++) {
				data.set(i,j, DefMatVal);
			}
		}
	}


	public Map<String, String> getPIDToName(){return this.pidToName;}

	public void write(BufferedWriter bw) throws IOException {
		//write("row");
		for(String columnName : columnIndexMap.keySet()) {
			bw.write("\t");
			bw.write(columnName);
		}
		bw.newLine();
		bw.flush();
		for(String rowName : rowIndexMap.keySet()) {
			bw.write(rowName);
			for(String colName : columnIndexMap.keySet()) {
				bw.write("\t");
				bw.write(String.valueOf(get(rowName, colName)));
			}
			bw.newLine();
		}
	}
	
	public void write(String fileName) throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(fileName));
		write(bw);
		bw.close();
		
	}
	
	public void writeGCT(BufferedWriter bw) throws IOException {
		bw.write("#1.2");
		bw.newLine();
		bw.write(data.getRowDimension() +"\t"+data.getColumnDimension());
		bw.newLine();
		
		bw.write("name\tdescription");
		for(String columnName : columnIndexMap.keySet()) {
			bw.write("\t");
			bw.write(columnName);
		}
		bw.newLine();
		for(String rowName : rowIndexMap.keySet()) {
			bw.write(rowName);
			bw.write("\t");
			if(this.pidToName==null || pidToName.isEmpty() || !pidToName.containsKey(rowName)){bw.write(rowName);}
			else{bw.write(pidToName.get(rowName));}
			for(String colName : columnIndexMap.keySet()) {
				bw.write("\t");
				bw.write(String.valueOf(get(rowName, colName)));
			}
			bw.newLine();
		}
	}
	
	public void writeGCT(String save) throws IOException {
		writeGCT(save, this.getColumnNames(), this.getRowNames());
	}
	
	public void writeGCT(String save, Collection<String> columns) throws IOException {
		writeGCT(save, columns, this.getRowNames());
	}
	
	
	public void writeGCT(String save, Map<String, Collection<String>> experimentInfo, boolean ordered) throws IOException {
		FileWriter writer=new FileWriter(save);
		writer.write("#1.2\n");
		
		int columnSize=0;
		for(String group: experimentInfo.keySet()){
			Collection<String> experiment=experimentInfo.get(group);
			for(String e: experiment){
				columnSize++;
			}
		}
		
		writer.write(rows.size() +"\t"+columnSize+"\n");
				
		writer.write("name\tdescription");
		
		
		for(String group: experimentInfo.keySet()){
			Collection<String> experiment=experimentInfo.get(group);
			for(String e: experiment){
				writer.write("\t"+e);
			}
		}
		writer.write("\n");
		
		
		for(String rowName : this.getRowNames()) {
			if(this.pidToName==null || pidToName.isEmpty() || !pidToName.containsKey(rowName)){writer.write(rowName+"\t"+rowName);}
			else{writer.write(rowName+"\t"+pidToName.get(rowName));}
			
			for(String group: experimentInfo.keySet()){
				Collection<String> experiment=experimentInfo.get(group);
				for(String e: experiment){
					writer.write("\t"+get(rowName, e));
				}
			}
			
			writer.write("\n");
		}
		writer.close();
	}
	
	public void writeGCT(String save, Map<String, Collection<String>> experimentInfo, Map<String, ExpressionExperimentInfo> fullInfo, boolean info) throws IOException {
		FileWriter writer=new FileWriter(save);
		writer.write("#1.2\n");
		
		int columnSize=0;
		for(String group: experimentInfo.keySet()){
			Collection<String> experiment=experimentInfo.get(group);
			for(String e: experiment){
				columnSize++;
			}
		}
		
		writer.write(rows.size() +"\t"+columnSize+"\n");
				
		writer.write("name\tdescription");
		
		
		for(String group: experimentInfo.keySet()){
			Collection<String> experiment=experimentInfo.get(group);
			for(String e: experiment){
				writer.write("\t"+fullInfo.get(e).getExperimentName());
			}
		}
		writer.write("\n");
		
		
		for(String rowName : this.getRowNames()) {
			if(this.pidToName==null || pidToName.isEmpty() || !pidToName.containsKey(rowName)){writer.write(rowName+"\t"+rowName);}
			else{writer.write(rowName+"\t"+pidToName.get(rowName));}
			
			for(String group: experimentInfo.keySet()){
				Collection<String> experiment=experimentInfo.get(group);
				for(String e: experiment){
					writer.write("\t"+get(rowName, e));
				}
			}
			
			writer.write("\n");
		}
		writer.close();
	}
	
	public void writeGCT(String save, Collection<String> columns, Collection<String> rows) throws IOException {
		FileWriter writer=new FileWriter(save);
		writer.write("#1.2\n");
		writer.write(rows.size() +"\t"+columns.size()+"\n");
				
		writer.write("name\tdescription");
		for(String columnName : columns) {
			writer.write("\t"+columnName);
		}
		writer.write("\n");
		for(String rowName : rows) {
			if(this.pidToName==null || pidToName.isEmpty() || !pidToName.containsKey(rowName)){writer.write(rowName+"\t"+rowName);}
			else{writer.write(rowName+"\t"+pidToName.get(rowName));}
			
			for(String colName : columns) {
				writer.write("\t"+get(rowName, colName));
			}
			writer.write("\n");
		}
		writer.close();
	}
	
	public void writeGCTWithHeaders(String save, Map<String, String> columnMapping) throws IOException {
		FileWriter writer=new FileWriter(save);
		writer.write("#1.2\n");
		writer.write(rows.size() +"\t"+columns.size()+"\n");
				
		writer.write("name\tdescription");
		for(String columnName : columns) {
			String info=columnMapping.get(columnName);
			writer.write("\t"+info);
		}
		writer.write("\n");
		for(String rowName : rows) {
			if(this.pidToName==null || pidToName.isEmpty() || !pidToName.containsKey(rowName)){writer.write(rowName+"\t"+rowName);}
			else{writer.write(rowName+"\t"+pidToName.get(rowName));}
			
			for(String colName : columns) {
				writer.write("\t"+get(rowName, colName));
			}
			writer.write("\n");
		}
		writer.close();
	}
	
	//SPECIFY THE PRECISION FOR PRINTING DOUBLE
	public void writeGCT(String save, Collection<String> columns, Collection<String> rows, int precision) throws IOException {
		FileWriter writer=new FileWriter(save);
		writer.write("#1.2\n");
		writer.write(rows.size() +"\t"+columns.size()+"\n");
				
		writer.write("name\tdescription");
		for(String columnName : columns) {
			writer.write("\t"+columnName);
		}
		writer.write("\n");
		for(String rowName : rows) {
			if(this.pidToName==null || pidToName.isEmpty() || !pidToName.containsKey(rowName)){writer.write(rowName+"\t"+rowName);}
			else{writer.write(rowName+"\t"+pidToName.get(rowName));}
			
			for(String colName : columns) {
				writer.write("\t");
				String s= new Double(get(rowName, colName)).toString();
				if (s.length()<2+precision) writer.write (s); 
				else writer.write (s,0,2+precision);
			}
			writer.write("\n");
		}
		writer.close();
	}
	
	
	public void writeGCT(String save, Map<String, String> columns, Map<String, String> rows) throws IOException {
		FileWriter writer=new FileWriter(save);
		writer.write("#1.2\n");
		writer.write(rows.size() +"\t"+columns.size()+"\n");
				
		writer.write("name\tdescription");
		for(String columnName : columns.keySet()) {
			writer.write("\t"+columns.get(columnName));
		}
		writer.write("\n");
		for(String rowName : rows.keySet()) {
			if(this.pidToName==null || pidToName.isEmpty() || !pidToName.containsKey(rowName)){writer.write(rowName+"\t"+rowName);}
			else{writer.write(rowName+"\t"+pidToName.get(rowName));}
			for(String colName : columns.keySet()) {
				writer.write("\t"+get(rowName, colName));
			}
			writer.write("\n");
		}
		writer.close();
	}
	
	public void writeGCT(String save, Map<String, String> columns, Collection<String> rows) throws IOException {
		FileWriter writer=new FileWriter(save);
		writer.write("#1.2\n");
		writer.write(rows.size() +"\t"+columns.size()+"\n");
				
		writer.write("name\tdescription");
		int counter=0;
		for(String columnName : columns.keySet()) {
			writer.write("\t"+columns.get(columnName)+"_"+counter);
			counter++;
		}
		writer.write("\n");
		for(String rowName : rows) {
			if(this.pidToName==null || pidToName.isEmpty() || !pidToName.containsKey(rowName)){writer.write(rowName+"\t"+rowName);}
			else{writer.write(rowName+"\t"+pidToName.get(rowName));}
			for(String colName : columns.keySet()) {
				writer.write("\t"+get(rowName, colName));
			}
			writer.write("\n");
		}
		writer.close();
	}
	
	public void writeGCT(String save, Map<String, String> columns) throws IOException {
		writeGCT(save, columns, this.getRowNames());
	}
	
	public void writeGCT(String save, File chipFile) throws IOException {
		Map<String, String> mapping=ParseGCTFile.parseChipFile(chipFile);
		FileWriter writer=new FileWriter(save);
		writer.write("#1.2\n");
		writer.write(data.getRowDimension() +"\t"+data.getColumnDimension()+"\n");
				
		writer.write("name\tdescription");
		for(String columnName : columnIndexMap.keySet()) {
			writer.write("\t"+columnName);
		}
		writer.write("\n");
		for(String rowName : rowIndexMap.keySet()) {
			writer.write(rowName+"\t"+mapping.get(rowName));
			for(String colName : columnIndexMap.keySet()) {
				writer.write("\t"+get(rowName, colName));
			}
			writer.write("\n");
		}
		writer.close();
	}
	
	
	/**
	 * Takes the inverse (or pseudo inverse if the underlying data is not a squared matrix)
	 * @return
	 */
	public Matrix dataInverse() {
		return data.inverse();		
	}
	
	public int columnDimension() { return data.getColumnDimension();}
	public int rowDimension() { return data.getRowDimension();}
	public double get(int i, int j) {return data.get(i,j);}
	
	public double[] getColumn(int j){
		double[] rtrn=new double[rowDimension()];
		for(int i=0; i<rtrn.length; i++){
			rtrn[i]=get(i,j);
		}
		return rtrn;
	}
	
	public double[] getColumn(String columnName){
		double[] rtrn=new double[this.rowDimension()];
		for(int i=0; i<rtrn.length; i++){
			rtrn[i]=get(i,columnName);
		}
		return rtrn;
	}
	
	public double[] getRow(String rowName){
		double[] rtrn=new double[this.columnDimension()];
		//System.err.println(rowName);
		
		for(int j=0; j<rtrn.length; j++){
			rtrn[j]=get(rowName, j);
		}
		
		return rtrn;
	}
	
	public void set(int i, int j, double val) {data.set(i, j, val);}
	public void set(String row, String column, double value) {
		if(!rowIndexMap.containsKey(row)) {
			throw new IllegalArgumentException("Trying to add a value using a row ("+row+") not in the matrix");
		}
		
		if(!columnIndexMap.containsKey(column)) {
			throw new IllegalArgumentException("Trying to add a value using a column ("+column+") not in the matrix");
		}
		data.set(rowIndexMap.get(row), columnIndexMap.get(column), value);
	}
	
	public void setRow(String row, double[] vals){
		for(int i=0; i<vals.length; i++){
			set(row, i, vals[i]);
		}
	}
	
	public void setColumn(double[] vals, String column){
		for(int i=0; i<vals.length; i++){
			set(i, column, vals[i]);
		}
	}
	
	public void setColumn(double[] vals, int column){
		for(int i=0; i<vals.length; i++){
			set(i, column, vals[i]);
		}
	}
	
	public void set(String row, int colIdx, double value) {
		data.set(rowIndexMap.get(row), colIdx, value);
		
	}
	
	public void set(int rowIdx, String column, double value) {
		data.set(rowIdx, columnIndexMap.get(column), value);
		
	}
	
	public double get(String row, String column) {
		if(!rowIndexMap.containsKey(row) ) {
			//System.err.println("Row " + row + " not found");
		}
		
		if(!columnIndexMap.containsKey(column) ) {
			//System.err.println("Column " + column + " not found");
		}
		
		return data.get(rowIndexMap.get(row), columnIndexMap.get(column));
	}
	public boolean containsColumn (String colName) {
		return columnIndexMap.containsKey(colName);
	}
	
	public boolean containsRow (String rowName) {
		return rowIndexMap.containsKey(rowName);
	}
	public double get(String rowName, int colIdx) {return data.get(rowIndexMap.get(rowName), colIdx);}
	public double get(int rowIdx, String colName) {return data.get(rowIdx, columnIndexMap.get(colName));}
	
	public List<String> getColumnNames() { return new ArrayList<String>(columnIndexMap.keySet());}
	public List<String> getRowNames() { return new ArrayList<String>(rowIndexMap.keySet());}
	public String getRowName(int i){return rows.get(i);}
	public String getColoumnName(int i){return columns.get(i);}
	public List<String> getRowDescriptions() { return new ArrayList<String>(rowDescrIndexMap.keySet());}
	public boolean hasColumn(String column) { return columnIndexMap.containsKey(column);}
	public boolean hasRow(String row) { return rowIndexMap.containsKey(row);}

	public List<Integer> getIndecesForRowDescription(String rowDescription) {
		return rowDescrIndexMap.get(rowDescription);
	}

	public void setRowDescription(int rowIdx, String rowDescription) {
		List<Integer> descrIdxs =  rowDescrIndexMap.get(rowDescription);
		if(descrIdxs == null) {
			descrIdxs = new ArrayList<Integer>();
			rowDescrIndexMap.put(rowDescription, descrIdxs);
		}
		descrIdxs.add(rowIndexMap.get(rowIdx));
	}
	
	public void setRowDescription(String row, String rowDescription) {
		setRowDescription(rowIndexMap.get(row), rowDescription);
	}
	
	public MatrixWithHeaders filterValuesLargerThanUsingColumn(int colIdx, double maxValue) {
		List<String> passingRows = new ArrayList<String>();
		Set<String> rows = rowIndexMap.keySet();
		for(String row : rows) {
			if(get(row,colIdx) < maxValue) {
				passingRows.add(row);
			}
		}
		
		return new MatrixWithHeaders(this, passingRows, getColumnNames());
	}
	
	public MatrixWithHeaders filterValuesSmallerThanUsingColumn(int colIdx, double minValue) {
		List<String> passingRows = new ArrayList<String>();
		Set<String> rows = rowIndexMap.keySet();
		for(String row : rows) {
			if(get(row,colIdx) > minValue) {
				passingRows.add(row);
			}
		}
		
		return new MatrixWithHeaders(this, passingRows, getColumnNames());
	}
	
	public MatrixWithHeaders filterValuesRangeUsingColumn(int colIdx, double minVal, double maxVal) {
		List<String> passingRows = new ArrayList<String>();
		Set<String> rows = rowIndexMap.keySet();
		for(String row : rows) {
			double cell = get(row,colIdx);
			if(cell > minVal && cell < maxVal ) {
				passingRows.add(row);
			}
		}
		
		return new MatrixWithHeaders(this, passingRows, getColumnNames());
	}

	
	/**
	 * returns true if the column has a nonzero value
	 */
	public boolean isColumnNonTrivial(String colName) {
		boolean trivial = false;
		
		for(int i = 0; i < data.getRowDimension(); i++) {
			if(get(i, colName) != 0) {
				trivial = true;
				break;
			}
		}
		return trivial;
	}
	
	/**
	 * Replaces each set of rows with same description with a row with column values the mean of the original row column values
	 * 
	 */
	public void compressMatrixByRowDescriptions() {
		List<List<Double>> newValList = new ArrayList<List<Double>>(rowDescrIndexMap.keySet().size());
		List<String> newRowNames = new ArrayList<String>(rowDescrIndexMap.keySet().size());
		for(String descr : rowDescrIndexMap.keySet()) {
			List<Integer> descrIdxs = rowDescrIndexMap.get(descr);
			List<Double> medianValues = new ArrayList<Double>(data.getColumnDimension());
			newValList.add(medianValues);
			newRowNames.add(descr);
			double [] columnValues = new double[descrIdxs.size()];
			for(int j = 0; j < data.getColumnDimension(); j++) {
				for(int i = 0 ; i < columnValues.length; i++) {
					columnValues[i] = data.get(descrIdxs.get(i), j);
				}
				medianValues.add(Statistics.median(columnValues));
			}
		}
		init(newValList, newRowNames, getColumnNames());
	}
	
	public void normalizeByMean() {
		double mean = Statistics.mean(data);
		normalizeBy(mean);
	}
	
	public void zColumnNormalize() {
		double [] tmpNormColVals = new double [rowDimension()];
		for (int j = 0; j < columnDimension(); j++) {
			double [] column = getColumn(j);
			for(int i = 0; i < rowDimension(); i++) {
				tmpNormColVals[i] = Statistics.zScore(get(i,j), column);
			}
			for(int i = 0; i < rowDimension(); i++) {
				set(i,j, tmpNormColVals[i]);
			}
		}
	}
	
	public void zRowNormalize() {
		double [] tmpNormRowVals = new double [columnDimension()];
		for (int i = 0; i < rowDimension(); i++) {
			double [] row = getRow(i);
			for(int j = 0; j < columnDimension(); j++) {
				tmpNormRowVals[j] = Statistics.zScore(get(i,j), row);
			}
			for(int j = 0; j < columnDimension(); j++) {
				set(i,j, tmpNormRowVals[j]);
			}
		}
	}
	
	
	public void normalizeBy(double constant) {
		for(int i = 0; i < data.getRowDimension(); i++){
			for(int j = 0; j < data.getColumnDimension(); j++) {
				data.set(i,j, data.get(i,j)/constant);
			}
		}
	}
	
	public void normalizeColumnsByMean() {
		for(int j = 0; j < columnDimension(); j++) {
			double [] colVals = getColumn(j);
			Double mean = Statistics.mean(colVals);
			for(int i = 0; i < rowDimension(); i++) {
				data.set(i, j, data.get(i,j)/mean);
			}
		}
	}
	
	public void normalizeColumnsByMedian() {
		for(int j = 0; j < columnDimension(); j++) {
			double [] colVals = getColumn(j);
			double median = Statistics.median(colVals);
			for(int i = 0; i < rowDimension(); i++) {
				data.set(i, j, data.get(i,j)/median);
			}
		}
	}
	
	
	public MatrixWithHeaders submatrixByRowDescriptions(Collection<String> rowDescriptions, List<String> newRowNames) {
		List<String> filteredMatrixRowDescriptions;
		List<String> rowNames;
		List<Integer> indices = new ArrayList<Integer>();
		
		int rowNum = 0;
		filteredMatrixRowDescriptions = new ArrayList<String>();
		//Get all indices found for the given row descriptions and add the row descriptions 
		for(String rowDescription : rowDescriptions) {
			List<Integer> positionsForDescription = getIndecesForRowDescription(rowDescription);
			if(positionsForDescription != null) {
				indices.addAll(positionsForDescription);
				for(int i = 0; i < positionsForDescription.size();i++) {
					filteredMatrixRowDescriptions.add(rowDescription);
				}
				rowNum += positionsForDescription.size();
			}
		}
		List<String> allRowNames = new ArrayList<String>(rowNum);
		Matrix dataSubset = new Matrix(rowNum, columnDimension());
		rowNames = new ArrayList<String>(rowNum);

		for(int i = 0; i <  indices.size(); i++) {
			int idx = indices.get(i);
			rowNames.add(newRowNames == null || newRowNames.isEmpty() ? allRowNames.get(idx) : newRowNames.get(i));
			for(int j = 0; j < columnDimension(); j++) {
				dataSubset.set(i, j, data.get(idx, j));
			}
		}
		MatrixWithHeaders filteredMWH = new MatrixWithHeaders(dataSubset, rowNames, getColumnNames(), filteredMatrixRowDescriptions);
		return filteredMWH;
	}
	
	public MatrixWithHeaders submatrixByRowNames(Collection<String> rowNames) {
		ArrayList<String> names=new ArrayList<String>();
		
		for(String rowName: rowNames){
			if(this.rowIndexMap.containsKey(rowName)){names.add(rowName);}
		}
		// return null if no names are found
		if (names.isEmpty())
			return null;
		
		MatrixWithHeaders rtrn=new MatrixWithHeaders(names, getColumnNames());
		if(this.pidToName!=null){rtrn.setPIDToName(this.pidToName);}
		
		for(String rowName: rowNames){
			if(this.rowIndexMap.containsKey(rowName)){rtrn.setRow(rowName, getRow(rowName));}
		}
		
		return rtrn;
	}
	
	public MatrixWithHeaders submatrixByRowNames(String rowName) {
		ArrayList<String> s= new ArrayList<String>();
		s.add(rowName);
		return (submatrixByRowNames(s));
	}
	public MatrixWithHeaders submatrixByRowNames(String [] rowNames) {
		List<String> rowList = new ArrayList<String>(rowNames.length);
		for(String row: rowNames) {
			rowList.add(row);
		}

		return submatrixByRowNames(rowList);
	}
	
	public MatrixWithHeaders submatrixByColumnNames(Collection<String> columnNames) {
		MatrixWithHeaders rtrn=new MatrixWithHeaders(getRowNames(), new ArrayList<String>(columnNames));
		if(this.pidToName!=null){rtrn.setPIDToName(this.pidToName);}
		
		for(String columnName: columnNames){
			//System.err.println(columnName);
			double[] array=getColumn(columnName);
			rtrn.setColumn(array, columnName);
		}
		
		return rtrn;
	}
	
	public MatrixWithHeaders submatrixByColumnIndex(Collection<Integer> columnIndex) {
		ArrayList<String> columnNames=new ArrayList();
		for(Integer index: columnIndex){columnNames.add(this.getColumnNames().get(index));}
		
		return this.submatrixByColumnNames(columnNames);
	}
	
	public MatrixWithHeaders submatrixByColumnNames(String [] columnNames) {
		List<String> columnList = new ArrayList<String>(columnNames.length);
		for(String col: columnNames) {
			columnList.add(col);
		}
		return submatrixByColumnNames(columnList);
	}
	
	
	
	public MatrixWithHeaders zScoresByColumn(List<String> controlRows) {
		Matrix zscoreMatrix = new Matrix(rowDimension(), columnDimension());
		for(int j = 0; j < columnDimension(); j++) {
			double [] columnControlValues = getColumnValues(controlRows, j);
			for(int i = 0; i < rowDimension(); i++) {
				zscoreMatrix.set(i, j, Statistics.zScore(get(i,j), columnControlValues));
			}
		}
		
		return new MatrixWithHeaders(zscoreMatrix, getRowNames(), getColumnNames());
	}
	
	public MatrixWithHeaders tScoresByRow(List<Integer> controlRows, List<Integer> sampleRows, String comparisonName, double fudgeFactor) {
		Map<String, Double> rtrn=new TreeMap();
		Matrix matrix=new Matrix(rowDimension(), 1);
		for(int j = 0; j < rowDimension(); j++) {
			double [] columnControlValues = getRowValuesByIndex(controlRows, j);
			double[] columnSampleValues=getRowValuesByIndex(sampleRows, j);
			matrix.set(j, 0, Statistics.tstat(columnControlValues, columnSampleValues, fudgeFactor));
			//rtrn.put(this.getRowNames().get(j), Statistics.tstat(columnControlValues, columnSampleValues, 0));
		}
		
		ArrayList columnNames=new ArrayList();
		columnNames.add(comparisonName);
		
		return new MatrixWithHeaders(matrix, getRowNames(), columnNames);
	}
	
	//returns a new MatrixWithHeaders that is the value fot he rank in each row
	public MatrixWithHeaders rank(){
		MatrixWithHeaders rtrn=new MatrixWithHeaders(getRowNames(), getColumnNames());
		
		//for each experiment get its corresponding values
		for(int i=0; i<columnDimension(); i++){
			double[] vals=getColumn(i);
			double[] rank=Statistics.rank(vals, true);
			rtrn.setColumn(rank, getColumnNames().get(i));
		}
		return rtrn;
	}
		
	
	public void quantileNormalizeColumns() {
		data.quantileNormalizeColumns();
	}
	
	public MatrixWithHeaders confidenceByColumn(List<String> controlRows) {
		MatrixWithHeaders confidence = zScoresByColumn(controlRows);

		for(int j = 0; j < columnDimension(); j++) {
			//Set up datastructures
			List<List<Double>> colPermutations = new ArrayList<List<Double>>(controlRows.size());
			for(int k = 0; k < controlRows.size(); k++) {
				colPermutations.add(new ArrayList<Double>());
			}
			List<Double> colObservations = new ArrayList<Double>(rowDimension());
			double [] columnControlValues = getColumnValues(controlRows, j);

			//Now fill them up
			for(int i = 0; i < rowDimension(); i++) {
				colObservations.add(confidence.get(i,j));
				double [] permutedVals = Statistics.computePermutedZScores(confidence.get(i,j), columnControlValues);
				for(int k = 0; k < controlRows.size(); k++) {
					colPermutations.get(k).add(permutedVals[k]);
				}
			}

			//sort
			Collections.sort(colObservations);
			for(List<Double> permutation : colPermutations) {
				Collections.sort(permutation);
			}
			//System.err.println("colObservations sorted: " + colObservations);
			//System.err.println("first col permutation sorted " + colPermutations.get(0));
			//Compute confidence values
			for(int i = 0; i < rowDimension(); i++) {
				double val = confidence.get(i,j);
				int k = val < 0 ? 0 : colObservations.size() - 1;
				double obsVal = 0;
				double fdr = 1;
				while(k >=0 && k < colObservations.size() && (colObservations.get(k) * obsVal >= 0)  ){// check there is no sign change {
					//System.err.println("obsVal: " + obsVal + " val " + val);
					obsVal = colObservations.get(k);
					if(Math.abs(obsVal)  <=  Math.abs(val)) {
						int observationCount = obsVal < 0 ? countValuesLessThan(obsVal, colObservations) : countValuesLargerThan(obsVal, colObservations);
						int permutedCount = 0;
						for(List<Double> permutation: colPermutations) {
							permutedCount += val < 0 ? countValuesLessThan(obsVal, permutation) : countValuesLargerThan(obsVal, permutation);
						}
						//System.err.println("position ("+i+","+j+") val " + val + " observed count " + observationCount + " permuted count " +( permutedCount/(double)controlRows.size()));
						fdr = Math.min(fdr, permutedCount/(double)(observationCount * controlRows.size()));
					}
					k = val < 0 ? k + 1 : k - 1;
				}
				confidence.set(i, j, 1 - fdr);
			}
		}
		
		return confidence;
	}
	
	//Sort matrix by values for a reference column index
	public MatrixWithHeaders sortList(int columnIndex){
		double[] vals=this.getColumn(columnIndex);
		double[] ranks=Statistics.rank(vals, false);
		
		Map<Double, ArrayList> temp=new TreeMap();
		for(int i=0; i<ranks.length; i++){
			ArrayList list=new ArrayList();
			if(temp.containsKey(ranks[i])){list=temp.get(ranks[i]);}
			list.add(i);
			temp.put(ranks[i], list);
		}
		
		List<String> rowNames= getRowNames();
		List<String> newRowOrder=new ArrayList<String>();
		
		for(Double rank:temp.keySet()){
			ArrayList<Integer> indexes=temp.get(rank);
			for(int index: indexes){newRowOrder.add(rowNames.get(index));}
		}
		
		MatrixWithHeaders rtrn=new MatrixWithHeaders(newRowOrder, this.getColumnNames());
		
		for(String name: newRowOrder){
			rtrn.setRow(name, this.getRow(name));
		}
		
		return rtrn;
	}
	
	public void append(MatrixWithHeaders other) {
		if(other == null || other.rowDimension() != rowDimension()) {
			throw new IllegalArgumentException ("To append a matrix both matrices must have same number of rows");
		}
		
		Matrix newData = new Matrix(rowDimension(), columnDimension() + other.columnDimension());
		List<String> columnNames = getColumnNames();
		List<String> rowNames    = getRowNames();
		columnNames.addAll(other.getColumnNames());
		for(int i = 0; i < rowDimension(); i++) {
			for(int j = 0; j < columnDimension(); j++) {
				newData.set(i, j, data.get(i,j));
			}
			for(int j = 0; j < other.columnDimension(); j++) {
				newData.set(i, j + columnDimension(), other.get(i,j));
			}
		}
		this.data = newData;
		initNameIndexMaps(rowNames, columnNames);
	}
	
	public void appendColumns(MatrixWithHeaders other){
		append(other);
	}
	
	public void appendRows(MatrixWithHeaders other) {
		if(other == null || other.columnDimension() != columnDimension()) {
			throw new IllegalArgumentException ("To append a matrix both matrices must have same number of columns");
		}
		
		Matrix newData = new Matrix(rowDimension() + other.rowDimension(), columnDimension());
		List<String> columnNames = getColumnNames();
		List<String> rowNames    = getRowNames();
		rowNames.addAll(other.getRowNames());
		for(int j = 0; j < columnDimension(); j++) {
			for(int i = 0; i < rowDimension(); i++) {
				newData.set(i, j, data.get(i,j));
			}
			for(int i = 0; i < other.rowDimension(); i++) {
				newData.set(i + rowDimension(), j, other.get(i,j));
			}
		}
		this.data = newData;
		initNameIndexMaps(rowNames, columnNames);
	}
	
	public void add(MatrixWithHeaders m2) {
		data.plusEquals(m2.data);
	}
	
	public void log2() {
		log2(0);
	}
	
	public void log2(double fudge) {
		for(int i = 0; i < data.getRowDimension(); i++) {
			for(int j = 0; j < data.getColumnDimension(); j++) {
				data.set(i, j, MathUtil.log2(data.get(i,j) + fudge));
			}
		}
	}
	
	public void log10() {
		log10(0);
	}
	
	public void log10(double fudge) {
		for(int i = 0; i < data.getRowDimension(); i++) {
			for(int j = 0; j < data.getColumnDimension(); j++) {
				data.set(i, j, Math.log10(data.get(i,j) + fudge));
			}
		}
	}
	
	public void pow() {
		for(int i = 0; i < data.getRowDimension(); i++) {
			for(int j = 0; j < data.getColumnDimension(); j++) {
				data.set(i, j, Math.pow(data.get(i,j),2));
			}
		}
	}
	
	public void round() {
		for(int i = 0; i < data.getRowDimension(); i++) {
			for(int j = 0; j < data.getColumnDimension(); j++) {
				data.set(i, j, Math.round(data.get(i,j)));
			}
		}
	}
	

	public MatrixWithHeaders times(MatrixWithHeaders other) {
		if(data.getColumnDimension() != other.rowDimension()) {
			throw new IllegalArgumentException ("Trying to multiply non compatible matrices  cols on left mat " + columnDimension() + " rows on right " + other.rowDimension());
		}
		
		Matrix resultData = data.times(other.data);
		MatrixWithHeaders resultMatrixWithHeaders =  new MatrixWithHeaders(resultData, getRowNames(), other.getColumnNames(), getRowDescriptions());
		System.err.println ("W row " + rowDimension() + " H col dim " + other.columnDimension() + " result dims: " + resultData.getRowDimension() + "x" + resultData.getColumnDimension() + " matrix with headers " +  resultMatrixWithHeaders.rowDimension() + "x" + resultMatrixWithHeaders.columnDimension());
		System.err.println("result column names " + resultMatrixWithHeaders.getColumnNames() + " H columns " + other.getColumnNames());
		return resultMatrixWithHeaders;
	}
	
	public void minus(MatrixWithHeaders other) {
		if(data.getColumnDimension() != other.columnDimension() || data.getRowDimension() != other.rowDimension()) {
			throw new IllegalArgumentException ("Trying to substract non compatible matrices  this has dimension " + rowDimension() + "x" + columnDimension() + "  other is " + other.rowDimension() + "x" +other.columnDimension());
		}
		
		data.minusEquals(other.data);
	}
	
	private int countValuesLargerThan(double val,List<Double> colObservations) {
		int count = colObservations.size() - 1;
		//System.err.print("val " + val + " obs " + colObservations);
		while( count >= 0 && val <= colObservations.get(count) ) {
			count--;
		}
		//System.err.println("# observed: " + (colObservations.size() -  count - 1));
		return colObservations.size() -  count - 1;
	}
	
	private int countValuesLessThan(double val,List<Double> colObservations) {
		int count = 0;
		//System.err.print("val " + val + " obs " + colObservations);
		while( count < colObservations.size() && val >= colObservations.get(count) ) {
			count++;
		}
		//System.err.println("# observed: " + count);
		return count;
	}
	
	private double [] getColumnValues(List<String> rowNames, int col) {
		double [] vals = new double[rowNames.size()];
		for(int i = 0; i < rowNames.size(); i++) {
			vals[i] = get(rowNames.get(i), col);
		}
		return vals;
	}
	
	private double [] getRowValues(List<String> columnNames, int row) {
		double [] vals = new double[columnNames.size()];
		for(int i = 0; i < columnNames.size(); i++) {
			vals[i] = get(row, columnNames.get(i));
		}
		return vals;
	}
	
	private double [] getRowValuesByIndex(List<Integer> columnNames, int row) {
		double [] vals = new double[columnNames.size()];
		for(int i = 0; i < columnNames.size(); i++) {
			vals[i] = get(row, columnNames.get(i));
		}
		return vals;
	}

	protected void initFromReader(BufferedReader br) throws IOException,
			ParseException {
		String header = br.readLine();
		if(header.startsWith("#1.2")) {
			//this is a GCT file
			initGCTFromReader(br);
		} 
		else if(header.startsWith("Code Class")){
			//this is a nanostring file
			initFromNanostring(br, header);
		}
		else {
			initFromRegularMatrix(br, header);
		}
	}
	
	public Matrix getData() { return data;}
	
	protected void initFromRegularMatrix(BufferedReader br, String header)
			throws IOException, ParseException {
		String [] columnNames = header.split("\t");
		List<String> columnNameList = new ArrayList<String> (columnNames.length);
		List<String> rowNameList = new ArrayList<String>();
		
		for(int i = 1; i < columnNames.length; i++){
			columnNameList.add(columnNames[i]);
		}
		
		String line = null;
		List<List<Double>> rawData = new ArrayList<List<Double>>();
		int lineNum = 1;
		while( (line = br.readLine()) != null) {
			String [] info = line.split("\t");
			if(info.length != columnNames.length) {
				throw new ParseException("Line " + lineNum + " has " + info.length + " columns but header had " + columnNames.length + " columns");
			}
			rowNameList.add(info[0]);
			List<Double> lineData = new ArrayList<Double>(info.length - 1);
			rawData.add(lineData);
			for(int i = 1 ; i < info.length; i++) {
				lineData.add(Double.parseDouble(info[i]));
			}
			lineNum++;
		}
		
		init(rawData, rowNameList, columnNameList);
	}
	
	protected void initGCTFromReader(BufferedReader br) throws IOException, ParseException {
		String line = br.readLine();
		if(line.startsWith("#1.2")) {
			line = br.readLine();
		}
		String [] dimensionsStr = line.split("\t");
		int expectedRowDimension = Integer.parseInt(dimensionsStr[0]);
		int expectedColumnDimension = Integer.parseInt(dimensionsStr[1]);
		
		String header = br.readLine();
		String [] columnNames = header.split("\t");
		List<String> columnNameList = new ArrayList<String> (columnNames.length - 2);
		if(columnNames.length - 2 != expectedColumnDimension) {
			System.err.println("WARNING: expected "+expectedColumnDimension+ " columns but read " +  (columnNames.length - 2) );
		}
		List<String> rowNameList = new ArrayList<String>(expectedRowDimension);
		List<String> rowDescrList = new ArrayList<String>(expectedRowDimension);
		
		for(int i = 2; i < columnNames.length; i++){
			columnNameList.add(columnNames[i]);
		}
		
		
		this.pidToName=new TreeMap<String, String>();
		List<List<Double>> rawData = new ArrayList<List<Double>>(expectedRowDimension);
		int lineNum = 0;
		while( (line = br.readLine()) != null) {
			String [] info = line.split("\t");
			if(info.length != columnNames.length) {
				throw new ParseException("Line " + lineNum + " has " + info.length + " columns but header had " + columnNames.length + " columns");
			}
			rowNameList.add(info[0].toUpperCase().intern());
			rowDescrList.add(info[1].toUpperCase().intern());
			this.pidToName.put(info[0].toUpperCase(), info[1].toUpperCase());
			List<Double> lineData = new ArrayList<Double>(info.length - 2);
			rawData.add(lineData);
			for(int i = 2 ; i < info.length; i++) {
				lineData.add(Double.parseDouble(info[i]));
			}
			if(lineNum % 1000 == 0) {System.out.println("Line  " +  lineNum + " Used Mem " + ((Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1024)  );}
			lineNum++;
		}
		if(lineNum != expectedRowDimension) {
			System.err.println("WARNING: Expected " + expectedRowDimension + " but  read " + lineNum);
		}
		init(rawData, rowNameList, columnNameList, rowDescrList);
	}
	
	
	public Map<String, Collection<Integer>> getNanostringReplicateMap(){return this.replicateMap;}
	
	protected void initFromNanostring(BufferedReader br, String header) throws IOException, ParseException {
		this.replicateMap=new TreeMap<String, Collection<Integer>>();
		
		String line = header;
		System.err.println(line);
		String [] columnNames = line.split("\t");
		List<String> columnNameList = new ArrayList<String> ();
		
		List<String> rowNameList = new ArrayList<String>();
		List<String> rowDescrList = new ArrayList<String>();
		
		for(int i = 3; i < columnNames.length; i++){
			columnNameList.add(columnNames[i]);
			Collection<Integer> indexes=new ArrayList<Integer>();
			if(replicateMap.containsKey(columnNames[i])){indexes=replicateMap.get(columnNames[i]);}
			indexes.add(i-3);
			this.replicateMap.put(columnNames[i], indexes);
		}
		
		
		this.pidToName=new TreeMap<String, String>();
		this.probeClass=new TreeMap<String, String>();
		List<List<Double>> rawData = new ArrayList<List<Double>>();
		int lineNum = 0;
		while( (line = br.readLine()) != null) {
			String [] info = line.split("\t");
			rowNameList.add(info[1]);
			rowDescrList.add(info[2]);
			this.pidToName.put(info[1], info[2]);
			this.probeClass.put(info[1], info[0]);
			
			List<Double> lineData = new ArrayList<Double>();
			rawData.add(lineData);
			for(int i = 3 ; i < info.length; i++) {
				lineData.add(Double.parseDouble(info[i]));
			}
			lineNum++;
		}
		
		init(rawData, rowNameList, columnNameList, rowDescrList);
	}
	
	public Map<String, String> getNanostringProbeClasses(){return this.probeClass;}

	private void init(List<List<Double>> rawData, List<String> rowNameList, List<String> columnNameList, List<String> rowDescriptionList ) {
		int m = rowNameList.size();
		int n = columnNameList.size();
		
		data = new Matrix(m,n);
		for(int i = 0; i < m; i++){
			for(int j = 0; j < n; j++) {
				data.set(i, j, rawData.get(i).get(j));
			}
		}
		initNameIndexMaps( rowNameList, columnNameList, rowDescriptionList);
	}
	
	private void init(List<List<Double>> rawData, List<String> rowNameList, List<String> columnNameList ) {
		init(rawData, rowNameList, columnNameList, new ArrayList<String>());
	}

	protected void initNameIndexMaps(List<String> rowNameList,List<String> columnNameList) {
		initNameIndexMaps(rowNameList, columnNameList, new ArrayList<String>());
	}
	
	protected void setData(Matrix data) {
		this.data = data;
	}
	
	protected void initRowIndexMaps(List<String> rowNameList) {
		int m = rowNameList.size();
		rowIndexMap    = new LinkedHashMap<String, Integer>(m);
		HashMap<String, Integer> occurrence = new HashMap<String, Integer>();
		occurrence = new HashMap<String, Integer>();
		for(int i = 0; i < rowNameList.size(); i ++) {
			String rowName = rowNameList.get(i);
			String uniqueRowName = getUniqueName(occurrence, rowName);
			rowIndexMap.put(uniqueRowName, i);
		}
	}
	
	protected void initColIndexMaps(List<String> columnNameList) {
		int n = columnNameList.size();
		columnIndexMap = new LinkedHashMap<String, Integer>(n);
		HashMap<String, Integer> occurrence = new HashMap<String, Integer>();
		for(int i = 0; i < columnNameList.size(); i ++) {
			String colName = columnNameList.get(i);
			String uniqueColName = getUniqueName(occurrence, colName);
			columnIndexMap.put(uniqueColName, i);
		}
	}
	
	protected void initRowDescrIndexMaps(List<String> rowDescriptionList) {
		rowDescrIndexMap = new LinkedHashMap<String, List<Integer>>();
	
		for(int i = 0; i < rowDescriptionList.size(); i ++) {
			String rowDescr = rowDescriptionList.get(i).toUpperCase();
			List<Integer> rowDescrIndeces = rowDescrIndexMap.get(rowDescr);
			if(rowDescrIndeces == null) {
				rowDescrIndeces = new ArrayList<Integer>();
				rowDescrIndexMap.put(rowDescr, rowDescrIndeces);
			}
			rowDescrIndeces.add(i);
		}
	}
	
	private void initNameIndexMaps(List<String> rowNameList,List<String> columnNameList, List<String> rowDescriptionList) {
		initRowIndexMaps(rowNameList);
		initColIndexMaps(columnNameList);
		initRowDescrIndexMaps(rowDescriptionList);
		rows = rowNameList;
		columns = columnNameList;
		
	}
	private String getUniqueName(HashMap<String, Integer> occurrence, String key) {
		Integer colNameOccurrence = occurrence.get(key) ;
		String uniqueKey = key;
		if(colNameOccurrence == null){
			occurrence.put(key,1);
		} else {
			occurrence.put(key,colNameOccurrence++);
			uniqueKey = getUniqueName(occurrence, key + "_" + colNameOccurrence);
		}
		
		return uniqueKey;
	}

	public int getNumberColumns(){return this.columnDimension();}
	public int getNumberRows(){return this.rowDimension();}
	
	public Map<String, double[]> toMap(){
		Map<String, double[]> rtrn=new TreeMap();
		
		for(String geneName: this.getRowNames()){
			double[] array=this.getRow(geneName);
			//System.err.println(geneName+" "+Statistics.average(array)+" "+array.length);
			rtrn.put(geneName, array);
		}
		
		return rtrn;
	}



	public void writeGMT(String save, double fold) throws IOException {
		FileWriter writer=new FileWriter(save);
		System.err.println("here");
		Collection<String> columnNames=getColumnNames();
		
		for(String name: columnNames){
			List<String> rowNames=getRowNames();
			System.err.println(name);
			double[] vals=getColumn(name);
			writer.write(name+"_UP"+"\t"+name+"_UP");
			for(int i=0; i<vals.length; i++){
				if(vals[i]>fold){writer.write("\t"+rowNames.get(i));}
			}
			writer.write("\n");
			
			writer.write(name+"_DOWN"+"\t"+name+"_DOWN");
			for(int i=0; i<vals.length; i++){
				if(vals[i]<(-fold)){writer.write("\t"+rowNames.get(i));}
			}
			writer.write("\n");
		}
		
		writer.close();
	}

	//Moran's version
	/*public MatrixWithHeaders medianNorm(){
		MatrixWithHeaders rtrn=new MatrixWithHeaders(this.rows, this.columns);
		
		for(int i=0; i<data.getColumnDimension(); i++){
			double median=Statistics.median(data.getColumn(i));
			double[] vals=medianNormColumn(data.getColumn(i), median);
			rtrn.setColumn(vals, i);
			System.err.println(this.columns.get(i)+" "+median+" "+this.getColumn(i).length);
		}
		
		return rtrn;
	}*/
	
	private double[] shiftArrayByConstant(double[] column, double median) {
		double[] rtrn=new double[column.length];
		for(int i=0; i<column.length; i++){
			rtrn[i]=column[i]-median;
		}
		return rtrn;
	}

	
	public MatrixWithHeaders quantileNormalize(){
		MatrixWithHeaders copy=copy();
		copy.setPIDToName(this.pidToName);
		copy.quantileNormalizeColumns();
		return copy;
	}
	
	public MatrixWithHeaders copy(){
		MatrixWithHeaders rtrn=new MatrixWithHeaders(this.rows, this.columns);
		
		for(int i=0; i<this.columnDimension(); i++){
			rtrn.setColumn(this.getColumn(i), i);
		}
		
		return rtrn;
	}

	
	public MatrixWithHeaders medianCenterColumns() {
		MatrixWithHeaders rtrn=new MatrixWithHeaders(this.rows, this.columns);
		rtrn.setPIDToName(this.pidToName);
		//double maxMAD=Statistics.maxStdev(data);
		//double maxMAD=Statistics.geometricMeanMAD(data);
		
		for(int i=0; i<data.getColumnDimension(); i++){
			double median=Statistics.median(data.getColumn(i));
			//System.err.println("Median :"+median);
			//double sampleMAD=Statistics.stdev(data.getColumn(i));
			//double sampleMAD=Statistics.MAD(data.getColumn(i));
			double[] vals=shiftArrayByConstant(data.getColumn(i), median);
			rtrn.setColumn(vals, i);
			//System.err.println(this.columns.get(i)+" "+median+" "+sampleMAD+" "+maxMAD+" "+(maxMAD/sampleMAD));
		}
		return rtrn;
	}
	
	public MatrixWithHeaders medianCenterRows() {
		MatrixWithHeaders rtrn=new MatrixWithHeaders(this.rows, this.columns);
		rtrn.setPIDToName(this.pidToName);
		//double maxMAD=Statistics.maxStdev(data);
		//double maxMAD=Statistics.geometricMeanMAD(data);
		
		for(int i=0; i<data.getRowDimension(); i++){
			double median=Statistics.median(data.getRow(i));
			//System.err.println("Median :"+median);
			//double sampleMAD=Statistics.stdev(data.getColumn(i));
			//double sampleMAD=Statistics.MAD(data.getColumn(i));
			double[] vals=shiftArrayByConstant(data.getRow(i), median);
			rtrn.setColumn(vals, i);
			//System.err.println(this.columns.get(i)+" "+median+" "+sampleMAD+" "+maxMAD+" "+(maxMAD/sampleMAD));
		}
		return rtrn;
	}

	/**
	 * Mean centers the columns of the current matrix rather than creating a new one.
	 *
	 */
	public void medianCenterColumnsThis() {
		for(int i=0; i<data.getColumnDimension(); i++){
			double median=Statistics.median(data.getColumn(i));
			for(int j =0; j < data.getRowDimension(); j++ ) {
				data.set(j, i, data.get(j,i) - median);
			}
		}
	}
	
	/**
	 * Mean centers the rows of the current matrix rather than creating a new one.
	 * 
	 */
	public void medianCenterRowsThis() {

		for(int i=0; i<data.getRowDimension(); i++){
			double median=Statistics.median(data.getRow(i));
			for(int j =0; j < data.getColumnDimension(); j++ ) {
				data.set(i, j, data.get(i,j) - median);
			}
		}
	}
	
	/**
	 * Substitutes each value in a column with the ratio of the value and the mean of column value
	 * Only applicable to positive value matrices.
	 *
	 */
	public void meanNormalizeColumns() throws IllegalStateException {
		for(int j=0; j<data.getColumnDimension(); j++){
			double mean=Statistics.mean(data.getColumn(j));
			if(mean > 0) {
				for(int i =0; i < data.getRowDimension(); i++ ) {
					if(data.get(i,j) < 0) {
						throw new IllegalStateException("mean normalization is only applicable to positive valued matrices, value at ("+i+","+j+") is " + data.get(i,j));
					}
					data.set(i, j, data.get(i,j)/mean);
				}
			}
		}
	}
	
	/**
	 * Substitutes each value in a column with the ratio of the value and the mean of column value
	 * Only applicable to positive value matrices.
	 *
	 */
	public void medianNormalizeColumns() throws IllegalStateException {
		for(int j=0; j<data.getColumnDimension(); j++){
			double median=Statistics.median(data.getColumn(j));
			if(median > 0) {
				for(int i =0; i < data.getRowDimension(); i++ ) {
					if(data.get(i,j) < 0) {
						throw new IllegalStateException("mean normalization is only applicable to positive valued matrices, value at ("+i+","+j+") is " + data.get(i,j));
					}
					data.set(i, j, data.get(i,j)/median);
				}
			}
		}
	}
	
	/**
	 * Substitutes each value in a row with the ratio of the value and the mean of row value
	 * Only applicable to positive value matrices.
	 */
	public void meanNormalizeRows() throws IllegalStateException {
		for(int i=0; i<data.getRowDimension(); i++){
			double mean=Statistics.mean(data.getRow(i));
			if(mean > 0) {
				for(int j =0; j < data.getColumnDimension(); j++ ) {
					if(data.get(i,j) < 0) {
						throw new IllegalStateException("mean normalization is only applicable to positive valued matrices, value at ("+i+","+j+") is " + data.get(i,j));
					}
					data.set(i, j, data.get(i,j) / mean);
				}
			}
		}
	}
	
	/**
	 * Substitutes each value in a row with the ratio of the value and the mean of row value
	 * Only applicable to positive value matrices.
	 */
	public void medianNormalizeRows() throws IllegalStateException {
		for(int i=0; i<data.getRowDimension(); i++){
			double median=Statistics.median(data.getRow(i));
			if(median > 0) {
				for(int j =0; j < data.getColumnDimension(); j++ ) {
					if(data.get(i,j) < 0) {
						throw new IllegalStateException("mean normalization is only applicable to positive valued matrices, value at ("+i+","+j+") is " + data.get(i,j));
					}
					data.set(i, j, data.get(i,j) / median);
				}
			}
		}
	}
	
	//Performs normalization on each column so that 1) medians are all the same 2) the variance across columns is the same
	//From Yang et al. 2002 NAR
	//Compute scale factor ai for each array i and then multiply each value by 1/ai
	//ai is estimated as the sample MAD divided by the geometric mean of all MADs
	public MatrixWithHeaders scaleNorm() {
		MatrixWithHeaders rtrn=new MatrixWithHeaders(this.rows, this.columns);
		rtrn.setPIDToName(this.pidToName);
		//double maxMAD=Statistics.maxStdev(data);
		double maxMAD=Statistics.geometricMeanMAD(data);
		
		for(int i=0; i<data.getColumnDimension(); i++){
			double median=Statistics.median(data.getColumn(i));
			//double sampleMAD=Statistics.stdev(data.getColumn(i));
			double sampleMAD=Statistics.MAD(data.getColumn(i));
			double[] vals=scaleNormColumn(data.getColumn(i), sampleMAD, maxMAD, median);
			rtrn.setColumn(vals, i);
			System.err.println(this.columns.get(i)+" "+median+" "+sampleMAD+" "+maxMAD+" "+(maxMAD/sampleMAD));
		}
		return rtrn;
	}

	private double[] scaleNormColumn(double[] vals, double sampleMAD, double maxMAD, double median){
		double scaleFactor=maxMAD/sampleMAD;
		
		double[] rtrn=new double[vals.length];
		
		for(int i=0; i<rtrn.length; i++){
			double MVal=vals[i]-median;
			double newMVal=scaleFactor*MVal;
			rtrn[i]=newMVal;
		}
		return rtrn;
	}
	
	/**
	 * Normalizes each column so it adds up to 1
	 * @throws IllegalStateException if a column contains negative numbers
	 */
	public void columnDensityNormalization() throws IllegalStateException{
		for(int j = 0 ; j < data.getColumnDimension(); j++) {
			double [] normalizedCol = normalizeToDensity(data.getColumn(j));
			for(int i = 0 ; i < data.getRowDimension(); i++) {
				data.set(i, j, normalizedCol[i]);
			}
		}
		
		
	}
	
	/**
	 * Normalizes each row so it adds up to 1
	 * @throws IllegalStateException if a column contains negative numbers
	 */
	public void rowDensityNormalization() throws IllegalStateException{
		for(int i = 0 ; i < data.getRowDimension(); i++) {
			double [] normalizedRow = normalizeToDensity(data.getRow(i));
			for(int j = 0 ; j < data.getColumnDimension(); j++) {
				data.set(i, j, normalizedRow[j]);
			}
		}	
	}	
	
	private double [] normalizeToDensity(double [] d) throws IllegalStateException {
		double [] rtrn = new double[d.length];
		double sum = Statistics.sum(d);
		if(sum  != 0  ) {
			for(int i = 0; i < d.length; i++) {
				if(d[i]>=0) {
					rtrn[i] = d[i]/sum;
				} else {
					throw new IllegalStateException("Negative number found in array. Only positive arrays can be normalized to density");
				}
			}
		}
		
		return rtrn;
	}
	
	public void setPIDToName(Map<String, String> rowDescriptions){
		if(rowDescriptions==null){
			System.err.println("Row description is null");
			return;
		}
		Map<String, String> allUpper=new TreeMap<String, String>();
		for(String pid: rowDescriptions.keySet()){
			String name=rowDescriptions.get(pid);
			allUpper.put(pid.toUpperCase(), name.toUpperCase());

			List<Integer> descriptionRowIdxs = rowDescrIndexMap.get(name);
			if(descriptionRowIdxs == null) {
				descriptionRowIdxs = new ArrayList<Integer>();
				rowDescrIndexMap.put(name, descriptionRowIdxs);
			}
			descriptionRowIdxs.add(rowIndexMap.get(pid.toUpperCase()));
		}
		
		this.pidToName=allUpper;
	}


	public void writeCLS(String string, Collection<String> subset, Collection<String> negatives) throws IOException {
		FileWriter writer=new FileWriter(string);
		writer.write((subset.size()+negatives.size())+"\t2\t1\n");
		writer.write("# sample\tcontrol\n");
		for(String sample: subset){writer.write("0\t");}
		writer.write("\n");
		for(String sample: negatives){writer.write("1\t");}
		writer.write("\n");
		writer.write("control\n");
		writer.close();
	}
	
	
	public void writeGCTAndCLS(String save, Map<String, ExpressionExperimentInfo> experimentInfo) throws IOException{
		FileWriter writerCLS=new FileWriter(save+".cls");
		FileWriter writerGCT=new FileWriter(save+".gct");
		
		Collection<String> passedQC=new TreeSet<String>();
		Map<String, Collection<String>> groups=new TreeMap<String, Collection<String>>();
		
		Collection<String> samples1=new TreeSet<String>();
		
		for(String sample: getColumnNames()){
			if(experimentInfo.containsKey(sample)){samples1.add(sample);}
		}
		
		
		for(String barcode: getColumnNames()){
			if(samples1.contains(barcode)){
			ExpressionExperimentInfo info=experimentInfo.get(barcode);
			if(info.passedQC()){
				passedQC.add(barcode);
				if(info.isControl()){
					Collection<String> samples=new TreeSet<String>();
					if(groups.containsKey("Controls")){samples=groups.get("Controls");}
					samples.add(barcode);
					groups.put("Controls", samples);
				}
				else{
					Collection<String> samples=new TreeSet<String>();
					if(groups.containsKey(info.getExperimentName())){samples=groups.get(info.getExperimentName());}
					samples.add(barcode);
					groups.put(info.getExperimentName(), samples);
				}
			}	
			}
		}
		
		writerGCT.write("#1.2\n"+rowDimension()+"\t"+passedQC.size()+"\n");
		writerCLS.write(passedQC.size()+"\t"+groups.keySet().size()+"\t1\n");
		writerCLS.write("#");
		for(String group: groups.keySet()){
			writerCLS.write("\t"+group);
		}
		writerCLS.write("\n");
		
		writerGCT.write("PID\tName");
		for(String group: groups.keySet()){
			Collection<String> samples=groups.get(group);
			for(String barcode: samples){
				ExpressionExperimentInfo info=experimentInfo.get(barcode);
				if(info.passedQC()){
					if(info.isControl()){writerCLS.write("Controls");}
					else{writerCLS.write(info.getExperimentName());}
					writerCLS.write("\t");
					writerGCT.write("\t"+barcode);
				}
			}
		}
		writerCLS.write("\n");
		writerGCT.write("\n");
		
		
		for(String gene: getRowNames()){
			writerGCT.write(gene+"\t"+gene);
			for(String group: groups.keySet()){
				Collection<String> samples=groups.get(group);
				for(String barcode: samples){
					ExpressionExperimentInfo info=experimentInfo.get(barcode);
					if(info.passedQC()){
						writerGCT.write("\t"+get(gene, barcode));
					}
				}	
			}
			writerGCT.write("\n");
		}
		
		writerCLS.close();
		writerGCT.close();
	}

   //randomly permutes each column of the matrix (independently)
	public void randomPermuteColumns(){
		for (int i=0;  i<this.getNumberColumns(); i++){
			double[] a= this.getColumn(i);
			int [] randperm=Statistics.randomPermutation(this.getNumberRows());
			for (int j=0; j<this.getNumberRows(); j++){this.set(j,i, a[randperm[j]-1]);}
		}
	
	}


	//Only works is the matrix was defines to add rows
	public void addRow(String row, String description, double[] vals) {
		
		int index=this.getNumberRows();
		rowIndexMap.put(row, index);
		
		List<Integer> rowDescrIndeces = rowDescrIndexMap.get(description);
		if(rowDescrIndeces == null) {
			rowDescrIndeces = new ArrayList<Integer>();
			rowDescrIndexMap.put(description, rowDescrIndeces);
		}
		rowDescrIndeces.add(index);
		
		for(int i=0; i<vals.length; i++){
			set(row, i, vals[i]);
		}
		
	}

	public void addColumn(String column) {
		
		Matrix newData = new Matrix(rowDimension(), columnDimension()+1);
		for(int i=0; i< rowDimension(); i++){
			for(int j = 0; j < columnDimension(); j++) {
				newData.set(i, j, get(i,j));
			}
		}
		
		setData(newData);	
		int index=columns.size();
		columnIndexMap.put(column, index);
		columns.add(column);
	}


	public MatrixWithHeaders filterInvariantGenes(MatrixWithHeaders data, double fold) {
		List<String> rows=new ArrayList<String>();
		
		for(String row: data.getRowNames()){
			double[] vals=data.getRow(row);
			double foldChange=Statistics.fold(vals, true);
			if(foldChange>fold){rows.add(row);}
		}
		
		return data.submatrixByRowNames(rows);
	}


	public void writeBox(String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		for(String rowName : rowIndexMap.keySet()) {
			for(String colName : columnIndexMap.keySet()) {
				writer.write(String.valueOf(get(rowName, colName))+"\t");
			}
			writer.write("\n");
		}
		
		writer.close();
	}


	public MatrixWithHeaders excludeByRowNames(Collection<String> flaggedGenes) {
		Collection<String> include=new TreeSet();
		
		
		for(String row: this.getRowNames()){
			if(!flaggedGenes.contains(row)){include.add(row);}
		}
		
		return this.submatrixByRowNames(include);
	}


	public double[] getRow(int i) {
		return this.data.getRow(i);
	}
	
	/**
	 * Resets all values below floor to floor.
	 * @param floor
	 */
	public void floor(double floor) {
		for(int i = 0; i < data.getRowDimension(); i++) {
			for(int j = 0; j < data.getColumnDimension(); j++) {
				if(get(i,j) < floor) {
					set(i,j, floor);
				}
			}
		}
	}


	public void writeRowNamesToFile(String outFile) throws IOException {
		
		FileWriter writer=new FileWriter(outFile);
		for (int i=0; i< this.getNumberRows(); i++ )
			{writer.write(this.getRowName(i)+"\n");}
		writer.close();
		
	}
    
 public void writeColoumnNamesToFile(String outFile) throws IOException {
		
		FileWriter writer=new FileWriter(outFile);
		for (int i=0; i< this.getNumberColumns(); i++ )
			{writer.write(this.getColoumnName(i)+"\n");}
		writer.close();
		
	}


public double[] getMedianOverAllRows() {
	
	double [] res =new double[this.columnDimension()];
	
	for (int c=0; c <this.columnDimension(); c++){
		ArrayList<Double> arr= new ArrayList<Double>();
		for(int r=0; r<this.rowDimension();r++){
		 arr.add(this.get(r,c));
		}
		res[c]= Statistics.median(arr);
	}
	
	return res;
}


public double[] getMeanOverAllRows() {
	
	double [] res =new double[this.columnDimension()];
	
	for (int c=0; c <this.columnDimension(); c++){
		ArrayList<Double> arr= new ArrayList<Double>();
		for(int r=0; r<this.rowDimension();r++){
		 arr.add(this.get(r,c));
		}
		res[c]= Statistics.mean(arr);
	}
	
	return res;
	
}

public double[] getRow(String gene, Collection<String> group1) {
	double[] rtrn=new double[group1.size()];
	
	int i=0;
	for(String column: group1){
		rtrn[i]=get(gene, column);
		i++;
	}
	
	return rtrn;
}


public double[] getValues(String gene, Collection<String> controls) {
	double[] rtrn=new double[controls.size()];
	
	int i=0;
	for(String control: controls){
		rtrn[i++]=get(gene, control);
	}
	
	return rtrn;
}


public boolean hasNanostringProbeClasses() {
	if(this.getNanostringProbeClasses()==null || this.getNanostringProbeClasses().isEmpty()){return false;}
	return true;
}

/**
 * @author skadri
 * Multiplies each column with a separate constant. That is, given a vector (dimension same as number of columns) multiplies all rows of column with same constant 
 * @return
 */
public MatrixWithHeaders multiplyColumnsWithConstants(double[] constants){
	
	MatrixWithHeaders resultMat = new MatrixWithHeaders(this.getRowNames(),this.getColumnNames());
	if(resultMat.columnDimension() != constants.length) {
		throw new IllegalArgumentException ("Trying to multiply non compatible matrix and vector. Columns on matrix " + resultMat.columnDimension() + " Vector Dimensions on right " + constants.length);
	}
	else{

		for(int j=0;j<resultMat.columnDimension();j++){
			for(int i=0;i<resultMat.rowDimension();i++){
				resultMat.set(i, j, (this.get(i,j)*constants[j]));
			}
		}
	}
	return resultMat;
}


/**
 * @author skadri
 * Multiplies each column with a separate constant. That is, given a vector (dimension same as number of columns) multiplies all rows of column with same constant 
 * @return
 */
public MatrixWithHeaders multiplyColumnsWithConstants(Map<String,Double> constants){
	
	MatrixWithHeaders resultMat = new MatrixWithHeaders(this.getRowNames(),this.getColumnNames());
	if(resultMat.columnDimension() != constants.keySet().size()) {
		throw new IllegalArgumentException ("Trying to multiply non compatible matrix and vector. Columns on matrix " + resultMat.columnDimension() + " Vector Dimensions on right " + constants.size());
	}
	else{

		for(int j=0;j<resultMat.columnDimension();j++){
			double constant = constants.get(resultMat.getColoumnName(j));
			for(int i=0;i<resultMat.rowDimension();i++){
				double value = constant*this.get(i, j);
				//System.out.println(constants.get(resultMat.getColoumnName(j))+" * "+ this.get(i, j));
				resultMat.set(i, j, value);
			}
		}
	}
	return resultMat;
}







	
	


	
	
	
}
