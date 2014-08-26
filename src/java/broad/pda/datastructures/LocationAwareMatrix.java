package broad.pda.datastructures;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import Jama.Matrix;
import broad.core.annotation.AnnotationHandler;
import broad.core.annotation.BasicGenomicAnnotation;
import broad.core.annotation.BasicLightweightAnnotation;
import broad.core.annotation.GenomicAnnotation;
import broad.core.annotation.LightweightGenomicAnnotation;
import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.MatrixWithHeaders;
import broad.core.datastructures.IntervalTree.Node;
import broad.core.error.ParseException;

public class LocationAwareMatrix extends MatrixWithHeaders implements AnnotationHandler {
	private LinkedHashMap<String,IntervalTree<String>> chrRowLocationTreeMap;
	private Map<String, LightweightGenomicAnnotation> rowLocationMap;

	public LocationAwareMatrix() {
		super();
		rowLocationMap = new LinkedHashMap<String, LightweightGenomicAnnotation>();
		chrRowLocationTreeMap = new LinkedHashMap<String, IntervalTree<String>>();
	}
	
	public LocationAwareMatrix (List<String> rows, List<String> columns, Map<String, LightweightGenomicAnnotation> rowLocations ) {
		super(rows, columns);
		this.rowLocationMap = rowLocations;
		chrRowLocationTreeMap = new LinkedHashMap<String, IntervalTree<String>>();
		for(String row : rowLocations.keySet()) {
			LightweightGenomicAnnotation location = rowLocations.get(row);
			IntervalTree<String> chrTree = chrRowLocationTreeMap.get(location.getChromosome());
			if(chrTree == null) {
				chrTree = new IntervalTree<String>();
				chrRowLocationTreeMap.put(location.getChromosome(), chrTree);
			}
			chrTree.put(location.getStart(),location.getEnd(), row);
		}
	}
	
	public LocationAwareMatrix(BufferedReader br, String chr) throws IOException,
			ParseException {
		String line = br.readLine();
		if(!line.startsWith("chromosome\t")) {
			throw new ParseException("LocationAwareMatrices must be instantiated from IGV formatted files");
		}
		
		initFromIGVMatrix(br, line, chr);
	}
	
	public LocationAwareMatrix(List<? extends LightweightGenomicAnnotation> probeData, List<String> samples) {
		chrRowLocationTreeMap = new LinkedHashMap<String, IntervalTree<String>>();
		rowLocationMap        = new LinkedHashMap<String, LightweightGenomicAnnotation>();
		List<String> rows = new ArrayList<String>(probeData.size());
		for(LightweightGenomicAnnotation probe : probeData) {
			rows.add(probe.toUCSC());
			IntervalTree<String> chrTree = chrRowLocationTreeMap.get(probe.getChromosome());
			if(chrTree == null) {
				chrTree  = new IntervalTree<String>();
				chrRowLocationTreeMap.put(probe.getChromosome(), chrTree);
			}
			chrTree.put(probe.getStart(),probe.getEnd(), probe.toUCSC());
			rowLocationMap.put(probe.toUCSC(), probe);
		}
		Matrix data = new Matrix(probeData.size(), samples.size());
		initNameIndexMaps(rows, samples);
		setData(data);
		
	}
	
	public LocationAwareMatrix(MatrixWithHeaders mwh) {
		chrRowLocationTreeMap = new LinkedHashMap<String, IntervalTree<String>>();
		rowLocationMap        = new LinkedHashMap<String, LightweightGenomicAnnotation>();
		List<String> rows = mwh.getRowNames();
		
		for(String row : rows) {
			LightweightGenomicAnnotation lga = BasicLightweightAnnotation.createFromUCSC(row);
			IntervalTree<String> chrTree = chrRowLocationTreeMap.get(lga.getChromosome());
			if(chrTree == null) {
				chrTree  = new IntervalTree<String>();
				chrRowLocationTreeMap.put(lga.getChromosome(), chrTree);
			}
			chrTree.put(lga.getStart(),lga.getEnd(), lga.toUCSC());
			rowLocationMap.put(lga.toUCSC(), lga);
		}
		Matrix data = mwh.getData();
		initNameIndexMaps(rows, mwh.getColumnNames());
		setData(data);
		
	}
	
	protected void initFromIGVMatrix(BufferedReader br, String header, String chromosomeToLoad)
	throws IOException, ParseException {
		String [] columnNames = header.split("\t");
		List<String> columnNameList = new ArrayList<String> (columnNames.length);
		List<String> rowNameList = new ArrayList<String>();
		
		for(int i = 4; i < columnNames.length; i++){
			columnNameList.add(columnNames[i]);
		}
		
		chrRowLocationTreeMap = new LinkedHashMap<String, IntervalTree<String>>();
		rowLocationMap        = new LinkedHashMap<String, LightweightGenomicAnnotation>();
		String line = null;
		List<List<Double>> rawData = new ArrayList<List<Double>>();
		int lineNum = 1;
		while( (line = br.readLine()) != null) {
			String [] info = line.split("\t");
			if(info.length != columnNames.length) {
				throw new ParseException("Line " + lineNum + " has " + info.length + " columns but header had " + columnNames.length + " columns");
			}
			String chr = info[0];//.replace("chr","");
			if(chromosomeToLoad == null || chromosomeToLoad.length()==0 || chromosomeToLoad.equals(chr)) {//Replaced to be compatible with 1.5
				int start = Integer.parseInt(info[1]);
				int end   = Integer.parseInt(info[2]);
				rowNameList.add(info[3]);
				List<Double> lineData = new ArrayList<Double>(info.length - 4);
				rawData.add(lineData);
				for(int i = 4 ; i < info.length; i++) {
					lineData.add(Double.parseDouble(info[i]));
				}
				IntervalTree<String> chrTree = chrRowLocationTreeMap.get(chr);
				if(chrTree == null) {
					chrTree  = new IntervalTree<String>();
					chrRowLocationTreeMap.put(chr, chrTree);
				}
				chrTree.put(start,end, info[3]);
				rowLocationMap.put(info[3], new BasicGenomicAnnotation(info[3],chr, start, end));
				lineNum++;
			}
		}
		super.initNameIndexMaps(rowNameList, columnNameList);
		Matrix data = new Matrix(rawData.size(), columnNameList.size());
		setData(data);
		for(int i =0 ; i < rawData.size(); i++) {
			List<Double> rowData = rawData.get(i);
			for(int j = 0; j < rowData.size(); j++) {
				data.set(i, j, rowData.get(j));
			}
		}
		
	}
	
	public Set<String> getChromosomes() {
		return chrRowLocationTreeMap.keySet();
	}
	
	public LightweightGenomicAnnotation getRowPosition(String row) { return rowLocationMap.get(row);}
	
	public void setColumns(List<String> columns) {
		initColIndexMaps(columns);
	}
	
	public void setDataDimensions(int rowDim, int columnDim) {
		setData(new Matrix(rowDim, columnDim));
	}
	
	public Collection<LightweightGenomicAnnotation> getRowAnnotations() {
		return rowLocationMap.values();
	}
	
	public void write(BufferedWriter bw) throws IOException {
		bw.write("chromosome\tstart\tend\tname");
		List<String> columns = getColumnNames();
		for(String columnName : columns) {
			bw.write("\t");
			bw.write(columnName);
		}
		bw.newLine();
		bw.flush();
		
		for(String chr : chrRowLocationTreeMap.keySet()) {
			Iterator<Node<String>> chrTreeNodeIt = chrRowLocationTreeMap.get(chr).iterator();
			while(chrTreeNodeIt.hasNext()) {
				Node<String> node = chrTreeNodeIt.next();
				String rowName = node.getValue();
				bw.write(chr);
				bw.write("\t" + node.getStart());
				bw.write("\t" + node.getEnd());
				bw.write("\t");
				bw.write(rowName);
				for(String colName : columns) {
					bw.write("\t");
					bw.write(String.valueOf(get(rowName, colName)));
				}
				bw.newLine();
			}
		}
	}

	public void annotation(GenomicAnnotation annotation) {
		IntervalTree<String> chrTree = chrRowLocationTreeMap.get(annotation.getChromosome());
		if(chrTree == null) {
			chrTree  = new IntervalTree<String>();
			chrRowLocationTreeMap.put(annotation.getChromosome(), chrTree);
		}
		chrTree.put(annotation.getStart(),annotation.getEnd(), annotation.getLocationString());	
	}
	
	public void eof() {
		List<String> rows = new ArrayList<String>();
		for(String  chr : chrRowLocationTreeMap.keySet()) {
			IntervalTree<String> chrTree = chrRowLocationTreeMap.get(chr);
			Iterator<Node<String>> nodeIt = chrTree.iterator();
			while(nodeIt.hasNext()) {
				rows.add(nodeIt.next().getValue());
			}
		}
		super.initRowIndexMaps(rows);
	}
	
	public void begin(){
		chrRowLocationTreeMap = new LinkedHashMap<String, IntervalTree<String>>();	
	}
	
	public IntervalTree<String> getRowTree(String chromosome) {return chrRowLocationTreeMap.get(chromosome);}
	
	public LocationAwareMatrix submatrixByColumns(Collection<String> columnNames) {
		MatrixWithHeaders tmp = super.submatrixByColumnNames(columnNames);
		LocationAwareMatrix rtn = new LocationAwareMatrix(tmp);
		return rtn;
	}
	
	public LocationAwareMatrix submatrixByColumns(String [] columnNames) {
		MatrixWithHeaders tmp = super.submatrixByColumnNames(columnNames);
		LocationAwareMatrix rtn = new LocationAwareMatrix(tmp);
		return rtn;
	}
	
	public List<LightweightGenomicAnnotation> getOverlappers(LightweightGenomicAnnotation annotation) {
		List<LightweightGenomicAnnotation> overlappers = new ArrayList<LightweightGenomicAnnotation>();
		//System.err.println("Trying to get overlappers for " + annotation.toUCSC() );
		if(chrRowLocationTreeMap.containsKey(annotation.getChromosome())) {
			//System.err.print("There is a tree for the annotation chr " + annotation.getChromosome());
			IntervalTree<String> tree = chrRowLocationTreeMap.get(annotation.getChromosome());
			Iterator<Node<String>> overlapperIt = tree.overlappers(annotation.getStart(), annotation.getEnd());
			//System.err.println(" Are there overlappers? " + overlapperIt.hasNext());
			while(overlapperIt.hasNext()) {
				String overlapperRow = overlapperIt.next().getValue();
				overlappers.add(rowLocationMap.get(overlapperRow));
			}
		}
		
		return overlappers;
	}
	

	public void browserLine(String line) {
		// Ignore
		
	}

	public void track(String line) {
		// Not clear what to do with different tracks, for now, load them all.
		
	}

	public LocationAwareMatrix removeOverlappingRows(LocationAwareMatrix other) {
		List<String> rowsToKeep = new ArrayList<String>();
		LinkedHashMap<String, LightweightGenomicAnnotation> trimmedRowLocationMap = new LinkedHashMap<String, LightweightGenomicAnnotation>();
		
		for(String row : getRowNames()) {
			if(other.getOverlappers(rowLocationMap.get(row)).isEmpty()) {
				rowsToKeep.add(row);
				trimmedRowLocationMap.put(row, rowLocationMap.get(row));
				//System.err.println("going to keep " + row);
			} 
		}
		
		MatrixWithHeaders trimmedMatrix = super.submatrixByRowNames(rowsToKeep);
		LocationAwareMatrix result = new LocationAwareMatrix(trimmedMatrix.getRowNames(), trimmedMatrix.getColumnNames(), trimmedRowLocationMap);
		for(int i = 0; i < trimmedMatrix.rowDimension();i++) {
			for(int j = 0; j < trimmedMatrix.columnDimension(); j++) {
				result.set(i, j, trimmedMatrix.get(i, j));
			}
		}
		
		return result;
	}

	public void appendRows(LocationAwareMatrix other) {
		super.appendRows(other);
		
		Iterator<String> chrIt = other.chrRowLocationTreeMap.keySet().iterator();
		while(chrIt.hasNext()) {
			String chr = chrIt.next();
			Iterator<String> otherChrAnnotationIt = other.chrRowLocationTreeMap.get(chr).valueIterator();
			while(otherChrAnnotationIt.hasNext()){
				String annotationName = otherChrAnnotationIt.next();
				LightweightGenomicAnnotation regionAnnotation = other.rowLocationMap.get(annotationName);
				
				IntervalTree<String> thisChrAnnotationTree = chrRowLocationTreeMap.get(chr);
				if(thisChrAnnotationTree == null) {
					thisChrAnnotationTree = new IntervalTree<String>();
					chrRowLocationTreeMap.put(chr, thisChrAnnotationTree);
				}
				thisChrAnnotationTree.put(regionAnnotation.getStart(), regionAnnotation.getEnd(), annotationName);
				rowLocationMap.put(annotationName, regionAnnotation);
			}
			
		}

	}

	public Map<String, LightweightGenomicAnnotation> getRowLocationMap() {return rowLocationMap;}

	public LightweightGenomicAnnotation getClosest(LightweightGenomicAnnotation annotation) {
		IntervalTree<String> chrAnnotationTree = chrRowLocationTreeMap.get(annotation.getChromosome());
		LightweightGenomicAnnotation closest = null;
		if(chrAnnotationTree != null) {
			Iterator<Node<String>> overlapperIt = chrAnnotationTree.overlappers(annotation.getStart(), annotation.getEnd());
			if(overlapperIt.hasNext()) {
				closest = rowLocationMap.get(overlapperIt.next().getValue());
				
			} else {
				Node<String> closestRowAfter = chrAnnotationTree.min(annotation.getStart(), annotation.getEnd());
				Node<String> closestRowBefore = chrAnnotationTree.max(annotation.getStart(), annotation.getEnd());
				int distToAfterAnnotation = closestRowAfter == null ? Integer.MAX_VALUE : annotation.getDistanceTo(rowLocationMap.get(closestRowAfter.getValue()));
				int distToBeforeAnnotation = closestRowBefore == null ? Integer.MAX_VALUE : annotation.getDistanceTo(rowLocationMap.get(closestRowBefore.getValue()));
				if(closestRowAfter != null && distToAfterAnnotation < distToBeforeAnnotation) {
					closest = rowLocationMap.get(closestRowAfter.getValue());
				} else if (closestRowBefore != null) {
					closest = rowLocationMap.get(closestRowBefore.getValue());
				}
			}			
		}
		
		return closest;
		
	}
	
	/**
	 * 
	 * @param annotation
	 * @return the annotation corresponding to the row that is before in chromosome order
	 */
	public LightweightGenomicAnnotation getPriorClosest(LightweightGenomicAnnotation annotation) {
		IntervalTree<String> chrAnnotationTree = chrRowLocationTreeMap.get(annotation.getChromosome());
		LightweightGenomicAnnotation closest = null;
		Node<String> closestRowBefore = chrAnnotationTree.max(annotation.getStart(), annotation.getStart()+1);
		closest = closestRowBefore != null ? rowLocationMap.get(closestRowBefore.getValue()) : null;
		return closest;
		
	}
	
	/**
	 * Return 
	 * @param annotation
	 * @return the annotation corresponding to the row that is after the given annotation in chromosome order
	 */
	public LightweightGenomicAnnotation getNextClosest(LightweightGenomicAnnotation annotation) {
		IntervalTree<String> chrAnnotationTree = chrRowLocationTreeMap.get(annotation.getChromosome());
		LightweightGenomicAnnotation closest = null;
		Node<String> closestRowAfter = chrAnnotationTree.min(annotation.getEnd()-1, annotation.getEnd());
		closest = closestRowAfter != null ? rowLocationMap.get(closestRowAfter.getValue()) : null;
		return closest;
		
	}

}
