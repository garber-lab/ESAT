package umms.core.utils;

import java.util.Collections;
import java.util.HashMap;
import java.util.ArrayList;
import java.util.Set;
import java.util.List;
import java.io.File;

import umms.core.utils.InDropPreprocess;

public class ExperimentMap {
	private HashMap<String, Integer> exp2index = new HashMap<String, Integer>();   // experiment-to-index map
	private HashMap<Integer, String> index2exp = new HashMap<Integer, String>();   // index-to-experiment map
	private boolean singleCell = false;                   // single-cell analysis flag
	private HashMap<String,ArrayList<File>> fileMap = new HashMap<String, ArrayList<File>>();   // experiment-to-bamfile map
	private int nEntries = 0;

	/* empty constructor */
	public ExperimentMap() {
	}
	
	/* non-single-cell constructor */
	public ExperimentMap(HashMap<String,ArrayList<File>> bamFiles) {
		// not single-cell data:
		singleCell = false;
		// Create a mapping of experiments to indices, with names sorted alphabetically:
		List<String> expList = new ArrayList<String> (bamFiles.keySet());
		Collections.sort(expList);
		// Fill the exp2index and index2exp maps:
		int i=0;
		for (String e:expList) {
			exp2index.put(e, i);
			index2exp.put(i, e);
			i++;
		}
		// log the number of entries in the maps:
		nEntries = i;
		// fill the file map:
		fileMap = bamFiles;
	}
	
	/* inDrop single-cell constructor */
	public ExperimentMap(HashMap<String,ArrayList<File>> bamFiles, InDropPreprocess idInfo) {
		
		// single-cell data:
		singleCell = true;
		// Create a mapping of exp:barcodes to indices, with names sorted alphabetically:
		List<String> bcList = new ArrayList<String> (idInfo.getBcCounts().keySet());
		Collections.sort(bcList);
		// Fill the exp2index and index2exp maps:
		int i=0;
		for (String bc:bcList) {
			String bcStr = bc.toString();
			exp2index.put(bcStr, i);
			index2exp.put(i, bcStr);
			i++;
		}
		// log the number of entries in the maps:
		nEntries = i;
		// fill the file map:
		fileMap = bamFiles;
	}	

	public String getName(int idx) {
		return index2exp.get(idx);
	}
	
	public Integer getIndex(String s) {
		if (exp2index.containsKey(s)) {
			return exp2index.get(s);
		} else { 
			return -1;     // indicate that this <experiment>:<barcode> should be skipped
		}
	}	
	
	public Integer getNexp() {
		return nEntries;
	}

	public HashMap<String,ArrayList<File>> getBamFiles() {
		return fileMap;
	}
	
	public boolean isSingleCell() {
		return singleCell;
	}
}
