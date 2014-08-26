package broad.pda.differentialExpression;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;

import broad.core.datastructures.MatrixWithHeaders;
import broad.core.error.ParseException;
import broad.core.math.MathUtil;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;

/**
 * For a set of several experiments and a fixed set of control experiments this class wrapps
 * differential expression tools to perform analysis of different subsets of experiments against
 * the controls.
 * @author mgarber
 *
 */
public class FixedControlDifferentialExpression {
	private static String USAGE = "Does differential expression on the given matrix, it scores and evaluates the significance of each cell (or group of cells) in the matrix as to whether they are \"differentially expressed\" when compared with the fixed control set."+
	"\nParameters \n\t-in <Matrix with headers with data> "+
	"\n-controls <comma separated experiment column names  to be used as controls>" +
	"\n\t-log <If the data is not in log space use this flag to convert internally DATA IS ASSUMED TO BE IN LOG SPACE make sure to set this flag if it is not.> " +
	"\n\t-out <Output base name>" +
	"\n\t-minFoldChange <Minimum fold change in order to call significant, default is 2>" +
	"\n";
	public static void main (String [] args) throws IOException, ParseException {
		ArgumentMap argMap = CLUtil.getParameters(args, USAGE,"default");
		
		String outPrefix = argMap.getOutput();
		BufferedReader br = argMap.getInputReader();
		double minFoldChange = argMap.containsKey("minFoldChange") ? argMap.getDouble("minFoldChange") : 2;
		MatrixWithHeaders data = new MatrixWithHeaders(br);
		br.close();
		
		String controlStr = argMap.getMandatory("controls");
		boolean doLog = argMap.containsKey("log");
		
		if(doLog) {
			data.log2();
		}
		
		List<String> groups = null;
		HashMap<String, Collection<String>> groupMap = new LinkedHashMap<String, Collection<String>>();
		
		if(argMap.isPresent("groups")) {
			String [] groupNameArr = argMap.getMandatory("groups").split(",");
			groups = CLUtil.listFromArray(groupNameArr);
			for(int i = 0; i < groupNameArr.length; i++) {
				Collection<String> groupList = CLUtil.listFromArray(argMap.getMandatory("group"+i).split(","));
				groupMap.put(groupNameArr[i], groupList);
			}
		} else {
			groups = data.getColumnNames();
			for(String group : groups) {
				Collection<String> groupList = new ArrayList<String>(1);
				groupList.add(group);
				groupMap.put(group, groupList);
			}
		}
		
		String [] controlArray = controlStr.split(",");
		List<String> controlGroup = CLUtil.listFromArray(controlArray);
		
		MatrixWithHeaders dataScored = new MatrixWithHeaders(data.getRowNames(), groups);
		MatrixWithHeaders dataSignificance = new MatrixWithHeaders(data.getRowNames(), groups);

		for(String group : groups) {
			System.out.println("Doing group " + group);
			DifferentialExpression de = new DifferentialExpression(data, controlGroup, groupMap.get(group), true, false);
			for(int i = 0; i < data.rowDimension(); i++) {
				dataScored.set(i, group, de.testStatistics.get(i,0));
				double logFold = Math.abs(de.testStatistics.get(i,1));
				//dataSignificance.set(i, group, logFold > minFoldChange ? de.fdrMatrix.get(i, 0) : 0);
				dataSignificance.set(i, group, de.fdrMatrix.get(i, 0));
			}
		}
		if(argMap.containsKey("gct")) {
			dataScored.writeGCT(outPrefix + ".teststat.gct");
			dataSignificance.writeGCT(outPrefix + ".confidence.gct");
		}else {
			dataScored.write(outPrefix + ".teststat");
			dataSignificance.write(outPrefix + ".confidence");
		}
	}
}
