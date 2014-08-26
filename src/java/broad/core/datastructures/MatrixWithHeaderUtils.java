package broad.core.datastructures;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;

import broad.core.error.ParseException;
import broad.core.hmm.BadModelException;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;

public class MatrixWithHeaderUtils {
	static final String USAGE = "\nTasks" +
	"\n quantileNormalizeCols  \n\t\t-in <input matrix (in GCT or regular format) standard input assumed if this parameter is not present> \n\t\t-out <Output file (or standard output if non is supplied, use .gct extension if such format is desired>"+
	"\n medianCenterColumns \n\t\t-in <input matrix (in GCT or regular format) standard input assumed if this parameter is not present> \n\t\t-out <Output file (or standard output if non is supplied, use .gct extension if such format is desired>"+
	"\n medianCenterRows \n\t\t-in <input matrix (in GCT or regular format) standard input assumed if this parameter is not present> \n\t\t-out <Output file (or standard output if non is supplied, use .gct extension if such format is desired>"+
	"\n meanNormalizeColumns \n\t\t-in <input matrix (in GCT or regular format) standard input assumed if this parameter is not present> \n\t\t-out <Output file (or standard output if non is supplied, use .gct extension if such format is desired>"+
	"\n meanNormalizeRows \n\t\t-in <input matrix (in GCT or regular format) standard input assumed if this parameter is not present> \n\t\t-out <Output file (or standard output if non is supplied, use .gct extension if such format is desired>"+
	"\n densityNormalizeColumns  \n\t\t-in <input matrix (in GCT or regular format) standard input assumed if this parameter is not present> \n\t\t-out <Output file (or standard output if non is supplied, use .gct extension if such format is desired>"+
	"\n densityNormalizeRows  \n\t\t-in <input matrix (in GCT or regular format) standard input assumed if this parameter is not present> \n\t\t-out <Output file (or standard output if non is supplied, use .gct extension if such format is desired>"+
	"\n zNormalizeColumns  \n\t\t-in <input matrix (in GCT or regular format) standard input assumed if this parameter is not present> \n\t\t-out <Output file (or standard output if non is supplied, use .gct extension if such format is desired>"+
	"\n zNormalizeRows  \n\t\t-in <input matrix (in GCT or regular format) standard input assumed if this parameter is not present> \n\t\t-out <Output file (or standard output if non is supplied, use .gct extension if such format is desired>"+
	"\n";
	
	public static void main (String [] args) throws IOException, ParseException, BadModelException {
		
		ArgumentMap argMap = CLUtil.getParameters(args,USAGE , "quantileNormalizeCols");
		 if ("quantileNormalizeCols".equalsIgnoreCase(argMap.getTask())) {
			 BufferedReader br = argMap.getInputReader();
			 MatrixWithHeaders in = new MatrixWithHeaders(br);
			 br.close();
			 
			 in.quantileNormalizeColumns();
			 BufferedWriter bw = argMap.getOutputWriter();
			 in.write(bw);
			 bw.close();
		 }else if ("mediancentercolumns".equalsIgnoreCase(argMap.getTask())) {
			 BufferedReader br = argMap.getInputReader();
			 MatrixWithHeaders in = new MatrixWithHeaders(br);
			 br.close();
			 
			 in.medianCenterColumnsThis();
			 BufferedWriter bw = argMap.getOutputWriter();
			 in.write(bw);
			 bw.close(); 
		 }else if ("mediancenterRows".equalsIgnoreCase(argMap.getTask())) {
			 BufferedReader br = argMap.getInputReader();
			 MatrixWithHeaders in = new MatrixWithHeaders(br);
			 br.close();
			 
			 in.medianCenterRowsThis();
			 BufferedWriter bw = argMap.getOutputWriter();
			 in.write(bw);
			 bw.close(); 
		 }else if ("meannormalizecolumns".equalsIgnoreCase(argMap.getTask())) {
			 BufferedReader br = argMap.getInputReader();
			 MatrixWithHeaders in = new MatrixWithHeaders(br);
			 br.close();
			 
			 in.meanNormalizeColumns();
			 BufferedWriter bw = argMap.getOutputWriter();
			 in.write(bw);
			 bw.close(); 
		 }else if ("meannormalizeRows".equalsIgnoreCase(argMap.getTask())) {
			 BufferedReader br = argMap.getInputReader();
			 MatrixWithHeaders in = new MatrixWithHeaders(br);
			 br.close();
			 
			 in.meanNormalizeRows();
			 BufferedWriter bw = argMap.getOutputWriter();
			 in.write(bw);
			 bw.close(); 
		 }else if ("densitynormalizecolumns".equalsIgnoreCase(argMap.getTask())) {
			 BufferedReader br = argMap.getInputReader();
			 MatrixWithHeaders in = new MatrixWithHeaders(br);
			 br.close();
			 
			 in.columnDensityNormalization();
			 BufferedWriter bw = argMap.getOutputWriter();
			 in.write(bw);
			 bw.close(); 
		 }else if ("densitynormalizeRows".equalsIgnoreCase(argMap.getTask())) {
			 BufferedReader br = argMap.getInputReader();
			 MatrixWithHeaders in = new MatrixWithHeaders(br);
			 br.close();
			 
			 in.rowDensityNormalization();
			 BufferedWriter bw = argMap.getOutputWriter();
			 in.write(bw);
			 bw.close(); 
		 }else if("znormalizeRows".equalsIgnoreCase(argMap.getTask())) {
			 boolean log = argMap.containsKey("log");
			 BufferedReader br = argMap.getInputReader();
			 MatrixWithHeaders in = new MatrixWithHeaders(br);
			 br.close();
			 
			 if(log) {
				 in.log2();
			 }
			 in.zRowNormalize();
			 BufferedWriter bw = argMap.getOutputWriter();
			 in.write(bw);
			 bw.close(); 
		 }else if("znormalizeColumns".equalsIgnoreCase(argMap.getTask())) {
			 boolean log = argMap.containsKey("log");
			 BufferedReader br = argMap.getInputReader();
			 MatrixWithHeaders in = new MatrixWithHeaders(br);
			 br.close();
			 
			 if(log) {
				 in.log2();
			 }
			 in.zColumnNormalize();
			 BufferedWriter bw = argMap.getOutputWriter();
			 in.write(bw);
			 bw.close(); 
		 }
		 
	}
}
