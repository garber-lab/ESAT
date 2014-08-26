package broad.core.util;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import broad.core.annotation.AnnotationReader;
import broad.core.annotation.AnnotationReaderFactory;
import broad.core.annotation.BasicGenomicAnnotation;
import broad.core.annotation.GenomicAnnotation;
import broad.core.error.ParseException;
import broad.core.math.MathUtil;

public class CLUtil {

	public CLUtil() {
		super();
	}
	
	public static ArgumentMap getParameters(String [] args, String usage, String defaultTask) {
		ArgumentMap argMap = new ArgumentMap(args.length, usage, defaultTask);
		for(int i = 0; i < args.length; i++) {
			if(args[i].startsWith("-")) {
				String key = args[i].substring(1);
				String val = "";
				if(i + 1 < args.length && !args[i + 1].startsWith("-")) {
					val = args[i + 1];
					i++;
				}
				//System.out.println("Processing " + key + ", " + val);
				if(argMap.containsKey(key)) { // hey this is multivalued
					List<String> values = new ArrayList<String>();
					values.add(argMap.get(key));
					values.add(val);
				} 
				argMap.put(key, val);
				//System.out.println(argMap.get(key));
			} else { //ONLY for backwards compatibility, before we implemented unix type parameter passing.
				String[] arg = args[i].split("=");
				//System.out.println("processing " +arg[i]);
				if(arg.length != 2){
					throw new IllegalArgumentException(usage);
				}
				argMap.put(arg[0],arg[1]);
		
			}
		}
		return argMap;
	}
	
	public static ArgumentMap getParameters(String [] args, String usage) {

		return getParameters(args, usage, null);
	}
	
	public static class ArgumentMap extends HashMap<String, List<String>> {
		private static final long serialVersionUID = 2312363L;
		private String usage;
		private String defaultTask;
		private String task;
		private String input;
		private String output;
		private String inputDir;
		private String outputDir;
		
		public ArgumentMap(int size, String usage, String defaultTask) {
			super(size);
			this.usage = usage;
			this.defaultTask = defaultTask;
		}
		
		public String get(String key) {
			List<String> result = super.get(key);
			return result == null ||  result.size() == 0 ? "" : result.get(0);
		}
		
		public String get(String key, String defaultValue) {
			return super.containsKey(key) ? get(key) : defaultValue;
		}
		
		public void put(String key, String value) {
			if(key.toLowerCase().equals("task")) {
				this.task = value;
			} else if (key.toLowerCase().equals("in")) {
				this.input = value;
			} else if (key.toLowerCase().equals("out")) {
				this.output = value;
			} else if (key.toLowerCase().equals("outdir")) {
				this.outputDir = value;
			} else if (key.toLowerCase().equals("indir")) {
				this.inputDir = value;
			} else {
				List<String> values = super.get(key);
				if(values == null) {
					values = new ArrayList<String>();
					super.put(key, values);
				}
				values.add(value);
			}
					
		}
		
		public String getTask() {
			if(task == null && defaultTask == null) {
				throw new IllegalArgumentException("Missing task\n"+ usage);
			}
			return task == null ? defaultTask : task;
		}
		

		public String getInput() {
			if(input == null) {
				throw new IllegalArgumentException("Must provide \"in\"" + usage);
			}
			return input;
		}
		
		public boolean hasInputFile() {
			return input != null;
		}
		
		public BufferedReader getInputReader() throws IOException {
			BufferedReader br = null;
			if(input != null) {
				br = new BufferedReader(new FileReader(input));
			} else {
				br = new BufferedReader(new InputStreamReader(System.in));
			}
			
			return br;
		}
		
		public InputStream getInputStream() throws IOException {
			InputStream is = null;
			if(input != null) {
				is = new FileInputStream(input);
			} else {
				is = System.in;
			}
			
			return is;
		}
		
		public String getOutput() {
			if(output == null) {
				throw new IllegalArgumentException("Must provide an \"out\"\n" + usage);
			}
			return output;
		}
		
		public BufferedWriter getOutputWriter() throws IOException {
			BufferedWriter bw = null;
			if(output != null) {
				bw = new BufferedWriter(new FileWriter(output));
			} else {
				bw = new BufferedWriter(new OutputStreamWriter(System.out));
			}
			
			return bw;			
		}
		
		public OutputStream getOutputStream() throws IOException {
			OutputStream os = null;
			if(output != null) {
				os = new FileOutputStream(output);
			} else {
				os = System.out;
			}
			
			return os;
		}
		
		public boolean isOutputSet() {
			return output != null;
		}
		
		public String getOutputDir() {
			if(outputDir == null) {
				outputDir = "./";
			}
				
			return outputDir;
		}
		
		public String getInputDir() {
			if(inputDir == null) {
				throw new IllegalArgumentException("Must provide an \"indir\"\n" + usage);
			}
			return inputDir;
		}
		
		public List<String> getAll(String key) {
			
			return super.get(key) == null ? new ArrayList<String>() : super.get(key);
		}
		
		public   Map<String, List<? extends GenomicAnnotation>>  getRegionMapFromParameters()
		throws ParseException, IOException {
			Map<String, List<? extends GenomicAnnotation>> regionChrMap = new LinkedHashMap<String, List<? extends GenomicAnnotation>>();
			if(containsKey("regions")) {
				String annotationFile = getMandatory("regions");
				String annotationFileFormat = containsKey("regionFormat") ?getMandatory("regionFormat") : "BED";

				AnnotationReader<? extends GenomicAnnotation> reader = AnnotationReaderFactory.create(annotationFile, annotationFileFormat);
				Iterator<String> chrIt = reader.getChromosomeIterator();
				while(chrIt.hasNext()) {
					String chr = chrIt.next();
					regionChrMap.put(chr, reader.getChromosomeBEDs(chr));	
				}
			} else {
				int start = getInteger("start");
				int end   = getInteger("end");
				String chr = getMandatory("chr").replace("chr", "");
				BasicGenomicAnnotation annotation = new BasicGenomicAnnotation("a", chr, start, end); 
				List<GenomicAnnotation> regionListTmp = new ArrayList<GenomicAnnotation>();
				regionListTmp.add(annotation);
				regionChrMap.put(chr, regionListTmp);			
			}
			return regionChrMap;
		}

		
		
		/**
		 * 
		 * @param key - the key whose value is presumably an integer
		 * @return An integer representing the value.
		 * @throws NumberFormatException - if the value could not be converted to an integer.
		 */
		public int getInteger(String key) throws NumberFormatException {
			String val = getMandatory(key);
			return Integer.parseInt(val);
		}
		
		public int getInteger(String key, int defaultValue) throws NumberFormatException {
			return super.containsKey(key) ? getInteger(key) : defaultValue;
		}
		
		public List<Integer> getIntegers(String key) throws NumberFormatException{
			List<String> stringVals = getAllMandatory(key);
			List<Integer> paramList = new ArrayList<Integer>(stringVals.size());
			
			Iterator<String> valIt = stringVals.iterator();
			while(valIt.hasNext()) {
				String val = valIt.next();
				paramList.add(Integer.parseInt(val));
			}
			
			return paramList;
		}
		
		public String getMandatory(String key)  throws IllegalArgumentException{
			List<String> parameter = super.get(key);
			if(parameter == null || parameter.size() == 0) {
				throw new IllegalArgumentException("Argument "+key+" is mandatory\n"+usage);
			}
			return parameter.get(0);
		}
		
		public List<String> getAllMandatory(String key) throws IllegalArgumentException{
			List<String> vals = super.get(key);
			if(vals == null || vals.size() == 0) {
				throw new IllegalArgumentException("Argument "+key+" is mandatory, at least one should be provided\n"+usage);
			}
			
			return vals;
		}
		public float getFloat(String key) throws NumberFormatException{
			String val = getMandatory(key);
			return Float.parseFloat(val);
		}
		
		public float getFloat(String key, float defaultValue) throws NumberFormatException {
			return super.containsKey(key) ? getFloat(key) : defaultValue;
		}
		
		public List<Float> getFloats(String key) throws NumberFormatException{
			List<String> stringVals = getAllMandatory(key);
			List<Float> paramList = new ArrayList<Float>(stringVals.size());
			
			Iterator<String> valIt = stringVals.iterator();
			while(valIt.hasNext()) {
				String val = valIt.next();
				paramList.add(Float.parseFloat(val));
			}
			
			return paramList;
		}

		public boolean isPresent(String key) {
			return super.get(key) != null && super.get(key).size() > 0;
		}

		public boolean isFlagTrue(String key) {
			return containsKey(key) && "TRUE".equalsIgnoreCase(get(key));
		}

		public double getDouble(String key) {
			String val = getMandatory(key);
			return Double.parseDouble(val);
		}
		
		public double getDouble(String key, double defaultValue) throws NumberFormatException {
			return super.containsKey(key) ? getDouble(key) : defaultValue;
		}
		
		public List<Double> getDoubles(String key) throws NumberFormatException{
			List<String> stringVals = getAllMandatory(key);
			List<Double> paramList = new ArrayList<Double>(stringVals.size());
			
			Iterator<String> valIt = stringVals.iterator();
			while(valIt.hasNext()) {
				String val = valIt.next();
				paramList.add(Double.parseDouble(val));
			}
			
			return paramList;
		}
		
		public String toArgString() {
			String result = "";
			for (String key : keySet()) {
				result = result + " -" + key + " " + get(key);
			}
			
			// for some reason these aren't stored in the hashmap
			if (this.task != null) {
				result = result + " -task " + this.task;
			}
			if (this.input != null) {
				result = result + " -in " + this.input;
			}
			if (this.output != null) {
				result = result + " -out " + this.output;
			}
			if (this.inputDir != null) {
				result = result + " -indir " + this.inputDir;
			}
			if (this.outputDir != null) {
				result = result + " -outdir " + this.outputDir;
			}
			
			return result;
		}
	}
	
	public static void writeLeftJustifiedField(BufferedWriter bw, String value, int fieldSize) throws IOException {
		if(value.length() > fieldSize) {
			bw.write(value.substring(0,fieldSize));
		} else {
			bw.write(value);
			for( int i = 0; i < (fieldSize - value.length()); i++) {
				bw.write(" ");
			}
		}
	}
	
	public static void writeRightJustifiedField(BufferedWriter bw, String value, int fieldSize) throws IOException {
		if(value.length() > fieldSize) {
			bw.write(value.substring(0,fieldSize));
		} else {
			for( int i = 0; i < (fieldSize - value.length()); i++) {
				bw.write(" ");
			}
			bw.write(value);
		}
	}
	
	public static String replaceFileExtension(String fileName, String newExtension) {
		String [] fileNameComponents = fileName.split("\\.");
		StringBuilder newFileName = new StringBuilder(fileNameComponents[0]); 
		
		for(int i = 1; i < fileNameComponents.length - 1; i++) {
			newFileName.append(".").append(fileNameComponents[i]);
		}
		
		newFileName.append(newExtension);
		return newFileName.toString();
	}

	public static double log2(double x) {
		return MathUtil.log2(x);
	}
	
	public static <T> List<T> listFromArray(T[] array) {
		List<T> list = new ArrayList<T>(array.length);
		for(int i = 0; i < array.length; i++) {
			list.add(array[i]);
		}
		return list;
	}
}
