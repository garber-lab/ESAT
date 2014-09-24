package umms.core.readers;

import java.io.File;
import java.io.IOException;
import java.io.BufferedReader;
import java.io.FileReader;
import java.util.Hashtable;

import org.apache.log4j.Logger;

public class MappingTableReader {
	
	public Logger logger;
	private BufferedReader reader;
	private String headerString;    // the file header (assumed to be the first line
	private Hashtable<String, Integer> headerFields;
	private String sep;
	
	// basic constructor
	public MappingTableReader(final File mapFile, String delimiter) throws IOException {
		reader = new BufferedReader(new FileReader(mapFile));
		sep = delimiter;
		headerString = reader.readLine();    // the header string is saved and the first valid line will be considered the line after the header
		// parse the header for later use:
		headerFields = new Hashtable<String, Integer>();   // save the header fields and their column position
		// Split the header string into components and fill the headerFields
		String[] fields = headerString.split(sep);
		for (int i=0; i<fields.length; i++) {
			headerFields.put(fields[i], i);
		}
    }

	// default method (tab-delimited file)
	public MappingTableReader(final File mapFile) throws IOException {
		this(mapFile, "\t");
	}
	
	// METHODS
	public boolean hasMandatoryFields(String[] requiredFields) {
		boolean result = true;
		
		// Test to make sure that all mandatory fields are present:
		for (String field:requiredFields) {
			if (!headerFields.containsKey(field)) {
				logger.warn("Gene mapping file header missing mandatory field: "+field);
				result = false;
			}
		}
		
		return result;
	}
	
	public Hashtable<String, Integer> getColumnMapping() {
		return headerFields;
	}
	
	public String readLine() throws IOException {
		return reader.readLine();
	}
	
	public String[] readOrderedFieldsFromLine(String[] fields) throws IOException {
		// returns the values from the columns specified in the "fields" array, in the same order
		String[] outStr = new String[fields.length];
		String line = reader.readLine();
		if (line==null) {
			outStr = new String[0];
		} else {
			String[] lineFields = line.split(sep);
			for (int i=0; i<fields.length; i++) {
				int idx = headerFields.get(fields[i]);    // extract the data from the correct column
				outStr[i] = lineFields[idx];              // copy it, in order
			}
		}
		
		return outStr;
	}
		
}
