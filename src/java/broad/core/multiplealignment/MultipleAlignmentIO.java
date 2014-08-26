package broad.core.multiplealignment;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.InputStream;
import java.util.List;

import broad.core.error.ParseException;

public interface MultipleAlignmentIO {
	MultipleAlignment load(String fileName) throws IOException, ParseException;
	MultipleAlignment load(InputStream in) throws IOException, ParseException;
	MultipleAlignment load(String fileName, List<String> sequencesToLoad) throws IOException, ParseException;
	MultipleAlignment load(InputStream in, List<String> sequencesToLoad) throws IOException, ParseException;
	void write(BufferedWriter bw, MultipleAlignment ma) throws IOException;
	void write(BufferedWriter bw, MultipleAlignment ma, List<String> orderOfSequences) throws IOException ;
	String getPreferredFileExtension();
}
