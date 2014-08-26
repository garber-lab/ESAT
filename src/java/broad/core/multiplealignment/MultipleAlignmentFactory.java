package broad.core.multiplealignment;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;

import broad.core.error.ParseException;

public class MultipleAlignmentFactory {
	
	public static MultipleAlignment create(String fileName, String format) throws IOException, ParseException {
		MultipleAlignmentIO maio = MultipleAlignmentIOFactory.create(format);
		MultipleAlignment ma = maio.load(fileName);
		ma.setIOHelper(maio);
		return ma;
	}
	
	public static MultipleAlignment create(InputStream in, String format) throws IOException, ParseException {
		MultipleAlignmentIO maio = MultipleAlignmentIOFactory.create(format);
		MultipleAlignment ma = maio.load(in);
		ma.setIOHelper(maio);
		return ma;
	}
	
	public static MultipleAlignment create(String format) {
		MultipleAlignmentIO maio = MultipleAlignmentIOFactory.create(format);
		MultipleAlignment ma = new MultipleAlignment();
		ma.setIOHelper(maio);
		return ma;
	}
	
	public static void main(String [] args) throws IOException, ParseException {
		String input = args[0];
		String output= args[1];

		
		MultipleAlignment ma = MultipleAlignmentFactory.create(input, "FASTA");
		
		BufferedWriter bw = new BufferedWriter(new FileWriter(output));
		ma.setIOHelper(MultipleAlignmentIOFactory.create("PHYLIP"));
		ma.write(bw);
		bw.close();
	}

}
