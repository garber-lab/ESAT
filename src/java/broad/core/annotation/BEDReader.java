package broad.core.annotation;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;

import broad.core.error.ParseException;

public class BEDReader extends AnnotationReader<BED>{
	
	public BEDReader() {
		super();
	}
	
	public BEDReader(BufferedReader br) throws ParseException, IOException {
		super();
		load(br, AnnotationFactoryFactory.bedFactory);
	}
	
	public BEDReader(String input) throws ParseException, IOException {
		super();
		load(new File(input), AnnotationFactoryFactory.bedFactory);
	}
	
	public void load(BufferedReader br) throws IOException, ParseException {
		super.load(br, AnnotationFactoryFactory.bedFactory);
	}
	public void load(BufferedReader br,GenomicAnnotationFilter<BED> filter) throws IOException, ParseException {
		super.load(br, AnnotationFactoryFactory.bedFactory, filter);
	}
	
	public void load(File  inFile ,GenomicAnnotationFilter<BED> filter) throws IOException, ParseException {
		super.load(inFile, AnnotationFactoryFactory.bedFactory, filter);
	}
	
	public int parse(BufferedReader br, GenomicAnnotationFilter<BED> filter, AnnotationHandler handler) throws ParseException {
		return super.parse(br, AnnotationFactoryFactory.bedFactory, filter, handler);
	}

	public BED createAnnotation(GenomicAnnotation a) {
		return new BED(a);
	}

	public int parse(String file, GenomicAnnotationFilter<BED> filter,
			AnnotationHandler handler) throws ParseException, IOException {
		return super.parse(file,AnnotationFactoryFactory.bedFactory, filter, handler);
		
	}

}
