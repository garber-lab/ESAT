package broad.core.annotation;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;

import broad.core.error.ParseException;

public class GFFReader extends AnnotationReader<GFF>{
	
	public GFFReader() {
		super();
	}
	
	public GFFReader(BufferedReader br) throws ParseException, IOException {
		super();
		load(br, AnnotationFactoryFactory.gffFactory);
	}
	
	public GFFReader(String input) throws IOException, ParseException {
		super();
		load(new File(input), AnnotationFactoryFactory.gffFactory);
	}
	
	public GFFReader(String input, GenomicAnnotationFilter<GFF> filter) throws IOException, ParseException {
		super();
		load(new File(input), AnnotationFactoryFactory.gffFactory, filter);
	}

	@Override public GFF createAnnotation(GenomicAnnotation a) {
		return new GFF(a);
	}
	
	public int parse(BufferedReader br, GenomicAnnotationFilter<GFF> filter, AnnotationHandler handler) throws ParseException {
		return super.parse(br, AnnotationFactoryFactory.gffFactory, filter, handler);
	}

	@Override
	public int  parse(String file, GenomicAnnotationFilter<GFF> filter, AnnotationHandler handler) throws ParseException, IOException {
		return super.parse(file, AnnotationFactoryFactory.gffFactory, filter, handler);
		
	}
	
	

}
