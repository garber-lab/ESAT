package broad.core.annotation;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;

import nextgen.core.annotation.Annotation;

import broad.core.error.ParseException;

public class BasicAnnotationReader extends AnnotationReader<BasicGenomicAnnotation> {
	public BasicAnnotationReader() {
		super();
	}
	
	public BasicAnnotationReader(BufferedReader br) throws ParseException, IOException {
		super();
		load(br, AnnotationFactoryFactory.basicGenomicAnnotationfactory);
	}
	
	public BasicAnnotationReader(String input) throws IOException, ParseException {
		super();
		load(new File(input), AnnotationFactoryFactory.basicGenomicAnnotationfactory);
	}

	@Override	public BasicGenomicAnnotation createAnnotation(GenomicAnnotation a) {
		return new BasicGenomicAnnotation(a);
	}
	
	public int parse(BufferedReader br, GenomicAnnotationFilter<BasicGenomicAnnotation> filter, AnnotationHandler handler) throws ParseException {
		return super.parse(br, AnnotationFactoryFactory.basicGenomicAnnotationfactory, filter, handler);
	}
	
	public int parse(String file, GenomicAnnotationFilter<BasicGenomicAnnotation> filter, AnnotationHandler handler) throws ParseException, IOException {
		return super.parse(file, AnnotationFactoryFactory.basicGenomicAnnotationfactory, filter, handler);
	}

	public BasicGenomicAnnotation createAnnotation(Annotation ac) {
		throw new UnsupportedOperationException("TODO");
	}
	
	
}
