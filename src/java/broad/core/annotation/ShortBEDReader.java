package broad.core.annotation;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;

import broad.core.error.ParseException;

public class ShortBEDReader extends AnnotationReader<ShortBED>{
	
	public ShortBEDReader() {
		super();
	}
	
	public ShortBEDReader(BufferedReader br) throws ParseException, IOException {
		super();
		load(br, AnnotationFactoryFactory.shortBEDFactory);
	}
	
	public ShortBEDReader(String input) throws ParseException, IOException {
		super();
		load(new File(input), AnnotationFactoryFactory.shortBEDFactory);
	}
	
	public void load(BufferedReader br) throws IOException, ParseException {
		super.load(br, AnnotationFactoryFactory.shortBEDFactory);
	}
	
	public int parse(BufferedReader br, GenomicAnnotationFilter<ShortBED> filter, AnnotationHandler handler) throws ParseException {
		return super.parse(br, AnnotationFactoryFactory.shortBEDFactory, filter, handler);
	}

	public ShortBED createAnnotation(GenomicAnnotation a) {
		return new ShortBED(a);
	}

	public int parse(String file, GenomicAnnotationFilter<ShortBED> filter,
			AnnotationHandler handler) throws ParseException, IOException {
		return super.parse(file,AnnotationFactoryFactory.shortBEDFactory, filter, handler);
		
	}

}
