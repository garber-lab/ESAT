package broad.core.annotation;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.Collection;

import nextgen.core.alignment.Alignment;
import nextgen.core.feature.Window;

import broad.core.error.ParseException;
import broad.pda.seq.segmentation.ReadFilter;

public class BEDGraphReader extends AnnotationReader<BEDGraph> {
	
	public BEDGraphReader() {
		super();
	}
	
	public BEDGraphReader(BufferedReader br) throws ParseException, IOException {
		super();
		load(br, AnnotationFactoryFactory.bedGraphFactory);
	}
	
	public BEDGraphReader(String input) throws ParseException, IOException {
		super();
		load(new File(input), AnnotationFactoryFactory.bedGraphFactory);
	}
	
	public void load(BufferedReader br) throws IOException, ParseException {
		super.load(br, AnnotationFactoryFactory.bedGraphFactory);
	}
	
	public int parse(BufferedReader br, GenomicAnnotationFilter<BEDGraph> filter, AnnotationHandler handler) throws ParseException {
		return super.parse(br, AnnotationFactoryFactory.bedGraphFactory, filter, handler);
	}

	public BEDGraph createAnnotation(GenomicAnnotation a) {
		return new BEDGraph(a);
	}

	public int parse(String file, GenomicAnnotationFilter<BEDGraph> filter,
			AnnotationHandler handler) throws ParseException, IOException {
		return super.parse(file,AnnotationFactoryFactory.bedGraphFactory, filter, handler);
		
	}

	
}
