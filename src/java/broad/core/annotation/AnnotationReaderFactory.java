package broad.core.annotation;

import java.io.IOException;

import broad.core.error.ParseException;

public class AnnotationReaderFactory {
	
	public static AnnotationReader<? extends GenomicAnnotation> create(String sourceFile, String annotationType) throws ParseException, IOException {
		AnnotationReader<? extends GenomicAnnotation> r;
		if("BED".equals(annotationType)) {
			r = new BEDReader(sourceFile);
		} else if("GFF".equals(annotationType)) {
			r = new GFFReader(sourceFile);
		} else if ("BEDGraph".equalsIgnoreCase(annotationType)) {
			r = new BEDGraphReader(sourceFile);
		} else {
			r = new BasicAnnotationReader(sourceFile);
		}
		
		return r;
	}
	
	public static AnnotationReader<? extends GenomicAnnotation> create(String annotationType) throws ParseException, IOException {
		AnnotationReader<? extends GenomicAnnotation> r;
		if("BED".equals(annotationType)) {
			r = new BEDReader();
		} else if("GFF".equals(annotationType)) {
			r = new GFFReader();
		} else {
			r = new BasicAnnotationReader();
		}
		
		return r;
	}
	
	public static AnnotationReader<? extends GenomicAnnotation> create(String sourceFile, String annotationType, final double minScore) throws ParseException, IOException {
		AnnotationReader<? extends GenomicAnnotation> r;
		if("BED".equals(annotationType)) {
			BEDReader bedR = new BEDReader();
			bedR.parse(sourceFile, new GenomicAnnotationFilter<BED>() {

				public boolean accept(BED annotation) { return annotation.getScore() > minScore;}

				public boolean isEnough(BED annotation) {return false;}
				
			}, new SimpleAnnotationHandler(bedR));
			r = bedR;
		} else if("GFF".equals(annotationType)) {
			GFFReader gffR = new GFFReader();
			gffR.parse(sourceFile, new GenomicAnnotationFilter<GFF>() {

				public boolean accept(GFF annotation) {	return annotation.getScore() > minScore;}

				public boolean isEnough(GFF annotation) {return false;}
				
			}, new SimpleAnnotationHandler(gffR));
			r = gffR;
		} else if("BEDGraph".equals(annotationType)) {
				BEDGraphReader bedR = new BEDGraphReader();
				bedR.parse(sourceFile, new GenomicAnnotationFilter<BEDGraph>() {

					public boolean accept(BEDGraph annotation) { return annotation.getScore() > minScore;}

					public boolean isEnough(BEDGraph annotation) {return false;}
					
				}, new SimpleAnnotationHandler(bedR));
				r = bedR;
		} else {
			BasicAnnotationReader basicR = new BasicAnnotationReader(); 
			basicR.parse(sourceFile, new GenomicAnnotationFilter<BasicGenomicAnnotation>() {

				public boolean accept(BasicGenomicAnnotation annotation) {	return annotation.getScore() > minScore;}

				public boolean isEnough(BasicGenomicAnnotation annotation) {return false;}
				
			}, new SimpleAnnotationHandler(basicR));
			r = basicR;
		}
		
		return r;
	}
	
	// Jesse June 13, 2012
	public static AnnotationReader<? extends GenomicAnnotation> create(String sourceFile, String annotationType, final String chromosome) throws ParseException, IOException {
		AnnotationReader<? extends GenomicAnnotation> r;
		
		if (chromosome == null || "all".equals(chromosome)) {
			return create(sourceFile, annotationType);
		}
		
		if("BED".equals(annotationType)) {
			BEDReader bedR = new BEDReader();
			bedR.parse(sourceFile, new GenomicAnnotationFilter<BED>() {

				public boolean accept(BED annotation) { return chromosome.equals(annotation.getChromosome());}

				public boolean isEnough(BED annotation) {return false;}
				
			}, new SimpleAnnotationHandler(bedR));
			r = bedR;
		} else if("GFF".equals(annotationType)) {
			GFFReader gffR = new GFFReader();
			gffR.parse(sourceFile, new GenomicAnnotationFilter<GFF>() {

				public boolean accept(GFF annotation) {	return chromosome.equals(annotation.getChromosome());}

				public boolean isEnough(GFF annotation) {return false;}
				
			}, new SimpleAnnotationHandler(gffR));
			r = gffR;
		} else {
			BasicAnnotationReader basicR = new BasicAnnotationReader(); 
			basicR.parse(sourceFile, new GenomicAnnotationFilter<BasicGenomicAnnotation>() {

				public boolean accept(BasicGenomicAnnotation annotation) {	return chromosome.equals(annotation.getChromosome());}

				public boolean isEnough(BasicGenomicAnnotation annotation) {return false;}
				
			}, new SimpleAnnotationHandler(basicR));
			r = basicR;
		}
		
		return r;
	}

}

class SimpleAnnotationHandler  implements AnnotationHandler {
	AnnotationReader reader;
	public SimpleAnnotationHandler(AnnotationReader r) {
		this.reader = r;
	}
	public void track(String line) {
		reader.startTrack(line);
		
	}

	public void browserLine(String line) {
		reader.addBrowserLine(line);
	}

	public void annotation(GenomicAnnotation annotation) {
		reader.addAnnotation(annotation);
		
	}

	@Override
	public void eof() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void begin() {
		// TODO Auto-generated method stub
		
	}
	
}
