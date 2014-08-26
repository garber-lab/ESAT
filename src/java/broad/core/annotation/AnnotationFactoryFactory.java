package broad.core.annotation;

import broad.core.error.ParseException;

public class AnnotationFactoryFactory  {
	
	public static BEDFactory bedFactory = new BEDFactory();
	public static GFFFactory gffFactory = new GFFFactory();
	public static ShortBEDFactory shortBEDFactory = new ShortBEDFactory();
	public static WIGFactory wigFactory = new WIGFactory();
	public static BasicGenomicFactory basicGenomicAnnotationfactory = new BasicGenomicFactory();
	public static RepeatMaskerAnnotationFactory repeatMaskerAnnotationFactory = new RepeatMaskerAnnotationFactory();
	public static BEDGraphFactory bedGraphFactory = new BEDGraphFactory();
	
	public static AnnotationFactory<? extends GenomicAnnotation> getFactory(String annotationType) {
		//System.out.print("getFactory for type: " + annotationType + " ");
		if("BED".equalsIgnoreCase(annotationType)) {
			return bedFactory;
		} else if("GFF".equalsIgnoreCase(annotationType)) {
			return gffFactory;
		} else if ("ShortBED".equalsIgnoreCase(annotationType)) {
			return shortBEDFactory;
		} else if ("BEDGraph".equalsIgnoreCase(annotationType)) {
			return bedGraphFactory;
		} else {
			return basicGenomicAnnotationfactory;
		}
		
	}
	public static class BEDFactory implements AnnotationFactory<BED> {

		public BED create(String[] rawFields) throws ParseException {
			return new BED(rawFields);
		}
		
		public BED create(GenomicAnnotation a) {
			return new BED(a);
		}
		
		public BED create(String name) {
			return new BED(name);
		}
		
	}
	
	public static class WIGFactory implements AnnotationFactory<BED> {

		public BED create(String[] rawFields) throws ParseException {
			String chr = rawFields[0];
			int start = Integer.parseInt(rawFields[1]);
			int end = Integer.parseInt(rawFields[2]);
			double score = Double.parseDouble(rawFields[2]);
			BED bed = bedFactory.create(chr+":"+start+"-"+end);
			bed.setChromosome(chr);
			bed.setStart(start);
			bed.setEnd(end);
			bed.setScore(score);
			return  bed;
		}
		
		public BED create(GenomicAnnotation a) {
			return new BED(a);
		}
		
		public BED create(String name) {
			return new BED(name);
		}
		
	}
	
	public static class ShortBEDFactory implements AnnotationFactory<ShortBED> {

		public ShortBED create(String[] rawFields) throws ParseException {
			return new ShortBED(rawFields);
		}
		
		public ShortBED create(GenomicAnnotation a) {
			return new ShortBED(a);
		}
		
		public ShortBED create(String name) {
			return new ShortBED(name);
		}
	}
	
	public static class BEDGraphFactory implements AnnotationFactory<BEDGraph> {

		public BEDGraph create(String[] rawFields) throws ParseException {
			return new BEDGraph(rawFields);
		}
		
		public BEDGraph create(GenomicAnnotation a) {
			return new BEDGraph(a);
		}
		
		public BEDGraph create(String name) {
			return new BEDGraph(name);
		}
	}
	
	public static class GFFFactory implements AnnotationFactory<GFF> {

		public GFF create(String[] rawFields) throws ParseException {
			return new GFF(rawFields);
		}
		
		public GFF create(GenomicAnnotation a) {
			return new GFF(a);
		}
		
		public GFF create(String name) {
			return new GFF(name);
		}
		
	}
	
	public static class BasicGenomicFactory implements AnnotationFactory<BasicGenomicAnnotation> {

		public BasicGenomicAnnotation create(String[] rawFields) throws ParseException {
			return new BasicGenomicAnnotation(rawFields);
		}
		
		public BasicGenomicAnnotation create(GenomicAnnotation a) {
			return new BasicGenomicAnnotation(a);
		}
		
		public BasicGenomicAnnotation create(String name) {
			return new BasicGenomicAnnotation(name);
		}
	}
	
	
	public static class RepeatMaskerAnnotationFactory implements AnnotationFactory<RepeatMaskerAnnotation> {
		public RepeatMaskerAnnotation create(String[] rawFields) throws ParseException {
			return new RepeatMaskerAnnotation(rawFields);
		}
		
		public RepeatMaskerAnnotation create(GenomicAnnotation a) {
			return new RepeatMaskerAnnotation(a);
		}
		
		public RepeatMaskerAnnotation create(String name) {
			return new RepeatMaskerAnnotation(name);
		}
	}
	
}
