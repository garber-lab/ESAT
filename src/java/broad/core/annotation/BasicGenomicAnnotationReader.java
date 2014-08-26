package broad.core.annotation;

import java.awt.Graphics2D;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import broad.core.annotation.AnnotationFactory;
import broad.core.annotation.AnnotationFactoryFactory;
import broad.core.annotation.AnnotationReader;
import broad.core.annotation.BasicGenomicAnnotation;
import broad.core.annotation.GenomicAnnotation;
import broad.core.annotation.GenomicAnnotationFilter;
import broad.core.error.ParseException;


public class BasicGenomicAnnotationReader extends AnnotationReader<GenomicAnnotation>{
	private List<GenomicAnnotation> annotationList;
	protected File source;
	public static int RIGHT_MARGIN = 50;
	public static int VERTICAL_SPACING = 100;
	public static int INTER_VERTICAL_SPACING = 40;
	public static final Pattern REG_PATTERN = Pattern.compile("(chr)(\\w{1,2})(:)([0-9]+)(-)([0-9]+)");
	public static final Pattern POINT_PATTERN = Pattern.compile("(chr)(\\w{1,2})(:)([0-9]+)");
	String annotationType;
	
	public BasicGenomicAnnotationReader() {
		super();
	}
	public BasicGenomicAnnotationReader(String annotationFile) throws FileNotFoundException {
		super();
		annotationList = new ArrayList<GenomicAnnotation>();
		source = new File(annotationFile);
		if(!source.exists()) {
			throw new FileNotFoundException("File "+ annotationFile + " does not exist");
		}	
		String [] fileExtensions = annotationFile.split("\\.");
		annotationType = fileExtensions[fileExtensions.length - 1].toUpperCase();
	}
	
	public BasicGenomicAnnotationReader(String annotationFile, String annotationType) throws FileNotFoundException{
		this(annotationFile);
		this.annotationType = annotationType; 
	}
	
	protected void setAnnotationList(List<GenomicAnnotation> list) {annotationList = list;}
	public List<GenomicAnnotation> getAnnotationList() { return annotationList;}
	
	public List<GenomicAnnotation> find(final String chr, final int startPos,final int endPos) 
	throws IOException, ParseException {
		AnnotationFactory<? extends GenomicAnnotation> factory = AnnotationFactoryFactory.getFactory(annotationType);
		final String cleanChr = chr.replace("chr", "");
		GenomicAnnotationFilter<GenomicAnnotation> filter = new GenomicAnnotationFilter<GenomicAnnotation>() {

			public boolean accept(GenomicAnnotation annotation) {
				
				return endPos >= annotation.getStart() && startPos <= annotation.getEnd() && cleanChr.equals(annotation.getChromosome());
			}

			public boolean isEnough(GenomicAnnotation annotation) {
				return false;
			}
			
		};
		
		load(source, factory, filter);
		

		return getAnnotationList();
	}
	
	public void load() throws IOException, ParseException  {
		AnnotationFactory<? extends GenomicAnnotation> factory = AnnotationFactoryFactory.getFactory(annotationType);
		load(source, factory);
	}
	
	
	public void plot(Graphics2D g2d, int start, int end, List<GenomicAnnotation> annots) {
		// FIRST plot reference
		int yStart = 100;
		int length = end - start;
		float scale = 1;
		int tickSpacing = 500;
		int ticks = (int) Math.floor(length/500);

		System.out.println("Plotting reference at of length "+ length);
		g2d.drawLine(RIGHT_MARGIN, yStart , length + RIGHT_MARGIN, yStart );		
		for(int i = 0 ; i <= ticks; i++) {
			int pos = i*tickSpacing;
			g2d.drawLine(RIGHT_MARGIN + pos, yStart - 20, RIGHT_MARGIN + pos, yStart);
			g2d.drawString(String.valueOf(pos / scale),RIGHT_MARGIN + pos - 50,yStart - 40);
		}
		
		Collections.sort(annots, new Comparator<GenomicAnnotation>() {
			public int compare(GenomicAnnotation arg0, GenomicAnnotation arg1) {				
				return (int) (Math.min(arg0.getStart(), arg0.getEnd()) - 
					Math.min(arg1.getStart(), arg1.getEnd()));
			} 
		});
		
		int yBase = yStart + 100;
		Iterator<GenomicAnnotation> annotIt = annots.iterator();
		int contigVerticalPos = yBase;
		long maxSubjectBase = 0;
		while(annotIt.hasNext()) {				
			LightweightGenomicAnnotation annot = annotIt.next();
			if(maxSubjectBase > annot.getStart() - start) {
				contigVerticalPos = contigVerticalPos +  (int)(INTER_VERTICAL_SPACING); 
			} else {
				contigVerticalPos = yBase;
			}
			//System.out.println("plotting " + annot);
			g2d.drawLine(RIGHT_MARGIN + (int) ((annot.getStart() - start) * scale), 
					contigVerticalPos, RIGHT_MARGIN + (int) ((annot.getEnd() - start) * scale), 
					contigVerticalPos);
			
			maxSubjectBase = Math.max(annot.getEnd() - start, maxSubjectBase);
		}
		
	}
	
	/**
	 * Parses an array of cannonical region descriptors of the form chrN:start-end
	 * @param rawRegions
	 * @return
	 */
	public static List<GenomicAnnotation> parseRawRegions(String[] rawRegions) {
		List<GenomicAnnotation> annotations = new ArrayList<GenomicAnnotation>();
		for(int i = 0; i < rawRegions.length; i++) {
			String regStr = rawRegions[i].trim();
			Matcher m = REG_PATTERN.matcher(regStr);
			if(m.find()) {
				String chr = m.group(2);
				Integer start = Integer.parseInt(m.group(4));
				Integer end   = Integer.parseInt(m.group(6));

				GenomicAnnotation ga  = new BasicGenomicAnnotation("annot" + i,chr,start, end);
				annotations.add(ga);
			} else {
				throw new IllegalArgumentException("region " + regStr + " is not in cannonical format");
			}
		}
		return annotations;
	}

	public static List<GenomicAnnotation> parseRawPointRegions(String[] rawRegions) {
		List<GenomicAnnotation> annotations = new ArrayList<GenomicAnnotation>();
		for(int i = 0; i < rawRegions.length; i++) {
			String regStr = rawRegions[i].trim();
			Matcher m = POINT_PATTERN.matcher(regStr);
			if(m.find()) {
				String chr = m.group(2);
				Integer base = Integer.parseInt(m.group(4));

				GenomicAnnotation ga  = new BasicGenomicAnnotation("annot" + i,chr,base, base + 1);
				annotations.add(ga);
			} else {
				throw new IllegalArgumentException("region " + regStr + " is not in cannonical format");
			}
		}
		return annotations;
	}

	@Override public GenomicAnnotation createAnnotation(GenomicAnnotation a) {
		return new BasicGenomicAnnotation(a);
	}

	@Override
	public int parse(String file,
			GenomicAnnotationFilter<GenomicAnnotation> filter,
			AnnotationHandler handler) throws ParseException, IOException {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public int parse(BufferedReader br,
			GenomicAnnotationFilter<GenomicAnnotation> filter,
			AnnotationHandler handler) throws ParseException, IOException {
		// TODO Auto-generated method stub
		return 0;
	}


}
