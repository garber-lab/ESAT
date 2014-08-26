package broad.core.annotation;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;

import java.util.Iterator;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.LinkedHashSet;

import broad.core.error.ParseException;

public class RepeatMaskerAnnotationReader extends AnnotationReader<RepeatMaskerAnnotation>{
	private enum RepeatTaxon { NAME, CLASS, FAMILY };
	Map<RepeatTaxon, Set<String>> uniqueNames;
	
	public RepeatMaskerAnnotationReader() {
		super();
		uniqueNames = new HashMap<RepeatTaxon, Set<String>>();
	}
	
	public RepeatMaskerAnnotationReader(BufferedReader br) throws ParseException, IOException {
		this();
		load(br, AnnotationFactoryFactory.repeatMaskerAnnotationFactory);
	}
	
	public RepeatMaskerAnnotationReader(String input) throws ParseException, IOException {
		this();
		load(new File(input), AnnotationFactoryFactory.repeatMaskerAnnotationFactory);
	}
	
	public void load(BufferedReader br) throws IOException, ParseException {
		super.load(br, AnnotationFactoryFactory.repeatMaskerAnnotationFactory);
	}
	
	public int parse(BufferedReader br, GenomicAnnotationFilter<RepeatMaskerAnnotation> filter, AnnotationHandler handler) throws ParseException {
		return super.parse(br, AnnotationFactoryFactory.repeatMaskerAnnotationFactory, filter, handler);
	}

	public RepeatMaskerAnnotation createAnnotation(GenomicAnnotation a) {
		return new RepeatMaskerAnnotation(a);
	}

	public int parse(String file, GenomicAnnotationFilter<RepeatMaskerAnnotation> filter,
			AnnotationHandler handler) throws ParseException, IOException {
		return super.parse(file,AnnotationFactoryFactory.repeatMaskerAnnotationFactory, filter, handler);
	}

	public Set<String> getUniqueRepeatNames() throws IOException {
		if (this.size() > 0 && !uniqueNames.containsKey(RepeatTaxon.NAME)) 
			collectUniqueNames();
		return uniqueNames.get(RepeatTaxon.NAME);
	}
	
	public Set<String> getUniqueRepeatClasses() throws IOException {
		if (this.size() > 0 && !uniqueNames.containsKey(RepeatTaxon.CLASS)) 
			collectUniqueNames();
		return uniqueNames.get(RepeatTaxon.CLASS);
	}
	
	public Set<String> getUniqueRepeatFamilies() throws IOException {
		if (this.size() > 0 && !uniqueNames.containsKey(RepeatTaxon.FAMILY)) 
			collectUniqueNames();
		return uniqueNames.get(RepeatTaxon.FAMILY);
	}
	
	private void collectUniqueNames() {
		LinkedHashSet<String> names = new LinkedHashSet<String>();
		LinkedHashSet<String> classes = new LinkedHashSet<String>();
		LinkedHashSet<String> families = new LinkedHashSet<String>();
		
		Iterator<RepeatMaskerAnnotation> itr = getAnnotationList().iterator();
		while (itr.hasNext()) {
			RepeatMaskerAnnotation curr = itr.next();
			names.add(curr.getRepeatName());
			classes.add(curr.getRepeatClass()); 
			families.add(curr.getRepeatFamily());
		}
		uniqueNames.put(RepeatTaxon.NAME, names);
		uniqueNames.put(RepeatTaxon.CLASS, classes);
		uniqueNames.put(RepeatTaxon.FAMILY, families);
	}
	
}
