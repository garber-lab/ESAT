package umms.core.annotation;

import java.io.*;
import java.util.Collection;

import umms.core.coordinatesystem.CoordinateSpace;
import umms.core.general.TabbedReader;
import umms.core.general.Predicates;
import org.apache.commons.collections15.Predicate;
import umms.core.general.CloseableFilterIterator;
import net.sf.samtools.util.CloseableIterator;

/**
 * @author engreitz
 * This class is responsible for reading annotation files and generating AnnotationLists from them.
 * Parsing logic is partially copied from broad.core.annotation.AnnotationReader
 */
public class AnnotationFileReader {
	
	public static <T extends Annotation> AnnotationList<T> load(File file, Class<T> clazz, TabbedReader.Factory<? extends T> factory) throws IOException {
		return load(file, clazz, factory, null);
	}
	

	public static <T extends Annotation> AnnotationList<T> load(File file, Class<T> clazz, TabbedReader.Factory<? extends T> factory, CoordinateSpace cs)  throws IOException {
		return load(file, clazz, factory, cs, Predicates.alwaysTrue());	
	}
	
	public static <T extends Annotation> AnnotationList<T> load(File file, Class<T> clazz, TabbedReader.Factory<? extends T> factory, CoordinateSpace cs, Collection<Predicate<? super T>> filters)  throws IOException {
		return load(file, clazz, factory, cs, Predicates.and(filters));
	}

	
	/**
	 * Reads an Annotation file line by line and return an AnnotationList containing the results
	 * @param <T>
	 * @param file File containing annotations
	 * @param clazz	The class of the desired AnnotationList
	 * @param cs CoordinateSpace with which to initialize the AnnotationList. Can be null.
	 * @param factory Factory for parsing and creating the annotation type
	 * @param filter Annotation filters to control the subset of annotations that are stored in the AnnotationList.
	 * @return AnnotationList containing all annotations from file that pass the filter
	 * @throws IOException 
	 */
	public static <T extends Annotation> AnnotationList<T> load(File file, Class<T> clazz, TabbedReader.Factory<? extends T> factory, CoordinateSpace cs, Predicate<? super T> filter) throws IOException {
		AnnotationList<T> annotations = new AnnotationList<T>(cs);
		CloseableIterator<T> itr = read(file, clazz, factory, filter);
		while (itr.hasNext()) {
			annotations.add(itr.next());
		}		
		itr.close();
		return annotations;
	}
	
	
	/**
	 * Read an Annotation file line by line in a buffered format and return an iterator over the results.
	 * Should be refactored with the "load" function as appropriate ("load" can call this function and then
	 * add it to its Annotation list)
	 * @param file
	 * @param clazz
	 * @param factory
	 * @param cs
	 * @param filter
	 * @return
	 */
	public static <T extends Annotation> CloseableIterator<T> read(File file, Class<T> clazz, TabbedReader.Factory<? extends T> factory, Predicate<? super T> filter) throws IOException {
		return new CloseableFilterIterator<T>(new AnnotationIterator<T>(file, factory), filter);
	}
	
	
	public static <T extends Annotation> CloseableIterator<T> read(File file, Class<T> clazz, TabbedReader.Factory<? extends T> factory) throws IOException  {
		return read(file, clazz, factory, Predicates.alwaysTrue());
	}
	
	
	private static class AnnotationIterator<T> extends umms.core.general.TabbedReader.TabbedIterator<T> {
		
		public AnnotationIterator(File file, TabbedReader.Factory<? extends T> factory) throws IOException {
			super(file, factory, 0);
		}

		@Override
		protected String getNextLine() {
			while (itr.hasNext()) {
				String line = itr.next();

				if (line.toLowerCase().startsWith("track")) {
					throw new IllegalArgumentException("AnnotationFileReader does not support files with track headers (TODO)");
				}
				
				line = line.trim();
				if (line.startsWith("#") || line.length() == 0) {
					continue;
				} else {
					return line;
				}
			}
			return null;
		}
	}
	
}
