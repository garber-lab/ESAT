package umms.core.general;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import broad.core.error.ParseException;

import org.apache.commons.collections15.Predicate;
import org.apache.commons.io.LineIterator;
import net.sf.samtools.util.CloseableIterator;

public class TabbedReader {

	/**
	 * Read an Annotation file line by line in a buffered format and return an iterator over the results.
	 * Should be refactored with the "load" function as appropriate ("load" can call this function and then
	 * add it to its Annotation list)
	 * @param file
	 * @param clazz
	 * @param factory
	 * @param cs
	 * @param filter
	 * @param skipRows number of rows to skip at the beginning of the file
	 * @return
	 */
	public static <T> CloseableIterator<T> read(File file, Class<T> clazz, TabbedReader.Factory<? extends T> factory, Predicate<? super T> filter, int skipRows) throws IOException {
		return new CloseableFilterIterator<T>(new TabbedIterator<T>(file, factory, skipRows), filter);
	}
	
	public static <T> CloseableIterator<T> read(File file, Class<T> clazz, TabbedReader.Factory<? extends T> factory, Predicate<? super T> filter) throws IOException {
		return read(file, clazz, factory, filter, 0);
	}
	
	public static <T> CloseableIterator<T> read(File file, Class<T> clazz, TabbedReader.Factory<? extends T> factory, int skipRows) throws IOException {
		return read(file, clazz, factory, Predicates.alwaysTrue(), skipRows);
	}
	
	public static <T> CloseableIterator<T> read(File file, Class<T> clazz, TabbedReader.Factory<? extends T> factory) throws IOException {
		return read(file, clazz, factory, Predicates.alwaysTrue());
	}
	
	public static class TabbedIterator<T> implements CloseableIterator<T> {
		protected LineIterator itr;
		private T curr;
		protected Factory<? extends T> factory;
		BufferedReader br;
		
		public TabbedIterator(File file, Factory<? extends T> factory, int skipRows) throws IOException {
			br = new BufferedReader(new FileReader(file));
			itr = new LineIterator(br);
			for (int i = 0; i < skipRows; i++) itr.next();
			this.factory = factory;
			advance();
		}
		
		@Override
		public void close() {
			try {
				br.close();
			} catch (Exception e) {
				throw new IllegalStateException("Could not close BufferedReader");
			}
		}
		
		@Override
		public boolean hasNext() {
			return (curr != null);
		}

		@Override
		public T next() {
			T result = curr;
			advance();
			return result;
		}

		private void advance() {
			curr = null;
			String nextLine = getNextLine();
			if (nextLine != null) {
				nextLine.trim();
				curr = factory.create(nextLine.split("\t"));
			}
		}
		
		/**
		 * Override this function to do any custom parsing (e.g., skip commented lines)
		 * @return
		 */
		protected String getNextLine() {
			if (itr.hasNext()) return itr.next();
			else return null;
		}
		
		@Override
		public void remove() {
			throw new UnsupportedOperationException("Remove not supported");
		}
		
	}
	
	
	public interface Factory<T> {
		T create(String[] rawFields) throws ParseException;
	}
}
