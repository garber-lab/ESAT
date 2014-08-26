package umms.core.general;

import java.util.Iterator;
import org.apache.commons.collections15.Predicate;
import org.apache.commons.collections15.iterators.FilterIterator;
import net.sf.samtools.util.CloseableIterator;

public class CloseableFilterIterator<T> implements CloseableIterator<T> {
	FilterIterator<T> filterIterator;
	boolean closeable;
	
	public CloseableFilterIterator(CloseableIterator<T> itr, Predicate<? super T> filter) {
		filterIterator = new FilterIterator<T>(itr, filter);
		closeable = true;
	}
	
	public CloseableFilterIterator(Iterator<T> itr, Predicate<? super T> filter) {
		filterIterator = new FilterIterator<T>(itr, filter);
		closeable = false;
	}
	
	public boolean hasNext() {
		return filterIterator.hasNext();
	}
	
	public T next() {
		return filterIterator.next();
	}
	
	public void remove() {
		filterIterator.remove();
	}
	
	public void close() {
		if (closeable) ((CloseableIterator<T>) filterIterator.getIterator()).close();
	}
}
