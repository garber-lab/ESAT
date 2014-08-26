package umms.core.coordinatesystem;

import net.sf.samtools.util.CloseableIterator;
import umms.core.annotation.Annotation;

public class ShuffledIterator<T extends Annotation> implements CloseableIterator<T> {
	protected CloseableIterator<T> itr;
	protected Annotation region;
	protected CoordinateSpace cs;
	
	public ShuffledIterator(CloseableIterator<T> itr, CoordinateSpace cs, Annotation region) {
		this.itr = itr;
		this.cs = cs;
		this.region = region;
	}
	public boolean hasNext() { return itr.hasNext(); }
	public T next() {
		T next = itr.next();
		cs.permuteAnnotation(next, region);
		return next;
	}
	public void remove() { throw new UnsupportedOperationException(); }
	public void close() { itr.close(); }
}
