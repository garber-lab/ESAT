package broad.core.datastructures;

import java.util.Iterator;
import java.util.ListIterator;

public class ReversibleIterator<E> implements Iterator<E> {
	private ListIterator<E> listIterator;
	private boolean traverseForward;

	public ReversibleIterator (ListIterator<E> listIt, boolean traverserseForward) {
		this.listIterator = listIt;
		this.traverseForward = traverserseForward;
	}
	public boolean hasNext() {
		return traverseForward ? listIterator.hasNext() : listIterator.hasPrevious();
	}

	public void remove() {
		listIterator.remove();
	}
	
	public E next() {
		return traverseForward ? listIterator.next() : listIterator.previous();
	}

}
