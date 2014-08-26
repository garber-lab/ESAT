package umms.core.general;

import java.util.Arrays;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import org.apache.commons.collections15.Predicate;

/**
 * @author engreitz
 * Inspired by Google Guava Predicates class
 */
public final class Predicates {
	
	
	public static <T> Predicate<T> alwaysTrue() {
		return new Predicate<T>() { public boolean evaluate(T t) { return true; } } ;
	}
	
	
	public static <T> Predicate<T> alwaysFalse() {
		return new Predicate<T>() { public boolean evaluate(T t) { return false; } } ;
	}
	
	
	/**
	 * Returns a predicate that evaluates to {@code true} if the object being tested
	 * passes all component predicates.  If {@components} is empty, the predicate will
	 * always return true.  It defensively copies the iterable passed in, so future
	 * changes to the iterable will not alter the behavior of this predicate.
	 * @param <T>
	 * @param components Collection of predicates to evaluate
	 * @return
	 */
	public static <T> Predicate<T> and(Iterable<? extends Predicate<? super T>> components) {
		return new AndPredicate<T>(defensiveCopy(components));
	}
	
	
	public static <T> Predicate<T> and(Predicate<? super T>... components) {
		return new AndPredicate<T>(defensiveCopy(components));
	}
	
	public static <T> Predicate<T> not(Predicate<T> predicate) {
		return new NotPredicate<T>(predicate);
	}
	
	public static <T> Predicate<T> or(Iterable<? extends Predicate<? super T>> components) {
		return new OrPredicate<T>(defensiveCopy(components));
	}
	
	
	public static <T> Predicate<T> or(Predicate<? super T>... components) {
		return new OrPredicate<T>(defensiveCopy(components));
	}

	
	private static class AndPredicate<T> implements Predicate<T> {
	    private final Iterable<? extends Predicate<? super T>> components;

	    private AndPredicate(Iterable<? extends Predicate<? super T>> components) {
	      this.components = components;
	    }
	    
	    @Override
	    public boolean evaluate(T t) {
	    	Iterator<? extends Predicate<? super T>> itr = components.iterator();
	    	while (itr.hasNext()) {
	    		if (!itr.next().evaluate(t)) return false; 
	    	}
	    	return true;
	    }
	}
	
	private static class OrPredicate<T> implements Predicate<T> {
		private final Iterable<? extends Predicate<? super T>> components;
		private OrPredicate(Iterable<? extends Predicate<? super T>> components) {
			this.components = components;
		}
		@Override
		public boolean evaluate(T t) {
			Iterator<? extends Predicate<? super T>> itr = components.iterator();
			while (itr.hasNext()) {
				if (itr.next().evaluate(t)) return true;
			}
			return false;
		}
	}
	
	private static class NotPredicate<T> implements Predicate<T> {
		Predicate<T> p;
		private NotPredicate(Predicate<T> p) {
			this.p = p;
		}
		@Override
		public boolean evaluate(T t) {
			return !p.evaluate(t);
		}
	}
	
	private static <T> List<T> defensiveCopy(T... array) {
		return defensiveCopy(Arrays.asList(array));
	}
	
	private static <T> List<T> defensiveCopy(Iterable<T> iterable) {
		ArrayList<T> list = new ArrayList<T>();
		for (T element : iterable) {
			if (element == null) {
				throw new IllegalArgumentException("Null predicate passed");
			}
			list.add(element);
		}
		return list;
	}
}
