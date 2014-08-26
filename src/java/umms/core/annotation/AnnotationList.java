package umms.core.annotation;

import java.util.*;

import org.apache.log4j.Logger;

import broad.core.datastructures.IntervalTree.Node;
import broad.core.datastructures.IntervalTree;

import umms.core.annotation.filter.*;
import umms.core.coordinatesystem.CoordinateSpace;
import umms.core.coordinatesystem.ShuffledIterator;
import umms.core.exception.IncompleteMethodImplementationException;

import umms.core.general.Predicates;
import umms.core.general.CloseableFilterIterator;
import umms.core.readFilters.SameOrientationFilter;
import net.sf.samtools.util.CloseableIterator;
import org.apache.commons.collections15.iterators.IteratorChain;
import org.apache.commons.collections15.Predicate;


/**
 * @author engreitz
 * Class for storing and oeprating on a set of annotations in memory 
 * Allows multiple identical annotations.
 * Iterators, counting, etc. always consider the filters present in the AnnotationList.
 * TODO: Add a different class that will dynamically read BED files, perhaps by using an indexing system
 * TODO: Add a global extension factor
 * TODO: Integrate with modified reader classes that can populate an AnnotationList from a file
 *       (this class should not be responsible for parsing or reading - only for storing and manipulating)
 * TODO: Alternative (faster?) implementation:  remove all annotations that don't pass filters so you don't have to test every annotation every time you iterate.
 * 		 In this implementation, addFilter would apply the filter to all annotations upon addition, removing those that do not pass.
 * @param <T>
 */
public class AnnotationList<T extends Annotation> extends AbstractAnnotationCollection<T> implements Iterable<T> {
	
	static Logger logger = Logger.getLogger(AnnotationList.class.getName());
	
	public static boolean ENFORCE_REFERENCE = false;
	
	protected final CoordinateSpace coordinateSpace;
	protected TreeMap<String, IntervalTree<T>> annotations = new TreeMap<String, IntervalTree<T>>();
	
	// Single predicate that contains all filters.  All filters must be true for an annotation to pass.
	private Predicate<T> allFilters = Predicates.alwaysTrue();
	

	public AnnotationList() {
		this(null);
	}
	
	
	/**
	 * Initializing with a coordinate space allows for some additional filtering and sanity checks
	 * @param cs
	 */
	public AnnotationList(CoordinateSpace cs) {
		//TODO
		//logger.warn("AnnotationList has not been completed/tested yet (TODO)");
		coordinateSpace = cs;
		if (cs != null) {
			for (String reference : coordinateSpace.getChromosomeNames()) {
				annotations.put(reference, new IntervalTree<T>());
			}
		}
	}
	
	
	public AnnotationList(CoordinateSpace cs, Collection<? extends T> annotations) {
		this(cs);
		addAll(annotations);
	}
	
	
	public AnnotationList(CoordinateSpace cs, T... annotations) {
		this(cs);
		for (int i = 0; i < annotations.length; i++) {
			add(annotations[i]);
		}
	}
	

	@Override
	public int size() {
		// Considers filters
		return getIteratorCount(iterator());
	}
	
	
	@Override
	public double getGlobalCount() {
		return (double) size();
	}
	
	
	@Override
	public double getCount(Annotation region, boolean fullyContained) {
		return getIteratorCount(getOverlappingAnnotations(region, fullyContained));
	}
	

	@Override
	public double getCountExcludingRegion(Annotation region, Annotation excluded) {
		Predicate<? super T> filter = Predicates.not(new OverlapFilter(excluded));
		return getIteratorCount(new CloseableFilterIterator<T>(getOverlappingAnnotations(region, false), filter));
	}
	
	
	/**
	 * Count the number of items in the iterator (exhausts the iterator)
	 * The iterators in this class are Closeable just to match the AnnotationCollection interface - 
	 * they don't really need to be closed.
	 * @param itr
	 * @return
	 */
	private int getIteratorCount(CloseableIterator<? extends T> itr) {
		int count = getIteratorCount((Iterator) itr);
		itr.close();
		return count;
	}
	
	
	/**
	 * Count the number of items in the iterator (exhausts the iterator)
	 * @param itr
	 * @return
	 */
	private int getIteratorCount(Iterator<? extends T> itr) {
		int count = 0;
		while (itr.hasNext()) {
			itr.next();
			++count;
		}
		return count;
	}
	
	
	@Override
	public int getBasesCovered(Annotation region, boolean fullyContained) {
		Annotation collapsed = collapse(region, fullyContained);
		return (collapsed == null) ? 0 : collapsed.length();
	}
	
	
	@Override
	public CloseableIterator<T> getOverlappingAnnotations(Annotation region, boolean fullyContained) {
		if (!validAnnotation(region)) {
			return new CloseableFilterIterator<T>(new IntervalTree<T>().valueIterator(), null);
		}
		
		IntervalTree<T> tree = annotations.get(region.getReferenceName());
		Iterator<T> itr = tree.overlappingValueIterator(region.getStart(), region.getEnd());
		
		Predicate<? super T> filter = fullyContained ? new FullyContainedFilter(region) : new OverlapFilter(region);
		filter = Predicates.and(allFilters, filter);
		return new CloseableFilterIterator<T>(itr, filter);
	}
	

	public int getNumOverlappingAnnotations(Annotation region) {
		int i = 0;
		CloseableIterator<T> itr = getOverlappingAnnotations(region);
		while (itr.hasNext()) {
			itr.next();
			i++;
		}
		itr.close();
		return i;
	}
	
	
	@Override
	public CloseableIterator<T> getPermutedAnnotations(Annotation region) {
		if (coordinateSpace == null) {
			throw new IllegalStateException("Must initialize AnnotationList with a coordinate space in order to permute annotations.");
		}
		return new ShuffledIterator<T>(getOverlappingAnnotations(region, true), coordinateSpace, region);
	}
	
	
	@Override
	public CloseableIterator<T> iterator() {
		Collection<Iterator<? extends T>> iterators = new ArrayList<Iterator<? extends T>>();
		for (IntervalTree<T> tree : annotations.values()) {
			iterators.add(tree.valueIterator());
		}
		IteratorChain<T> chain = new IteratorChain<T>(iterators);
		return new CloseableFilterIterator<T>(chain, allFilters);
	}


	
	@Override
	public void addFilter(Predicate<T> filter) {
		allFilters = Predicates.and(allFilters, filter);
	}
		
	@Override
	public void addFilters(Collection<Predicate<T>> filters) {
		allFilters = Predicates.and(allFilters, Predicates.and(filters));
	}
	
	/**
	 * Add a collection of annotations
	 * @param annotations
	 */
	public void addAll(Iterable<? extends T> annotations) {
		for (T a : annotations) {
			add(a);
			//logger.info(a.getReferenceName() + ":" + this.annotations.get(a.getReferenceName()).size());
		}
	}
	
	/**
	 * Add another AnnotationList.  Duplicates allowed
	 * @param otherSet
	 */
	public void addAll(AnnotationList<? extends T> otherSet) {
		CloseableIterator<? extends T> itr = otherSet.iterator();
		while (itr.hasNext()) add(itr.next());
	}
	
	/**
	 * Add an annotation
	 * @param annotation
	 */
	public void add(T annotation) {
		if (validAnnotation(annotation)) {
			if (coordinateSpace == null && !annotations.containsKey(annotation.getReferenceName())) {
				annotations.put(annotation.getReferenceName(), new IntervalTree<T>());
			}
			annotations.get(annotation.getReferenceName()).put(annotation.getStart(), annotation.getEnd(), annotation);
		}
	}
	
	public void removeAnnotation(T annotation) {
		if (annotations.containsKey(annotation.getReferenceName())) {
			IntervalTree<T> tree = annotations.get(annotation.getReferenceName());
			tree.remove(annotation.getStart(), annotation.getEnd(), annotation);
		}
	}
	
	
	/**
	 * Make sure the annotation has an acceptable reference name
	 * @param annotation
	 */
	private boolean validAnnotation(Annotation annotation) {
		return true; // TODO Pam changed so will work with transcriptome space
		/*if (coordinateSpace != null && !annotations.containsKey(annotation.getReferenceName())) { 
			if (ENFORCE_REFERENCE) {
				throw new IllegalArgumentException("Region reference name " + annotation.getReferenceName() + " not found.");
			} else {
				logger.info("Skipping annotation " + annotation.toUCSC() + " because reference name " + annotation.getReferenceName() + " is not recognized by the AnnotationList.");
				return false;
			}
		}
		return true;*/
	}
	
	
	/**
	 *  Collapse overlapping annotations.
	 */
	public void collapse() {
		throw new UnsupportedOperationException();
	}

	
	/**
	 *  Returns the intersection of this AnnotationList with another.
	 *  Behaves like bedIntersect
	 *  Different than getOverlappers
	 * @param other
	 * @return
	 */
	public AnnotationList<T> intersect(AnnotationList<? extends Annotation> other) {
		// TODO: add getIntersection to IntervalTree?  might be able to create a faster algorithm
		throw new UnsupportedOperationException();
	}

	
	/**
	 * Returns a new AnnotationList that does not contain any annotations that overlap annotations in other
	 * @param other
	 * @return
	 */
	public AnnotationList<T> minus(AnnotationList<? extends Annotation> other) {
		AnnotationList<T> nonOverlappers = new AnnotationList<T>(coordinateSpace);
		
		CloseableIterator<T> itr = iterator();
		while (itr.hasNext()) {
			T next = itr.next();
			if (!other.hasAnnotationThatOverlaps(next)) {
				nonOverlappers.add(next);
			}
		}
		itr.close();
		return nonOverlappers;
	}


	/**
	 * @param other
	 * @return AnnotationList containing all annotations that overlap an annotation in other.
	 */
	public AnnotationList<T> getOverlappingAnnotationList(AnnotationList<? extends Annotation> other) {
		return getOverlappingAnnotationList(other, false);
	}
	
	
	/**
	 * @param region
	 * @return true if this has any Annotations that overlap region
	 */
	public boolean hasAnnotationThatOverlaps(Annotation region) {
		int count = (int) getCount(region, false);
		return (count > 0);
	}
	
	
	/**
	 * @param region
	 * @return true if this has any Annotations that fully contain region
	 */
	public boolean hasAnnotationThatContains(Annotation region) {
		CloseableIterator<T> itr = getOverlappingAnnotations(region);
		boolean result = false;
		while (itr.hasNext()) {
			if (itr.next().contains(region)) {
				result = true;
				break;
			}
		}
		return result;
	}
	
	
	/**
	 * @param other
	 * @param fullyContained
	 * @return AnnotationList containing all annotations that overlap (or contained by) an annotation in other.
	 */
	public AnnotationList<T> getOverlappingAnnotationList(AnnotationList<? extends Annotation> other, boolean fullyContained) {
		// TODO: add getOverlapping to IntervalTree?  might be able to create a faster algorithm
		AnnotationList<T> overlappers = new AnnotationList<T>(coordinateSpace);
		
		CloseableIterator<T> itr = iterator();
		while (itr.hasNext()) {
			T next = itr.next();
			if ((fullyContained && other.hasAnnotationThatContains(next)) || (!fullyContained && other.hasAnnotationThatOverlaps(next))) {
				overlappers.add(next);
			}
		}
		itr.close();
		
		return overlappers;
	}
	
	
	/**
	 * @param query
	 * @return Closest annotation in the set, or null if there are no annotations on the same chromosome
	 */
	public T getClosest(final Annotation query) {
		return getClosest(query, false);
	}
	
	
	/**
	 * @param query
	 * @return Closest non-overlapping annotation in the set, or null if there are no non-overlapping annotations on the same chromosome
	 */
	public T getClosestNonOverlapping(final Annotation query) {
		return getClosest(query, true);
	}
	
	
	/**
	 * @param query
	 * @return Closest (non-overlapping) annotation in the set, or null if there are no matches on the same chromosome
	 */
	public T getClosest(final Annotation query, boolean nonOverlapping) {
		String ref = query.getReferenceName();
		if (!annotations.containsKey(ref) || annotations.get(ref).size() == 0) return null;
		T closest = null;
		
		CloseableIterator<T> itr = getOverlappingAnnotations(query);
		if (itr.hasNext() && !nonOverlapping) {
			closest = itr.next();
		} else {
			T closestAfter = getClosestHelper(query, true);
			T closestBefore = getClosestHelper(query, false);
			int distToAfter = closestAfter == null ? Integer.MAX_VALUE : query.getDistanceTo(closestAfter);
			int distToBefore = closestBefore == null ? Integer.MAX_VALUE : query.getDistanceTo(closestBefore);
			if (closestAfter != null && distToAfter < distToBefore) {
				closest = closestAfter;
			} else if (closestBefore != null) {
				closest = closestBefore;
			}
		}
	
		return closest;
	}
	
	
	private T getClosestHelper(final Annotation query, boolean after) {
		IntervalTree<T> tree = annotations.get(query.getReferenceName());
		Node<T> closest = null;
		T result = null;
		if (after) {
			closest = tree.min(query.getEnd(), query.getEnd()+1);
		} else {
			closest = tree.max(query.getStart()-1, query.getStart());
		}
		if (closest != null) result = closest.getValue();
		return result;
	}
	
	
	/**
	 * Returns the closest non-overlapping annotation 5-prime of the query, considering strand.
	 * @param query
	 * @return
	 */
	public T getClosestUpstream(final Annotation query) {
		throw new UnsupportedOperationException();
	}
	
	
	/**
	 * Returns the closest non-overlapping annotation 3-prime of the query, considering strand.
	 * @param query
	 * @return
	 */
	public T getClosestDownstream(final Annotation query) {
		throw new UnsupportedOperationException();
	}
	
	
	/**
	 * @param query
	 * @return distance in coordinate space between the query and the closest element in the set, or null if there are no others on the same chromosome
	 */
	public Integer getDistanceToClosest(final Annotation query) {
		T closest = getClosest(query);
		return query.getDistanceTo(closest);
	}
	
	
	/**
	 * @return A list of all annotations passing the filter
	 */
	public List<T> toList() {
		List<T> list = new ArrayList<T>();
		CloseableIterator<T> itr = iterator();
		while (itr.hasNext()) list.add(itr.next());
		return list;
	}
	
	
	public CoordinateSpace getCoordinateSpace() {
		return coordinateSpace;
	}

	@Override
	public double getRefSequenceLambda(String refname) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double getCountStrandedExcludingRegion(Annotation region,
			Annotation excluded) {
		try {
			throw new IncompleteMethodImplementationException();
		} catch (IncompleteMethodImplementationException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return 0;
	}
}
