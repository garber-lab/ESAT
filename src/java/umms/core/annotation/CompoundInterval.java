package umms.core.annotation;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.TreeMap;

import umms.core.coordinatesystem.TranscriptomeSpace;

import org.apache.log4j.Logger;

import com.sleepycat.persist.model.Persistent;

/**
 * @author engreitz
 * This class contains multiple non-overlapping intervals on an arbitrary integer coordinate space.
 * Intervals are stored using red-black trees so O(log n) performance is guaranteed.
 */
@Persistent
public class CompoundInterval implements Comparable<CompoundInterval>, java.io.Serializable {

	static Logger logger = Logger.getLogger(CompoundInterval.class.getName());
	TreeSet<SingleInterval> blockTree= new TreeSet<SingleInterval>();

	// Caches for mapping between absolute and relative coordinates
	TreeMap<SingleInterval,SingleInterval> absoluteToRelative;  
	TreeMap<SingleInterval,SingleInterval> relativeToAbsolute;
	
	// Boolean to keep track of whether the caches should be recalculated or not
	boolean modified = false; 
	boolean startEndCalculated = false;
	private int start, end;	
	
	public CompoundInterval() {}
	
	public CompoundInterval(int start, int end) {
		addInterval(start, end);
	}
	
	public CompoundInterval(CompoundInterval other) {
		for (SingleInterval interval : other.getBlocks()) {
			blockTree.add(interval);  // do not need to make copies of SingleInterval because its variables are immutable
		}
		modified = true;
	}
	
	public int getStart() {
		if (blockTree.size() == 0) throw new IllegalStateException("CompoundInterval does not have any blocks.");
		if (!startEndCalculated) calculateStartEnd();
		return this.start;
	}
	
	
	public int getEnd() {
		if (blockTree.size() == 0) throw new IllegalStateException("CompoundInterval does not have any blocks.");
		if (!startEndCalculated) calculateStartEnd();
		return this.end;
	}
	
	
	private void calculateStartEnd() {
		this.start = blockTree.first().getStart();
		this.end = blockTree.last().getEnd();
		startEndCalculated = true;
	}
	
	public int numBlocks() {
		return blockTree.size();
	}
	
	public int length() {
		int length = 0;
		for (SingleInterval block : blockTree) {
			length += block.getLength();
		}
		return length;
	}
	
	public int getSpan() {
		return blockTree.last().getEnd() - blockTree.first().getStart();
	}
	
	
	private void generateAbsoluteToRelativeCache() {
		int size = 0;
		absoluteToRelative = new TreeMap<SingleInterval,SingleInterval>();  
		relativeToAbsolute = new TreeMap<SingleInterval,SingleInterval>();	
		for (SingleInterval absoluteBlock : blockTree) {
			SingleInterval relativeBlock = new SingleInterval(size, size + absoluteBlock.getLength());
			size += absoluteBlock.getLength();
			absoluteToRelative.put(absoluteBlock, relativeBlock);
			relativeToAbsolute.put(relativeBlock, absoluteBlock);
		}
		modified = false;
	}
	
	
	/**
	 * @param positionInInterval 0-based position in the compound interval
	 * @return reference coordinate at the given interval position.  e.g. for positionInInterval=0, returns getStart()
	 */
	public int getCoordinateAtPosition(int positionInInterval) {
		if (positionInInterval < 0) throw new IllegalArgumentException("PositionInInterval must be non-negative.");
		if (positionInInterval >length()) throw new IllegalArgumentException("PositionInInterval is longer than the entire interval.");
		
		if (modified) generateAbsoluteToRelativeCache();
		
		SingleInterval query = new SingleInterval(positionInInterval, positionInInterval + 1);
		SingleInterval relativeBlock = relativeToAbsolute.floorKey(query);
		
		return relativeToAbsolute.get(relativeBlock).getStart() + positionInInterval - relativeBlock.getStart();
		/*  uncached implementation
		Iterator<SingleInterval> itr = blockTree.iterator();
		while (itr.hasNext()) {
			SingleInterval curr = itr.next();
			if (curr.getLength() >= positionInInterval) {
				return curr.getStart() + positionInInterval;
			} else {
				positionInInterval -= curr.getLength();
			}
		}
		
		throw new IllegalStateException("should never reach this point");
		*/
	}
	
	
	/**
	 * @param coordinateInReference
	 * @return 0-based interval position at given reference coordinate
	 */
	public int getPositionAtCoordinate(int coordinateInReference) {
		if (coordinateInReference < getStart() | coordinateInReference >= getEnd()) {
			throw new IllegalArgumentException("coordinateInReference outside bounds of interval");
		}
		
		if (modified) generateAbsoluteToRelativeCache();
		
		SingleInterval query = new SingleInterval(coordinateInReference, coordinateInReference + 1);
		SingleInterval absoluteBlock = absoluteToRelative.floorKey(query);
		
		return absoluteToRelative.get(absoluteBlock).getStart() + coordinateInReference - absoluteBlock.getStart();
		
		/* uncached implementation
		SingleInterval query = new SingleInterval(coordinateInReference, coordinateInReference+1);
		SingleInterval potentialOverlapper = blockTree.floor(query);
		if (potentialOverlapper != null && query.overlaps(potentialOverlapper)) {
			return coordinateInReference - potentialOverlapper.getStart();
		}
		return position;
		*/
	}
	
	
	/**
	 * Extends or trims the start coordinate.  If trimming, removes any blocks that
	 * fall before the new start position and truncates the block that overlaps the
	 * start position, if it exists.  If expanding, extends the first block to the new start position.
	 * @param start new start coordinate
	 */
	public void setStart(int start) {
		if (start > getEnd()) throw new IllegalArgumentException("Interval start cannot be greater than end.");
		SingleInterval newInterval = null;
		if (start <= getStart()) {
			newInterval = new SingleInterval(start, blockTree.first().getEnd());
			blockTree.remove(blockTree.first());
		} else {
			Iterator<SingleInterval> itr = blockTree.iterator();
			while (itr.hasNext()) {
				SingleInterval curr = itr.next();
				if (start >= curr.getEnd()) {
					// interval is entirely after the new start point - remove it.
					itr.remove();
				} else {
					// extend the new first interval to the new start point
					newInterval = new SingleInterval(start, curr.getEnd());
					itr.remove();
					break;  // we've passed all relevant intervals
				}
			}
		}
		blockTree.add(newInterval);
		modified = true;
		startEndCalculated = false;
	}
	
	
	/**
	 * Extends or trims the end coordinate. If trimming, removes blocks that fall
	 * after the new end coordinate.  If expanding, extends the last block to the
	 * new end coordinate.
	 * @param end new end coordinate
	 */
	public void setEnd(int end) {
		if (end < getStart()) throw new IllegalArgumentException("Interval end cannot be less than its start.");
		SingleInterval newInterval = null;
		if (end >= getEnd()) {
			// just extend the last interval
			newInterval = new SingleInterval(blockTree.last().getStart(), end);
			blockTree.remove(blockTree.last());
		} else {
			Iterator<SingleInterval> itr = blockTree.descendingIterator();
			while (itr.hasNext()) {
				SingleInterval curr = itr.next();
				if (end <= curr.getStart()) { //MG: changed <= to < since our intervals are closed/open so if end is equal to start the interval [s,s) is valid, it is a single point.
					// interval is entirely after the new end point - remove it.
					logger.debug("Trying to set a short end "+end+" to interval which will remove the current block: " + curr);
					itr.remove();
				} else {
					// extend the new last interval to the new end point
					newInterval = new SingleInterval(curr.getStart(), end);  
					itr.remove();
					break;  // we've passed all relevant intervals
				}
			}
		}
		
		blockTree.add(newInterval);
		modified = true;
		startEndCalculated = false;
	}
	
	
	/**
	 * Adds an interval, updating the compound interval to reflect the intersection of its
	 * previous state and the new interval.
	 * @author engreitz
	 * @param newInterval
	 */
	public void addInterval(SingleInterval newInterval) {
		SingleInterval start = blockTree.floor(newInterval);
		if (start == null) start = newInterval;
		
		Iterator<SingleInterval> itr = blockTree.tailSet(start, true).iterator();
		
		while (itr.hasNext()) {
			SingleInterval i = itr.next();
			if (i.getStart() > newInterval.getEnd()) break;  // we've passed all relevant intervals
			
			if (i.isAdjacent(newInterval) || i.overlaps(newInterval)) {	
				newInterval = i.union(newInterval);
				itr.remove();
			}
		}
		
		blockTree.add(newInterval);
		modified = true;
		startEndCalculated = false;
	}
	
	
	public void addInterval(int start, int end) {
		addInterval(new SingleInterval(start, end));
	}

	/**
	 * Shift all blocks in the compound interval by a given delta
	 * @param delta
	 */
	public void shift(int delta) {
		TreeSet<SingleInterval> newTree = new TreeSet<SingleInterval>();
		Iterator<SingleInterval> itr = blockTree.iterator();
		while (itr.hasNext()) {
			SingleInterval curr = itr.next();
			newTree.add(new SingleInterval(curr.getStart() + delta, curr.getEnd() + delta));
		}
		blockTree = newTree;
		startEndCalculated = false;
		modified = true;
	}
	
	/**
	 * Move an annotation, preserving the relationships between its blocks, to a new coordinate.
	 * @param coordinateInReference
	 */
	public void moveToCoordinate(int coordinateInReference) {
		int delta = coordinateInReference - getStart();
		shift(delta);
	}
	
	/**
	 * @param i
	 * @return true if an interval with these exact boundaries is present in the CompoundInterval
	 */
	public boolean containsExactInterval(SingleInterval i) {
		return blockTree.contains(i);
	}
	
	/**
	 * @param i
	 * @return true if i is entirely contained within a block (i can be smaller than the block)
	 */
	public boolean containsInterval(SingleInterval i) {
		SingleInterval floor = blockTree.floor(i);
		return ((floor != null) && floor.contains(i));
	}
	
	/**
	 * @param other
	 * @return returns true if all blocks in "other" are contained in blocks of this object
	 */
	public boolean contains(CompoundInterval other) {
		for (SingleInterval block : other.getBlocks()) {
			if (!containsInterval(block)) return false;
		}
		return true;
	}
	
	/**
	 * @param other
	 * @return returns true if any block in "other" overlaps with any block of this object
	 */
	public boolean overlaps(CompoundInterval other) {
		return overlaps(other, false);
	}
		
	/**
	 * @param other
	 * @param ignoreBlocks if true, then overlap will consider only the start and end boundaries, ignoring blocks
	 * @return true if the two intervals overlap
	 */
	public boolean overlaps(CompoundInterval other, int buffer, boolean ignoreBlocks) {
		if (!ignoreBlocks) {
			for (SingleInterval interval : other.getBlocks()) {
				if (overlaps(interval, buffer)) return true;
			}
			return false;
		} else {
			return ((getStart() < other.getEnd() + buffer) && (getEnd() > other.getStart() - buffer));
		}
	}
	
	
	/**
	 * @param other
	 * @param ignoreBlocks
	 * @return true if the boundaries of the intervals overlap
	 */
	public boolean overlaps(CompoundInterval other, boolean ignoreBlocks) {
		return overlaps(other, 0, ignoreBlocks);
	}
	
	
	/**
	 * @param other
	 * @param buffer
	 * @return true if any of the blocks come within the buffer distance of each other
	 */
	public boolean overlaps(CompoundInterval other, int buffer) {
		return overlaps(other, buffer, false);
	}
	
	
	/**
	 * @param interval
	 * @return returns true if "interval" overlaps with any block of this object
	 */
	public boolean overlaps(SingleInterval interval) {
		return (this.intersect(interval).size() > 0) ? true : false;
	}
	
	
	/**
	 * @param interval
	 * @param buffer
	 * @return true if any block comes within the buffer distance from the given interval
	 */
	public boolean overlaps(SingleInterval interval, int buffer) {
		SingleInterval expanded = new SingleInterval(interval.getStart() - buffer, interval.getEnd() + buffer);
		return overlaps(expanded);
	}
	
	
	/**
	 * Remove an interval.  Throws exception if the interval does not exist - so check first.
	 * @param i
	 */
	public void removeInterval(SingleInterval i) {
		if (!containsExactInterval(i)) 
			throw new IllegalStateException("Tried to remove interval that didn't exist");
		else
			blockTree.remove(i);
		
		modified = true;
		startEndCalculated = false;
	}
	

	/**
	 * @param other
	 * @return new CompoundInterval containing intersections of the blocks of the two objects
	 */
	public CompoundInterval intersect(CompoundInterval other) {
		CompoundInterval result = new CompoundInterval();
		for (SingleInterval block : other.getBlocks()) {
			List<SingleInterval> intersections = intersect(block);
			for (SingleInterval intersection : intersections) {
				result.addInterval(intersection);
			}
		}
		return result;
	}
	

	/**
	 * @param other
	 * @return new CompoundInterval containing the (blocked) union of the two objects
	 */
	public CompoundInterval union(CompoundInterval other) {
		CompoundInterval i = new CompoundInterval(this);
		for (SingleInterval newInterval : other.getBlocks()) {
			i.addInterval(newInterval);
		}
		return i;
	}
	
	
	/**
	 * @return a new CompoundInterval containing the gaps between the blocks in this interval
	 */
	public CompoundInterval complement() {
		return complement(getStart(), getEnd());
	}
	
	/**
	 * @return complement of this CompoundInterval using the provided bounds
	 */
	public CompoundInterval complement(int start, int end) {
		CompoundInterval bounds = new CompoundInterval(start, end);
		if (this.numBlocks() == 0) {
			return bounds;
		}
		
		CompoundInterval result = new CompoundInterval();
		
		if (start < getStart()) result.addInterval(start, getStart());
		
		Iterator<SingleInterval> itr = blockTree.iterator();
		SingleInterval interval2 = itr.next();
		
		while (itr.hasNext()) {
			SingleInterval interval1 = interval2;
			interval2 = itr.next();
			result.addInterval(interval1.getEnd(), interval2.getStart());
		}
		
		if (end > getEnd()) result.addInterval(getEnd(), end);
		
		// Trim if necessary
		result = result.intersect(bounds);
		
		return result;
	}
	
	
	public CompoundInterval minus(CompoundInterval other) {
		return this.intersect(other.complement(this.getStart(), this.getEnd()));
	}
	
	/**
	 * @return individual SingleIntervals.  should not be modified
	 */
	public final SortedSet<SingleInterval> getBlocks() { return blockTree; }

	
	private List<SingleInterval> intersect(SingleInterval interval) {
		List<SingleInterval> intersections = new ArrayList<SingleInterval>();
		
		SingleInterval start = blockTree.floor(interval);
		if (start == null) start = interval;
		
		SortedSet<SingleInterval> subset = blockTree.tailSet(start, true);
		for (SingleInterval block : subset) {
			if (block.getStart() >= interval.getEnd()) break;  // we've passed all relevant intervals
			SingleInterval intersection = interval.intersect(block);
			if (intersection != null) intersections.add(intersection);
		}
		return intersections;
	}
	
	
	/**
	 * @param other
	 * @return true if all contained intervals are equal
	 */
	public boolean equals(CompoundInterval other) {
		return blockTree.equals(other.getBlocks());
	}
	
	
	public int compareTo(CompoundInterval other) {
		if (this.equals(other)) return 0;
		if (getStart() == other.getStart()) {
			return getEnd() - other.getEnd();
		} else {
			return getStart() - other.getStart();
		}
	}
	
	
	// For testing purposes
	public String toString() {
		StringBuilder b = new StringBuilder("CompoundInterval: ");
		for (SingleInterval i : blockTree) {
			b.append(" ");
			b.append(i);
		}
		return b.toString();
	}
}
