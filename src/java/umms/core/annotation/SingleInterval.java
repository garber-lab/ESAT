package umms.core.annotation;

import com.sleepycat.persist.model.Persistent;

import broad.core.datastructures.Interval;

/**
 * @author engreitz
 * This class represents a single interval on a integer coordinate space.
 */
@Persistent
public final class SingleInterval extends Interval.Impl implements Comparable<SingleInterval>, java.io.Serializable {
	
	//TODO: Mock Constructor
	public SingleInterval(){
		super(0,0);
	}
	
	public SingleInterval(int start, int end) {
		super(start, end);
	}
	
	public SingleInterval(Interval other) {
		super(other.getStart(), other.getEnd());
	}

	public boolean overlaps(SingleInterval other) {
		int relationship = super.getRelationship(other);
		return ((relationship & Interval.HAS_OVERLAPPING_PART) == Interval.HAS_OVERLAPPING_PART);
	}
	
	public boolean contains(SingleInterval other) {
		return (getStart() <= other.getStart() && getEnd() >= other.getEnd());
	}
	
	
	/**
	 * @param other
	 * @return 0 if overlapping, otherwise the minimum absolute distance from one to the other.
	 */
	public int getDistanceTo(SingleInterval other) {
		int result = 0;
		if (!overlaps(other)) {
			if (this.getStart() < other.getStart()) {
				result = other.getStart() - this.getEnd();
			} else {
				result = this.getStart() - other.getEnd();
			}
		}
		return result;
	}
	
	
	/**
	 * @param other
	 * @return New interval of intersection; if they don't intersect, then returns null
	 */
	public SingleInterval intersect(SingleInterval other) {
		if (!overlaps(other))
			return null;
		else
			return new SingleInterval(Math.max(this.getStart(), other.getStart()), Math.min(this.getEnd(), other.getEnd()));
	}
	
	/**
	 * @param other
	 * @return New interval spanning both intervals, including the space between them if they don't overlap
	 */
	public SingleInterval union(SingleInterval other) {
		return new SingleInterval(Math.min(this.getStart(), other.getStart()), Math.max(this.getEnd(), other.getEnd()));
	}
	
	
	/* (non-Javadoc)
	 * @see java.lang.Comparable#compareTo(java.lang.Object)
	 * Compares by start and end coordinates.
	 */
	public int compareTo(SingleInterval other) {
		if (getStart() < other.getStart()) {
			return -1;
		} else if (getStart() > other.getStart()) {
			return 1;
		} else {
			return new Integer(other.getEnd()).compareTo(new Integer(getEnd()));
		}
	}
	
	public boolean equals(SingleInterval other) {
		return (compareTo(other) == 0);
	}
	
	public String toString() {
		return "[" + getStart() + "," + getEnd() + ")"; //Intervals are open/closed by default
	}
}
