package umms.core.coordinatesystem;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;
import java.util.Random;

import org.apache.log4j.Logger;

import umms.core.annotation.*;
import umms.core.annotation.filter.FullyContainedFilter;
import umms.core.feature.GenomeWindow;
import umms.core.feature.Window;
import broad.core.annotation.ShortBEDReader;
import broad.pda.annotation.BEDFileParser;

import org.apache.commons.collections15.Predicate;
import org.apache.commons.collections15.iterators.FilterIterator;


/**
 * @author engreitz
 * Class representing a genome with defined linear chromosomes.  Masking behavior can be controlled
 * through setting masked regions and percent masked allowed parameters, which will filter out windows
 * that overlap masked regions.
 */
public class GenomicSpace implements CoordinateSpace {

	Map<String, Integer> chromosomeSizes;
	AnnotationList<Annotation> maskedRegions = new AnnotationList<Annotation>();
	double pctMaskedAllowed = 0;
	boolean overlapAllowed = false;
	
	static Logger logger = Logger.getLogger(GenomicSpace.class.getName());
	
	static final public int PERMUTATION_ATTEMPTS = 10;
	static private Random generator = new Random();
	

	/**
	 * 
	 * @param chromosomeSizes
	 */
	public GenomicSpace(Map<String, Integer> chromosomeSizes){
		this.chromosomeSizes=chromosomeSizes;
	}
	
	/**
	 * Chromosome name and accompanying sizes
	 * @param chromosomeSizes
	 */
	public GenomicSpace(String chromosomeSizes){
		this(BEDFileParser.loadChrSizes(chromosomeSizes));
	}
	
	public GenomicSpace(String chromosomeSizes, String maskedRegionFile) {
		this(chromosomeSizes);
		if (maskedRegionFile != null) {
			try {
				maskedRegions = AnnotationFileReader.load(new File(maskedRegionFile), Annotation.class, new BasicAnnotation.Factory());
			} catch (IOException e) {
				throw new IllegalArgumentException(maskedRegionFile + " not found: " + e.getMessage());
			}
		}
	}
	
	public GenomicSpace(Map<String, Integer> chromosomeSizes, String maskedRegionFile) {
		this(chromosomeSizes);
		if (maskedRegionFile != null) {
			try {
				maskedRegions = AnnotationFileReader.load(new File(maskedRegionFile), Annotation.class, new BasicAnnotation.Factory());
			} catch (IOException e) {
				throw new IllegalArgumentException(maskedRegionFile + " not found: " + e.getMessage());
			}
		}
	}
	
	public GenomicSpace(String chromosomeSizes, String maskedRegionFile, double pctMaskedAllowed) {
		this(chromosomeSizes, maskedRegionFile);
		setPercentMaskedAllowed(pctMaskedAllowed);
	}

	
	public void setPercentMaskedAllowed(double pct) {
		if (pct < 0 || pct > 100) {
			throw new IllegalArgumentException("pct must be in the range of [0,100]");
		}
		pctMaskedAllowed = pct;
		this.overlapAllowed = Math.abs(pct - 0.0) > 0.0001;
	}
	
	
	/**
	 * Returns an interval tree over the genomic region
	 */
	public Collection<? extends Window> getOverlappingRegion(String chr,int start,int end) {
		Collection<Window> rtrn=new TreeSet<Window>();
		/*
		 * Check that the fragment coordinates are contained within the coordinate space 
		 * If not, adjust the fragment ends. 
		 */
		if(start<0)
			start = 0;
		if(end>this.chromosomeSizes.get(chr))
			end = this.chromosomeSizes.get(chr);
		
		Window window = new GenomeWindow(chr,start,end);
		
		if (window!= null & this.chromosomeSizes.containsKey(window.getReferenceName())){
			if(window.getEnd()<this.chromosomeSizes.size())
				rtrn.add(window);
		}
		return rtrn;
	}

	/**
	 * Returns an iterator on chromosome chr starting at start, ending at end, of size, windowSize with overlap between consecutive windows 
	 */
	public Iterator<Window> getWindowIterator(int windowSize, String chr,int start, int end,int overlap) {
		return getWindowIterator(windowSize, chr, start, end, overlap, false);
	}
	
	public Iterator<Window> getWindowIterator(int windowSize, String chr,int start, int end, int overlap, boolean backward) {
		Iterator<Window> rtrn = new WindowsIterator(chr, windowSize, start, end, overlap, backward);
		return new FilterIterator<Window>(rtrn, new MaskFilter<Window>());
	}
	
	@Override
	public Iterator<? extends Window> getWindowIterator(Annotation window, int windowSize, int overlap) {
		return getWindowIterator(windowSize, window.getChr(), window.getStart(), window.getEnd(), overlap);
	}
	
	public Iterator<? extends Window> getWindowIterator(Annotation window, int windowSize, int overlap, boolean considerStrand) {
		return getWindowIterator(windowSize, window.getChr(), window.getStart(), window.getEnd(), overlap, considerStrand && window.getStrand() == Annotation.Strand.NEGATIVE);
	}

	/**
	 * Get a window iterator from start to finish on chromosome chr
	 */
	public Iterator<Window> getWindowIterator(String chr, int windowSize, int overlap) {
		return getWindowIterator(windowSize, chr,0, this.chromosomeSizes.get(chr), overlap);
	}
	
	/**
	 * Get a window iterator from start to finish in the entire space
	 */
	public Iterator<Window> getWindowIterator(int windowSize, int overlap) {
		Iterator<Window> rtrn = new ChromosomeIterator(chromosomeSizes, windowSize, overlap);
		rtrn = new FilterIterator<Window>(rtrn, new MaskFilter<Window>());
		return rtrn;
	}
	
	
	public boolean isValidWindow(Annotation window) {
		if (!chromosomeSizes.containsKey(window.getReferenceName())) return false;
		if (!getReferenceAnnotation(window.getReferenceName()).contains(window)) return false;
		if (! new MaskFilter<Annotation>().evaluate(window)) return false;
		return true;
	}

	/**
	 * Returns true if this coordinate space contains data for the specified chromosome
	 */
	public boolean containsDataFor(String chromosome){
		return this.chromosomeSizes.containsKey(chromosome);
	}
	
	
	private static class ChromosomeIterator implements Iterator<Window> {
		WindowsIterator currentIterator;
		Iterator<String> chromosomeIterator;
		Map<String, Integer> chromosomeSizes;
		int windowSize;
		int overlap;
		
		ChromosomeIterator(Map<String, Integer> chromosomeSizes, int windowSize, int overlap){
			this.chromosomeIterator=chromosomeSizes.keySet().iterator();
			this.chromosomeSizes=chromosomeSizes;
			this.windowSize=windowSize;
			this.overlap=overlap;
		}
		
		@Override
		public boolean hasNext() {
			if(currentIterator!=null && currentIterator.hasNext()){
				return true;
			}
			else{
				return chromosomeIterator.hasNext();
			}
		}

		@Override
		public Window next() {
			if(currentIterator!=null && currentIterator.hasNext()){
				return currentIterator.next();
			}
			else{
				//make new windows iterator for chromosome
				currentIterator=makeNewIterator(chromosomeIterator.next());
				return next();
			}
		}

		private WindowsIterator makeNewIterator(String chr) {
			return new WindowsIterator(chr, windowSize, 0, chromosomeSizes.get(chr), overlap);
		}

		@Override
		public void remove() {
			throw new UnsupportedOperationException();
		}
		
	}
	
	
	/**
	 * An iterator over the windows on chromosome
	 * 
	 */
	private static class WindowsIterator implements Iterator<Window> {

		private String chr;
		private int currPosition;
		private int start;
		private int windowSize;
		private int step;	
		private int end;
		private boolean backward; 
		
		/**
		 * Constructs an iterator on chromosome chr starting at start, ending at end, of size, windowSize with overlap between consecutive windows 
		 * @param chr
		 * @param windowSize
		 * @param start
		 * @param overlap
		 * @param end
		 */
		public WindowsIterator(String chr, int windowSize, int start, int end, int overlap, boolean backward) {
			if (overlap > windowSize) throw new IllegalArgumentException("Overlap cannot be greater than window size.");
			if (overlap < 0) throw new IllegalArgumentException("Overlap cannot be less than 0");
			
			this.chr = chr;
			this.windowSize = windowSize;
			this.step = windowSize - overlap;
			this.start = start;
			this.end = end;
			this.backward = backward;
			
			if (backward) {
				this.currPosition = end;
			} else {
				this.currPosition = start;
			}
		}
		
		public WindowsIterator(String chr, int windowSize, int start, int end, int overlap) {
			this(chr, windowSize, start, end, overlap, false);
		}
		
		/**
		 * Returns true if there is another window
		 */
		public boolean hasNext() {
			return (!backward && end >= (currPosition + windowSize)) || (backward && start <= (currPosition - windowSize));
		}

		/**
		 * Returns the next window
		 */
		public Window next() {
			Window w;
			if (backward) {
				w = new GenomeWindow(chr, currPosition - windowSize, currPosition);
				currPosition = currPosition - step;
			} else {
				w = new GenomeWindow(chr, currPosition, currPosition + windowSize);
				currPosition = currPosition + step;
			}
			return w;
		}

		/**
		 * N/A here
		 */
		public void remove() {
			throw new UnsupportedOperationException(); 
		}
		
		
	}


	@Override
	public long getLength() {
		long rtrn=0;
		for(String chr: this.chromosomeSizes.keySet()){
			rtrn+=this.chromosomeSizes.get(chr);
		}
		return rtrn;
	}

	@Override
	public long getLength(String chr) {
		return this.chromosomeSizes.get(chr);
	}

	public long getUnmaskedLength(String chr){
		return (this.getLength(chr) - maskedRegions.getBasesCovered(this.getReferenceAnnotation(chr)));
	}
	@Override
	public Annotation permuteAnnotation(Annotation a) {
		return permuteAnnotation(a, new BasicAnnotation(a.getChr(), 0, chromosomeSizes.get(a.getChr())));
	}
	
	public boolean hasChromosome(String chr) {
		return chromosomeSizes.containsKey(chr);
	}
		
	@Override
	public Annotation permuteAnnotation(Annotation a, Annotation bounds) {
		MaskFilter<Annotation> filter = new MaskFilter<Annotation>();
		
		int permutationSpace = bounds.size() - a.getLengthOnReference();
		Annotation newAnnotation = a.copy();
		newAnnotation.setReferenceName(bounds.getReferenceName());

		boolean found = true;
		for (int i = 0; i < PERMUTATION_ATTEMPTS; i++) {
			int newStart = generator.nextInt(permutationSpace) + bounds.getStart();
			newAnnotation.moveToCoordinate(newStart);
			if (filter.evaluate(newAnnotation)) {
				found = true;
				a.moveToCoordinate(newStart);
				break;
			}
		}
		if (!found) throw new PermutationNotFoundException(PERMUTATION_ATTEMPTS);
		
		return newAnnotation;
	}
	

	private class MaskFilter<T extends Annotation> implements Predicate<T> {
		public boolean evaluate(T w) {
			//System.out.println("MaskFilter for " + w.toUCSC());
			boolean reject = false;
			if (maskedRegions.size() > 0) {
				if (!overlapAllowed) {
					reject = (maskedRegions.getNumOverlappingAnnotations(w) > 0);
				} else {
					try {
						double basesCovered = (double) maskedRegions.getBasesCovered(w);
						reject = (basesCovered / (double) w.length())*100.0 > pctMaskedAllowed;
					} catch (Exception e) {
						throw new RuntimeException(e.getMessage());
					}
				}
			}
			return !reject;
		}
	}
	
	
	@Override
	public final Collection<String> getReferenceNames() {
		return chromosomeSizes.keySet();
	}
	
	@Override
	public final Collection<? extends Annotation> getReferenceAnnotations() {
		List<Annotation> annotations = new ArrayList<Annotation>();
		for (String chr : chromosomeSizes.keySet()) {
			annotations.add(new BasicAnnotation(chr, 0, (int) getLength(chr)));
		}
		return annotations;
	}
	
	/**
	 * Get an annotation spanning an entire chromosome
	 * @param chrName The chromosome name
	 * @return Annotation consisting of the entire chromosome
	 */
	public Annotation getEntireChromosome(String chrName) {
		if(!chromosomeSizes.containsKey(chrName)) {
			throw new IllegalArgumentException("Chromosome name " + chrName + " not recognized.");
		}
		return new BasicAnnotation(chrName, 0, chromosomeSizes.get(chrName).intValue(), chrName);
	}
	
	@Override
	public Annotation getReferenceAnnotation(String name) {
		if(!chromosomeSizes.containsKey(name)) {
			throw new IllegalArgumentException("Chromosome name " + name + " not recognized.");
		}
		return new BasicAnnotation(name, 0, (int) getLength(name));
	}

	@Override
	public Iterator<Window> getWindowIterator(Collection<Gene> baseGenes, int windowSize, int overlap) {
		throw new UnsupportedOperationException("TODO");
	}

	/**
	 * Set the seed for the random generator (used for permuting regions).
	 */
	public static void setSeed(long seed) {
		generator.setSeed(seed);
	}


	/**
	 * Returns a fragment of type Window between start and end on chromosome
	 * @param 
	 */
	public Collection<? extends Window> getFragment(String chr, int start, int end) {
		Collection<Window> rtrn=new TreeSet<Window>();
		if(!this.chromosomeSizes.containsKey(chr)){
			logger.warn(chr+" is not in the genomic space it is possible that the fragment will extend past the chromosome end");
			rtrn.add(new GenomeWindow(chr,start,end));
			return rtrn;
		}
		/*
		 * Check that the fragment coordinates are contained within the coordinate space.
		 * If fragment is entirely outside the coordinates, throw a warning
		 * If one end of the fragment is outside the coordinates, adjust the fragment end.
		 */
		
		if (start > chromosomeSizes.get(chr) || end < 0) {
			logger.warn("Requested " + chr + ":" + start + "-" + end + " is entirely out of bounds. \n" + 
					"Chromosome bounds: " + chr + ":0-" + chromosomeSizes.get(chr));
			throw new AnnotationOutOfBoundsException("Requested " + chr + ":" + start + "-" + end + " is entirely out of bounds. \n" + 
					"Chromosome bounds: " + chr + ":0-" + chromosomeSizes.get(chr));
		}
		
		if(start<0)
			start = 0;
		if(end>this.chromosomeSizes.get(chr))
			end = this.chromosomeSizes.get(chr);
		
		rtrn.add(new GenomeWindow(chr,start,end));
		return rtrn;
	}

	
	@Override
	public Collection<? extends Window> getFragment(Annotation annotation) {
		
		return getFragment(annotation.getChr(), annotation.getStart(), annotation.getEnd());
	}

	@Override
	public int getSize(Annotation region) {
		Collection<? extends Window> fragments = getFragment(region);
		int size = fragments.iterator().next().size();
		for(Window window : fragments) {
			if(window.size() != size) {
				throw new IllegalStateException("Fragment set for annotation consists of multiple regions with different sizes.");
			}
		}
		return size;
	}

	@Override
	public Collection<String> getChromosomeNames() {
		return getReferenceNames();
	}

}

