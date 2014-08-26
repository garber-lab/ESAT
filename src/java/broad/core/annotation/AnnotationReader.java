package broad.core.annotation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.core.error.ParseException;

import broad.pda.datastructures.Alignments;

import nextgen.core.annotation.Annotation;
import nextgen.core.feature.Window;

import org.apache.commons.collections15.Predicate;

import org.apache.commons.math3.stat.*;

public abstract class AnnotationReader <T extends GenomicAnnotation> {
	private List<String> browserLines;
	//private List<String> trackLines;
	private List< AnnotationSet<T>> annotationSetMap; 
	//private List<Map<String, List<T>>> trackSequenceMapList;
	
	public AnnotationReader() {
		super();
		//trackSequenceMapList = new ArrayList<Map<String,List<T>>>();
		init();
	}
	
	public void init() {
		browserLines = new ArrayList<String>();
		//trackLines = new ArrayList<String>();
		annotationSetMap = new ArrayList<AnnotationSet<T>>();
		AnnotationSet<T> first = new AnnotationSet<T>();
		annotationSetMap.add(first);
	}
	
	/**
	 * Jesse Engreitz
	 * September 26, 2012
	 * Allow for creating a reader from a collection of annotations
	 * @param annotations
	 */
	public AnnotationReader(Collection<T> annotations) {
		this();
		for (T annot : annotations) 
			addAnnotation(annot);
	}

	
	public void clear() {
		browserLines.clear();
		annotationSetMap.clear();
	}
	
	public void load(File source, AnnotationFactory<? extends T> factory) throws IOException, ParseException{
		load(new BufferedReader(new FileReader(source)), factory);
	}
	
	public void load(BufferedReader br, AnnotationFactory<? extends T> factory) throws IOException, ParseException{
		load(br, factory, new GenomicAnnotationFilter<T> () {

			public boolean accept(T annotation) {
				return true;
			}

			public boolean isEnough(T annotation) {
				return false;
			}
			
		});
	}
	
	public void load(File source, AnnotationFactory<? extends T> factory, GenomicAnnotationFilter<T> filter) throws IOException, ParseException{
		load(new BufferedReader(new FileReader(source)), factory, filter);
	}
	
	public void load(BufferedReader br, AnnotationFactory<? extends T> factory, GenomicAnnotationFilter<T> filter) throws IOException, ParseException{
		String line;
		try {
			while((line = br.readLine()) != null) {
				line = line.trim();
				if(line.startsWith("#") || line.length() == 0) {
					continue;
				}
				
				if( line.toLowerCase().startsWith("browser")){
					addBrowserLine(line);
					continue;
				}
				
				if( line.toLowerCase().startsWith("track")){ //Assuming only one track per annotation set, TODO: check!
					if(annotationSetMap.size() == 1) {// This is a hack resulting from a bad choice made before to actually initialize the annotationSetMap in the constructor TODO: FIXIT
						annotationSetMap.get(0).setInfoFromRawData(line);
					} else  {
						startTrack(line);
					}
					continue;
				}
				//System.out.println(line.replace("\t", "-TAB-"));
				String [] lineSplit = line.split("\t");
				//System.err.println("Creating " + line);
				T annotation = factory.create(lineSplit);
				
				if(! filter.accept(annotation)) {
					continue;
				}
				
				if(filter.isEnough(annotation)) {
					break;
				}
				
				AnnotationSet<T> currentSet = null;
				if(annotationSetMap.isEmpty()) {
					currentSet = new AnnotationSet<T>(); 
					annotationSetMap.add(currentSet);
				} else {
					currentSet = annotationSetMap.get(annotationSetMap.size() - 1);
				}
				
				currentSet.addAnnotation(annotation);
			}
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			try {
				//System.out.print("Closing "+input);
				br.close();
				//System.out.print(" ..... Closed\n");
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}
	
	public abstract int parse (String file, GenomicAnnotationFilter<T> filter, AnnotationHandler handler ) throws ParseException, IOException;
	protected int parse(String file, AnnotationFactory<? extends T> factory, GenomicAnnotationFilter<T> filter, AnnotationHandler handler) throws ParseException, IOException {
		BufferedReader br = new BufferedReader(new FileReader(file));
		int numRead = 0;
		try {
			numRead = parse(br, factory, filter, handler);
		} finally {
			if(br != null) {
				br.close();
			}
		}
		return numRead;
	}
	
	public abstract int parse (BufferedReader br, GenomicAnnotationFilter<T> filter, AnnotationHandler handler ) throws ParseException, IOException;
	protected int parse(BufferedReader br, AnnotationFactory<? extends T> factory, GenomicAnnotationFilter<T> filter, AnnotationHandler handler) throws ParseException {
		String line;
		handler.begin();
		int numRead = 0;
		try {
			while((line = br.readLine()) != null) {
				line = line.trim();
				if(line.startsWith("#") || line.length() == 0) {
					//comment
				} else 	if( line.toLowerCase().startsWith("browser")){
					
					handler.track(line);
					
				} else if( line.toLowerCase().startsWith("track")){ //Assuming only one track per annotation set, TODO: check!
					handler.browserLine(line);
				} else {
					//System.out.println(line.replace("\t", "-TAB-"));
					String [] lineSplit = line.split("\t");
					T annotation = factory.create(lineSplit);
					
					if(filter == null ||  (filter.accept(annotation) && !filter.isEnough(annotation))) {
						handler.annotation(annotation);
						numRead++;
					}
				}
			}
			handler.eof();
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			try {
				//System.out.print("Closing "+input);
				br.close();
				//System.out.print(" ..... Closed\n");
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}		
		return numRead;
	}
	
	public Map<String, List<T>> getChromosomeAnnotationMap() {
		return getChromosomeAnnotationMap(0);
	}
	
	public Iterator<String> getChromosomeIterator() { return annotationSetMap.get(0).getAnnotationTreeMap().keySet().iterator();}
	
	public ShortBEDReader cloneRegions() {
		ShortBEDReader newReader = new ShortBEDReader();
		for (T annot : getAnnotationList()) {
			ShortBED newBed = new ShortBED(annot);
			newReader.addAnnotation(newBed);
		}
		return newReader;
	}
	
	public Map<String, List<T>> getChromosomeAnnotationMap(int i) {
		return annotationSetMap.size() > i ? annotationSetMap.get(i).getAsChromosomeAnnotationMap() : null;
	}
	
	public int getBaseCoverage(List<T> annotations) {
		int cov = 0;
		Iterator<T> it = annotations.iterator();
		while(it.hasNext()) {
			T annotation = it.next();
			cov = cov +  annotation.length();
		}
		return cov;
	}
	
	public int numberCovered(List<? extends GenomicAnnotation> annotations) {
		int num = 0;
		
		Iterator<? extends GenomicAnnotation> it = annotations.iterator();
		while(it.hasNext()) {
			LightweightGenomicAnnotation annot = it.next();
			IntervalTree<T> tree = getChromosomeTree(annot.getChromosome());
			if(tree != null) {
				Iterator<Node<T>> overlapIt = tree.overlappers(annot.getStart(), annot.getEnd());
				if(overlapIt.hasNext()){
					num++;
				}
			}
		}
		return num;
	}
	
	
	public void intersect(AnnotationReader <T> other) {
		
		HashMap<String, IntervalTree<T>> intersectTree = new HashMap<String, IntervalTree<T>>();
		Iterator<String> chrIt = other.getChromosomeIterator();
		while(chrIt.hasNext()) {
			String chr = chrIt.next();
			Iterator<T> annotationIt = other.getChromosomeTree(chr).valueIterator();
			IntervalTree<T> tree = getChromosomeTree(chr);
			if(tree == null) {
				continue;
			}

			while(annotationIt.hasNext()) {
				T annot = annotationIt.next();
				Iterator<Node<T>> overlaperIt = tree.overlappers(annot.getStart(), annot.getEnd());
				//System.err.println("Overlappers for " + annot.getName() + " got one? " + overlaperIt.hasNext());
				int intersectNum = 0; 
				while(overlaperIt.hasNext()) {
					Node<T> overlapperNode = overlaperIt.next();
					T overlapper = overlapperNode.getValue();
					T intersect = createAnnotation(annot);
					//intersect.setId(overlapper.getId());
					intersect.setName(overlapper.getName()+"_" + annot.getName());
					intersect.setScore(overlapper.getScore());
					
					IntervalTree<T> newTree = intersectTree.get(annot.getChromosome());
					if(newTree == null) {
						newTree = new IntervalTree<T>();
						intersectTree.put(annot.getChromosome(),newTree);
					}
					intersect.takeIntersection(overlapper);
					
					newTree.put(intersect.getStart(), intersect.getEnd(),intersect);
					intersectNum++;
				}
			}
		}
		//return intersectTree;
		setChromosomeTreeMap(0, intersectTree);
		
	}
	
	public void intersect(List<? extends GenomicAnnotation> annotations) {
		
		HashMap<String, IntervalTree<T>> intersectTree = new HashMap<String, IntervalTree<T>>();
		for (GenomicAnnotation annot : annotations) {
			IntervalTree<T> tree = getChromosomeTree(annot.getChromosome());
			if(tree == null) {
				continue;
			}
			Iterator<Node<T>> overlaperIt = tree.overlappers(annot.getStart(), annot.getEnd());
			//System.err.println("Overlappers for " + annot.getName() + " got one? " + overlaperIt.hasNext());
			int intersectNum = 0; 
			while(overlaperIt.hasNext()) {
				Node<T> overlapperNode = overlaperIt.next();
				T overlapper = overlapperNode.getValue();
				T intersect = createAnnotation(annot);
				//intersect.setId(overlapper.getId());
				intersect.setName(overlapper.getName()+"_" + annot.getName());
				intersect.setScore(overlapper.getScore());
				
				IntervalTree<T> newTree = intersectTree.get(annot.getChromosome());
				if(newTree == null) {
					newTree = new IntervalTree<T>();
					intersectTree.put(annot.getChromosome(),newTree);
				}
				intersect.takeIntersection(overlapper);
				
				newTree.put(intersect.getStart(), intersect.getEnd(),intersect);
				intersectNum++;
			}
		}
		//return intersectTree;
		setChromosomeTreeMap(0, intersectTree);
		
	}

	public Collection<T> takeIntersection(List<? extends GenomicAnnotation> annotations) {
		Set<T> intersection = new TreeSet<T>();
		
		for (GenomicAnnotation annot : annotations) {
			IntervalTree<T> tree = getChromosomeTree(annot.getChromosome());
			if(tree == null) {
				continue;
			}
			
			Iterator<Node<T>> overlaperIt = tree.overlappers(annot.getStart(), annot.getEnd());
			//System.err.println("Overlappers for " + annot.getName() + " got one? " + overlaperIt.hasNext());
			int intersectNum = 0; 
			while(overlaperIt.hasNext()) {
				Node<T> overlapperNode = overlaperIt.next();
				T overlapper = overlapperNode.getValue();
				T intersect = createAnnotation(annot);
				//intersect.setId(overlapper.getId());
				intersect.setName(overlapper.getName()+"_" + annot.getName());
				intersect.setScore(overlapper.getScore());
				intersect.takeIntersection(overlapper);
				intersection.add(intersect);
				intersectNum++;
			}
		}

		return intersection;
	}
	
	
	// Find the closest annotation (overlapping ok, 0 distance)
	public T findClosest(final LightweightGenomicAnnotation annotation) {
		return findClosest(annotation, false);
	}
	
	// Find the closest non-overlapping annotation
	public T findClosestNonOverlapping(final LightweightGenomicAnnotation annotation) {
		return findClosest(annotation, true);
	}
	
	/**
	 *  Added "nonOverlapping" option June 13, 2012
	 *  @param nonOverlapping		Set to true to get the closest non-overlapping annotation
	 */
	public T findClosest(final LightweightGenomicAnnotation annotation, boolean nonOverlapping) {
		IntervalTree<T> chrTree = getChromosomeTree(annotation.getChromosome());
		T closest = null;
		if(chrTree != null) {
			Iterator<Node<T>> overlapperIt = chrTree.overlappers(annotation.getStart(), annotation.getEnd());
			if(overlapperIt.hasNext() && !nonOverlapping) {
				closest = overlapperIt.next().getValue();
			} else {
				Node<T> closestAfter = chrTree.min(annotation.getEnd(), annotation.getEnd()+1);
				Node<T> closestBefore = chrTree.max(annotation.getStart()-1, annotation.getStart());
				int distToAfter = closestAfter == null ? Integer.MAX_VALUE : annotation.getDistanceTo(closestAfter.getValue());
				int distToBefore = closestBefore == null ? Integer.MAX_VALUE : annotation.getDistanceTo(closestBefore.getValue());
				if(closestAfter != null && distToAfter < distToBefore) {
					closest = closestAfter.getValue();
				} else if (closestBefore != null) {
					closest = closestBefore.getValue();
				}
			}
		}
		return closest;
	}
	

	/**
	 * Jesse Engreitz
	 * June 24, 2012
	 * Return a map of elements in this reader to the closest elements in another reader
	 * @param otherReader		second set of annotations
	 * @param nonOverlapping	set to true to get the closest non-overlapping elements
	 * @return
	 */
	public Map<T,LightweightGenomicAnnotation> findClosestForAll(final AnnotationReader<? extends GenomicAnnotation> otherReader, boolean nonOverlapping) {
		Map<T,LightweightGenomicAnnotation> map = new TreeMap<T,LightweightGenomicAnnotation>();
		
		Iterator<String> chrIt = this.getChromosomeAnnotationMap().keySet().iterator();
		while (chrIt.hasNext()) {
			String chr = chrIt.next();
			IntervalTree<T> tree = this.getChromosomeTree(chr);
			Iterator<T> annotIt = tree.valueIterator();
			while (annotIt.hasNext()) {
				T annot = annotIt.next();
				LightweightGenomicAnnotation closest = otherReader.findClosest(annot, nonOverlapping);
				// Note: closest might be null if there is no element in otherReader on the same chromosome (e.g. chrM)
				map.put(annot, closest);
			}
		}
		return map;
	}
	public Map<T,LightweightGenomicAnnotation> findClosestForAll(final AnnotationReader<? extends GenomicAnnotation> otherReader) { return findClosestForAll(otherReader, false); }
	public Map<T,LightweightGenomicAnnotation> findClosestNonoverlappingForAll(final AnnotationReader<? extends GenomicAnnotation> otherReader) { return findClosestForAll(otherReader, true); }	
	
	
	/**
	 * Jesse - June 13, 2012
	 * @param annotation
	 * @return
	 */
	public int getDistanceToClosest(final LightweightGenomicAnnotation annotation) {
		T closest = findClosest(annotation);
		int rtrn = Integer.MAX_VALUE;
		if (closest != null) {
			rtrn = annotation.getDistanceTo(closest);
		}
		return rtrn;
	}
	
	
	public Map<T, Integer> getDistancesToClosest(final AnnotationReader<? extends GenomicAnnotation> otherReader) {
		Map<T,Integer> map = new TreeMap<T,Integer>();
		
		Iterator<String> chrIt = this.getChromosomeAnnotationMap().keySet().iterator();
		while (chrIt.hasNext()) {
			String chr = chrIt.next();
			IntervalTree<T> tree = this.getChromosomeTree(chr);
			Iterator<T> annotIt = tree.valueIterator();
			while (annotIt.hasNext()) {
				T annot = annotIt.next();
				int distance = otherReader.getDistanceToClosest(annot);
				map.put(annot, distance);
			}
		}
		return map;
	}
	
	
	/**
	 * Jesse Engreitz
	 * June 24, 2012
	 * Takes the annotations in the reader and permutes them on their respective chromosomes.
	 * Note that this alters the positions of annotations in the reader, although not their sizes;
	 * @param sizes		map of chromosome sizes
	 */
	public void permuteAnnotationsOnChromosome(final Map<String, Integer> sizes) {
		Collection<T> newReads = new ArrayList<T>();
		Iterator<String> chrIt = this.getChromosomeAnnotationMap().keySet().iterator();
		while (chrIt.hasNext()) {
			String chr = chrIt.next();
			if (!sizes.containsKey(chr)) continue; 
			int size = sizes.get(chr);
			IntervalTree<T> tree = this.getChromosomeTree(chr);
			Iterator<T> itr = tree.valueIterator();
			while (itr.hasNext()) {
				T curr = itr.next();
				int start = 0 + new Double(Math.random()*(size-curr.length())).intValue();
				int end = start + curr.length();
				curr.setStart(start);
				curr.setEnd(end);
				newReads.add(curr);
			}
		}
		
		init();
		for (T t : newReads) addAnnotation(t);
	}
	
	
	/**
	 * Takes the complement in the given list in this Reader. That is, if an annotation in the list
	 * reader does not overlap any in the this reader, the annotation is kept intact.
	 * If it is contained in any annotation in the reader then it is removed. If it incompletely overlapped, then 
	 * the annotation resulting from taking out overlapping regions replaces the overlapped one.
	 */
	public BasicAnnotationReader takeComplement(List<? extends Annotation> annotations) {
		
		BasicAnnotationReader complement = new BasicAnnotationReader();
		for (Annotation a: annotations) {
			List<Annotation> overlapers = getOverlappers(a);
			//List<Annotation> annotationComplement = a.minus(overlapers);
			//for(Annotation ac : annotationComplement) {
				complement.addAnnotation(complement.createAnnotation(a.minus(overlapers)));
			//}
		}
		complement.merge();
		return complement;
	}
	
	public List<Annotation> getOverlappers(Annotation a) {
		List<Annotation> list=new ArrayList<Annotation>();
		list.add(a);
		return getOverlappers(list);
	}

	/**
	 * Removes from this reader all elements that overlap the given list.
	 * @param annotations
	 */
	public void minus(List<? extends GenomicAnnotation> annotations) {
		Iterator<? extends GenomicAnnotation> it = annotations.iterator();
		List<T> toRemove = new ArrayList<T>();
		while(it.hasNext()) {
			LightweightGenomicAnnotation annot = it.next();
			IntervalTree<T> tree = getChromosomeTree(annot.getChromosome());
			if(tree == null) {
				continue;
			}
			Iterator<T> overlaperIt = new IntervalTree.ValuesIterator<T>(tree.overlappers(annot.getStart(), annot.getEnd()));
			while(overlaperIt.hasNext()) {
				T overlapper = overlaperIt.next();
				//System.out.println("Annotation " + annot.toString() + " has overlapper: " + overlapper.toString());
				if(!toRemove.contains(overlapper) ) {
					toRemove.add(overlapper);
				}
			}
		}
		
		Iterator<T> toRemoveIt = toRemove.iterator();
		while(toRemoveIt.hasNext()) {
			T t = toRemoveIt.next();
			IntervalTree<T> tree = getChromosomeTree(t.getChromosome());
			tree.remove(t.getStart(), t.getEnd());
		}
	}
	
	/**
	 * Filters this reader by removing all elements that do not overlap the given list of annotations
	 * @param annotations a list of Annotations to filter
	 * 
	 */
	public void filterByOverlap(List<? extends Annotation> annotations) {
		filterByOverlap(annotations, 0);
		/*
		HashMap<String, IntervalTree<T>> intersectTree = new HashMap<String, IntervalTree<T>>();
		Iterator<LightweightGenomicAnnotation> it = annotations.iterator();
		while(it.hasNext()) {
			LightweightGenomicAnnotation annot = it.next();
			IntervalTree<T> tree = getChromosomeTree(annot.getChromosome());
			if(tree == null) {
				continue;
			}
			
			Iterator<Node<T>> overlaperIt = tree.overlappers(annot.getStart(), annot.getEnd());
			while(overlaperIt.hasNext()) {
				Node<T> overlapperNode = overlaperIt.next();
				T overlapper = overlapperNode.getValue();
				
				IntervalTree<T> newTree = intersectTree.get(annot.getChromosome());
				if(newTree == null) {
					newTree = new IntervalTree<T>();
					intersectTree.put(annot.getChromosome(),newTree);
				}
				if (newTree.find(overlapperNode.getStart(), overlapperNode.getEnd()) == null) {
					newTree.put(overlapper.getStart(), overlapper.getEnd(), overlapper);
				}
			}
		}
		//return intersectTree;
		setChromosomeTreeMap(0, intersectTree);
		*/
	}
	
	public void filterByOverlap(List<? extends Annotation> annotations, int extensionFactor) {
		HashMap<String, IntervalTree<T>> intersectTree = new HashMap<String, IntervalTree<T>>();
		Iterator<? extends Annotation> it = annotations.iterator();

		
		while(it.hasNext()) {
			Annotation annot = it.next();
			IntervalTree<T> tree = getChromosomeTree(annot.getChr());
			if(tree == null) {
				continue;
			}
			
			Iterator<Node<T>> overlaperIt = tree.overlappers(annot.getStart() - extensionFactor, annot.getEnd() + extensionFactor);
			while(overlaperIt.hasNext()) {
				Node<T> overlapperNode = overlaperIt.next();
				T overlapper = overlapperNode.getValue();
				
				IntervalTree<T> newTree = intersectTree.get(annot.getChr());
				if(newTree == null) {
					newTree = new IntervalTree<T>();
					intersectTree.put(annot.getChr(),newTree);
				}
				if (newTree.find(overlapperNode.getStart(), overlapperNode.getEnd()) == null) {
					newTree.put(overlapper.getStart(), overlapper.getEnd(), overlapper);
				}
			}
		}
		//return intersectTree;
		setChromosomeTreeMap(0, intersectTree);
	}
	
	public void filterByOverlapOfLastExonOrAfterLastExon(List<? extends GenomicAnnotation> annotations, int extensionFactor) {
		HashMap<String, IntervalTree<T>> intersectTree = new HashMap<String, IntervalTree<T>>();
		Iterator<? extends GenomicAnnotation> it = annotations.iterator();
		while(it.hasNext()) {
			GenomicAnnotation annot = it.next();
			IntervalTree<T> tree = getChromosomeTree(annot.getChromosome());
			if(tree == null) {
				continue;
			}
			
			if (annot.getOrientation().equals("*"))
				continue;
			
			Annotation lastExon = null;
			int annotEnd = annot.getOrientedEnd();
			List<? extends Annotation> annotBlocks = annot.getBlocks();
			Iterator<? extends Annotation> annotBlocksIt = annotBlocks.iterator();
			while (annotBlocksIt.hasNext()) {
				Annotation block = annotBlocksIt.next();
				if (block.getStart() == annotEnd || block.getEnd() == annotEnd) {
					lastExon = block;
					break;
				}
			}
			
			int start;
			int end;
			if (lastExon.getStart() == annotEnd) {
				start = annotEnd - extensionFactor;
				end = lastExon.getEnd();
			} else {
				start = lastExon.getStart();
				end = annotEnd + extensionFactor;
			}
			
			Iterator<Node<T>> overlaperIt = tree.overlappers(start, end);
			while(overlaperIt.hasNext()) {
				Node<T> overlapperNode = overlaperIt.next();
				T overlapper = overlapperNode.getValue();
				
				IntervalTree<T> newTree = intersectTree.get(annot.getChromosome());
				if(newTree == null) {
					newTree = new IntervalTree<T>();
					intersectTree.put(annot.getChromosome(),newTree);
				}
				if (newTree.find(overlapperNode.getStart(), overlapperNode.getEnd()) == null) {
					newTree.put(overlapper.getStart(), overlapper.getEnd(), overlapper);
				}
			}
		}
		//return intersectTree;
		setChromosomeTreeMap(0, intersectTree);
	}
	
	public void filterByOverlapOfLastExon(List<? extends GenomicAnnotation> annotations, int extensionFactor) {
		HashMap<String, IntervalTree<T>> intersectTree = new HashMap<String, IntervalTree<T>>();
		Iterator<? extends GenomicAnnotation> it = annotations.iterator();
		while(it.hasNext()) {
			GenomicAnnotation annot = it.next();
			IntervalTree<T> tree = getChromosomeTree(annot.getChromosome());
			if(tree == null) {
				continue;
			}
			
			Annotation lastExon = null;
			int annotEnd = annot.getOrientedEnd();
			List<? extends Annotation> annotBlocks = annot.getBlocks();
			Iterator<? extends Annotation> annotBlocksIt = annotBlocks.iterator();
			while (annotBlocksIt.hasNext()) {
				Annotation block = annotBlocksIt.next();
				if (block.getStart() == annotEnd || block.getEnd() == annotEnd) {
					lastExon = block;
					break;
				}
			}
			
			Iterator<Node<T>> overlaperIt = tree.overlappers(lastExon.getStart() - extensionFactor, lastExon.getEnd() + extensionFactor);
			while(overlaperIt.hasNext()) {
				Node<T> overlapperNode = overlaperIt.next();
				T overlapper = overlapperNode.getValue();
				
				IntervalTree<T> newTree = intersectTree.get(annot.getChromosome());
				if(newTree == null) {
					newTree = new IntervalTree<T>();
					intersectTree.put(annot.getChromosome(),newTree);
				}
				if (newTree.find(overlapperNode.getStart(), overlapperNode.getEnd()) == null) {
					newTree.put(overlapper.getStart(), overlapper.getEnd(), overlapper);
				}
			}
		}
		//return intersectTree;
		setChromosomeTreeMap(0, intersectTree);
	}
	
	/**
	 * ADDED BY MORAN Feb 8th 2010
	 * Filters this reader by removing all elements that do not overlap AND HAVE THE SAME ORIENTATION as the given list of annotations 
	 * @param annotations a list of Annotations to filter
	 * 
	 */
	public void filterByOverlapAndOrientation(List<? extends GenomicAnnotation> annotations) {
		HashMap<String, IntervalTree<T>> intersectTree = new HashMap<String, IntervalTree<T>>();
		Iterator<? extends GenomicAnnotation> it = annotations.iterator();
		while(it.hasNext()) {
			LightweightGenomicAnnotation annot = it.next();
			
			
			IntervalTree<T> tree = getChromosomeTree(annot.getChromosome());
			if(tree == null) {
				continue;
			}
			//Get every annotation in this tree (the current object, this)
			//that is between the start-end of the annotation in set2 
			Iterator<Node<T>> overlaperIt = tree.overlappers(annot.getStart(), annot.getEnd());
			while(overlaperIt.hasNext()) {
				Node<T> overlapperNode = overlaperIt.next();
				T overlapper = overlapperNode.getValue();
				
				IntervalTree<T> newTree = intersectTree.get(annot.getChromosome());
				if(newTree == null) {
					newTree = new IntervalTree<T>();
					intersectTree.put(annot.getChromosome(),newTree);
				}
				if (newTree.find(overlapperNode.getStart(), overlapperNode.getEnd()) == null && overlapper.getOrientation().equals(annot.getOrientation())) {
					newTree.put(overlapper.getStart(), overlapper.getEnd(), overlapper);
				}
			}
		}
		//return intersectTree;
		setChromosomeTreeMap(0, intersectTree);
	}
	
	public static Map<Annotation, Integer> getLastExons(List<? extends Annotation> annotations) {
		Map<Annotation, Integer> lastExons = new HashMap<Annotation, Integer>();
		
		Iterator<? extends Annotation> it = annotations.iterator();
		while(it.hasNext()) {
			Annotation annot = it.next();
			
			Annotation lastExon = null;
			int annotEnd = annot.getOrientedEnd();
			
			List<? extends Annotation> annotBlocks = annot.getBlocks();
			Iterator<? extends Annotation> annotBlocksIt = annotBlocks.iterator();
			while (annotBlocksIt.hasNext()) {
				Annotation block = annotBlocksIt.next();
				if (block.getStart() == annotEnd || block.getEnd() == annotEnd) {
					lastExon = block;
					break;
				}
			}
			lastExons.put(lastExon, annot.getOrientedEnd());
		}
		
		return lastExons;
	}
	
	/**
	 * Moran FEB23 2010
	 * Increments by one the score of the annotations in the current object if they overlap the input
	 * set. Set to zero flag will first initialize all scores to zero
	 * @return
	 */
	public void IncrementScoreIfOverlap(AnnotationReader<? extends GenomicAnnotation> annotReader,int setNumber, Boolean setToZero) {
		
		Map<String, IntervalTree<T>> originalTreeMap = getChromosomeTreeMap(setNumber);
		Iterator<String> chrIt = originalTreeMap.keySet().iterator();
		while(chrIt.hasNext()) {
			String chr = chrIt.next();
			IntervalTree<T> oldTree = originalTreeMap.get(chr);
			Iterator<T> valueIt = new IntervalTree.ValuesIterator<T>(oldTree.iterator());
			while(valueIt.hasNext()) {
				T currentElement = valueIt.next();
				if (setToZero)
				{
					currentElement.setScore(0);
				}
				if (! annotReader.getOverlappers(currentElement).isEmpty())
				{
					currentElement.setScore(currentElement.getScore()+1);
				}
			}
		}
	}
	
	
	
	public int getBaseCoverage() { return getBaseCoverage(getAnnotationList());}
	
	public abstract T createAnnotation(GenomicAnnotation a);
	
	public void shift (int amountToShift) {
		Iterator<T> it = getAnnotationList().iterator();
		while(it.hasNext()) {
			T annotation = it.next();
			IntervalTree<T> tree = getChromosomeTree(annotation.getChromosome());
			tree.remove(annotation.getStart(), annotation.getEnd());
			annotation.setStart(annotation.getStart() + amountToShift);
			annotation.setEnd(annotation.getEnd() + amountToShift);
			tree.put(annotation.getStart(), annotation.getStart(),annotation);
			
		}
		
	}

	/**
	 * Merge list of annotations by combining overlapping annotations into one.
	 */
	public void merge () { merge(0);}
	
	/**
	 * Gets a list of overlapping annotations
	 */
	public List<Annotation> getOverlappers(List<? extends Annotation> annotationList) {
		return getOverlappers(annotationList, false);
	}
		
		
	public List<Annotation> getOverlappers(List<? extends Annotation> annotationList, boolean fullyContained) {
		List<Annotation> overlappers = new ArrayList<Annotation>();
		Iterator<? extends Annotation> annotationIt = annotationList.iterator();
		while(annotationIt.hasNext()) {
			Annotation annot = annotationIt.next();
			IntervalTree<T> tree = getChromosomeTree(annot.getChr());
			if(tree == null) {
				continue;
			}
			Iterator<T> overlapperIt = new IntervalTree.ValuesIterator<T>(tree.overlappers(annot.getStart(), annot.getEnd()));
			while(overlapperIt.hasNext()) {
				T overlapper = overlapperIt.next();
				if(!overlappers.contains(overlapper)) {
					if (!fullyContained || (fullyContained && annot.contains(overlapper))) {
						overlappers.add(overlapper);
					}
				}
			}
		}
		return overlappers;
	} 
	
	public List<Annotation> getOverlappers(AnnotationReader<? extends GenomicAnnotation> other) {
		return getOverlappers(other.getAnnotationList());
	}

	
	/**
	 * Adds a 0 fudge fact
	 * @param region
	 * @return
	 * @throws IOException
	 */
	public int getBasesCovered(Annotation region) {
		return getBasesCovered(region, 0);
	}
	
	public int getBasesCovered(Annotation region, boolean fullyContained) {
		return getBasesCovered(region, 0, fullyContained);
	}
	
	public int getBasesCovered(Annotation region, int fudgeFactor) {
		return getBasesCovered(region, fudgeFactor, false);
	}
	
	public int getBasesCovered(Annotation region, int fudgeFactor, boolean fullyContained) {
		// initialize:
		int sum=0;
		if(region.getBlocks().size()>1){
			for(Annotation annotation: region.getBlocks()){
				//This is by defnition a single interval
				sum+=getBasesCovered(annotation);
			}
		}
		else{
		
			boolean covered[] = new boolean[region.length()];
			for (int i=0;i<region.length();i++) covered[i] = false;
	
			Iterator<Annotation> iter = getOverlappers(region.getBlocks(), fullyContained).iterator();
			while (iter.hasNext()) {
				Annotation record = iter.next();
				record.setStart(record.getStart() - fudgeFactor);
				record.setEnd(record.getEnd() + fudgeFactor);
				
				int startIndex = Math.max(0, record.getStart() - region.getStart());
				int endIndex = Math.min(region.length(), record.getEnd() - region.getStart());
	
				for (int i = startIndex; i < endIndex; i++) {
					covered[i] = true;
				}
			}
			
			int counter = 0;
			for (int i = 0; i < covered.length; i++) counter += covered[i] ? 1 : 0;
			sum+= counter;
		}
		
		return sum;
	}
	
	
	
	/**
	 * Gets a list of overlapping annotations
	 *
	public Collection<T> getOverlappers(AnnotationReader<T> other) {
		Set<T> overlappers = new TreeSet<T>();
		Iterator<String> chrIt = other.getChromosomeIterator();
		while(chrIt.hasNext()) {
			String chr = chrIt.next();
			Iterator<T> annotationIt = other.getChromosomeTree(chr).valueIterator();
			IntervalTree<T> tree = getChromosomeTree(chr);
			if(tree == null) {
				continue;
			}
			while(annotationIt.hasNext()) {
				LightweightGenomicAnnotation annot = annotationIt.next();
				Iterator<T> overlapperIt = new IntervalTree.ValuesIterator<T>(tree.overlappers(annot.getStart(), annot.getEnd()));
				while(overlapperIt.hasNext()) {
					T overlapper = overlapperIt.next();
					overlappers.add(overlapper);
				}
			}
		}
		return overlappers;
	} */
	
	public List<T> getOverlappers(LightweightGenomicAnnotation annotation) {
		Set<T> overlappers = new TreeSet<T>();

		IntervalTree<T> tree = getChromosomeTree(annotation.getChromosome());
		if(tree != null) {
			Iterator<T> overlapperIt = new IntervalTree.ValuesIterator<T>(tree.overlappers(annotation.getStart(), annotation.getEnd()));
			while(overlapperIt.hasNext()) {
				T overlapper = overlapperIt.next();
				if(!overlappers.contains(overlapper)) {
					overlappers.add(overlapper);
				}
			}	
		}
		
		return new ArrayList<T>(overlappers);
	}
	
	
	public double getScoreSum(LightweightGenomicAnnotation region) {
		List<T> overlappers = getOverlappers(region);
		double sum = 0;
		for (int i = 0; i < overlappers.size(); i++) {
			sum += overlappers.get(i).getScore();
		}
		return sum;
	}
	
	public double getAverageScore(LightweightGenomicAnnotation region) {
		List<T> overlappers = getOverlappers(region);
		double[] values = new double[overlappers.size()];
		for (int i = 0; i < values.length; i++) {
			values[i] = overlappers.get(i).getScore();
		}
		double result = StatUtils.mean(values);
		return result;
	}
	
	/**
	 * Merge list of annotations by combining overlapping annotations into one.
	 */
	public void merge (int setNumber) {
		Map<String, IntervalTree<T>> mergedTreeMap = new HashMap<String, IntervalTree<T>>();
		Map<String, IntervalTree<T>> originalTreeMap = getChromosomeTreeMap(setNumber);
		
		Iterator<String> chrIt = originalTreeMap.keySet().iterator();
		while(chrIt.hasNext()) {
			String chr = chrIt.next();
			IntervalTree<T> newTree = new IntervalTree<T>();
			mergedTreeMap.put(chr, newTree);
			IntervalTree<T> oldTree = originalTreeMap.get(chr);
			Iterator<T> valueIt = new IntervalTree.ValuesIterator<T>(oldTree.iterator());
			while(valueIt.hasNext()) {
				T currentElement = valueIt.next();
				Iterator<T> overlapperIt = new IntervalTree.ValuesIterator<T>(newTree.overlappers(currentElement.getStart(), currentElement.getEnd()));
				if(overlapperIt.hasNext()) {
					T overlapper = overlapperIt.next();
					if(overlapperIt.hasNext()) {
						throw new IllegalStateException("Two elements overlapped the current element in the new merged tree, this should never happen");
					}
					newTree.remove(overlapper.getStart(),overlapper.getEnd());
					overlapper.takeUnion(currentElement);
					newTree.put(overlapper.getStart(), overlapper.getEnd(), overlapper);
				} else {
					newTree.put(currentElement.getStart(), currentElement.getEnd(), currentElement);
				}
				
				
			}
		}
		
		setChromosomeTreeMap(setNumber, mergedTreeMap);
	}
	
	
	
	
	/**
	 * Moran, Feb22 2010
	 * concat  an annotation set to an exciting one (without merging)
	 * WARNING: assumes that same types are concatenated- uses casting to current 
	 * objets annotation type
	 */
	
	@SuppressWarnings("unchecked")
	public void concatAnnotationSet (int setNumber, List<? extends GenomicAnnotation> list){
		
		Iterator<? extends GenomicAnnotation> it = list.iterator();
		while(it.hasNext()) {
			T annot = (T) it.next();
			addAnnotation(annot,setNumber);
		}
	}
	
	public void addAnnotation(T annotation, int track) {
		AnnotationSet<T> set = annotationSetMap.get(track);
		if(set != null) {
			set.addAnnotation(annotation);
		} else {
			throw new IndexOutOfBoundsException("No track for given index " + track);
		}
	}
	
	public void addAnnotation(T annotation) {addAnnotation(annotation, 0);}
	
	public void addAnnotationToLastTrack(T annotation) {addAnnotation(annotation, annotationSetMap.size() - 1);}
	
	public List<T> getAnnotationsForSequence(String sequenceName) {
		return getChromosomeAnnotationMap(0).get(sequenceName);
	}
	
	public List<T> getAnnotationsForSequence(String sequenceName, int setNumber) {
		return getChromosomeAnnotationMap(setNumber).get(sequenceName);
	}
	
	public Map<String, List<T>> getSequenceAnnotationMap() {
		return getChromosomeAnnotationMap(0);
	}
	
	public Map<String, List<T>> getSequenceAnnotationMap(int setNumber) {
		return getChromosomeAnnotationMap(setNumber);
	}
	
	public void addBrowserLine(String line) {
		browserLines.add(line);
	}
	
	public void startTrack(String line) {
		AnnotationSet<T> track = new AnnotationSet<T>(); 
		track.setInfoFromRawData(line);
		annotationSetMap.add( track);		
	}
	
	public void startTrack(String line, int trackNum) {
		AnnotationSet<T> track = new AnnotationSet<T>(); 
		track.setInfoFromRawData(line);
		annotationSetMap.set(trackNum, track);
	}
	
	public IntervalTree<T> getChromosomeTree(String chr, int i) {
		return annotationSetMap.get(i).getAnnotationTree(chr);
	}
	
	public IntervalTree<T> getChromosomeTree(String chr) {
		return annotationSetMap.get(0).getAnnotationTree(chr);
	}
	
	private void setChromosomeTreeMap(int i, Map<String, IntervalTree<T>> treeMap) {
		annotationSetMap.get(i).setAnnotationTreeMap(treeMap);
	}
	
	public Map<String,IntervalTree<T>> getChromosomeTreeMap(int i) {
		return annotationSetMap.get(i).getAnnotationTreeMap();
	}
	
	public Map<String,IntervalTree<T>> getChromosomeTreeMap() {
		return annotationSetMap.get(0).getAnnotationTreeMap();
	}
	
	
	public String getTrackInfo(int i) {
		return annotationSetMap.get(i).getInfoAsString();
	}
	
	public List<T> getChromosomeBEDs(String chr) {
		return getChromosomeAnnotationMap().get(chr);
	}
	
	protected void addAnnotationSet(AnnotationSet<T> track) {
		annotationSetMap.add(track);
	}
	
	public List<T> getChromosomeBEDs(String chr, int i) {
		return getChromosomeAnnotationMap(i).get(chr);
	}	
	
	public List<T> getAnnotationList(int setNumber) {
		Map<String, List<T>> map = getSequenceAnnotationMap(setNumber);
		List<T> all = new ArrayList<T>();
		Iterator<List<T>> sequenceAnnotationIt = map.values().iterator();
		while(sequenceAnnotationIt.hasNext()) {
			all.addAll(sequenceAnnotationIt.next());
		}
		
		return all;
	}
	
	public T getAnnotation(int setNumber, String name ) {
		AnnotationSet<T> set = annotationSetMap.get(setNumber);
		return set.getAnnotation(name);
	}
	
	public boolean containsAnnotation(int setNumber, String name ) {
		AnnotationSet<T> set = annotationSetMap.get(setNumber);
		return set.containsAnnotation(name);
	}
	
	public T getAnnotation(String name ) {
		return getAnnotation(0, name);
	}
	
	public boolean containsAnnotation(String name ) {
		return containsAnnotation(0, name);
	}
	
	public void shuffleScores(String chr, int setNum) {
		Random r = new Random();
		List<T> original = getChromosomeBEDs(chr, setNum);
		int size = original.size();
		for(int i = 0; i < size; i++) {
			int rIdx = i + r.nextInt(size - i);
			double newScore = original.get(rIdx).getScore();
			original.get(rIdx).setScore(original.get(i).getScore());
			original.get(i).setScore(newScore);
		}
	}
	
	public void shuffleScores(String chr) {
		shuffleScores(chr,0);
	}
	
	public List<T> getAnnotationList() { return getAnnotationList(0);};
	
	public static class AnnotationSet <T extends LightweightGenomicAnnotation>{
		static final Pattern trackFreeFormAttributes = Pattern.compile("([^ ]+)= *\"(.+)\"");
		static final Pattern trackDeterminedAttributes =  Pattern.compile("([^ ]+)= *(.+)");
		String name;
		HashMap<String,String> setFreeFormInfo;
		HashMap<String, String> setDeterminedInfo;
		
		Map<String, IntervalTree<T>> chromosomeAnnotationTree;
		Map<String, T> nameAnnotationMap;
		
		public AnnotationSet() {
			super();
			setFreeFormInfo = new HashMap<String, String>();
			setDeterminedInfo = new HashMap<String, String>();
			chromosomeAnnotationTree = new HashMap<String, IntervalTree<T>>();
			nameAnnotationMap = new HashMap<String, T>();
		}

		public String getInfoAsString() {
			StringBuilder buf = new StringBuilder("track ");
			if(name != null) {
				buf.append("name=").append(name).append("  ");
			}
			Iterator<String> attributeIt = setFreeFormInfo.keySet().iterator();
			while(attributeIt.hasNext()) {
				String attribute = attributeIt.next();
				buf.append(attribute).append("=\"").append(setFreeFormInfo.get(attribute)).append("\"  ");
			}
			
			attributeIt = setDeterminedInfo.keySet().iterator();
			while(attributeIt.hasNext()) {
				String attribute = attributeIt.next();
				buf.append(attribute).append("=").append(setDeterminedInfo.get(attribute)).append("  ");
			}
			return buf.toString();
		}

		public void setInfoFromRawData(String rawInfo) {
			Matcher m = trackFreeFormAttributes.matcher(rawInfo);
			while(m.find()) {
				String attribute = m.group(0);
				String value     = m.group(1);
				if("name".equals(attribute)) {
					this.name = value;
				} else {
					setFreeFormInfo.put(attribute, value);
				}
					
			}
			
			m = trackDeterminedAttributes.matcher(rawInfo);
			while(m.find()) {
				String attribute = m.group(0);
				String value     = m.group(1);
				if("name".equals(attribute)) {
					this.name = value;
				} else {
					setDeterminedInfo.put(attribute, value);
				}
					
			}
		}
		
		public void setName(String name) { this.name = name;}
		public String getName() {return name;} 
		
		public void addAnnotation(T annotation) {
			String chr = annotation.getChromosome();
			IntervalTree<T> chrAnnotationTree = chromosomeAnnotationTree.get(chr);
			if(chrAnnotationTree == null) {
				chrAnnotationTree = new IntervalTree<T>();
				chromosomeAnnotationTree.put(chr, chrAnnotationTree);
			}
			if(annotation.getStart() > annotation.getEnd()) {
				System.err.println("Annotation " + annotation + " start is after end. Skipping");
			} else {
				chrAnnotationTree.put(annotation.getStart(), annotation.getEnd(), annotation);
				nameAnnotationMap.put(annotation.getName(), annotation);
			}
		}
		
		public IntervalTree<T> getAnnotationTree(String chr) {
			return chromosomeAnnotationTree.get(chr);
		}
		
		boolean containsAnnotation(String name ) { return nameAnnotationMap.containsKey(name);}
		
		public T getAnnotation(String name) {
			return nameAnnotationMap.get(name);
		}
		
		public Map<String, IntervalTree<T>> getAnnotationTreeMap() {
			return chromosomeAnnotationTree;
		}
		
		public void setAnnotationTreeMap(Map<String, IntervalTree<T>> map) {
			this.chromosomeAnnotationTree = map;
		}
		
		public Map<String, List<T>> getAsChromosomeAnnotationMap() {
			HashMap<String, List<T>> chrAnnotationMap = new HashMap<String, List<T>>(chromosomeAnnotationTree.size());
			Iterator<String> chrIt = chromosomeAnnotationTree.keySet().iterator();
			while(chrIt.hasNext()) {
				String chr = chrIt.next();
				IntervalTree<T> chrTree = chromosomeAnnotationTree.get(chr);
				List<T> chrAnnotations = new ArrayList<T>(chrTree.size());
				chrAnnotationMap.put(chr, chrAnnotations);
				Iterator<T> annotationIt = new IntervalTree.ValuesIterator<T>(chrTree.iterator());
				while(annotationIt.hasNext()) {
					chrAnnotations.add(annotationIt.next());
				}
			}

			return chrAnnotationMap;
		}
		
		public int size(){
			int res=0;
			Iterator<String> chrIt = chromosomeAnnotationTree.keySet().iterator();
			while(chrIt.hasNext()) {
				String chr = chrIt.next();
				IntervalTree<T> chrTree = chromosomeAnnotationTree.get(chr);
				res+=chrTree.size();
			}
			return res;
		}
		
	}
	
	public static class trivialFilter implements GenomicAnnotationFilter< GenomicAnnotation> {

	
		public boolean accept(GenomicAnnotation annotation) {

			return true;
		}


		public boolean isEnough(GenomicAnnotation annotation) {

			return false;
		}
		
	}

	public int numOfAnnotation(int setNum) {
		return(this.annotationSetMap.get(setNum).size());
	}
	
	/**
	 * Jesse June 23, 2012
	 * @return 		total number of annotations 
	 */
	//@Override
	public int size() {
		int size = 0;
		Iterator<AnnotationSet<T>> itr = annotationSetMap.iterator();
		while (itr.hasNext()) size = size + itr.next().size();
		return size;
	}

	
	/**
	 * Jesse June 23, 2012
	 * @return 		string representation of all annotations, one per line
	 */
	public String toStringForWriting() {
		// TODO: Are there multiple different ways to convert this object to a String?
		StringBuffer buffer = new StringBuffer();
		Iterator<String> chrIt = getChromosomeAnnotationMap().keySet().iterator();
		while (chrIt.hasNext()) {
			String chr = chrIt.next();
			IntervalTree<T> tree = getChromosomeTree(chr);
			Iterator<T> annotIt = tree.valueIterator();
			while (annotIt.hasNext()) {
				T annot = annotIt.next();
				buffer.append(annot.toString());
				buffer.append("\n");
			}
		}
		return buffer.toString();
	}
	
	/**
	 * @author mgarber
	 * @param BufferedWriter - buffer to write annotations to
	 * @return 		writes all annotations to the provided BufferedWriter
	 * @throws IOException 
	 */
	public void write(BufferedWriter bw) throws IOException {
		Iterator<String> chrIt = getChromosomeAnnotationMap().keySet().iterator();
		while (chrIt.hasNext()) {
			String chr = chrIt.next();
			IntervalTree<T> tree = getChromosomeTree(chr);
			Iterator<T> annotIt = tree.valueIterator();
			while (annotIt.hasNext()) {
				T annot = annotIt.next();
				bw.write(annot.toString());
				bw.newLine();
			}
		}
	}

	/**
	 * @author mgarber
	 * @param BufferedWriter - buffer to write annotations to
	 * @return 		writes all annotations to the provided BufferedWriter
	 * @throws IOException 
	 */
	public void write(String outputFile) throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(outputFile));
		write(bw);
		bw.close();
	}

	public void filterByScore(double set1FilterScore) {
		// TODO Auto-generated method stub		
	}
	

	public void addFilter(Predicate<T> filter) {
		throw new UnsupportedOperationException("not implemented yet");
	}

	
	public void addFilters(Collection<Predicate<T>> filters) {
		throw new UnsupportedOperationException("not implemented yet");
	}
}
