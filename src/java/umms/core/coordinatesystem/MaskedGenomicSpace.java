package umms.core.coordinatesystem;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.log4j.Logger;

import broad.core.annotation.ShortBEDReader;
import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.core.error.ParseException;
import broad.pda.datastructures.Alignments;

import umms.core.alignment.Alignment;
import umms.core.annotation.Annotation;
import umms.core.annotation.BasicAnnotation;
import umms.core.annotation.Gene;
import umms.core.feature.GeneWindow;
import umms.core.feature.GenomeWindow;
import umms.core.feature.Window;

public class MaskedGenomicSpace implements CoordinateSpace{

	static Logger logger = Logger.getLogger(MaskedGenomicSpace.class.getName());
	private boolean constantWidth; //This is whether to consider all windows of constant size
	private TranscriptomeSpace gappedSpace;
	private GenomicSpace continuousSpace;
	private Map<String, Collection<Gene>> genes;
	//private Map<String, Annotation> excluded;
	//private Map<String, IntervalTree<Annotation>> excludedRegions;
	
	/**
	 * This is meant to represent a coordinate space that excludes masked regions
	 * This way we wont count reads in the global stats in these regions
	 * Wont use them in estimating parameters or length calculations
	 * @param maskedRegions a collection of masked region
	 * @param chrSizes the lengths of the associated chromosomes
	 */
	//TODO Make constructor take  extra option to account for fixed or relative position
	public MaskedGenomicSpace(Map<String, Integer> chrSizes, Map<String, Collection<Annotation>> maskedRegions, int distance) {
		this(chrSizes, maskedRegions, distance, true);
	}

	/**
	 * This is meant to represent a coordinate space that excludes masked regions
	 * This way we wont count reads in the global stats in these regions
	 * Wont use them in estimating parameters or length calculations
	 * @param maskedRegions a collection of masked region
	 * @param chrSizes the lengths of the associated chromosomes
	 * @param constantWidth a boolean (true= ensure that windows have the same final length even if represent different genomic space, false= all windows will have same genomic space even if shorter window length because of masked regions)
	 */
	public MaskedGenomicSpace(Map<String, Integer> chrSizes, Map<String, Collection<Annotation>> maskedRegions, int distance, boolean constantWidth){
		//this.excludedRegions=makeTree(maskedRegions);
		//this.excluded=new TreeMap<String, Annotation>();
		this.genes=makeGenes(chrSizes, maskedRegions, distance);
		this.gappedSpace=new TranscriptomeSpace(genes);
		this.continuousSpace=new GenomicSpace(chrSizes);
		this.constantWidth=constantWidth;
	}
	
	
	private Map<String, Collection<Gene>> makeGenes(Map<String, Integer> chrSizes, Map<String, Collection<Annotation>> excludedRegions, int distance) {
		Map<String, Collection<Gene>> rtrn=new TreeMap<String, Collection<Gene>>();
		for(String chr: chrSizes.keySet()){
			Collection<Gene> set=new TreeSet<Gene>();
			int size=chrSizes.get(chr);
			Collection<Annotation> regions=excludedRegions.get(chr);
			set.addAll(makeGene(regions, size, chr));
			rtrn.put(chr, set);
		}	
		return rtrn;
	}


	private Collection<Gene> makeGene(Collection<Annotation> excludedRegions, Integer chrSize, String chr) {
		Annotation chrAnn=new BasicAnnotation(chr, 0, chrSize);
		chrAnn.setName(chr);
		if(excludedRegions!=null && excludedRegions.size()>1){ //TODO The size >1 is a hack to make it work with some bug in compound interval
			BasicAnnotation annotation=new BasicAnnotation(excludedRegions);
			//this.excluded.put(chr, annotation);
			Annotation diff=chrAnn.minus(annotation);
			chrAnn=diff;
			//Split diff by distance
			/*Collection<Gene> genes =splitByDistance(diff, distance);
			return genes;*/
		}
		Collection<Gene> genes=new TreeSet<Gene>();
		genes.add(new Gene(chrAnn));
		return genes;
	}


	private Collection<Gene> splitByDistance(Annotation diff, int distance) {
		//Iterate over the annotation
		//If block seperated by greater than distance then split
		Collection<Gene> rtrn=new TreeSet<Gene>();
		Annotation previousBlock=null;
		Collection<Annotation> list=new TreeSet<Annotation>();
		
		for(Annotation block: diff.getBlocks()){
			if(previousBlock!=null){
				//check distance
				int newDistance=block.getStart()-previousBlock.getEnd();
				if(newDistance>distance){
					//split--> make annotation from running collection
					Gene g=new Gene(list);
					rtrn.add(g);
					//make new collection and add this
					list=new TreeSet<Annotation>();
				}
			}
			list.add(block);
			previousBlock=block;
		}
		
		if(!list.isEmpty()){
			Gene g=new Gene(list);
			rtrn.add(g);
		}
		
		logger.info(rtrn.size());
		
		return rtrn;
	}

	/**
	 * Will return the continuous fragment
	 */
	@Override
	public Collection<? extends Window> getFragment(String chr, int start, int end) {
		//This should always return the fragment defined in the continuousSpace
		return this.continuousSpace.getFragment(chr, start, end);
	}

	/**
	 * Will return the overlapping regions in gapped space
	 */
	@Override
	public Collection<? extends Window> getOverlappingRegion(String chr, int start, int end) {
		//This should always return the gapped space
		return this.gappedSpace.getOverlappingRegion(chr, start, end);
	}
	
	/**
	 * Interprets the windowSize based on how it was initialized
	 */
	@Override
	public Iterator<? extends Window> getWindowIterator(int windowSize,	String chr, int start, int end, int overlap) {
		//if continuousWidth is true, use gapped
		if(this.constantWidth){
			return this.gappedSpace.getWindowIterator(windowSize, chr, start, end, overlap);
		}
		else{
			//return new MaskedIterator(this.continuousSpace.getWindowIterator(windowSize, chr, start, end, overlap));
			return new MaskedIterator(new Alignments(chr, start,end), windowSize, overlap);
		}
	}

	@Override
	public Iterator<? extends Window> getWindowIterator(String chr,	int windowSize, int overlap) {
		//if continuousWidth is true, use gapped
		if(this.constantWidth){
			return this.gappedSpace.getWindowIterator(chr, windowSize, overlap);
		}
		else{
			//return new MaskedIterator(this.continuousSpace.getWindowIterator(chr, windowSize, overlap));
			return new MaskedIterator(this.continuousSpace.getEntireChromosome(chr), windowSize, overlap);
		}
	}

	@Override
	public Iterator<? extends Window> getWindowIterator(int windowSize,	int overlap) {
		//if continuousWidth is true, use gapped
		if(this.constantWidth){
			return this.gappedSpace.getWindowIterator(windowSize, overlap);
		}
		else{
			throw new UnsupportedOperationException();
			//return new MaskedIterator(this.continuousSpace.getWindowIterator(windowSize, overlap));
		}
	}

	@Override
	public Iterator<? extends Window> getWindowIterator(Annotation window,	int windowSize, int overlap) {
		return getWindowIterator(windowSize, window.getChr(), window.getStart(), window.getEnd(), overlap);
	}


	/**
	 * Always returns the gapped space length
	 */
	@Override
	public long getLength() {
		return this.gappedSpace.getLength();
	}
	
	@Override
	public long getLength(String chr) {
		return this.gappedSpace.getLength(chr);
	}

	
	/**
	 * Wrap the genomic iterator so it splices out the gap from the interval
	 * @author mguttman
	 */
	private class MaskedIterator implements Iterator<Window> {
		String chr;
		private int windowSize;
		private int step;	
		GeneWindow scanRegion;
		Iterator<? extends Annotation> blockIter;
		Iterator<Window> currentWindows;
		boolean finishedUngapped=false;
		
		/**
		 * Constructs an iterator on chromosome chr starting at start, ending at end, of size, windowSize with overlap between consecutive windows 
		 * @param chr
		 * @param windowSize
		 * @param start
		 * @param overlap
		 * @param end
		 */
		public MaskedIterator(Annotation region, int windowSize, int overlap){
			this.scanRegion=genes.get(region.getChr()).iterator().next().trimAbsolute(region.getStart(), region.getEnd()); //TODO This needs to go through the full list
			this.blockIter=this.scanRegion.getBlocks().iterator();
			this.chr=region.getChr();
			this.windowSize = windowSize;
			this.step = windowSize - overlap;
		}
		
		
		/**
		 * Returns true if there is another window
		 */
		public boolean hasNext() {
			if(currentWindows!=null && currentWindows.hasNext()){
				return true;
			}
			else if(blockIter.hasNext() && !finishedUngapped){
				this.currentWindows=getBlockWindows(blockIter.next());
				return hasNext();
			}
			else if(!blockIter.hasNext()&& !finishedUngapped){
				//starting gapped windows
				finishedUngapped=true;
				//reset blockIter
				blockIter=this.scanRegion.getBlocks().iterator();
				return hasNext();
			}
			else if(blockIter.hasNext() && finishedUngapped){
				this.currentWindows=getNonBlocked(blockIter.next());
				return hasNext();
			}
			return false;
		}

		private Iterator<Window> getNonBlocked(Annotation block) {
			//TODO Should do a dynamic programming thing here
			logger.info("get gapped blocks");
									
			Collection<Window> rtrn=new ArrayList<Window>();
			
			//Get the windows spanning gaps
			int start=block.getEnd()-windowSize;
			int end=block.getEnd();
				
			//TODO Get the block from block.getEnd
			GeneWindow fullRegion=this.scanRegion.trimAbsolute(start, end+windowSize);
			
			logger.info(block.getStart()+" "+block.getEnd()+" "+start+" "+end);
			
			long newTime=0;
			long oldTime=0;
			
			for(int i=start; i<end; i+=step){
				long startTime=System.currentTimeMillis();
				Window w=new GeneWindow(fullRegion);
				w.setStart(i);
				w.setEnd(i+windowSize);
				long endTime=System.currentTimeMillis();
				newTime+=(endTime-startTime);
				
				//startTime=System.currentTimeMillis();
				//w=fullRegion.trimAbsolute(i, i+windowSize);
				//endTime=System.currentTimeMillis();
				//oldTime+=(endTime-startTime);
				
				rtrn.add(w);
			}
			
			logger.info("got "+rtrn.size()+" blocks in "+(oldTime)+" vs "+newTime);
			
			return rtrn.iterator();
		}

		private Iterator<Window> getBlockWindows(Annotation annotation) {
			//create window over block from start--> end-window size
			int start=annotation.getStart();
			int end=annotation.getEnd()-windowSize;
			
			Collection<Window> rtrn=new ArrayList<Window>();
			
			for(int i=start; i<end; i+=step){
				Window w=new GenomeWindow(chr,i, i+windowSize);
				rtrn.add(w);
			}
			
			return rtrn.iterator();
		}

		/**
		 * Returns the next window
		 */
		public Window next() {
			if(hasNext()){
				return this.currentWindows.next();
			}
			else{throw new IllegalStateException("Called next but does not have next. Call hasNext() first");}
		}
		
		

		/**
		 * N/A here
		 */
		public void remove() {
			// 
		}
		
		
	}
	
	/**
	 * Wrap the genomic iterator so it splices out the gap from the interval
	 * @author mguttman
	 */
	private class MaskedIterator3 implements Iterator<Window> {

		Iterator<? extends Window> iter;
		Window next;
		
		MaskedIterator3(Iterator<? extends Window> iter){
			this.iter=iter;
			hasNext(); //To initialize
		}
		
		@Override
		public boolean hasNext() {
			//we will iterate through until we find a Window that passes filter
			//if we hit the end return false
			//else save the alignment and return true (cache fact that we have an unreturned element saved)
			if(next!=null){return true;}
			else{
				//test for a next element and cache result
				while(iter.hasNext()){
					Window alignment=getNext();
					if(alignment!=null){
						next=alignment;
						return true;
					}
				}
				return false;
			}
		}
		
		private Window getNext(){
			Window alignment=iter.next();
			Window trim=trim(alignment);
			return trim;
		}

		private Window trim(Window w) {
			Collection<? extends Window> list=gappedSpace.getOverlappingRegion(w.getChr(), w.getStart(), w.getEnd());
			if(list!=null && !list.isEmpty()){
				Window rtrn=list.iterator().next();
				if(rtrn!=null && rtrn.getBlocks().size()>0){
					return rtrn;
				}
			}
			return null;
			
			//intersect with gapped map
			/*Gene nonGaps=nonMaskedRegions.get(w.getChr()).iterator().next();
			GeneWindow rtrn=new GeneWindow(nonGaps.intersect(w));
			rtrn.addSourceAnnotation(w);
						
			if(rtrn.getBlocks()!=null && rtrn.getBlocks().size()>0){
				System.out.println(rtrn.toBED());
				return rtrn;
			}
			return null;*/
		}

		@Override
		public Window next() {
			//return next
			Window previousNext=next;
			next=null;
			//call hasNext
			hasNext();
			return previousNext;
		}

		@Override
		public void remove() {
			// TODO Auto-generated method stub
			
		}
		
	}

	@Override
	public Annotation permuteAnnotation(Annotation a) {
		throw new UnsupportedOperationException("TODO");
	}
	
	@Override
	public Annotation permuteAnnotation(Annotation a, Annotation bounds) {
		throw new UnsupportedOperationException("TODO");
	}
	
	@Override
	public Collection<String> getReferenceNames() {
		return gappedSpace.getReferenceNames();
	}
	
	@Override
	public Collection<? extends Annotation> getReferenceAnnotations() {
		throw new UnsupportedOperationException("TODO");
	}
	
	@Override
	public Annotation getReferenceAnnotation(String name) {
		//System.out.println(gappedSpace.getReferenceAnnotation(name).size());
		return gappedSpace.getReferenceAnnotation(name);
		//TODO: Confirm this is correct
	}
	
	@Override
	public Iterator<Window> getWindowIterator(Collection<Gene> baseGenes, int windowSize, int overlap) {
		throw new UnsupportedOperationException("TODO");
	}

	@Override
	public Collection<? extends Window> getFragment(Annotation annotation) {
		return getFragment(annotation.getChr(), annotation.getStart(), annotation.getEnd());
	}

	@Override
	public int getSize(Annotation region) {
		throw new UnsupportedOperationException("TODO");
	}

	@Override
	public Collection<String> getChromosomeNames() {
		throw new UnsupportedOperationException("TODO");
	}

	public boolean hasChromosome(String chr) {
		return genes.containsKey(chr);
	}
	
	@Override
	public Annotation getEntireChromosome(String chrName) {
		throw new UnsupportedOperationException("TODO");
	}
	
	@Override
	public boolean isValidWindow(Annotation window) {
		return gappedSpace.isValidWindow(window);
	}
}
