package broad.core.util;

import java.io.*;
import java.util.*;

import umms.core.annotation.Annotation;
import umms.core.annotation.Annotation.Strand;
import umms.core.annotation.BEDFileParser;
import umms.core.annotation.BasicAnnotation;
import umms.core.annotation.Gene;
import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;



public class CollapseByIntersection {
	
	static int numOverlap;
	static int numOverlapGenes;
	private static final int MAX_RECURSION = 50;
	
	public static Set<Annotation> collapseByIntersection(Collection<? extends Annotation> BasicAnnotation, boolean intersection){
		numOverlap=10; //just to start
		
		Set<Annotation> rtrn=new TreeSet<Annotation>(BasicAnnotation);
		
		int i=0;
		while(numOverlap>1){ //TODO: do not use a static variable for this. 
			numOverlap=0;
			Set<Annotation> set=new TreeSet<Annotation>();
			Map<String, IntervalTree<Annotation>> trees=makeIntervalTree(rtrn);
			//for each region get all overlappign and collapse by intersection
			for(Annotation align: rtrn){
				Iterator<Node<Annotation>> regions=getOverlapping(align, trees);
				Annotation collapsed=collapse(regions, intersection);
				if(collapsed.getSize()>0 && collapsed.getStart()<collapsed.getEnd()){set.add(collapsed);}
			}
			i++;
			rtrn=set;
			//System.err.println("Iteration "+i+" Num Overlap "+numOverlap);
		}
		
		return rtrn;
	}
	
	
	/**
	 * This method returns only exons that are trimmed down by intron location. If its fully within intron it just returns it. If bases overlap intron it trims down
	 * @param exons
	 * @param introns
	 * @param strand
	 * @return
	 */
	public static Collection<Annotation> DecollapseByIntronLocation(Collection<Annotation> exons, Collection<Annotation> introns, String strand){
		if(introns==null || introns.isEmpty()){return exons;}
		
		introns=filter(introns, strand);
		
		TreeSet<Annotation> decollapsedExons=new TreeSet<Annotation>();
		
		Map<String, IntervalTree<Annotation>> intronTree=makeIntervalTree(introns);
		
		for(Annotation exon: exons){			
			//There is a reason for the -1 and +1. If an intron abuts but does not overlap the exon overlaps an intron that truly overlaps the exon,
			//the result of the decollapse will not contain the segment that abuts the original exon. The -1 and +1 ensure that butting introns are
			//included and thus exons that abut them are reported.
			Iterator<Node<Annotation>> iter=intronTree.get(exon.getChr()).overlappers(exon.getStart() - 1, exon.getEnd() + 1); 
			if(!iter.hasNext()) {
				decollapsedExons.add(exon);
			} else {
				LinkedList<Annotation> overlappingIntrons = new LinkedList<Annotation>(); 
				while(iter.hasNext()){
					overlappingIntrons.add(iter.next().getValue());
				}
				//Collections.sort(overlappingIntrons);
				//TreeSet<BasicAnnotation> diff = new TreeSet<BasicAnnotation>();
				//getRecursiveDifference(exon, overlappingIntrons, rtrn, 1);
				carveExon(exon, overlappingIntrons, decollapsedExons);
				//rtrn.addAll(diff);				
			}
		}
		
		// We are now going to filter out any exon that fully contains any exon-intron-exon combination. These affect highly expressed genes producing large numbers of intronic reads.
		Map<String, IntervalTree<Annotation>> exonTree = makeIntervalTree(decollapsedExons);
		Set<Annotation> filteredExons = new TreeSet<Annotation>();
		
		for(Annotation exon : decollapsedExons) {
			Iterator<Node<Annotation>> intronIter=intronTree.get(exon.getChr()).overlappers(exon.getStart(), exon.getEnd());
			boolean keep = true;
			while(intronIter.hasNext() && keep) {
				Annotation intron = intronIter.next().getValue();
				if(exon.fullyContains(intron)) {
					BasicAnnotation intronLeftPoint = new BasicAnnotation(intron.getChr(), intron.getStart()-1, intron.getStart());
					BasicAnnotation intronRightPoint = new BasicAnnotation(intron.getChr(), intron.getEnd(), intron.getEnd()+1);						
					Iterator<Node<Annotation>> leftExons=exonTree.get(exon.getChr()).overlappers(intronLeftPoint.getStart(), intronLeftPoint.getEnd());
					Iterator<Node<Annotation>> rightExons=exonTree.get(exon.getChr()).overlappers(intronRightPoint.getStart(), intronRightPoint.getEnd());
					while(leftExons.hasNext() && keep) {
						Annotation leftExon = leftExons.next().getValue();
						if(leftExon.equals(exon)) {
							continue;
						}
						while(rightExons.hasNext() && keep) {
							Annotation rightExon = rightExons.next().getValue();
							if(rightExon.equals(exon)) {
								continue;
							}
							keep = exon.fullyContains(rightExon) && exon.fullyContains(leftExon);
							
						}
					}
				}
			}
	
			if(keep) {
				filteredExons.add(exon);
			}
		}
		
		return filteredExons;
	}
	
	public static Collection<Annotation> filter(Collection<Annotation> introns, String strand) {
		if(strand.equalsIgnoreCase("+")){
			Collection<Annotation> rtrn=new TreeSet<Annotation>();
			for(Annotation intron: introns){if(intron.getOrientation().equals(Strand.POSITIVE)){rtrn.add(intron);}}
			return rtrn;
		}
		else if(strand.equalsIgnoreCase("-")){
			Collection<Annotation> rtrn=new TreeSet<Annotation>();
			for(Annotation intron: introns){if(intron.getOrientation().equals(Strand.NEGATIVE)){rtrn.add(intron);}}
			return rtrn;
		}
		else{return introns;}
	}

	public static Collection<Annotation> DecollapseByIntronLocation(Collection<Annotation> exons, Collection<Annotation> introns){
		return DecollapseByIntronLocation(exons, introns, "*");
	}
	
	
	private static void carveExon(Annotation exon, LinkedList<Annotation> overlappingIntrons, TreeSet<Annotation> difference) {
		if( difference.contains(exon)) {
			return;
		}
		//System.err.println("Carving " + exon.toUCSC());
		if(overlappingIntrons.size() == 0 ) {
			difference.add(exon);
		} else {
			TreeSet<BasicAnnotation> intronLeftEnds = new TreeSet<BasicAnnotation>();
			TreeSet<BasicAnnotation> intronRightEnds = new TreeSet<BasicAnnotation>();
			BasicAnnotation enlargedExon = new BasicAnnotation(exon.getChr(), exon.getStart()-1,exon.getEnd()+1);
		
			for(Annotation intron : overlappingIntrons){
				/*if(exon.fullyContained(intron)){
					BasicAnnotation leftEnd = new BasicAnnotation(intron.getChr(), intron.getStart()-1, intron.getStart());
					BasicAnnotation rightEnd = new BasicAnnotation(intron.getChr(), intron.getEnd(), intron.getEnd()+1);
					intronLeftEnds.add(leftEnd);
					intronRightEnds.add(rightEnd);
				} else {
					BasicAnnotation leftEnd = new BasicAnnotation(intron.getChr(), intron.getStart(), intron.getStart()+1);
					BasicAnnotation rightEnd = new BasicAnnotation(intron.getChr(), intron.getEnd() - 1, intron.getEnd());
					if(enlargedExon.overlaps(leftEnd)) {
						intronLeftEnds.add(leftEnd);
					}
					if(enlargedExon.overlaps(rightEnd)) {
						intronRightEnds.add(rightEnd);
					} 
				}*/
				
				BasicAnnotation leftEnd = new BasicAnnotation(intron.getChr(), intron.getStart(), intron.getStart()+1);
				BasicAnnotation rightEnd = new BasicAnnotation(intron.getChr(), intron.getEnd() - 1, intron.getEnd());
				if(enlargedExon.overlaps(leftEnd)) {
					intronLeftEnds.add(leftEnd);
				}
				if(enlargedExon.overlaps(rightEnd)) {
					intronRightEnds.add(rightEnd);
				}
			}
				
			
			if(intronRightEnds.isEmpty()) {
				intronRightEnds.add(new BasicAnnotation(exon.getChr(), exon.getStart()-1, exon.getStart()));
			}
			if(intronLeftEnds.isEmpty()) {
				intronLeftEnds.add(new BasicAnnotation(exon.getChr(), exon.getEnd(), exon.getEnd()+1));
			}

			for(BasicAnnotation re : intronRightEnds) {
				for(BasicAnnotation le : intronLeftEnds) {
					if(re.getEnd() < le.getStart()) {
						Annotation carving = new BasicAnnotation(exon.getChr(), re.getEnd(),le.getStart());
						difference.add(carving);
					} else {
						// carve out the middle of the exon
						if(le.getEnd() > exon.getStart()) {
							BasicAnnotation leftPiece = new BasicAnnotation (exon.getChr(), exon.getStart(), Math.min(exon.getEnd(),le.getStart()));//new BasicAnnotation (exon.getChr(), exon.getStart(), Math.min(exon.getEnd(),le.getEnd()+1));
							difference.add(leftPiece);
						}
						if(re.getStart() < exon.getEnd()) {
							BasicAnnotation rightPiece = new BasicAnnotation (exon.getChr(), re.getEnd(), exon.getEnd());//new BasicAnnotation (exon.getChr(), re.getStart()+1, exon.getEnd());
							difference.add(rightPiece);
						}
					}
				}
			}
			
		}
		

	}
	
	private static void getRecursiveDifference(Annotation exon, LinkedList<Annotation> overlappingIntrons, TreeSet<Annotation> difference, int recursionNum) {
		//TreeSet<BasicAnnotation> rtrn = new TreeSet<BasicAnnotation>();
		//System.err.println("Recursion: " + recursionNum + ", Exon " + exon.toUCSC() + ", #introns: " + overlappingIntrons.size());
		Annotation lastIntron = null;
		if( difference.contains(exon)) {
			return;
		}
		if(overlappingIntrons.size() == 0 || recursionNum >= MAX_RECURSION ) {
			difference.add(exon);
			if(recursionNum >= MAX_RECURSION) {
				//System.err.println("Maximum number of recursions hit (exon "+exon.toUCSC() +"). Graph will get too complex, this is usually the result of alingment artifacts.");
			}
		} else {
			LinkedList<Annotation> usedIntrons = new LinkedList<Annotation>();
			while(overlappingIntrons.size() > 0) {
				Annotation intron = overlappingIntrons.pop();
				//System.err.println("\tIntron " + intron.toUCSC());
				if(lastIntron != null && !intron.overlaps(lastIntron)) {
					break;
				}
				Collection<Annotation> exonIntronDiff = takeDifference(exon, intron);
				//System.err.println("\texonIntrondiff = " + exonIntronDiff);
				
				for(Annotation exonPart : exonIntronDiff) {
					LinkedList<Annotation> nonOverlappingIntrons = new LinkedList<Annotation>();
					for(Annotation otherIntron : overlappingIntrons) {
						if(!otherIntron.overlaps(intron) && otherIntron.overlaps(exonPart)) {
							nonOverlappingIntrons.add(otherIntron);
						}
					}
					for(Annotation usedIntron : usedIntrons) {
						if(!usedIntron.overlaps(intron) && intron.overlaps(exonPart)) {
							nonOverlappingIntrons.add(usedIntron);
						}
					}
					Collections.sort(nonOverlappingIntrons);
					//System.err.println("\t\tnonoverlapping introns" + nonOverlappingIntrons);
					if(!difference.contains(exonPart)){
						getRecursiveDifference(exonPart, nonOverlappingIntrons, difference , recursionNum++);
					}
				}
				lastIntron = intron;
				usedIntrons.push(intron);
			}
		}
		//System.err.println("Done with " + exon.toUCSC());
		
		//return rtrn;
	}

	/*
	 * MG: Fixed what must lead to a bug in the cases where a differences returns more than one annotation 7/24/17
	 */

	private static Collection<Annotation> takeDifference(Annotation exon,  Annotation intron) {
		Collection<Annotation> rtrn = new ArrayList<Annotation>();
		if(!exon.overlaps(intron)){
			rtrn.add(exon);
		}else if(!exon.fullyContains(intron)) {
			Annotation newExons=exon.minus(intron);
			rtrn.addAll(newExons.getBlocks());
		}
		
		return rtrn;
	}


	public static Collection<Annotation> DecollapseByIntronLocationOld(Collection<Annotation> exons, Collection<Annotation> introns){
		Collection<Annotation> rtrn=new TreeSet<Annotation>();
		
		Map<String, IntervalTree<Annotation>> intronTree=makeIntervalTree(introns);
		
		for(Annotation exon: exons){
			rtrn.add(exon);
			Iterator<Node<Annotation>> iter=intronTree.get(exon.getChr()).overlappers(exon.getStart(), exon.getEnd());
			while(iter.hasNext()){
				Annotation intron=iter.next().getValue();
				rtrn.addAll(takeDifference(exon, intron));	
			}
		}
		
		return rtrn;
	}
	
	public static Collection<Annotation> DecollapseByIntronLocationWorking(Collection<Annotation> exons, Collection<Annotation> introns){
		if(introns==null || introns.isEmpty()){return exons;}
		Collection<Annotation> rtrn=new TreeSet<Annotation>();
		
		Map<String, IntervalTree<Annotation>> intronTree=makeIntervalTree(introns);
		
		for(Annotation exon: exons){
			rtrn.add(exon);
			Iterator<Node<Annotation>> iter=intronTree.get(exon.getChr()).overlappers(exon.getStart(), exon.getEnd());
			while(iter.hasNext()){
				Annotation intron=iter.next().getValue();
				//MG: Added this case because crashed when an exon overlapped an intron completely
				// Removed this hack 7/24/17
				if(exon.fullyContains(intron)){
					//System.err.println("Odd case Exon "+exon.toUCSC()+" Intron: "+intron.toUCSC());
					Annotation newExons=exon.minus(intron);
					rtrn.addAll(newExons.getBlocks());
				}
				
			}
		}
		
		return rtrn;
	}
	
	
	/*public static Set<BasicAnnotation> CollapseByUnion(Collection<BasicAnnotation> BasicAnnotation){
		
		Collection<BasicAnnotation> accountedFor=new TreeSet();
		Set set=new TreeSet();
		Map<String, IntervalTree> trees=makeIntervalTree(BasicAnnotation);
		//for each region get all overlappign and collapse by intersection
		for(BasicAnnotation align: BasicAnnotation){
			if(!accountedFor.contains(align)){
			Collection<BasicAnnotation> regions=getOverlappingAll(align, trees);
			accountedFor.addAll(regions);
			BasicAnnotation collapsed=collapse(regions);
			if(collapsed.getSize()>0){set.add(collapsed);}
			}
		}
		
		
		return set;
	}*/
	
	
	public static Collection<Gene> updateBoundariesWithoutCrossingIntrons(Collection<Gene> genes){
		numOverlapGenes=10; //just to start
		
		Collection<Gene> rtrn=genes;
		
		int i=0;
		while(numOverlapGenes>1){
			numOverlapGenes=0;
			Set set=new TreeSet();
			Map<String, IntervalTree<Gene>> trees=makeIntervalTreeForGenes(rtrn);
			for(Gene align: rtrn){
				Iterator<Node<Gene>> regions=getOverlapping(align, trees);
				Gene collapsed=collapseGenes(regions);
				if(collapsed!=null){set.add(collapsed);}
			}
			rtrn=set;
			i++;
			//System.err.println("Iteration "+i+" Num Overlap "+numOverlapGenes);
		}
		
		return rtrn;
	}
	
		
	private static Gene collapseGenes(Iterator<Node<Gene>> regions){
		//these are all regions overlapping a given interval
		
		//now take all exons and collapse
		int i=0;
		Collection<Annotation> exons=new TreeSet();
		while(regions.hasNext()){
			i++;
			Gene gene=regions.next().getValue();
			exons.addAll(gene.getExonSet());
		}
		numOverlapGenes=Math.max(i, numOverlapGenes);
		
		exons=collapseByIntersection(exons, false);
		
		//then make into transcript
		Gene gene=null;
		if(exons!=null && !exons.isEmpty() && exons.size()>0){
			gene=new Gene(exons);
			//System.err.println(gene.getAlignment().toUCSC());
		}
		
		
		
		return gene;
	}
	
	private static BasicAnnotation trimExon(Annotation exon, Annotation intron){
		//if exon starts at or after intron and ends at or before intron -> return null
		if(intron.fullyContains(exon)){return null;}
		
		//if exon starts before intron
		if(exon.getStart()<intron.getStart()){
			return new BasicAnnotation(exon.getChr(), exon.getStart(), Math.min(exon.getEnd(), intron.getStart()));
		}
		
		//if exon starts after intron
		if(exon.getEnd()>intron.getEnd()){
			return new BasicAnnotation (exon.getChr(), Math.max(exon.getStart(), intron.getEnd()), exon.getEnd());
		}
		
		return null;
		
	}
	
	private static Iterator<Node<Annotation>> getOverlapping(Annotation align, Map<String, IntervalTree<Annotation>> trees){
		IntervalTree<Annotation> tree=trees.get(align.getChr());
		Iterator<Node<Annotation>> iter=tree.overlappers(align.getStart(), align.getEnd());
		return iter;
	}
	
	//can loop through everything in the overlapping list and get all of their overlaps as well
	private static Collection<Annotation> getOverlappingAll(Annotation align, Map<String, IntervalTree> trees){
		IntervalTree<Annotation> tree=trees.get(align.getChr());
		
		Collection<Annotation> rtrn=new TreeSet<Annotation>();
		
		//go through and get overlapping
		Iterator<Node<Annotation>> iter=tree.overlappers(align.getStart(), align.getEnd());
		rtrn=addAll(rtrn, iter);
		
		while(iter.hasNext()){
			Node<Annotation> node=iter.next();
			Iterator overlappers=tree.overlappers(node.getStart(), node.getEnd());
			rtrn=addAll(rtrn, overlappers);
		}
		
		return rtrn;
	}
	
	private static Collection<Annotation> addAll(Collection<Annotation> set, Iterator<Node<Annotation>> overlappers){
		Collection<Annotation> rtrn=set;
		while(overlappers.hasNext()){
			rtrn.add(overlappers.next().getValue());
		}
		return rtrn;
	}
	
	private static Iterator<Node<Gene>> getOverlapping(Gene gene, Map<String, IntervalTree<Gene>> trees){
		IntervalTree<Gene> tree=trees.get(gene.getChr());
		Iterator<Node<Gene>> iter=tree.overlappers(gene.getAlignment().getStart(), gene.getAlignment().getEnd());
		return iter;
	}
	
	
	
	
	private static Annotation collapse(Iterator<Node<Annotation>> iter, boolean intersection){
		int start=-1;
		int end=Integer.MAX_VALUE;
		String chr="";
		
		int i=0;
		while(iter.hasNext()){
			Annotation align=iter.next().getValue();
			if(i==0){start=align.getStart(); end=align.getEnd();}
			if(intersection){
				start=Math.max(start, align.getStart());
				end=Math.min(end,  align.getEnd());
			}
			else{
				start=Math.min(start, align.getStart());
				end=Math.max(end,  align.getEnd());
			}
			chr=align.getChr();
			i++;
		}
		
		numOverlap=Math.max(i, numOverlap); //TODO: This is not good. Why modify this class variable here? It makes it horribly non-thread safe for no reason.
		Annotation align=new BasicAnnotation(chr, start, end);
		return align;
	}
	
	/*private static BasicAnnotation collapse(Collection<BasicAnnotation> set){
		int start=-1;
		int end=Integer.MAX_VALUE;
		String chr="";
		
		int i=0;
		for(BasicAnnotation align: set){
			if(i==0){start=align.getStart(); end=align.getEnd();}
			start=Math.min(start, align.getStart());
			end=Math.max(end,  align.getEnd());
			chr=align.getChr();
			i++;
		}
		
		numOverlap=Math.max(i, numOverlap);
		BasicAnnotation align=new BasicAnnotation(chr, start, end);
		return align;
	}*/
	
	private static void write(String save, Collection set)throws IOException{
		FileWriter writer=new FileWriter(save);
		
		for(Object k: set){writer.write(k.toString()+"\n");}
		
		writer.close();
	}
	
	public static Map<String, IntervalTree<Annotation>> makeIntervalTree(Collection<? extends Annotation> BasicAnnotation){
		Map<String, IntervalTree<Annotation>> rtrn=new TreeMap<String, IntervalTree<Annotation>> ();
		
		for(Annotation align: BasicAnnotation){
			IntervalTree<Annotation> tree=new IntervalTree<Annotation>();
			String chr=align.getChr();
			if(rtrn.containsKey(chr)){tree=rtrn.get(chr);}
			//if(align.getStart()>=align.getEnd()){System.err.println("ERROR: " +align.toUCSC());}
			//System.err.println(align);
			tree.put(align.getStart(), align.getEnd(), align);
			rtrn.put(align.getChr(), tree);
		}
		
		return rtrn;
	}
	
	
	public static Map<String, IntervalTree<Gene>> makeIntervalTreeForGenes(Collection<Gene> BasicAnnotation){
		Map<String, IntervalTree<Gene>> rtrn=new TreeMap<String, IntervalTree<Gene>>();
		
		for(Gene align: BasicAnnotation){
			IntervalTree<Gene> tree=new IntervalTree<Gene>();
			String chr=align.getChr();
			if(rtrn.containsKey(chr)){tree=rtrn.get(chr);}
			Node<Gene> node = tree.find(align.getAlignment().getStart(), align.getAlignment().getEnd()+1);
			if (node != null)
				node.incrementCount();
			else
				tree.put(align.getAlignment().getStart(), align.getAlignment().getEnd()+1, align);
			rtrn.put(align.getChr(), tree);
		}
		
		return rtrn;
	}
	
	

	public static Collection<Annotation> CollapseGenesByIntersection(Collection<Gene> temp, boolean intersection) {
		Collection<Annotation> BasicAnnotation=convert(temp);
		
		return collapseByIntersection(BasicAnnotation, intersection);
	}


	private static Collection<Annotation> convert(Collection<Gene> temp) {
		Collection<Annotation> rtrn=new TreeSet();
		
		for(Gene gene: temp){rtrn.add(gene.getAlignment());}
		
		return rtrn;
	}

	public static Map<String, IntervalTree<Annotation>> makeIntervalTreeForGeneExons(Collection<Gene> genes) {
		Collection<Annotation> exons=new TreeSet();
		
		for(Gene gene: genes){
			exons.addAll(gene.getSortedAndUniqueExons());
		}
		
		return makeIntervalTree(exons);
	}




	
	
}
