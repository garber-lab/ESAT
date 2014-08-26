package broad.pda.seq.segmentation;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Gene;

import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.core.util.CollapseByIntersection;
import broad.pda.datastructures.Alignments;

//TODO integrate with TranscriptomeSpace
@Deprecated
public class GenomeWithGaps {

	IntervalTree<Annotation> cDNASpace;
	IntervalTree<Annotation> gapSpace;
	Gene currentInterval;
	Gene nextInterval;
	int windowSize;
	Collection<Annotation> gaps;
	
	public GenomeWithGaps(Collection<Alignments> gaps, String chr, int genomeSize, int windowSize, int startPosition){
		System.err.println(chr);
		
		//first need to collapse gaps into largest overlapping gap
		this.gaps=CollapseByIntersection.collapseByIntersection(gaps, true);
		
		
		cDNASpace=convertToCovered(this.gaps, chr, genomeSize);
		this.windowSize=windowSize;
		Alignments align=new Alignments(chr, startPosition, startPosition);
		this.currentInterval=new Gene(align);
		this.nextInterval=peakNextWindow();
		System.err.println("First tested interval "+this.nextInterval.getAlignment().toUCSC());
		
		
	}
	
	private void write(String save, Collection<Alignments> gaps)throws IOException{
		FileWriter writer=new FileWriter(save);
		
		for(Alignments gap: gaps){writer.write(gap+"\n");}
		
		writer.close();
	}
	
	public Collection<Annotation> getGaps(){return this.gaps;}
	public IntervalTree<Annotation> getCDNAs(){return this.cDNASpace;}
		
	//is there an interval of size windowSize thats still available in the current genome?
	public boolean hasNext(){
		if(nextInterval.getGappedSize()<windowSize){return false;}
		return true;
	}
	
	private IntervalTree<Annotation> convertToCovered(Collection<Annotation> gaps, String chr, int genomeSize){
		IntervalTree<Annotation> rtrn=new IntervalTree();
		IntervalTree<Annotation> gapTree=makeIntervalTree(gaps).get(chr);
		
		int currentBP=0;
		
		Iterator<Node<Annotation>> iter=gapTree.iterator();
		while(iter.hasNext()){
			Annotation gap=iter.next().getValue();
			Annotation block=new Alignments(gap.getChr(), currentBP, gap.getStart());
			currentBP=gap.getEnd();
			if(block.getStart()>block.getEnd()){System.err.println(block);}
			rtrn.put(block.getStart(), block.getEnd(), block);
			
		}
		
		Alignments end=new Alignments(chr, currentBP, genomeSize);
		if(end.getSize()>0){rtrn.put(end.getStart(), end.getEnd(), end);}
		
		return rtrn;
	}
	
	
	private Map<String, IntervalTree> makeIntervalTree(Collection<Annotation> alignments){
		Map<String, IntervalTree> rtrn=new TreeMap();
		
		for(Annotation align: alignments){
			IntervalTree tree=new IntervalTree();
			String chr=align.getChr();
			if(rtrn.containsKey(chr)){tree=rtrn.get(chr);}
			tree.put(align.getStart(), align.getEnd(), align);
			rtrn.put(align.getChr(), tree);
		}
		
		return rtrn;
	}
	
	public Gene getNextWindow(){
		currentInterval=nextInterval;
		nextInterval=peakNextWindow();
		return currentInterval;
	}
	
	private Gene peakNextWindow(){
		
		int start=currentInterval.getStart()+1;
		int originalStart=start;
		Iterator<Node<Annotation>> iter=cDNASpace.iterator();
		if(currentInterval.getStart()>0){iter=cDNASpace.iterator(start, start);}
		
		
		//start with region overlapping it
		Iterator<Node<Annotation>> overlappers=cDNASpace.overlappers(start, start+1);
		
		int width=windowSize;
		boolean done=false;
		Collection<Annotation> exonBlocks=new TreeSet();
		while((overlappers.hasNext() || iter.hasNext()) && !done){
			Annotation align=null;
			if(overlappers.hasNext()){align=overlappers.next().getValue();}
			else{align=iter.next().getValue();}
			int startPos=Math.max(start, align.getStart());
			int endPos =Math.min(align.getEnd(), startPos+width);
			Alignments block=new Alignments(align.getChr(), startPos, endPos);
			width=width-block.getSize();
			start=block.getEnd();
			exonBlocks.add(block);
			if(width<=0){done=true;}
			start=endPos;
		}
		
		Gene rtrn=new Gene(exonBlocks);
		
		return rtrn;
	}
	
	public Gene getNextWindow(int position){
		
		int start=position;
		Iterator<Node<Annotation>> iter=cDNASpace.iterator(start, start);
		
		
		//start with region overlapping it
		Iterator<Node<Annotation>> overlappers=cDNASpace.overlappers(start, start+this.windowSize);
		
		int width=windowSize;
		boolean done=false;
		Collection<Annotation> exonBlocks=new TreeSet();
		while((overlappers.hasNext() || iter.hasNext()) && !done){
			Annotation align=null;
			if(overlappers.hasNext()){align=overlappers.next().getValue();}
			else{align=iter.next().getValue();}
			int startPos=Math.max(start, align.getStart());
			int endPos =Math.min(align.getEnd(), startPos+width);
			Alignments block=new Alignments(align.getChr(), startPos, endPos);
			if(startPos>endPos){System.err.println("ERROR");}
			width=width-block.getSize();
			start=block.getEnd();
			exonBlocks.add(block);
			if(width<=0){done=true;}
			start=endPos;
		}
		
		Gene rtrn=new Gene(exonBlocks);
		
		return rtrn;
	}
	
	
}
