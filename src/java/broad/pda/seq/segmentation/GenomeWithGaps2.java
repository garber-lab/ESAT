package broad.pda.seq.segmentation;

import java.util.Collection;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Gene;

import org.apache.log4j.Logger;

import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.core.util.CollapseByIntersection;
import broad.pda.datastructures.Alignments;


// TODO integrate with TranscriptomeSpace
public class GenomeWithGaps2 {

	private IntervalTree<Annotation> cDNASpace; //indexed by relative position
	private int genomeSize;
	private boolean directOrientation;
	static Logger logger = Logger.getLogger(GenomeWithGaps2.class.getName());
	
	public GenomeWithGaps2(Collection<Annotation> gaps, String chr, int genomeSize, String orientation){
		gaps=CollapseByIntersection.collapseByIntersection(gaps, true);
		cDNASpace=convertToCovered(gaps, chr, genomeSize);
		this.genomeSize = genomeSize;
		this.directOrientation = "+".equals(orientation);
	}
	
	
	//Need to convert this relative space indexes
	private IntervalTree<Annotation> convertToCovered(Collection<Annotation> gaps, String chr, int genomeSize){
		IntervalTree<Annotation> rtrn=new IntervalTree();
		
		Map<String, IntervalTree<Annotation>> trees=makeIntervalTree(gaps);
		IntervalTree<Annotation> gapTree=trees.get(chr);
		
		int currentBP=0;
		if(gapTree != null) {
			Iterator<Node<Annotation>> iter=gapTree.iterator();
			
			int counter=0;
			while(iter.hasNext()){
				Annotation gap=iter.next().getValue();
				Alignments block=new Alignments(gap.getChr(), currentBP, gap.getStart());
				currentBP=gap.getEnd();
				rtrn.put(counter, counter+block.getSize(), block);
				counter=counter+block.getSize();
			}
			
			Alignments end=new Alignments(chr, currentBP, genomeSize);
			if(end.getSize()>0){rtrn.put(counter, counter+end.getSize(), end);}
		} else {
			logger.debug("No gaps tree for chr " + chr);
		}
		return rtrn;
	}
	
	private Map<String, IntervalTree<Annotation>> makeIntervalTree(Collection<Annotation> alignments){
		Map<String, IntervalTree<Annotation>> rtrn=new TreeMap<String, IntervalTree<Annotation>>();
		
		for(Annotation align: alignments){
			IntervalTree<Annotation> tree=new IntervalTree<Annotation>();
			String chr=align.getChr();
			if(rtrn.containsKey(chr)){tree=rtrn.get(chr);}
			tree.put(align.getStart(), align.getEnd(), align);
			rtrn.put(align.getChr(), tree);
		}
		
		return rtrn;
	}
	
	public Gene getRelativeWindow(int relativeStart, int relativeEnd){
		//System.err.println("Calculating relative window for  " + relativeStart  + " - " + relativeEnd + " region length " + genomeSize);
		if(relativeEnd > genomeSize) {
			System.err.println("ERROR " + relativeStart  + " - " + relativeEnd + " have no overlappers in this CDNA of size " + (genomeSize));
			return null;  //TODO: Should throw an illegal argument exception since this means that coordinates are outside of range
		}
	
		Iterator<Node<Annotation>> relativeMapping=cDNASpace.overlappers(relativeStart, relativeEnd);
		Collection<Annotation> exons=new TreeSet<Annotation>();
		//if(!relativeMapping.hasNext()) { System.out.println("start-end: "+ relativeStart+"-"+relativeEnd+" had no exons mapping");}
		while(relativeMapping.hasNext()){
			Node<Annotation> node=relativeMapping.next();
			
			//if rs and re in 1 block then done. Might be able to tell this case from num overlappers
			if(node.getStart()<relativeStart && node.getEnd()>relativeEnd){
				Annotation absAlignment=node.getValue();
				int absStart=relativeStart-node.getStart();
				int absEnd=relativeEnd-node.getStart();
				Alignments firstExon=new Alignments(absAlignment.getChr(), absAlignment.getStart()+absStart, absAlignment.getStart()+absEnd); //end unless fully contained
				exons.add(firstExon);
			}
			
			//if rs is within the block start at the appropriate point
			else if(node.getStart()<relativeStart){
				Annotation absAlignment=node.getValue();
				int absStart=relativeStart-node.getStart();
				Annotation firstExon=new Alignments(absAlignment.getChr(), absAlignment.getStart()+absStart, absAlignment.getEnd()); //end unless fully contained
				exons.add(firstExon);
			}
			
			else if(node.getEnd()>relativeEnd){
				Annotation absAlignment=node.getValue();
				int absEnd=relativeEnd-node.getStart();
				Annotation lastExon=new Alignments(absAlignment.getChr(), absAlignment.getStart(), absAlignment.getStart()+absEnd); //end unless fully contained
				exons.add(lastExon);
			}
			else{
				Annotation absAlignment=node.getValue();
				exons.add(absAlignment);
			}
		}
		
		Gene gene=null;
		if(!exons.isEmpty()){gene=new Gene(exons);}
		
		return gene;
	}
	
	public int getRelativeGenomeLength(){
		return this.cDNASpace.max().getEnd();
	}
	
}
