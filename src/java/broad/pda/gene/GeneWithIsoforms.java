package broad.pda.gene;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import nextgen.core.alignment.Alignment;
import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Gene;

import org.omg.CORBA.portable.RemarshalException;

import broad.core.annotation.GFF;
import broad.core.annotation.LightweightGenomicAnnotation;
import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.core.util.CollapseByIntersection;
import broad.pda.datastructures.Alignments;


/* This class was generated to resolve the problem of treating isoforms with the same exact
stat-end coordinates while storing genes in a sorted manner in an intervalTree.
The class extends RefSeqGene.
The fields inherited by RefSeqGene will represent the isoform with the longest transcript.
Isoforms may be with different orientation
 */
public class GeneWithIsoforms extends Gene{

	Collection <Gene> isoforms;
	int numOfIsoforms;
	
	
	//Ctor:
	/*public RefSeqGeneWithIsoforms(Alignment align) {
		super(align);
		this.numOfIsoforms=1;  this.isoforms=new ArrayList <Gene>();
	}*/
	public GeneWithIsoforms(Collection<? extends LightweightGenomicAnnotation> exons){
		super(exons);
		this.numOfIsoforms=1; this.isoforms=new ArrayList <Gene>();
	}
	public GeneWithIsoforms(LightweightGenomicAnnotation align){
		super(align);
		this.numOfIsoforms=1; this.isoforms=new ArrayList <Gene>();
	}
	public GeneWithIsoforms(String pslString){
		super(pslString);
		this.numOfIsoforms=1; this.isoforms=new ArrayList <Gene>();
	}
	public GeneWithIsoforms(String chrString, int start, int end){
		super(chrString, start,end);
		this.numOfIsoforms=1; this.isoforms=new ArrayList <Gene>();
	}
	
	@Deprecated  // start/end and exons are redundant
	public GeneWithIsoforms(String chr, int start, int end, String name, Strand orientation, Collection<Annotation> exons){
		super(chr,name,orientation,exons);
		this.numOfIsoforms=1; this.isoforms=new ArrayList <Gene>();
	}
	/*public GeneWithIsoforms(String chr, int start, int end, String name, String orientation, int[] exonsStart, int[] exonsEnd){
		super(chr,start,end,name,orientation,exonsStart,exonsEnd);
		this.numOfIsoforms=1; this.isoforms=new ArrayList <Gene>();
	}*/
	@Deprecated  // start/end and exons are redundant
	public GeneWithIsoforms(String chr, int start, int end, String name, String orientation, List<Integer> exonsStart, List<Integer> exonsEnd){
		super(chr,start,end,name,orientation,exonsStart,exonsEnd);
		this.numOfIsoforms=1; this.isoforms=new ArrayList <Gene>();
	}
	@Deprecated  // start/end and exons are redundant
	public GeneWithIsoforms(String chr, int start, int end, String name, String orientation, List<Integer> exonsStart, List<Integer> exonsEnd, int blockStart, int blockEnd){
		super(chr,start,end,name,orientation,exonsStart,exonsEnd,blockStart,blockEnd);
		this.numOfIsoforms=1; this.isoforms=new ArrayList <Gene>();
	}
	@Deprecated  // start/end and exons are redundant
	public GeneWithIsoforms(String chr, int start, int end, String name, String orientation, List<Integer> exonsStart, List<Integer> exonsEnd, int blockStart, int blockEnd,String [] extraData){
		super(chr,start,end,name,orientation,exonsStart,exonsEnd,blockStart,blockEnd, extraData);
		this.numOfIsoforms=1; this.isoforms=new ArrayList <Gene>();
	}
	@Deprecated  // start/end and exons are redundant
	public GeneWithIsoforms(String chr, int start, int end, String name, String orientation, List<Integer> exonsStart, List<Integer> exonsEnd, String [] extraData){
		super(chr,start,end,name,orientation,exonsStart,exonsEnd, extraData);
		this.numOfIsoforms=1; this.isoforms=new ArrayList <Gene>();
	}
		
	
	public GeneWithIsoforms(Gene gene) {
		//super (gene.toBED(),false);
		super(gene.getChr(), gene.getName(), gene.getOrientation(), gene.getExonSet(), gene.getCDSStart(), gene.getCDSEnd());
		
		this.numOfIsoforms=1; 
		isoforms=new ArrayList <Gene>();
		initFromGene(gene);
	}
	
	public boolean IsExactOverlappingIsoform(Gene gene){
		return gene.getChr().equalsIgnoreCase(this.getChr()) && gene.getStart()== this.getStart() && gene.getEnd()==this.getEnd();
	}
	
	//Add an isoform
	public boolean AddIsoform (Gene gene) {
		return AddIsoform(gene, true);
	}
	
	public boolean AddIsoform (Gene gene, boolean updateMainIsoform){
		
		if (this.IsExactOverlappingIsoform(gene)){
			if (!updateMainIsoform || gene.getTranscriptLength() <= this.getTranscriptLength()){
				this.isoforms.add(gene);
				this.numOfIsoforms++;
			}
			else{
				Gene currGene= this.makeRefSeqGeneInstance();
				this.isoforms.add(currGene);
				this.numOfIsoforms++;
				//this.updateMainIsoform(gene);
				
			}
			return true;
		} else {
			System.err.println("Could not add isoform " + gene.toUCSC() + " to iso set with main iso: " + toUCSC());
		}
		return false;
	}
	
	
	public boolean addContainedIsoform (Gene isoform){
		boolean okToAdd = geneSpanContains(isoform) && overlaps(isoform);
		if(okToAdd) {

			isoforms.add(isoform);
			numOfIsoforms++;
		}
		
		return okToAdd;
	}
	
	public boolean AddAllIsoforms(GeneWithIsoforms overlapper) {
		return AddAllIsoforms(overlapper, true);
	}
	
	public boolean AddAllIsoforms(GeneWithIsoforms overlapper, boolean updateMainIsoform) {
		boolean res=true;
		Collection<Gene> set=overlapper.getAllIsoforms();
		if (set.isEmpty()) return false;
		for (Gene G : set){
			res=!res ? false : this.AddIsoform(G, updateMainIsoform);
			//if(!res) return false;
		}
		return res;
		
	}
	
	public boolean AddAllIsoNotIncludedforms(GeneWithIsoforms other) {
		boolean res=true;
		Collection<Gene> otherIsos=other.getAllIsoforms();
		Collection<Gene> thisIsos = getAllIsoforms();
		if (otherIsos.isEmpty()) return false;
		for (Gene G : otherIsos){
			boolean contained = false;
			Iterator<Gene> thisIsoIt = thisIsos.iterator();
			while(!contained && thisIsoIt.hasNext()) {
				Gene iso = thisIsoIt.next();
				contained = iso.compareTo(G) == 0;
			}
			if(!contained) {
				res=!res ? false : this.AddIsoform(G, false);
			}
		}
		return res;
	}
	
	public Collection<Gene> getAllIsoforms() {

		return getAllIsoforms(true);
	}
	
	/**
	 * In cases were the containg gene is artificial, like when building gene loci and adding
	 * all isoforms contain in the loci. The containing gene is not desired. This methods allows
	 * the caller to exclude the containing isoform from the list.
	 * @param includeContainingIsoform - boolean flag indicating whether or not to include the containing isoform
	 * @return Collection<RefSeqGene> of all the isoforms for this gene
	 */
	public Collection<Gene> getAllIsoforms(boolean includeContainingIsoform) {
		
		//RefSeqGene currGene= makeRefSeqGeneInstance();
		Collection <Gene> rtrn= new ArrayList<Gene>();
		if (! this.isoforms.isEmpty())
			rtrn.addAll(this.isoforms);
		if(includeContainingIsoform) {
			rtrn.add(this);
		}
		return rtrn;
	}
	
	//This function assigns the same bed score to all isoforms
	public void setBedScore(double scr){
		//throw new UnsupportedOperationException("Broken");
		super.setBedScore(scr);
		if (! this.isoforms.isEmpty()){
			for (Gene iso:this.isoforms)
				iso.setBedScore(scr);
		}
	}
	
	public void setExtraFields(double [] scores) {
		super.setExtraFields(scores);
		/**
		 * Commented out by skadri on 04/02/12
		 * Replaces extra fields for all children of a gene. 
		 */
		/*if (! this.isoforms.isEmpty()){
			for (RefSeqGene iso: isoforms)
				iso.setExtraFields(scores);
		}*/
	}
	
	//This function assigns the update the count score to equal the bed score in all isoforms
	public void updateScrToBedScore(double score){
		super.setBedScore(score);
		for (Gene Iso:this.isoforms)
			Iso.setBedScore(score);
	}
	
	public List<Annotation> getScoredExons(){
		List<Annotation> allExons= new ArrayList<Annotation>();
		for (Gene Iso:this.getAllIsoforms()){
			allExons.addAll(Iso.getExonSet());
		}
		return allExons;
	}
	//Get Merged Transcript
	
	
	/*private void updateMainIsoform(Gene gene) {

		chr=gene.getChr();
		start=gene.getStart();
		stop=gene.getEnd();
		name=gene.getName();
		orientation=gene.getOrientation();
		sequence=gene.getSequence();
		bedScore=gene.getBedScore(); //the transcript score as it appears in a bed file
		extraFields=gene.getExtraFields();
		setExons(gene.getExons(),gene.getExonsScores());
		
		
	}*/
	
	private Gene makeRefSeqGeneInstance() {
		
		//return new RefSeqGene(this.getChr(), this.getStart(), this.getEnd(), this.getName(),this.getCountScore(),this.getBedScore(), this.getOrientation(), this.exonStarts, this.exonEnds,this.exonScores,this.sequence);
		//RefSeqGene g=new RefSeqGene(toBED(),false);
		Gene g=new Gene(this);
		//throw new UnsupportedOperationException("Broken");
		g.setSequence(getSequence());
		g.setBedScore(getBedScore());
		g.setCountScore(getCountScore());
		return  g;
	}
	
	public boolean overlapsByGenomicRegion(
			IntervalTree<GeneWithIsoforms> otherTree, boolean considerOrientation) {
			
			Iterator<Node< GeneWithIsoforms>> geneIt=otherTree.overlappers(this.getStart(),this.getEnd());
			while (geneIt.hasNext()){
				Gene otherGene=geneIt.next().getValue();
				if (considerOrientation && this.overlapsGene(otherGene))
					return true;
				if (!considerOrientation && this.overlapsGeneInAnyOrientation(otherGene))
					return true;
			}
			return false;
	}

	public void dedup() {
		isoforms = new TreeSet<Gene>(isoforms);
	}

	public Gene findCompatibleGenes(IntervalTree<GeneWithIsoforms> overlapTree, int[] numIntrons) {
		
		Gene rtrn=null;
		int bestIntronNum=0;
		Iterator <GeneWithIsoforms> gIt=overlapTree.valueIterator();
		while(gIt.hasNext()){
			Collection<Gene> genes=gIt.next().getAllIsoforms();
			for (Gene g: genes){
				Collection<Gene> myIso=this.getAllIsoforms();
				for (Gene iso:myIso){
					int intronNum=iso.numOfCompatibleIntrons(g);
					if (intronNum > bestIntronNum){
						bestIntronNum=intronNum;
						rtrn=g;
					}
				}
			}
		}
		
		numIntrons[0]=bestIntronNum;
		return rtrn;
	}
	public void expandUtrs(Integer utr1, Integer utr2) {
		super.expandUtrs(utr1,utr2);
		for(Gene g: this.isoforms){
			g.expandUtrs(utr1,utr2);
		}
		
	}
	
	/*public void set3PrimeEnd(int updated3Prime) {
		super.set3PrimeEnd(updated3Prime);
		for(Gene g: this.isoforms){
			g.set3PrimeEnd(updated3Prime);
		}
	}*/
	
	/*public void set5PrimeEnd(int updated5Prime) {
		super.set5PrimeEnd(updated5Prime);
		for(Gene g: this.isoforms){
			g.set5PrimeEnd(updated5Prime);
		}
	}*/
	
	
	//adds a suffix to all the isoforms
	public void addSuffixToName(String refName) {

		this.setName(this.getName()+refName);
		for (Gene g:this.isoforms){
			g.setName(g.getName()+refName);
		}
		
	}
	
	
	
	public boolean overlapsExon(LightweightGenomicAnnotation exon){
		boolean overlaps = super.overlapsExon(exon);
		
		Iterator<Gene> isoformIt = isoforms.iterator();
		while(!overlaps && isoformIt.hasNext()) {
			overlaps = isoformIt.next().overlapsExon(exon);
		}
		
		return overlaps;
	}
	
	
	
	public boolean overlapsExon(Gene gene) {
		boolean overlaps = super.overlaps(gene);
		
		Iterator<Gene> isoformIt = isoforms.iterator();
		while(!overlaps && isoformIt.hasNext()) {
			overlaps = isoformIt.next().overlapsExon(gene);
		}
		
		return overlaps;
	}
	
	
	public boolean overlaps(Gene other) {
		boolean overlaps = super.overlaps(other);
		
		Iterator<Gene> isoformIt = isoforms.iterator();
		while(!overlaps && isoformIt.hasNext()) {
			Gene isoform = isoformIt.next(); 
			if(!isoform.equals(this)) {
				overlaps = isoform.overlaps(other);
			}
		}
		
		return overlaps;
	}
	
	public Gene getMerged() {

		Gene curr=this;
		for (Gene iso: this.isoforms)
				curr=curr.takeUnion(iso);
		return curr;
	}

	
public void cleanIsoforms() {
		isoforms = new ArrayList<Gene>();
		numOfIsoforms = 1;
	}
	
	
	
	public Collection<Annotation> getAllIsoformExonSet(){
		Collection<Annotation> rtrn=new TreeSet<Annotation>(); 
		
		for(Gene isoform : isoforms) {
			rtrn.addAll(isoform.getExonSet());
		}

		return rtrn;
	}
	
	public Collection<Annotation> getAllIsoformIntronSet(){
		Collection<Annotation> rtrn=new TreeSet<Annotation>(); 
		
		for(Gene isoform : isoforms) {
			rtrn.addAll(isoform.getIntronSet());
		}

		return rtrn;
	}
	
	public Gene constituentIsoform() {
		Collection<? extends Annotation> exons = getAllIsoformExonSet();
		exons = CollapseByIntersection.collapseByIntersection(exons, true);
		Set<Alignments> constituentExons = new TreeSet<Alignments>();
		Collection<Gene> isoforms = getAllIsoforms();
		for(Annotation exon : exons) {
			Iterator<Gene> isoformIt = isoforms.iterator();
			boolean overlaps = true;
			while(overlaps &&  isoformIt.hasNext()) {
				Gene isoform = isoformIt.next();
				overlaps = isoform.overlapsExon(exon);
				//System.err.println("exon " + exon.toUCSC() +" overlaps " + isoform.toUCSC() + "("+isoform.getName()+")? " + overlaps);
			}
			if(overlaps) {
				//System.err.println("Adding exon " + exon.toUCSC());
				constituentExons.add(new Alignments(exon));
			}
		}
		
		
		Gene constituentIsoform = null;
		
		if(constituentExons.size() >0) {
			constituentIsoform = new Gene(constituentExons);
			constituentIsoform.setName(getName());
			constituentIsoform.setBedScore(getBedScore());
			constituentIsoform.setOrientation(getOrientation());
		}
		return constituentIsoform;
	}
	
	/**
	 * Constructs an artificial gene where the constituent introns (regions that are intronic in all isoforms) are the exons of the gene.
	 * @return
	 */
	public Gene constituentIntrons() {
		Collection<? extends Annotation> introns = getAllIsoformIntronSet();
		introns = CollapseByIntersection.collapseByIntersection(introns, true);
		Set<Alignments> constituentIntrons = new TreeSet<Alignments>();
		Collection<Gene> isoforms = getAllIsoforms();
		for(Annotation intron : introns) {
			Iterator<Gene> isoformIt = isoforms.iterator();
			boolean isContained = true;
			while(isContained &&  isoformIt.hasNext()) {
				Gene isoform = isoformIt.next();
				Iterator<? extends Annotation> isoIntronIt = isoform.getIntronSet().iterator();
				boolean containedInIsoformIntron = false;
				while(!containedInIsoformIntron && isoIntronIt.hasNext()) {
					Annotation isoIntron = isoIntronIt.next();
					containedInIsoformIntron = isoIntron.contains(intron);
				}
				isContained = containedInIsoformIntron;
				//System.err.println("exon " + exon.toUCSC() +" overlaps " + isoform.toUCSC() + "("+isoform.getName()+")? " + overlaps);
			}
			if(isContained) {
				constituentIntrons.add(new Alignments(intron));
			}
		}
		
		Gene constituentIntroform = new Gene(constituentIntrons);
		constituentIntroform.setName(getName());
		constituentIntroform.setBedScore(getBedScore());
		constituentIntroform.setOrientation(getOrientation());
		
		return constituentIntroform;
	}
	
	
	public void colapse(int exonFudgeFactor) {
		List<Gene> collapsedIsoforms = new ArrayList<Gene>();
		Iterator<Gene> isoIt = isoforms.iterator(); 
		while ( isoIt.hasNext()) {
			Gene iso = isoIt.next();
			boolean foundAlmostEqual = this.almostEqual(iso, exonFudgeFactor);
			Iterator<Gene> collapsedIsoIt = collapsedIsoforms.iterator();
			while(!foundAlmostEqual && collapsedIsoIt.hasNext()) {
				foundAlmostEqual = collapsedIsoIt.next().almostEqual(iso, exonFudgeFactor);
			}
			if(!foundAlmostEqual) {
				collapsedIsoforms.add(iso);
			}
		}
		
		isoforms = collapsedIsoforms;
		/*if(numOfIsoforms != isoforms.size()+1) {
			System.err.println("Collapsed isoforms for " + getName() + " from  " + numOfIsoforms + " to " + (isoforms.size()+1));
		}*/
		numOfIsoforms = collapsedIsoforms.size() + 1; //TODO: why do we have this field? it is completely redundant.
		
	}
	
	public boolean almostContains(Gene other, int exonFudgeFactor) {
		boolean isAlmostContained = false;
		//System.err.println("Checking whether refseqGene " + other.toBED() + " is contained in refseqgenewithisos " + toBED() );
		Iterator<Gene> isoIt = getAllIsoforms().iterator();
		while(!isAlmostContained && isoIt.hasNext()) {
			Gene iso = isoIt.next();
			throw new UnsupportedOperationException("Broken");
			//isAlmostContained = iso.almostContainsRefSeq(other, exonFudgeFactor);
		}
		
		return isAlmostContained;
	
	}
	
	public boolean almostContains(GeneWithIsoforms other, int exonFudgeFactor) {
		boolean isAlmostContained = false;
		Iterator<Gene> otherIsoIt = other.getAllIsoforms().iterator();
		while(!isAlmostContained && otherIsoIt.hasNext()) {
			Gene otherIso = otherIsoIt.next();
			isAlmostContained = almostContains(otherIso, exonFudgeFactor);
		}
		
		return isAlmostContained;
	
	}
	
	
	
	
	public boolean hasOverlapping3PUTR(GeneWithIsoforms other) {
		boolean overlaps = false;
		if("*".equals(getOrientation()) || "*".equals(other.getOrientation())) {
			overlaps = false;
		} else {
			overlaps =  get3PrimeExon().overlaps(other.get3PrimeExon());
		}
		return overlaps;
	}
	
	public boolean hasOverlapping5PUTR(GeneWithIsoforms other) {
		return get5PrimeExon().overlaps(other.get5PrimeExon());
	}

	public String allIsoformsToBED() {
		StringBuilder sb = new StringBuilder();
		for(Gene iso : getAllIsoforms()) {
			sb.append(iso.toBED());
			sb.append("\n");
		}
		return sb.toString();
	}
	public void standardizeOreintation() {
		int plus = 0;
		int minus = 0;
		for(Gene iso : getAllIsoforms()) {
			if ("+".equals(iso.getOrientation())) {
				plus++;
			} else if ("-".equals(iso.getOrientation())) {
				minus++;
			} 
		}
		
		if(plus > minus) {
			for(Gene iso : getAllIsoforms()) {
				iso.setOrientation("+");
			}
		} else 	if(minus > plus) {
			for(Gene iso : getAllIsoforms()) {
				iso.setOrientation("-");
			}
		}else{
			for(Gene iso : getAllIsoforms()) {
				iso.setOrientation(Strand.UNKNOWN);
			}
		}
		
	}
	
 	public String toString(String name){
		String rtrn="";
		Collection<Gene> genes=this.getAllIsoforms();
		for(Gene gene: genes){
			gene.setName(name);
			rtrn+=gene.toBED()+"\n";
		}
		return rtrn;
 	}
	public Integer getNumberOfIsoforms() {
		
		return this.isoforms.size();
	}
	
	
	public Collection<? extends Gene> selectRandIsoSubset(
			Integer maxIsoPerLoci) {
		
		Collection<Gene> rtrn = new TreeSet <Gene>();
		int i=0;
		for (Gene g: this.getAllIsoforms()){
			if (i< maxIsoPerLoci){
				rtrn.add(g); i++;
			}
			else
				break;
		}
			
		return rtrn;
	}


}
