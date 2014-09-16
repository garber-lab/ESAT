package umms.core.annotation;

import broad.core.parser.StringParser;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import umms.core.feature.GeneWindow;

import org.apache.log4j.Logger;

import broad.core.annotation.BasicGenomicAnnotation;
import broad.core.annotation.BasicLightweightAnnotation;
import broad.core.annotation.GFF;
import broad.core.annotation.LightweightGenomicAnnotation;
import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.core.sequence.Sequence;
import broad.pda.datastructures.Alignments;
import broad.pda.rnai.ExtractSequence;

public class Gene extends BasicAnnotation {
	static Logger logger = Logger.getLogger(Gene.class.getName());
	ArrayList<Double> scores;
	Map<String, String> attributes;
	int cdsStart; // beginning of CDS
	int cdsEnd; // end of CDS
	double[] exonScores;
	private String sequence;
	private String samRecord;
	private double countScore=0; //Moran -added this as init value 
	double bedScore; //the transcript score as it appears in a bed file
	String [] extraFields;
	private Collection<Gene> isoforms;
	
	public static final Pattern START_CODON_PATRN= Pattern.compile("ATG",  Pattern.CASE_INSENSITIVE);
	public static final String [] STOP_CODONS =  {"TAG", "TAA", "TGA"};
	
	
	public Gene(String chr, String name, Strand orientation, Collection<? extends Annotation> exons, int cdsStart, int cdsEnd){
		super(chr, orientation, name, exons);
		init(cdsStart, cdsEnd);
	}
	
	public Gene (String chr, int start, int end, String name, Strand orientation, int cdsStart, int cdsEnd) {
		super(chr, start, end, orientation, name);
		init(cdsStart, cdsEnd);
	}
	
	public Gene(Collection<? extends Annotation> exons, String name, Strand orientation, int cdsStart, int cdsEnd){
		super(exons, orientation, name);
		init(cdsStart, cdsEnd);
	}
	

	private void init(int cdsStart, int cdsEnd) {
		this.exonScores = new double[numBlocks()];
		this.cdsStart=cdsStart;
		this.cdsEnd=cdsEnd;
	}
	

	public Gene(String chr, String name, Strand orientation, Collection<? extends Annotation> exons){
		this(chr, name, orientation, exons, 0,0);
	}
	
	public Gene(Annotation annotation, String name, Strand orientation, Collection<? extends Annotation> exons){
		this(annotation.getChr(), name, orientation, exons, 0,0);
	}
	
	public Gene(Annotation annotation, Collection<? extends Annotation> exons){
		this(annotation.getChr(), annotation.getName(), annotation.getOrientation(), exons, 0,0);
	}
	
	public Gene(Annotation annotation, Collection<? extends Annotation> exons, String name, int cdsStart, int cdsEnd){
		this(annotation.getChr(), name, annotation.getOrientation(), exons, cdsStart, cdsEnd);
	}
	
	public Gene(Annotation annotation, Collection<? extends Annotation> exons, String name){
		this(annotation, exons, name, annotation.getStart(), annotation.getStart());
	}
		
	public Gene(Collection<? extends Annotation> exons, String name, Strand orientation){
		this(exons, name, orientation, 0, 0);
	}
	
	public Gene (String chr, int start, int end, String name, Strand orientation) {
		this(chr, start, end, name, orientation, 0, 0);
	}
	
	public Gene (String chr, int start, int end, String name){
		this(chr, start, end, name, Strand.UNKNOWN);
	}
	
	public Gene(Collection<? extends Annotation> exons, String name) {
		this(exons, name, exons.iterator().next().getOrientation());
	}
	
	public Gene(Collection<? extends Annotation> exons) {
		this(exons, "");
	}
	
	public Gene (String chr, int start, int end){
		this(chr, start, end, "");
	}
	
	@Deprecated
	public Gene (Annotation annotation) {
		this(annotation.getChr(), annotation.getName(), annotation.getOrientation(), annotation.getBlocks());
	}
	
	public Gene(String pslString){
		this(pslString, true);
	}
	
	public Gene(String pslString, boolean isPSLFormat){
		this(makeGene(pslString, isPSLFormat));
	}
	
	public Gene(String chr2, int start2, int end, String name2,	double bedScore2, String orientation2, String sequence2,String[] extraFields2) {
		this(chr2,start2,end,name2,bedScore2,orientation2,null,null,null,sequence2,extraFields2);
	}
	
	public Gene(String chr2, int start2, int end, String name2,	double bedScore2, String orientation2, String sequence2) {
		this(chr2,start2,end,name2,bedScore2,orientation2,null,null,null,sequence2);
	}
	
	public Gene(String chr, int start, int end, String name,double mybedScore, String orientation, int[] exonsStart, int[] exonsEnd){
		this(chr, name, AbstractAnnotation.getStrand(orientation), makeExons(chr, exonsStart, exonsEnd));
		this.bedScore=mybedScore;	
		//logger.error("Caution: Gene constructor specifies both start/end and exons; exons will override start/end.  If you can remove constructor, do it. -JE");
	}

	public Gene(String chr, int start, int end, String name, String orientation, List<Integer> exonsStart, List<Integer> exonsEnd){
		this(chr, start, end, name, orientation, exonsStart, exonsEnd, null);
		//logger.error("Caution: Gene constructor specifies both start/end and exons; exons will override start/end.  If you can remove constructor, do it. -JE");
	}
	
	public Gene(String chr, int start, int end, String name, String orientation, List<Integer> exonsStart, List<Integer> exonsEnd, int blockStart, int blockEnd){
		this(chr, name, AbstractAnnotation.getStrand(orientation), makeExons(chr, exonsStart, exonsEnd));
		this.cdsStart=blockStart;
		this.cdsEnd=blockEnd;
	}
	
	//Bug fix Dec 28th 2010 , was not intializing extraFields
	public Gene(String chr, int start, int end, String name, String orientation, List<Integer> exonsStart, List<Integer> exonsEnd, String [] extraData){
		this(chr, name, AbstractAnnotation.getStrand(orientation), makeExons(chr, exonsStart, exonsEnd));
		//logger.error("Caution: Gene constructor specifies both start/end and exons; exons will override start/end.  If you can remove constructor, do it. -JE");
		if(extraData != null) 
			setExtraFields(extraData);
	}
	
	public Gene(String chr2, int start2, int end, String name2,	double bedScore2, String orientation2, int[] exonStarts2,int[] exonEnds2, double[] exonScores2,String sequence2) {
		this(chr2, name2, AbstractAnnotation.getStrand(orientation2), makeExons(chr2, exonStarts2, exonEnds2));	
		this.exonScores=exonScores2;
		this.bedScore=bedScore2;
		this.sequence=sequence2;
		//logger.error("Caution: Gene constructor specifies both start/end and exons; exons will override start/end.  If you can remove constructor, do it. -JE");
	}
	
	public Gene(String chr2, int start2, int end, String name2,	double bedScore2, String orientation2, int[] exonStarts2,int[] exonEnds2, double[] exonScores2,String sequence2,String[] extraFields2) {
		this(chr2, start2, end, name2, bedScore2, orientation2, exonStarts2, exonEnds2, exonScores2, sequence2);	
		this.extraFields=extraFields2;
	}
	
	public Gene(String chr, int start, int end, String name,double score, String strand, List<Integer> exonsStart,	List<Integer> exonsEnd, int blockStart2, int blockEnd2) {
		this( chr,  start,  end,  name, strand,  exonsStart,  exonsEnd,  blockStart2, blockEnd2);
		this.bedScore=score;
		//logger.error("Caution: Gene constructor specifies both start/end and exons; exons will override start/end.  If you can remove constructor, do it. -JE");
	}
	
	public Gene(String chr, int start, int end, String name, String orientation, List<Integer> exonsStart, List<Integer> exonsEnd, int blockStart, int blockEnd, String [] extraData){
		this(chr, name, AbstractAnnotation.getStrand(orientation), makeExons(chr, exonsStart, exonsEnd), blockStart, blockEnd);
		if(extraData != null) {
			this.extraFields = extraData;
		}
		//logger.error("Caution: Gene constructor specifies both start/end and exons; exons will override start/end.  If you can remove constructor, do it. -JE");
	}
	
	public Gene(String chr, int start, int end, String name,double score, String strand, List<Integer> exonsStart,List<Integer> exonsEnd, int blockStart2, int blockEnd2, String [] extraColumns) {
		this( chr,  start,  end,  name, strand,  exonsStart,  exonsEnd,  blockStart2, blockEnd2, extraColumns);
		this.bedScore=score;
		//logger.error("Caution: Gene constructor specifies both start/end and exons; exons will override start/end.  If you can remove constructor, do it. -JE");
	}
	
	
	
	public Gene(Gene gene) {
		this(gene.getChr(), gene.getName(), gene.getOrientation(), gene.getExonSet(), gene.getCDSStart(), gene.getCDSEnd());
		this.bedScore=gene.getScore();
		initFromGene(gene);
			
	}

	protected void initFromGene(Gene gene) {
		if (gene.extraFields != null)
			setExtraFields(gene.getExtraFields());
		if (gene.scores !=null)
			setScores(gene.getScores());		
		if (gene.attributes !=null)
			setAttributes(gene.getAttributes());
		if (gene.samRecord!=null)
			this.samRecord=gene.getSAMString();
		if (gene.countScore!=0)  
			this.countScore=gene.getCountScore();
		if(gene.sequence!=null)
			setSequence(gene.sequence);
		this.setBedScore(gene.getBedScore());
	}
	
	

	
	private static Collection<? extends Annotation> makeExons(String chr, int[] exonsStart, int[] exonsEnd) {
		Collection<Annotation> rtrn=new ArrayList<Annotation>();
		
		for(int i=0; i<exonsStart.length; i++){
			Annotation exon=new BasicAnnotation(chr, exonsStart[i], exonsEnd[i]);
			rtrn.add(exon);
		}
		
		return rtrn;
	}
	
	private static Collection<? extends Annotation> makeExons(String chr, List<Integer> exonsStart, List<Integer> exonsEnd) {
		Collection<Annotation> rtrn=new ArrayList<Annotation>();
		
		for(int i=0; i<exonsStart.size(); i++){
			Alignments align=new Alignments(chr, exonsStart.get(i), exonsEnd.get(i));
			rtrn.add(align);
		}
		
		return rtrn;
	}
	
	public void setSequence(String seq){this.sequence=seq;}
	
	public void setSequenceFromChromosome(Sequence chrSequence) {
		Sequence  geneSequence = new Sequence(getName(),this.length());
		Set<? extends Annotation> exons = getExonSet();
		for(Annotation exon : exons) {
			Sequence exonSeq = chrSequence.getSubSequence(getName(), exon.getStart(), exon.getEnd());
			geneSequence.append(exonSeq.getSequenceBases());
		}
		
		if(this.isNegativeStrand()) {
			geneSequence.reverse();
		}
		this.sequence = geneSequence.getSequenceBases();
	}
	
	public void setOrientation(String orientation){
		setOrientation(AbstractAnnotation.getStrand(orientation));
	}
	
	public void setSAMString(String sam){this.samRecord=sam;}
	
	public void setExtraFields(double[] d) {
		extraFields = new String[d.length];
		for (int i = 0; i < d.length; i++) {
			String s = (new Double(d[i])).toString();
			extraFields[i] = s;
		}
	}
	
	public void setExtraFields(double val, int i) {
		if (i>0 && i<=this.extraFields.length)
			this.extraFields[i-1]=String.valueOf(val);
		
	}
	
	public void setExtraFields(String[] extraData){
		String[] extraDataCpy= new String[extraData.length];
		for (int i=0; i<extraData.length; i++)
			extraDataCpy[i]=extraData[i];
		this.extraFields = extraDataCpy;
	}
	
	public int getCDSEnd() {
		return this.cdsEnd;
	}

	public int getCDSStart() {
		return this.cdsStart;
	}
	
	private void setScores(ArrayList<Double> scores2) {
		ArrayList<Double> newscores = new ArrayList<Double>();
		for (int i=0; i< scores2.size(); i++)
			newscores.add(new Double(scores2.get(i)));
		this.scores=newscores;
		
	}
	
	private static List<Integer> [] setBlockStartsAndEnds(String[] blockStarts, String[] blockSizes, int size, int start){
		List<Integer> starts=new ArrayList<Integer> ();
		List<Integer>  end=new ArrayList<Integer> ();
		for(int i=0; i<size; i++){
			starts.add(start+new Integer(blockStarts[i].replaceAll("\"", "").trim()));
			end.add((Integer)starts.get(i)+new Integer(blockSizes[i].replaceAll("\"", "").trim()));
		}
		List [] rtrn={starts, end};
		return rtrn;
	}
	
	private void setAttributes(Map<String, String> attrs) {
		Map<String, String> m= new HashMap<String, String>();
		for (String key: attrs.keySet())
			m.put(key, attrs.get(key));
		this.attributes=m;
	}
	
	/*public void set5PrimeEnd(int updated5Prime) {
		if("-".equals(orientation)) {
			for(int i = 0 ; i < exonEnds.length; i++){
				if(exonEnds[i] ==  stop) {
					exonEnds[i] = updated5Prime;
				}
			}
			this.stop = updated5Prime;
		} else {
			for(int i = 0 ; i < exonEnds.length; i++){
				if(exonStarts[i] ==  start) {
					exonStarts[i] = updated5Prime;
				}
			}
			this.start = updated5Prime;
			
		}
	}
	
	public void set3PrimeEnd(int updated3Prime) {
		if("-".equals(orientation)) {
			for(int i = 0 ; i < exonEnds.length; i++){
				if(exonStarts[i] ==  start) {
					exonStarts[i] = updated3Prime;
				}
			}
			this.start = updated3Prime;
		} else {
			for(int i = 0 ; i < exonEnds.length; i++){
				if(exonEnds[i] ==  stop) {
					exonEnds[i] = updated3Prime;
				}
			}
			this.stop = updated3Prime;
			
		}
	}
	
	public void set5PrimeEndForGeneOnly(int updated5Prime) {
		if("-".equals(orientation)) {
			for(int i = 0 ; i < exonEnds.length; i++){
				if(exonEnds[i] ==  stop) {
					exonEnds[i] = updated5Prime;
				}
			}
			this.stop = updated5Prime;
		} else {
			for(int i = 0 ; i < exonEnds.length; i++){
				if(exonStarts[i] ==  start) {
					exonStarts[i] = updated5Prime;
				}
			}
			this.start = updated5Prime;
			
		}
	}
	
	public void set3PrimeEndForGeneOnly(int updated3Prime) {
		if("-".equals(orientation)) {
			for(int i = 0 ; i < exonEnds.length; i++){
				if(exonStarts[i] ==  start) {
					exonStarts[i] = updated3Prime;
				}
			}
			this.start = updated3Prime;
		} else {
			for(int i = 0 ; i < exonEnds.length; i++){
				if(exonEnds[i] ==  stop) {
					exonEnds[i] = updated3Prime;
				}
			}
			this.stop = updated3Prime;
			
		}
	}*/
	
	/**
	 * Set coding region coordinates (includes start codon but not stop codon)
	 * @param start start position
	 * @param end end position
	 */
	public void setCDSRegion(int start, int end){
		if(start > end) {
			throw new IllegalArgumentException("Start position " + start + " is greater than end position " + end);
		}
		if((start < getStart() && start < getEnd()) || (start > getStart() && start > getEnd())) {
			throw new IllegalArgumentException("Start position " + start + " is out of range of gene start and end " + getStart() + "-" + getEnd());
		}
		if((end < getStart() && end < getEnd()) || (end > getStart() && end > getEnd())) {
			throw new IllegalArgumentException("End position " + end + " is out of range of gene start and end " + getStart() + "-" + getEnd());
		}
		this.cdsStart = start;
		this.cdsEnd = end;
	}
	
	//TODO I dont know how to deal with this
	/*protected void setExons(Alignments[] exons,	double[] giveExonsScores) {
		
		this.exonStarts=new int[exons.length] ;
		this.exonEnds=new int[exons.length];
		this.exonScores=new double[exons.length];
		for (int i=0; i<exons.length; i++) {
			this.exonStarts[i]=exons[i].getStart();
			this.exonEnds[i]=exons[i].getEnd();
			this.exonScores[i]= giveExonsScores != null ? giveExonsScores [i] : 0;
		}
		
	}*/
	
	
	public void setBedScore (double d){bedScore=d;}
	
	public void setCountScore(double score){this.countScore=score;}
	
	private String setN(int len){
		String rtrn="";
		for(int i=0; i<len; i++){rtrn+="N";}
		return rtrn;
	}
	
	public int getSize(){
		return length();
	}
	
	public int getGenomicLength() { return getLengthOnReference(); }
	
	public boolean isNegativeStrand(){return getOrientation().equals(Strand.NEGATIVE);}
	
	public String getSequence(){return this.sequence;}
	
	public Sequence getSequence(Sequence chrSequence) {
		Sequence seq=new Sequence(getName());
		Collection<? extends Annotation> exons=this.getSortedAndUniqueExons();
		
		for(Annotation exon: exons){
			Sequence sequence=chrSequence.getSubSequence("", exon.getStart(), exon.getEnd());
			seq.append(sequence.getSequenceBases());
		}
		
		return seq;
	}
	
	public String getSequence(Sequence chr2, boolean repeatMask, boolean stranded) throws Exception {
		return ExtractSequence.getSequenceForGene(this, chr2, repeatMask, new TreeMap(), stranded);
	}
	
	public Sequence getSequenceObject(){
		Sequence seq=new Sequence(getName());
		seq.setSequenceBases(this.sequence);
		return seq;
	}
	
	public String getString(){
		return getChr()+":"+getStart()+"-"+getEnd();
	}
	
	public String getSAMString(){return this.samRecord;}
	
	/**
	 * 
	 * @param name
	 * @return attribute value or 0 if the attribute is not defined.
	 */
	public String getAttribute(String name) {
		return attributes == null || !attributes.containsKey(name) ? null : attributes.get(name);
	}
	
	public Map getAttributes() { return attributes;} //TODO: Return a copy, attributes should not be modified outside of accessor methods.
	
	public String[] getExtraFields() {
		if (this.extraFields!=null)
			return extraFields;
		else 
			return null;
	}
		
	public String getExtraFields(int i) {
		if (this.extraFields!=null &&  this.extraFields.length>=(i+1))
			return this.extraFields[i];
		else
		{	
			//System.err.println("no extra fields");
			return "";
		}
	}
	
	private ArrayList<Double> getScores() {
		return this.scores;
	}
	
	public double getNormalizedScore(){
		double score=this.getScore();
		Set<? extends Annotation> exons=this.getExonSet();
		
		if(exons.isEmpty()) return 0;
		
		return score/length();
	}
	
	public double getRPKM() {
		if (extraFields != null)
			return Double.parseDouble(extraFields[4]);
		else
			return -1;
	}
	
	public double getPValue() {
		if (extraFields != null)
			return Double.parseDouble(extraFields[0]);
		else
			return -1;
	}
	
	public double getCountScore(){return this.countScore;}
	
	public double[] getExonsScores(){
		return this.exonScores;
	}
	
	@Deprecated // JE - length() returns the blocked length.  getLengthOnReference() returns getEnd() - getStart()
	public int getTranscriptLength(){
		return length();
	}
	
	/**
	 * Get the set of exons comprising the 5' UTR
	 * @return the set of exons comprising the 5' UTR as a collection of Alignments objects
	 */
	public Collection<? extends Annotation> get5Utr() {
		return this.intersect(get5UtrIncludingIntrons()).getBlocks();
	}
	
	/**
	 * Get the set of exons comprising the 3' UTR
	 * @return the set of exons comprising the 3' UTR as a collection of Alignments objects
	 */
	public Collection<? extends Annotation> get3Utr() {
		return this.intersect(get3UtrIncludingIntrons()).getBlocks();
	}

	// TODO JE - I'm confused why this doesn't return a blocked UTR 
	public Gene get3UTRGene() {
		if(!hasCDS()){return new Gene(this);}
		Annotation UTRRegion=get3UTR();
		Gene rtrn=this.trimAbsolute(UTRRegion.getStart(), UTRRegion.getEnd());	
		if(rtrn == null) {
			return rtrn;
		}
		String geneName = "";
		if(getName() == null) {
			geneName += getChr() + "_" + getStart() + "_" + getEnd();
		} else {
			geneName += getName();
		}
		rtrn.setName(geneName + "_3UTR");
		return rtrn;
	}

	// TODO JE - I'm confused why this doesn't return a blocked UTR 
	public Gene get5UTRGene() {
		if(!hasCDS()){return new Gene(this);}
		Annotation UTRRegion=get5UTR();
		Gene rtrn=this.trimAbsolute(UTRRegion.getStart(), UTRRegion.getEnd());		
		if(rtrn == null) {
			return rtrn;
		}
		String geneName = "";
		if(getName() == null) {
			geneName += getChr() + "_" + getStart() + "_" + getEnd();
		} else {
			geneName += getName();
		}
		rtrn.setName(geneName + "_5UTR");
		return rtrn;
	}
	
	// TODO JE - I'm confused why this doesn't return a blocked UTR 
	public Annotation get5UTR(){
		if(getOrientation() == Strand.POSITIVE){
			Annotation rtrn=new BasicAnnotation(getChr(), getStart(), this.getCDSRegion().getStart(), getStrand());
			return rtrn;
		}
		Annotation rtrn=new BasicAnnotation(getChr(), this.getCDSRegion().getEnd(), getEnd(), getStrand());
		return rtrn;
	}
	
	// TODO JE - I'm confused why this doesn't return a blocked UTR 
	public Annotation get3UTR(){
		if(getOrientation() == Strand.NEGATIVE){
			Annotation rtrn=new BasicAnnotation(getChr(), getStart(), this.getCDSRegion().getStart(), getStrand());
			return rtrn;
		}
		Annotation rtrn=new BasicAnnotation(getChr(), this.getCDSRegion().getEnd(), getEnd(), getStrand());
		return rtrn;
	}
	
	public boolean hasCDS(){
		if(this.cdsStart==this.cdsEnd){return false;}
		return true;
	}
	
	/**
	 * Get an iterator over all reference positions making up the gene
	 * @return Iterator over all positions in the gene
	 */
	public Iterator<Integer> getAllPositions() {
		TreeSet<Integer> positions = new TreeSet<Integer>();
		Collection<GeneWindow> windows = getWindows(1);
		for(GeneWindow window : windows) {
			positions.add(window.getStart());
		}
		return positions.iterator();
	}

	public Collection<GeneWindow> getWindows(int windowSize) {
		return getWindows(windowSize, 1, 0);
	}
	
	public Collection<GeneWindow> getWindows(int windowSize, int stepSize, int start) {
		if(stepSize < 1) {
			throw new IllegalArgumentException("Step size must be >= 1");
		}
		Collection<GeneWindow> subGenes=new TreeSet<GeneWindow>();
		if (windowSize > length()){
			GeneWindow window = new GeneWindow(this);
			window.addSourceAnnotation(this);
			subGenes.add(window);
		} else {
			for(int i=start; i< (length()+1)-windowSize; i=i+stepSize){
				GeneWindow subGene=trimGene(i, i+windowSize);
				if(subGene!=null){
					subGenes.add(subGene);
				}
		}
		}
		return subGenes;
	}

	/*public Collection<GeneWindow> getWindows(int windowSize, int stepSize, int start) {
		// Account for erroneous input by user
		if(stepSize < 1) {
			throw new IllegalArgumentException("Step size must be >= 1");
		}
		
		// Initializes a collection of GeneWindow (or Genes) called subGenes which is a TreeSet
		Collection<GeneWindow> subGenes = new TreeSet<GeneWindow>();
		
		// Our collection exons which is a generic implementation of Annotation is a collection of all exons
		
		Collection<? extends Annotation > exons = getBlocks();
		
		// Iterates through all exons in our collection exons
		for(Annotation exon: exons) {
			//logger.info("Exon: " + exon.toUCSC());
			// If the exon is larger than the input window size
			if(exon.getSize() >= windowSize) { //TODO should be >= ?
				logger.debug(this.toUCSC() + "\texon_size\t" + exon.getSize() + "\twindow_size\t" + windowSize);
				for(int i = exon.getStart(); i <= exon.getEnd() - windowSize; i++) {
					logger.debug(this.toUCSC() + "\ti\t" + i);
					GeneWindow subGene = new GeneWindow(new Gene(getChr(), i, i + windowSize));
					subGene.addSourceAnnotation(this);
					subGene.setOrientation(this.getOrientation());
					// GeneWindow subGene = new Basic(getChr(), i, i + windowSize);
					subGenes.add(subGene);
					logger.debug(this.toUCSC() + "\tadded\t" + subGene.toUCSC());
				}
			}
			
			for(int i = Math.max(exon.getEnd() - windowSize, exon.getStart()); i < exon.getEnd(); i++) {
				//logger.info("Size: " + this.getSize() + " ExonSize: " + exon.getSize() + " " + i + " start: " + exon.getStart() + " end: " + exon.getEnd());
				int relativeStart = getPositionAtReferenceCoordinate(i, true);
				int size = relativeStart + windowSize;
				logger.debug(this.toUCSC() + "\ti\t" + i + "\trel_start\t" + relativeStart + "\tsize\t" + size);
				if(relativeStart + windowSize < this.getSize()) { //TODO should be <= ?
					//Gene exonGene =  new Gene(exon);
					//GeneWindow subGene = this.trimAbsolute(i, i+windowSize);
					GeneWindow subGene = this.trimGene(relativeStart, relativeStart + windowSize);
					//GeneWindow subGene = this.trimGeneNew(exonGene, relativeStart, relativeStart+windowSize);
					subGene.setOrientation(this.getOrientation());
					subGenes.add(subGene);
					//logger.info("" + subGene.toBED());
					logger.debug(this.toUCSC() + "\tsubgene\t" + subGene.toUCSC());
				}
			}	
		}*/
		
		// Need full gene sequence, not just exons
		/*for(Annotation exon: exons) {
			// do intron windows
			for(int i = exon.getEnd() - windowSize + 1; i < exon.getEnd() - windowSize; i++) {
				 subGene = new Basic(getChr(), i, i + windowSize);
				subGenes.add((GeneWindow) subGene);
			}
		}
			
		return subGenes;
	}*/
	
	public Gene getStartCodon() {
		if(this.getCDS() == null) return null;
		Gene gene=this.getCDS();
		//System.err.println(gene);
		if(getOrientation() == null) return null;
		if(gene.getOrientation() == Strand.POSITIVE){
			return gene.trimGene(0, 3);
		}
		else if(gene.getOrientation() == Strand.NEGATIVE){
			return gene.trimGene(gene.length()-3, gene.length());
		}
		else{
			System.err.println("UNKNOWN Orientation... assuming +");
			return gene.trimAbsolute(gene.getStart(), gene.getStart()+3);
		}
	}
	

	
	/**
	 * Get the region of the gene span that lies 5' of start codon, including introns and the start codon itself
	 * @return the region of the gene span that lies 5' of start codon as an Alignments object
	 */
	public Annotation get5UtrIncludingIntrons(){
		if(getOrientation() == Strand.NEGATIVE){
			Annotation rtrn=new BasicAnnotation(getChr(), this.getCDSRegion().getEnd() + 1, getEnd(), getOrientation());
			return rtrn;
		}
		Annotation rtrn=new BasicAnnotation(getChr(), getStart(), this.getCDSRegion().getStart() - 1, getOrientation());
		return rtrn;
	}
	
	/**
	 * Return the chromosome number without the chr
	 * @return the chromosome number without the chr
	 */
	public String getChrNum() {
		String chr=getChr();
		return chr.replaceAll("chr", "");
	}
	
	/**
	 * Get the region of the gene span that lies 3' of the stop codon, including introns but not the stop codon itself
	 * @return the region of the gene span that lies 3' of the stop codon as an Alignments object
	 */
	public Alignments get3UtrIncludingIntrons(){
		if(getOrientation() == Strand.NEGATIVE){
			Alignments rtrn=new Alignments(getChr(), getStart(), this.getCDSRegion().getStart() - 1, getOrientation());
			return rtrn;
		}
		Alignments rtrn=new Alignments(getChr(), this.getCDSRegion().getEnd() + 1, getEnd(), getOrientation());
		return rtrn;
	}
	
	public int getNumExons(){
		return numBlocks();
	}
	
	public Alignments[] getExons(){
		Collection<? extends Annotation> blocks=getBlocks();
		
		Alignments[] exons=new Alignments[blocks.size()];
		
		int count=0;
		for(Annotation block: blocks){
			exons[count++]=new Alignments(block);
		}
		
		return exons;
	}
	
	public Set<? extends Annotation> getExonSet() {
		Set<Annotation> rtrn=new TreeSet<Annotation>();
		
		for(Annotation block: getBlocks()){
			rtrn.add(block);
		}
		return rtrn;
	}
	
	public Collection<? extends Annotation> getSortedAndUniqueExons(){
		return getExonSet();
	}
	
	public IntervalTree<Annotation> getExonTree() {
		IntervalTree<Annotation> rtrn=new IntervalTree();
		
		for(Annotation exon: this.getExonSet()){rtrn.put(exon.getStart(), exon.getEnd(), exon);}
		
		return rtrn;
	}
	
	/**
	 * @return first exon in coordinate order
	 */
	public Annotation getFirstExon(){
		return getBlocks().get(0);
	}
	
	/**
	 * @return last exon in coordinate order
	 */
	public Annotation getLastExon(){
		return getBlocks().get(numBlocks()-1);
	}
	
	public Annotation getOrientedLastExon() {
		if (this.getOrientation() == Strand.NEGATIVE)
			return getFirstExon();
		else
			return getLastExon();
	}
	
	public Annotation get5PrimeExon() {
		if (this.getOrientation() == Strand.NEGATIVE)
			return getLastExon();
		else
			return getFirstExon();
	}
	
	public Annotation get3PrimeExon(int size) {
		//if + return the last exon
		if(getOrientation() == Strand.POSITIVE){
			Annotation rtrn = getLastExon();
			if(rtrn.getSize()>size){
				rtrn=new BasicAnnotation(rtrn.getChr(), rtrn.getEnd()-size, rtrn.getEnd(),Strand.POSITIVE);
			}
			return rtrn;
		}
		//if minus return the first exon
		if(getOrientation() == Strand.NEGATIVE){
			Annotation rtrn = getFirstExon();
			if(rtrn.getSize()>size){
				rtrn=new BasicAnnotation(rtrn.getChr(), rtrn.getStart(), rtrn.getStart()+size, Strand.NEGATIVE);
			}
			return rtrn;
		}
		//else return null;
		return null;
	}
	
	public Annotation get3PrimeExon() {
		if (this.getOrientation() == Strand.NEGATIVE)
			return getFirstExon();
		else
			return getLastExon();
	}

	/**
	 * Get introns
	 * @return Gene object whose blocks are the introns, or null if no introns
	 */
	public Gene getIntrons(){
		Collection<? extends Annotation> introns = this.getIntronSet();
		if(introns.size() == 0) return null;
		Gene rtrn = new Gene(introns);
		rtrn.setOrientation(getOrientation());
		return rtrn;
	}
	
	public Annotation getLastIntron() {
		Gene introns = getIntrons();
		if(introns.getNumExons() == 0) throw new IllegalArgumentException("Gene has no introns");
		return introns.get3PrimeExon();
	}
	
	/**
	 * Returns the first intron of this gene
	 * @return
	 */
	public Annotation getFirstIntron() {
		Gene introns = getIntrons();
		if(introns == null) return introns;
		return introns.get5PrimeExon();
	}
	
	public Alignments[] getIntronsBlocks(){
		
		Object[] sortedUniqueExons=this.getSortedAndUniqueExons().toArray();
		if(sortedUniqueExons.length==0){return new Alignments[0];}
		Alignments[] rtrn=new Alignments[sortedUniqueExons.length-1];
		
		
		for(int i=0; i<sortedUniqueExons.length-1; i++){
			Alignments current=(Alignments)sortedUniqueExons[i];
			Alignments next=(Alignments) sortedUniqueExons[i+1];
			Alignments align=new Alignments(getChr(), current.getEnd(), next.getStart());
			rtrn[i]=align;
		}
		return rtrn;
	}
	
	public Collection<? extends Annotation> getIntronSet(){
		Collection<BasicAnnotation> rtrn=new TreeSet<BasicAnnotation>(); 
		
		List<? extends Annotation> sortedUniqueExons=this.getBlocks();
		
		for(int i=0; i< sortedUniqueExons.size()-1; i++){
			Annotation current = sortedUniqueExons.get(i);
			Annotation next = sortedUniqueExons.get(i+1);
			//Alignments align=new Alignments(this.chr, current.getEnd()+1, next.getStart()-1); // If we assume Exons are in [start,end) format there should be no need add one to the end.  
			BasicAnnotation intron=new BasicAnnotation(getChr(), current.getEnd(), next.getStart());   
			//if(align.getSize()<0){System.err.println(align.toUCSC()+" "+current.toUCSC()+" "+next.toUCSC()+" "+toBED());}
			intron.setOrientation(getOrientation());
			intron.setName(getName() + "_intron_" + i);
			rtrn.add(intron);
		}
		return rtrn;
	}
	
	public Gene getIntronTranscript() {
	
		if (this.getNumExons()<=1)
			return null;

		Gene g =new Gene(getChr(), this.getName(), this.getOrientation(), this.getIntronSet());
		
		/*
		Alignments[] iarr=this.getIntronsBlocks();
		int g_st=iarr[0].getStart();
		int g_end=iarr[iarr.length -1].getEnd();
		Gene g =new Gene(getChr(),g_st,g_end,this.getName(),this.getOrientation(),this.getIntronSet());
		*/
		
		return g;
	}
	
	/*public RefSeqGene getCDS(){
		Alignments cds=getCDSRegion();
		Collection<Alignments> exons=new TreeSet<Alignments>();
		
		
		for(Alignments exon: this.getExonSet()){
			if(exon.overlapsAtAll(cds)){exons.add(exon);}	
		}
		
		Alignments firstExon=exons.iterator().next();
		Alignments lastExon=(Alignments)exons.toArray()[exons.size()-1];
		
		exons.remove(firstExon);
		exons.remove(lastExon);
		
		firstExon=trimFirst(firstExon, cds);
		lastExon=trimLast(lastExon, cds);
		
		exons.add(firstExon);
		exons.add(lastExon);
		
		return new RefSeqGene(exons);
	}*/
	
	/**
	 * Whether the gene is coding or noncoding
	 * @return True if the cds exists and has size at least 3
	 */
	public boolean isCodingGene() {
		Gene cds = getCDS();
		if(cds == null) return false;
		if(cds.size() < 3) return false;
		return true;
	}
	
	public Gene getCDS(){
		Alignments cds=getCDSRegion();
		Gene rtrn=this.trimAbsolute(cds.getStart(), cds.getEnd());		
		if(rtrn == null) {
			return null;
		}
		String geneName = "";
		if(getName() == null) {
			geneName += getChr() + "_" + getStart() + "_" + getEnd();
		} else {
			geneName += getName();
		}
		rtrn.setName(geneName + "_CDS");
		return rtrn;
	}
	
	public Alignments getCDSRegion(){
		Alignments align=new Alignments(getChr(), this.cdsStart, this.cdsEnd);
		return align;
	}
	
	public Annotation getAlignment(){
		Alignments align=new Alignments(getChr(), getStart(), getEnd(), getOrientation(), getName());
		return align;
	}
	
	public Alignments getPromoter(int fudgeFactor){ return getPromoter(fudgeFactor,fudgeFactor);}
	
	public Alignments getPromoter(int upstream,int downstream){
		Alignments align=null;
		if(getOrientation() == Strand.POSITIVE){
			align=new Alignments(getChr(), getStart()-upstream, getStart()+downstream);
		}
		else{
			align=new Alignments(getChr(), getEnd()-downstream, getEnd()+upstream);
		}
		align.setName(getName());
		align.setOrientation(getOrientation().toString());
		return align;
	}
	
	public double getNumberPositions(int rate){
		return ((getEnd()- getStart())/rate)+1;
		
	}
	
	public int getGappedSize(){
		int sum=0;
		Alignments[] exons=this.getExons();
		if(exons!=null){
			for(int i=0; i<exons.length; i++){
				sum+=exons[i].getSize();
			}
		}
		return sum;
	}
	
	public double getBedScore(){return bedScore;}
	
	public int get3PrimeGenomicDistance(LightweightGenomicAnnotation annot) {
		
		int res=0;
		if (getOrientation() == Strand.NEGATIVE)
			res=annot.getEnd()-this.getStart();
		else
			res=this.getEnd()-annot.getStart();
			
	
		return res;
	}
	
		
	public Gene getOverlap(Gene second) {
		Set<? extends Annotation> otherExons = second.getExonSet();
		List<? extends Annotation> thisExons  = new ArrayList<Annotation>(getExonSet());
		Collection<Annotation> overlappingExons = new TreeSet<Annotation>();
		for(Annotation exon: otherExons) {
			Annotation exonClone = new Alignments(exon);
			List<Annotation> tmp = exonClone.intersect(thisExons);
			overlappingExons.addAll(exonClone.intersect(thisExons));
		}
		
		Gene rtrn = null;
		if(!overlappingExons.isEmpty()) {
			rtrn = new Gene(overlappingExons);
			rtrn.setOrientation(getOrientation());
		}
		
		return rtrn;
	}
	
	/**
	 * Get all overlaps between this gene and any gene in the collection
	 * @param others The collection of genes to check for overlap
	 * @return A RefSeqGene object whose "exons" consist of overlapping intervals between this gene and the collection
	 */
	public Gene getOverlap(Collection<Gene> others) {
		Set<Annotation> otherExons = new TreeSet<Annotation>();
		for(Gene other : others) otherExons.addAll(other.getExonSet());
		List<? extends Annotation> thisExons  = new ArrayList<Annotation>(getExonSet());
		Collection<Annotation> overlappingExons = new TreeSet<Annotation>();
		for(Annotation exon: otherExons) {
			Annotation exonClone = new BasicAnnotation(exon);
			overlappingExons.addAll(exonClone.intersect(thisExons));
		}
		
		Gene rtrn = null;
		if(!overlappingExons.isEmpty()) {
			rtrn = new Gene(overlappingExons);
			rtrn.setOrientation(getOrientation());
		}
		
		return rtrn;		
	}
	
	private int getOverlap(Annotation exon, Annotation align) {
		Annotation intersect=exon.intersect(align);
		return intersect.length();
	}
	
	
	//position of splice junction relative to full length sequence (0 index)
	public ArrayList<Integer> getSpliceJunctionCoordinates()throws Exception{
		ArrayList<Integer> coordinates=new ArrayList();
		
		Collection<? extends Annotation> exons=this.getExonSet();
		
		int length=0;
		int counter=0;
		for(Annotation exon: exons){
			length+=exon.getSize();
			counter++;
			if(counter<exons.size()){coordinates.add(length);}
		}
		
		
		return coordinates;
	}
	
	public IntervalTree<Integer> getSpliceJunctionCoordinatesTree() throws Exception {
		ArrayList<Integer> junctions=getSpliceJunctionCoordinates();
		IntervalTree<Integer> tree=new IntervalTree();
		
		for(Integer pos: junctions){
			tree.put(pos, pos+1, pos);
		}
		
		return tree;
	}
	
	
	public String toString(){
		return this.toBED();
		//return this.name;
	}
	
	public String toRefSeq(){
		Set<? extends Annotation> exons = getExonSet();
		String rtrn=getName()+"\t"+(getName() == null ? toUCSC() : getName())+"\t"+getChr()+"\t"+getOrientation()+"\t"+getStart()+"\t"+getEnd()+"\t"+getStart()+"\t"+getEnd()+"\t"+exons.size();
		String starts="";
		String ends="";
		for(Annotation exon : exons) {
			starts+=exon.getStart();
			ends+=exon.getEnd();
		}
		rtrn=rtrn+"\t"+starts+"\t"+ends;
		return rtrn;
	}
	
	public String toUCSC() {
		return getChr()+":" +getStart() + "-" + getEnd();
	}
	
	public String toBED (){
		return toBED(true);
	}
	
	public String toBED(int r, int g, int b) {
		return toBED(true, r, g, b);
	}
	
	public String toBED(boolean useExtraFields){
		return toBED(true, 0, 0, 0);
	}
	
	/*
	 * Can't use method in AbstractAnnotation because need CDS start and end. -PR
	 */
	public String toBED(boolean useExtraFields, int r, int g, int b){
		if(r < 0 || r > 255 || g < 0 || g > 255 || b < 0 || b > 255) {
			throw new IllegalArgumentException("RGB values must be between 0 and 255");
		}
		String rgb = r + "," + g + "," + b;
		List<? extends Annotation> exons = getBlocks();
		String rtrn=getReferenceName()+"\t"+getStart()+"\t"+getEnd()+"\t"+(getName() == null ? toUCSC() : getName())+"\t"+getBedScore()+"\t"+getOrientation()+"\t"+getCDSStart()+"\t"+getCDSEnd()+"\t"+rgb+"\t"+exons.size();
		String sizes="";
		String starts="";
		for(Annotation exon : exons){
			sizes=sizes+(exon.length())+",";
			starts=starts+(exon.getStart()-getStart())+",";
		}
		rtrn=rtrn+"\t"+sizes+"\t"+starts;
		if(extraFields != null & useExtraFields) {
			for(String field : extraFields) {
				rtrn = rtrn + "\t" + field;
			}
		}
		return rtrn;
	}
	
	public String toBEDwithBedScore(){
		return toBED();
	}
		
	public String toGTF(String source) {
		return this.toGTF(source,this.getName(),getName() +".0");
	
	}
	
	public String toGTF(String source,String geneName ,String transcriptName) {
		
		StringBuilder rtrn;
		rtrn = new StringBuilder();
			Set<? extends Annotation> exons = getExonSet();
			GFF transcriptGFF = new GFF(new Gene(getChr(), getStart(), getEnd(), getName()));
			transcriptGFF.setFeature("transcript");
			transcriptGFF.setSource(source);
			transcriptGFF.setOrientation(this.getOrientation().toString());
			transcriptGFF.addAttribute("gene_id", transcriptName);
			if(attributes != null && !attributes.isEmpty()) {
				for(String attr : attributes.keySet()) {
					transcriptGFF.addAttribute(attr, String.valueOf(attributes.get(attr)));
				}
			}
			//Add +1 when we print GFF as it is 1 based rather than 0 based
			transcriptGFF.setStart(getStart());  // BUG FIX : MORAN AUG 17TH, ADD +1 ; 2nd bug fix - 11.29.10 only during print we add +1 to start pos
			transcriptGFF.setEnd(getEnd()); //BUG FIX : MORAN AUG 17TH, END IS CORRECT AS GFF IS INCLUSIVE BUT ALSO USES 1 BASE CORRDINATES (WE USE 0 BASED)
			transcriptGFF.setChromosome(getChr());
			rtrn.append(transcriptGFF.toString(true));
			rtrn.append("\n");
			int i = 0;
			for(Annotation exon : exons) {
				GFF exonGFF = new GFF(exon);
				exonGFF.setFeature("exon");
				exonGFF.setSource(source);
				exonGFF.setName(exon.getName());
				exonGFF.addAttribute("gene_id", exon.getName());
				exonGFF.addAttribute("transcript_id", transcriptName);
				exonGFF.addAttribute("Parent",   transcriptName);
				if(attributes != null && !attributes.isEmpty()) {
					for(String attr : attributes.keySet()) {
						exonGFF.addAttribute(attr, String.valueOf(attributes.get(attr)));
					}
				}
				rtrn.append(exonGFF.toString(true));
				rtrn.append("\n");
			}
	
		return rtrn.toString();
	}
		
	
	
	//MORAN: added GTF in cufflinks format without parent
	//Note: To use the attributes transcript_id or gene_id  as is, set geneName and transcriptNme to empy strings
	//else the name will be concatenated to the attributes values
	public String toCufflinksGTF(String source,String geneName ,String transcriptName,String nonParsedAttrs){
		StringBuilder rtrn = new StringBuilder();
		Set<? extends Annotation> exons = getExonSet();
		GFF transcriptGFF = new GFF(getName());
		transcriptGFF.setSource(source);
		transcriptGFF.setOrientation(this.getOrientation().toString());
		transcriptGFF.setFeature("transcript");//As in version 1353 
		
		if ( nonParsedAttrs.isEmpty()){
			transcriptGFF.clearAttribute("gene_id");
			if (! geneName.isEmpty())
				transcriptGFF.addAttribute("gene_id", geneName);//As in version 1353 to compile with cufflinks
			transcriptGFF.clearAttribute("transcript_id");
			if (! transcriptName.isEmpty())
				transcriptGFF.addAttribute("transcript_id", transcriptName);//
			if(attributes != null && !attributes.isEmpty()) {
				for(String attr : attributes.keySet()) {
					transcriptGFF.addAttribute(attr, String.valueOf(attributes.get(attr)));
				}
			}
		}
		else {
			transcriptGFF.setAttributes(nonParsedAttrs);
		}
		//Add +1 when we print GFF as it is 1 based rather than 0 based
		transcriptGFF.setStart(getStart());  //BUG FIX : MORAN AUG 17TH, ADD +1; 2nd bug fix - 11.29.10 only during print we add +1 to start pos
		transcriptGFF.setEnd(getEnd()); //BUG FIX : MORAN AUG 17TH, END IS CORRECT AS GFF IS INCLUSIVE BUT ALSO USES 1 BASE CORRDINATES (WE USE 0 BASED)
		transcriptGFF.setChromosome(getChr());
		rtrn.append(transcriptGFF.toCufflinksString(true));
		rtrn.append("\n");
		
		for(Annotation exon : exons) {
			GFF exonGFF = new GFF(exon);
			exonGFF.setFeature("exon");
			exonGFF.setSource(source);
			
			if ( nonParsedAttrs.isEmpty()){
				exonGFF.setName(exon.getName());
				exonGFF.clearAttribute("gene_id");
				if (! geneName.isEmpty())
					exonGFF.addAttribute("gene_id", geneName); //As in version 1353 to compile with cufflinks
				exonGFF.clearAttribute("transcript_id");
				if (! transcriptName.isEmpty())
					exonGFF.addAttribute("transcript_id", transcriptName);
				//exonGFF.addAttribute("Parent",   transcriptName);
				if(attributes != null && !attributes.isEmpty()) {
					for(String attr : attributes.keySet()) {
						exonGFF.addAttribute(attr, String.valueOf(attributes.get(attr)));
					}
				}
			}
			else
				exonGFF.setAttributes(nonParsedAttrs);
			
			rtrn.append(exonGFF.toCufflinksString(true));
			rtrn.append("\n");
		}
	
		return rtrn.toString();
	}
	
	//Fixed off by one error
	public String toSAM(){
		String qualityString="255";
		String qualityLetters="";
		String mateReference="*";
		String matePosition="0";
		String insertSize="0";
		
		int strandNum=0;
		if(getOrientation() == Strand.NEGATIVE){strandNum=16;}
		if(this.sequence==null){this.sequence="";}
		
		for(int i=0; i<sequence.toCharArray().length; i++){qualityLetters=qualityLetters+"I";}
		
		String cigar=toCigar();
		
		String rtrn=getName()+"\t"+strandNum+"\t"+getChr()+"\t"+getSAMStart()+"\t"+qualityString+"\t"+cigar+"\t"+mateReference+"\t"+matePosition+"\t"+insertSize+"\t"+sequence+"\t"+qualityLetters;
		
		return rtrn;
	}
	
	public String toSAM(String originalSAMLine) {
		return toSAM(originalSAMLine, false);
	}
	
	public String toSAM(String originalSAMLine, boolean flip) {
		String[] tokens=originalSAMLine.split("\t");
		//all we want to do is update the relevant columns which are the chromosome position and cigar string
		//orientation is [1]
		//chr is [2]
		//position is [3]
		//cigar is [5]
		//sequence is [9]
		//quality is [10]
		
		String sequence=tokens[9];
		//String quality=tokens[10];
		String orientation=tokens[1];
		if(flip){
			if(orientation.equalsIgnoreCase("0")){orientation="16";}
			else if(orientation.equalsIgnoreCase("16")){orientation="0";}
			else{System.err.println(orientation);}
			sequence=Sequence.reverseSequence(sequence);
			//quality=reverse(quality);
		}
		 
		String rtrn="";
		for(int i=0; i<tokens.length; i++){
			if(i!=1 && i!=2 && i!=3 && i!=5 && i!=9){rtrn+=tokens[i]+"\t";}
			else if(i==1){rtrn+=orientation+"\t";}
			else if(i==2){rtrn+=this.getChr()+"\t";}
			else if(i==3){rtrn+=(getSAMStart())+"\t";}
			else if(i==5){String cigar=toCigar(); rtrn+=cigar+"\t";}
			else if(i==9){rtrn+=sequence+"\t";}
			//else if(i==10){rtrn+=quality+"\t";}
		}
		return rtrn;
	}
	
	public void setCDS(Annotation cdsRegion){
		this.cdsStart=cdsRegion.getStart();
		this.cdsEnd=cdsRegion.getEnd();
	}
	
	@Deprecated  // should use getPositionAtReferenceCoordinate instead
	public int absoluteToRelativePosition(int referenceCoordinate) {
		return getPositionAtReferenceCoordinate(referenceCoordinate);
		/*
		Collection<? extends Annotation> exons=this.getExonSet();
		
		int count=0;
		for(Annotation exon: exons){
			if(exon.fullyContained(new Alignments(getChr(), position, position))){
				count+=(position-exon.getStart()-1);
			}
			else if(exon.getEnd()<position){count+=exon.length();}
		}
		
		return count;
		*/
	}
	
	public boolean isInExon(int snp) {
		Collection<? extends Annotation> exons=this.getExonSet();
		
		
		for(Annotation exon: exons){
			if(exon.getStart()<snp && exon.getEnd()>snp){return true;}
		}
		return false;
	}
	
	/**
	 * Get all ORFs as gene objects
	 * @param trim3UTRs whether to trim 3'UTRs to end at the beginning of next ORF
	 * @return the collection of ORFs
	 */
	public Collection<Gene> findAllORFs(boolean trim3UTRs) {
		
		if(trim3UTRs) return findAllOrfsConservative3UTRs();
		
		Collection<int[]> orfs = findAllORFs(sequence);
		String cdsOrientation = getOrientation().toString();
		if(Strand.UNKNOWN.equals(getOrientation())) {
			System.err.println("Not finding ORFs in gene " + getName() + " because strand is unknown.");
		}
		
		Collection<Gene> rtrn=new TreeSet<Gene>();
		for(int[] orf: orfs){
			int orfGenomicStart = "-".equals(cdsOrientation) ? mapToGenomic(orf[1]):  mapToGenomic(orf[0])  ; 
			int orfGenomicEnd   = "-".equals(cdsOrientation) ? mapToGenomic(orf[0]):  mapToGenomic(orf[1]);
			
			//System.err.println("ORF genomic start end " + orfGenomicStart +"-" + orfGenomicEnd);
			
			LightweightGenomicAnnotation geneCDS = new BasicLightweightAnnotation(getChr(), orfGenomicStart, orfGenomicEnd);
			Set<Annotation> cdsExons = new TreeSet<Annotation>();
			Set<? extends Annotation> exonSet = getExonSet();
			for(Annotation e : exonSet) {
				if(geneCDS.overlaps(e)) {
					cdsExons.add(e.intersect(geneCDS));
				}
			}
			
			Gene cds = null;
			if(cdsExons.size() > 0 ) {
				cds = new Gene(cdsExons);
				cds.setOrientation(cdsOrientation);
				cds.setName(getName());
				cds.setCountScore((orf[1]- orf[0])/(double)(sequence.length()));
			}
			
			//TODO Add 5' UTR and 3' UTR to each ORF
			//cds.setCDS(cds.getAlignment());
			Gene testCDS=this.copy();
			testCDS.setCDS(cds.getAlignment());
						
			//rtrn.add(cds);
			rtrn.add(testCDS);
		}
		return rtrn;
	}
	
	/**
	 * Get all ORFs as genes, with 3'UTRs trimmed so they end at the beginning of next downstream ORF
	 * @return The collection of all ORFs
	 */
	private Collection<Gene> findAllOrfsConservative3UTRs() {
		
		if(Strand.UNKNOWN.equals(getOrientation())) {
			System.err.println("Not finding ORFs in gene " + getName() + " because strand is unknown.");
		}
		
		//System.err.println("Finding ORFs with shortened 3'UTRs for gene " + name + " " + chr + ":" + start + "-" + stop);
		
		TreeSet<Gene> orfs = (TreeSet<Gene>) findAllORFs(false);
		
		Collection<Gene> rtrn = new TreeSet<Gene>();
		
		if(getOrientation().equals(Strand.POSITIVE)) {
			
			Iterator<Gene> iter = orfs.iterator();
			if(!iter.hasNext()) return rtrn;
			
			while(iter.hasNext()) {
				Gene orf = iter.next();
				int cdsEnd = orf.cdsEnd;
				int endPos = orf.getEnd();
				Iterator<Gene> after = orfs.tailSet(orf).iterator();
				while(after.hasNext()) {
					int nextStart = after.next().cdsStart;
					if(nextStart > cdsEnd) {
						endPos = nextStart;
						break;
					}
				}
				//System.err.println("Setting end from " + orf.stop + " to " + endPos);
				Gene trimmed = orf.trimAbsolute(orf.getStart(), endPos);			
				//System.err.println("New end is " + trimmed.stop);
				if(trimmed != null) {
					rtrn.add(trimmed);
					//System.err.println("Added ORF with CDS " + trimmed.blockStart + "-" + trimmed.blockEnd);
				}
			}
		}
		
		if(getOrientation().equals(Strand.NEGATIVE)) {
			Iterator<Gene> iter = orfs.descendingIterator();
			if(!iter.hasNext()) return rtrn;
			while(iter.hasNext()) {
				Gene orf = iter.next();
				int cdsStart = orf.cdsStart;
				int startPos = orf.getStart();
				Iterator<Gene> after = ((TreeSet) orfs.headSet(orf)).descendingIterator();
				while(after.hasNext()) {
					int nextEnd = after.next().cdsEnd;
					if(nextEnd < cdsStart) {
						startPos = nextEnd;
						break;
					}
				}
				//System.err.println("Setting start from " + orf.start + " to " + startPos);
				Gene trimmed = orf.trimAbsolute(startPos, orf.getEnd());			
				//System.err.println("New start is " + trimmed.start);
				if(trimmed != null) {
					rtrn.add(trimmed);
					//System.err.println("Added ORF with CDS " + orf.blockStart + "-" + orf.blockEnd);
				}
			}
		}
		
		return rtrn;
		
	}
	
	/**
	 * Find start codons
	 * @return collection of all start codon genomic coordinates
	 */
	public TreeSet<BasicLightweightAnnotation> findAllStartCodons() {
		Collection<int[]> orfs = findAllORFs(sequence);
		String cdsOrientation = getOrientation().toString();
		if(Strand.UNKNOWN.equals(getOrientation())) {
			System.err.println("Skipping because we dont know strand");
			return null;
		}
		
		TreeSet<BasicLightweightAnnotation> rtrn = new TreeSet<BasicLightweightAnnotation>();
		for(int[] orf: orfs){
			int startCodonGenomicStart = "-".equals(cdsOrientation) ? mapToGenomic(orf[1]):  mapToGenomic(orf[0])  ; 
			int startCodonGenomicEnd   = "-".equals(cdsOrientation) ? mapToGenomic(orf[1] + 2):  mapToGenomic(orf[0] + 2);
			
			BasicLightweightAnnotation startCodon = new BasicLightweightAnnotation(getChr(), startCodonGenomicStart, startCodonGenomicEnd);
			rtrn.add(startCodon);
			System.err.println("Found start codon at " + getChr() + " " + startCodonGenomicStart + "-" + startCodonGenomicEnd);
		}
		return rtrn;
		
	}
	
	
	public static Collection<int[]> findAllORFs(String sequence) {
		Collection<int[]> allORFs=new ArrayList<int[]>();
		int lastStart = 0;
		int lastEnd = 0;
		Matcher m = START_CODON_PATRN.matcher(sequence);
		while(m.find()) {
			int startCodonPos = m.start(); //start codon
			int thisORFEnd = startCodonPos;
			boolean foundStopCodon = false;
			//System.err.println("new orf start: " + startCodonPos);
			while(thisORFEnd < sequence.length() - 3  && !foundStopCodon) {
				String currentCodon = sequence.substring(thisORFEnd, thisORFEnd+3);
				//System.err.print(" "+currentCodon+ " ");
				thisORFEnd += 3;
				foundStopCodon = isStopCodon(currentCodon);
				//int[] pos={startCodonPos, thisORFEnd};
				//allORFs.add(pos);
			}
			if(foundStopCodon){
				int[] pos={startCodonPos, thisORFEnd};
				allORFs.add(pos);
			}
			//System.err.println("Orf length: " + (thisORFEnd - startCodonPos) );
			if(lastEnd - lastStart < thisORFEnd - startCodonPos) {
				lastStart = startCodonPos;
				lastEnd   = thisORFEnd;
				//System.err.println("It was a winner");
			}
		}
		int [] startEnd = {lastStart, lastEnd};
		return allORFs;
	}
	
	
	
	public Collection<Gene> findAllORFs(Sequence chrSeq, boolean trim3UTRs) {
		this.setSequenceFromChromosome(chrSeq);
		return findAllORFs(trim3UTRs);
	}
	
	
	public void addSequence(Sequence seq){
		this.sequence=seq.getSequenceBases();
	}
	
	public void addSuffixToName(String refName) {
		this.setName(this.getName()+refName);
		
	}
	
	public void addAttribute(String name, String val){
		if(attributes == null) {
			attributes = new HashMap<String, String>();
		}
		attributes.put(name, val);
	}
	
	public void addExtraField (String value) {
		String [] newExtraFields = null;
		if(extraFields != null) {
			newExtraFields = new String[extraFields.length + 1];
			for(int i = 0; i < extraFields.length; i++) {
				newExtraFields[i] = extraFields[i];
			}
		} else {
			newExtraFields = new String [1];
		}
		
		newExtraFields[newExtraFields.length - 1] = value;
		
		this.extraFields = newExtraFields;
	}
	
	public void expandUtrs(Integer utr1, Integer utr2) {
		this.setStart(getStart()-utr1);
		this.setEnd(getEnd()+utr2);
	}
		
	
	public boolean geneSpanContains(Gene gene){
		return gene.getStart() >= getStart() && gene.getEnd() <= getEnd();
	}
	
	public boolean containsAttribute(String key) {
		return attributes != null && attributes.containsKey(key);
	}
	
	public static boolean isStopCodon(String currentCodon) {
		boolean isStopCodon = false;
		for( String stop : STOP_CODONS) {
			if(currentCodon.equalsIgnoreCase(stop)) {
				isStopCodon = true;
				break;
			}
		}
		return isStopCodon;
	}
	
	
	public boolean isFullyCompatible(Gene iso) {
	
		boolean res=false;
		Gene myIntron = getIntronTranscript();
		Gene isoIntron=iso.getIntronTranscript();
		if (myIntron==null | isoIntron==null)
			return false;
		if (myIntron.isEqualStructure(isoIntron))
			res=true;
		return res;
	}
	
	private boolean isEqualStructure(Gene iso) {
	
		boolean res=true;
		if (this.getNumExons()!=iso.getNumExons())
			return false;
		Alignments[] myEx= this.getExons();
		Alignments[] isoEx= iso.getExons();
		for (int i=0; i<this.getNumExons();i++){
			if (!(myEx[i].equals(isoEx[i])))
				res=false;
				
		}
		
		return res;
	}
	
	public boolean overlapsExon(Annotation gen) {
		Collection<? extends Annotation> exonSet=this.getExonSet();
		for(Annotation align: exonSet){
			if(align.overlaps(new BasicGenomicAnnotation(gen.getName(),gen.getChr(),gen.getStart(),gen.getEnd())) ) {return true;}
		}
		return false;
	}
	
	public boolean overlapsExon(LightweightGenomicAnnotation exon){
		Collection<? extends Annotation> exonSet=this.getExonSet();
		
		for(Annotation align: exonSet){
			if(align.overlaps(exon)){return true;}
		}
		
		return false;
	}
	
	/**
	 * Whether an exon of this gene overlaps an exon of other gene and the two genes have the same orientation
	 * @param other Other gene
	 * @return Whether genes have same orientation and at least one overlapping exon
	 */
	public boolean overlaps(Gene other) {
		return overlaps(other, false);
	}
	
	/**
	 * Get all genes overlapping a window
	 * @param genes The genes
	 * @param chr Window chr
	 * @param start Window start
	 * @param end Window end
	 * @return All genes overlapping the window, not necessarily fully contained
	 */
	public static Collection<Gene> getOverlappers(Collection<Gene> genes, String chr, int start, int end, Strand strand, boolean ignoreOrientation) {
		Gene window = new Gene(new BasicAnnotation(chr, start, end, strand));
		Collection<Gene> rtrn = new TreeSet<Gene>();
		for(Gene gene : genes) {
			if(gene.overlaps(window, ignoreOrientation))	{
				rtrn.add(gene);
			}
		}
		return rtrn;
	}
	
	/**
	 * Get inner distance between this gene and its nearest neighbor in a collection
	 * Ignore genes that overlap this gene in any orientation
	 * Looks at entire gene span, not just exons
	 * @param others Other genes to check
	 * @return Min distance to a non-overlapper in the collection, or -1 if no other genes on same chromosome
	 */
	public int distanceToNearestNonOverlapper(Collection<Annotation> others) {
		TreeSet<Gene> thisChr = new TreeSet<Gene>();
		for(Annotation g : others) {
			if(g.getChr().equals(getChr())) {
				thisChr.add(new Gene(g.getChr(), g.getStart(), g.getEnd(), g.getName()));
			}
		}
		if(thisChr.isEmpty()) {
			return -1;
		}
		int minDist = Integer.MAX_VALUE;
		Iterator<Gene> headSetDescending = ((TreeSet<Gene>) thisChr.headSet(this)).descendingIterator();
		while(headSetDescending.hasNext()) {
			Gene g = headSetDescending.next();
			if(overlaps(g, true)) {
				continue;
			}
			int dist = getStart() - g.getEnd();
			if(dist < minDist) {
				minDist = Math.max(0, dist);
				break;
			}
		}
		Iterator<Gene> tailSetAscending = ((TreeSet<Gene>) thisChr.tailSet(this)).iterator();
		while(tailSetAscending.hasNext()) {
			Gene g = tailSetAscending.next();
			if(overlaps(g, true)) {
				continue;
			}
			int dist = g.getStart() - getEnd();
			if(dist < minDist) {
				minDist = dist;
				break;
			}
		}
		if(minDist == Integer.MAX_VALUE) {
			throw new IllegalStateException("No min distance found");
		}
		return minDist;
	}
	
	/**
	 * Whether an exon of this gene overlaps an exon of other gene
	 * @param other Other gene
	 * @param ignoreOrientation Ignore orientation. If set to false, orientations must be equal or at least one orientation must be unknown.
	 * @return Whether genes overlap at the exon level
	 */
	public boolean overlaps(Gene other, boolean ignoreOrientation) {
		//System.err.println("Findig overlap betweein this " + toBED() + "\nand\n"+other.toBED());
		boolean overlaps = false;
		//System.err.println("\t\t this gene's orientation [" + orientation + "] others [" + other.getOrientation()+ "]");
		if(ignoreOrientation || getOrientation().equals(other.getOrientation()) || Strand.UNKNOWN.equals(getOrientation()) || Strand.UNKNOWN.equals(other.getOrientation())) {
			//System.err.println("\t\tOrientation is compatible");
			Collection<? extends Annotation> exons = getSortedAndUniqueExons();
			Collection<? extends Annotation> otherExons = other.getSortedAndUniqueExons();
			Iterator<? extends Annotation> otherExonsIt  = otherExons.iterator();
			
			while (!overlaps && otherExonsIt.hasNext()) {
				Annotation otherExon = otherExonsIt.next();
				Iterator<? extends Annotation> exonIt  = exons.iterator();
				while(!overlaps && exonIt.hasNext()) {
					Annotation exon = exonIt.next();
					overlaps = exon.overlaps(otherExon);
					if(exon.compareTo(otherExon) > 0) {
						break;
					}
					//System.err.println("exon " + exon.toUCSC() + " otherExon " + otherExon.toUCSC() + " overlaps? " + overlaps + " exons size " + exons.size() + " other exons size " + otherExons.size());
				}
			}
		}
		return overlaps;
	}
	
	/**
	 * Returns true if this and the other gene overlap but at least minPctOverlap
	 * @param other The other gene with which overlap is checked
	 * @param minPctOverlap minimum percent of overlap
	 * @return
	 */
	public boolean overlaps(Gene other, double minPctOverlap) {

		boolean compatible = overlaps(other);
		
		if(compatible) {
			double pctOverlap = Math.min(percentOverlapping(other), other.percentOverlapping(this));
			compatible = pctOverlap >= minPctOverlap;
		}
		return compatible;
	}
	
	/**
	 * Returns true if this and the other gene overlap but at least minPctOverlap
	 * @param other The other gene with which overlap is checked
	 * @param minPctOverlap minimum percent of overlap
	 * @return
	 */
	public boolean overlapsStranded(Gene other, double minPctOverlap) {

		boolean compatible = (overlaps(other) && this.getOrientation().equals(other.getOrientation()));
		
		if(compatible) {
			double pctOverlap = Math.min(percentOverlapping(other), other.percentOverlapping(this));
			compatible = pctOverlap >= minPctOverlap;
		}
		return compatible;
	}
	/**
	 * Whether this gene overlaps any gene in the collection at the exon level
	 * @param others The genes to check for overlap
	 * @return True iff this gene overlaps at least one gene in the collection at the exon level
	 */
	public boolean overlapsGene(Collection<Gene> others) {
		for(Gene other : others) {
			if(overlaps(other)) return true;
		}
		return false;
	}
	
	/*private String makeCigar(){
		String rtrn="";
		Alignments[] exons=this.getExons();
		for(int i=0; i<exons.length; i++){
			Alignments exon=(Alignments)exons[i];
			if(i!=0){
				Alignments previous=(Alignments)exons[i-1];
				Alignments intron=new Alignments(exon.getChr(), previous.getEnd(), exon.getStart());
				if(intron.getSize()>0){rtrn=rtrn+(intron.getSize())+"N";}
			}
			rtrn=rtrn+(exon.getSize())+"M";
		}
		return rtrn;
	}*/
	
	//asks if the whole genomic regions of the two genes overlap in the same orientation
	//if one gene does have a specified orientation, assume they 2 genes are in the same orientation
	public boolean overlapsGene(Gene gen) {
		if (  this.getOrientation().equals(gen.getOrientation()) || this.getOrientation() == Strand.UNKNOWN || gen.getOrientation() == Strand.UNKNOWN) 
		{	
		   Alignments g1 =new Alignments(this.getChr(),this.getStart(),this.getEnd());
		   Alignments g2 =new Alignments(gen.getChr(),gen.getStart(),gen.getEnd());
		   return g1.overlaps(g2);
		}
	   else
		   return false;
	}
	
		
	public boolean overlapsGeneInAnyOrientation(Gene otherGene) {
		   Alignments g1 =new Alignments(this.getChr(),this.getStart(),this.getEnd());
		   Alignments g2 =new Alignments(otherGene.getChr(),otherGene.getStart(),otherGene.getEnd());
		   return g1.overlaps(g2);
	}
	
	public boolean overlaps2(LightweightGenomicAnnotation next){
		if(next==null){return false;}
		if(!getChr().equals(next.getChromosome())){return false;}
		if(next.getStart()>=getStart() && next.getStart()<=getEnd()){return true;}
		if(next.getEnd()>=getStart() && next.getEnd()<=getEnd()){return true;}
		return false;
	}
	
	
	public boolean almostEqual(Gene other, int exonFudgeFactor) {
		List<Annotation> otherExons = new ArrayList<Annotation>(other.getExonSet());
		List<Annotation> thisExons = new ArrayList<Annotation>(getExonSet());
		boolean areAlmostEqual = otherExons.size() == thisExons.size();
		
		int i = 0;
		while( areAlmostEqual && i < thisExons.size()) {
			Annotation thisExon = thisExons.get(i);
			Annotation otherExon = otherExons.get(i);
			areAlmostEqual = Math.abs(thisExon.getStart() - otherExon.getStart()) < exonFudgeFactor && Math.abs(thisExon.getEnd() - otherExon.getEnd()) < exonFudgeFactor;
			i++;
		}
		
		return areAlmostEqual;
	}
	
	public boolean almostSameStructure(Gene other, int intronFudgeFactor) {
		List<Annotation> otherIntrons = new ArrayList<Annotation>(other.getIntronSet());
		List<Annotation> thisIntrons = new ArrayList<Annotation>(getIntronSet());
		boolean almostSameStructure = otherIntrons.size() == thisIntrons.size();
		
		int i = 0;
		while( almostSameStructure && i < thisIntrons.size()) {
			Annotation thisIntron = thisIntrons.get(i);
			Annotation otherIntron = otherIntrons.get(i);
			almostSameStructure = Math.abs(thisIntron.getStart() - otherIntron.getStart()) < intronFudgeFactor && Math.abs(thisIntron.getEnd() - otherIntron.getEnd()) < intronFudgeFactor;
			i++;
		}
		
		return almostSameStructure;
	}
	
	public boolean almostContains(Gene other, int exonFudgeFactor) {
		return almostContainsRefSeq(other, exonFudgeFactor);
	}
	
	protected boolean almostContainsRefSeq(Gene other, int exonFudgeFactor) {
		if(!other.getOrientation().equals(getOrientation())) {
			return false;
		}
		List<Annotation> otherExons = new ArrayList<Annotation>(other.getExonSet());
		List<Annotation> thisExons = new ArrayList<Annotation>(getExonSet());
		boolean isAlmostContained = otherExons.size() <= thisExons.size();
		
		if(getOrientation() == Strand.NEGATIVE) {
			int i = otherExons.size() - 1;
			int j = thisExons.size() - 1;
			while( isAlmostContained && i >=0) {
				Annotation thisExon = thisExons.get(j);
				Annotation otherExon = otherExons.get(i);
				isAlmostContained = Math.abs(thisExon.getStart() - otherExon.getStart()) < exonFudgeFactor && Math.abs(thisExon.getEnd() - otherExon.getEnd()) < exonFudgeFactor;
				i--;
				j--;
			}
		} else {
			int i = 0;
			while( isAlmostContained && i < otherExons.size()) {
				Annotation thisExon = thisExons.get(i);
				Annotation otherExon = otherExons.get(i);
				isAlmostContained = Math.abs(thisExon.getStart() - otherExon.getStart()) < exonFudgeFactor && Math.abs(thisExon.getEnd() - otherExon.getEnd()) < exonFudgeFactor;
				i++;
			}
		} 
		
		return isAlmostContained;
	}
	
	
	public int hashCode() {
		return toBED().hashCode();
	}
		
	//sort by midpoint
	@Override
	public int compareTo (Annotation other) {
		int result = super.compareTo(other);
		//Changed by MG (12/24/12) to make it
		if (result==0 && other.getClass().isInstance(Gene.class)) {
			Gene b = (Gene) other;
			List<Annotation> aExons =new ArrayList<Annotation> (getExonSet());
			List<Annotation> bExons =new ArrayList<Annotation> (b.getExonSet());
			
			int minLength = Math.min(aExons.size(), bExons.size());
			int idx = 0;
			while(idx < minLength && result == 0) {
				result = aExons.get(idx).compareTo(bExons.get(idx));
				idx++;
			}
			if(result == 0) {
				result = aExons.size() - bExons.size();
			}
			//check for location of the cds
			if(result ==0){
				result=getCDSRegion().compareTo(b.getCDSRegion());
			}
			result=getCDSRegion().compareTo(b.getCDSRegion()); 
		}
		return result;
	}
	
	public Gene copy(){
		return new Gene(this);
	}
	
	
	/**
	 * Trims the gene in relative space
	 * Returns a GeneWindow with pointer to the original
	 */
	public GeneWindow trimGene(int relativeStart, int relativeEnd){
		//on first call, cache the relative to absolute coordinates
		int absoluteStart=this.getReferenceCoordinateAtPosition(relativeStart, true);
		int absoluteEnd=this.getReferenceCoordinateAtPosition(relativeEnd, true);
		
		//logger.info(relativeStart+"-"+relativeEnd+" "+absoluteStart+"-"+absoluteEnd);
		
		return trimAbsolute(absoluteStart, absoluteEnd);
		
		/*GeneWindow window;
		
		if(getSortedAndUniqueExons().size()>1){
		
			Collection<Annotation> gaps=new TreeSet<Annotation>(); //introns in relative space
			for(Annotation intron: getIntronSet()){
				Annotation relativeIntron=new Alignments(getChr(), intron.getStart()-this.getStart(), intron.getEnd()-this.getStart());
				gaps.add(relativeIntron);
			}
				
			//System.err.println(gaps);
			GenomeWithGaps2 gwg=new GenomeWithGaps2(gaps, getChr(), getAlignment().getSize(), getOrientation().toString());
			
			Gene relativeWindow = gwg.getRelativeWindow(relativeStart, relativeEnd);
			Gene relativeGene = null; 
			if(relativeWindow != null) {	//TODO: REMOVE THIS CHECK ONLY TEMPORARY	
				relativeGene= relativeWindow.addConstantFactor(this.getStart());
			}
			
			window=new GeneWindow(relativeGene);
		}
		
		else{
			Alignments exon=new Alignments(getChr(), this.getStart()+(relativeStart), this.getStart()+relativeEnd, getOrientation(), getName());
			window=new GeneWindow(new Gene(exon));
		}
		
		window.addSourceAnnotation(this);
		window.setOrientation(getOrientation());
		window.setCDSStart(getCDSStart());
		window.setCDSEnd(getCDSEnd());
		window.setName(getName()+"_"+relativeStart+"_"+relativeEnd);
		return window;
		*/
	}
	
		
	

	public void setCDSStart(int start){
		this.cdsStart=start;
	}
	
	public void setCDSEnd(int end){
		this.cdsEnd=end;
	}
	
	public Gene findLongestORF() {
		int [] orf = findLongestORF(sequence);
		String cdsOrientation = getOrientation().toString();
		if(Strand.UNKNOWN.equals(getOrientation())) {
			int [] reverseOrf = findLongestORF(Sequence.reverseSequence(sequence));
			cdsOrientation = "+";
			if(reverseOrf[1] - reverseOrf[0] > orf[1]- orf[0]) {
				orf[0] = reverseOrf[0];
				orf[1] = reverseOrf[1];
				cdsOrientation = "-";
			}
		}
		
		int orfGenomicStart = "-".equals(cdsOrientation) ? mapToGenomic(orf[1]):  mapToGenomic(orf[0])  ; 
		int orfGenomicEnd   = "-".equals(cdsOrientation) ? mapToGenomic(orf[0]):  mapToGenomic(orf[1]);
		
		//System.err.println("ORF genomic start end " + orfGenomicStart +"-" + orfGenomicEnd);
		
		LightweightGenomicAnnotation geneCDS = new BasicLightweightAnnotation(getChr(), orfGenomicStart, orfGenomicEnd);
		Set<Annotation> cdsExons = new TreeSet<Annotation>();
		Set<? extends Annotation> exonSet = getExonSet();
		for(Annotation e : exonSet) {
			if(geneCDS.overlaps(e)) {
				cdsExons.add(e.intersect(geneCDS));
			}
		}
		
		Gene cds = null;
		if(cdsExons.size() > 0 ) {
			cds = new Gene(cdsExons);
			cds.setOrientation(cdsOrientation);
			cds.setName(getName());
			cds.setCountScore((orf[1]- orf[0])/(double)(sequence.length()));
		}
		return cds;
	}
	
	
	public int numOverlappingExons(Gene iso) {
		int res=0;
		for (Annotation myEx: this.getExonSet()){
			for (Annotation isoEx: iso.getExonSet()){
				if (myEx.overlaps(isoEx)){
					res++; break;
				}
			}
		}
		return res;
	}
	
	public int numOfCompatibleIntrons(Gene g) {
		int rtrn=0;
		Alignments[] myIntrons= this.getIntronsBlocks();
		Alignments[] gIntrons= g.getIntronsBlocks();
		int myNan=-10;
		int match=0;
		int[]ix=new int[gIntrons.length];
		for(int i=0; i<gIntrons.length; i++){
			ix[i]=myNan;
			for(int j=0; j<myIntrons.length; j++){
				if (gIntrons[i].equals(myIntrons[j])){
					ix[i]=j;
					match++;
				}
			}
		}
		if (match==0)
			return 0;
		int currLength=0;
		double Prev=-1;
		for (int i=0;i<ix.length;i++){
			if(ix[i]-Prev==1 | (Prev==-1 & ix[i]!= myNan)){ //continue chain, or start chain 
				Prev=ix[i];
				currLength++;
			}
			else
			{
				if (currLength > rtrn)
					rtrn=currLength;
				currLength=0;
				Prev=-1;
			}
		}
		if (currLength > rtrn)
			rtrn=currLength;
		return rtrn;
	}

	//TODO: MG Testing
	//@return A RefSeqGene that represent the trimmed piece if the region overlaps any exonic regions
	public GeneWindow trimAbsolute(int alignmentStart, int alignmentEnd) {
		Annotation rtrn=this.intersect(new Alignments(this.getChr(), alignmentStart, alignmentEnd));
		if (rtrn == null) return null;
		
		/*
		//If ends requested exceed the gene then trim
		if(alignmentStart<this.getStart()){alignmentStart=this.getStart();}
		if(alignmentEnd>this.getEnd()){alignmentEnd=this.getEnd();}
		
		Collection<Annotation> rtrn=new TreeSet<Annotation>();
		//get exons overlapping from start to end
		Alignments region=new Alignments(getChr(), alignmentStart, alignmentEnd);
		
		//go through each exon
		Collection<? extends Annotation> exons=this.getSortedAndUniqueExons();
		for(Annotation exon: exons){
			if(exon.overlaps(region)){
				rtrn.add(exon);
			}
		}

		if(rtrn.isEmpty()){return null;}

		//now trim by ends
		//Get first and last exons
		Annotation first=(rtrn.iterator().next());
		Annotation last=(Annotation)rtrn.toArray()[rtrn.size()-1];
				
		//remove them from the current collection
		rtrn.remove(first);
		rtrn.remove(last);
		
		//If alignment start and end exceeds the exon then reset the bounds to avoid extending
		if(alignmentStart<first.getStart()){alignmentStart=first.getStart();}
		if(alignmentEnd>last.getEnd()){alignmentEnd=last.getEnd();}
		
		if(first.equals(last)){
			//trim first starting at alignmentStart
			Alignments newFirst=new Alignments(first.getChr(), alignmentStart, alignmentEnd);
			rtrn.add(newFirst);
		}
		else{
			//trim first starting at alignmentStart
			Alignments newFirst=new Alignments(first.getChr(), alignmentStart, first.getEnd());
			//trim last ending alignmentEnd
			Alignments newEnd=new Alignments(last.getChr(), last.getStart(), alignmentEnd);
			
			//add new first and last to the collection 
			rtrn.add(newFirst);
			rtrn.add(newEnd);
		}
		
		//make and return RefSeqGene
		//System.err.println(rtrn);
		*/
		
		GeneWindow rtrnGene = new GeneWindow(rtrn);
		rtrnGene.setCDS(this.getCDSRegion());
		rtrnGene.addSourceAnnotation(this);
		return rtrnGene;
	}
	
	
	//computes the percentage of this that overlaps with gene2
	public double percentOverlapping(Gene gene2){
		if(!getChr().equalsIgnoreCase(gene2.getChr())){return 0.0;}
		IntervalTree<Annotation> test=makeExonTree();
		
		Collection<? extends Annotation> exons=gene2.getExonSet();
		int totalOverlap=0;
		for(Annotation exon: exons){
			Iterator<Node<Annotation>> overlappers=test.overlappers(exon.getStart(), exon.getEnd());
			while(overlappers.hasNext()){
				Annotation align=overlappers.next().getValue();
				int overlap=getOverlap(exon, align);
				totalOverlap+=overlap;
			}
		}
		return (double)totalOverlap/getSize();
	}
	
	//computes the percentage of the genomic region spanned by this object that overlaps with that of gene2
	public Double percentGenomeOverlapping(Gene gene2) {
		
		if(!getChr().equalsIgnoreCase(gene2.getChr())){return 0.0;}
		Alignments g1 = new Alignments(this.getChr(),this.getStart(),this.getEnd());
		Alignments g2 = new Alignments(gene2.getChr(),gene2.getStart(),gene2.getEnd());
		int overlap=getOverlap(g1, g2);
		return (double)overlap/this.getGenomicLength();
		
	}
	
	private IntervalTree<Annotation> makeExonTree() {
		IntervalTree<Annotation> tree=new IntervalTree<Annotation>();
		
		Collection<? extends Annotation> exons= getExonSet();
		for(Annotation exon: exons){
			tree.put(exon.getStart(), exon.getEnd(), exon);
		}
		
		return tree;
	}
	
		
	private Gene addConstantFactor(int factor){
		Collection<Annotation> set=new TreeSet();
		for(Annotation exon: this.getExonSet()){
			Alignments abs=new Alignments(exon.getChr(), exon.getStart()+factor, exon.getEnd()+factor);
			set.add(abs);
		}
		return new Gene(set, getName(), getOrientation(), getCDSStart(), getCDSEnd());
	}
	
	/**
	 * Maps cDNA coordinates to genomic coordinates
	 */
	private int mapToGenomic(int position) {
		return getReferenceCoordinateAtPosition(position);
	}
	
	public int transcriptToGenomicPosition(int transcriptPosition) {
		int relativePosition= transcriptPosition;
		int genomicPosition = -1;
		if(isNegativeStrand()){
			relativePosition=(getSize())-(transcriptPosition) - 1; //Recall that the transcript is [0,L) closed open from 0 to the length L
		}
		Gene samRecord= trimGene(relativePosition, relativePosition+1 ); //TODO: Should we do relativePosition-1, relativePosition when transcript in reversed orientation?
		if(samRecord != null) {
			genomicPosition = samRecord.getStart() ;
		} else {
			Exception t = new Exception ("n");
			logger.warn("Could not map back to genome position relative position was: "+ relativePosition+ " this transicript: " +toBED() + " stack trace: ", t );
		}
		return genomicPosition;
	}
	
	public int getMidpointGenomicCoords() {
		return transcriptToGenomicPosition(getSize() / 2);
	}
	
	/**
	 * Convert an interval in transcriptome space to genome space
	 * @param startPosOnTranscript Start position in transcriptome space
	 * @param endPosOnTranscript End position in transcriptome space
	 * @return The interval (possibly spliced) in genome space
	 */
	public Gene transcriptToGenomicPosition(int startPosOnTranscript, int endPosOnTranscript) {
		if(startPosOnTranscript > endPosOnTranscript) {
			throw new IllegalArgumentException("Start position in transcriptome space cannot be greater than end position.");
		}
		int genomicStart = isNegativeStrand() ? transcriptToGenomicPosition(endPosOnTranscript) : transcriptToGenomicPosition(startPosOnTranscript);
		int genomicEnd = isNegativeStrand() ? transcriptToGenomicPosition(startPosOnTranscript) : transcriptToGenomicPosition(endPosOnTranscript);
		Gene genomicInterval = new Gene(getChr(), genomicStart, genomicEnd);
		return getOverlap(genomicInterval);
	}
	
	/**
	 * Maps genomic coordinates to cDNA coordinates
	 * @param genomicPosition within transcript
	 * @return 0-based position within oriented start of transcript
	 */
	public int genomicToTranscriptPosition(int genomicPosition) {
		if(genomicPosition < getStart() || genomicPosition > getEnd() ) {
			return -1;
		}
		
		List<Annotation> exons = new ArrayList<Annotation>(getExonSet());
				
		int position = 0;

		for(int i = 0; i < exons.size() ; i++) {
			Annotation e = exons.get(i);

			if(genomicPosition > e.getEnd()) {
				position += e.length();
			} else if (e.getStart() <= genomicPosition && genomicPosition < e.getEnd()) {
				position +=  genomicPosition - e.getStart();
				break;
			} else {
				return -1;
			}
		}
	
		if(getOrientation().equals(Strand.NEGATIVE)) {
			position = position > -1 ? length() - 1 - position : position;
		}
		return position;
	}
	
	/**
	 * Get the genomic coordinate of the given position plus offset along the transcript
	 * @param genomicPosition Genomic position
	 * @param offset Offset (positive gives result in 3' direction; negative gives result in 5' direction)
	 * @return The genomic position at the given transcript distance from the given position
	 */
	public int genomicToGenomicPositionWithOffset(int genomicPosition, int offset) {
		int transcriptPos = genomicToTranscriptPosition(genomicPosition);
		if(transcriptPos + offset < 0) {
			throw new IllegalArgumentException(getName() + " " + genomicPosition + " is already within " + offset + " positions of 5' end of transcript");
		}
		if(transcriptPos + offset >= getSize()) {
			throw new IllegalArgumentException(getName() + " " + genomicPosition + " is already within " + offset + " positions of 3' end of transcript");
		}		
		int nextPos = transcriptPos + offset;
		return transcriptToGenomicPosition(nextPos);
	}
	
	public static int [] findLongestORF(String sequence) {
		int lastStart = 0;
		int lastEnd = 0;
		Matcher m = START_CODON_PATRN.matcher(sequence);
		while(m.find()) {
			int startCodonPos = m.start(); 
			int thisORFEnd = startCodonPos;
			boolean foundStopCodon = false;
			//System.err.println("new orf start: " + startCodonPos);
			while(thisORFEnd < sequence.length() - 3  && !foundStopCodon) {
				String currentCodon = sequence.substring(thisORFEnd, thisORFEnd+3);
				//System.err.print(" "+currentCodon+ " ");
				thisORFEnd += 3;
				foundStopCodon = isStopCodon(currentCodon);
			}
			//System.err.println("Orf length: " + (thisORFEnd - startCodonPos) );
			if(lastEnd - lastStart < thisORFEnd - startCodonPos) {
				lastStart = startCodonPos;
				lastEnd   = thisORFEnd;
				//System.err.println("It was a winner");
			}
		}
		int [] startEnd = {lastStart, lastEnd};
		return startEnd;
	}
	
	/***************************************************************************/
	
	public Gene takeUnion(Gene other) {
		TreeSet<Annotation> combinedExons = new TreeSet<Annotation>();
		combinedExons.addAll(getExonSet());
		combinedExons.addAll(other.getExonSet());
		List<Annotation> mergedExons = BasicLightweightAnnotation.stitchList(combinedExons, 0);
		Gene union = new Gene(mergedExons);
		union.setName(getName());
		union.setOrientation(getOrientation());
		return union;
		
	}

	public double getPercentCDS() {
		return (double)this.getCDS().getSize()/(double)this.getSize();
	}

	public double getPercent3UTR() {
		if(this.get3UTRGene()!=null){
			return (double)this.get3UTRGene().getSize()/(double)this.getSize();
		}
		return 0.0;
	}
	

	/**
	 * Get new gene extended in either direction
	 * @param deltaStart Distance to extend at begin position
	 * @param deltaEnd Distance to extend at end position
	 * @return New extended gene
	 */
	public Gene extendAnnotation(int deltaStart, int deltaEnd) { 
		Gene subAnnotation = new Gene(this.copy());
		subAnnotation.expand(deltaStart, deltaEnd);
		return subAnnotation;
	}

	
	private static Gene makeGene(String rawData, boolean isPSLFormat) {
		if(isPSLFormat) {
			return makeFromPSL(rawData);
			
		} else { 
			//Assume regular BED format
			return makeFromBED(rawData);
		}
	}
	
	
	private static Gene makeFromBED(String rawData) {
		StringParser s = new StringParser();
		s.parse(rawData);
		String[] tokens = s.getStringArray();
		String chr=(tokens[0]);
		int start=new Integer(tokens[1]);
		int end=new Integer(tokens[2]);
		Strand orientation=Strand.UNKNOWN;
		if(tokens.length > 3) {
			String name=tokens[3];
			if(tokens.length > 4) {
				double bedScore = new Double(tokens[4]);
				if(tokens.length > 5){
					orientation= AbstractAnnotation.getStrand(tokens[5]);
					
					if(tokens.length > 6) {
						int cdsStart=Integer.parseInt(tokens[6]);
						int cdsEnd=Integer.parseInt(tokens[7]);
						
						String[] blockSizes=tokens[10].split(",");
						String[] blockStarts=tokens[11].split(",");
						List<Integer>[] exonStartEnd=setBlockStartsAndEnds(blockStarts, blockSizes, new Integer(tokens[9]), start);
						Collection<Annotation> exons=new ArrayList<Annotation>();
						
						for(int i = 0; i < blockSizes.length; i++ ) {
							Annotation exon=new BasicAnnotation(chr, exonStartEnd[0].get(i), exonStartEnd[1].get(i), orientation, name);
							exons.add(exon);
						}
						
						String [] extraColumns = null;
						
						Gene g = new Gene(chr, name, orientation, exons, cdsStart, cdsEnd);
						assert(g.getStart() == start && g.getEnd() == end); // JE implementation check
						assert(g.getCDSStart() == cdsStart && g.getEnd() == cdsEnd);
  						
						g.setBedScore(bedScore);
						g.setScore(bedScore);
						
						if(tokens.length > 12) {
							extraColumns = new String[tokens.length - 12];
							for(int j = 12; j < tokens.length; j++) {
								extraColumns[j - 12] = tokens[j];
							}
							g.setExtraFields(extraColumns);
						}
						return g;
						
					}
					else{
						Gene g=new Gene(chr, start, end, name, orientation);
						g.setBedScore(bedScore);
						g.setScore(bedScore);
						return g;
					}
				}
				else{
					Gene g=new Gene(chr, start, end, name);
					g.setBedScore(bedScore);
					g.setScore(bedScore);
					return g;
				}
			}
			else{
				Gene g = new Gene(chr, start, end, name);
				return g;
			}
		}
		else{
			Gene g = new Gene(chr, start, end);
			return g;
		}
		
	}

	private static Gene makeFromPSL(String rawData) {
		String[] tokens=rawData.split("\t");
		
		String chr=tokens[13];
		int start=new Integer(tokens[15]);
		int end=new Integer(tokens[16]);
		String name=tokens[9];
		Strand orientation=AbstractAnnotation.getStrand(tokens[8]);
		int numExons=new Integer(tokens[17]);
		
		Collection<Annotation> exons=new ArrayList<Annotation>();
		
		int[] exonStart=new int[numExons];
		int[] exonEnd=new int[numExons];
		
		String[] exonStarts=tokens[20].replaceAll("\"", "").trim().split(",");
		String[] exonSizes=tokens[18].replaceAll("\"", "").trim().split(",");
		
		for(int i=0; i<exonStarts.length; i++){
			exonStart[i]=new Integer(exonStarts[i].toString());
			exonEnd[i]=exonStart[i]+new Integer(exonSizes[i].toString());
			exons.add(new Alignments(chr, exonStart[i], exonEnd[i], orientation, name));
		}
		
		Sequence seq=new Sequence(name);
		seq.setSequenceBases(tokens[21].replaceAll(",", "").trim().toUpperCase());
		
		Gene g=new Gene(chr, name, orientation, exons);
		assert(g.getStart() == start && g.getEnd() == end);  // JE implementation check
		
		g.addSequence(seq);
		
		
		return g;
	}
	
	public void addIsoform(Gene g){
		if(this.isoforms==null || this.isoforms.isEmpty()){
			this.isoforms=new TreeSet<Gene>();
			this.isoforms.add(new Gene(this));
		}
		this.isoforms.add(g);
	}
	
	public Collection<Gene> getIsoforms() {
		if(this.isoforms==null || this.isoforms.isEmpty()){
			this.isoforms=new TreeSet<Gene>();
			this.isoforms.add(new Gene(this));
		}
		
		return this.isoforms;
	}
	
	/*public List<BasicAnnotation> getBlocks() {
		
		List<BasicAnnotation> rtrn = new ArrayList<BasicAnnotation>();
		
		
		if (exonStarts != null){
			for(int i=0; i<exonStarts.length; i++){
				BasicAnnotation exon = new BasicAnnotation(chr, exonStarts[i], exonEnds[i],orientation);
				exon.setName(getName()+"_"+i);
				rtrn.add(exon);
			}
		} else {
			BasicAnnotation exon = new BasicAnnotation(chr, getStart(), getEnd(),orientation);
			rtrn.add(exon);
		}

		return rtrn;
	}*/
	
}
