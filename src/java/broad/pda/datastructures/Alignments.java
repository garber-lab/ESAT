package broad.pda.datastructures;

import java.io.FileWriter;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.*;

import org.broad.igv.sam.Alignment;

import broad.core.annotation.BasicGenomicAnnotation;
import broad.core.annotation.LightweightGenomicAnnotation;
import broad.core.math.Statistics;
import broad.core.sequence.Extractor;
import broad.core.sequence.Sequence;
import broad.pda.rnai.ExtractSequence;

import nextgen.core.annotation.AbstractAnnotation;
import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Gene;

//import broad.pda.rnai.ExtractSequence;


public class Alignments extends BasicGenomicAnnotation {

	String line;
	List<Double> scores; //What are these?
	private double countScore; //What are these?
	
	
	public Alignments(String ucsc){
		this(ucsc.split(":")[0], new Integer(ucsc.split(":")[1].split("-")[0]), new Integer(ucsc.split(":")[1].split("-")[1]));
	}

	public Alignments(String chr, int start, int end){
		this(chr, start, end, Strand.UNKNOWN);
	}
	
	public Alignments(String chr, int start, int end, Strand strand){
		this(chr, start, end, strand, "");
	}
	
	public Alignments(String chr, int start, int end, String strand){
		this(chr, start, end, AbstractAnnotation.getStrand(strand));
	}
	
	public Alignments(String chr, int start, int end, Strand strand, String name){
		super(name, chr, start, end, strand);
	}
	
	public Alignments(String chr, String start, String end){
		this(chr, Integer.parseInt(start), Integer.parseInt(end));
	}
	
	public Alignments(Annotation a) {
		this(a.getReferenceName(), a.getStart(), a.getEnd(), a.getOrientation(), a.getName());
	}

	
	public Alignments(String fileName, boolean fromFileName){
		this(fileName.split("_")[0], new Integer(fileName.split("_")[1].split("-")[0]), new Integer(fileName.split("_")[1].split("-")[1]));
	}
	
	public Alignments(LightweightGenomicAnnotation a) {
		this(a.getReferenceName(), a.getStart(), a.getEnd(), a.getOrientation(), a.getName());
	}

	
	public Alignments(Alignment record) {
		this(record.getChromosome(), record.getStart(), record.getEnd(), strand(record.isNegativeStrand()), record.getReadName());
	}
	
	private static Strand strand(boolean strand) {
		Strand s=Strand.POSITIVE;
		if(strand){s=Strand.NEGATIVE;}
		return s;
	}

	public Alignments(nextgen.core.alignment.Alignment record) {
		this(record.getReferenceName(), record.getStart(), record.getEnd(), record.getOrientation(), record.getName());
	}
	
	
	
	/*************Might be bad************/
	public double getCountScore(){return this.countScore;}
	
	public List<Double> getScores() {return this.scores;}

	public String getSequence(Sequence chr, boolean repeatMask) throws Exception {
		return ExtractSequence.getSequenceUnoriented(this, chr, repeatMask, new TreeMap());
	}
	
	public String getSequence(String genomeDir) throws Exception {
		return ExtractSequence.getSequenceUnoriented(this, genomeDir, false, new TreeMap());
	}

	public void setCountScore(double score){this.countScore=score;}

	public void setLine(String line){this.line=line;}

	public void setScores(List<Double> values) { this.scores = values;}

	public void addScore(double score){
		if(this.scores==null){this.scores=new ArrayList<Double>();}
		this.scores.add(score);
	}
	
	public double getScore(){
		return countScore;
	}
	
	
	/*public boolean fullyContained(Annotation g){
		if(!g.getChr().equals(getChr())){return false;}
		if(g.getStart()>=getStart() && g.getEnd()<=getEnd()){return true;}
		return false;
	}*/
	
	
	public String toString(){
		if(this.line==null || this.line.isEmpty()){return getChr()+"\t"+getStart()+"\t"+getEnd()+"\t"+getName()+"\t"+getScore()+"\t"+getOrientation();}
		else{return this.line;}
	}
	
	public String toFileName(){return getChr()+"_"+getStart()+"-"+getEnd();}
	public String toUCSC(){return this.getChr()+":"+getStart()+"-"+getEnd();}
	
	public String toStringNum(){
		return getChr().replaceAll("chr", "")+"\t"+getStart()+"\t"+getEnd();
	}	
	

	public boolean overlapsCollection(Collection<Alignments> c){
		for(Alignments align: c){
			if(overlaps(align)){return true;}
		}
		return false;
	}

	//This is closed-closed
	public boolean overlapsAtAll(Alignments align){
		if(!align.getChr().equalsIgnoreCase(getChr())){return false;}
		if(align.getStart()>=getStart() && align.getStart()<=getEnd()){return true;}
		if(getStart()>=align.getStart() && getStart()<=align.getEnd()){return true;}
		return false;
	}
	
	
	public boolean overlapsAtAll(Set<Alignments> set){
		
		for(Alignments align: set){
			if(overlapsAtAll(align)){return true;}
		}
		return false;
	}
	
	public boolean within(Alignments bigger){
		if(!bigger.getChr().equalsIgnoreCase(getChr())){return false;}
		if(getStart()>=bigger.getStart() && getStart()<=bigger.getEnd()){return true;}
		if(getEnd()>=bigger.getStart() && getEnd()<=bigger.getEnd()){return true;}
		return false;
	}
	
	public int length(){
		return getEnd()-this.getStart();
	}
	
	public boolean containedWithin(Alignments align){
		if(!align.getChr().equalsIgnoreCase(getChr())){return false;}
		
		if(align.getStart()>getStart() && align.getEnd()<getEnd()){return true;}
		return false;
	}
	
	public boolean containedWithin(Alignments interval, Alignments gene){
		if(!gene.getChr().equalsIgnoreCase(interval.getChr())){return false;}
		
		if(gene.getStart()>=interval.getStart() && gene.getEnd()<=interval.getEnd()){return true;}
		return false;
	}
	
	public static int distanceBetween(Alignments align1, Alignments align2){
		int start=Math.max(align1.getStart(), align2.getStart());
		int end=Math.min(align1.getEnd(), align2.getEnd());
		return start-end;
	}
	
	public Alignments merge(Alignments a1, Alignments a2){
		Alignments temp=new Alignments(getChr(), Math.min(a1.getStart(), a2.getStart()), Math.max(a1.getEnd(), a2.getEnd()));
		temp=new Alignments(getChr(), Math.min(temp.getStart(), getStart()), Math.max(temp.getEnd(), getEnd()));
		//System.err.println(a1.toUCSC()+" "+a2.toUCSC()+" "+this.toUCSC()+temp.toUCSC());
		
		return temp;
	}
	
	public Alignments merge(Alignments a1){
		Alignments temp=new Alignments(getChr(), Math.min(a1.getStart(), getStart()), Math.max(a1.getEnd(), getEnd()));
		
		return temp;
	}
		
	
	public static void write(String save, Collection<Alignments> alignments) throws IOException {
		write(save, alignments, true);
	}
	
	/**
	 * Write alignments 
	 * @param save
	 * @param alignments
	 * @param sort				sort by genomic position
	 * @throws IOException
	 */
	public static void write(String save, Collection<Alignments> alignments, boolean sort) throws IOException{
		BufferedWriter writer=new BufferedWriter(new FileWriter(save));
		write(writer, alignments, sort);
		writer.close();
	}
	
	
	public int getMidPoint() {
		return (getStart()+getEnd())/2;
	}
	
	/**
	 * Write alignments sorted by score i
	 * @param args
	 */
	public static void write(String save, Collection<Alignments> alignments, final int sortOnScore) throws IOException {
		List<Alignments> keys = new LinkedList<Alignments>(alignments);
		Collections.sort(keys, new Comparator<Alignments>() {
	         @Override
		        public int compare(Alignments o1, Alignments o2) {
		            return o1.getScores().get(sortOnScore).compareTo(o2.getScores().get(sortOnScore));
		         }
		 });
		write(save, alignments, false);
	}
	
	public static void write(BufferedWriter writer, Collection<Alignments> alignments) throws IOException {
		write(writer, alignments, false);
	}
	
	public static void write(BufferedWriter writer, Collection<Alignments> alignments, boolean sort) throws IOException {
		List<Alignments> keys = new LinkedList<Alignments>(alignments);
		if (sort) {
			Collections.sort(keys, new Comparator<Alignments>() {
		         @Override
		         public int compare(Alignments o1, Alignments o2) {
		             return o1.compareTo(o2);
		         }
		     });
		}
		
		for(Alignments align: keys) {
			writer.write(align.toStringWithScores());
			writer.write("\n");
		}
	}
	
	
	public String toStringWithScores() {
		StringBuffer sb = new StringBuffer(toString() + "\t" + countScore);
		for (Double score: getScores()) {
			sb.append("\t" + score);
		}
		return sb.toString();
	}
	

}
