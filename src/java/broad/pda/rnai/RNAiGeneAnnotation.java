package broad.pda.rnai;

import nextgen.core.annotation.Gene;
import broad.core.sequence.Sequence;
import broad.pda.datastructures.Alignments;

public class RNAiGeneAnnotation {

	public static String nanostringHeader="GENE.SOURCEID\tGENE.TAXONID\tGENE.NCBI_BUILDID\tGENE.SOURCECHR\tGENE.SOURCESTART\tGENE.SOURCEEND\tGENE.SOURCESTRAND\tTRANS.SOURCEID\tTRANS.SEQ\tTRANS.LENGTH";
	String geneSourceID;
	String geneSourceVersion;
	String geneSymbol;
	String geneEntrezGeneId;
	String geneTaxonId;
	String geneChr;
	String geneMapLocation;
	String geneNcbiBuild;
	String geneSourceContig;
	String geneSourceStart;
	String geneSourceEnd;
	String geneSourceStrand;
	String geneCreateDate;
	String geneUpdateDate;
	String transcriptSourceId;
	String transcriptSourceVersion;
	String transcriptSeq;
	String transcriptLength;
	String transcriptSeqCreateDate;
	String transcriptRefSeq;
	String transcriptGenBankId;
	Gene gene;
	
	public RNAiGeneAnnotation(Gene gene, Alignments linc, String geneName, int versionNumber, String seq, String species, String ncbiBuild){
		this.gene=gene;
		geneSourceID=geneName;
		geneSourceVersion="1";
		geneSymbol=geneName;
		geneEntrezGeneId="NA";
		geneTaxonId=species;
		geneChr=gene.getChr();
		geneMapLocation="NA";
		geneNcbiBuild=ncbiBuild;
		geneSourceContig=linc.getChr();
		geneSourceStart=new Integer(linc.getStart()).toString();
		geneSourceEnd=new Integer(linc.getEnd()).toString();
		geneSourceStrand=gene.getOrientation().toString();
		geneCreateDate="NA";
		geneUpdateDate="NA";
		transcriptSourceId=geneName+"_"+versionNumber;
		transcriptSourceVersion="1";
		transcriptSeq=seq;
		transcriptLength=new Integer(seq.toCharArray().length).toString();
		transcriptSeqCreateDate="NA";
		transcriptRefSeq="NA";
		transcriptGenBankId="NA";
	}
	
	public RNAiGeneAnnotation(String rnaiString){
		String[] tokens=rnaiString.split("\t");
		geneSourceID=tokens[0];
		geneSourceVersion=tokens[1];
		geneSymbol=tokens[2];
		geneEntrezGeneId=tokens[3];
		geneTaxonId=tokens[4];
		geneChr=tokens[5];
		geneMapLocation=tokens[6];
		geneNcbiBuild=tokens[7];
		geneSourceContig=tokens[8];
		geneSourceStart=tokens[9];
		geneSourceEnd=tokens[10];
		geneSourceStrand=tokens[11];
		geneCreateDate=tokens[12];
		geneUpdateDate=tokens[13];
		transcriptSourceId=tokens[14];
		transcriptSourceVersion=tokens[15];
		transcriptSeq=tokens[16];
		transcriptLength=tokens[17];
		transcriptSeqCreateDate=tokens[18];
		transcriptRefSeq=tokens[19];
		transcriptGenBankId=tokens[20];
	}
	
	public String toString(){
		return geneSourceID+"\t"+geneSourceVersion+"\t"+geneSymbol+"\t"+geneEntrezGeneId+"\t"+geneTaxonId+"\t"+geneChr+"\t"+geneMapLocation+"\t"+geneNcbiBuild+"\t"+geneSourceContig+"\t"+geneSourceStart+"\t"+geneSourceEnd+"\t"+geneSourceStrand+"\t"+geneCreateDate+"\t"+geneUpdateDate+"\t"+transcriptSourceId+"\t"+transcriptSourceVersion+"\t"+transcriptSeq+"\t"+transcriptLength+"\t"+transcriptSeqCreateDate+"\t"+transcriptRefSeq+"\t"+transcriptGenBankId;
	}
	
	public String toRNAi(){return toString();}

	public String getTranscriptSequence() {
		return this.transcriptSeq;
	}

	public String getTranscriptName() {
		return this.transcriptSourceId;
	}

	public Alignments getRegion() {
		Alignments rtrn= new Alignments(geneChr, geneSourceStart, geneSourceEnd);
		rtrn.setOrientation(geneSourceStrand);
		return rtrn;
	}

	public String getGeneName() {
		return this.geneSourceID;
	}

	public Sequence getSequence() {
		Sequence rtrn=new Sequence(this.transcriptSourceId);
		rtrn.setSequenceBases(this.transcriptSeq);
		return rtrn;
	}

	public String toNanostring() {
		return geneSourceID+"\t"+geneTaxonId+"\t"+geneNcbiBuild+"\t"+geneSourceContig+"\t"+geneSourceStart+"\t"+geneSourceEnd+"\t"+geneSourceStrand+"\t"+transcriptSourceId+"\t"+transcriptSeq+"\t"+transcriptLength;
	}
	
}
