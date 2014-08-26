package broad.core.multiplealignment;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Stack;


import broad.core.annotation.GFF;
import broad.core.annotation.GFFReader;
import broad.core.error.ParseException;
import broad.core.multiplealignment.MultipleAlignment.AlignedSequence;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;

public class AlignmentExtractor {
	public static ArrayList<String> degenerateCodonStarts;
	public static ArrayList<String> twoFoldDegenerateCodons;
	public static ArrayList<String> threeFoldDegenerateCodons;
	public static ArrayList<String> nonDegenerateCodons;
	public static HashMap<String, String> codonToAminoacid;
	
	public static String USAGE = "Usage: AlignmentExtractor TASK=<task_num> <task_args>\n" +
	"\tTasks:\n" +
	"\t\t1. Extract Regions: IN=<Multiple alignment file> OUT=<output file> ANNOTGFF=<GFF annotation file, sequence name should match one of the aligned sequence IDs> INFORMAT=<[FASTA]/EXON/PHYLIP/SEQPHYLIP> OUTFORMAT<[FASTA]/EXON/PHYLIP> SEQS=<comma separated list of sequence ids>\n" +
	"\t\t2. Extract 4D sites from exon alignment: IN=<Exon multiple alignment> OUT=<output file> INFORMAT=<input format (default FASTA)> OUTFORMAT=<output format (default FASTA)>\n " +
	"\t\t3. Extract Conserved sites from exon alignment: IN=<Exon multiple alignment> OUT=<output alignment file> INFORMAT=<input format (default FASTA)> OUTFORMAT=<output format (default FASTA)>\n"+
	"\t\t4. Extract region IN=<Multiple alignment file> OUT=<output file> REGION=<start..end> REF=<reference> INFORMAT=<[FASTA]/EXON/PHYLIP> OUTFORMAT<[FASTA]/EXON/PHYLIP/SEQPHYLIP> \n" +
	"\t\t5. Convert alignment file IN=<Multiple Alignment File> OUT=<output file> INFORMAT=<[FASTA]/EXON/PHYLIP> OUTFORMAT<FASTA/EXON/[PHYLIP]/SEQPHYLIP> [-compress <if the output alignment should be reference gap free>] \n" +
	"\t\t6. Introduce Random Gaps -in <Multiple Alignment File default is standard input> -out <output gapped alignment defalut is standard out> -informat <[FASTA]/EXON/PHYLIP> -outformat <FASTA/EXON/[PHYLIP]/SEQPHYLIP> -maxGaps <maximum number of gaps> [-ref <reference if no gaps should be introduced in reference>]\n" +
	"\t\t7. Permute Columns -in <Multiple Alignment File default is standard input> -out <output gapped alignment defalut is standard out> -informat <[FASTA]/EXON/PHYLIP> -outformat <FASTA/EXON/[PHYLIP]/SEQPHYLIP> " +
	"\n\t8. Sample Columns from Alignment -in <Multiple Alignment File default is standard input> -out <output gapped alignment defalut is standard out> -informat <[FASTA]/EXON/PHYLIP> -outformat <[FASTA]/EXON/PHYLIP/SEQPHYLIP> -cols <Number of columns to sample> [-consecutiveCols <Sample consecutive columns rather than single ones>]" +
	"\n\t9. Combined alignments -informat <[FASTA], PHYLIP, SEQPHYLIP> -outformat <[FASTA], PHYLIP, SEQPHYLIP> -in <Standard input of file with a list of alignment files> -out <standard output or file name> [-maxColumns <If specified, and the combined alignment is larger than this value, a sampled alignmnet will be generated instead of the combined one>]" +
	"\n";
	
	
	static {
		codonToAminoacid = new HashMap<String, String>();
		
		degenerateCodonStarts = new ArrayList<String>(8);
		degenerateCodonStarts.add("CT"); //Leu
		degenerateCodonStarts.add("GT"); //Val
		degenerateCodonStarts.add("TC"); //Ser
		degenerateCodonStarts.add("CC"); //Pro
		degenerateCodonStarts.add("AC"); //Thr
		degenerateCodonStarts.add("GC"); //Ala
		degenerateCodonStarts.add("CG"); //Arg
		degenerateCodonStarts.add("GG"); //Gly
		
		addDegenerateCodons("Leu", "CT");
		addDegenerateCodons("Val", "GT");
		addDegenerateCodons("Ser", "TC");
		addDegenerateCodons("Pro", "CC");
		addDegenerateCodons("Thr", "AC");
		addDegenerateCodons("Ala", "GC");
		addDegenerateCodons("Arg", "CG");
		addDegenerateCodons("Gly", "GG");
		
		twoFoldDegenerateCodons = new ArrayList<String>(18);
		twoFoldDegenerateCodons.add("TAT"); //Tyr
		codonToAminoacid.put("TAT", "Tyr");
		twoFoldDegenerateCodons.add("TAC"); //Tyr
		codonToAminoacid.put("TAC", "Tyr");
		twoFoldDegenerateCodons.add("GAT"); //Asp
		codonToAminoacid.put("GAT", "Asp");
		twoFoldDegenerateCodons.add("GAC"); //Asp
		codonToAminoacid.put("GAC", "Asp");
		twoFoldDegenerateCodons.add("GAA"); //Glu
		codonToAminoacid.put("GAA", "Glu");
		twoFoldDegenerateCodons.add("GAG"); //Glu
		codonToAminoacid.put("GAG", "Glu");
		twoFoldDegenerateCodons.add("TTT"); //Phe
		codonToAminoacid.put("TTT", "Phe");
		twoFoldDegenerateCodons.add("TTC"); //Phe
		codonToAminoacid.put("TTC", "Phe");
		twoFoldDegenerateCodons.add("TGT"); //Cys
		codonToAminoacid.put("TGT", "Cys");
		twoFoldDegenerateCodons.add("TGC"); //Cys
		codonToAminoacid.put("TGC", "Cys");
		twoFoldDegenerateCodons.add("CAT"); //His
		codonToAminoacid.put("CAT", "His");
		twoFoldDegenerateCodons.add("CAC"); //His
		codonToAminoacid.put("CAC", "His");
		twoFoldDegenerateCodons.add("CAA"); //Gln
		codonToAminoacid.put("CAA", "Gln");
		twoFoldDegenerateCodons.add("CAG"); //Gln
		codonToAminoacid.put("CAG", "Gln");
		twoFoldDegenerateCodons.add("AAT"); //Lys
		codonToAminoacid.put("AAT", "Lys");
		twoFoldDegenerateCodons.add("AAC"); //Lys
		codonToAminoacid.put("AAC", "Lys");
		twoFoldDegenerateCodons.add("AAG"); //Asn
		codonToAminoacid.put("AAG", "Asn");
		twoFoldDegenerateCodons.add("AAA"); //Asn
		codonToAminoacid.put("AAA", "Asn");
		twoFoldDegenerateCodons.add("AGT"); //Ser
		codonToAminoacid.put("AGT", "Ser");
		twoFoldDegenerateCodons.add("AGC"); //Ser
		codonToAminoacid.put("AGC", "Ser");
		
		threeFoldDegenerateCodons = new ArrayList<String>(3);
		threeFoldDegenerateCodons.add("ATT"); //Ile
		threeFoldDegenerateCodons.add("ATC"); //Ile
		threeFoldDegenerateCodons.add("ATA"); //Ile
		codonToAminoacid.put("ATT", "Ile");
		codonToAminoacid.put("ATC", "Ile");
		codonToAminoacid.put("ATA", "Ile");
		
		
		nonDegenerateCodons = new ArrayList<String>(2);
		nonDegenerateCodons.add("ATG"); //Met
		codonToAminoacid.put("ATG", "Met");
		nonDegenerateCodons.add("TGG"); //Trp
		codonToAminoacid.put("TGG", "Trp");
	}
	
	public static void main(String [] args) throws IOException, IllegalAccessException, ParseException {
		ArgumentMap argMap = CLUtil.getParameters(args, USAGE);
		
		
		if("1".equals(argMap.getTask())) {
			String alignmentFile = argMap.getInput();
			String outputFile    = argMap.getOutput();
			String gff 			 = argMap.getMandatory("ANNOTGFF");
			String inFormat      = argMap.containsKey("INFORMAT") ? argMap.get("INFORMAT") : "FASTA";
			String outFormat     = argMap.containsKey("OUTFORMAT") ? argMap.get("OUTFORMAT") : "FASTA" ;
			String [] sequences     = argMap.getMandatory("SEQS").split(",");
			extractGFFRegions(alignmentFile, outputFile, gff, inFormat, outFormat, sequences);
		} else if("2".equals(argMap.getTask() )) {
			String exonMultialignment = argMap.getInput();
			String output             = argMap.getOutput();
			String reference     = argMap.getMandatory("REF");
			String inFormat      = argMap.containsKey("INFORMAT") ? argMap.get("INFORMAT") : "FASTA";
			String outFormat     = argMap.containsKey("OUTFORMAT") ? argMap.get("OUTFORMAT") : "FASTA" ;
			String [] sequences     = argMap.getMandatory("OTHERSEQS").split(",");
			extractAFourFoldDegenerateBases(exonMultialignment, output, reference, sequences, inFormat, outFormat);
		} else if("3".equals(argMap.getTask() )) {
			String exonMultialignment = argMap.getInput();
			String output             = argMap.getOutput();
			String reference     = argMap.getMandatory("REF");
			String inFormat      = argMap.containsKey("INFORMAT") ? argMap.get("INFORMAT") : "FASTA";
			String outFormat     = argMap.containsKey("OUTFORMAT") ? argMap.get("OUTFORMAT") : "FASTA" ;
			String [] sequences     = argMap.getMandatory("OTHERSEQS").split(",");
			extractConservedBases(exonMultialignment, output, reference, sequences, inFormat, outFormat);
		} else if ("4".equals(argMap.getTask())) {
			String [] startEnd = argMap.getMandatory("REGION").split("\\.\\.");
			int start = Integer.parseInt(startEnd[0]);
			int end   = Integer.parseInt(startEnd[1]);
			String multiAlignment = argMap.getInput();
			String output        = argMap.getOutput();
			String reference     = argMap.getMandatory("REF");
			String inFormat      = argMap.containsKey("INFORMAT") ? argMap.get("INFORMAT") : "FASTA";
			String outFormat     = argMap.containsKey("OUTFORMAT") ? argMap.get("OUTFORMAT") : "FASTA" ;
			extractSubAlignment(multiAlignment, output, reference, start, end, inFormat, outFormat);
		} else if ("5".equals(argMap.getTask())) {
			String inFormat       = argMap.containsKey("INFORMAT") ? argMap.get("INFORMAT") : "FASTA";
			String outFormat      = argMap.containsKey("OUTFORMAT") ? argMap.get("OUTFORMAT") : "PHYLIP";
			InputStream is = argMap.getInputStream();
			MultipleAlignment ma = MultipleAlignmentFactory.create(is, inFormat);
			is.close();
			
			if(argMap.containsKey("compress")) {
				if(argMap.containsKey("ref")) {
					ma.setReferenceId(argMap.get("ref"));
				}
				ma.compress();
			}
			
			BufferedWriter bw = argMap.getOutputWriter();
			ma.setIOHelper(MultipleAlignmentIOFactory.create(outFormat));
			ma.write(bw);
			bw.close();
		} else if ("6".equals(argMap.getTask())) {
			String inFormat       = argMap.containsKey("informat") ? argMap.get("informat") : "FASTA";
			String outFormat      = argMap.containsKey("outformat") ? argMap.get("outformat") : "PHYLIP";
			String reference      = argMap.get("ref");
			int maxGaps           = argMap.getInteger("maxGaps");
			
			InputStream is = argMap.getInputStream();
			MultipleAlignment ma = MultipleAlignmentFactory.create(is, inFormat);
			is.close();
			
			ma.setIOHelper(MultipleAlignmentIOFactory.create(outFormat));
			if(reference != null) {
				ma.setReferenceId(reference);
			}
			ma.introduceGaps(reference == null, maxGaps);
			
			BufferedWriter bw = argMap.getOutputWriter();
			ma.write(bw);
			bw.close();
		} else if("7".equals(argMap.getTask())) {
			String inFormat       = argMap.containsKey("informat") ? argMap.get("informat") : "FASTA";
			String outFormat      = argMap.containsKey("outformat") ? argMap.get("outformat") : "FASTA";
			InputStream is = argMap.getInputStream();
			MultipleAlignment ma = MultipleAlignmentFactory.create(is, inFormat);
			is.close();
			
			ma.setIOHelper(MultipleAlignmentIOFactory.create(outFormat));
			ma.permuteColumns();
			
			BufferedWriter bw = argMap.getOutputWriter();
			ma.write(bw);
			bw.close();
		}  else if("8".equals(argMap.getTask())) {
			String inFormat       = argMap.containsKey("informat") ? argMap.get("informat") : "FASTA";
			String outFormat      = argMap.containsKey("outformat") ? argMap.get("outformat") : "FASTA";
			int cols 			  = argMap.getInteger("cols");
			int numOfConsecutiveCols = argMap.containsKey("consecutiveCols") ? argMap.getInteger("consecutiveCols") : 1;
			InputStream is = argMap.getInputStream();
			MultipleAlignment ma = MultipleAlignmentFactory.create(is, inFormat);
			is.close();
			
			MultipleAlignment sampledMA = ma.sampleColumns(cols, numOfConsecutiveCols);
			sampledMA.setIOHelper(MultipleAlignmentIOFactory.create(outFormat));
			
			BufferedWriter bw = argMap.getOutputWriter();
			sampledMA.write(bw);
			bw.close();
		} else if("9".equals(argMap.getTask())) {
			String inFormat       = argMap.containsKey("informat") ? argMap.get("informat") : "FASTA";
			String outFormat      = argMap.containsKey("outformat") ? argMap.get("outformat") : "FASTA";
			int numOfConsecutiveCols = argMap.containsKey("consecutiveCols") ? argMap.getInteger("consecutiveCols") : 1;
			
			MultipleAlignment combinedAln = new MultipleAlignment();
			BufferedReader br = argMap.getInputReader();
			String line = null;
			while((line = br.readLine()) != null) {
				line.trim();
				//System.out.println(line);
				if(!line.startsWith("#") && line.length() > 0) {
					//System.out.println("\tloading...");
					MultipleAlignment ma = MultipleAlignmentFactory.create(line, inFormat);
					if(!ma.isEmpty()) {
						combinedAln.append(ma);
					}
				}
			}
			br.close();
			
			if(argMap.containsKey("maxColumns") ) {
				int maxColumns = argMap.getInteger("maxColumns");
				if(combinedAln.length() > maxColumns) {
					MultipleAlignment sampled = combinedAln.sampleColumns(maxColumns, numOfConsecutiveCols);
					combinedAln = sampled;
				}
			}
			BufferedWriter bw = argMap.getOutputWriter();
			combinedAln.setIOHelper(MultipleAlignmentIOFactory.create(outFormat));
			combinedAln.write(bw);
			bw.close();
		} else {
			System.err.println(USAGE);
		} 
	}

	private static void addDegenerateCodons(String aminoacid, String firstTwoPositions) {
		codonToAminoacid.put(firstTwoPositions + "A", aminoacid);
		codonToAminoacid.put(firstTwoPositions + "C", aminoacid);
		codonToAminoacid.put(firstTwoPositions + "G", aminoacid);
		codonToAminoacid.put(firstTwoPositions + "T", aminoacid);
	}

	private static void extractSubAlignment(String multiAlignment, String output, String reference, int start, int end, String inFormat, String outFormat) throws IOException, ParseException {
		MultipleAlignment source = MultipleAlignmentFactory.create(multiAlignment, inFormat);
		MultipleAlignment subAln = source.getSubAlignment(reference, start, end, false);
		subAln.setIOHelper(MultipleAlignmentIOFactory.create(outFormat));
		BufferedWriter bw = new BufferedWriter(new FileWriter(output));
		subAln.write(bw);
		bw.close();
	}

	private static void extractAFourFoldDegenerateBases(String exonMultialignment, String output, String reference, String [] seqIds, String inFormat, String outFormat) 
	throws IOException, ParseException {
		MultipleAlignment ma = MultipleAlignmentFactory.create(exonMultialignment, inFormat);
		MultipleAlignment fourDAln =  MultipleAlignmentFactory.create(outFormat);
		fourDAln.setReferenceId(reference);
		//System.out.println("Aligned Sequences: " + ma.getAlignedSequenceIds() );
		AlignedSequence ref4D = new AlignedSequence(reference);
		ref4D.setName(reference);
		fourDAln.addSequence(ref4D);
		String ref = ma.getAlignedSequence(reference).getSequenceBases();
		String [] nonRefSeqs = new String[seqIds.length];
		AlignedSequence [] nonRef4Ds  = new AlignedSequence[seqIds.length];
		for(int i = 0; i < seqIds.length; i++) {
			//System.out.println("Processing sequence " + seqIds[i] + " has alignment? " + (ma.getAlignedSequence(seqIds[i]) != null));
			nonRefSeqs[i] = ma.getAlignedSequence(seqIds[i]).getSequenceBases();
			nonRef4Ds[i]  = new AlignedSequence(seqIds[i]);
			nonRef4Ds[i].setName(seqIds[i]);
			fourDAln.addSequence(nonRef4Ds[i]);
		}
		
		
		int codon = 0;
		int pos = 0;
		while(pos < ref.length()-3) {			
			String refCodon = ref.substring(pos, pos + 3);
			//System.out.print("Ref codon " + refCodon);
			if(degenerateCodonStarts.contains(refCodon.substring(0, 2))) {
				//System.out.print(" is 4D ");
				ref4D.append(refCodon.charAt(2));
				for(int i = 0; i < nonRefSeqs.length; i++) {
					nonRef4Ds[i].append(nonRefSeqs[i].charAt(pos+2));
				}
			}
			//System.out.println("");
			codon++;
			pos = codon * 3;
		}
		BufferedWriter bw = new BufferedWriter(new FileWriter(output));
		fourDAln.write(bw);
		bw.close();
	}
		
	private static void extractConservedBases(String exonMultialignment, String output, String reference, String [] seqIds, String inFormat, String outFormat) 
		throws IOException, ParseException {
			MultipleAlignment ma = MultipleAlignmentFactory.create(exonMultialignment, inFormat);
			MultipleAlignment consAln =  MultipleAlignmentFactory.create(outFormat);
			consAln.setReferenceId(reference);
			AlignedSequence refCons = new AlignedSequence(reference);
			refCons.setName(reference);
			consAln.addSequence(refCons);
			String ref = ma.getAlignedSequence(reference).getSequenceBases();
			String [] nonRefSeqs = new String[seqIds.length];
			AlignedSequence [] nonRefCons  = new AlignedSequence[seqIds.length];
			for(int i = 0; i < seqIds.length; i++) {
				nonRefSeqs[i] = ma.getAlignedSequence(seqIds[i]).getSequenceBases();
				nonRefCons[i]  = new AlignedSequence(seqIds[i]);
				nonRefCons[i].setName(seqIds[i]);
				consAln.addSequence(nonRefCons[i]);
			}
			
			
			int codon = 0;
			int pos = 0;
			while(pos < ref.length() - 3) {			
				String refCodon = ref.substring(pos, pos + 3);
				//System.out.print("Ref codon " + refCodon);

					//System.out.print(" is Conserved ");
				refCons.append(refCodon.substring(0,2));
				for(int i = 0; i < nonRefSeqs.length; i++) {
					nonRefCons[i].append(nonRefSeqs[i].substring(pos, pos+2));
				}
				//System.out.println("");
				codon++;
				pos = codon * 3;
			}
		
			BufferedWriter bw = new BufferedWriter(new FileWriter(output));
			consAln.write(bw);
			bw.close();
	}

	public static void extractGFFRegions(String alignmentFile, String outputFile, String gff, String inFormat, String outFormat, String [] seqs) throws IOException, ParseException {
		Stack<GeneAlignmentAnalysis> problemGenes = new Stack<GeneAlignmentAnalysis>();
		MultipleAlignment ma = MultipleAlignmentFactory.create(alignmentFile, inFormat);
		if(!inFormat.equals(outFormat)) {
			MultipleAlignmentIO maio = MultipleAlignmentIOFactory.create(outFormat);
			ma.setIOHelper(maio);
		}
		GFFReader gffR = new GFFReader(gff);
		Map<String, List<GFF>> seqGFFMap = gffR.getSequenceAnnotationMap();
		if (seqGFFMap.keySet().size() != 1) {
			throw new IllegalStateException("The loaded GFF file contains annotation for more than one sequence: " + seqGFFMap.keySet());
		}
		
		String refId = seqGFFMap.keySet().iterator().next();
		Iterator<GFF> annotationIt = seqGFFMap.get(refId).iterator();
		if (!annotationIt.hasNext() ) {
			throw new RuntimeException ("The GFF list loaded from " + gff + " was empty");
		}
		
		GFF annotation = annotationIt.next();
		MultipleAlignment subAlignment = ma.getSubAlignment(refId, annotation.getStart(), annotation.getEnd(), annotation.inReversedOrientation());
		GeneAlignmentAnalysis firstAnalysis = new GeneAlignmentAnalysis(annotation.getFirstValue("refSeqId"));

		firstAnalysis.addExonAlignment(ma, Integer.valueOf(annotation.getFirstValue("ExonNum")));
		problemGenes.push(firstAnalysis);
		
		MultipleAlignment tmpAlignment = null;
		while (annotationIt.hasNext() ) {
			annotation = annotationIt.next();
			String geneId = annotation.getFirstValue("refSeqId");
			int exonNumber =  Integer.parseInt(annotation.getFirstValue("ExonNum"));
			//GeneAlignmentAnalysis currentAnalysis = problemGenes.peek();
			tmpAlignment = ma.getSubAlignment(refId, annotation.getStart(), annotation.getEnd(), annotation.inReversedOrientation());
			/*
			if (geneId.equals(currentAnalysis.getGeneId())) {
				currentAnalysis.addExonAlignment(ma,exonNumber);
			}else {
				currentAnalysis = new GeneAlignmentAnalysis(geneId);
				currentAnalysis.addExonAlignment(ma, exonNumber);
				problemGenes.add(currentAnalysis);
			}
			*/
			subAlignment.append(tmpAlignment);
		}
		
		BufferedWriter bw = new BufferedWriter(new FileWriter(outputFile));
		subAlignment.write(bw);
		bw.close();
		
		/*
		bw = new BufferedWriter(new FileWriter("exonAlignmentStats.txt", true));
		for(int i = 0; i < problemGenes.size(); i ++) {
			problemGenes.get(i).write(bw, seqs);
			bw.newLine();
		}
		bw.close();
		*/
	}

	private static class GeneAlignmentAnalysis {
		private String geneId;
		private HashMap<String, List<Integer>> speciesGappedExons;
		private HashMap<String, Integer> speciesFirstMissenseExons;
		
		public GeneAlignmentAnalysis(String geneId) {
			this.geneId = geneId;
			speciesGappedExons = new HashMap<String, List<Integer>>();
			speciesFirstMissenseExons = new HashMap<String, Integer>();
		}
		
		public boolean hasIssues() { return !(speciesGappedExons.isEmpty() && speciesFirstMissenseExons.isEmpty()); } 
		
		public boolean hasMissenseExons() {return !speciesFirstMissenseExons.isEmpty();}
		
		public String getGeneId() { return geneId;}
		
		public void addExonAlignment(MultipleAlignment ma, int exonNumber) {
			Iterator<String> it = ma.getAlignedSequenceIds().iterator();
			while(it.hasNext()) {
				String id = it.next();
				AlignedSequence seq = ma.getAlignedSequence(id);
				float percentGaps = seq.getPercentGaps();
				if(percentGaps > 0.5) {
					List<Integer> gappedExons = speciesGappedExons.get(id);
					if(gappedExons == null) {
						gappedExons = new ArrayList<Integer>();
					}
					gappedExons.add(exonNumber);
					continue;
				}
				List<Integer> gapSizes = seq.getGapSizes();
				if (!speciesFirstMissenseExons.containsKey(id) && !gapSizes.isEmpty() && (gapSizes.size() > 1 || !seq.getSequenceBases().endsWith("-"))) {
					 int totalGaps = 0;
					 for(int i = 0; i < gapSizes.size(); i++) {
						 totalGaps =+ gapSizes.get(i);
					 }
					 if (totalGaps % 3 > 0) {
						 speciesFirstMissenseExons.put(id, exonNumber);
					 }
				}
				
			}
		}
		
		public void write(BufferedWriter bw, String [] seqIds) throws IOException {
			bw.write(geneId);
			bw.write("\t");
			for(int j = 0; j < seqIds.length; j++) {
				String id = seqIds[j];
				if(speciesGappedExons.containsKey(id)) {
					List<Integer> unalignedExons = speciesGappedExons.get(id);
					for(int i = 0; i < unalignedExons.size(); i ++) {
						bw.write(String.valueOf(unalignedExons.get(i)));
						if(i < unalignedExons.size() - 1) {
							bw.write(",");
						}
					}
				} else {
					bw.write("0");
				}
				bw.write("\t");
				
				if(speciesFirstMissenseExons.containsKey(id)) {
					if (speciesFirstMissenseExons.containsKey(id)) {
						bw.write(String.valueOf(speciesFirstMissenseExons.get(id)));
					} else {
						bw.write(0);
					}
				}
			}
		}
	}
	
	
}
