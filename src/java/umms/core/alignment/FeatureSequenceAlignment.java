package umms.core.alignment;

import broad.core.parser.CommandLineParser;
import jaligner.matrix.Matrix;
import jaligner.matrix.MatrixGenerator;
import jaligner.SmithWatermanGotoh;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.log4j.Logger;

import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.pda.annotation.BEDFileParser;

import umms.core.annotation.Annotation;
import umms.core.annotation.BasicAnnotation;
import umms.core.annotation.Gene;

/**
 * @author prussell
 *
 */
public class FeatureSequenceAlignment {
	
	private Map<String, broad.core.sequence.Sequence> chromosomes;
	private Map<Gene, jaligner.Sequence> sequences;
	private TreeSet<UnorderedGenePair> featurePairs;
	private Matrix scoringMatrix;
	private float gapOpen;
	private float gapExtend;
	private static Logger logger = Logger.getLogger(FeatureSequenceAlignment.class.getName());
	private TreeMap<UnorderedGenePair, jaligner.Alignment> senseAlignments;
	private TreeMap<UnorderedGenePair, jaligner.Alignment> antisenseAlignments;
	
	/**
	 * Default match score for Smith Waterman
	 */
	public static float DEFAULT_MATCH_SCORE = 5;
	
	/**
	 * Default mismatch score for Smith Waterman
	 */
	public static float DEFAULT_MISMATCH_SCORE = -4;
	
	/**
	 * Default gap open penalty for Smith Waterman
	 */
	public static float DEFAULT_GAP_OPEN_PENALTY = 8;
	
	/**
	 * Default gap extend penalty for Smith Waterman
	 */
	public static float DEFAULT_GAP_EXTEND_PENALTY = 2;
	
	/**
	 * @param featureBedFile Bed file of features
	 * @param genomeFasta Genome fasta file
	 * @throws IOException
	 */
	public FeatureSequenceAlignment(String featureBedFile, String genomeFasta) throws IOException {
		this(featureBedFile, genomeFasta, DEFAULT_MATCH_SCORE, DEFAULT_MISMATCH_SCORE, DEFAULT_GAP_OPEN_PENALTY, DEFAULT_GAP_EXTEND_PENALTY);
	}
	
	/**
	 * @param genes Set of features
	 * @param genomeFasta Genome fasta file
	 * @throws IOException
	 */
	public FeatureSequenceAlignment(Collection<Gene> genes, String genomeFasta) throws IOException {
		this(genes, genomeFasta, DEFAULT_MATCH_SCORE, DEFAULT_MISMATCH_SCORE, DEFAULT_GAP_OPEN_PENALTY, DEFAULT_GAP_EXTEND_PENALTY);
	}
	
	/**
	 * @param featureBedFile Bed file of features
	 * @param genomeFasta Genome fasta file
	 * @param matchScore Match score for Smith Waterman
	 * @param mismatchScore Mismatch score for Smith Waterman
	 * @param gapOpenPenalty Gap open penalty for Smith Waterman
	 * @param gapExtendPenalty Gap extend penalty for Smith Waterman
	 * @throws IOException
	 */
	public FeatureSequenceAlignment(String featureBedFile, String genomeFasta, float matchScore, float mismatchScore, float gapOpenPenalty, float gapExtendPenalty) throws IOException {
		this(BEDFileParser.loadData(new File(featureBedFile)), genomeFasta, matchScore, mismatchScore, gapOpenPenalty, gapExtendPenalty);
	}

	/**
	 * @param genes The features
	 * @param chrs Map of chromosome names to sequences
	 * @param matchScore Match score for Smith Waterman
	 * @param mismatchScore Mismatch score for Smith Waterman
	 * @param gapOpenPenalty Gap open penalty for Smith Waterman
	 * @param gapExtendPenalty Gap extend penalty for Smith Waterman
	 * @throws IOException
	 */
	public FeatureSequenceAlignment(Collection<Gene> genes, Map<String, broad.core.sequence.Sequence> chrs, float matchScore, float mismatchScore, float gapOpenPenalty, float gapExtendPenalty) throws IOException {
		
		chromosomes = new TreeMap<String, broad.core.sequence.Sequence>();
		chromosomes.putAll(chrs);
		
		// Store genes and sequences
		sequences = new TreeMap<Gene, jaligner.Sequence>();
		for(Gene gene : genes) {
			broad.core.sequence.Sequence broadSeq = chromosomes.get(gene.getChr()).getSubsequence(gene);
			jaligner.Sequence jSeq = new jaligner.Sequence(broadSeq.getSequenceBases());
			jSeq.setId(gene.getName());
			sequences.put(gene, jSeq);
		}
		logger.info("Loaded " + sequences.size() + " features with sequences.");
		
		// Store unordered pairs of genes
		featurePairs = new TreeSet<UnorderedGenePair>();
		for(Gene gene1 : genes) {
			for(Gene gene2 : genes) {
				if(gene1.equals(gene2)) continue;
				featurePairs.add(new UnorderedGenePair(gene1, gene2));
			}
		}
		logger.info("There are " + featurePairs.size() + " unordered pairs of features.");
		
		gapOpen = gapOpenPenalty;
		gapExtend = gapExtendPenalty;
		scoringMatrix = MatrixGenerator.generate(matchScore, mismatchScore);
		senseAlignments = new TreeMap<UnorderedGenePair, jaligner.Alignment>();
		antisenseAlignments = new TreeMap<UnorderedGenePair, jaligner.Alignment>();
		
		logger.info("Match score: " + matchScore);
		logger.info("Mismatch score: " + mismatchScore);
		logger.info("Gap open penalty: " + gapOpen);
		logger.info("Gap extend penalty: " + gapExtend);
	
	}
	
	/**
	 * @param genes Set of features
	 * @param genomeFasta Genome fasta file
	 * @param matchScore Match score for Smith Waterman
	 * @param mismatchScore Mismatch score for Smith Waterman
	 * @param gapOpenPenalty Gap open penalty for Smith Waterman
	 * @param gapExtendPenalty Gap extend penalty for Smith Waterman
	 * @throws IOException
	 */
	public FeatureSequenceAlignment(Collection<Gene> genes, String genomeFasta, float matchScore, float mismatchScore, float gapOpenPenalty, float gapExtendPenalty) throws IOException {
		this(genes, FastaSequenceIO.getChrSequencesFromFasta(genomeFasta), matchScore, mismatchScore, gapOpenPenalty, gapExtendPenalty);
	}
	
	
	/**
	 * Get the alignment between two genes
	 * @param genes Gene pair
	 * @return The alignment of transcript sequences
	 */
	private jaligner.Alignment getSenseAlignment(UnorderedGenePair genes) {
		if(!senseAlignments.containsKey(genes)) {
			alignSense(genes);
		}
		return senseAlignments.get(genes);
	}
	
	/**
	 * Get the antisense alignment between two genes
	 * @param genes Gene pair
	 * @return The alignment of transcript sequences
	 */
	private jaligner.Alignment getAntisenseAlignment(UnorderedGenePair genes) {
		if(!antisenseAlignments.containsKey(genes)) {
			alignAntisense(genes);
		}
		return antisenseAlignments.get(genes);
	}
	
	/**
	 * Align the sequences of two genes in sense direction; store alignment
	 * @param gene1 Gene 1
	 * @param gene2 Gene 2
	 */
	private void alignSense(UnorderedGenePair genes) {
		Gene gene1 = genes.getFirstGene();
		Gene gene2 = genes.getSecondGene();
		jaligner.Sequence seq1 = sequences.get(gene1);
		jaligner.Sequence seq2 = sequences.get(gene2);
		jaligner.Alignment senseAlign = SmithWatermanGotoh.align(seq1, seq2, scoringMatrix, gapOpen, gapExtend);
		senseAlignments.put(genes, senseAlign);
	}
		
	/**
	 * Align the sequences of two genes in antisense direction; store alignment
	 * @param gene1 Gene 1
	 * @param gene2 Gene 2
	 */
	private void alignAntisense(UnorderedGenePair genes) {
		Gene gene1 = genes.getFirstGene();
		Gene gene2 = genes.getSecondGene();
		jaligner.Sequence seq1 = sequences.get(gene1);
		jaligner.Sequence seq2 = sequences.get(gene2);
		String seq2bases = seq2.getSequence();
		String seq2antisenseBases = Sequence.reverseSequence(seq2bases);
		jaligner.Sequence seq2antisense = new jaligner.Sequence(seq2antisenseBases);
		seq2antisense.setId(seq2.getId() + "_antisense");
		jaligner.Alignment antisenseAlign = SmithWatermanGotoh.align(seq1, seq2antisense, scoringMatrix, gapOpen, gapExtend);
		int correctedStart2 = seq2antisense.length() - antisenseAlign.getStart2() - antisenseAlign.getLength() + antisenseAlign.getGaps2();
		antisenseAlign.setStart2(correctedStart2);
		antisenseAlignments.put(genes, antisenseAlign);
	}

	/**
	 * String of information about the alignment
	 * @param align The alignment
	 * @return Formatted string for printing
	 */
	private static String formattedInfoString(jaligner.Alignment align) {
		String seq1 = new String(align.getSequence1());
		String seq2 = new String(align.getSequence2());
		String markup = new String(align.getMarkupLine());
		String rtrn = align.getSummary();
		rtrn += seq1 + "\n";
		rtrn += markup + "\n";
		rtrn += seq2 + "\n";
		return rtrn;
	}
	
	
	
	/**
	 * Get all sense direction alignments passing thresholds
	 * @param minAlignLength Min alignment length
	 * @param minPctIdentity Min percent identity
	 * @return All sense alignments of unordered gene pair passing criteria
	 */
	public Map<UnorderedGenePair, jaligner.Alignment> getAllPairwiseSenseAlignments(float minAlignLength, float minPctIdentity) {
		Map<UnorderedGenePair, jaligner.Alignment> rtrn = new TreeMap<UnorderedGenePair, jaligner.Alignment>();
		for(UnorderedGenePair genes : featurePairs) {
			jaligner.Alignment align = getSenseAlignment(genes);
			if(align.getLength() >= minAlignLength && align.getPercentIdentity() >= minPctIdentity) {
				rtrn.put(genes, align);
			}
		}
		return rtrn;
	}
	
	/**
	 * Get all antisense direction alignments passing thresholds
	 * @param minAlignLength Min alignment length
	 * @param minPctIdentity Min percent identity
	 * @return All antisense alignments of unordered gene pair passing criteria
	 */
	public Map<UnorderedGenePair, jaligner.Alignment> getAllPairwiseAntisenseAlignments(float minAlignLength, float minPctIdentity) {
		Map<UnorderedGenePair, jaligner.Alignment> rtrn = new TreeMap<UnorderedGenePair, jaligner.Alignment>();
		for(UnorderedGenePair genes : featurePairs) {
			jaligner.Alignment align = getAntisenseAlignment(genes);
			if(align.getLength() >= minAlignLength && align.getPercentIdentity() >= minPctIdentity) {
				rtrn.put(genes, align);
			}
		}
		return rtrn;
	}
	
	/**
	 * Write all pairwise alignments (sense direction only) to bed file in genome coordinates
	 * @param outBedFile Output bed file
	 * @param minAlignLength Min alignment length to keep
	 * @param minPctIdentity Min percent identity to keep
	 * @throws IOException
	 */
	public void writeAllSenseAlignmentsToBed(String outBedFile, int minAlignLength, float minPctIdentity) throws IOException {
		writeAllSenseAlignmentsToBed(outBedFile, false, minAlignLength, minPctIdentity);
	}
	
	/**
	 * Write all pairwise alignments (antisense direction only) to bed file in genome coordinates
	 * @param outBedFile Output bed file
	 * @param minAlignLength Min alignment length to keep
	 * @param minPctIdentity Min percent identity to keep
	 * @throws IOException
	 */
	public void writeAllAntisenseAlignmentsToBed(String outBedFile, int minAlignLength, float minPctIdentity) throws IOException {
		writeAllAntisenseAlignmentsToBed(outBedFile, false, minAlignLength, minPctIdentity);
	}

	/**
	 * Write all pairwise alignments (sense direction only) to bed file in genome coordinates
	 * @param outBedFile Output bed file
	 * @param append Write to end of file rather than beginning
	 * @param minAlignLength Min alignment length to keep
	 * @param minPctIdentity Min percent identity to keep
	 * @throws IOException
	 */
	public void writeAllSenseAlignmentsToBed(String outBedFile, boolean append, int minAlignLength, float minPctIdentity) throws IOException {
		logger.info("Writing all sense alignments to bed file " + outBedFile);
		logger.info("Min alignment length: " + minAlignLength);
		logger.info("Min percent identity: " + minPctIdentity);
		Map<UnorderedGenePair, jaligner.Alignment> senseAligns = getAllPairwiseSenseAlignments(minAlignLength, minPctIdentity);
		FileWriter w = new FileWriter(outBedFile, append);
		for(UnorderedGenePair genes : senseAligns.keySet()) {
			jaligner.Alignment align = senseAligns.get(genes);
			try {
				String print = senseAlignmentAsBed(genes);
				w.write(print + "\n");
			} catch (Exception e) {
				logger.error("Caught exception when getting alignment for regions " + genes.getFirstGene().getName() + " and " + genes.getSecondGene().getName());
				logger.error("Alignment:");
				logger.error(formattedInfoString(align));
				e.printStackTrace();
			}				
		}
		w.close();
		logger.info("Done writing sense alignments.");
	}

	/**
	 * Write all pairwise alignments (antisense direction only) to bed file in genome coordinates
	 * @param outBedFile Output bed file
	 * @param append Write to end of file rather than beginning
	 * @param minAlignLength Min alignment length to keep
	 * @param minPctIdentity Min percent identity to keep
	 * @throws IOException
	 */
	public void writeAllAntisenseAlignmentsToBed(String outBedFile, boolean append, int minAlignLength, float minPctIdentity) throws IOException {
		logger.info("Writing all antisense alignments to bed file " + outBedFile);
		logger.info("Min alignment length: " + minAlignLength);
		logger.info("Min percent identity: " + minPctIdentity);
		Map<UnorderedGenePair, jaligner.Alignment> antisenseAligns = getAllPairwiseAntisenseAlignments(minAlignLength, minPctIdentity);
		FileWriter w = new FileWriter(outBedFile, append);
		for(UnorderedGenePair genes : antisenseAligns.keySet()) {
			jaligner.Alignment align = antisenseAligns.get(genes);
			try {
				String print = antisenseAlignmentAsBed(genes);
				w.write(print + "\n");
			} catch (Exception e) {
				logger.error("Caught exception when getting alignment for regions " + genes.getFirstGene().getName() + " and " + genes.getSecondGene().getName());
				logger.error("Alignment:");
				logger.error(formattedInfoString(align));
				e.printStackTrace();
			}				
		}
		w.close();
		logger.info("Done writing antisense alignments.");
	}


	
	/**
	 * Get the aligned region as a feature in genomic coordinates
	 * @param genes The pair of genes
	 * @return The aligned region
	 */
	private Annotation senseAlignmentAsFeature(UnorderedGenePair genes) {
		return asFeature(getSenseAlignment(genes), genes);
	}
	
	/**
	 * Get the aligned region of the antisense alignment as a feature in genomic coordinates
	 * @param genes The pair of genes
	 * @return The aligned region
	 */
	private Annotation antisenseAlignmentAsFeature(UnorderedGenePair genes) {
		return asFeature(getAntisenseAlignment(genes), genes);
	}

	/**
	 * The alignment as a feature in genomic coordinates
	 * @param align The alignment
	 * @param genes The pair of genes that were aligned
	 * @return The aligned region
	 */
	private static Annotation asFeature(jaligner.Alignment align, UnorderedGenePair genes) {
		Gene gene1 = genes.getFirstGene();
		int gene1start = align.getStart1();
		int gene1end = gene1start + align.getLength() - align.getGaps1() - 1;
		Gene alignedRegion1 = gene1.transcriptToGenomicPosition(gene1start, gene1end);
		int newStart1 = alignedRegion1.isNegativeStrand() ? alignedRegion1.getStart() + 1 : alignedRegion1.getStart() - 1;
		alignedRegion1.setStart(newStart1);
		Gene gene2 = genes.getSecondGene();
		int gene2start = align.getStart2();
		int gene2end = gene2start + align.getLength() - align.getGaps2() - 1;
		Gene alignedRegion2 = gene2.transcriptToGenomicPosition(gene2start, gene2end);
		int newStart2 = alignedRegion2.isNegativeStrand() ? alignedRegion2.getStart() + 1 : alignedRegion2.getStart() - 1;
		alignedRegion2.setStart(newStart2);
		BasicAnnotation rtrn = new BasicAnnotation(alignedRegion1);
		rtrn.addBlocks(alignedRegion2);
		rtrn.setName(gene1.getChr() + ":" + gene1.getStart() + "-" + gene1.getEnd() + ":" + gene2.getStart() + "-" + gene2.getEnd());
		return rtrn;
	}
	
	/**
	 * Get the aligned region as a bed line in genomic coordinates
	 * @param genes The pair of genes
	 * @return
	 */
	private String senseAlignmentAsBed(UnorderedGenePair genes) {
		Annotation alignedRegion = senseAlignmentAsFeature(genes);
		return alignedRegion.toBED(41,144,41);
	}
	
	/**
	 * Get the aligned region of the antisense alignment as a bed line in genomic coordinates
	 * @param genes The pair of genes
	 * @return
	 */
	private String antisenseAlignmentAsBed(UnorderedGenePair genes) {
		Annotation alignedRegion = antisenseAlignmentAsFeature(genes);
		return alignedRegion.toBED(144,59,144);
	}
	
	/**
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-b", "Bed file of features", true);
		p.addStringArg("-g", "Genome fasta file", true);
		p.addStringArg("-os", "Output bed file for sense alignments", false, null);
		p.addStringArg("-oa", "Output bed file for antisense alignments", false, null);
		p.addFloatArg("-ma", "Match score", false, DEFAULT_MATCH_SCORE);
		p.addFloatArg("-mi", "Mismatch score", false, DEFAULT_MISMATCH_SCORE);
		p.addFloatArg("-go", "Gap open penalty", false, DEFAULT_GAP_OPEN_PENALTY);
		p.addFloatArg("-ge", "Gap extend penalty", false, DEFAULT_GAP_EXTEND_PENALTY);
		p.addIntArg("-l", "Min alignment length", false, 0);
		p.addFloatArg("-id", "Min alignment idenity", false, 0);
		p.parse(args);
		String bedFile = p.getStringArg("-b");
		String fastaFile = p.getStringArg("-g");
		String outBedFileSense = p.getStringArg("-os");
		String outBedFileAntisense = p.getStringArg("-oa");
		float matchScore = p.getFloatArg("-ma");
		float mismatchScore = p.getFloatArg("-mi");
		float gapOpen = p.getFloatArg("-go");
		float gapExtend = p.getFloatArg("-ge");
		int minLength = p.getIntArg("-l");
		float minIdentity = p.getFloatArg("-id");
		
		
		FeatureSequenceAlignment fsa = new FeatureSequenceAlignment(bedFile, fastaFile, matchScore, mismatchScore, gapOpen, gapExtend);
		
		if(outBedFileSense != null) {
			fsa.writeAllSenseAlignmentsToBed(outBedFileSense, minLength, minIdentity);
		}
		
		if(outBedFileAntisense != null) {
			fsa.writeAllAntisenseAlignmentsToBed(outBedFileAntisense, minLength, minIdentity);
		}
		
	}
	
	
	/**
	 * An unordered pair of genes
	 * @author prussell
	 *
	 */
	public class UnorderedGenePair implements Comparable<UnorderedGenePair> {
		
		private Gene firstGene;
		private Gene secondGene;
		
		/**
		 * Determine which gene is first and which is second
		 * @param gene1 Gene 1
		 * @param gene2 Gene 2
		 */
		public UnorderedGenePair(Gene gene1, Gene gene2) {
			if(gene1.compareTo(gene2) < 0) {
				firstGene = gene1;
				secondGene = gene2;
			}
			else {
				firstGene = gene2;
				secondGene = gene1;
			}
		}
		
		/**
		 * @return The first gene
		 */
		public Gene getFirstGene() {
			return firstGene;
		}
		
		/**
		 * @return The second gene
		 */
		public Gene getSecondGene() {
			return secondGene;
		}

		@Override
		public int compareTo(UnorderedGenePair o) {
			if(firstGene.compareTo(o.getFirstGene()) == 0) {
				return secondGene.compareTo(o.getSecondGene());
			}
			return firstGene.compareTo(o.getFirstGene());
		}
		
		
		@Override
		public boolean equals(Object o) {
			UnorderedGenePair p = (UnorderedGenePair) o;
			return firstGene.equals(p.getFirstGene()) && secondGene.equals(p.getSecondGene());
		}

		@Override
		public int hashCode() {
			String str = firstGene.toString() + secondGene.toString();
			return str.hashCode();
		}
		
	}
	
	
}
