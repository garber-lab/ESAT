/**
 * 
 */
package broad.core.sequence;

import broad.core.parser.CommandLineParser;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.log4j.Logger;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.BasicAnnotation;
import nextgen.core.annotation.Gene;
import nextgen.core.annotation.Annotation.Strand;
import nextgen.core.utils.CountLogger;
import broad.pda.annotation.BEDFileParser;

/**
 * @author prussell
 *
 */
public class TranscribedRegions {
	
	private Map<String, Sequence> chromosomes;
	private Map<String, Collection<Gene>> genes;
	private Map<String, TreeMap<HalfOpenInterval, Strand>> geneStrands;
	private static Logger logger = Logger.getLogger(TranscribedRegions.class.getName());
	private String genomeFile;
	
	/**
	 * @param genomeFasta Genome fasta file
	 * @param bedFile Bed gene annotation
	 * @throws IOException
	 */
	public TranscribedRegions(String genomeFasta, String bedFile) throws IOException {
		genomeFile = genomeFasta;
		chromosomes = new TreeMap<String, Sequence>();
		// Load genes
		logger.info("Loading genes from file " + bedFile + "...");
		genes = BEDFileParser.loadDataByChr(new File(bedFile));
		int numGenes = 0;
		for(String chr : genes.keySet()) {
			numGenes += genes.get(chr).size();
		}
		logger.info("Loaded " + numGenes + " genes.");
		// Collapse overlapping exons and determine strands
		computeStrands();
	}
	
	private void loadGenome() throws IOException {
		chromosomes.clear();
		FastaSequenceIO fsio = new FastaSequenceIO(genomeFile);
		List<Sequence> chrs = fsio.loadAll();
		for(Sequence chr : chrs) {
			chromosomes.put(chr.getId(), chr);
		}
	}

	/**
	 * Get the genes
	 * @return Genes by chromosome
	 */
	public Map<String, Collection<Gene>> getGenes() {
		return genes;
	}
	
	/**
	 * Shift a position along the direction of transcription
	 * If the shifted position is not in the same block as the original position, this method currently throws an exception
	 * Throws an exception if the position is not in any transcribed interval
	 * @param chr Chromosome
	 * @param pos Original position
	 * @param offset Offset along direction of transcription (can be negative)
	 * @return Shifted position within the transcript
	 */
	public int shiftPosition(String chr, int pos, int offset) {
		for(HalfOpenInterval interval : geneStrands.get(chr).keySet()) {
			if(interval.contains(pos)) {
				Strand posStrand = geneStrands.get(chr).get(interval);
				int offsetPos = 0;
				if(posStrand.equals(Strand.UNKNOWN)) {
					throw new IllegalArgumentException("Strand unknown.");
				}
				if(posStrand.equals(Strand.POSITIVE)) {
					offsetPos = pos + offset;
					logger.debug("Strand " + posStrand.toString() + " orig " + pos + " offset " + offsetPos);
				} else if (posStrand.equals(Strand.NEGATIVE)) {
					offsetPos = pos - offset;
					logger.debug("Strand " + posStrand.toString() + " orig " + pos + " offset " + offsetPos);
				} 
				if(interval.contains(offsetPos)) {
					logger.debug("Returning " + offsetPos);
					return offsetPos;
				}
				throw new UnsupportedOperationException("Method does not support shifting position to another block.");
			}
		}
		throw new IllegalArgumentException("No interval in data structure contains the position.");
	}
	
	/**
	 * Get the transcribed orientation of the position based on the gene set
	 * Returns unknown if the position does not overlap a gene or strand is ambiguous
	 * @param chr Chromosome
	 * @param pos Position
	 * @return The orientation
	 */
	public Strand getOrientation(String chr, int pos) {
		for(HalfOpenInterval interval : geneStrands.get(chr).keySet()) {
			if(interval.contains(pos)) {
				return geneStrands.get(chr).get(interval);
			}
		}
		return Strand.UNKNOWN;
	}

	/**
	 * Get the oriented transcribed base at the position based on orientation provided
	 * or based on orientation of overlapping exon if orientation provided is null
	 * @param chr Chromosome
	 * @param pos Position
	 * @param strand Strand or null if get from annotation
	 * @return The transcribed base
	 * @throws IOException 
	 */
	public char getTranscribedBase(String chr, int pos, Strand strand) throws IOException {
		if(chromosomes.isEmpty()) loadGenome();
		Strand s = strand == null ? getOrientation(chr, pos) : strand;
		if(s.equals(Strand.UNKNOWN)) {
			return 'N';
		}
		Sequence subSeq = chromosomes.get(chr).getSubSequence("", pos, pos+1);
		if(s.equals(Strand.NEGATIVE)) {
			subSeq.reverse();
		}
		return subSeq.getSequenceBases().toUpperCase().charAt(0);		
	}
	
	/**
	 * Get the oriented transcribed base at the position based on orientation of overlapping exon
	 * @param chr Chromosome
	 * @param pos Position
	 * @return The transcribed base
	 * @throws IOException 
	 */
	public char getTranscribedBase(String chr, int pos) throws IOException {
		return getTranscribedBase(chr, pos, null);
	}
	
	/**
	 * Get the total count of each nucleotide
	 * @return Map of nucleotide as string to count
	 * @throws IOException
	 */
	private Map<String, Integer> getOverallNucleotideCounts() throws IOException {
		logger.info("Getting overall nucleotide counts for all transcribed exons...");
		int countA = 0;
		int countC = 0;
		int countG = 0;
		int countT = 0;
		int countN = 0;
		for(String chr : genes.keySet()) {
			logger.info(chr);
			// First list all postions contained in any exon
			Collection<Integer> allPositions = new TreeSet<Integer>();
			for(Gene gene : genes.get(chr)) {
				for(Annotation exon : gene.getBlocks()) {
					for(int i = exon.getStart(); i < exon.getEnd(); i++) {
						allPositions.add(Integer.valueOf(i));
					}
				}
			}
			// Now get transcribed nucleotide for each position and add to tracks
			int numPositions = allPositions.size();
			CountLogger c = new CountLogger(numPositions,10);
			for(Integer pos : allPositions) {
				c.advance();
				char base = getTranscribedBase(chr, pos.intValue());
				switch(base) {
				case 'A':
					countA++;
					break;
				case 'a':
					countA++;
					break;
				case 'C':
					countC++;
					break;
				case 'c':
					countC++;
					break;
				case 'G':
					countG++;
					break;
				case 'g':
					countG++;
					break;
				case 'T':
					countT++;
					break;
				case 't':
					countT++;
					break;
				case 'N':
					countN++;
					break;
				case 'n':
					countN++;
					break;
				default:
					logger.warn("Nucleotide " + base + " not recognized.");
					break;
				}
			}
		}
		Map<String, Integer> rtrn = new TreeMap<String, Integer>();
		rtrn.put("A", Integer.valueOf(countA));
		rtrn.put("C", Integer.valueOf(countC));
		rtrn.put("G", Integer.valueOf(countG));
		rtrn.put("T", Integer.valueOf(countT));
		rtrn.put("N", Integer.valueOf(countN));
		logger.info("Done getting nucleotide counts.");
		return rtrn;
	}
	
	
	/**
	 * Write overlapping intervals to bed file for QC
	 * @param bedFile Output bed file
	 * @throws IOException
	 */
	@SuppressWarnings("unused")
	private void writeOverlappingIntervalsToBed(String bedFile) throws IOException {
		FileWriter w = new FileWriter(bedFile);
		for(String chr : geneStrands.keySet()) {
			for(HalfOpenInterval interval : geneStrands.get(chr).keySet()) {
				Strand strand = geneStrands.get(chr).get(interval);
				Annotation a = new BasicAnnotation(chr, (int)interval.getStart(), (int)interval.getEnd());
				a.setOrientation(strand);
				w.write(a.toBED() + "\n");
			}
		}
		w.close();
	}

	
	/**
	 * Collapse overlapping exons and save strand for each collapsed interval
	 */
	private void computeStrands() {
		logger.info("Collapsing overlapping exons and identifying strand...");
		geneStrands = new TreeMap<String, TreeMap<HalfOpenInterval,Strand>>();
		for(String chr : genes.keySet()) {
			logger.info(chr);
			TreeMap<HalfOpenInterval, Strand> chrStrands = new TreeMap<HalfOpenInterval, Strand>();
			for(Gene gene : genes.get(chr)) {
				Strand strand = gene.getOrientation();
				for(HalfOpenInterval block : exonBlocks(gene)) {
					boolean merged = false;
					for(HalfOpenInterval existingInterval : chrStrands.keySet()) {
						if(block.overlaps(existingInterval)) {
							existingInterval.merge(block);
							if(!chrStrands.get(existingInterval).equals(strand)) {
								chrStrands.put(existingInterval, Strand.UNKNOWN);
							}
							merged = true;
						}
					}
					if(!merged) {
						chrStrands.put(block, strand);
					}
				}
			}
			geneStrands.put(chr, chrStrands);
		}
		logger.info("Done identifying intervals and strands.");
	}

	private Collection<HalfOpenInterval> exonBlocks(Gene gene) {
		Collection<HalfOpenInterval> rtrn = new TreeSet<HalfOpenInterval>();
		for(Annotation exon : gene.getBlocks()) {
			rtrn.add(new HalfOpenInterval(exon.getStart(), exon.getEnd()));
		}
		return rtrn;
	}
	
	
	
	/**
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-g", "Genome bed file", true);
		p.addStringArg("-b", "Annotation bed file", true);
		p.addStringArg("-oc", "For each nucleotide get total number of transcribed positions and write table to this file", false, null);
		p.parse(args);
		String genomeFasta = p.getStringArg("-g");
		String bedFile = p.getStringArg("-b");
		String outCountsTable = p.getStringArg("-oc");
		if(p.getFlagsAndValues().isEmpty()) {
			p.printHelpMessage();
			System.exit(-1);
		}
		
		TranscribedRegions t = new TranscribedRegions(genomeFasta, bedFile);
		
		if(outCountsTable != null) {
			Map<String, Integer> countsByBase = t.getOverallNucleotideCounts();
			logger.info("Writing nucleotide counts to file " + outCountsTable + ".");
			FileWriter w = new FileWriter(outCountsTable);
			for(String base : countsByBase.keySet()) {
				w.write(base + "\t" + countsByBase.get(base).toString() + "\n");
			}
			w.close();
			logger.info("Done writing file.");
		}
		
		logger.info("All done.");
		
	}
	
	
	
	
	private class HalfOpenInterval implements Comparable<HalfOpenInterval> {

		private double start;
		private double end;
		
		/**
		 * @param startPos Start position (inclusive)
		 * @param endPos End position (exclusive)
		 */
		public HalfOpenInterval(double startPos, double endPos) {
			start = startPos;
			end = endPos;
		}
		
		/**
		 * Get start position
		 * @return Start position (inclusive)
		 */
		public double getStart() {
			return start;
		}
		
		/**
		 * Get end position
		 * @return End position (exclusive)
		 */
		public double getEnd() {
			return end;
		}
		
		/**
		 * Set start
		 * @param s New start
		 */
		public void setStart(double s) {
			start = s;
		}
		
		/**
		 * Set end
		 * @param e New end
		 */
		public void setEnd(double e) {
			end = e;
		}
		
		/**
		 * Whether this interval overlaps other interval
		 * @param other Other interval
		 * @return True iff the intervals share an inclusive position
		 */
		public boolean overlaps(HalfOpenInterval other) {
			if(start < other.getStart()) {
				return end > other.getStart();
			}
			return start < other.getEnd();
		}
		
		/**
		 * Merge with the other interval if they overlap
		 * Reset start and end positions to the union of the intervals
		 * If they do not overlap, do nothing
		 * @param other Other interval to merge
		 */
		public void merge(HalfOpenInterval other) {
			if(!overlaps(other)) {
				return;
			}
			setStart(Math.min(start, other.getStart()));
			setEnd(Math.max(end, other.getEnd()));
		}
		
		/**
		 * Whether the interval contains the number
		 * @param i Number to check
		 * @return True iff i is strictly contained in the interval
		 */
		public boolean contains(int i) {
			return i >= start && i < end;
		}
		
		/**
		 * Whether the interval contains the number
		 * @param d Number to check
		 * @return True iff d is strictly contained in the interval
		 */
		@SuppressWarnings("unused")
		public boolean contains(double d) {
			return d >= start && d < end;
		}
		
		@Override
		public boolean equals(Object o) {
			HalfOpenInterval oi = (HalfOpenInterval)o;
			return start == oi.getStart() && end == oi.getEnd();
		}
		
		@Override
		public String toString() {
			String s = Double.valueOf(start).toString();
			String e = Double.valueOf(end).toString();
			return s + "-" + e;
		}
		
		@Override
		public int hashCode() {
			return toString().hashCode();
		}
		
		@Override
		public int compareTo(HalfOpenInterval other) {
			if(start == other.getStart()) {
				return 0;
			}
			if(start < other.getStart()) {
				return -1;
			}
			return 1;
		}

	}

}
