/**
 * 
 */
package broad.core.sequence;

import broad.core.parser.CommandLineParser;
import broad.core.parser.StringParser;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;


import broad.core.math.MannWhitney;
import broad.core.math.Statistics;
import broad.pda.annotation.BEDFileParser;

import net.sf.samtools.util.CloseableIterator;
import nextgen.core.alignment.AbstractPairedEndAlignment.TranscriptionRead;
import nextgen.core.alignment.Alignment;
import nextgen.core.annotation.Annotation.Strand;
import nextgen.core.annotation.Gene;
import nextgen.core.coordinatesystem.TranscriptomeSpace;
import nextgen.core.model.TranscriptomeSpaceAlignmentModel;
import nextgen.core.utils.AnnotationUtils;
import nextgen.core.utils.FileUtil;

/**
 * @author prussell
 *
 */
public class TranscribedSequenceComposition {

	private Map<String, Collection<Gene>> genes;
	private Map<String, Collection<Gene>> regions;
	private TranscribedRegions transcribedRegions;
	private TranscriptomeSpaceAlignmentModel data;
	private Map<Gene, Gene> regionToGeneMap;
	private static Logger logger = Logger.getLogger(TranscribedSequenceComposition.class.getName());
	protected Collection<String> allNucleotides;
	private Map<Integer, Map<Gene, Map<String, Collection<NucleotidePosition>>>> regionBeginningPositionsByOffsetAndNucleotide;
	
	
	/**
	 * @param bamFile Bam file
	 * @param geneBedFile Bed file of genes
	 * @param genomeFasta Genome fasta file
	 * @param transcriptionRead Transcription read
	 * @throws IOException
	 */
	public TranscribedSequenceComposition(String bamFile, String geneBedFile, String genomeFasta, TranscriptionRead transcriptionRead) throws IOException {
		this(bamFile, geneBedFile, geneBedFile, genomeFasta, transcriptionRead);
	}
	
	/**
	 * @param bamFile Bam file
	 * @param geneBedFile Bed file of genes
	 * @param bedFileRegionsOfInterest Bed file of regions within genes
	 * @param genomeFasta Genome fasta file
	 * @param transcriptionRead Transcription read
	 * @throws IOException
	 */
	public TranscribedSequenceComposition(String bamFile, String geneBedFile, String bedFileRegionsOfInterest, String genomeFasta, TranscriptionRead transcriptionRead) throws IOException {
		logger.info("");
		logger.info("Loading genes and regons...");
		genes = BEDFileParser.loadDataByChr(new File(geneBedFile));
		regions = BEDFileParser.loadDataByChr(new File(bedFileRegionsOfInterest));
		logger.info("");
		logger.info("Loading alignment data...");
		data = new TranscriptomeSpaceAlignmentModel(bamFile, new TranscriptomeSpace(genes), transcriptionRead);
		logger.info("");
		logger.info("Loading sequences...");
		transcribedRegions = new TranscribedRegions(genomeFasta, geneBedFile);
		logger.info("");
		logger.info("Identifying regions with parent genes...");
		regionToGeneMap = AnnotationUtils.mapChildToLargestParent(regions, genes);
		allNucleotides = new TreeSet<String>();
		allNucleotides.add("A");
		allNucleotides.add("C");
		allNucleotides.add("G");
		allNucleotides.add("T");
		regionBeginningPositionsByOffsetAndNucleotide = new TreeMap<Integer, Map<Gene, Map<String, Collection<NucleotidePosition>>>>();
	}
	
	
	/**
	 * Get map of nucleotide bases to beginning positions of reads starting with the nucleotide
	 * @param region The region
	 * @param offset Offset along the parent transcript from the beginning position of each read, or null if region does not have parent
	 * @throws IOException 
	 */
	private Map<String, Collection<NucleotidePosition>> getBeginningPositionsByNucleotide(Gene region, int offset) throws IOException {
		if(!regionBeginningPositionsByOffsetAndNucleotide.containsKey(Integer.valueOf(offset))) {
			Map<Gene, Map<String, Collection<NucleotidePosition>>> mapThisOffset = new TreeMap<Gene, Map<String, Collection<NucleotidePosition>>>();
			regionBeginningPositionsByOffsetAndNucleotide.put(Integer.valueOf(offset), mapThisOffset);
		}
		if(!regionBeginningPositionsByOffsetAndNucleotide.get(Integer.valueOf(offset)).containsKey(region)) {
			calculateBeginningPositionsByNucleotide(region, offset);
		}
		return regionBeginningPositionsByOffsetAndNucleotide.get(Integer.valueOf(offset)).get(region);
	}
	
	/**
	 * Get transcribed base of parent gene at first position of alignment
	 * @param parent Parent gene
	 * @param read Alignment 
	 * @param offset Offset along the parent transcript from the beginning position
	 * @return The transcribed base or N if read orientation does not match parent gene strand
	 * @throws IOException 
	 */
	public char getTranscribedBaseAtBeginningPosition(Gene parent, Alignment read, int offset) throws IOException {
		NucleotidePosition n = getBeginningPositionAndBase(parent, read, offset);
		if(n == null) {
			return 'N';
		}
		return n.getBase();
	}
	
	
	/**
	 * Get transcribed base and position of parent gene at first position of alignment
	 * @param parent Parent gene
	 * @param read Alignment 
	 * @param offset Offset along the parent transcript from the beginning position
	 * @return The transcribed position and nucleotide
	 * @throws IOException 
	 */
	private NucleotidePosition getBeginningPositionAndBase(Gene parent, Alignment read, int offset) throws IOException {
		Strand parentStrand = parent.getOrientation();
		Strand readStrand = read.getFragmentStrand();
		String chr = parent.getChr();
		logger.debug("Read " + read.toBED());
		if(!readStrand.equals(parentStrand)) {
			logger.debug("Returning null for read " + read.getName() + " because strand does not match gene " + parent.getName());
			return null;
		}
		try {
			int beginningPos = readStrand.equals(Strand.POSITIVE) ? read.getStart() : parent.genomicToGenomicPositionWithOffset(read.getEnd(), 1);
			int newPos = parent.genomicToGenomicPositionWithOffset(beginningPos, offset);
			NucleotidePosition n = new NucleotidePosition(chr, newPos, transcribedRegions.getTranscribedBase(chr, newPos, parentStrand));
			logger.debug("fragment_strand=" + readStrand.toString() + "\tnuc=" + n.getBase() + "\tbeginning_pos=" + beginningPos + "\toffset_pos=" + newPos);
			return n;
			
		} catch (IllegalArgumentException e) {
			logger.debug("Returning null for read " + read.getName() + " because caught exception: " + e.getMessage());
			return null;
		}
	}
	
	/**
	 * Calculate map of nucleotide bases to beginning positions of reads starting with the nucleotide
	 * @param region The region
	 * @param offset Offset along the parent transcript from the beginning position of each read
	 * @throws IOException 
	 */
	private void calculateBeginningPositionsByNucleotide(Gene region, int offset) throws IOException {
		logger.debug("Calculating beginning positions by nucleotide (offset = " + offset + ") for region " + region.toBED());
		Strand regionStrand = region.getStrand();
		if(regionStrand.equals(Strand.UNKNOWN)) {
			throw new IllegalArgumentException("Strand must be known: " + region.toBED());
		}
		if(!regionToGeneMap.containsKey(region)) {
			logger.warn("Region has no child-parent mappings: " + region.getName());
			regionBeginningPositionsByOffsetAndNucleotide.get(Integer.valueOf(offset)).put(region, null);
			return;
		}
		Gene parent = regionToGeneMap.get(region);
		logger.debug("Parent gene is " + parent.toBED());
		Strand parentStrand = parent.getOrientation();
		if(parentStrand.equals(Strand.UNKNOWN)) {
			throw new IllegalArgumentException("Strand must be known " + parent.toBED());
		}
		if(!regionStrand.equals(parentStrand)) {
			throw new IllegalStateException("Region and parent gene strands must match.");
		}
		CloseableIterator<Alignment> reads = data.getOverlappingReads(region, false);
		Map<String, Collection<NucleotidePosition>> thisRegion = new TreeMap<String, Collection<NucleotidePosition>>();
		while(reads.hasNext()) {
			Alignment read = reads.next();
			logger.debug("Read " + read.toBED());
			try {
				NucleotidePosition n = getBeginningPositionAndBase(parent, read, offset);
				if(n == null) {
					continue;
				}
				String baseString = n.getBaseString();
				if(!thisRegion.containsKey(baseString)) {
					Collection<NucleotidePosition> cn = new ArrayList<NucleotidePosition>();
					thisRegion.put(baseString, cn);
				}
				thisRegion.get(baseString).add(n);
			} catch (IllegalArgumentException e) {
				logger.debug("Caught exception: " + e.getMessage());
				continue;
			}
		}
		reads.close();
		regionBeginningPositionsByOffsetAndNucleotide.get(Integer.valueOf(offset)).put(region, thisRegion);
	}
	
	/**
	 * Get median beginning position of reads starting with any nucleotide in the set
	 * @param region The region
	 * @param offset Offset along the parent transcript from the beginning position of each read
	 * @param nucleotideSet The set of nucleotides to use
	 * @param convertToTranscriptCoordinates Convert beginning positions to transcript coordinates with respect to parent gene
	 * @return Map of base to median beginning position of positions having the base, or -1 if region does not have parent
	 * @throws IOException 
	 */
	private double medianBeginningPositionByNucleotideSet(Gene region, int offset, TreeSet<String> nucleotideSet, boolean convertToTranscriptCoordinates) throws IOException {
		Map<String, Collection<NucleotidePosition>> allPositions = getBeginningPositionsByNucleotide(region, offset);
		if(allPositions == null) {
			return -1;
		}
		Gene parent = regionToGeneMap.get(region);
		Collection<Double> positions = new ArrayList<Double>();
		for(String nuc : nucleotideSet) {
			if(!allPositions.containsKey(nuc)) {
				continue;
			}
			for(NucleotidePosition pos : allPositions.get(nuc)) {
				int p = convertToTranscriptCoordinates ? parent.genomicToTranscriptPosition(pos.getPos()) : pos.getPos();
				positions.add(Double.valueOf(p));
			}
		}
		return Statistics.medianCollection(positions);
	}

	
	/**
	 * Get median beginning position of reads starting with each nucleotide
	 * @param region The region
	 * @param offset Offset along the parent transcript from the beginning position of each read
	 * @param convertToTranscriptCoordinates Convert beginning positions to transcript coordinates with respect to parent gene
	 * @return Map of base to median beginning position of positions having the base
	 * @throws IOException 
	 */
	private Map<String, Double> medianBeginningPositionByNucleotide(Gene region, int offset, boolean convertToTranscriptCoordinates) throws IOException {
		Map<String, Double> rtrn = new TreeMap<String, Double>();
		for(String nuc : allNucleotides) {
			TreeSet<String> thisNuc = new TreeSet<String>();
			thisNuc.add(nuc);
			rtrn.put(nuc, Double.valueOf(medianBeginningPositionByNucleotideSet(region, offset, thisNuc, convertToTranscriptCoordinates)));
		}
		return rtrn;
	}
	
	/**
	 * Get the number of fragments overlapping the region with beginning position at the nucleotide
	 * @param region The region
	 * @param offset Offset along the parent transcript from the beginning position of each read
	 * @return Map of base to number of overlapping fragments whose first position is the base, or null if region does not have parent
	 * @throws IOException 
	 */
	private Map<String, Integer> countBeginningPositionByNucleotide(Gene region, int offset) throws IOException {
		Map<String, Collection<NucleotidePosition>> allPositions = getBeginningPositionsByNucleotide(region, offset);
		if(allPositions == null) {
			return null;
		}
		Map<String, Integer> rtrn = new TreeMap<String, Integer>();
		for(String s : allPositions.keySet()) {
			rtrn.put(s, Integer.valueOf(allPositions.get(s).size()));
		}
		return rtrn;
	}
	
	/**
	 * Get the total number of fragments overlapping all regions with beginning position at the nucleotide
	 * @param offset Offset along the parent transcript from the beginning position of each read
	 * @return Map of base to total number of fragments whose first position is the base
	 * @throws IOException 
	 */
	@SuppressWarnings("unused")
	private Map<String, Integer> totalBeginningPositionByNucleotide(int offset) throws IOException {
		return totalBeginningPositionByNucleotide(offset, null);
	}
	
	/**
	 * Get the total number of fragments (possibly on one chromosome only) overlapping all regions with beginning position at the nucleotide
	 * @param offset Offset along the parent transcript from the beginning position of each read
	 * @param chrListFile Single chromosome to use or null if all chromosomes
	 * @return Map of base to total number of fragments whose first position is the base
	 * @throws IOException 
	 */
	private Map<String, Integer> totalBeginningPositionByNucleotide(int offset, String chrListFile) throws IOException {
		logger.info("");
		logger.info("Calculating total number of fragments with beginning position at each nucleotide...");
		Collection<String> chrs = new TreeSet<String>();
		if(chrListFile == null) {
			chrs.addAll(regions.keySet());
		} else {
			chrs.addAll(FileUtil.fileLinesAsList(chrListFile));
		}
		Map<String, Integer> rtrn = new TreeMap<String, Integer>();
		for(String n : allNucleotides) {
			rtrn.put(n, Integer.valueOf(0));
		}
		for(String c : chrs) {
			logger.info(c);
			if(!regions.containsKey(c)) continue;
			for(Gene region : regions.get(c)) {
				Map<String, Integer> counts = countBeginningPositionByNucleotide(region, offset);
				if(counts == null) {
					continue;
				}
				for(String s : counts.keySet()) {
					int nextVal = counts.containsKey(s) ? Math.max(0, counts.get(s).intValue()) : 0;
					rtrn.put(s, Integer.valueOf(rtrn.get(s).intValue() + nextVal));
				}
			}
		}
		logger.info("Done calculating totals.");
		return rtrn;
	}
	
	/**
	 * Write to a file the total number of fragments overlapping all regions with beginning position at the nucleotide
	 * @param outFile Output file
	 * @param offset Offset along the parent transcript from the beginning position of each read
	 */
	@SuppressWarnings("unused")
	private void writeTotalCountsBeginningPositionByNucleotide(String outFile, int offset) throws IOException {
		writeTotalCountsBeginningPositionByNucleotide(outFile, offset, null);
	}
	
	/**
	 * Write to a file the total number of fragments (possibly on certain chromosomes only) overlapping all regions with beginning position at the nucleotide
	 * @param outFile Output file
	 * @param offset Offset along the parent transcript from the beginning position of each read
	 * @param chrListFile File containing list of chromosomes to use or null if all chromosomes
	 */
	private void writeTotalCountsBeginningPositionByNucleotide(String outFile, int offset, String chrListFile) throws IOException {
		logger.info("");
		logger.info("Writing total number of fragments with beginning position at each nucleotide to file " + outFile + "...");
		Map<String, Integer> counts = totalBeginningPositionByNucleotide(offset, chrListFile);
		FileWriter w = new FileWriter(outFile);
		for(String s : counts.keySet()) {
			w.write(s + "\t" + counts.get(s) + "\n");
		}
		w.close();
	}
	
	/**
	 * For each region write median beginning position of reads starting with each nucleotide
	 * @param outFile Output table
	 * @param offset Offset along the parent transcript from the beginning position of each read
	 * @param convertToTranscriptCoordinates Convert beginning positions to transcript coordinates with respect to parent gene
	 * @throws IOException
	 */
	@SuppressWarnings("unused")
	private void writeMedianBeginningPositionByNucleotide(String outFile, int offset, boolean convertToTranscriptCoordinates) throws IOException {
		writeMedianBeginningPositionByNucleotide(outFile, offset, convertToTranscriptCoordinates, null);
	}
		
	/**
	 * Get start positions of fragments beginning with any nucleotide in a set
	 * @param region The region
	 * @param nucleotideSet The set of nucleotides to use
	 * @param offset Offset along the parent transcript from the beginning position of each read
	 * @return Start positions of all fragments beginning with any nucleotide in the set, or null if region does not have parent
	 * @throws IOException 
	 */
	private ArrayList<Double> getFragmentStartPositionsInNucleotideSet(Gene region, TreeSet<String> nucleotideSet, int offset) throws IOException {
		ArrayList<Double> rtrn = new ArrayList<Double>();
		Map<String, Collection<NucleotidePosition>> positionsByNuc = getBeginningPositionsByNucleotide(region, offset);
		if(positionsByNuc == null) {
			return null;
		}
		for(String nuc : nucleotideSet) {
			if(!positionsByNuc.containsKey(nuc)) {
				continue;
			}
			for(NucleotidePosition pos : positionsByNuc.get(nuc)) {
				rtrn.add(Double.valueOf(pos.getPos()));
			}
		}
		return rtrn;
	}
	
	
	/**
	 * Get Mann Whitney P value of start positions of fragments beginning with two sets of nucleotides, and with a minimum number of fragments applying to each nucleotide set
	 * @param region Region for overlapping fragments
	 * @param offset Offset along the parent transcript from the beginning position of each read
	 * @param nucleotideSet1 One set of nucleotides
	 * @param nucleotideSet2 Other set of nucleotides
	 * @return Mann Whitney P value or -1 if not enough of each nucleotide or region does not have parent transcript
	 * @throws IOException 
	 */
	private double mannWhitneyPvalBeginningPositions(Gene region, int offset, TreeSet<String> nucleotideSet1, TreeSet<String> nucleotideSet2) throws IOException {
		boolean nonemptyIntersection = nucleotideSet1.removeAll(nucleotideSet2);
		if(nonemptyIntersection) {
			throw new IllegalArgumentException("The two sets of nucleotides must not intersect.");
		}
		ArrayList<Double> positionSet1 = getFragmentStartPositionsInNucleotideSet(region, nucleotideSet1, offset);
		ArrayList<Double> positionSet2 = getFragmentStartPositionsInNucleotideSet(region, nucleotideSet2, offset);
		if(positionSet1 == null || positionSet2 == null) {
			return -1;
		}
		MannWhitney mw = new MannWhitney(positionSet1, positionSet2);
		return mw.getPvalue();
	}
	
	/**
	 * For each region write median start positions and Mann Whitney P value for hypothesis that fragments starting with nucleotides in two sets start at different places
	 * @param outFile Output file
	 * @param offset Offset along the parent transcript from the beginning position of each read
	 * @param nucleotideSet1 One set of nucleotides
	 * @param nucleotideSet2 Other set of nucleotides
	 * @throws IOException
	 */
	@SuppressWarnings("unused")
	private void writeBeginningPositionSeparation(String outFile, int offset, TreeSet<String> nucleotideSet1, TreeSet<String> nucleotideSet2) throws IOException {
		writeBeginningPositionSeparation(outFile, offset, null, nucleotideSet1, nucleotideSet2, 1);
	}
	
	/**
	 * For each region write median start positions and Mann Whitney P value for hypothesis that fragments starting with nucleotides in two sets start at different places
	 * @param outFile Output file
	 * @param offset Offset along the parent transcript from the beginning position of each read
	 * @param chrListFile File containing list of chromosomes to write
	 * @param nucleotideSet1 One set of nucleotides
	 * @param nucleotideSet2 Other set of nucleotides
	 * @throws IOException
	 */
	@SuppressWarnings("unused")
	private void writeBeginningPositionSeparation(String outFile, int offset, String chrListFile, TreeSet<String> nucleotideSet1, TreeSet<String> nucleotideSet2) throws IOException {
		writeBeginningPositionSeparation(outFile, offset, chrListFile, nucleotideSet1, nucleotideSet2, 1);
	}
	
	/**
	 * For each region write median start positions and Mann Whitney P value for hypothesis that fragments starting with nucleotides in two sets start at different places
	 * @param outFile Output file
	 * @param offset Offset along the parent transcript from the beginning position of each read
	 * @param nucleotideSet1 One set of nucleotides
	 * @param nucleotideSet2 Other set of nucleotides
	 * @param minNumberOfEach For Mann Whitney, minimum number of fragments starting with nucleotides in each set
	 * @throws IOException
	 */
	@SuppressWarnings("unused")
	private void writeBeginningPositionSeparation(String outFile, int offset, TreeSet<String> nucleotideSet1, TreeSet<String> nucleotideSet2, int minNumberOfEach) throws IOException {
		writeBeginningPositionSeparation(outFile, offset, null, nucleotideSet1, nucleotideSet2, minNumberOfEach);
	}
	
	/**
	 * For each region write median start positions and Mann Whitney P value for hypothesis that fragments starting with nucleotides in two sets start at different places
	 * @param outFile Output file
	 * @param offset Offset along the parent transcript from the beginning position of each read
	 * @param chrListFile File containing list of chromosomes to write
	 * @param nucleotideSet1 One set of nucleotides
	 * @param nucleotideSet2 Other set of nucleotides
	 * @param minNumberOfEach For Mann Whitney, minimum number of fragments starting with nucleotides in each set
	 * @throws IOException
	 */
	private void writeBeginningPositionSeparation(String outFile, int offset, String chrListFile, TreeSet<String> nucleotideSet1, TreeSet<String> nucleotideSet2, int minNumberOfEach) throws IOException {
		logger.info("");
		logger.info("Writing median start position in transcript coordinates and Mann-Whitney P value of start position in two sets of nucleotides for each region to file " + outFile + "...");
		if(minNumberOfEach < 1) {
			throw new IllegalArgumentException("Min number of each must be at least 1.");
		}
		Collection<String> chrs = new TreeSet<String>();
		if(chrListFile == null) {
			chrs.addAll(regions.keySet());
		} else {
			chrs.addAll(FileUtil.fileLinesAsList(chrListFile));
		}
		FileWriter w = new FileWriter(outFile);
		String header = "parent_gene";
		header += "\tregion";
		header += "\tnum_fragments_" + commaSeparatedString(nucleotideSet1);
		header += "\tnum_fragments_" + commaSeparatedString(nucleotideSet2);
		header += "\tmedian_beginning_position_transcript_coord_" + commaSeparatedString(nucleotideSet1);
		header += "\tmedian_beginning_position_transcript_coord_" + commaSeparatedString(nucleotideSet2);
		header += "\tmedian_" + commaSeparatedString(nucleotideSet1) + "_minus_median_" + commaSeparatedString(nucleotideSet2);
		header += "\tmann_whitney_pval";
		w.write(header + "\n");
		for(String c : chrs) {
			logger.info(c);
			if(!regions.containsKey(c)) continue;
			for(Gene region : regions.get(c)) {
				ArrayList<Double> positionSet1 = getFragmentStartPositionsInNucleotideSet(region, nucleotideSet1, offset);
				ArrayList<Double> positionSet2 = getFragmentStartPositionsInNucleotideSet(region, nucleotideSet2, offset);
				if(positionSet1 == null || positionSet2 == null) {
					continue;
				}
				int numSet1 = positionSet1.size();
				int numSet2 = positionSet2.size();
				if(numSet1 < minNumberOfEach || numSet2 < minNumberOfEach) {
					continue;
				}
				double median1 = medianBeginningPositionByNucleotideSet(region, offset, nucleotideSet1, true);
				double median2 = medianBeginningPositionByNucleotideSet(region, offset, nucleotideSet2, true);
				double pval = mannWhitneyPvalBeginningPositions(region, offset, nucleotideSet1, nucleotideSet2);
				String line = regionToGeneMap.get(region).getName();
				line += "\t" + region.getName();
				line += "\t" + numSet1;
				line += "\t" + numSet2;
				line += "\t" + median1;
				line += "\t" + median2;
				double diff = median1 - median2;
				line += "\t" + diff;
				line += "\t" + pval;
				w.write(line + "\n");
			}
		}
		w.close();
		logger.info("Done writing file.");
	}
	
	/**
	 * For each region write median beginning position of reads starting with each nucleotide, possible only for certain chromosomes
	 * @param outFile Output table
	 * @param offset Offset along the parent transcript from the beginning position of each read
	 * @param convertToTranscriptCoordinates Convert beginning positions to transcript coordinates with respect to parent gene
	 * @param chrListFile File containing list of chromosomes to write or null if all chromosomes
	 * @throws IOException
	 */
	private void writeMedianBeginningPositionByNucleotide(String outFile, int offset, boolean convertToTranscriptCoordinates, String chrListFile) throws IOException {
		logger.info("");
		logger.info("Writing fragment median start positions by nucleotide for each region to file " + outFile + "...");
		Collection<String> chrs = new TreeSet<String>();
		if(chrListFile == null) {
			chrs.addAll(regions.keySet());
		} else {
			chrs.addAll(FileUtil.fileLinesAsList(chrListFile));
		}
		FileWriter w = new FileWriter(outFile);
		String header = "parent_gene\tregion";
		
		for(String nuc : allNucleotides) {
			header += "\t" + nuc;
		}
		w.write(header + "\n");
		for(String c : chrs) {
			logger.info(c);
			if(!regions.containsKey(c)) continue;
			for(Gene region : regions.get(c)) {
				Map<String, Double> medianPositions = medianBeginningPositionByNucleotide(region, offset, convertToTranscriptCoordinates);
				String line = regionToGeneMap.get(region).getName() + "\t" + region.getName();
				for(String nuc : allNucleotides) {
					if(medianPositions.containsKey(nuc)) {
						line += "\t" + medianPositions.get(nuc);
					} else {
						line += "\t-";
					}
				}
				w.write(line + "\n");
			}
		}
		w.close();
		logger.info("Done writing file.");
	}
	
	/**
	 * For each region write number of reads starting with each nucleotide
	 * @param outFile Output table
	 * @param offset Offset along the parent transcript from the beginning position of each read
	 * @throws IOException
	 */
	@SuppressWarnings("unused")
	private void writeCountBeginningPositionByNucleotide(String outFile, int offset) throws IOException {
		writeCountBeginningPositionByNucleotide(outFile, offset, null);
	}
	
	/**
	 * For each region write number of reads starting with each nucleotide, possible only for certain chromosomes
	 * @param outFile Output table
	 * @param offset Offset along the parent transcript from the beginning position of each read
	 * @param chrListFile File containing list of chromosomes to write or null if all chromosomes
	 * @throws IOException
	 */
	private void writeCountBeginningPositionByNucleotide(String outFile, int offset, String chrListFile) throws IOException {
		logger.info("");
		logger.info("Writing number of overlapping fragments starting with each nucleotide for each region to file " + outFile + "...");
		Collection<String> chrs = new TreeSet<String>();
		if(chrListFile == null) {
			chrs.addAll(regions.keySet());
		} else {
			chrs.addAll(FileUtil.fileLinesAsList(chrListFile));
		}
		FileWriter w = new FileWriter(outFile);
		String header = "parent_gene\tregion";
		
		for(String nuc : allNucleotides) {
			header += "\t" + nuc;
		}
		w.write(header + "\n");
		for(String c : chrs) {
			logger.info(c);
			if(!regions.containsKey(c)) continue;
			for(Gene region : regions.get(c)) {
				Map<String, Integer> counts = countBeginningPositionByNucleotide(region, offset);
				if(counts == null) {
					continue;
				}
				String line = regionToGeneMap.get(region).getName() + "\t" + region.getName();
				for(String nuc : allNucleotides) {
					if(counts.containsKey(nuc)) {
						line += "\t" + counts.get(nuc);
					} else {
						line += "\t-";
					}
				}
				w.write(line + "\n");
			}
		}
		w.close();
		logger.info("Done writing file.");
	}

	private static String commaSeparatedString(TreeSet<String> set) {
		String rtrn = "";
		for(String s : set) {
			rtrn += s + ",";
		}
		if(rtrn.length() > 0) {
			return rtrn.substring(0, rtrn.length()-1);
		}
		return rtrn;
	}
	
	private static TreeSet<String> commaSeparatedTokens(String commaSeparatedList) {
		StringParser p = new StringParser();
		p.parse(commaSeparatedList, ",");
		String[] tokens = p.getStringArray();
		TreeSet<String> rtrn = new TreeSet<String>();
		for(int i=0; i < tokens.length; i++) {
			String t = tokens[i];
			if(t.length() != 1) {
				throw new IllegalArgumentException("Must provided comma separated list of single nucleotides. String " + t + " is not valid.");
			}
			rtrn.add(t);
		}
		return rtrn;
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-b", "Bam file", true);
		p.addStringArg("-g", "Gene bed file", true);
		p.addStringArg("-r", "Region bed file", true);
		p.addStringArg("-f", "Genome fasta file", true);
		p.addIntArg("-s", "Offset along the parent transcript from the beginning position of each read (positive = 3' direction; negative = 5' direction", false, 0);
		p.addStringArg("-otc", "Output table of number of fragments beginning with each nucleotide for each region", false, null);
		p.addStringArg("-ots", "Output table of difference in fragment beginnings between two nucleotide sets", false, null);
		p.addStringArg("-oc", "Output counts of number of fragments beginning with each nucleotide", false, null);
		p.addBooleanArg("-d", "Debug logging", false, false);
		p.addStringArg("-c", "File containing list of chromosomes to write", false, null);
		p.addBooleanArg("-ft", "First read transcription strand", false, false);
		p.addIntArg("-mn", "Min number of fragments with start nucleotide in each set for difference table", false, 1);
		p.addStringArg("-n1", "One comma separated set of nucleotides for difference table", false, "A,G");
		p.addStringArg("-n2", "Other comma separated set of nucleotides for difference table", false, "C,T");
		p.parse(args);
		String bamFile = p.getStringArg("-b");
		String geneBed = p.getStringArg("-g");
		String regionBed = p.getStringArg("-r");
		String genomeFasta = p.getStringArg("-f");
		int offset = p.getIntArg("-s");
		String outCounts = p.getStringArg("-oc");
		boolean debug = p.getBooleanArg("-d");
		String chrListFile = p.getStringArg("-c");
		boolean firstReadTranscriptionStrand = p.getBooleanArg("-ft");
		TranscriptionRead tr = firstReadTranscriptionStrand ? TranscriptionRead.FIRST_OF_PAIR : TranscriptionRead.SECOND_OF_PAIR;
		String outTableCount = p.getStringArg("-otc");
		if(debug) logger.setLevel(Level.DEBUG);
		String outTableSeparation = p.getStringArg("-ots");
		int minNumberOfEach = p.getIntArg("-mn");
		TreeSet<String> nucleotideSet1 = commaSeparatedTokens(p.getStringArg("-n1"));
		TreeSet<String> nucleotideSet2 = commaSeparatedTokens(p.getStringArg("-n2"));
		
		TranscribedSequenceComposition t = new TranscribedSequenceComposition(bamFile, geneBed, regionBed, genomeFasta, tr);
		if(outTableCount != null) t.writeCountBeginningPositionByNucleotide(outTableCount, offset, chrListFile);
		if(outCounts != null) t.writeTotalCountsBeginningPositionByNucleotide(outCounts, offset, chrListFile);
		if(outTableSeparation != null) t.writeBeginningPositionSeparation(outTableSeparation, offset, chrListFile, nucleotideSet1, nucleotideSet2, minNumberOfEach);
		
		logger.info("");
		logger.info("All done.");
		
	}

	private class NucleotidePosition {
		
		private String chr;
		private int pos;
		private char nuc;
		
		@SuppressWarnings("unused")
		protected NucleotidePosition(String chrName, int position) {
			chr = chrName;
			pos = position;
		}
		
		protected NucleotidePosition(String chrName, int position, char base) {
			chr = chrName;
			pos = position;
			nuc = base;
		}
		
		@SuppressWarnings("unused")
		protected String getChr() {
			return chr;
		}
		
		protected int getPos() {
			return pos;
		}
		
		protected String getBaseString() {
			char[] b = new char[1];
			b[0] = nuc;
			return new String(b);
		}
		
		protected char getBase() {
			return nuc;
		}
		
		@SuppressWarnings("unused")
		protected void setBase(char b) {
			nuc = b;
		}
		
	}
	
}
