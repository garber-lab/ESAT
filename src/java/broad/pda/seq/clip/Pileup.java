/**
 * 
 */
package broad.pda.seq.clip;

import broad.core.parser.CommandLineParser;
import broad.core.parser.StringParser;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

import org.apache.commons.collections15.Predicate;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import broad.core.math.EmpiricalDistribution;
import broad.core.math.MathUtil;
import broad.core.math.UnfairDie;
import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.pda.annotation.BEDFileParser;
import net.sf.samtools.util.CloseableIterator;
import nextgen.core.alignment.Alignment;
import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Annotation.Strand;
import nextgen.core.annotation.BasicAnnotation;
import nextgen.core.annotation.Gene;
import nextgen.core.coordinatesystem.TranscriptomeSpace;
import nextgen.core.model.AlignmentModel;
import nextgen.core.readFilters.FragmentLengthFilter;
import nextgen.core.readFilters.PairedAndProperFilter;
import nextgen.core.readers.WigReader;
import nextgen.core.utils.AnnotationUtils;

/**
 * @author prussell
 *
 */
public class Pileup {
	
	private AlignmentModel data;
	private Map<String, Collection<Gene>> annotation;
	private Map<Annotation, Double> overlapperCounts;
	private Map<Annotation, Annotation> peakRegions;
	private Map<String, Collection<Gene>> regions;
	private Map<Gene, Gene> regionToParent;
	static Logger logger = Logger.getLogger(Pileup.class.getName());
	private Map<String, Sequence> genome;
	private Map<Gene, Sequence> geneSequences;
	private WigReader wigFragmentEnds;
	private static int NUM_BINS = 200;
	
	/**
	 * @param bamFile Bam file of alignments
	 * @param annotationBedFile Bed transcriptome annotation
	 * @param genomeFasta Fasta file of genome
	 * @param readFilters Read filters to use
	 * @throws IOException
	 */
	public Pileup(String bamFile, String annotationBedFile, String subregionBedFile, String genomeFasta, Collection<Predicate<Alignment>> readFilters) throws IOException {
		logger.info("Loading gene annotation from " + annotationBedFile + "...");
		annotation = BEDFileParser.loadDataByChr(annotationBedFile);
		logger.info("Loading regions from " + subregionBedFile + "...");
		regions = BEDFileParser.loadDataByChr(subregionBedFile);
		logger.info("Mapping subregion to parent genes...");
		regionToParent = AnnotationUtils.mapChildToLargestParent(regions, annotation);
		logger.info("Loading alignment data from " + bamFile + "...");
		data = new AlignmentModel(bamFile, new TranscriptomeSpace(annotation), readFilters);
		logger.info("Loading genome from " + genomeFasta + "...");
		genome = FastaSequenceIO.getChrSequencesFromFasta(genomeFasta);
		geneSequences = new TreeMap<Gene, Sequence>();
		overlapperCounts = new TreeMap<Annotation, Double>();
		peakRegions = new TreeMap<Annotation, Annotation>();
		wigFragmentEnds = null;
		logger.info("Done constructing object.");
	}
	
	/**
	 * @param bamFile Bam file of alignments
	 * @param annotationBedFile Bed transcriptome annotation
	 * @param genomeFasta Fasta file of genome
	 * @param regionOverlapperCountFile File of overlapper counts for regions (can be incomplete)
	 * @param peakRegionBedFile Bed file of peak region (subregion of the parent gene enclosing all reads overlapping the region) corresponding to each region. Bed name is the same as region name.
	 * @param readFilters Read filters to use
	 * @throws IOException
	 */
	public Pileup(String bamFile, String fragmentEndWig, String annotationBedFile, String subregionBedFile, String genomeFasta, String regionOverlapperCountFile, String peakRegionBedFile, Collection<Predicate<Alignment>> readFilters) throws IOException {
		this(bamFile, annotationBedFile, subregionBedFile, genomeFasta, readFilters);
		if(fragmentEndWig != null) {
			logger.info("Reading fragment end pileups from file " + fragmentEndWig);
			wigFragmentEnds = new WigReader(fragmentEndWig);
		}
		if(regionOverlapperCountFile != null) {
			overlapperCounts.putAll(readCountsFromFile(regionOverlapperCountFile));
		}
		if(peakRegionBedFile != null) {
			peakRegions.putAll(readPeaksPerRegionFromFile(peakRegionBedFile));
		}
	}
	
	/**
	 * Get a region by name from stored regions
	 * @param regionName Region name
	 * @return The annotation object for the region
	 */
	private Annotation regionByName(String regionName) {
		for(String chr : regions.keySet()) {
			for(Annotation region : regions.get(chr)) {
				if(region.getName().equals(regionName)) {
					return region;
				}
			}
		}
		throw new IllegalArgumentException("Region " + regionName + " not recognized.");
	}
	
	/**
	 * Read peak region (span of all reads overlapping region) from a bed file where the bed name matches the original region name
	 * @param file Bed file
	 * @return Map of original region to peak region
	 * @throws IOException
	 */
	private Map<Annotation, Annotation> readPeaksPerRegionFromFile(String file) throws IOException {
		logger.info("Reading peak regions from file " + file);
		Map<String, Gene> byName = BEDFileParser.loadDataByName(new File(file));
		Map<Annotation, Annotation> rtrn = new TreeMap<Annotation, Annotation>();
		for(String name : byName.keySet()) {
			rtrn.put(regionByName(name), byName.get(name));
		}
		return rtrn;
	}
	
	/**
	 * Read overlapper counts for regions from a table file (format: region_name\tcount) instead of getting from bam file
	 * @param file Table file
	 * @return Map of region to overlapper count
	 * @throws IOException
	 */
	private Map<Annotation, Double> readCountsFromFile(String file) throws IOException {
		logger.info("Reading region counts from file " + file);
		Map<Annotation, Double> rtrn = new TreeMap<Annotation, Double>();
		FileReader r = new FileReader(file);
		BufferedReader b = new BufferedReader(r);
		StringParser s = new StringParser();
		while(b.ready()) {
			s.parse(b.readLine());
			if(s.getFieldCount() == 0) continue;
			if(s.getFieldCount() != 2) {
				r.close();
				b.close();
				throw new IllegalArgumentException("Line format: region_name\tcount");
			}
			rtrn.put(regionByName(s.asString(0)), Double.valueOf(s.asDouble(1)));
		}
		r.close();
		b.close();
		return rtrn;
	}

	
	/**
	 * Get reads overlapping a region
	 * @param region The region
	 * @param fullyContained Include fully contained reads only
	 * @return Iterator over reads overlapping the region
	 */
	private CloseableIterator<Alignment> getPileup(Annotation region, boolean fullyContained) {
		logger.debug("Getting pileup over region " + region.getName());
		return data.getOverlappingAnnotations(region, fullyContained);
	}
	
	/**
	 * Get count of fragments overlapping a region from the wig file instead of bam file
	 * @param region Region
	 * @return Sum of all wig values in the exons of the region
	 */
	private int getCountFromWig(Annotation region) {
		
		double rtrn = 0;
		Map<Integer, Double> vals = wigFragmentEnds.getValues(region);
		for(Integer pos : vals.keySet()) {
			rtrn += vals.get(pos).doubleValue();
		}
		
		return (int) rtrn;
		
	}
	
	/**
	 * Count the number of reads overlapping the region
	 * @param region Region
	 * @param fullyContained Include fully contained reads only
	 * @return Count of reads overlapping the region
	 */
	private int getCountFromBam(Annotation region, boolean fullyContained) {
		
		if(!overlapperCounts.containsKey(region)) {
			double d = data.getCount(region, fullyContained);
			overlapperCounts.put(region, Double.valueOf(d));
		}
		
		return (int) overlapperCounts.get(region).doubleValue();
		
	}
	
	/**
	 * Getting peak region (span of all overlapping reads) from the cache or compute from bam file
	 * @param region The region
	 * @param fullyContained Include fully contained reads only
	 * @return The peak region
	 */
	private Annotation getPeak(Annotation region, boolean fullyContained) {
		
		if(!peakRegions.containsKey(region)) {
			logger.debug("Getting peak region for " + region.toUCSC());
			peakRegions.put(region, data.getPeak(regionToParent.get(region), region, fullyContained));
		}
		
		return peakRegions.get(region);
		
	}
	
	/**
	 * Get sequence of a gene using a cache
	 * @param gene The gene
	 * @return The nucleotide sequence of the gene
	 */
	private Sequence getGeneSequence(Gene gene) {
		if(!geneSequences.containsKey(gene)) {
			//logger.debug("Getting sequence of gene " + gene.getName());
			Sequence seq = genome.get(gene.getChr()).getSubsequence(gene);
			geneSequences.put(gene, seq);
		}
		return geneSequences.get(gene);
	}
	
	/**
	 * Check that a nucleotide probability map is valid
	 * @param nucProbs Map of nucleotide to background probability
	 */
	private static void validateNucProbabilityMap(Map<String, Double> nucProbs) {
		double totalProb = 0;
		for(String nuc : nucProbs.keySet()) {
			if(nuc.length() != 1) {
				throw new IllegalArgumentException("Not a valid nucleotide: " + nuc);
			}
			if(!Character.isUpperCase(nuc.charAt(0))) {
				throw new IllegalArgumentException("All nucleotides in map must be uppercase. Queries are converted to uppercase.");
			}
			double prob = nucProbs.get(nuc).doubleValue();
			if(prob < 0 || prob > 1) {
				throw new IllegalArgumentException("Not a valid probability: " + prob);
			}
			totalProb += prob;
		}
		if(!MathUtil.closeTo1(totalProb)) {
			throw new IllegalArgumentException("Probabilities must add up to 1. Sum: " + totalProb);
		}
	}
	
	/**
	 * Get string of nucleotide probabilities for printing
	 * @param nucProbs Map of nucleotide to probability
	 * @return String representation for printing
	 */
	@SuppressWarnings("unused")
	private static String nucleotideProbString(Map<String, Double> nucProbs) {
		String nucProbString = "";
		for(String nuc : nucProbs.keySet()) {
			nucProbString += nuc + ":" + nucProbs.get(nuc) + "\t";
		}
		return nucProbString;
	}
	
	/**
	 * Get the probability of randomly choosing each position in the sequence, biased by background nucleotide probabilities
	 * @param nucProbs Background nucleotide probabilities
	 * @param sequence The sequence
	 * @return Probability of each position in sequence; adds up to 1
	 */
	private static double[] nucleotideBiasedPositionProbabilities(Map<String, Double> nucProbs, Sequence sequence) {
	
		//logger.debug("Getting biased position probabilities for sequence " + sequence.getSequenceBases());
		//logger.debug("Nucleotide probabilities\t" + nucleotideProbString(nucProbs));
		
		validateNucProbabilityMap(nucProbs);
		String bases = sequence.getSequenceBases();
		int seqLen = bases.length();
		
		// Count the number of occurrences of each nucleotide in the sequence
		// Only count nucleotides that actually appear
		Map<String, Integer> numOccurrences = new TreeMap<String, Integer>();
		for(int i = 0; i < seqLen; i++) {
			String nuc = bases.substring(i, i+1).toUpperCase();
			if(!nucProbs.containsKey(nuc)) {
				throw new IllegalArgumentException("Nucleotide probability map does not contain value for nucleotide " + nuc + ".");
			}
			if(!numOccurrences.containsKey(nuc)) {
				numOccurrences.put(nuc, Integer.valueOf(1));
				continue;
			} else {
				numOccurrences.put(nuc, Integer.valueOf(numOccurrences.get(nuc).intValue() + 1));
			}
		}
		
		// Get the scale factor to divide individual nucleotide probabilities by to make the whole sequence add up to 1
		double scaleFactor = 0;
		for(String nuc : numOccurrences.keySet()) {
			double d = nucProbs.get(nuc).doubleValue() * numOccurrences.get(nuc).doubleValue();
			scaleFactor += d;
		}
		
		// The probability of each position is the background probability of that nucleotide divided by the scale factor
		double[] rtrn = new double[seqLen];
		for(int i = 0; i < seqLen; i++) {
			String nuc = bases.substring(i, i+1).toUpperCase();
			rtrn[i] = nucProbs.get(nuc).doubleValue() / scaleFactor;
		}
		
		// Make sure adds up to 1
		double rtrnTotal = 0;
		for(int i = 0; i < rtrn.length; i++) {
			rtrnTotal += rtrn[i];
		}
		if(!MathUtil.closeTo1(rtrnTotal)) {
			throw new IllegalStateException("Assigned probabilities add up to " + rtrnTotal + ".");
		}
				
		return rtrn;
		
	}
	
	/**
	 * Get string of array values for printing
	 * @param vals Array
	 * @return String of values on single line
	 */
	@SuppressWarnings("unused")
	private static String arrayString(double[] vals) {
		String rtrnString = "";
		for(int i = 0; i < vals.length; i++) {
			rtrnString += vals[i] + "\t";
		}
		return rtrnString;
	}
	
	/**
	 * Read nucleotide probabilities from a file. Include uppercase nucleotide letters only. Format: nucleotide\tprobability
	 * @param file File to read from
	 * @return Map of nucleotide to probability
	 * @throws IOException
	 */
	@SuppressWarnings("unused")
	private static Map<String, Double> readNucProbsFromFile(String file) throws IOException{
		return readNucProbsFromFile(file, false);
	}
	
	/**
	 * Read nucleotide probabilities from a file. Include uppercase nucleotide letters only. Format: nucleotide\tprobability
	 * @param file File to read from
	 * @param scaleProbsToOne If the probabilities in the file do not add up to 1, scale them so they do
	 * @return Map of nucleotide to probability
	 * @throws IOException
	 */
	private static Map<String, Double> readNucProbsFromFile(String file, boolean scaleProbsToOne) throws IOException{
		FileReader r = new FileReader(file);
		BufferedReader b = new BufferedReader(r);
		StringParser s = new StringParser();
		Map<String, Double> rtrn = new TreeMap<String, Double>();
		while(b.ready()) {
			String line = b.readLine();
			s.parse(line);
			if(s.getFieldCount() == 0) {
				continue;
			}
			if(s.getFieldCount() != 2) {
				r.close();
				b.close();
				throw new IllegalArgumentException("Line format: nucleotide\tprob");
			}
			rtrn.put(s.asString(0), Double.valueOf(s.asDouble(1)));
		}
		r.close();
		b.close();
		
		if(scaleProbsToOne) {
			double total = 0;
			for(Double d : rtrn.values()) {
				total += d.doubleValue();
			}
			double factor = 1 / total;
			for(String str : rtrn.keySet()) {
				double scaled = rtrn.get(str).doubleValue() * factor;
				rtrn.put(str, Double.valueOf(scaled));
			}
		}
		
		validateNucProbabilityMap(rtrn);
		return rtrn;
	}
	
	/**
	 * For each region get several shufflings of the positions
	 * @param fullyContained When computing number of fragments and peak region, include fully contained fragments only
	 * @param nucProbs Background nucleotide probabilities
	 * @param numPermutations Number of permutations per region
	 * @return Map of region to collection of random permutations
	 */
	public Map<Annotation, Collection<EmpiricalDistribution>> getShuffledPositionsAllRegions(boolean fullyContained, Map<String, Double> nucProbs, int numPermutations) {
		
		logger.info("Getting " + numPermutations + " shuffled sets of positions for each region");
		Map<Annotation, Collection<EmpiricalDistribution>> rtrn = new TreeMap<Annotation, Collection<EmpiricalDistribution>>();
		for(String chr : regions.keySet()) {
			logger.info(chr);
			for(Annotation region : regions.get(chr)) {
				rtrn.put(region, getShuffledPositions(region, fullyContained, nucProbs, numPermutations));
			}
		}
		
		return rtrn;
		
	}
	
	/**
	 * For each region, find the highest fragment end pileup from the bam file.
	 * Randomly permute the fragments in a biased way according to background nucleotide probabilities.
	 * Get a scan P value for the real highest pileup in the region.
	 * @param fullyContained When computing number of fragments and peak region, include fully contained fragments only
	 * @param nucProbs Background nucleotide probabilities
	 * @param numPermutations Number of permutations per region
	 * @param fivePrimeEnd Use the 5' end of fragments (if false, use 3' end)
	 * @return Map of region to P value
	 */
	public Map<Annotation, Double> getScanPvalFragmentEndsMaxPositionAllRegionsBam(boolean fullyContained, Map<String, Double> nucProbs, int numPermutations, boolean fivePrimeEnd) {
		logger.info("Getting scan p-values of fragment end pileups for all regions");
		Map<Annotation, Double> rtrn = new TreeMap<Annotation, Double>();
		for(String chr : regions.keySet()) {
			for(Annotation region : regions.get(chr)) {
				rtrn.put(region, Double.valueOf(getScanPvalFragmentEndsMaxPositionBam(region, fullyContained, nucProbs, numPermutations, fivePrimeEnd)));
			}
		}
		return rtrn;
	}
	
	/**
	 * Randomly permute the fragments in a biased way according to background nucleotide probabilities.
	 * Get a scan P value for pileup height of each position in the region as read from wig file.
	 * @param region The region
	 * @param fullyContained If computing peak region from the bam file, include fully contained fragments only
	 * @param nucProbs Background nucleotide probabilities
	 * @param numPermutations Number of permutations
	 * @return Scan P value for each position in region
	 */
	public Map<Integer, Double> getScanPvalsFragmentEndPileupsWig(Annotation region, boolean fullyContained, Map<String, Double> nucProbs, int numPermutations) {
				
		// Make empirical scan distribution for shuffled fragments
		Collection<EmpiricalDistribution> shuffledFragmentEnds = getShuffledPositions(region, fullyContained, nucProbs, numPermutations);
		// TODO add option for pileup height / number of overlappers
		Collection<EmpiricalDistribution> shuffledPileupHeights = new ArrayList<EmpiricalDistribution>();
		for(EmpiricalDistribution shuffled : shuffledFragmentEnds) {
			EmpiricalDistribution d = shuffled.getDistributionOfHistogramValues();
			logger.debug("Max pileup for shuffled fragment ends is " + d.getMax());
			shuffledPileupHeights.add(d);
		}
		EmpiricalDistribution shuffleScanDistribution = EmpiricalDistribution.getEmpiricalMaxStatisticDistribution(shuffledPileupHeights);

		// Get real pileup heights from wig file
		Map<Integer, Double> rtrn = new TreeMap<Integer, Double>();
		Map<Integer, Double> realPileupHeights = wigFragmentEnds.getValues(region);
		for(Integer pos : realPileupHeights.keySet()) {
			Double val = realPileupHeights.get(pos);
			double pval = shuffleScanDistribution.getPValue(val.doubleValue());
			logger.debug("P value for " + region.getChr() + " " + pos.toString() + " (height " + val.toString() + ") is " + pval);
			rtrn.put(pos, Double.valueOf(pval));
		}
	
		return rtrn;
		
	}

	/**
	 * Randomly permute the fragments in a biased way according to background nucleotide probabilities.
	 * For each region get a scan P value for pileup height of each position in the region as read from wig file.
	 * @param fullyContained If computing peak regions from the bam file, include fully contained fragments only
	 * @param nucProbs Background nucleotide probabilities
	 * @param numPermutations Number of permutations per region
	 * @return Map of region to scan P value for each position in region
	 */
	public Map<Annotation, Map<Integer, Double>> getScanPvalsFragmentEndPileupsAllRegionsWig(boolean fullyContained, Map<String, Double> nucProbs, int numPermutations) {
		logger.info("Getting scan p-values of fragment end pileups for all regions");
		Map<Annotation, Map<Integer, Double>> rtrn = new TreeMap<Annotation, Map<Integer, Double>>();
		for(String chr : regions.keySet()) {
			for(Annotation region : regions.get(chr)) {
				logger.info(region.getName());
				rtrn.put(region, getScanPvalsFragmentEndPileupsWig(region, fullyContained, nucProbs, numPermutations));
				for(Integer pos : rtrn.get(region).keySet()) {
					logger.debug(region.toUCSC() + "\t" + pos.toString() + "\t" + rtrn.get(region).get(pos));
				}
			}
		}
		return rtrn;
	}

	/**
	 * Randomly permute the fragments in a biased way according to background nucleotide probabilities.
	 * For each region get a scan P value for pileup height of each position in the region as read from wig file.
	 * Write information to a table.
	 * @param outTable Output table file
	 * @param outBed Bed file for significant positions
	 * @param pvalCutoffForBed P value cutoff for significant position
	 * @param fullyContained If computing peak regions from the bam file, include fully contained fragments only
	 * @param nucProbs Background nucleotide probabilities
	 * @param numPermutations Number of permutations per region
	 * @return Map of region to scan P value for each position in region
	 */
	public void writeScanPvalsFragmentEndPileupsAllRegionsWig(String outTable, String outBed, double pvalCutoffForBed, boolean fullyContained, Map<String, Double> nucProbs, int numPermutations) throws IOException {
		
		logger.info("Writing P values based on wig file to table " + outTable);
		logger.info("Writing positions with pval < " + pvalCutoffForBed + " to file " + outBed);
		FileWriter wt = new FileWriter(outTable);
		FileWriter wb = new FileWriter(outBed);
		
		String header = "parent_gene\t";
		header += "region_name\t";
		header += "region_coord\t";
		header += "pos\t";
		header += "pos_wig_value\t";
		header += "pval\t";
		wt.write(header + "\n");
		
		Map<Annotation, Map<Integer, Double>> pvals = getScanPvalsFragmentEndPileupsAllRegionsWig(fullyContained, nucProbs, numPermutations);
		for(Annotation region : pvals.keySet()) {
			for(Integer pos : pvals.get(region).keySet()) {
				String line = regionToParent.get(region).getName() + "\t";
				line += region.getName() + "\t";
				line += region.toUCSC() + "\t";
				line += pos.toString() + "\t";
				line += wigFragmentEnds.getValue(region.getChr(), pos.intValue()) + "\t";
				double pval = pvals.get(region).get(pos).doubleValue();
				line += pval + "\t";
				// If significant write to bed file
				if(pval < pvalCutoffForBed) {
					BasicAnnotation ba = new BasicAnnotation(region.getChr(), pos.intValue(), pos.intValue(), region.getOrientation(), region.toUCSC() + ":" + pos.toString() + ":" + pval);
					ba.setScore(pval);
					wb.write(ba.toBED(255,0,0) + "\n");
				}
				wt.write(line + "\n");
			}
		}
		
		wt.close();
		wb.close();
		logger.info("Done writing files.");
		
	}
	
	/**
	 * Find the highest fragment end pileup in the region from the bam file.
	 * Randomly permute the fragments in a biased way according to background nucleotide probabilities.
	 * Get a scan P value for the real highest pileup in the region.
	 * @param region The region
	 * @param fullyContained When computing number of fragments and peak region, include fully contained fragments only
	 * @param nucProbs Background nucleotide probabilities
	 * @param numPermutations Number of permutations per region
	 * @param fivePrimeEnd Use the 5' end of fragments (if false, use 3' end)
	 * @return Map of region to P value
	 */
	public double getScanPvalFragmentEndsMaxPositionBam(Annotation region, boolean fullyContained, Map<String, Double> nucProbs, int numPermutations, boolean fivePrimeEnd) {
		
		EmpiricalDistribution realFragmentEnds = getFragmentEndsFromBamFile(region, fullyContained, fivePrimeEnd);
		EmpiricalDistribution pileupHeights = realFragmentEnds.getDistributionOfHistogramValues();
		double maxRealPileup = pileupHeights.getMax();
		logger.debug("Max real pileup for region " + region.toUCSC() + " is " + maxRealPileup + ".");
		
		Collection<EmpiricalDistribution> shuffledFragmentEnds = getShuffledPositions(region, fullyContained, nucProbs, numPermutations);
		Collection<EmpiricalDistribution> shuffledPileupHeights = new ArrayList<EmpiricalDistribution>();
		for(EmpiricalDistribution shuffled : shuffledFragmentEnds) {
			if(realFragmentEnds.numBins() != shuffled.numBins()) {
				throw new IllegalArgumentException("All distributions must have same number of bins");
			}
			EmpiricalDistribution d = shuffled.getDistributionOfHistogramValues();
			logger.debug("Max pileup for shuffled fragment ends is " + d.getMax());
			shuffledPileupHeights.add(d);
		}
		EmpiricalDistribution shuffleScanDistribution = EmpiricalDistribution.getEmpiricalMaxStatisticDistribution(shuffledPileupHeights);
		
		double pVal = shuffleScanDistribution.getPValue(maxRealPileup);
		logger.debug("P value is " + pVal);
		
		return pVal;
		
	}
		
	
	/**
	 * Get shuffled sets of positions within the span covered by fragments overlapping the region
	 * Each set contains a number of positions equal to the number of fragments overlapping the region
	 * Random shuffled positions are biased by nucleotide
	 * @param region Region within parent gene
	 * @param fullyContained When getting the span covered by fragments overlapping the region, include fully contained fragments only
	 * @param nucProbs Background probabilities of each nucleotide; must add up to 1
	 * @param numPermutations Number of random sets to get
	 * @return Collection of random position sets; each set is presented in an empirical distribution
	 */
	public Collection<EmpiricalDistribution> getShuffledPositions(Annotation region, boolean fullyContained, Map<String, Double> nucProbs, int numPermutations) {
		
		if(!regionToParent.containsKey(region)) {
			throw new IllegalArgumentException("Region " + region.getName() + " not in stored regions");
		}
		Gene parent = regionToParent.get(region);
		
		//logger.debug("");
		logger.debug("Getting " + numPermutations + " sets of permuted positions for parent gene " + parent.getName() + " and region " + region.getName());
		//logger.debug("Nucleotide probabilities\t" + nucleotideProbString(nucProbs));

		// Count the number of fragments overlapping the region
		int numOverlappers = -1;
		if(wigFragmentEnds != null) {
			numOverlappers = getCountFromWig(region);
		} else {
			numOverlappers = getCountFromBam(region, fullyContained);
		}
		logger.debug("Region has " + numOverlappers + " overlappers.");
		
		// Get parent gene sequence from genome
		Sequence parentSequence = getGeneSequence(parent);
		
		// Get span covered by overlappers
		//logger.debug("Getting peak region");
		Annotation span = getPeak(region, fullyContained);
		logger.debug("Parent\t" + parent.getName() + "\tregion\t" + region.getName() + "\tpeak\t" + span.toUCSC());
		
		// First position of parent covered by a read overlapping the transcript, in transcript coordinates
		int start = parent.genomicToTranscriptPosition(span.getStart());
		// Last position of parent covered by a read overlapping the transcript, in transcript coordinates
		int end = parent.genomicToTranscriptPosition(span.getEnd() - 1);
		int regionPeakFirstPos = Math.min(start, end);
		int regionPeakLastPos = Math.max(start, end);
		//logger.debug("Peak coordinates on parent transcript: " + regionPeakFirstPos + "-" + regionPeakLastPos);
		
		// Get subsequence of gene comprising the covered region
		Sequence regionPeakSubsequence = parentSequence.getSubSequence(parentSequence.getId() + ":" + regionPeakFirstPos + "-" + Integer.valueOf(regionPeakLastPos+1).toString(), regionPeakFirstPos, regionPeakLastPos + 1);
		//logger.debug("Peak sequence " + regionPeakSubsequence.getSequenceBases());
		
		// Get position probabilities for subsequence
		double[] positionProbs = nucleotideBiasedPositionProbabilities(nucProbs, regionPeakSubsequence);
		//logger.debug("Position probabilities\t" + arrayString(positionProbs));
		
		// Make the shuffled positions
		Collection<EmpiricalDistribution> rtrn = new ArrayList<EmpiricalDistribution>();
		UnfairDie unfairDie = new UnfairDie(positionProbs);
		for(int i = 0; i < numPermutations; i++) {
			double[] shuffledPositions = new double[numOverlappers];
			for(int j = 0; j < numOverlappers; j++) {
				shuffledPositions[j] = unfairDie.roll();
			}
			EmpiricalDistribution dist = new EmpiricalDistribution(shuffledPositions, NUM_BINS);
			Collection<Double> vals = dist.getAllDataValues();
			String valsString = "";
			for(Double val : vals) {
				valsString += val.toString() + "\t";
			}
			logger.debug("Random permutation " + i + ":\t" + valsString);
			rtrn.add(dist);
		}
		return rtrn;
		
	}
		
	/**
	 * Get an empirical distribution storing the fragment ends of all fragments overlapping the region
	 * Fragment ends are in transcript coordinates with respect to the parent transcript
	 * @param parent Parent gene
	 * @param region Region within parent to get pileup over
	 * @param fullyContained Only include fragments that are fully contained in the subregion
	 * @param fivePrimeEnd Report the 5' end of fragments (if false, reports 3' end)
	 * @return Distribution of fragment endpoints in transcript coordinates
	 */
	public EmpiricalDistribution getFragmentEndsFromBamFile(Annotation region, boolean fullyContained, boolean fivePrimeEnd) {
		Gene parent = regionToParent.get(region);
		//logger.debug("");
		//logger.debug("Getting fragment ends for region " + region.getName() + " in parent " + parent.getName() + ".");
		CloseableIterator<Alignment> overlappers = getPileup(region, fullyContained);
		Collection<Double> endpoints = new ArrayList<Double>();
		Strand parentStrand = parent.getOrientation();
		if(parentStrand.equals(Strand.UNKNOWN)) {
			throw new IllegalArgumentException("Parent gene orientation must be known");
		}
		while(overlappers.hasNext()) {
			Alignment alignment = overlappers.next();
			int endpointGenomicCoord = -1;
			int firstGenomicPosition = alignment.getAlignmentStart();
			int lastGenomicPosition = alignment.getAlignmentEnd() - 1;
			@SuppressWarnings("unused")
			String message = "";
			if(fivePrimeEnd) {
				endpointGenomicCoord = parentStrand.equals(Strand.POSITIVE) ? firstGenomicPosition : lastGenomicPosition;
				message += "Fragment " + alignment.getName() + " 5' end genomic coordinate is " + endpointGenomicCoord + "."; 
			} else {
				endpointGenomicCoord = parentStrand.equals(Strand.NEGATIVE) ? firstGenomicPosition : lastGenomicPosition;
				message += "Fragment " + alignment.getName() + " 3' end genomic coordinate is " + endpointGenomicCoord + "."; 
			}
			Double endpointTranscriptCoord = Double.valueOf(parent.genomicToTranscriptPosition(endpointGenomicCoord));
			if(endpointTranscriptCoord.doubleValue() < 0) {
				String m = "";
				if(fivePrimeEnd) m += "Skipping 5' end of fragment " + alignment.getName() + " because end is not within mature transcript";
				else m += "Skipping 5' end of fragment " + alignment.getName() + " because end is not within mature transcript";
				logger.warn(m);
				continue;
			}
			//logger.debug(message + " Transcript coordinate is " + endpointTranscriptCoord.toString() + ".");
			endpoints.add(endpointTranscriptCoord);
		}
		return new EmpiricalDistribution(endpoints, NUM_BINS);
	}
	
	/**
	 * Get the transcriptome annotation
	 * @return The annotation
	 */
	public Map<String, Collection<Gene>> getAnnotation() {
		return annotation;
	}
	
	private static void writeFragmentEndStats(Pileup pileup, String regionBed, String outFile, boolean fullyContained) throws IOException {
		logger.info("Writing information about pileups over regions in " + regionBed + " to file " + outFile + "...");
		Map<String, Collection<Gene>> regionsByChr = BEDFileParser.loadDataByChr(regionBed);
		Map<Gene, Gene> regionToParentAnnotation = AnnotationUtils.mapChildToLargestParent(regionsByChr, pileup.getAnnotation());
		FileWriter w = new FileWriter(outFile);
		String header = "Region\t";
		header += "Parent_gene\t";
		header += "Num_overlapping_fragments\t";
		header += "Fragment_5prime_end_median_pos\t";
		header += "Fragment_3prime_end_median_pos\t";
		w.write(header + "\n");
		for(Gene region : regionToParentAnnotation.keySet()) {
			Gene parent = regionToParentAnnotation.get(region);
			EmpiricalDistribution fragment5pEnds = pileup.getFragmentEndsFromBamFile(region, fullyContained, true);
			EmpiricalDistribution fragment3pEnds = pileup.getFragmentEndsFromBamFile(region, fullyContained, false);
			int numOverlappingFragments = fragment5pEnds.getTotalObservations();
			double median5p = fragment5pEnds.getMedianOfAllDataValues();
			double median3p = fragment3pEnds.getMedianOfAllDataValues();
			String line = region.getName() + "\t";
			line += parent.getName() + "\t";
			line += numOverlappingFragments + "\t";
			line += median5p + "\t";
			line += median3p + "\t";
			w.write(line + "\n");
		}
		w.close();
		logger.info("Done writing file.");
	}

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-b", "Bam file", true);
		p.addStringArg("-a", "Bed annotation for transcriptome space", true);
		p.addStringArg("-r", "Bed file of regions", true);
		p.addBooleanArg("-pp", "When using bam file, include paired reads and proper pairs only", false, true);
		p.addIntArg("-f", "When using bam file, max fragment length", false, 1000);
		p.addBooleanArg("-d", "Debug logging on", false, false);
		p.addStringArg("-oe", "Output table for fragment end stats", false, null);
		p.addStringArg("-g", "Genome fasta", true);
		p.addStringArg("-np", "Nucleotide probability file", false, null);
		p.addBooleanArg("-sp", "When providing nucleotide probability file, scale all probabilitie to add up to one if they don't already", false, false);
		p.addIntArg("-nr", "Number of random position permutations per region", false, 100);
		p.addBooleanArg("-fc", "When using bam file to get fragments overlapping regions,  use fully contained fragments only", false, true);
		p.addStringArg("-oc", "Region overlapper count file to use instead of bam file or wig file (format: region_name\tcount", false, null);
		p.addStringArg("-pb", "When using bam file, bed file of region to use for randomized positions per region. If not provided or incomplete for some regions, program will compute from bam file and use subregion of the parent gene enclosing all reads overlapping the region. Bed name is the same as region name.", false, null);
		p.addStringArg("-wf", "Wig file of fragment end pileups to use instead of getting overlappers from bam file", false, null);
		p.addStringArg("-otw", "Output table for position P values based on wig counts (requires -np and -obs)", false, null);
		p.addStringArg("-obs", "Output bed file of significant pileup positions (requires -np and -otw)", false, null);
		p.addDoubleArg("-p", "P value cutoff for significant pileup", false, 0.001);
		p.parse(args);
		if(p.getBooleanArg("-d")) {
			logger.setLevel(Level.DEBUG);
		}
		
		String genomeFasta = p.getStringArg("-g");
		String bamFile = p.getStringArg("-b");
		String bedFile = p.getStringArg("-a");
		Collection<Predicate<Alignment>> readFilters = new ArrayList<Predicate<Alignment>>();
		if(p.getBooleanArg("-pp")) {
			readFilters.add(new PairedAndProperFilter());
		}
		int maxFragmentLength = p.getIntArg("-f");
		String regionBed = p.getStringArg("-r");
		String outputFragmentEndInfo = p.getStringArg("-oe");
		String nucProbFile = p.getStringArg("-np");
		boolean scaleProbsToOne = p.getBooleanArg("-sp");
		int numRand = p.getIntArg("-nr");
		boolean fullyContained = p.getBooleanArg("-fc");
		String overlapperCountFile = p.getStringArg("-oc");
		String peakRegionBed = p.getStringArg("-pb");
		String wigFile = p.getStringArg("-wf");
		String outTableWigPval = p.getStringArg("-otw");
		String outBedSigPos = p.getStringArg("-obs");
		double pvalCutoff = p.getDoubleArg("-p");
		
		Pileup pileup = new Pileup(bamFile, wigFile, bedFile, regionBed, genomeFasta, overlapperCountFile, peakRegionBed, readFilters);
		pileup.data.addFilter(new FragmentLengthFilter(pileup.data.getCoordinateSpace(), maxFragmentLength));
		
		if(outputFragmentEndInfo != null) {
			writeFragmentEndStats(pileup, regionBed, outputFragmentEndInfo, false);
		}
		
		if(nucProbFile != null && outTableWigPval != null && outBedSigPos != null) {
			pileup.writeScanPvalsFragmentEndPileupsAllRegionsWig(outTableWigPval, outBedSigPos, pvalCutoff, fullyContained, readNucProbsFromFile(nucProbFile, scaleProbsToOne), numRand);
		}
		
		logger.info("");
		logger.info("All done.");
		
	}

}
