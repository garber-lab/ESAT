/**
 * 
 */
package broad.core.overlaputils;

import broad.core.parser.CommandLineParser;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;

import nextgen.core.annotation.Gene;
import nextgen.core.utils.AnnotationUtils;

import broad.pda.annotation.BEDFileParser;



/**
 * @author prussell
 *
 */
public class GeneSetIntersect {

	private static Logger logger = Logger.getLogger(GeneSetIntersect.class.getName());
	private static DecimalFormat decimalFormat = new DecimalFormat("#.##");

	
	/**
	 * Map each gene to the collection of overlappers in another set
	 * @param genes The genes by chromosome
	 * @param intersectSet The other set by chromosome
	 * @return By chromosome: map of each gene that has overlappers in the other set to its overlappers
	 */
	public static Map<String, Map<Gene, Collection<Gene>>> mapGenesToOverlappers(Map<String, Collection<Gene>> genes, Map<String, Collection<Gene>> intersectSet) {
		Map<String, Map<Gene, Collection<Gene>>> rtrn = new TreeMap<String, Map<Gene, Collection<Gene>>>();
		for(String chr : genes.keySet()) {
			if(!intersectSet.containsKey(chr)) {
				continue;
			}
			Map<Gene, Collection<Gene>> overlappers = new TreeMap<Gene, Collection<Gene>>();
			for(Gene gene : genes.get(chr)) {
				Collection<Gene> overlappersThisGene = new TreeSet<Gene>();
				for(Gene other : intersectSet.get(chr)) {
					if(other.getEnd() < gene.getStart() || other.getStart() > gene.getEnd()) {
						continue;
					}
					if(gene.overlaps(other)) {
						overlappersThisGene.add(other);
					}
				}
				if(!overlappersThisGene.isEmpty()) {
					overlappers.put(gene, overlappersThisGene);
				}
			}
			if(!overlappers.isEmpty()) {
				rtrn.put(chr, overlappers);
			}
		}
		return rtrn;
	}
	
	/**
	 * Get overlap between two gene sets
	 * @param genes Genes by chromosome
	 * @param intersectSet Set to intersect with by chromosome
	 * @return Set of Gene objects representing the overlaps
	 */
	public static Map<String, Collection<Gene>> getOverlap(Map<String, Collection<Gene>> genes, Map<String, Collection<Gene>> intersectSet) {
		return getOverlap(genes, intersectSet, false);
	}
	
	/**
	 * Get overlap between two gene sets
	 * @param genes Genes by chromosome
	 * @param intersectSet Set to intersect with by chromosome
	 * @param verbose Print messages
	 * @return Set of Gene objects representing the overlaps
	 */
	public static Map<String, Collection<Gene>> getOverlap(Map<String, Collection<Gene>> genes, Map<String, Collection<Gene>> intersectSet, boolean verbose) {
		Map<String, Collection<Gene>> excludeGenes = new TreeMap<String, Collection<Gene>>();
		return getOverlap(genes, intersectSet, excludeGenes, verbose);
	}
	
	/**
	 * Get overlap between two gene sets and automatically exclude all genes in a set
	 * @param genes Genes by chromosome
	 * @param intersectSet Set to intersect with by chromosome
	 * @param excludeSet Exclude genes if they overlap a region from this set
	 * @return Set of Gene objects representing the overlaps
	 */
	public static Map<String, Collection<Gene>> getOverlap(Map<String, Collection<Gene>> genes, Map<String, Collection<Gene>> intersectSet, Map<String, Collection<Gene>> excludeSet) {
		return getOverlap(genes, intersectSet, excludeSet, false);
	}
	
	/**
	 * Get overlap between two gene sets and automatically exclude all genes in a set
	 * @param genes Genes by chromosome
	 * @param intersectSet Set to intersect with by chromosome
	 * @param excludeSet Exclude genes from first gene set if they overlap a gene from this set
	 * @param verbose Print messages
	 * @return Set of Gene objects representing the overlaps
	 */
	public static Map<String, Collection<Gene>> getOverlap(Map<String, Collection<Gene>> genes, Map<String, Collection<Gene>> intersectSet, Map<String, Collection<Gene>> excludeSet, boolean verbose) {
		if(excludeSet.isEmpty()) {
			for(String chr : genes.keySet()) {
				Collection<Gene> empty = new TreeSet<Gene>();
				excludeSet.put(chr, empty);
			}
		}
		Map<String, Collection<Gene>> rtrn = new TreeMap<String, Collection<Gene>>();
		for(String chr : genes.keySet()) {
			if(verbose) {
				logger.info(chr);
			}
			Collection<Gene> thisChrOverlaps = new TreeSet<Gene>();
			for(Gene gene : genes.get(chr)) {
				for(Gene other : excludeSet.get(chr)) {
					if(gene.overlaps(other)) continue;
				}
				if(!intersectSet.containsKey(chr)) continue;
				Gene overlap = gene.getOverlap(intersectSet.get(chr));
				if(overlap == null) continue;
				overlap.setName(gene.getName());
				thisChrOverlaps.add(overlap);
			}
			rtrn.put(chr, thisChrOverlaps);
		}
		return rtrn;
	}
	
	
	/**
	 * Count number of regions that overlap at least one region in other set
	 * @param regionsToCount Regions to count
	 * @param regionsToOverlapWith Potential overlappers
	 * @return The number of regions in the first set that overlap at least one region in other set
	 */
	public static int countOverlappers(Map<String, Collection<Gene>> regionsToCount, Map<String, Collection<Gene>> regionsToOverlapWith) {
		return totalGenes(getOverlap(regionsToCount, regionsToOverlapWith));
	}
	
	/**
	 * Count number of regions that overlap at least one region in other set
	 * @param regionsToCount Regions to count
	 * @param regionsToOverlapWith Potential overlappers
	 * @param excludeSet Exclude regions from first set if they overlap a region in this set
	 * @return The number of regions in the first set that overlap at least one region in other set
	 */
	public static int countOverlappers(Map<String, Collection<Gene>> regionsToCount, Map<String, Collection<Gene>> regionsToOverlapWith, Map<String, Collection<Gene>> excludeSet) {
		return totalGenes(getOverlap(regionsToCount, regionsToOverlapWith, excludeSet));
	}
	
	/**
	 * Intersect a set of regions with another set of genes and randomize positions of the regions within each gene
	 * @param regions The regions of interest
	 * @param intersectGenes The genes to intersect with and randomize within
	 * @return Regions that have been intersected with the gene set and then randomized within each intersect gene
	 */
	public static Map<String, Collection<Gene>> getOverlapWithRandomizedPositions(Map<String, Collection<Gene>> regions, Map<String, Collection<Gene>> intersectGenes) {
		Map<String, Collection<Gene>> excludeGenes = new TreeMap<String, Collection<Gene>>();
		for(String chr : regions.keySet()) {
			Collection<Gene> empty = new TreeSet<Gene>();
			excludeGenes.put(chr, empty);
		}
		return getOverlapWithRandomizedPositions(regions, intersectGenes, excludeGenes);
	}
	
	/**
	 * Intersect a set of regions with another set of genes and randomize positions of the regions within each gene
	 * @param regions The regions of interest
	 * @param intersectGenes The genes to intersect with and randomize within
	 * @param excludeFeatures Do not consider features that overlap one of these
	 * @return Regions that have been intersected with the gene set and then randomized within each intersect gene
	 */
	public static Map<String, Collection<Gene>> getOverlapWithRandomizedPositions(Map<String, Collection<Gene>> regions, Map<String, Collection<Gene>> intersectGenes, Map<String, Collection<Gene>> excludeFeatures) {
		Map<String, Collection<Gene>> rtrn = new TreeMap<String, Collection<Gene>>();
		logger.info("Getting overlaps with randomized positions...");
		Random rand = new Random();
		for(String chr : intersectGenes.keySet()) {
			Collection<Gene> randomizedOverlapsThisChr = new TreeSet<Gene>();
			if(!regions.containsKey(chr)) continue;
			if(regions.get(chr).isEmpty()) continue;
			logger.info(chr);
			Map<String, Collection<Gene>> genesThisChr = new HashMap<String, Collection<Gene>>();
			genesThisChr.put(chr, regions.get(chr));
			for(Gene feature : intersectGenes.get(chr)) {
				for(Gene gene : genesThisChr.get(chr)) {
					if(feature.getStart() > gene.getEnd()) continue;
					if(feature.getEnd() < gene.getStart()) break;
					if(feature.overlaps(gene) && feature.getStart() >= gene.getStart() && feature.getEnd() <= gene.getEnd()) {
						logger.debug("Region " + gene.getName() + " contains region " + feature.getName());
						int overlapperStart = gene.genomicToTranscriptPosition(feature.getStart());
						int overlapperEndInclusive = gene.genomicToTranscriptPosition(feature.getEnd() - 1);
						int overlapperStartOnTranscript = Math.min(overlapperStart,overlapperEndInclusive);
						int overlapperEndOnTranscript = Math.max(overlapperStart, overlapperEndInclusive) + 1;
						int overlapperLengthOnTranscript = overlapperEndOnTranscript - overlapperStartOnTranscript;
						int lastPossibleRandomStartPos = gene.getSize() - overlapperLengthOnTranscript;
						int randomStartOnTranscript = rand.nextInt(lastPossibleRandomStartPos + 1);
						int randomEndOnTranscript = randomStartOnTranscript + overlapperLengthOnTranscript + 1;
						Gene randomOverlap = gene.copy();
						try {
							randomOverlap = new Gene(randomOverlap.trim(randomStartOnTranscript, gene.getSize() - randomEndOnTranscript));
						} catch (IllegalArgumentException e) {
							logger.warn("Caught exception when randomizing position of feature " + feature.getName() + " within gene " + gene.getName() + ". Skipping.");
							continue;
						}
						randomOverlap.setName(gene.getName() + "_" + feature.getName() + "_random");
						randomOverlap.setOrientation(feature.getOrientation());
						randomizedOverlapsThisChr.add(randomOverlap);
					}
				}
			}
			rtrn.put(chr, randomizedOverlapsThisChr);
		}
		return rtrn;
	}
	
	/**
	 * Get sum of gene sizes
	 * @param genes The gene set by chromosome
	 * @return Total transcribed bases
	 */
	private static int totalBases(Map<String, Collection<Gene>> genes) {
		int rtrn = 0;
		Map<String, Collection<Gene>> collapsed = AnnotationUtils.collapseOverlappers(genes, false);
		for(String chr : collapsed.keySet()) {
			for(Gene gene : collapsed.get(chr)) {
				rtrn += gene.getSize();
			}
		}
		return rtrn;
	}
	
	/**
	 * Get total size of gene set
	 * @param genes The gene set by chromosome
	 * @return Total number of genes
	 */
	private static int totalGenes(Map<String, Collection<Gene>> genes) {
		int rtrn = 0;
		for(String chr : genes.keySet()) rtrn += genes.get(chr).size();
		return rtrn;
	}


	/**
	 * Get number of overlapping bases between two gene sets
	 * @param genes Genes by chromosome
	 * @param intersectSet Set to intersect with by chromosome
	 * @return The total number of overlapping bases
	 */
	public static int numOverlappingBases(Map<String, Collection<Gene>> genes, Map<String, Collection<Gene>> intersectSet) {
		Map<String, Collection<Gene>> excludeGenes = new TreeMap<String, Collection<Gene>>();
		return numOverlappingBases(genes, intersectSet, excludeGenes);
	}
	
	/**
	 * Get number of overlapping bases between two gene sets and automatically exclude all genes in a set
	 * @param genes Genes by chromosome
	 * @param intersectSet Set to intersect with by chromosome
	 * @param excludeSet Exclude genes if they overlap a gene from this set
	 * @return The total number of overlapping bases
	 */
	public static int numOverlappingBases(Map<String, Collection<Gene>> genes, Map<String, Collection<Gene>> intersectSet, Map<String, Collection<Gene>> excludeSet) {
		Map<String, Collection<Gene>> overlap = getOverlap(genes, intersectSet, excludeSet);
		return totalBases(overlap);
	}

	/**
	 * Get number of overlapping genes between two gene sets
	 * @param genes Genes by chromosome
	 * @param intersectSet Set to intersect with by chromosome
	 * @return The total number of overlapping bases
	 */
	public static int numOverlappingGenes(Map<String, Collection<Gene>> genes, Map<String, Collection<Gene>> intersectSet) {
		Map<String, Collection<Gene>> excludeGenes = new TreeMap<String, Collection<Gene>>();
		return numOverlappingGenes(genes, intersectSet, excludeGenes);
	}
	
	/**
	 * Get number of overlapping genes between two gene sets and automatically exclude all genes in a set
	 * @param genes Genes by chromosome
	 * @param intersectSet Set to intersect with by chromosome
	 * @param excludeGenes Exclude genes from first gene set if they overlap a gene from this set
	 * @return The total number of overlapping bases
	 */
	public static int numOverlappingGenes(Map<String, Collection<Gene>> genes, Map<String, Collection<Gene>> intersectSet, Map<String, Collection<Gene>> excludeGenes) {
		Map<String, Collection<Gene>> overlap = getOverlap(genes, intersectSet, excludeGenes);
		return totalGenes(overlap);
	}


	/**
	 * Get overlap between two gene sets specified in bed files and write to a file
	 * @param geneFile Bed file of genes
	 * @param intersectFile Bed file of regions to intersect with
	 * @param outFile Output bed file to write genes representing the overlaps
	 * @param randomizePositionsWithinEachGene Whether to randomize the position of each overlap within the gene
	 * @throws IOException
	 */
	public static void writeOverlap(String geneFile, String intersectFile, String outFile, boolean randomizePositionsWithinEachGene) throws IOException {
		Map<String, Collection<Gene>> genes = BEDFileParser.loadDataByChr(new File(geneFile));
		Map<String, Collection<Gene>> intersectGenes = BEDFileParser.loadDataByChr(new File(intersectFile));
		writeOverlap(genes, intersectGenes, outFile, randomizePositionsWithinEachGene);
	}

	/**
	 * Get overlap between two gene sets specified in bed files, automatically excluding all genes in a set, and write to a file
	 * @param geneFile Bed file of genes
	 * @param intersectFile Bed file of regions to intersect with
	 * @param excludeFile Exclude genes if they overlap a gene from this set
	 * @param outFile Output bed file to write Genes representing the overlaps
	 * @param randomizePositionsWithinEachGene Whether to randomize the position of each overlap within the gene
	 * @throws IOException
	 */
	public static void writeOverlap(String geneFile, String intersectFile, String excludeFile, String outFile, boolean randomizePositionsWithinEachGene) throws IOException {
		logger.info("Writing overlap between " + geneFile + " and " + intersectFile + " to file " + outFile + "...");
		logger.info("Loading genes...");
		Map<String, Collection<Gene>> genes = BEDFileParser.loadDataByChr(new File(geneFile));
		logger.info("Loading annotations to intersect with...");
		Map<String, Collection<Gene>> intersectGenes = BEDFileParser.loadDataByChr(new File(intersectFile));
		logger.info("Loading annotations to exclude...");
		Map<String, Collection<Gene>> excludeGenes = new TreeMap<String, Collection<Gene>>();
		if(excludeFile != null) excludeGenes.putAll(BEDFileParser.loadDataByChr(new File(excludeFile)));
		logger.info("Done loading all annotations.");
		writeOverlap(genes, intersectGenes, excludeGenes, outFile, randomizePositionsWithinEachGene);
	}

	/**
	 * Get overlap between two gene sets and write to a file
	 * @param genes Genes by chromosome
	 * @param intersectSet Set to intersect with by chromosome
	 * @param outFile Output bed file to write Genes representing the overlaps
	 * @param randomizePositionsWithinEachGene Whether to randomize the position of each overlap within the gene
	 * @throws IOException
	 */
	public static void writeOverlap(Map<String, Collection<Gene>> genes, Map<String, Collection<Gene>> intersectSet, String outFile, boolean randomizePositionsWithinEachGene) throws IOException { 
		Map<String,Collection<Gene>> excludeGenes = new TreeMap<String,Collection<Gene>>();
		writeOverlap(genes, intersectSet, excludeGenes, outFile, randomizePositionsWithinEachGene);
	}
	
	/**
	 * Get overlap between two gene sets, automatically excluding all genes in a set, and write to a file
	 * @param genes Genes by chromosome
	 * @param intersectSet Set to intersect with by chromosome
	 * @param excludeSet Exclude genes if they overlap a gene from this set
	 * @param outFile Output bed file to write Genes representing the overlaps
	 * @param randomizePositionsWithinEachGene Whether to randomize the position of each overlap within the gene
	 * @throws IOException
	 */
	public static void writeOverlap(Map<String, Collection<Gene>> genes, Map<String, Collection<Gene>> intersectSet, Map<String, Collection<Gene>> excludeSet, String outFile, boolean randomizePositionsWithinEachGene) throws IOException {
		FileWriter w = new FileWriter(outFile);
		Map<String, Collection<Gene>> overlaps = new TreeMap<String, Collection<Gene>>();
		if(!randomizePositionsWithinEachGene) {
			logger.info("Getting overlaps...");
			overlaps.putAll(getOverlap(genes, intersectSet, excludeSet));
		}
		else overlaps.putAll(getOverlapWithRandomizedPositions(genes, intersectSet, excludeSet));
		logger.info("Writing overlaps to file " + outFile + "...");
		for(String chr : overlaps.keySet()) {
			for(Gene gene : overlaps.get(chr)) {
				w.write(gene.toBED() + "\n");
			}
		}
		w.close();
		logger.info("Done writing file.");
	}
	
	/**
	 * Write basic overlap stats to a file
	 * @param geneFile Gene bed file
	 * @param intersectFile1 Intersect set 1 bed file
	 * @param intersectFile2 Intersect set 2 bed file
	 * @param excludeFile Exclude genes if they overlap a gene in this bed file
	 * @param outFile Output file of stats
	 * @throws IOException
	 */
	public static void writeGeneralStats(String geneFile, String intersectFile1, String intersectFile2, String excludeFile, String outFile) throws IOException {
		logger.info("Writing general overlap stats to file " + outFile + "...");
		FileWriter w = new FileWriter(outFile);
		logger.info("Reading from bed files...");
		Map<String, Collection<Gene>> genes = BEDFileParser.loadDataByChr(geneFile);
		Map<String, Collection<Gene>> set1 = BEDFileParser.loadDataByChr(intersectFile1);
		Map<String, Collection<Gene>> set2 = new TreeMap<String, Collection<Gene>>();
		if(intersectFile2 != null) {
			set2.putAll(BEDFileParser.loadDataByChr(intersectFile2));
		}
		Map<String, Collection<Gene>> exclude = new TreeMap<String, Collection<Gene>>();
		if(excludeFile != null) {
			set2.putAll(BEDFileParser.loadDataByChr(excludeFile));
		}

		if(excludeFile != null) {
			logger.info("Counting number of genes excluded...");
			int numExcluded = numOverlappingGenes(genes, exclude);
			w.write("Excluded " + numExcluded + " genes in " + geneFile + " that overlap a gene in " + excludeFile + ".\n");
			w.write("\n");
		}
		
		logger.info("Computing overlap between genes and set1...");
		Map<String, Collection<Gene>> overlapGenesSet1 = getOverlap(genes, set1, exclude);
		int totalBasesGenes = totalBases(genes);
		int totalBasesSet1 = totalBases(set1);
		int totalBasesOverlapGenesSet1 = totalBases(overlapGenesSet1);
		int totalGenes = totalGenes(genes);
		int totalSet1 = totalGenes(set1);
		int totalSet2 = totalGenes(set2);
		double pctGeneBasesOverlapSet1 = 100 * (double) totalBasesOverlapGenesSet1 / totalBasesGenes;
		double pctGenesOverlapSet1 = 100 * (double) countOverlappers(genes, set1) / totalGenes;
		double pctSet1BasesOverlapGenes = 100 * (double) totalBasesOverlapGenesSet1 / totalBasesSet1;
		double pctSet1GenesOverlapGenes = 100 * (double) countOverlappers(set1, genes) / totalSet1;
		w.write("Overlap between " + geneFile + " and " + intersectFile1 + "\n");
		w.write(totalBasesOverlapGenesSet1 + " bases\n");
		w.write(decimalFormat.format(pctGeneBasesOverlapSet1) + "% of bases in " + geneFile + " (" + totalBasesGenes + " total bases)\n");
		w.write(decimalFormat.format(pctSet1BasesOverlapGenes) + "% of bases in " + intersectFile1 + " (" + totalBasesSet1 + " total bases)\n");
		w.write(decimalFormat.format(pctGenesOverlapSet1) + "% of genes in " + geneFile + " (" + totalGenes + " total regions)\n");
		w.write(decimalFormat.format(pctSet1GenesOverlapGenes) + "% of genes in " + intersectFile1 + " (" + totalSet1 + " total regions)\n");
		w.write("\n");

		if (intersectFile2 != null) {
			
			logger.info("Computing overlap between genes and set2...");
			Map<String, Collection<Gene>> overlapGenesSet2 = getOverlap(genes, set2, exclude);
			int totalBasesSet2 = totalBases(set2);
			int totalBasesOverlapGenesSet2 = totalBases(overlapGenesSet2);
			double pctGeneBasesOverlapSet2 = 100 * (double) totalBasesOverlapGenesSet2 / totalBasesGenes;
			double pctGenesOverlapSet2 = 100 * (double) countOverlappers(genes, set2) / totalGenes;
			double pctSetBases2OverlapGenes = 100 * (double) totalBasesOverlapGenesSet2 / totalBasesSet2;
			double pctSet2GenesOverlapGenes = 100 * (double) countOverlappers(set2, genes) / totalSet2;
			w.write("Overlap between " + geneFile + " and " + intersectFile2 + "\n");
			w.write(totalBasesOverlapGenesSet2 + " bases\n");
			w.write(decimalFormat.format(pctGeneBasesOverlapSet2) + "% of bases in " + geneFile + " (" + totalBasesGenes + " total bases)\n");
			w.write(decimalFormat.format(pctSetBases2OverlapGenes) + "% of bases in " + intersectFile2 + " (" + totalBasesSet2 + " total bases)\n");
			w.write(decimalFormat.format(pctGenesOverlapSet2) + "% of genes in " + geneFile + " (" + totalGenes + " total regions)\n");
			w.write(decimalFormat.format(pctSet2GenesOverlapGenes) + "% of genes in " + intersectFile2 + " (" + totalSet2 + " total regions)\n");
			w.write("\n");			
			
			logger.info("Computing overlap between genes, set1 and set2...");
			Map<String, Collection<Gene>> overlapGenesSet1Set2 = getOverlap(overlapGenesSet1, set2, exclude);
			int totalBasesOverlapGenesSet1Set2 = totalBases(overlapGenesSet1Set2);
			double pctGeneBasesOverlapSet1Set2 = 100 * (double) totalBasesOverlapGenesSet1Set2 / totalBasesGenes;
			double pctGenesOverlapSet1Set2 = 100 * (double) countOverlappers(genes, overlapGenesSet1Set2) / totalGenes;
			double pctSet1BasesOverlapGenesSet2 = 100 * (double) totalBasesOverlapGenesSet1Set2 / totalBasesSet1;
			double pctSet1GenesOverlapGenesSet2 = 100 * (double) countOverlappers(set1, overlapGenesSet1Set2) / totalSet1;
			double pctSet2BasesOverlapGenesSet1 = 100 * (double) totalBasesOverlapGenesSet1Set2 / totalBasesSet2;
			double pctSet2GenesOverlapGenesSet1 = 100 * (double) countOverlappers(set2, overlapGenesSet1Set2) / totalSet2;
			w.write("Overlap between " + geneFile + " and " + intersectFile1 + " and " + intersectFile2 + "\n");
			w.write(totalBasesOverlapGenesSet1Set2 + " bases\n");
			w.write(decimalFormat.format(pctGeneBasesOverlapSet1Set2) + "% of bases in " + geneFile + " (" + totalBasesGenes + " total bases)\n");
			w.write(decimalFormat.format(pctSet1BasesOverlapGenesSet2) + "% of bases in " + intersectFile1 + " (" + totalBasesSet1 + " total bases)\n");
			w.write(decimalFormat.format(pctSet2BasesOverlapGenesSet1) + "% of bases in " + intersectFile2 + " (" + totalBasesSet2 + " total bases)\n");
			w.write(decimalFormat.format(pctGenesOverlapSet1Set2) + "% of genes in " + geneFile + " (" + totalGenes + " total regions)\n");
			w.write(decimalFormat.format(pctSet1GenesOverlapGenesSet2) + "% of genes in " + intersectFile1 + " (" + totalSet1 + " total regions)\n");
			w.write(decimalFormat.format(pctSet2GenesOverlapGenesSet1) + "% of genes in " + intersectFile2 + " (" + totalSet2 + " total regions)\n");
			w.write("\n");
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
		p.addStringArg("-b", "Bed file of genes", true);
		p.addStringArg("-i1", "Bed file containing regions to intersect with first (set1)", true);
		p.addStringArg("-i2", "Bed file containing regions to intersect with second (set2)", false, null);
		p.addBooleanArg("-r", "Only has effect if -ob is provided. Intersect genes with set1 and write bed file. Ignore set2. If true, randomize position of set1 overlap within gene. ", false, false);
		p.addIntArg("-nr", "Randomize overlap positions of set1 within genes this many times, intersect with set2, and print counts only.", false, 0);
		p.addStringArg("-e", "Bed file containing regions to exclude - exclude any gene that overlaps one of these", false, null);
		p.addStringArg("-ob", "Output bed file for intersection with set1 (randomized or not)", false, null);
		p.addStringArg("-oc", "Output counts file for intersections with set1, randomized and intersected with set2. Requires -i2 and -nr.", false, null);
		p.addBooleanArg("-d", "Debug logging", false, false);
		p.addStringArg("-og", "Output file for general overlap stats", false, null);
		p.parse(args);
		boolean debug = p.getBooleanArg("-d");
		if(debug) {
			logger.setLevel(Level.DEBUG);
		}
		String geneFile = p.getStringArg("-b");
		String intersectFile1 = p.getStringArg("-i1");
		String intersectFile2 = p.getStringArg("-i2");
		String excludeFile = p.getStringArg("-e");
		int numRandom = p.getIntArg("-nr");
		boolean randomize = p.getBooleanArg("-r");
		String outBed = p.getStringArg("-ob");
		String outRandCounts = p.getStringArg("-oc");
		String outGeneralStats = p.getStringArg("-og");
		
		if(outGeneralStats != null) {
			writeGeneralStats(geneFile, intersectFile1, intersectFile2, excludeFile, outGeneralStats);
		}
		
		if(outBed != null) {
			writeOverlap(geneFile, intersectFile1, excludeFile, outBed, randomize);
		}
		
		if(outRandCounts != null) {
			
			if(intersectFile2 == null) {
				throw new IllegalArgumentException("-i2 is required for randomized overlap counts.");
			}
			
			FileWriter w = new FileWriter(outRandCounts);
			w.write("test\ttotal_bases_overlap\n");
			Map<String, Collection<Gene>> genes = BEDFileParser.loadDataByChr(new File(geneFile));
			Map<String, Collection<Gene>> intersectFeatures1 = BEDFileParser.loadDataByChr(new File(intersectFile1));
			Map<String, Collection<Gene>> intersectFeatures2 = BEDFileParser.loadDataByChr(new File(intersectFile2));
			Map<String, Collection<Gene>> excludeFeatures = new TreeMap<String, Collection<Gene>>();
			if(excludeFile != null) excludeFeatures.putAll(BEDFileParser.loadDataByChr(new File(excludeFile)));
			
			logger.info("");
			logger.info("Getting overlap between regions and set1...");
			Map<String, Collection<Gene>> overlapSet1 = getOverlap(genes, intersectFeatures1, excludeFeatures, true);
			logger.info("Overlap is " + totalGenes(overlapSet1) + " genes and " + totalBases(overlapSet1) + " bases.");
			
			logger.info("");
			logger.info("Getting overlap between regions and set1+set2...");
			Map<String, Collection<Gene>> overlapSet2 = getOverlap(overlapSet1, intersectFeatures2, true);
			logger.info("Overlap is " + totalGenes(overlapSet2) + " genes and " + totalBases(overlapSet2) + " bases.");
			logger.info("COUNT\treal_overlap\t" + totalBases(overlapSet2));
			w.write("real_overlap\t" + totalBases(overlapSet2) + "\n");
			
			logger.info("");
			logger.info("Randomizing overlap positions of set1 within regions...");
			for(int i=0; i < numRandom; i++) {
				logger.info("Random iteration " + Integer.valueOf(i+1).toString());
				Map<String, Collection<Gene>> randomizedOverlap = getOverlapWithRandomizedPositions(genes, intersectFeatures1);
				Map<String, Collection<Gene>> randomizedWithSet2 = getOverlap(randomizedOverlap, intersectFeatures2, true);
				//writeOverlap(randomizedOverlap, intersectFeatures2, "test_" + i + "_" + outRandCounts + ".bed", false);
				logger.info("COUNT\toverlap_randomized_" + i + "\t" + totalBases(randomizedWithSet2));
				w.write("overlap_randomized_" + i + "\t" + totalBases(randomizedWithSet2) + "\n");
			}
			
			logger.info("");
			logger.info("All done. Wrote counts to file " + outRandCounts + ".");
			
			w.close();
		}

	}

}
