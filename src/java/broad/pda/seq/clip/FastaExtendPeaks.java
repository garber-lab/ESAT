package broad.pda.seq.clip;

import broad.core.parser.CommandLineParser;

import java.io.BufferedWriter;
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

import nextgen.core.annotation.Annotation.Strand;
import nextgen.core.annotation.Gene;

/**
 * @author shari
 * Extends peaks to standardized size (for probe design, etc.)
 */
public class FastaExtendPeaks {

	private Map<String, Collection<Gene>> genes;
	private Map<String, Gene> genesByName;
	private Map<String, Collection<Gene>> peaks;
	private Map<String, Collection<Gene>> extendedPeaks;
	private Map<String, Sequence> chromosomes;
	private int windowSize;
	private static Logger logger = Logger.getLogger(FastaExtendPeaks.class.getName());
	private static int DEFAULT_WINDOW_SIZE = 120;
	
	/**
	 * @param genomeFasta Fasta file of sequences
	 * @param bedFile Bed file
	 * @throws IOException
	 */
	public FastaExtendPeaks(String genomeFasta, String bedFile, String peakFile, int windowSize) throws IOException {
		this(FastaSequenceIO.getChrSequencesFromFasta(genomeFasta), bedFile, peakFile, windowSize);
	}
	
	/**
	 * Load sequences and annotations
	 * @param chrsByName Chromosomes by name
	 * @param bedFile Genes in bed format
	 * @param windowSize size of return sequences
	 * @throws IOException
	 */
	public FastaExtendPeaks(Map<String, Sequence> chrsByName, String bedFile, String peakFile, int windowSize) throws IOException {
		
		logger.info("Loading genes from file " + bedFile + "...");
		
		genes = BEDFileParser.loadDataByChr(new File(bedFile));
		peaks = BEDFileParser.loadDataByChr(new File(peakFile));
		int numGenes = 0;
		for(String chr : genes.keySet()) {
			numGenes += genes.get(chr).size();
		}
		logger.info("Loaded " + numGenes + " genes.");
		int numPeaks = 0;
		for(String chr : peaks.keySet()) {
			numPeaks += peaks.get(chr).size();
		}
		logger.info("Loaded " + numPeaks + " peaks.");
		
		genesByName = new TreeMap<String, Gene>();
		for (String chr : genes.keySet()){
			Collection<Gene> chrGenes = genes.get(chr);
			if (chrGenes.size() > 0){
				for (Gene gene : chrGenes){
					genesByName.put(gene.getName(), gene);
				}
			}
		}
		
		chromosomes = chrsByName;
		setWindowSize(windowSize);
		//setWindowSize(DEFAULT_WINDOW_SIZE);
		
	}
	
	/**
	 * Load sequences and annotations
	 * @param chrsByName Chromosomes by name
	 * @param bedFile Genes in bed format
	 * @throws IOException
	 */
	public FastaExtendPeaks(Map<String, Sequence> chrsByName, String bedFile) throws IOException {
		
		logger.info("Loading genes from file " + bedFile + "...");
		
		genes = BEDFileParser.loadDataByChr(new File(bedFile));
		int numGenes = 0;
		for(String chr : genes.keySet()) {
			numGenes += genes.get(chr).size();
		}
		logger.info("Loaded " + numGenes + " genes.");
		
		chromosomes = chrsByName;
		
	}
	
	public FastaExtendPeaks(Map<String, Sequence> chrsByName, String bedFile, String peakFile) throws IOException {
		this(chrsByName,bedFile,peakFile,DEFAULT_WINDOW_SIZE);
	}
	
	public FastaExtendPeaks(String genomeFasta, String bedFile) throws IOException {
		this(FastaSequenceIO.getChrSequencesFromFasta(genomeFasta), bedFile);
	}

	/**
	 * Extend peaks to length windowSize
	 */
	private Map<String,Collection<Gene>> getExtendedPeaks(){
		
		Map<String, Collection<Gene>> rtrn=new TreeMap<String, Collection<Gene>>();
		
		for(String chr : peaks.keySet()){
			
			for (Gene peak : peaks.get(chr)){
				String[] parts = peak.getName().split(":");
				String geneName = parts[0];
				if (!genesByName.containsKey(geneName)){
					throw new IllegalArgumentException("Peak " + geneName + " not in gene file.");
				}
				Gene gene = genesByName.get(geneName);
				Gene newPeak = null;
				int midPoint = gene.genomicToTranscriptPosition(peak.getMidpointGenomicCoords());
				if (gene.getOrientation().equals(Strand.NEGATIVE)) {
					midPoint = gene.length() - midPoint;
				}
				if (midPoint > windowSize/2 & midPoint < gene.length()-windowSize/2){
					newPeak = gene.trimGene(midPoint-windowSize/2, midPoint+windowSize/2);
					logger.info("1");
				} else if (midPoint > windowSize/2 & gene.length() > windowSize){
					newPeak = gene.trimGene(gene.length()-windowSize, gene.length());
					logger.info("2");
				} else if (midPoint < gene.length()-windowSize/2 & gene.length() > windowSize){
					newPeak = gene.trimGene(0, windowSize);
					logger.info("3");
				} 
				logger.info("MIDPOINT\t" + midPoint);
				logger.info("STRAND\t" + gene.getOrientation());
				logger.info("ORIGINAL_PEAK\tcoords=" + peak.toBED());
				logger.info("GENE\tcoords=" + gene.toBED());
				logger.info("NEW_PEAK\tcoords=" + newPeak.toBED());
				logger.info("OVERLAP\t" + gene.getOverlap(newPeak));
				
				if (newPeak != null){
					newPeak.setName(peak.getName());
					newPeak.setOrientation(peak.getOrientation());
					Collection<Gene> data=new TreeSet<Gene>();
					if(rtrn.containsKey(newPeak.getChr())){
						data.addAll(rtrn.get(newPeak.getChr()));
					}
					data.add(newPeak);
					rtrn.put(newPeak.getChr(), data);
				}
			}
		}
		return rtrn;
	}
	
	
	/**
	 * Get the spliced sequence of the gene
	 * @param gene The gene
	 * @return The sequence in 5' to 3' orientation
	 */
	private Sequence getGeneSequence(Gene gene) {
		String chr = gene.getReferenceName();
		return chromosomes.get(chr).getSubsequence(gene);
		
	}
	
	/**
	 * 
	 * @param outFile
	 * @throws IOException
	 */
	public void writeFasta(Map<String,Collection<Gene>> genesToWrite, String outFile) throws IOException {
		logger.info("Writing gene sequences to file " + outFile + "...");
		FileWriter w = new FileWriter(outFile);
		BufferedWriter b = new BufferedWriter(w);
		FastaSequenceIO fsio = new FastaSequenceIO();
		for(String chr : genesToWrite.keySet()) {
			for(Gene gene : genesToWrite.get(chr)) {
				try {
					Sequence seq = getGeneSequence(gene);
					seq.setId(gene.getName() + ":" + gene.getScore());
					seq.uppercase();
					fsio.write(seq, b);
				} catch (NullPointerException e) {
					logger.warn("Null pointer: " + gene.toBED());
				}
			}
		}
		b.close();
		logger.info("Done writing sequences.");
	}
	
	private void setWindowSize(int size){
		this.windowSize = size;
	}
	
	private void setExtendedPeaks(){
		this.extendedPeaks = getExtendedPeaks();
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-b", "Bed file of genes",false,null);
		p.addStringArg("-p", "Bed file of peaks", true);
		p.addStringArg("-f", "Fasta file of chromosomes", true);
		p.addIntArg("-w", "Fixed peak size", false, 120);
		p.addStringArg("-o", "Output fasta file of sequences", true);
		p.parse(args);
		String bedFile = p.getStringArg("-b");
		String peakFile = p.getStringArg("-p");
		String fastaFile = p.getStringArg("-f");
		String outFile = p.getStringArg("-o");
		int windowSize = p.getIntArg("-w");
		
		if (bedFile != null) {
			FastaExtendPeaks f = new FastaExtendPeaks(fastaFile, bedFile, peakFile, windowSize);
			f.setExtendedPeaks();
			f.writeFasta(f.extendedPeaks,outFile);
		} else {
			FastaExtendPeaks f = new FastaExtendPeaks(fastaFile,peakFile);
			f.writeFasta(f.genes, outFile);
		}
	}
	
}
