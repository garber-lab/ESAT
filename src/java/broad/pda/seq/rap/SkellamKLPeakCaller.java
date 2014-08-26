package broad.pda.seq.rap;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import broad.pda.annotation.BEDFileParser;

import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.util.Log;
import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Gene;
import nextgen.core.model.AlignmentModel;
import nextgen.core.model.score.WindowProcessor;
import nextgen.core.model.score.WindowScoreIterator;
import nextgen.core.coordinatesystem.GenomicSpace;

public class SkellamKLPeakCaller extends GenomeCommandLineProgram {

	private static final Log logger = Log.getInstance(SkellamKLPeakCaller.class);
	
    @Usage
    public String USAGE = "Computes significant peaks in RAP data from the input sample using the skellam statistic and KL-divergence." +
    		" Run RAPPeakCaller -task analyze";
    
	@Option(doc="Window size")
	public int WINDOW;
	
	@Option(doc="Overlap between windows")
	public int OVERLAP;

	@Option(doc="Alignment (mapped to genome) of RAP sample SAM or BAM file")
	public File TARGET;

	@Option(doc="Alignment (mapped to genome) of input/AS sample SAM or BAM file for normalization.")
	public File BACKGROUND;
	
	@Option(doc="Output file basename")
	public File OUTPUT;
	
	@Option(doc="Percentile cut-off (e.g. 0.001)", optional=true)
	public Double CUTOFF = 0.99;
	
	@Option(doc="Whether to look for regions BELOW the given cutoff, rather than above.")
	public Boolean DEPLETED=false;
	
	/**
	 * 
	 * @param args
	 */
	public static void main(final String[] args) {
		System.exit(new SkellamKLPeakCaller().instanceMain(args));
	}
	
	/**
	 * 
	 */
	protected int doWork(){
		
		try {
		Map<String, Collection<Gene>> maskedRegions = BEDFileParser.loadDataByChr(new File(MASK_FILE));
		
		Map<String,Integer> chrSizes = BEDFileParser.loadChrSizes(SIZES);
		//Get all the regions to be scores
		Map<String,Annotation> regions = getRegionsMap();
		
		// Get the coordinate space. Fails if not genomic space. To do: allow transcriptome space.
		((GenomicSpace) getCoordinateSpace()).setPercentMaskedAllowed(PCT_MASKED_ALLOWED);
		AlignmentModel rap = loadAlignmentModel(TARGET);
		AlignmentModel background = loadAlignmentModel(BACKGROUND);

		BufferedWriter bw = new BufferedWriter(new FileWriter(OUTPUT,true));
		BufferedWriter bwAll = new BufferedWriter(new FileWriter(new File(OUTPUT.getAbsolutePath() + ".allWindowDistributions.bed")));
		//For each chromosome
		for (String chr : regions.keySet()) {
			Annotation region = regions.get(chr);
			
			//Calculate lambdas
			double backgroundLambda = background.getRefSequenceLambda(chr);
			double length = (double)chrSizes.get(chr) - (double)getMaskedLength(maskedRegions,chr);
			backgroundLambda = (backgroundLambda*background.getRefSequenceLength(chr))/length;
			double rapLambda = rap.getRefSequenceLambda(chr);
			rapLambda = (rapLambda*rap.getRefSequenceLength(chr))/length;
			
			if(rapLambda ==0.0){
				logger.error(chr +" is not expressed in RAP file");
			}
			else{
				logger.info("Starting: " + region.toUCSC());
				SkellamScore.Processor processor = new SkellamScore.Processor(rap, background,rapLambda,backgroundLambda);
				List<Annotation> sigWindows;
				
				sigWindows = getSignificantWindows(rap, background, CUTOFF, region, processor,backgroundLambda,rapLambda,bwAll);
				
				List<SkellamScore> scoredWindows = scoreWindows(processor, sigWindows);
				writeResults(scoredWindows, bw);
				
			}
		}
		
			bwAll.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return 0;
	}
	
	private void writeResults(List<SkellamScore> scored, BufferedWriter bw) throws IOException {
		for (SkellamScore score : scored) {
			bw.write(score.toString());
		}
	}
	private List<SkellamScore> scoreWindows(WindowProcessor<SkellamScore> processor, List<? extends Annotation> windows) {
		// Using the same processor initialized above preserves the region counts
		
		LinkedList<SkellamScore> scored = new LinkedList<SkellamScore>();
		for (Annotation window : windows) {
			scored.addLast(processor.processWindow(window));
		}
		return scored;
	}
	/**
	 * Returns the length on this chromosome that is masked
	 * @param chr
	 * @return
	 */
	private int getMaskedLength(Map<String, Collection<Gene>> maskedRegions,String chr){
		
		int len = 0;
		if(maskedRegions.containsKey(chr)){
			for(Gene g:maskedRegions.get(chr)){
				len = g.getSize();
			}
		}
		return len;
	}

	private List<Annotation> getSignificantWindows(AlignmentModel rap,
			AlignmentModel background2, Double cutoff, Annotation region, WindowProcessor<SkellamScore> processor,
			double backgroundLambda, double rapLambda,BufferedWriter bw) throws IOException {
			
		LinkedList<Annotation> sigWindows = new LinkedList<Annotation>();
		int count = 0;
		WindowScoreIterator<SkellamScore> itr = rap.scan(region, WINDOW, OVERLAP, processor);
		while (itr.hasNext()) {
			SkellamScore curr = itr.next();
			
			if(curr.countLessThanLambda()){
				continue; 
			}
			else{
				double skellam = curr.getSkellamPValue();
				bw.write(curr.toString());
				if(skellam<CUTOFF){
					count++;
					
					// If this window overlaps the previous one, combine into one window
					boolean overlapping = false;
					if (sigWindows.size() > 0) {
						Annotation previous = sigWindows.getLast();
						if (curr.getAnnotation().overlaps(previous)) {
							sigWindows.removeLast();
							sigWindows.addLast(previous.union(curr.getAnnotation()));
							overlapping = true;
						}
					} 

					if (!overlapping) sigWindows.addLast(curr.getAnnotation());
				}
			}				
		}
		itr.close();
				
		logger.info("Found " + count + " significant windows.");
		logger.info("Collapsed to " + sigWindows.size() + " non-overlapping regions.");
		return sigWindows;
	}
}
