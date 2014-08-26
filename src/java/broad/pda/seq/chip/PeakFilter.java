package broad.pda.seq.chip;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Iterator;
import java.util.List;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;

import broad.core.annotation.BED;
import broad.core.annotation.BEDReader;
import broad.core.error.ParseException;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.datastructures.Alignments;
import broad.pda.seq.alignment.AlignmentUtils;
import broad.pda.seq.segmentation.ContinuousDataAlignmentModel;

public class PeakFilter {
	private static final Logger logger =  Logger.getLogger(PeakFilter.class.getName());
	static final String USAGE = "\nTasks" +
	"\n filterByOverlap -in <Set of peaks in BED format> -null <Null set of peaks> -maxOverlap <Peaks that are covered by more than the specified overlap are filtered> -out <Name of the filtered output file> "+
	"\n filterByEnrichment -in <Set of peaks in BED format>  -null <Null alignment>  -libAlignment <Alignment for the library> -minFoldOverNull <Minimum fold score over null peak when overlapping a null peal> -out <Name of the filtered output file>" +
	"\n";
	
	
	public static void main(String [] args) throws IOException, ParseException {
		ArgumentMap argMap = CLUtil.getParameters(args,USAGE , "filterByOverlap");
		Globals.setHeadless(true);
		System.out.println("PeakFilter V0.1");
		if ("filterByOverlap".equalsIgnoreCase(argMap.getTask())) {
			String inputFile = argMap.getInput();
			double maxOverlap = argMap.getDouble("maxOverlap");
			String nullPeaksFile =  argMap.getMandatory("null");

			BEDReader peaks = new BEDReader(inputFile);
			logger.info("Loaded peaks");
			BEDReader nullPeaks = new BEDReader(nullPeaksFile);
			logger.info("Loaded null peaks");
			nullPeaks.merge();
			List<BED> peakList = peaks.getAnnotationList();
			nullPeaks.intersect(peaks);
			logger.info("Prepared null peaks");
			BufferedWriter bw = argMap.getOutputWriter();
			int total = 0;
			int passing = 0;
			for(BED peak : peakList) {
				List<BED> nullOverlap = nullPeaks.getOverlappers(peak);
				int coverage = 0;
				total++;
				for(BED nullPeak : nullOverlap) {
					coverage += nullPeak.length();
				}
				//System.err.println("Peak " + peak.toUCSC() + " has overlap " + (peak.length()/(double)coverage) + " maxoverlap is " + maxOverlap);
				if(coverage == 0 || peak.length()/(double)coverage < maxOverlap) {
					bw.write(peak.toShortString());
					bw.newLine();
					passing++;
				}
			}
			bw.close();
			logger.info("Done original peaks: " + total + ", passing peaks " + passing);
			
		}else if ("filterByNullEnrichment".equalsIgnoreCase(argMap.getTask())) {
			double minFold    = argMap.getDouble("minFoldOverNull");
			String inputFile = argMap.getInput();
			String nullAlignmentFile =  argMap.getMandatory("null");
			String libraryAlignmentFile  = argMap.getMandatory("libAlignment");
			int window = argMap.containsKey("window") ? argMap.getInteger("window") : 0;
			BEDReader peaks = new BEDReader(inputFile);
			logger.info("Loaded peaks");
			boolean loadPairsAsFragments = argMap.containsKey("loadPairsAsFragments") || argMap.containsKey("pairedEnd");
			
			ContinuousDataAlignmentModel nullData = AlignmentUtils.loadAlignmentData(nullAlignmentFile, true, 0, true, false, null, loadPairsAsFragments);
			ContinuousDataAlignmentModel libData = AlignmentUtils.loadAlignmentData(libraryAlignmentFile, true, 0, true, false, null, loadPairsAsFragments);			

			logger.info("Done loading " + libraryAlignmentFile);
			
			if (window > 0) {
				scorePeaks(peaks, libData, nullData, window, minFold);
			} else {
				scorePeaksWholeRegion(peaks, libData, nullData, minFold);
			}
			
			List<BED> peakList = peaks.getAnnotationList();
			BufferedWriter bw = argMap.getOutputWriter();
			int total = 0;
			int passing = 0;
			for(BED peak : peakList) {
				double fold = peak.getExtraScore(peak.getExtraScores().size() - 1);
				//System.err.println("peak " + peak.toUCSC() + " fold " + peak.getExtraScore(peak.getExtraScores().size() - 1));
				if(fold > minFold) {
					//peak.setScore(fold);
					//bw.write(peak.toShortString());
					bw.write(peak.toBED());
					logger.debug("Doing full bed");
					bw.newLine();
					passing++;
				}
				total++;
			}
			bw.close();
			logger.info("Done original peaks: " + total + ", passing peaks " + passing);
			
		}
	}
	

	private static void scorePeaks(BEDReader annotationReader,	ContinuousDataAlignmentModel libData, ContinuousDataAlignmentModel nullData, int window, double minFold) throws IOException {
		Iterator<String> chrIt = annotationReader.getChromosomeIterator();
		while(chrIt.hasNext()) {
			String chr = chrIt.next();
			List<BED> chrAnnotations = annotationReader.getChromosomeBEDs(chr);
			for(BED chrAnnotation: chrAnnotations) {
				if(chrAnnotation.length() < window) {
					int toAdd = window - chrAnnotation.length() +2;
					chrAnnotation.setStart(chrAnnotation.getStart() - toAdd/2);
					chrAnnotation.setEnd(chrAnnotation.getEnd() + toAdd/2);
				}
				//System.err.print("Processing " + chrAnnotation.toUCSC());
				int s = chrAnnotation.getStart();
				double maxScore = 0;
				int maxScoreStart = s;
				double maxScorePval = 1;
				double maxScoreNullPval = 1;
				double fold = 0;
				do {
					Alignments tmp = new Alignments(chr, s, Math.min(s + window, chrAnnotation.getEnd()));
					double [] libScores = libData.scoreSegment(tmp);
					double [] nullScores = nullData.scoreSegment(tmp);
					//System.err.println("lib score " + libScores[1] + " null score " + nullScores[1]);
					double foldOverNull = libScores[1]/(nullScores[1]+0.01);  //0 pvalue, 1 enrichment, 2 read number, 3 reads/base
					if(libScores[1] > maxScore && foldOverNull > minFold) {
						maxScoreStart = s;
						maxScore = libScores[1];
						fold = foldOverNull;
						maxScorePval = libScores[0];
						maxScoreNullPval = nullScores[0];
					}
					//System.err.println("fold over null " + fold);
					s += 10;
					
				} while(s + window < chrAnnotation.getEnd());
				chrAnnotation.addExtraScore(maxScorePval);
				chrAnnotation.addExtraScore(maxScoreNullPval);
				chrAnnotation.addExtraScore(chrAnnotation.getScore());
				chrAnnotation.addExtraScore(maxScore);
				chrAnnotation.addExtraScore(fold);
				
				chrAnnotation.setScore(maxScore);
				if(maxScore > 0 ) {
					chrAnnotation.setThickEnd(Math.min(maxScoreStart + window, chrAnnotation.getEnd()));
					chrAnnotation.setThickStart(maxScoreStart);
				}
			}
			
			nullData.resetTreeCache();
			libData.resetTreeCache();
		}
	}
	
	
	/**
	 * Jesse Engreitz
	 * June 29, 2012
	 * This function allows for scoring enrichment based on read counts for the entire peak, rather than picking the 
	 * window with the best fold enrichment.
	 * @param annotationReader
	 * @param libData
	 * @param nullData
	 * @param minFold
	 * @throws IOException
	 */
	private static void scorePeaksWholeRegion(BEDReader annotationReader, ContinuousDataAlignmentModel libData, ContinuousDataAlignmentModel nullData, double minFold) throws IOException {
		Iterator<String> chrIt = annotationReader.getChromosomeIterator();
		while(chrIt.hasNext()) {
			String chr = chrIt.next();
			List<BED> chrAnnotations = annotationReader.getChromosomeBEDs(chr);
			for(BED chrAnnotation: chrAnnotations) {
				Alignments tmp = new Alignments(chr, chrAnnotation.getStart(), chrAnnotation.getEnd());			
				double [] libScores = libData.scoreSegment(tmp);
				double [] nullScores = nullData.scoreSegment(tmp);
					//System.err.println("lib score " + libScores[1] + " null score " + nullScores[1]);
				double foldOverNull = libScores[1]/(nullScores[1]+0.01);  //0 pvalue, 1 enrichment, 2 read number, 3 reads/base
				chrAnnotation.addExtraScore(libScores[0]);
				chrAnnotation.addExtraScore(nullScores[0]);
				chrAnnotation.addExtraScore(chrAnnotation.getScore());
				chrAnnotation.addExtraScore(libScores[1]);
				chrAnnotation.addExtraScore(foldOverNull);
				chrAnnotation.setScore(libScores[1]);
			}
			nullData.resetTreeCache();
			libData.resetTreeCache();
		}
	}
	
}
