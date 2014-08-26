package broad.core.motif;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import nextgen.core.annotation.Annotation;

import org.apache.log4j.Logger;

import broad.core.annotation.AnnotationReader;
import broad.core.annotation.BED;
import broad.core.annotation.GenomicAnnotation;
import broad.core.annotation.LightweightGenomicAnnotation;
import broad.core.math.ComputeFDR;
import broad.core.math.EmpiricalDistribution;
import broad.core.math.Statistics;
import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.core.sequence.SequenceRegion;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;

public class PWMUtils {
	static Logger logger = Logger.getLogger(PWMUtils.class.getName());
	private static final int SAMPLE_NUM = 2000000;
	private static final NumberFormat DEFAULT_FORMATTER = NumberFormat.getNumberInstance();
	private static final int DEFAULT_SHUFFLES = 100;
	public static final String USAGE = "Usage: ArrayDesignUtilities TASK=<task_num> <task_args>\n" +
	"\tTasks:\n" +
	"\t\t1. Compute log likelohood score thresholds for a given PWM set. This task will compute the distribution of scores then output thresholds for commonly used quantiles. -pwm <File containing one or more PWM matrices> -pA <new equilibrium probability for A>  -pC <new equilibrium probability for C>  -pG <new equilibrium probability for G>  -pT <new equilibrium probability for T> -out <file name for score quantile cutoff per pwm file or standard out if none is specified> -outdir <output directory or current directory if this parameter is not specified> [-prefix <A prefix for each pwm distribution file>]> " +
	"\n\t\tscan. Scan a set a set of genomic regions for PWM ocurrencens -pwm <File containing one or more PWM matrices> -regions <Genomic Annotations> -seqdir <Sequence directory> -maxOnly <If only the maximal hit per region should be reported> -outdir <Output directory to write output file> [-minScoreToWrite <Filter output by minimum score, reudces output file size> -prefix <Prefix to use for each PWM hit file> -shuffles <Number of pwm shuffles to asses significance>]"+
	"\n\t\tpermute. Permute a PWM matrix -pwm <File contining one or more PWM matrics> -out <Output file with permuted matrices or standard out in not specified>"+
	"\n\t\treverse. Reverse-compliment a PWM -pwm <File contining one or more PWM matrics> -out <Output file with permuted matrices or standard out in not specified>" +
	"\n\t\tscanGenome. Scan the genome for PWM ocurrencens -pwm <File containing one or more PWM matrices> -seqdir <Sequence directory> -outdir <Output directory to write output file> -minScoreToWrite <Filter output by minimum score, reudces output file size> [-prefix <Prefix to use for each PWM hit file> -shuffles <Number of pwm shuffles to asses significance>]"+
	"\n";
	
	public static void main (String [] args) throws IllegalArgumentException, Exception {
		ArgumentMap argMap = CLUtil.getParameters(args, USAGE,"scan");
		
		if("scan".equalsIgnoreCase(argMap.getTask())) {
			String seqDir = argMap.getMandatory("seqdir");			
			Map<String, Sequence> go = FastaSequenceIO.loadSequencesByNameFromDirectory(new File(seqDir));
			int shuffles = argMap.containsKey("shuffles") ? argMap.getInteger("shuffles") : DEFAULT_SHUFFLES;
			String prefix = argMap.containsKey("prefix") ? argMap.getMandatory("prefix") : "";
			double minScoreToWrite = argMap.containsKey("minScoreToWrite") ? argMap.getDouble("minScoreToWrite") : Double.NEGATIVE_INFINITY;
			boolean maxOnly = argMap.containsKey("maxOnly");
			String outdir = argMap.getOutputDir();
			
			Map<String, List<? extends GenomicAnnotation>> regionChrMap = argMap.getRegionMapFromParameters();
			
			String pwmFile = argMap.getMandatory("pwm");
			PositionWeightMatrixIO pwmIO = new PositionWeightMatrixIO();
			FileInputStream fis = new FileInputStream(pwmFile);
			pwmIO.load(fis);
			fis.close();
			pwmIO.addPseudoCounts();	
			List<PositionWeightMatrix> pwms = pwmIO.getMatrices();
			
			double []  bgNucleotideFreqs = computeSequenceBGSequences(go, regionChrMap);
			logger.info("Background nuclotide frequencies: (" + bgNucleotideFreqs[0] +", " +bgNucleotideFreqs[1] +", " + bgNucleotideFreqs[1] +", " + bgNucleotideFreqs[3] +")");
			
			Iterator<String> chrIt = regionChrMap.keySet().iterator();
			while(chrIt.hasNext()) {
				String chr = chrIt.next();
				Iterator<? extends GenomicAnnotation> annotIt = regionChrMap.get(chr).iterator();
				Sequence chrSeq = go.get(chr);
				short [] encChrSeq = chrSeq.encodeSequenceIgnoreCase();
				System.err.print(" -sequence encoded(" + (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1000 +")- ");
				System.err.println(" -sequence unloaded(" + (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1000 + ")- starting region scanning." );
				while(annotIt.hasNext()) {
					LightweightGenomicAnnotation annot = annotIt.next();
					scanRegion(shuffles, prefix, minScoreToWrite, maxOnly,
							outdir, pwms, chr, annot, bgNucleotideFreqs,
							encChrSeq);
					//System.err.println("\t\tFinished region " + annot.toUCSC() + ", memory: " +  (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1000 );
				}
				System.err.println("Finished chromosome " + chr + ", memory: " +  (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1000 );
			}		
		} else if ("1".equals(argMap.getTask())) {	
			String prefix = argMap.containsKey("prefix") ? argMap.get("prefix") : "";
			String outdir = argMap.getOutputDir();
			float pA = argMap.getFloat("pA");
			float pC = argMap.getFloat("pC");
			float pG = argMap.getFloat("pG");
			float pT = argMap.getFloat("pT");
			String pwmFile = argMap.getMandatory("pwm");
			int binNum = 90;
			
			if(pA + pC + pG + pT != 1) {
				System.err.println("ERROR: pA + pC + pG + pT = " + (pA + pC + pG + pT) + " it is not a distribution");
				return;
			}

			PositionWeightMatrixIO pwmIO = new PositionWeightMatrixIO();
			FileInputStream fis = new FileInputStream(pwmFile);
			pwmIO.load(fis);
			fis.close();
			
			List<PositionWeightMatrix> pwmList = pwmIO.getMatrices();
			BufferedWriter oBW = argMap.getOutputWriter();
			double [] quantiles = {0.5, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99, 0.999, 0.9999};
			oBW.write("PWM\t50%\t75%\t80%\t85%\t90%\t95%\t97.5%\t99%\t99.9%\t99.99\n");
			for( PositionWeightMatrix pwm : pwmList) {
				List<Double> scores = pwm.computeScoreDistribution(pA, pC, pG, pT, SAMPLE_NUM);
				Collections.sort(scores);
				EmpiricalDistribution dist = new EmpiricalDistribution(binNum, scores.get(0), scores.get(scores.size() - 1));
				for (double s : scores) {
					dist.add(s);
				}
				
				oBW.write(pwm.getName());
				for(double q : quantiles) {
					oBW.write("\t");
					oBW.write(String.valueOf(Statistics.quantile(scores, q)));
				}
				oBW.newLine();
				
				BufferedWriter bw = new BufferedWriter(new FileWriter(outdir + "/" + prefix + pwm.getName() + ".dist",true));
				dist.write(bw);
				bw.close();
			}
			oBW.close();
			
		} else if("permute".equalsIgnoreCase(argMap.getTask())) {
			String pwmFile = argMap.getMandatory("pwm");
			PositionWeightMatrixIO pwmIO = new PositionWeightMatrixIO();
			FileInputStream fis = new FileInputStream(pwmFile);
			pwmIO.load(fis);
			fis.close();
			
			List<PositionWeightMatrix> pwmList = pwmIO.getMatrices();
			BufferedWriter bw = argMap.getOutputWriter();
			for( PositionWeightMatrix pwm : pwmList) {
				PositionWeightMatrix pPWM = pwm.permuteColumns(true);
				pPWM.write(bw,DEFAULT_FORMATTER);
				
			}

			bw.close();
		} else if("reverse".equalsIgnoreCase(argMap.getTask())) {
			String pwmFile = argMap.getMandatory("pwm");
			PositionWeightMatrixIO pwmIO = new PositionWeightMatrixIO();
			FileInputStream fis = new FileInputStream(pwmFile);
			pwmIO.load(fis);
			fis.close();
			
			List<PositionWeightMatrix> pwmList = pwmIO.getMatrices();
			BufferedWriter bw = argMap.getOutputWriter();
			for( PositionWeightMatrix pwm : pwmList) {
				PositionWeightMatrix pPWM = pwm.reverseComplement();
				pPWM.write(bw,DEFAULT_FORMATTER);
				
			}

			bw.close();
		}	else if("genomescan".equalsIgnoreCase(argMap.getTask())) {
			String seqDir = argMap.getMandatory("seqdir");			
			int shuffles = argMap.containsKey("shuffles") ? argMap.getInteger("shuffles") : DEFAULT_SHUFFLES;
			String prefix = argMap.containsKey("prefix") ? argMap.getMandatory("prefix") : "";
			double minScoreToWrite = argMap.getDouble("minScoreToWrite");
			String outdir = argMap.getOutputDir();
			
			String pwmFile = argMap.getMandatory("pwm");
			PositionWeightMatrixIO pwmIO = new PositionWeightMatrixIO();
			FileInputStream fis = new FileInputStream(pwmFile);
			pwmIO.load(fis);
			fis.close();
			pwmIO.addPseudoCounts();	
			List<PositionWeightMatrix> pwms = pwmIO.getMatrices();
			
			Map<String, Sequence> nonRandomChromsomes = FastaSequenceIO.loadSequencesByNameFromDirectory(new File(seqDir), true);
			
			for(Sequence chrSeq : nonRandomChromsomes.values()) {
				List<SequenceRegion> regions = chrSeq.chunk(100000, 0);
				String chr = chrSeq.getId();
				System.err.print("Scanning chromosome " + chr + " memory " + (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1000);
				System.err.print(" -sequence lodaded(" + (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1000 + ")- ");

				float gcPct = chrSeq.gcContent();
				double []  bgNucleotideFreqs = new double[4];
				bgNucleotideFreqs[0] = (1-gcPct)/(double)2;
				bgNucleotideFreqs[1] = gcPct/(double)2;
				bgNucleotideFreqs[2] = gcPct/(double)2;
				bgNucleotideFreqs[3] = (1-gcPct)/(double)2;
				
				short [] encChrSeq = chrSeq.encodeSequenceIgnoreCase();
				System.err.print(" -sequence encoded(" + (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1000 +")- ");
				System.err.println(" -sequence unloaded(" + (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1000 + ")- starting region scanning." );

				
				for(SequenceRegion annot : regions) {
					//System.err.println("\tScanning region - " + annot.toUCSC() + " memory: " +  (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1000);
					
						//System.err.println("\t\t\tAfter writting data, memory: " +  (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1000 );
					double t = System.currentTimeMillis();
					scanRegion(shuffles, prefix, minScoreToWrite, false, outdir, pwms, chr, annot, bgNucleotideFreqs, encChrSeq);
					System.err.println("Done with chunk " + annot.toUCSC() + " took: " + (System.currentTimeMillis() - t));
					//System.err.println("\t\tFinished region " + annot.toUCSC() + ", memory: " +  (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1000 );
				}
				System.err.println("Finished chromosome " + chr + ", memory: " +  (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1000 );
			}		
		}
	}

	private static double[] computeSequenceBGSequences(Map<String, Sequence> org, Map<String, List<? extends GenomicAnnotation>> regionChrMap) throws IOException {
		double [] bgFreqs = new double [4];
		int totalBases = 0;
		logger.debug("Computing sequences background nucleotide dsitribution");
		Iterator<String> chrIt = regionChrMap.keySet().iterator();
		while(chrIt.hasNext()) {
			String chr = chrIt.next();
			Sequence chrSeq = org.get(chr);
			for (GenomicAnnotation a : regionChrMap.get(chr)) {
				SequenceRegion aSR = new SequenceRegion(a.getName());
				aSR.setChromosome(chr);
				aSR.setStart(a.getStart());
				aSR.setEnd(a.getEnd());
				Sequence aSRSeq = chrSeq.getSubSequence(aSR, 0);
				float gc = aSRSeq.gcContent();
				int length = aSRSeq.getLength();
				totalBases += length;
				double totalAs = length * (1-gc)/2.0;
				double totalCs = length * gc/2.0;
				bgFreqs[0] += totalAs;
				bgFreqs[1] += totalCs;
				bgFreqs[2] += totalCs;
				bgFreqs[3] += totalAs;
				
			}
		}
		
		for(int i = 0; i < bgFreqs.length; i++) {
			bgFreqs[i] = bgFreqs[i]/(double)totalBases;
		}
		return bgFreqs;
	}
	
	public static double[] computeSequenceBGSequences(Map<String, Sequence> chrsByName, AnnotationReader<? extends GenomicAnnotation> reader) throws IOException {
		double [] bgFreqs = new double [4];
		int totalBases = 0;
		logger.debug("Computing sequences background nucleotide dsitribution");
		Iterator<String> chrIt = reader.getChromosomeIterator();
		while(chrIt.hasNext()) {
			String chr = chrIt.next();
			Sequence chrSeq = chrsByName.get(chr);
			Iterator<? extends GenomicAnnotation> aIt = reader.getChromosomeTree(chr).valueIterator();
			while (aIt.hasNext()) {
				GenomicAnnotation a = aIt.next();
				SequenceRegion aSR = new SequenceRegion(a.getName());
				aSR.setChromosome(chr);
				aSR.setStart(a.getStart());
				aSR.setEnd(a.getEnd());
				Sequence aSRSeq = chrSeq.getSubSequence(aSR, 0);
				float gc = aSRSeq.gcContent();
				int length = aSRSeq.getLength();
				totalBases += length;
				double totalAs = length * (1-gc)/2.0;
				double totalCs = length * gc/2.0;
				bgFreqs[0] += totalAs;
				bgFreqs[1] += totalCs;
				bgFreqs[2] += totalCs;
				bgFreqs[3] += totalAs;
				
			}
		}
		
		for(int i = 0; i < bgFreqs.length; i++) {
			bgFreqs[i] = bgFreqs[i]/(double)totalBases;
		}
		return bgFreqs;
	}

	public  static Map<PositionWeightMatrix, LightweightGenomicAnnotation > scanRegion(int shuffles, String prefix,
			double minScoreToWrite, boolean maxOnly, String outdir,
			List<PositionWeightMatrix> pwms, String chr,
			Annotation annot,
			double[] bgNucleotideFreqs, short[] encChrSeq) throws IOException {

		//System.err.println("\tScanning region - " + annot.toUCSC() + " memory: " +  (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1000);
		short [] annotSeq = getAnnotationSeq(encChrSeq, annot);
		Map<PositionWeightMatrix, LightweightGenomicAnnotation> bestHits = new LinkedHashMap<PositionWeightMatrix, LightweightGenomicAnnotation>(pwms.size());
		Iterator<PositionWeightMatrix> pwmIt = pwms.iterator();
		while(pwmIt.hasNext()) {
			PositionWeightMatrix pwm = pwmIt.next();
			logger.debug("\t\tScanning pwm - " + pwm.name);
			EmpiricalDistribution shuffledScoreDist = new EmpiricalDistribution(500, -50, 20);
			List<Double> maxPermutationVals = new ArrayList<Double>(shuffles);
			EmpiricalDistribution []  shuffledScoreDistArray = new EmpiricalDistribution[shuffles + 1];
			double t = System.currentTimeMillis();
			//List<Double> shuffledScores = new ArrayList<Double>(annot.length() * shuffles);
		//	System.err.println("\t\t\tBefore permutting, memory: " +  (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1000 );
			for (int i = 0; i < shuffles; i++) {
				PositionWeightMatrix shuffledPWM = pwm.permuteColumns(true);
				shuffledScoreDistArray[i + 1] = new EmpiricalDistribution(500, -50, 20);
				List<BED> hits = slidePWM(shuffledPWM, chr, annotSeq, bgNucleotideFreqs);
				double permMax = Double.NEGATIVE_INFINITY;
				for(int j = 0; j < hits.size(); j++) {
					BED hit = hits.get(j);
					double score = hit.getScore();
					//shuffledScores.add(hit.getScore());
					shuffledScoreDist.add(score);
					shuffledScoreDistArray[i + 1].add(score);
					permMax = permMax > score ? permMax : score;
				}
				maxPermutationVals.add(permMax);
			}
			logger.info("\t\tpermutations took: " + (System.currentTimeMillis() - t) +", memory: " +  (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1000 );

			Collections.sort(maxPermutationVals);
			logger.debug("\t\t\tand sorting permutation scores: " + (System.currentTimeMillis() - t) );
			//System.out.println("pwm " + pwm.getName() + " cutoff " + pwmCutoffs.get(pwm.getName()));
			t = System.currentTimeMillis(); 
			List<BED> hits = slidePWM(pwm, chr, annotSeq, bgNucleotideFreqs);
			logger.debug("\t\tSliding PWM took " + (System.currentTimeMillis() - t));
			//Add to permutated distributions
			shuffledScoreDistArray[0] = new EmpiricalDistribution(500, -50,20);
			for(BED hit : hits) {
				shuffledScoreDistArray[0].add(hit.getScore());
				shuffledScoreDist.add(hit.getScore());
				hit.setStart(hit.getStart() + annot.getStart());
				hit.setEnd(hit.getEnd() + annot.getStart());
			}
			//System.err.println("\t\t\tAfter scanning, memory: " +  (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1000 );
			BufferedWriter bw = null;
			if(outdir != null) {
				bw = new BufferedWriter(new FileWriter(outdir + "/" + prefix + pwm.getName() + ".bed",true));
			}
			BufferedWriter significanceBW = null;
			if(shuffles > 0 && outdir != null) {
				significanceBW = new BufferedWriter(new FileWriter(outdir + "/" + prefix + pwm.getName() + ".pvals",true));
				significanceBW.write("Location\tScore\tpvalue\tFWER\tFDR");
				significanceBW.newLine();
			}
			//System.err.println("\t\t\tAfter writting significance, memory: " +  (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1000 );
			try {
				BED bestHit = null;
				for(int i = 0; i <  hits.size(); i++) {
					BED hit = hits.get(i);
					//System.err.println(hit.toUCSC() + " <s="+hit.getScore());
					if(bestHit == null || hit.getScore() > bestHit.getScore()) {
						bestHit = hit;
					}
					//System.out.println("\tminScoreToReport " + minScoreToWrite + " hit " + hit.toString() + " maxOnly? " + maxOnly + " score high enough? " + (hit.getScore()  > minScoreToWrite));
					if(hit.getScore()  > minScoreToWrite && !maxOnly ) {
						if(shuffles > 0) {
							//hit.setScore(Statistics.pvalue(maxPermutationVals, hit.getScore()));
							double pval = 1 -shuffledScoreDist.getCumulativeProbability(hit.getScore());
							double fwer = Statistics.pvalue(maxPermutationVals, hit.getScore(), false);
							double fdr  = ComputeFDR.FDR(shuffledScoreDistArray[0], shuffledScoreDistArray, hit.getScore());
							if(significanceBW != null) {
								significanceBW.write(hit.toUCSC() +"\t" +hit.getScore() + "\t" + (pval) + "\t" + fwer + "\t" + fdr);
								significanceBW.newLine();
							}
							//hit.setScore(shuffledScoreDist.getCummulativeProbability(hit.getScore()));
						} 
						if( !maxOnly && bw != null) {
							bw.write(hit.toShortString());										
							bw.newLine();
						}
					}
				}
				bestHit.setName(annot.getName());
				bestHits.put(pwm,  bestHit);
				if(bw != null && maxOnly) {
					bw.write(bestHit.toShortString());										
					bw.newLine();
				}
			} finally {
				if(bw != null) {bw.close();}
				if(significanceBW != null) {significanceBW.close();}
			}
			logger.debug("\t\t\tAfter writting data, took(from start of PWM scan): " + (System.currentTimeMillis()-t)+"  memory: " +  (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory())/1000 );
		}
		return bestHits;
	}

	private static short[] getAnnotationSeq(short[] encChrSeq, Annotation annot) {
		short [] annotSeq = new short[annot.length()];
		
		for(int i =  annot.getStart() ; i < annot.getEnd(); i++) {
			annotSeq[i-annot.getStart()] = encChrSeq[i];
		}
		
		return annotSeq;
	}

	private static List<BED> slidePWM(PositionWeightMatrix pwm, String sequenceName, short[] encSeq, double [] bgFrequencies) {
		PositionWeightMatrix rpwm = pwm.reverseComplement();
		PositionWeightMatrix bg = pwm.createIsoPWM(bgFrequencies, "bg");
		int L = pwm.size();
		List<BED> scoredKmers = new ArrayList<BED>(encSeq.length - L);
		//double t = System.currentTimeMillis();
		for(int i = 0; i< encSeq.length - L; i++) {
			if(!containsGap(encSeq, i, L)){
				double neutralLod = bg.getLogLikelihood(encSeq, i);
				double directScore = pwm.getLogLikelihood(encSeq, i) - neutralLod;
				double reverseScore = rpwm.getLogLikelihood(encSeq, i) - neutralLod;
				//System.out.println("\tGot scores: " + (System.currentTimeMillis()-t));
				boolean directMatch = directScore > reverseScore;
				double score = directMatch ? directScore : reverseScore;
				BED match = new BED(null, sequenceName, i, i + L);
				//System.err.println("match " + match.toUCSC() + " neutral	" + neutralLod + " direct lod " + directScore);
				match.setScore(score);
				match.setChromosome(sequenceName);
				match.setOrientation(directMatch);
				scoredKmers.add(match);
			}
		}
		//System.out.println("\tslidePWM took: " + (System.currentTimeMillis()-t));
		return scoredKmers;
	}

	private static boolean containsGap(short[] encSeq, int i, int l) {
		for(int k = 0; k < l; k++) {
			if(encSeq[i+k] >3) {
				return true;
			}
		}
		return false;
	}
}
