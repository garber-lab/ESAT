package broad.core.siphy;

import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.apache.log4j.Logger;

import Jama.Matrix;
import broad.core.annotation.AnnotationReader;
import broad.core.annotation.BED;
import broad.core.annotation.GenomicAnnotation;
import broad.core.annotation.LightweightGenomicAnnotation;
import broad.core.error.ParseException;
import broad.core.math.Statistics;
import broad.core.motif.PWMUtils;
import broad.core.motif.PositionWeightMatrix;
import broad.core.motif.PositionWeightMatrixIO;
import broad.core.multiplealignment.MAFIO;
import broad.core.multiplealignment.MultipleAlignment;
import broad.core.multiplealignment.MultipleAlignment.AlignedSequence;
import broad.core.sequence.Sequence;
import broad.pda.seq.segmentation.ContinuousDataAlignmentModel;

public class ConservedPWMScanner {
	private List<PositionWeightMatrix> pwms;
	private int maxPWMLength;
	private MultipleAlignment currentAlignmentChunk;
	private double [] bgDistribution;
	private EvolutionaryModel model;
	private int alignmentChunkSize =  TreeScaler.MAF_CHUNK_SIZE;
	private List<String> ignoreList = new ArrayList<String>();
	static Logger logger = Logger.getLogger(ConservedPWMScanner.class.getName());
	
	public ConservedPWMScanner(String pwmFile) throws IOException, ParseException {
		PositionWeightMatrixIO pwmIO = new PositionWeightMatrixIO();
		FileInputStream fis = new FileInputStream(pwmFile);
		pwmIO.load(fis);
		fis.close();
		pwmIO.addPseudoCounts();
		//pwmIO.write(new FileOutputStream("pwm.with.pseudocounts.pwm"), NumberFormat.getNumberInstance());
		this.pwms = pwmIO.getMatrices();

		//Pre-compute things that are widely used
		Iterator<PositionWeightMatrix> pwmIt = pwms.iterator();
		while(pwmIt.hasNext()) {
			PositionWeightMatrix pwm = pwmIt.next();
			maxPWMLength = pwm.size() > maxPWMLength ? pwm.size() : maxPWMLength;
		}
	}

	public ConservedPWMScanner(List<PositionWeightMatrix> matrices) {
		this.pwms = matrices;

		//Pre-compute things that are widely used
		Iterator<PositionWeightMatrix> pwmIt = pwms.iterator();
		while(pwmIt.hasNext()) {
			PositionWeightMatrix pwm = pwmIt.next();
			maxPWMLength = pwm.size() > maxPWMLength ? pwm.size() : maxPWMLength;
		}
	}

	public Iterator<PositionWeightMatrix> getPWMiterator() {
		return pwms.iterator();
	}
	
	public void setBackgroundModel(AnnotationReader<? extends GenomicAnnotation> reader, Map<String, Sequence> org) throws IOException {
		this.bgDistribution = PWMUtils.computeSequenceBGSequences(org, reader);
	}
	
	protected List<BED> slidePWM(PositionWeightMatrix pwm, double minSeedScore, short[] encodedRef, List<int []> ungappedChunks, int numPermutations) {
		List<BED> scoredKmers = new ArrayList<BED>();
		PositionWeightMatrix rpwm = pwm.reverseComplement();
		List<PositionWeightMatrix> permPWMs = new ArrayList<PositionWeightMatrix>(numPermutations+1);
		List<PositionWeightMatrix> reversedPermPWMs = new ArrayList<PositionWeightMatrix>(numPermutations+1);
		PositionWeightMatrix bg = pwm.createIsoPWM(bgDistribution != null ? bgDistribution : model.getParameters().getBackgroundNucleotideFreqs(), "bg");
		permPWMs.add(bg);
		reversedPermPWMs.add(bg);

		PositionWeightMatrixModel pwmm = new PositionWeightMatrixModel(pwm, model);
		PositionWeightMatrixModel rpwmm = new PositionWeightMatrixModel(rpwm, model);
		AlignedSequence reference = currentAlignmentChunk.getReference();
		if(reference == null || reference.getLength() == 0) {
			return scoredKmers;
		}
		
		List<PositionWeightMatrixModel> permPWMMs = new ArrayList<PositionWeightMatrixModel>(numPermutations);
		List<PositionWeightMatrixModel> reversedPermPWMMs = new ArrayList<PositionWeightMatrixModel>(numPermutations);
		
		for(int i = 0; i < numPermutations; i++) {
			permPWMs.add(pwm.permuteColumns(false));
			reversedPermPWMs.add(permPWMs.get(i+1).reverseComplement());
			permPWMMs.add( new PositionWeightMatrixModel(permPWMs.get(i+1), model));
			reversedPermPWMMs.add( new PositionWeightMatrixModel(reversedPermPWMs.get(i+1), model));
		}
		
		
		int i = 0;
		int L = pwm.size();
		logger.trace("unggapped chunks " + ungappedChunks.size());
		for(int [] chunk : ungappedChunks) {
			for(i = chunk[0]; i< chunk[1] - L; i++) {
			//System.out.println("Start of loop: " + System.currentTimeMillis());
			//System.out.println("\tGot next sliding window: " + System.currentTimeMillis());
			//System.out.println("\tGot sequence bases: " + System.currentTimeMillis());
				double minDirectScore = Double.POSITIVE_INFINITY;
				double minReverseScore = Double.POSITIVE_INFINITY;
				double directLikelihood = pwm.getLogLikelihood(encodedRef, i);
				double reverseLikelihood = rpwm.getLogLikelihood(encodedRef, i) ;
				for(int j = 0; j < permPWMs.size(); j++) {
					minDirectScore = Math.min(minDirectScore, directLikelihood - permPWMs.get(j).getLogLikelihood(encodedRef, i));
					minReverseScore = Math.min(minReverseScore, reverseLikelihood - reversedPermPWMs.get(j).getLogLikelihood(encodedRef, i));
					
				}
				logger.trace("\tGot scores: " + System.currentTimeMillis());
				boolean directMatch = minDirectScore > minReverseScore;
				double score = directMatch ? minDirectScore : minReverseScore;
				logger.trace("\taffine score " + score + " min seed score " + minSeedScore);
				if(score >= minSeedScore) {
					int start = reference.getStart() + i;
					//System.out.println("\t\tgood hit, computing Probabilisty score " + start);
					//long startTime = System.currentTimeMillis();
					Map<String, Matrix> regionAlignment = currentAlignmentChunk.getColumnsAsVector(start, L);
					//long next = System.currentTimeMillis();
					//System.out.println("\tGot as column " + (next - startTime));
					
					double conservedScore = directMatch ? pwmm.score(regionAlignment) : rpwmm.score(regionAlignment);
					
					//long upDownTime = System.currentTimeMillis();
					//System.out.println("\tUpped & Downed: " + (upDownTime - next));
					BED match = new BED(null, reference.getChromosome(), start, start + L);
					match.setScore(conservedScore);
					match.setChromosome(reference.getChromosome());
					match.setOrientation(directMatch);
					for(int j = 0; j < permPWMMs.size(); j++) {
						//System.err.println("idx " + j + " permPWMs size " + permPWMMs.size() + " reversedPermPWMMs size " + reversedPermPWMMs.size());
						match.addExtraScore(conservedScore  - (directMatch ? permPWMMs.get(j).score(regionAlignment) : reversedPermPWMMs.get(j).score(regionAlignment)));
					}
					scoredKmers.add(match);
					//System.out.println("\tDone with monkey: " + (System.currentTimeMillis() - upDownTime));
				}
			}
			//System.out.println("\tEnd of loop: " + System.currentTimeMillis());
		}
			
		return scoredKmers;
	}

	public void setModel(EvolutionaryModel model) {
		this.model = model;
		
	}

	public int getAlignmentChunkSize() {
		return alignmentChunkSize;
	}

	public void setAlignmentChunkSize(int alignmentChunkSize) {
		this.alignmentChunkSize = alignmentChunkSize;
	}

	public Map<PositionWeightMatrix, GenomicAnnotation> scan(MAFIO chrMafIO, List<BED> annotations, int shuffles, float seedMinScore) throws IOException, ParseException {

		LinkedHashMap<PositionWeightMatrix, GenomicAnnotation> rtrnMap = new LinkedHashMap<PositionWeightMatrix, GenomicAnnotation>(); 

		for (GenomicAnnotation annot : annotations) {

			int chunkStart = annot.getStart();
			while(chunkStart < annot.getEnd()) {
	
				int chunkEnd =  shuffles == 0 ? 
						Math.min(chunkStart + getAlignmentChunkSize() + maxPWMLength - 1,annot.getEnd()) : //Overlap so one can report a hit at the end of the chunk.
							annot.getEnd();
				this.currentAlignmentChunk = ConservationUtils.setUpMAF(chrMafIO,ignoreList, model, chunkStart, chunkEnd); // could improve so that index file gets loaded only once.
				if(!currentAlignmentChunk.isEmpty()) {
					AlignedSequence reference = currentAlignmentChunk.getReference();
					logger.debug("Aligned ref " + reference.getSequenceBases());
					List<int []> ungappedChunks = reference.findUngappedSequenceChunks();
	
					short [] encodedReference = Sequence.encodeSequenceIgnoreCase(reference.getSequenceBuilder()); 
					currentAlignmentChunk.encodeAsMatrix();
					Iterator<PositionWeightMatrix> pwmIt = pwms.iterator();
					while(pwmIt.hasNext()) {
						PositionWeightMatrix pwm = pwmIt.next();
						pwm.write(new BufferedWriter(new PrintWriter(System.out)), NumberFormat.getNumberInstance());
						System.out.flush();
						long start = System.currentTimeMillis();
						double [] maxPermutationVals = new double [shuffles];
						//List<Double> shuffledScores = new ArrayList<Double>(annot.length() * shuffles);
						long permStart = System.currentTimeMillis();
						//System.out.println("pwm " + pwm.getName() + " cutoff " + pwmCutoffs.get(pwm.getName()));
						logger.debug("Sliding " + pwm.getName() + ", "+ pwm.getNumCol()+" on " + annot.toUCSC());
						List<BED> hits = slidePWM(pwm, seedMinScore, encodedReference, ungappedChunks, shuffles);
						logger.debug("got " + hits.size() + " hits");
						//Add to permutated distributions
						BED maxHit = null;
	
						for(BED hit : hits) {
							List<Double> hitShuffles = hit.getExtraScores();
							if(maxHit == null || maxHit.getScore() < hit.getScore()) {
								maxHit = new BED(hit);
							}
							for(int i = 0; i < hitShuffles.size(); i++) {
								maxPermutationVals[i] = Math.max(hitShuffles.get(i), maxPermutationVals[i]);
								//System.err.println("maxPermutation " + i + " is now " +maxPermutationVals[i] );
							}
						}
	
						if(maxHit != null) {
							Arrays.sort(maxPermutationVals);
							double fwer =  Statistics.pvalue(maxPermutationVals, maxHit.getScore(), false) ;
							maxHit.addExtraScore(fwer);
							rtrnMap.put(pwm, maxHit);
						} else {
							rtrnMap.put(pwm, null);
						}
					}
				}
				chunkStart = chunkStart + getAlignmentChunkSize();
			}

		}
		return rtrnMap;	
	}

	public void setIgnoreList(List<String> ignoreList) {
		this.ignoreList  = ignoreList;
		
	}
}
