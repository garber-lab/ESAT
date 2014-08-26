package broad.core.motif.misc;

import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.apache.log4j.Logger;

import broad.core.annotation.BED;
import broad.core.annotation.LightweightGenomicAnnotation;
import broad.core.error.ParseException;
import broad.core.math.ComputeFDR;
import broad.core.math.EmpiricalDistribution;
import broad.core.math.Statistics;
import broad.core.motif.PWMUtils;
import broad.core.motif.PositionWeightMatrix;
import broad.core.motif.PositionWeightMatrixIO;
import broad.pda.datastructures.Alignments;

public class ScoreMotifInstances {
	static Logger logger = Logger.getLogger(ScoreMotifInstances.class.getName());
	
	
	public ScoreMotifInstances(String pwmFile, String backgroundPWM, Alignments annotation, String genomeDirectory, String save) throws Exception{
		List<PositionWeightMatrix> pwms=getPWMs(pwmFile);
		PositionWeightMatrix bg=getPWMs(backgroundPWM).iterator().next();
		
		logger.info("Got PWM");
		
		String seq=annotation.getSequence(genomeDirectory);
		
		logger.info("Got Sequence");
		
		for(PositionWeightMatrix pwm: pwms){
			logger.info("scanning ...");
			
			List<BED> list=slidePWM(pwm, bg, annotation, seq);
			write(save, list);
		}
		
	}
	
	private List<PositionWeightMatrix> getPWMs(String pwmFile) throws ParseException, IOException {
		PositionWeightMatrixIO pwmIO = new PositionWeightMatrixIO();
		FileInputStream fis = new FileInputStream(pwmFile);
		pwmIO.load(fis);
		fis.close();
		pwmIO.addPseudoCounts();	
		List<PositionWeightMatrix> pwms = pwmIO.getMatrices();
		return pwms;
	}

	private void write(String save, List<BED> list) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(BED bed: list){
			writer.write(bed.getChr()+"\t"+bed.getStart()+"\t"+bed.getEnd()+"\t"+bed.getScore()+"\n");
		}
		
		writer.close();
	}

	private static List<BED> slidePWM(PositionWeightMatrix pwm, PositionWeightMatrix bg, Alignments annotation, String sequence) {
		PositionWeightMatrix rpwm = pwm.reverseComplement();
		int L = pwm.size();
		char[] seq=sequence.toCharArray();
		
		List<BED> scoredKmers = new ArrayList<BED>(seq.length - L);
		//double t = System.currentTimeMillis();
		for(int i = 0; i< seq.length - L; i++) {
			double neutralLod = bg.getLogLikelihood(getSeq(seq, i, i+L));
			double directScore = pwm.getLogLikelihood(getSeq(seq, i, i+L))- neutralLod;
			//double reverseScore = rpwm.getLogLikelihood(getSeq(seq, i, i+L))- neutralLod;
			//boolean directMatch = directScore > reverseScore;
			double score = Math.exp(directScore);
			BED match = new BED(null, annotation.getChr(), (annotation.getStart()+i), annotation.getStart()+(i + L));
			match.setScore(score);
			match.setChromosome(annotation.getChr());
			//match.setOrientation(directMatch);
			scoredKmers.add(match);
		}
		return scoredKmers;
	}
	
	
	private static char[] getSeq(char[] seq, int i, int j) {
		char[] rtrn=new char[j-i];
		
		for(int index=i; index<j; index++){
			rtrn[index-i]=seq[index];
		}
		
		return rtrn;
	}


	public static void main(String[] args) throws Exception{
		if(args.length>4){
			String pwm=args[0];
			String backgrund=args[1];
			Alignments region=new Alignments(args[2]);
			String genomeDirectory=args[3];
			String save=args[4];
			new ScoreMotifInstances(pwm, backgrund, region, genomeDirectory, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=pwm \n args[1]=background \n args[2]=region \n args[3]=genomeDirectory \n args[4]=save";
	
}
