package broad.core.siphy;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.io.RandomAccessFile;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Stack;
import java.util.TreeMap;
import org.apache.log4j.Logger;

import nextgen.core.annotation.Annotation;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;

import Jama.Matrix;
import broad.core.annotation.BED;
import broad.core.annotation.BEDReader;
import broad.core.annotation.BasicGenomicAnnotation;
import broad.core.annotation.GenomicAnnotation;
import broad.core.annotation.LightweightGenomicAnnotation;
import broad.core.error.ParseException;
import broad.core.math.ComputeFDR;
import broad.core.math.EmpiricalDistribution;
import broad.core.math.Statistics;
import broad.core.motif.PositionWeightMatrix;
import broad.core.motif.PositionWeightMatrixIO;
import broad.core.multiplealignment.MAFAlignment;
import broad.core.multiplealignment.MAFIO;
import broad.core.multiplealignment.MultipleAlignment;
import broad.core.multiplealignment.MultipleAlignment.AlignedSequence;
import broad.core.multiplealignment.MultipleAlignmentFactory;
import broad.core.multiplealignment.MultipleAlignmentIOFactory;
import broad.core.sequence.Sequence;
import broad.core.siphy.EvolutionaryModel.NodeLikelihoodParameters;
import broad.core.siphy.EvolutionaryModel.OmegaFit;
import broad.core.siphy.EvolutionaryModel.PiFit;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;

public class TreeScaler {
	
	private static final Logger logger = Logger.getLogger(TreeScaler.class.getName());
	public static final String USAGE = "Usage: TreeScaler TASK=<task_num> <task_args>\n" +
	"\tTasks:\n" +
	"\t\t1. Compute scaling of tree for each site in a multiple alignment. \n\t\tParameters:\n\t\t  -in <multiple alignment file> "+
		"\n\t\t  -mod <Neutral Evolutionary model consisting of aminoacid background distribution, mutation matrix and neutral phylogenetic tree>" +
		"\n\t\t  -format <Alignment format default is FASTA is default>" +
		"\n\t\t  -outdir <output directory> or -out <if a specific file name is desired>" +
		"\n\n\t  -ref <reference sequence id, necessary if the alignment is not in MAF format>" +
		"\n\t\t  -ignore <comma separated species to ignore>" +
		"\n\t\t  [-minTreeLength <minimum tree length to actually compute omega, sites with lesser total tree length are ignored. Default is 1>" +
		"\n\t\t  -window <Size of sliding window where tree is locally scaled, default is 1 > -windowOverlap <sliding window overlap, default is window length - 1>" +
		"\n\t\t  -neutralOmegaDist <A neutral sequence generated omega distribution by tree length,"+
		"\n\t\t\t you can use task 3 to generate this file if provided it will be use to compute p-values," +
		"\n\t\t\t the distribution should be generated with same window intended to run in this task>]" +
		"\n\t\t\t alternative, if you are sampling (see -withSampling option) you can specify: " +
		"\n\t\t [-neutralOmegas <a previously generated omega estimation run in neutral sequence with respect to the same model using full tree only>]" +
		"\n\t\t  sample evolutionary model file:" +
		"\n\t\t  ALPHABET: C,T,G,A  (if the data in the model is not using the default ordering of aminoacids)" +
		"\n\t\t  BACKGROUND: 0.214713 0.307324 0.248610 0.229353" + 
		"\n\t\t  RATE_MAT:"+
		"\n\t\t  -1.113121    0.267171    0.611653    0.234297" + 
		"\n\t\t  0.186659   -0.888845    0.199678    0.502508 " +
		"\n\t\t  0.528254    0.246835   -0.956902    0.181813 " + 
		"\n\t\t  0.219341    0.673340    0.197079   -1.089760 " + 
		"\n\t\t  TREE: ((((mm8:0.085233,rn4:0.098462):0.262242,hg18:0.128359):0.025266,canFam2:0.171487):0.308235,monDom4:0.308235);" +
		"\n\t\t  [-withSampling <Sample missing data from neutral model> -numSamplings <Number of times to run estimation to use in averaging omega. More than 10 iterations are redundant>" +
		"\n\t\t  [-start <If MAF file, you may specify the reference start coordinate> -end <If MAF file, you may specify the reference end coordinate>]" +
	"\n\t\t3. Estimate omega distribution in neutral sequence -in <Alignment file> " +
		"\n\t\t  -format <default is FASTA, MAF is also supported> " +
		"\n\t\t  -mod <Neutral Evolutionary model consisting of aminoacid background distribution, mutation matrix and neutral phylogenetic tree>" +
		"\n\t\t  -window <Window size for omega estimation>" +
		"\n\t\t  -ref <reference sequence> \n" +
		"\n\t\t  -ignore <comma separated regions to ignore>" +
	"\n\t\t4. Write model -kappa <HKY kappa parameter> - background <comma separated background nucleotide frequencies (A,C,G,T)> [-mu <Scaling factor default is 1>]" +
	"\n\t\t5. Generate sampled alignment -mod <Neutral Evolutionary model consisting of nucleotide background distribution, mutation matrix and neutral phylogenetic tree>" +
		"\n\t\t  -format <Alignment format default is FASTA >" +
		"\n\t\t  -out <output file>" +
		"\n\t\t  -colNum <Number of columns to generate>" +
	"\n\t\t6. Fill alignment gaps and missing data with neutral model generated probabilities." +
	"\n\t\t -in <Alignment file>" +
	"\n\t\t -out <Output file name> " +
	"\n\t\t -mod <Neutral model as defined in task 1>" +
	"\n\t\t [-informat <input alignment format, default is FASTA> -ignore <comma separated species to ignore> -outformat <output format default is FASTA]" +
	"\n\t\t7. Fit base frequency pi. Parameters are the same as task 1 except for:" +
	"\n\t\t one off -priorOmega <a prior omega to use when fitting rhow> "+ 
	"\n\t\t\t -siteOmegas <An output of the task 1 or 2 for the same region omega will then be avaraged in windows> -omegaWindow <windo in which to fit omega, default is 1kb> " +
	"\n\t\t8. Integrate omegas over windows, this assumes that the omegas provided have been computed on the same tree length (e.g. via sampling) " +
	"\n\t\t\t -in <A site by site omega calculation output file> -out <output file> -shift <A coordinate shift in case the original omega file is not in genomic coordinates> " +
	"\n\t\t\t -chr <chromosome the data belongs to> -window <window size to integrate over> -neturalOmegas <A site by site omega calculation output on neutral sequence>"  +
	"\n\t\t9. Define conserved elements based on a desired score, this works for any format provided you specify the column containing the score" +
	"\n\t\t\t -in <estimation file (any file generated from running tasks 1,2,7,8 or similar), default is standard input> -out <output file name defalut is standard out> [-separator <Column separator character, default is tab>]" +
	"\n\t\t\t -pvalcol <column of site p-value (first column is 1)> -positioncol <column of position, position is assumed to be the start of the window(first column is 1)> -minscore <Minimum score used to join a site> -maxscore <Maximum score to join> -maxgap <maximum gap of sites with less than pval to paste through>" +
	"\n\t\t\t [-chr <chromosome default is 'C'> -window <window size default is 1> " +
	"\n\t\t10. Score regions -alignment <alignment including regions, regions should be in coordinates that are consistent with the alignment> -mod <model file> -ref <reference sequence id default is the first sequence in the alignment> [-window <If you want to tile each region with a fixed window rather than fitting omega to the full region> -overlap <By default a scan of windows with window size - 1 overlap is done> -minTreeLength <If no omega should be computed if the minimum branch length of the kmer is below this threshold>]" +
	"\n\t\t\t -in <Input file BED annotation file, default is standard input> -out <Output file, default is standard output> [-shift <amount to shift position> -bedIsOneBased <Add this flad is the positions in BED file start at 1 rather than 0 ]" +
	"\n\t\t11. Integrate Stationary distribution in windows -window <window size> [-windowOverlap <sliding window overlap, default is window length - 1> -dist <Neutral distribution in empirical format see??? if pvalues are desired>]" +
	"\n\t\t\t -in <estimation file (any file generated from running tasks 1,2,7,8 or similar), default is standard input> -out <output file name defalut is standard out>" +
	"\n\t\t12. Add an empirical pvalue to computation based on the spcified column value -col <column number with values to use, first column is 0> -dist <Neutral distribution of same values to use in empirical pValue computations> " +
	"\n\t\t\t-in <File (standard input is default)> -out <new file with the extra column of > [-separator <character used to separate columns, default is tab '\\t>' -rightTail <add this flag if the pvalue should be 1- p(X>=x) rather than p(X<= x)>]" +
	"\n\t\t13. Integrate Omega log odds score distribution in windows -window <window size> [-dist <Neutral distribution in empirical format see??? if pvalues are desired>]" +
	"\n\t\t14. Given a PWM, an alignment and a neutral model slide PWM and compute log odds likelihood (of the window being generated by the neutral or PWM models " +
	"\n\t\t\t-in <Alignment file (or standard input) if MAF this is required> -pwm <File with PWM description> -mod <Neutral model> [-seedMinScore <Minimum affinity score in order to incurr in the expense of the phylogenetic computation> -outprefix <If a prefix to the automatically generated output file is desired>]" +
	"\n\t\t\t[-ref <if Alignment not in MAF format> -alignmentFormat <one of [FASTA],PHYLIP. Fasta is default> -start <If alignment is in MAF format> -end <If alignment is in MAF format> -outdir <Output directory if other than current dir> ]" +
	"\n\t\t15. Compute sample scores for a set of PWMs in <Alignment file (or standard input) if MAF this is required NOTE: Only MAF supported at the moment> -pwm <File with PWM description> -mod <Neutral model> -numOfSamples <Number of samples in which to compute scores> " +
	"\n\t\t\t[-ref <if Alignment not in MAF format> -alignmentFormat <one of [FASTA],PHYLIP. Fasta is default> -start <If alignment is in MAF format> -end <If alignment is in MAF format> -outdir <Output directory if other than current dir> ]" +
	"\n\t\t16. Given a PWM, an MAF alignment and a neutral model slide PWM and compute log odds likelihood (of the window being generated by the neutral or PWM models " +
	"\n\t\t\t-indir <Alignment directory of chromosome  MAF alignments> -pwm <File with PWM description> -mod <Neutral model> [-seedMinScore <Minimum affinity score in order to incurr in the expense of the phylogenetic computation> -ignore <comma separated species to ignore>]" +
	"\n\t\t\t-start <start of region > -end <end of region > -chr <chromosome region> -regions <alternatively you can specify an annotation file > -regionFormat <[BED], GFF or generic> " +
	"\n\t\t\t[-mafSuffix <A suffix for maf alignment files default is .maf> -outdir <Output directory if other than current dir> -outprefix <If a prefix to the automatically generated output file is desired> -minScoreToReport <Do not report scores less than this> -shuffles <Number of shuffles to do if suffling then a pvalue is reported>>]" +
	"\n\t\t17. Change equilibrium distribution in model. -mod <Model file to change (see task 1 for a description> -pA <new equilibrium probability for A>  -pC <new equilibrium probability for C>  -pG <new equilibrium probability for G>  -pT <new equilibrium probability for T> -out <Output file or standard out if none is specified>"+
	"\n\t\t18. Evolve sequence according to tree -ancestralSequence <A nucleotide sequence to evolve> -mod <Model to use> -ignore <optional -- comma separated species to ignore in the given model tree> -numColumns <Number of columns to sample> -bg <Optional -- new background distribution, as a comma separated list of the A,C,G,T frequencies> -" +
	"\n\t\tbayesian. Estimate posterior P(omega | Data). Basic data (Alignment and model should be specified per in task 1) specific parameters: " +
	"\n\t\t\t-out <Name of output file containig P(0.25 | data) and P(1 | data) for each window, another file will also be created with the aggregated empiric distribution of P(W | Data)> " +
	"\n\t\t\t-printFullDistribution <If set the program will print the posterior probability for the sampled values of omega for each position of the alignment. THIS GENERATES A HUGE FILE>" +
	"\n\t\t\t-likelihood  Computed the probability of an alignment given a model.\n\t\t -mod <Neutral model as defined in task 1> \n\t\t -in <Alignment file> \n\n\t -ref <reference sequence id, necessary if the alignment is not in MAF format> \n\t\t -ignore <comma separated species to ignore>"+
	"\n\t\tmaximalPWM Given a PWM, an MAF alignment and a neutral model slide PWM and compute the maximum log odds likelihood (of the window being generated by the neutral or PWM models " +
	"\n\t\t\t-indir <Alignment directory of chromosome  MAF alignments> -pwm <File with PWM description> -mod <Neutral model> [-seedMinScore <Minimum affinity score in order to incurr in the expense of the phylogenetic computation> -ignore <comma separated species to ignore>]" +
	"\n\t\t\t-regions <annotation file to score > -regionFormat <[BED], GFF or generic> " +
	"\n\t\t\t[-mafSuffix <A suffix for maf alignment files default is .maf> -out <Output file or standard out if non is specified>  -shuffles <Number of shuffles to do if suffling then a pvalue is reported>>]" +
	"\n";	

	EvolutionaryModel model;
	private MultipleAlignment alignment;
	private double minimumTreeLength;
	private HashMap<String, Phylogeny> prunnedTrees = new HashMap<String, Phylogeny>();
	private List<String> ignoreSequences;
	
	private static final double MIN_TREE_LENGTH = 1;
	private static final int INF = 1000000000;
	private static final int DEFAULT_OMEGA_WIN_FOR_RHO = 2000;
	private static final int DEFAULT_SAMPLINGS = 8;
	private static final double MIN_INTERESTING_OMEGA = 0.4;
	private static DecimalFormat numberFormat = new  DecimalFormat("##0.####");
	private static DecimalFormat tinnyNumberFormat = new  DecimalFormat("##0.#########");
	static int MAF_CHUNK_SIZE = 100000;
	private String chr;
	
	public TreeScaler() {
		super();
		minimumTreeLength = MIN_TREE_LENGTH;
	}
	
	
	public TreeScaler(Annotation region, File modelFile, String alnFile, String alnFileFormat)throws Exception{
		super();
		minimumTreeLength = MIN_TREE_LENGTH;
		setNeutralModel(modelFile);
		MultipleAlignment alignment=ConservationUtils.setUpMAF(alnFile, new ArrayList<String>(), getModel(), region.getStart(), region.getEnd());
		//System.err.println("Neutral model");
		setAlignment(alignment);
		encodeAsMatrix();
		this.chr=region.getChr();
	}
	
	public void setAlignment(String alnFile, String alnFileFormat)throws Exception{
		alignment = MultipleAlignmentFactory.create(alnFile, alnFileFormat);
	}
	
	public void encodeAsMatrix(){
		alignment.encodeAsMatrix();
	}
	
	public static void main(String[] args) throws Exception {
		logger.debug("Logger level set to DEBUG or less");
		
		ArgumentMap argMap = CLUtil.getParameters(args, USAGE);
		
		TreeScaler scaler = new TreeScaler();
		if ("1".equals(argMap.getTask())) {	
			File modelFile = new File(argMap.getMandatory("mod"));
			String alnFile = argMap.getInput();
			String alnFileFormat = argMap.containsKey("format") ? argMap.get("format") : "FASTA";
			int window = argMap.containsKey("window") ? argMap.getInteger("window") : 1;
			int overlap = argMap.containsKey("windowOverlap") ? argMap.getInteger("windowOverlap") : window - 1;
			String out = null;
			String outdir = null;
			if(argMap.isOutputSet()) {
				out = argMap.getOutput();
			} else {
				outdir = argMap.getOutputDir();
			}
			String ignoreListStr = argMap.get("ignore");
			boolean sample = argMap.containsKey("withSampling");
			int numSamplings = argMap.containsKey("numSamplings") ? argMap.getInteger("numSamplings") : DEFAULT_SAMPLINGS;
			double minTreeLength = argMap.containsKey("minTreeLength") ? argMap.getDouble("minTreeLength") : MIN_TREE_LENGTH;

			scaler.setNeutralModel(modelFile);
			scaler.setMinimumTreeLength(minTreeLength);	
			if(argMap.isPresent("neutralOmegaDist")) {
				scaler.model.setOmegaDistByTreeLength(argMap.get("neutralOmegaDist"));
			} if(argMap.isPresent("neutralOmegas")) {
				System.out.println("Setting neutral omega sitribution ");
				scaler.model.setNeutralOmegaDistribution(argMap.get("neutralOmegas"));
				System.out.println("done setting neutral omega stats: ");// + scaler.model.binStatistics);
			}
			
			List<String> ignoreList = processIgnoreListString(ignoreListStr);
			
			scaler.setUpAlignment(argMap, alnFile, alnFileFormat, ignoreList); 

			//scaler.alignment = MultipleAlignmentFactory.create(alnFile, alnFileFormat);
			//scaler.alignment.remove(ignoreList);
			//scaler.alignment.setReferenceId(ref);			
			//scaler.ignoreSequences = ignoreList; 
			scaler.alignment.encodeAsMatrix();
			
			String [] alnFilePath = alnFile.split("/");
			
			if(out == null) {
				out = outdir + "/" + alnFilePath[alnFilePath.length - 1].replaceFirst("\\..+$", ".omegas");
			}
			
			BufferedWriter bw = new BufferedWriter(new FileWriter(out));
			if(sample) {
				System.out.println("Doing " + numSamplings + " omega samplings per window");
				scaler.scaleTreeWithSampling(window, bw, ignoreList, numSamplings);
			} else {
				scaler.scaleTree(window, bw, ignoreList, overlap);
			}
			bw.close();
			/*
			scaler.alignment.setIOHelper(MultipleAlignmentIOFactory.create("PHYLIP"));
			bw = new BufferedWriter(new FileWriter(alnFile + ".sampled"));
			scaler.alignment.write(bw);
			bw.close();
			*/
		} else if("3".equals(argMap.getTask())) {
			File modelFile = new File(argMap.getMandatory("mod"));
			String alnFile = argMap.getInput();
			String alnFileFormat = argMap.containsKey("format") ? argMap.get("format") : "FASTA";
			int window = argMap.containsKey("window") ? argMap.getInteger("window") : 1;
			String ref = argMap.getMandatory("ref");
			String ignoreListStr = argMap.get("ignore");

			scaler.setNeutralModel(modelFile);
			
			String [] alnFilePath = alnFile.split("/");
			String omegaDistFile = alnFilePath[alnFilePath.length - 1].replaceFirst("\\..+$", "_omegas.dist");
			String omegasByPos   = alnFilePath[alnFilePath.length - 1].replaceFirst("\\..+$", "_omegas_by_pos.dist");
			
			List<String> ignoreList = processIgnoreListString(ignoreListStr);			
			
			Iterator<? extends MultipleAlignment> it = null;
			if("MAF".equals(alnFileFormat)) {
				String mafIndex = alnFile + ".index";
				MAFAlignment mafAln = new MAFAlignment( mafIndex);
				RandomAccessFile alnRaf = new RandomAccessFile(alnFile , "r");
				mafAln.load(alnRaf,  ignoreList);
				alnRaf.close();
				throw new Exception("Need to reimplement after reimplementing the MAF alignment class");
				//it = mafAln.getAlignmentBlocks().iterator();

			} else {
				List<MultipleAlignment> oneAlignmentList = new ArrayList<MultipleAlignment>(1);
				MultipleAlignment ma =  MultipleAlignmentFactory.create(alnFile, alnFileFormat);
				ma.remove(ignoreList);
				ma.encode();
				ma.setReferenceId(ref);
				oneAlignmentList.add(ma);
				it = oneAlignmentList.iterator();
			}
			

			BufferedWriter pbw = new BufferedWriter(new FileWriter(omegasByPos));
			Map<Double, List<Double>> treeLengthsOmegaDists = new HashMap<Double, List<Double>>();
			try {
				while(it.hasNext()) {
					MultipleAlignment aln = it.next();
					aln.encode();
					scaler.alignment = aln;
					//aln.setIOHelper(MultipleAlignmentIOFactory.create("PHYLIP"));
					//aln.write(bw);
					//bw.flush();
					Map<Integer, Map<Double, Double>> positionTreeScalings = 
						scaler.scaleTreeWithRandomPrunnings(window);
					Iterator<Integer> posIt = positionTreeScalings.keySet().iterator();
					while(posIt.hasNext()) {
						int pos = posIt.next();
						pbw.write(String.valueOf(pos));
						Map<Double, Double>  posScalings= positionTreeScalings.get(pos);
						Iterator<Double> treeLengthIt = posScalings.keySet().iterator();
						while(treeLengthIt.hasNext()) {
							double treeLength = treeLengthIt.next();
							double omega      = posScalings.get(treeLength);
							pbw.write("\t");
							pbw.write(String.valueOf(treeLength));
							pbw.write("-");
							pbw.write(String.valueOf(omega));
							List<Double> treeLengthOmegaDist = treeLengthsOmegaDists.get(treeLength);
							if(treeLengthOmegaDist == null) {
								treeLengthOmegaDist = new ArrayList<Double>();
								treeLengthsOmegaDists.put(treeLength, treeLengthOmegaDist);
							}
							
							treeLengthOmegaDist.add(omega);
						}
						pbw.newLine();
					}
				}

			} finally {
				pbw.close();
			}
			

			List<Double> treeLengthDist = new ArrayList<Double>(treeLengthsOmegaDists.keySet());
			Collections.sort(treeLengthDist);
			
			BufferedWriter dbw = new BufferedWriter(new FileWriter(omegaDistFile));
			try {
				Iterator<Double> treeLengthIt = treeLengthDist.iterator();
				while(treeLengthIt.hasNext()) {
					double treeLength = treeLengthIt.next();
					dbw.write(String.valueOf(treeLength));
					dbw.write("\t");
					List<Double> omegaDist = treeLengthsOmegaDists.get(treeLength);
					Collections.sort(omegaDist);
					Iterator<Double> omegaDistIt = omegaDist.iterator();
					while(omegaDistIt.hasNext()) {
						dbw.write(String.valueOf(omegaDistIt.next()));
						if(omegaDistIt.hasNext()) {
							dbw.write(",");
						}
					}
					dbw.newLine();
				}
			} finally {
				dbw.close();
			}
		} else if("4".equals(argMap.getTask())) {
			double kappa = argMap.getDouble("kappa");
			
			String [] bgStr = argMap.getMandatory("background").split(",");
			if(bgStr.length != 4) {
				System.err.println("Brackground frequencies must be 4! in the form of Afreq,Cfreq,Gfreq,Tfreq");
				return;
			}
			double [] bg = new double[bgStr.length];
			for(int i = 0; i < bgStr.length; i++) {
				bg[i] = Double.parseDouble(bgStr[i]);
			}
			double mu = argMap.containsKey("mu") ? argMap.getDouble("mu") : 1;
			
			EvolutionaryModelParameters emp = new EvolutionaryModelParameters(kappa, bg, mu);
			
			BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(System.out));
			emp.write(bw);
			bw.flush();
		}else if("5".equals(argMap.getTask())) {
			File modelFile = new File(argMap.getMandatory("mod"));
			String out = argMap.getOutput();
			String format = argMap.containsKey("format") ? argMap.get("format") : "FASTA";
			int colNums = argMap.getInteger("colNum");
			String ignoreListStr = argMap.get("ignore");
			List<String> ignoreList = processIgnoreListString(ignoreListStr);
			scaler.setNeutralModel(modelFile);
			Phylogeny prunedTree = ConservationUtils.pruneTree(ignoreList, scaler.model.getTree());
			scaler.model.getModelParameters().setTree(prunedTree);
			scaler.generateNeutralAlignment(colNums, format, out);
			
		} else if("6".equals(argMap.getTask())) {
			File modelFile = new File(argMap.getMandatory("mod"));
			String out = argMap.getOutput();
			String informat = argMap.containsKey("format") ? argMap.get("format") : "FASTA";
			String outformat = argMap.containsKey("outformat") ? argMap.get("outformat") : "FASTA";
			String in = argMap.getInput();
			double minTreeLength = argMap.containsKey("minTreeLength") ? argMap.getDouble("minTreeLength") : MIN_TREE_LENGTH;
			
			List<String> ignoreList = processIgnoreListString(argMap.get("ignore"));

			scaler.setNeutralModel(modelFile);
			
			scaler.alignment = MultipleAlignmentFactory.create(in, informat);
			scaler.alignment.remove(ignoreList);
			scaler.ignoreSequences = ignoreList;
			
			scaler.fillInAlignmentMissingData(minTreeLength, true);
			
			scaler.alignment.setIOHelper(MultipleAlignmentIOFactory.create(outformat));
			
			BufferedWriter bw = new BufferedWriter(new FileWriter(out));
			scaler.alignment.write(bw);
			bw.close();
			
		} else if ("7".equals(argMap.getTask())) {	
			File modelFile = new File(argMap.getMandatory("mod"));
			String alnFile = argMap.getInput();
			String alnFileFormat = argMap.containsKey("format") ? argMap.get("format") : "FASTA";
			String out = null;
			String outdir = null;
			if(argMap.isOutputSet()) {
				out = argMap.getOutput();
			} else {
				outdir = argMap.getOutputDir();
			}
			String ignoreListStr = argMap.get("ignore");
			boolean sample = argMap.containsKey("withSampling");
			double minTreeLength = argMap.containsKey("minTreeLength") ? argMap.getDouble("minTreeLength") : MIN_TREE_LENGTH;
			double priorOmega = argMap.containsKey("priorOmega") ? argMap.getDouble("priorOmega") : 1;
			
			scaler.setNeutralModel(modelFile);
			scaler.setMinimumTreeLength(minTreeLength);	
			scaler.model.setOmega(priorOmega);

			List<String> ignoreList = processIgnoreListString(ignoreListStr);
			
			scaler.setUpAlignment(argMap, alnFile, alnFileFormat, ignoreList); 
			
			if(sample) {
				System.out.print("Filling for missing data ... ");
				scaler.fillInAlignmentMissingData(minTreeLength, true);
				System.out.println("done");
			} else {
				scaler.alignment.encodeAsMatrix();
			}
			
			String [] alnFilePath = alnFile.split("/");
			
			if(out == null) {
				out = outdir + "/" + alnFilePath[alnFilePath.length - 1].replaceFirst("\\..+$", ".pi");
			}
			
			BufferedWriter bw = new BufferedWriter(new FileWriter(out));
			scaler.fitPI(bw, ignoreList);
			bw.close();
			/*
			scaler.alignment.setIOHelper(MultipleAlignmentIOFactory.create("PHYLIP"));
			bw = new BufferedWriter(new FileWriter(alnFile + ".sampled"));
			scaler.alignment.write(bw);
			bw.close();
			*/
		} else if("8".equals(argMap.getTask())) {
			int window = argMap.getInteger("window");
			String neutralDistFile = argMap.get("neutralDist");
			String chr = argMap.get("chr") == null ? "C" : argMap.get("chr");
			int overlap = argMap.containsKey("windowOverlap") ? argMap.getInteger("windowOverlap") : window - 1;
			
			EmpiricalDistribution neutralDist = null;
			if(neutralDist != null) {
				neutralDist = new EmpiricalDistribution(new File(neutralDistFile));
			}
			
			TreeScalerIO omegaio = new TreeScalerIO();
			
			InputStream is = argMap.getInputStream();
			omegaio.load(is);
			is.close();
			
			List<BED> windows = integrateStationaryDistributions(omegaio.getScaledWindows(), window, chr, neutralDist, overlap);
			Iterator<BED> windowIt = windows.iterator();
			BufferedWriter bw = argMap.getOutputWriter();
			while(windowIt.hasNext()) {
				bw.write(windowIt.next().toString());
				bw.newLine();
			}
			bw.close();
			
		} else if ("9".equals(argMap.getTask())) {
			int scorecol      = argMap.getInteger("scorecol");
			int positioncol  = argMap.getInteger("positioncol");
			double maxScore   = argMap.containsKey("maxscore") ? argMap.getDouble("maxscore") : INF;
			double minScore   = argMap.containsKey("minscore") ? argMap.getDouble("minscore") : 0;
			int maxgap       = argMap.getInteger("maxgap");

			String separator = argMap.containsKey("separator") ? argMap.get("separator") : "\t";
			String chr       = argMap.containsKey("chr") ? argMap.get("chr") : "C";
			//int shift        = argMap.containsKey("shift") ? argMap.getInteger("shift") : 0;
			int window       = argMap.containsKey("window") ? argMap.getInteger("window") : 1;
 			
			BufferedReader br = argMap.getInputReader();
			LinkedHashMap<Integer, Double> siteScoreMap = new LinkedHashMap<Integer, Double>();
			String line = null;
			while ( (line = br.readLine()) != null) {
				if(line.trim().length() == 0 || line.startsWith("#")) {
					continue;
				}
				
				String [] lineInfo = line.split(separator);
				int position = Integer.parseInt(lineInfo[positioncol - 1]);
				double score  = Double.parseDouble(lineInfo[scorecol - 1]);
				siteScoreMap.put(position, score);
			}
			br.close();
			
			List<BED> elements = threadElements(siteScoreMap, maxScore, minScore, maxgap, chr, window);  
			BufferedWriter bw = argMap.getOutputWriter();
			Iterator<BED> elementIt = elements.iterator();
			while(elementIt.hasNext()) {
				bw.write(elementIt.next().toString());
				bw.newLine();
			}
			bw.close();
			
		} else if ("10".equals(argMap.getTask()) ){
			String alnFile = argMap.getMandatory("alignment");
			File modelFile = new File(argMap.getMandatory("mod"));
			String alnFileFormat = argMap.containsKey("format") ? argMap.get("format") : "FASTA";
			List<String> ignoreList = processIgnoreListString(argMap.get("ignore"));
			int shift = argMap.containsKey("shift") ? argMap.getInteger("shift") : 0;
			boolean oneBased = argMap.containsKey("bedIsOneBased");
			boolean intervalsClosedClosed = argMap.containsKey("intervalsClosedClosed");
			String chr = argMap.get("chr");
			
			scaler.setNeutralModel(modelFile);
			MAFIO  chrMafIO = null;
			if("MAF".equalsIgnoreCase(alnFileFormat)){ 
				chrMafIO = new MAFIO(alnFile, true);
			} else {
				scaler.alignment = ConservationUtils.setUpAlignment(argMap, alnFile, alnFileFormat, ignoreList, scaler.model);
				scaler.alignment.encodeAsMatrix();
				scaler.alignment.remove(ignoreList); //TODO: It is not elegant/clear to remove sequences here for non MAF alignmets while doing so later for MAFs
			}

			BufferedReader ir = argMap.getInputReader();
			BEDReader reader = new BEDReader(ir);
			ir.close();
			Iterator<BED> regionIt = null;
			int regionNum = 0;
			if(chr == null || chr.trim().length() == 0) {
				List<BED> annotationList = reader.getAnnotationList();
				regionNum = annotationList.size();
				regionIt = annotationList.iterator();
			} else {
				List<BED> chrAnnotationList = reader.getChromosomeBEDs(chr);
				regionIt = chrAnnotationList.iterator();
				regionNum = chrAnnotationList.size();
			}
			System.err.println("Estimating  " + regionNum + " regions");
			BufferedWriter bw = argMap.getOutputWriter();
			AlignedSequence refSequence = null;
			if(!"MAF".equalsIgnoreCase(alnFileFormat)) {
				refSequence = scaler.alignment.getReference();
			}
			//System.err.println("Reference: " + refSequence.getId() + " "  + refSequence.getStart() + "-" + refSequence.getEnd());
			while(regionIt.hasNext()) {
				BED region = regionIt.next();
				//System.err.print("Region " + region.getName() + " chr" + region.getChromosome() + ":" + region.getStart() + "-" + region.getEnd() );
				if(oneBased) {
					region.setStart(region.getStart() - 1);
					region.setEnd(region.getEnd() - 1);
				}
				if(intervalsClosedClosed) {
					region.setEnd(region.getEnd() + 1);
				}
				if("MAF".equalsIgnoreCase(alnFileFormat)) {
					scaler.alignment = ConservationUtils.setUpMAF(chrMafIO,ignoreList, scaler.getModel(), region.getStart(), region.getEnd()); // could improve so that index file gets loaded only once.
					scaler.alignment.remove(ignoreList); 
					scaler.alignment.encodeAsMatrix();
				} else if(! (refSequence.getStart() <= region.getStart() && refSequence.getEnd() > region.getEnd())) {
					//System.err.print( " was not within alignment boundaries.");
					region = new BED(region.intersect(refSequence)); //Changed by Manuel to account for new Annotation logic
					//System.err.print(" -- Interesected it with reference: chr" + region.getChromosome() + ":" + region.getStart() + "-" + region.getEnd() );
					if(region.getLength() == 0) {
						//System.err.print(" but got empty intersection => Skipped");
						continue;
					}
				
				}
				
				if( argMap.containsKey("window")) {
					int window = argMap.getInteger("window");
					double minTreeLength = argMap.containsKey("minTreeLength") ? argMap.getDouble("minTreeLength") : MIN_TREE_LENGTH;
					scaler.setMinimumTreeLength(minTreeLength);
					int overlap = argMap.containsKey("overlap") ? argMap.getInteger("overlap") : window - 1;
					scaler.scaleTree(window, bw, ignoreList, overlap);
				} else {
					OmegaFit regionFit = scaler.scaleRegion(ignoreList, region);
					
					bw.write(region.getChromosomeString() + "\t" + (region.getStart() + shift) + "\t" + (region.getEnd() + shift) + "\t" + region.getName() + "\t" + region.getOrientation() + "\t" + (regionFit == null ? 1 : regionFit.getOmega() )+ "\t" + (regionFit == null ? 0 : regionFit.getLogOddsScore() )+ "\t" + (regionFit == null ? 1 :regionFit.getPVal() )+ "\t" + (regionFit == null ? 0 :regionFit.getTreeLength()));
					bw.newLine();
					//System.err.print(" wrote fit");
				}
				//System.err.println("");
				//scaler.writeSiteOmegaInfo(bw, region.getStart(), regionFit);
			}
			bw.close();
		} else if ("11".equals(argMap.getTask())) {
			int window = argMap.getInteger("window");
			String neutralDistFile = argMap.get("neutralDist");
			String chr = argMap.get("chr") == null ? "C" : argMap.get("chr");
			int overlap = argMap.containsKey("windowOverlap") ? argMap.getInteger("windowOverlap") : window - 1;
			
			EmpiricalDistribution neutralDist = null;
			if(neutralDist != null) {
				neutralDist = new EmpiricalDistribution(new File(neutralDistFile));
			}
			
			StationaryDistributionIO sdio = new StationaryDistributionIO();
			
			InputStream is = argMap.getInputStream();
			sdio.load(is);
			is.close();
			
			BufferedWriter bw = argMap.getOutputWriter();
			//List<BED> windows = integrateStationaryDistributions(sdio.getFits(), window, chr, neutralDist, overlap);
			integrateStationaryDistributions(sdio.getFits(), window, chr, neutralDist, overlap, bw);
			/*Iterator<BED> windowIt = windows.iterator();
			BufferedWriter bw = argMap.getOutputWriter();
			while(windowIt.hasNext()) {
				bw.write(windowIt.next().toString());
				bw.newLine();
			}
			*/
			bw.close();
		} else if ("12".equals(argMap.getTask())) {
			String distFile = argMap.getMandatory("dist");
			int col = argMap.getInteger("col");
			String separator = argMap.containsKey("separator") ? argMap.get("separator") : "\t";
			boolean rightTail = argMap.containsKey("rightTail");
			
			
			EmpiricalDistribution dist = new EmpiricalDistribution(new File(distFile));
			
			BufferedReader br = argMap.getInputReader();
			BufferedWriter bw = argMap.getOutputWriter();
			
			String line = null;
			while((line = br.readLine()) != null) {
				if(line.startsWith("#") || line.startsWith("track")) {
					bw.write(line);
					bw.newLine();
					continue;
				}
				
				String [] lineInfo = line.split(separator);
				double value = Double.parseDouble(lineInfo[col]);
				double pVal = rightTail ? 1 - dist.getCumulativeProbability(value) : dist.getCumulativeProbability(value);
				bw.write(line.trim());
				bw.write(separator);
				bw.write(tinnyNumberFormat.format(pVal));
				bw.newLine();
			}
			br.close();
			bw.close();
		} else if ("13".equals(argMap.getTask())) {
			int window = argMap.getInteger("window");
			String neutralDistFile = argMap.get("neutralDist");
			String chr = argMap.get("chr") == null ? "C" : argMap.get("chr");
			int overlap = argMap.containsKey("windowOverlap") ? argMap.getInteger("windowOverlap") : window - 1;
			EmpiricalDistribution neutralDist = null;
			if(neutralDist != null) {
				neutralDist = new EmpiricalDistribution(new File(neutralDistFile));
			}
			
			TreeScalerIO tsio = new TreeScalerIO();
			
			InputStream is = argMap.getInputStream();
			tsio.load(is);
			is.close();
			
			List<BED> windows = integrateStationaryDistributions(tsio.getScaledWindows(), window, chr, neutralDist, overlap);
			Iterator<BED> windowIt = windows.iterator();
			BufferedWriter bw = argMap.getOutputWriter();
			while(windowIt.hasNext()) {
				bw.write(windowIt.next().toString());
				bw.newLine();
			}
			bw.close();
		} else if ("14".equals(argMap.getTask())) {
			File modelFile = new File(argMap.getMandatory("mod"));
			scaler.setNeutralModel(modelFile);
			String alnFileFormat = argMap.containsKey("alignmentFormat") ? argMap.get("alignmentFormat") : "FASTA";
			String alnFile = argMap.getInput();
			String outdir = argMap.getOutputDir();
			float seedMinScore = argMap.containsKey("seedMinScore") ? argMap.getFloat("seedMinScore") : 0;
			int numPermutations = argMap.containsKey("numPermutations") ? argMap.getInteger("numPermutations") : 5;
			String prefix = argMap.containsKey("outprefix") ? argMap.get("outprefix") : "";
			System.err.println("... loading alignment");
			scaler.setUpAlignment(argMap, alnFile, alnFileFormat,new ArrayList<String>());
			System.err.println("Finished loading Alignment, length: " + scaler.alignment.length());
			String pwmFile = argMap.getMandatory("pwm");
			PositionWeightMatrixIO pwmIO = new PositionWeightMatrixIO();
			FileInputStream fis = new FileInputStream(pwmFile);
			pwmIO.load(fis);
			fis.close();
			List<PositionWeightMatrix> pwms = pwmIO.getMatrices();
			AlignedSequence reference = scaler.alignment.getReference();
			List<int []> ungappedChunks = reference.findUngappedSequenceChunks();
			short [] encodedReference = Sequence.encodeSequenceIgnoreCase(reference.getSequenceBuilder());
			scaler.alignment.encodeAsMatrix();
			Iterator<PositionWeightMatrix> pwmIt = pwms.iterator();
			while(pwmIt.hasNext()) {
				PositionWeightMatrix pwm = pwmIt.next();
				List<BED> hits = scaler.slidePWM(pwm, seedMinScore, encodedReference, ungappedChunks, numPermutations);
				Iterator<BED> hitIt = hits.iterator();
				BufferedWriter bw = new BufferedWriter(new FileWriter(outdir + "/" + prefix  + pwm.getName() + ".bed"));
				try {
					while(hitIt.hasNext()) {
						BED hit = hitIt.next();
						bw.write(hit.toString(false));
						bw.newLine();
					}
				} finally {
					if(bw != null) {bw.close();}
				}
			}

			
		}		
		else if ("15".equals(argMap.getTask())) {
			File modelFile = new File(argMap.getMandatory("mod"));
			scaler.setNeutralModel(modelFile);
			//String alnFileFormat = argMap.containsKey("alignmentFormat") ? argMap.get("alignmentFormat") : "FASTA";
			String alnFile = argMap.getInput();
			int start = argMap.getInteger("start");
			int end   = argMap.getInteger("end");
			float numOfSamples = argMap.getInteger("numOfSamples");
			
			//scaler.setUpAlignment(argMap, alnFile, alnFileFormat,new ArrayList<String>());
			
			String pwmFile = argMap.getMandatory("pwm");
			PositionWeightMatrixIO pwmIO = new PositionWeightMatrixIO();
			FileInputStream fis = new FileInputStream(pwmFile);
			pwmIO.load(fis);
			fis.close();
			List<PositionWeightMatrix> pwms = pwmIO.getMatrices();
			
			Iterator<PositionWeightMatrix> pwmIt = pwms.iterator();
			int maxPWMLength = 0;
			HashMap<PositionWeightMatrix, BufferedWriter> writerMap = new HashMap<PositionWeightMatrix, BufferedWriter>();
			List<PositionWeightMatrix> rpwms = new ArrayList<PositionWeightMatrix>(pwms.size());
			List<PositionWeightMatrixModel> pwmms = new ArrayList<PositionWeightMatrixModel>(pwms.size());
			List<PositionWeightMatrixModel> rpwmms = new ArrayList<PositionWeightMatrixModel>(pwms.size());
			while(pwmIt.hasNext()) {
				PositionWeightMatrix pwm = pwmIt.next();
				maxPWMLength = pwm.size() > maxPWMLength ? pwm.size() : maxPWMLength;
				writerMap.put(pwm, new BufferedWriter(new FileWriter(argMap.getOutputDir() + "/" + pwm.getName() + "_samples.bed")));
				PositionWeightMatrix rpwm = pwm.reverseComplement();
				rpwms.add(rpwm);
				PositionWeightMatrixModel pwmm = new PositionWeightMatrixModel(pwm, scaler.getModel());
				PositionWeightMatrixModel rpwmm = new PositionWeightMatrixModel(rpwm, scaler.getModel());
				pwmms.add(pwmm);
				rpwmms.add(rpwmm);
			}
			
			PositionWeightMatrix bg = pwms.get(0).createIsoPWM(scaler.getModel().getParameters().getBackgroundNucleotideFreqs(), "bg");
			
			MAFIO mafio = new MAFIO(alnFile, true);
			
			Random r = new Random();
			//BufferedWriter outBW = argMap.getOutputWriter();
			int samplesSoFar = 0;
			while(samplesSoFar < numOfSamples) {
				int startSample = start + r.nextInt(end - start - maxPWMLength );
				MultipleAlignment alnWindow = ConservationUtils.setUpMAF(mafio, new ArrayList<String>(), scaler.getModel(), startSample, startSample + maxPWMLength);
				if(alnWindow.getReference() == null || alnWindow.length() < maxPWMLength) {continue;}
				//alnWindow.setIOHelper(MultipleAlignmentIOFactory.create("PHYLIP"));

				 //alnWindow.write(outBW);
				 //outBW.flush();
				char [] windowRefChrs = alnWindow.getReference().getSequenceBases().toCharArray();
				//alnWindow.encodeAsMatrix();
				BED sample = new BED(alnWindow.getReference());
				sample.setChromosome(alnWindow.getReference().getChromosome());
				for(int j = 0; j < pwms.size(); j++) {
					PositionWeightMatrix pwm = pwms.get(j);
					PositionWeightMatrix rpwm = rpwms.get(j);
					PositionWeightMatrixModel pwmm = pwmms.get(j);
					PositionWeightMatrixModel rpwmm = rpwmms.get(j);
					MultipleAlignment subAlignment = alnWindow.getSubAlignment(startSample, startSample + pwm.size(), false);
					subAlignment.encodeAsMatrix();
					double bgScore = bg.getLogLikelihood(windowRefChrs);
					double directScore = pwm.getLogLikelihood(windowRefChrs) - bgScore;
					double reverseScore = rpwm.getLogLikelihood(windowRefChrs) - bgScore;
					//System.out.println("\tGot scores: " + System.currentTimeMillis());
					boolean directMatch = directScore > reverseScore;
					double score = directMatch ? directScore : reverseScore;

					Map<String, Matrix> windowAlignment = subAlignment.getColumnsAsVector(alnWindow.getReferenceStart(), pwm.size());
					//windowAlignment.get(subAlignment.getReferenceId()).print(7, 7);
					double conservedScore = directMatch ? pwmm.score(windowAlignment) : rpwmm.score(windowAlignment);
					sample.setScore(conservedScore);
					sample.setOrientation(directMatch);
					BufferedWriter bw = writerMap.get(pwm);
					bw.write(sample.toString(false)  + "\t" + score);
					bw.newLine();
				}
				samplesSoFar++;
			}
			mafio.destroyFileHandle();
			//outBW.close();
			pwmIt = pwms.iterator();
			while(pwmIt.hasNext()) {
				PositionWeightMatrix pwm = pwmIt.next();
				BufferedWriter bw = writerMap.get(pwm);
				if(bw != null) {bw.close();}
			}

			
		} else if ("16".equals(argMap.getTask())) {
			File modelFile = new File(argMap.getMandatory("mod"));
			scaler.setNeutralModel(modelFile);
			String alnDir = argMap.getInputDir();
			String outdir = argMap.getOutputDir();
			float seedMinScore = argMap.containsKey("seedMinScore") ? argMap.getFloat("seedMinScore") : 0;
			float seedQuantile = argMap.containsKey("seedQuantile") ? argMap.getFloat("seedQuantile") : 0f;
			boolean filterScores = argMap.containsKey("minScoreToReport");
			double minScoreToReport = filterScores ? minScoreToReport = argMap.getDouble("minScoreToReport") : -1000;
			
			boolean useQuantileForSeeding = argMap.containsKey("seedQuantile");
			HashMap<String, Double> pwmCutoffs = new HashMap<String, Double>();
			String prefix = argMap.containsKey("outprefix") ? argMap.get("outprefix") : "";
			String mafSuffix  = argMap.containsKey("mafSuffix") ? argMap.getMandatory("mafSuffix") : ".maf";
			String ignoreListStr = argMap.get("ignore");
			List<String> ignoreList = processIgnoreListString(ignoreListStr);
			scaler.model.pruneTree(ignoreList);
			
			int shuffles = argMap.containsKey("shuffles") ? argMap.getInteger("shuffles") : 5;
			
			String pwmFile = argMap.getMandatory("pwm");
			PositionWeightMatrixIO pwmIO = new PositionWeightMatrixIO();
			FileInputStream fis = new FileInputStream(pwmFile);
			pwmIO.load(fis);
			fis.close();
			pwmIO.addPseudoCounts();
			
			List<PositionWeightMatrix> pwms = pwmIO.getMatrices();
			
			Iterator<PositionWeightMatrix> pwmIt = pwms.iterator();
			int maxPWMLength = 0;
			while(pwmIt.hasNext()) {
				PositionWeightMatrix pwm = pwmIt.next();
				maxPWMLength = pwm.size() > maxPWMLength ? pwm.size() : maxPWMLength;
				if(useQuantileForSeeding) {
					List<Double> scoreDist = pwm.computeScoreDistribution((float)scaler.model.getPi().get(0,0),(float) scaler.model.getPi().get(1,1), (float)scaler.model.getPi().get(2,2),(float) scaler.model.getPi().get(3,3), 500000);
					Collections.sort(scoreDist);
					double cutoff = Statistics.quantile(scoreDist, seedQuantile);
					System.err.println("Using cutoff " + cutoff + " for " + pwm.getName() + " min score " + scoreDist.get(0) + " max " + scoreDist.get(scoreDist.size() - 1));
					pwmCutoffs.put(pwm.getName(), cutoff);
				} else {
					pwmCutoffs.put(pwm.getName(), (double)seedMinScore);
				}
			}

			
			Map<String, List<? extends GenomicAnnotation>> regionChrMap = getRegionMapFromParameters(argMap);
			
			Iterator<String> chrIt = regionChrMap.keySet().iterator();
			while(chrIt.hasNext()) {
				String chr = chrIt.next();
				Iterator<? extends GenomicAnnotation> annotIt = regionChrMap.get(chr).iterator();
				String mafAlnName = "chr" + chr + mafSuffix;
				MAFIO  chrMafIO = new MAFIO(alnDir+"/"+mafAlnName, true);
				while(annotIt.hasNext()) {
					LightweightGenomicAnnotation annot = annotIt.next();
					int chunkStart = annot.getStart();
					while(chunkStart < annot.getEnd()) {
						
						int chunkEnd =  shuffles == 0 ? 
								Math.min(chunkStart + MAF_CHUNK_SIZE + maxPWMLength - 1,annot.getEnd()) : //Overlap so one can report a hit at the end of the chunk.
									annot.getEnd();
						scaler.alignment = ConservationUtils.setUpMAF(chrMafIO,ignoreList, scaler.getModel(), chunkStart, chunkEnd); // could improve so that index file gets loaded only once.
						if(!scaler.alignment.isEmpty()) {
							AlignedSequence reference = scaler.alignment.getReference();
							//System.out.println("Aligned ref " + reference.getSequenceBases());
							List<int []> ungappedChunks = reference.findUngappedSequenceChunks();
	
							short [] encodedReference = Sequence.encodeSequenceIgnoreCase(reference.getSequenceBuilder()); 
							scaler.alignment.encodeAsMatrix();
							pwmIt = pwms.iterator();
							while(pwmIt.hasNext()) {
								PositionWeightMatrix pwm = pwmIt.next();
								long start = System.currentTimeMillis();
								EmpiricalDistribution shuffledScoreDist = new EmpiricalDistribution(500, -50, 20);
								List<Double> maxPermutationVals = new ArrayList<Double>(shuffles);
								EmpiricalDistribution []  shuffledScoreDistArray = new EmpiricalDistribution[shuffles + 1];
								//List<Double> shuffledScores = new ArrayList<Double>(annot.length() * shuffles);
								long permStart = System.currentTimeMillis();
								
								long last = System.currentTimeMillis();
								System.err.println("Permutations took " + ((last - permStart)/1000.0));
								Collections.sort(maxPermutationVals);
								long now = System.currentTimeMillis();								
								System.err.println("Sorting took: " + ((now - last)/1000.0) + " #items: " + maxPermutationVals.size());
								last = now;
								//System.out.println("pwm " + pwm.getName() + " cutoff " + pwmCutoffs.get(pwm.getName()));
								List<BED> hits = scaler.slidePWM(pwm, pwmCutoffs.get(pwm.getName()), encodedReference, ungappedChunks, shuffles);
								//Add to permutated distributions
								shuffledScoreDistArray[0] = new EmpiricalDistribution(500, -50,20);
								for(BED hit : hits) {
									shuffledScoreDistArray[0].add(hit.getScore());
									shuffledScoreDist.add(hit.getScore());
								}
								//BufferedWriter err = new BufferedWriter(new PrintWriter(System.err));
								//shuffledScoreDistArray[0].write(err);
								//err.flush();
								//err.close();
								now = System.currentTimeMillis();
								System.err.println("PWM scan took: " + ((now - last)/1000.0) );
								
								//System.out.println("#Hits: " +hits.size());
								//System.out.println("done sliding, "+pwm.getName()+ " took: " + (System.currentTimeMillis()-start)/1000 + " seconds");
								BufferedWriter bw = new BufferedWriter(new FileWriter(outdir + "/" + prefix + pwm.getName() + ".bed",true));
								BufferedWriter significanceBW = null;
								if(shuffles > 0) {
									significanceBW = new BufferedWriter(new FileWriter(outdir + "/" + prefix + pwm.getName() + ".pvals",true));
									significanceBW.write("Location\tScore\tpvalue\tFWER\tFDR");
									significanceBW.newLine();
								}
								try {
									long startWrite = System.currentTimeMillis();
									for(BED hit : hits) {
										//System.out.println("\tminScoreToReport " + minScoreToReport + " hit " + hit.toString());
										if(!filterScores || hit.getScore() > minScoreToReport) {
											if(shuffles > 0) {
												//hit.setScore(Statistics.pvalue(maxPermutationVals, hit.getScore()));
												double pval = 1 -shuffledScoreDist.getCumulativeProbability(hit.getScore());
												double fwer = Statistics.pvalue(maxPermutationVals, hit.getScore(), false);
												double fdr  = ComputeFDR.FDR(shuffledScoreDistArray[0], shuffledScoreDistArray, hit.getScore());
												significanceBW.write(hit.toUCSC() +"\t" +hit.getScore() + "\t" + (pval) + "\t" + fwer + "\t" + fdr);
												significanceBW.newLine();
												//hit.setScore(shuffledScoreDist.getCummulativeProbability(hit.getScore()));
										
											} 
											bw.write(hit.toString(false));										
											bw.newLine();
										}
									}
									now = System.currentTimeMillis();
									System.err.println("wrote all: " + ((now - startWrite)/1000.0) );
								} finally {
									if(bw != null) {bw.close();}
									if(significanceBW != null) {significanceBW.close();}
								}
							}
						}
						chunkStart = chunkStart + MAF_CHUNK_SIZE;
					}
				}
				chrMafIO.destroyFileHandle();
			}
			
		} else if ("17".equals(argMap.getTask())) {
			File modelFile = new File(argMap.getMandatory("mod"));
			double pA = argMap.getDouble("pA");
			double pC = argMap.getDouble("pC");
			double pG = argMap.getDouble("pG");
			double pT = argMap.getDouble("pT");
			
			if(pA + pC + pG + pT != 1) {
				System.err.println("ERROR: pA + pC + pG + pT = " + (pA + pC + pG + pT) + " it is not a distribution");
				return;
			}
			
			scaler.setNeutralModel(modelFile);
			EvolutionaryModel model = scaler.getModel();
			model.changeBackground(pA, pC, pG, pT);		
			EvolutionaryModelParameters modelParams = model.getModelParameters();
			BufferedWriter bw = argMap.getOutputWriter();
			bw.write(modelParams.toString());
			bw.close();
		} else if ("18".equals(argMap.getTask())) {
			File modelFile = new File(argMap.getMandatory("mod"));
			String ignoreListStr = argMap.get("ignore");
			List<String> ignoreList = processIgnoreListString(ignoreListStr);			
			scaler.setNeutralModel(modelFile);			
			scaler.model.pruneTree(ignoreList);
			String sequence = argMap.getMandatory("ancestralSequence");
			if(argMap.containsKey("bg")) {
				String [] bgs = argMap.get("bg").split(",");
				double pA = Double.parseDouble(bgs[0]);
				double pC = Double.parseDouble(bgs[1]);
				double pG = Double.parseDouble(bgs[2]);
				double pT = Double.parseDouble(bgs[3]);
				if(pA + pC + pG + pT != 1) {
					System.err.println("ERROR: pA + pC + pG + pT = " + (pA + pC + pG + pT) + " it is not a distribution");
					return;
				}
				scaler.model.changeBackground(pA, pC, pG, pT);	
			}
			Phylogeny evolvedSequenceTree = scaler.model.evolve(sequence);
			
			System.out.println(evolvedSequenceTree.toNewHampshire(true));

			
		} else if ("bayesian".equalsIgnoreCase(argMap.getTask())) {
			int samplingNumber = 100;
			boolean printFullDist = argMap.containsKey("printFullDistribution") || argMap.containsKey("printFullDist");
			if(printFullDist) {System.err.println("WARNING: Printing full distribution data, this can take a huge amount of space");}
			File modelFile = new File(argMap.getMandatory("mod"));
			String alnFile = argMap.getInput();
			int window = argMap.containsKey("window") ? argMap.getInteger("window") : 1;
			String alnFileFormat = argMap.containsKey("format") ? argMap.get("format") : "FASTA";
			String out = argMap.getOutput();
			String outdist = out + ".dist";

			String ignoreListStr = argMap.get("ignore");
			double minTreeLength = argMap.containsKey("minTreeLength") ? argMap.getDouble("minTreeLength") : MIN_TREE_LENGTH;
			
			scaler.setNeutralModel(modelFile);
			scaler.setMinimumTreeLength(minTreeLength);	

			List<String> ignoreList = processIgnoreListString(ignoreListStr);
			
			scaler.setUpAlignment(argMap, alnFile, alnFileFormat, ignoreList); 
			scaler.alignment.encodeAsMatrix();

			NormalDistribution priorOmega = new NormalDistribution(1,0.5);
			scaler.estimatePosterior(out, outdist, ignoreList, samplingNumber, priorOmega, window, printFullDist);
		
		
		}else if ("likelihood".equalsIgnoreCase(argMap.getTask())) {
			File modelFile = new File(argMap.getMandatory("mod"));
			String ignoreListStr = argMap.get("ignore");
			String alnFile = argMap.getInput();
			String alnFileFormat = argMap.containsKey("format") ? argMap.get("format") : "FASTA";
			List<String> ignoreList = processIgnoreListString(ignoreListStr);			
			scaler.setNeutralModel(modelFile);			
			scaler.model.pruneTree(ignoreList);
			scaler.minimumTreeLength = 0;
			scaler.setUpAlignment(argMap, alnFile, alnFileFormat, ignoreList);
			scaler.alignment.encodeAsMatrix();
			
			System.out.println(scaler.getLogLikelihood(1, scaler.model.getTree(), 0, 1));

			
		}else if ("maximalPWM".equalsIgnoreCase(argMap.getTask())) {
			File modelFile = new File(argMap.getMandatory("mod"));
			scaler.setNeutralModel(modelFile);
			String alnDir = argMap.getInputDir();
			String annotationFile = argMap.getMandatory("regions");
			float seedMinScore = argMap.containsKey("seedMinScore") ? argMap.getFloat("seedMinScore") : 0;
			
			String mafSuffix  = argMap.containsKey("mafSuffix") ? argMap.getMandatory("mafSuffix") : ".maf";
			String ignoreListStr = argMap.get("ignore");
			List<String> ignoreList = processIgnoreListString(ignoreListStr);
			scaler.model.pruneTree(ignoreList);
			BEDReader reader = new BEDReader(annotationFile);
			
			int shuffles = argMap.containsKey("shuffles") ? argMap.getInteger("shuffles") : 5;
			System.err.println("Filtering windows with affine score less than " + seedMinScore);
			String pwmFile = argMap.getMandatory("pwm");
			
			ConservedPWMScanner cpwms = new ConservedPWMScanner(pwmFile);
			cpwms.setModel(scaler.getModel());

			Iterator<PositionWeightMatrix> pwmIt = cpwms.getPWMiterator();
			BufferedWriter significanceBW = argMap.getOutputWriter();
			significanceBW.write("Region\tRegionUCSC\t");
			while(pwmIt.hasNext()) {
				PositionWeightMatrix pwm = pwmIt.next();
				significanceBW.write(pwm.getName()+"_MaxScore\t"+pwm.getName()+"_FWER\t"+pwm.getName()+"_hitStart\t"+pwm.getName()+"_hitEnd" );
			}
			significanceBW.newLine();
			
			cpwms.setAlignmentChunkSize(MAF_CHUNK_SIZE);
			cpwms.setIgnoreList(ignoreList);
			Iterator<String> chrIt = reader.getChromosomeIterator();

			while(chrIt.hasNext()) {
				String chr = chrIt.next();
				String mafAlnName = chr + mafSuffix;
				MAFIO  chrMafIO = new MAFIO(alnDir+"/"+mafAlnName, true);
				
				Map<PositionWeightMatrix, GenomicAnnotation> bestHitsForRegion = cpwms.scan(chrMafIO, reader.getChromosomeBEDs(chr), shuffles, seedMinScore);
				
				pwmIt = bestHitsForRegion.keySet().iterator();
				while(pwmIt.hasNext()) {
					PositionWeightMatrix pwm = pwmIt.next();
					pwm.write(new BufferedWriter(new PrintWriter(System.out)), NumberFormat.getNumberInstance());
					System.out.flush();
					long start = System.currentTimeMillis();
					GenomicAnnotation hit = bestHitsForRegion.get(pwm);
					
					if(hit != null) {
						significanceBW.write("\t"+hit.getScore()+"\t"+hit.getExtraScore(0)+"\t"+hit.getStart() + "\t" + hit.getEnd());
					} else {
						significanceBW.write("\t"+0+"\t"+1+"\t"+0 + "\t" + 0);
					}

				}
				significanceBW.newLine();
				chrMafIO.destroyFileHandle();
			}
			significanceBW.close();
		} else {
			System.out.println(USAGE);
		}
		
	}
	


	private void estimatePosterior(String outFile, String outDistFile, List<String> ignoreList, int samplingNumber, NormalDistribution priorOmega, int window, boolean printFullDist) throws IOException {
		double minOmega = 0.05;
		double maxOmega = 2.5;
		
		double samplingStep = (maxOmega - minOmega)/(double)samplingNumber;
		
		double [] samplingOmegas = new double[samplingNumber ]; // The first 2 are the putative conserved and neutral omegas.
		for(int i = 0; i < samplingNumber; i++) {
			samplingOmegas[i] = i*samplingStep + minOmega;
		}
		Phylogeny alnTree = ConservationUtils.pruneTree(ignoreList, model.getTree());
		System.out.println(alnTree.toNewHampshire(false));

		double alignmentTreeLength = getTotalDistanceFromNode(alnTree.getRoot());
		System.out.println("base tree total length " + alignmentTreeLength);
		if(alignmentTreeLength < minimumTreeLength) {
			return;
		}
		//System.out.println("reference start " + alignment.getReferenceStart());
		List<int[]> ungappedIslands = alignment.getUngappedReferenceIslands();
		Iterator<int []> ungappedRegionIt = ungappedIslands.iterator();

		System.out.println("Alignment starts at " + alignment.getReferenceStart() + " alignment length " + alignment.length());
		model.RPIDecomposition();
		//System.out.println("R: ");
		model.R.print(7, 5);
		
		
		BufferedWriter bw = new BufferedWriter(new FileWriter(outFile));
		BufferedWriter fullDistBW = null;
		if(printFullDist) {
			fullDistBW = new BufferedWriter(new FileWriter(outDistFile+".full"));
		}
		TreeMap<Double, Double> distribution = new TreeMap<Double, Double>();
		for(int i = 0; i < samplingOmegas.length; i++) {
			distribution.put(samplingOmegas[i], 0.0);
			if(printFullDist) {
				fullDistBW.write(String.valueOf(samplingOmegas[i]));
				if(i < samplingOmegas.length - 1) {fullDistBW.write("\t");}
			}
		}
		
		if(printFullDist){fullDistBW.newLine();}
		
		while(ungappedRegionIt.hasNext()) {
			int [] region = ungappedRegionIt.next();
				//System.out.println("\tgood ungapped island " + region[0]+"-"+region[1]);
			for(int i = region[0]; i < region[1] - window + 1; i++) { //Go through the sites within ungapped ref
				int refPosition = i + alignment.getReferenceStart();
				double [] weightedlogLikelihoods = new double[samplingOmegas.length]; 
				double PX = 0; //The denumerator (probability of data P(X)
				for(int j = 0; j < samplingOmegas.length; j++) {
					double prior = priorOmega.cumulativeProbability(samplingOmegas[j]-samplingStep/2.0, samplingOmegas[j]+samplingStep/2.0);
					if(prior == 0) {
						System.err.println("Ooops prior was 0 omega " + samplingOmegas[j]);
					}
					weightedlogLikelihoods[j] = getLogLikelihood(samplingOmegas[j], alnTree, refPosition, window) + Math.log(prior);
					PX += Math.exp(weightedlogLikelihoods[j]);
				}	
				if(Math.log(PX) == Double.NaN) {System.err.println(":-( PX was NAN at " + refPosition + " PX=" + PX);}
				double logPX = Math.log(PX);
				bw.write(refPosition+"\t" +(weightedlogLikelihoods[(int)((0.25 - minOmega)/samplingStep)] - logPX ) + "\t"+ (weightedlogLikelihoods[(int)((1 - minOmega)/samplingStep)] - logPX));
				bw.newLine();
				if(printFullDist) {fullDistBW.write(String.valueOf(refPosition));}
				for(int j = 0; j < samplingOmegas.length; j++) {
					double omegaPosterior = Math.exp(weightedlogLikelihoods[j] - logPX) ;
					if(printFullDist) {fullDistBW.write("\t"+omegaPosterior);}
					distribution.put(samplingOmegas[j], omegaPosterior + distribution.get(samplingOmegas[j]));
				}
				if(printFullDist) {fullDistBW.newLine();}
			}
		}
		bw.close();
		if(printFullDist) {
			fullDistBW.close();
		}
		
		BufferedWriter distBW = new BufferedWriter(new FileWriter(outDistFile));
		for(double omega : distribution.keySet()) {
			distBW.write(omega+"\t"+distribution.get(omega));
			distBW.newLine();
		}
		distBW.close();
	}
	
	
	
	public double getLogLikelihood(double omega, Phylogeny alnTree,  int startPosition, int window) {
		Map<String, Matrix> alignmentWindow = alignment.getColumnsAsVector(startPosition, window);
		alignmentWindow.get("hg18").print(3,3);
		double meanTreeLength = 0;
		for(int j = 0; j < window; j++) {
			List<String> gappedSeqs = ConservationUtils.getGappedSeqsInWindowMatrix(1, alignmentWindow, j);	
			System.err.println(gappedSeqs);
			ConservationUtils.setUninformativeNodes(alignmentWindow, gappedSeqs, j);
			Phylogeny columnTree = ConservationUtils.pruneTree(gappedSeqs, alnTree);
			meanTreeLength += getTotalDistanceFromNode(columnTree.getRoot());
		}
		meanTreeLength = meanTreeLength/(double)window;
		System.err.println("Mean tree length: " + meanTreeLength);
		EvolutionaryModel changedModel = model.copy();
		changedModel.setOmega(omega);
		
		double logLikelihood = 0;
		if(meanTreeLength  > minimumTreeLength) {
			for(int i = 0; i < window; i++){
				double likelihood = changedModel.computeLikelihood(alignmentWindow, alnTree.getRoot(), i);
				if(likelihood == 0){System.err.println("Ups 0 likelihood ... omega " + omega + " iteration " + i +" start position " + startPosition);}
				logLikelihood += Math.log(likelihood);
			}
		}
		
		return logLikelihood;
	}


	private static Map<String, List<? extends GenomicAnnotation>>  getRegionMapFromParameters(ArgumentMap argMap)
			throws ParseException, IOException {
		return argMap.getRegionMapFromParameters();
	}

	private void setUpAlignment(ArgumentMap argMap, String alnFile, String alnFileFormat, List<String> ignoreList)  throws IOException, ParseException,	FileNotFoundException {
		alignment = ConservationUtils.setUpAlignment(argMap, alnFile, alnFileFormat, ignoreList, model);
		ignoreSequences = ignoreList;
	}
	
	private static List<BED> integrateStationaryDistributions(List<? extends Fit> fits, int windowSize, String chr, EmpiricalDistribution neutralDist, int overlap) {
		Stack<BED> windows = new Stack<BED>();
		List<Fit> windowMembers = new ArrayList<Fit>(windowSize);
		Iterator<? extends Fit> fitIt = fits.iterator();
		int num = 0;
		BED lastWindow = null;
		while(fitIt.hasNext()) {
			Fit fit = fitIt.next();
		
			if(windowMembers.size() > 0 &&  fit.getPosition() - windowMembers.get(windowMembers.size() - 1).getPosition() > 1) {
					windowMembers.clear();
			}
			
			if(windowMembers.size() == windowSize) {
				BED window = new BED("w" + num);
				window.setChromosome(chr);
				window.setStart(windowMembers.get(0).getPosition());
				window.setEnd(windowMembers.get(windowSize - 1).getPosition() );
				double score = 0;
				if(lastWindow == null || (lastWindow.getEnd() - window.getStart() <= overlap)) {
					Iterator<Fit> windowMemberIt = windowMembers.iterator();
					while(windowMemberIt.hasNext()) {
						Fit member = windowMemberIt.next();
						score += member.getLogLikelihoodRatio();
					}
					if(neutralDist != null) {
						score = - Math.log(1 - neutralDist.getCumulativeProbability(score));
					}
					window.setScore(score);

					windows.add(window);
					num++;
					lastWindow = window;
				}

				windowMembers.remove(0);
			}
			
			windowMembers.add(fit);
			
			if(windowMembers.size() > windowSize) {
				throw new RuntimeException("ERROR: bug handling windows, window memebers exceed window size");
			}
		}
		
		return windows;
	}

	private static void integrateStationaryDistributions(List<? extends Fit> fits, int windowSize, String chr, EmpiricalDistribution neutralDist, int overlap, BufferedWriter bw) throws IOException {
		List<Fit> windowMembers = new ArrayList<Fit>(windowSize);
		Iterator<? extends Fit> fitIt = fits.iterator();
		int num = 0;
		BED lastWindow = null;
		while(fitIt.hasNext()) {
			Fit fit = fitIt.next();
		
			if(windowMembers.size() > 0 &&  fit.getPosition() - windowMembers.get(windowMembers.size() - 1).getPosition() > 1) {
					windowMembers.clear();
			}
			
			if(windowMembers.size() == windowSize) {
				BED window = new BED("w" + num);
				window.setChromosome(chr);
				window.setStart(windowMembers.get(0).getPosition());
				window.setEnd(windowMembers.get(windowSize - 1).getPosition() );
				double score = 0;
				if(lastWindow == null || (lastWindow.getEnd() - window.getStart() <= overlap)) {
					Iterator<Fit> windowMemberIt = windowMembers.iterator();
					while(windowMemberIt.hasNext()) {
						Fit member = windowMemberIt.next();
						score += member.getLogLikelihoodRatio();
					}
					if(neutralDist != null) {
						score = - Math.log(1 - neutralDist.getCumulativeProbability(score));
					}
					window.setScore(score);
					
					num++;
					lastWindow = window;
				}

				windowMembers.remove(0);
				bw.write(window.toString());
				bw.newLine();
			}
			windowMembers.add(fit);
			
			if(windowMembers.size() > windowSize) {
				throw new RuntimeException("ERROR: bug handling windows, window memebers exceed window size");
			}
		}
		
	}

	private static List<BED> threadElements(LinkedHashMap<Integer, Double> siteScoreMap, double maxScore, double minScore, int maxgap, String chr, int window) {
		Stack<BED> elements = new Stack<BED>();
		
		Iterator<Integer> posIt = siteScoreMap.keySet().iterator();
		while(posIt.hasNext()) {
			int pos = posIt.next();
			double score = siteScoreMap.get(pos);
			if(score < maxScore && score > minScore ) {
				if(elements.empty() || (pos - elements.peek().getEnd()) > maxgap) {
					//System.out.println("pos " + pos + " element end " + (elements.size() > 0 ? elements.peek().getEnd() : "null") );
					BED newElement = new BED("e" + elements.size());
					newElement.setStart(pos);
					newElement.setEnd(pos + window);
					newElement.setChromosome(chr);
					elements.push(newElement);
				} else {
					elements.peek().setEnd(pos + window);
				}
			}
		}
		return elements;
	}

	private void fitPI(BufferedWriter bw, List<String> ignoreList) throws IOException{
		Phylogeny alnTree = ConservationUtils.pruneTree(ignoreList, model.getTree());
		System.out.println(alnTree.toNewHampshire(false));

		double alignmentTreeLength = getTotalDistanceFromNode(alnTree.getRoot());
		System.out.println("base tree total length " + alignmentTreeLength);
		if(alignmentTreeLength < minimumTreeLength) {
			//System.out.println("To few species aligned, alignment tree is too short " + alignmentTreeLength);
			return;
		}
		//System.out.println("reference start " + alignment.getReferenceStart());
		List<int[]> ungappedIslands = alignment.getUngappedReferenceIslands();
		Iterator<int []> ungappedRegionIt = ungappedIslands.iterator();
		//int priorUngappedRegEnd = 0;
		//int lastSiteWithOmega = 0;
		System.out.println("Alignment starts at " + alignment.getReferenceStart() + " alignment length " + alignment.length());
		model.RPIDecomposition();
		//System.out.println("R: ");
		model.R.print(7, 5);
		while(ungappedRegionIt.hasNext()) {
			int [] region = ungappedRegionIt.next();
				//System.out.println("\tgood ungapped island " + region[0]+"-"+region[1]);
			for(int i = region[0]; i < region[1]; i++) { //Go through the sites within ungapped ref
				int refPosition = i + alignment.getReferenceStart();
				PiFit fit = fitPI(alnTree, refPosition);
				if(fit != null) {
					writePIFitInfo(bw, refPosition, fit);
				}
			}
		}
		
	}

	public PiFit fitPI(Phylogeny alnTree,  int i) {
		PiFit fit = null;
		//System.out.println("\tref start  " + refPosition + " num gaps so far " + refGaps);
		//Map<String, short[]> column = alignment.getColumns(i, window);
		//TODO: Update to encode alignment as Matrix and avoid the getColumnsAsVector call.
		Map<String, Matrix> column = alignment.getColumnsAsVector(i, 1);

		Phylogeny siteTree = removeGappedSequences(1, alnTree, column);

		String [] siteLeaves = siteTree.getAllExternalSeqNames();
		String [] originalLeaves = model.getTree().getAllExternalSeqNames();
		int treeBitVal = 0;
		for(int j = 0; j < originalLeaves.length; j++) {
			String leaf = originalLeaves[j];
			for(int k = 0; k < siteLeaves.length;k++) {
				if(siteLeaves[k].equals(leaf)) {
					treeBitVal += (1 << j);
				}
			}
		}

		//TODO: Cache this values for each tree, computation is not cheap.
		double treeDist = getTotalDistanceFromNode(siteTree.getRoot());
		//System.out.println("\t site " + refPosition + " window tree " + siteTree.toNewHampshire(true) + " total dist " + treeDist);
		if(treeDist  > minimumTreeLength) {

			//double omega =  model.fitParameters(column, siteTree, window);
			//TODO: update fitParameters to avoid 
			//System.out.println("Start of fit\n");
			//try {
			fit = model.fitPI(column, siteTree);
			if(fit == null) {
				System.err.println("Could not fit column " + i);
			} else {
				fit.setTreeLength(treeDist);
				fit.setTreeBit(treeBitVal);
			}
		
			//} catch (UnableToFitException e) {
			//	System.out.println("Could not fit tree at ref position " + refPosition);
			//}
			//System.out.println ("omega = " + omega);
			//
			//System.out.println("\nEnd of Fit\n\n");
			//lastSiteWithOmega = i;
		} else {
			//System.out.println("total site tree branch length is too short " + treeDist + " skipping....");
		}
		
		return fit;
	}
	
	
	public List<BED> slidePWM(PositionWeightMatrix pwm, double minSeedScore, short[] encodedRef, List<int []> ungappedChunks, int numPermutations) {
		List<BED> scoredKmers = new ArrayList<BED>();
		PositionWeightMatrix rpwm = pwm.reverseComplement();
		List<PositionWeightMatrix> permPWMs = new ArrayList<PositionWeightMatrix>(numPermutations+1);
		List<PositionWeightMatrix> reversedPermPWMs = new ArrayList<PositionWeightMatrix>(numPermutations+1);
		PositionWeightMatrix bg = pwm.createIsoPWM(model.getParameters().getBackgroundNucleotideFreqs(), "bg");
		permPWMs.add(bg);
		reversedPermPWMs.add(bg);

		PositionWeightMatrixModel pwmm = new PositionWeightMatrixModel(pwm, model);
		PositionWeightMatrixModel rpwmm = new PositionWeightMatrixModel(rpwm, model);
		AlignedSequence reference = alignment.getReference();
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
		for(int [] chunk : ungappedChunks) {
			for(i = chunk[0]; i< chunk[1] - L; i++) {
			//System.out.println("Start of loop: " + System.currentTimeMillis());
			//System.out.println("\tGot next sliding window: " + System.currentTimeMillis());
			//System.out.println("\twindow st " + w.getStart() + " end " + w.getEnd());
			//System.out.println("\tGot sequence bases: " + System.currentTimeMillis());
				double minDirectScore = Double.POSITIVE_INFINITY;
				double minReverseScore = Double.POSITIVE_INFINITY;
				double directLikelihood = pwm.getLogLikelihood(encodedRef, i);
				double reverseLikelihood = rpwm.getLogLikelihood(encodedRef, i) ;
				for(int j = 0; j < permPWMs.size(); j++) {
					minDirectScore = Math.min(minDirectScore, directLikelihood - permPWMs.get(j).getLogLikelihood(encodedRef, i));
					minReverseScore = Math.min(minReverseScore, reverseLikelihood - reversedPermPWMs.get(j).getLogLikelihood(encodedRef, i));
					
				}
				//System.out.println("\tGot scores: " + System.currentTimeMillis());
				boolean directMatch = minDirectScore > minReverseScore;
				double score = directMatch ? minDirectScore : minReverseScore;
				//System.err.println("\taffine score " + score + " min seed score " + minSeedScore);
				if(score >= minSeedScore) {
					int start = reference.getStart() + i;
					//System.out.println("\t\tgood hit, computing Probabilisty score " + start);
					//long startTime = System.currentTimeMillis();
					Map<String, Matrix> regionAlignment = alignment.getColumnsAsVector(start, L);
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
	
	
	public List<BED> slidePWM(PositionWeightMatrix pwm, double minSeedScore, int numPermutations) {
		List<int []> ungappedChunks = this.alignment.getReference().findUngappedSequenceChunks();
		short [] encodedReference = Sequence.encodeSequenceIgnoreCase(this.alignment.getReference().getSequenceBuilder());
		return slidePWM(pwm, minSeedScore, encodedReference, ungappedChunks, numPermutations);
	}

	/**
	 * Filles gaps (missing data) by using neutral model
	 * @param ignoreSeqList List of sequences in tree to ignore.
	 * @param minTreeLength Minimum tree length to consider.
	 * @param sample If true, it will sample a based from the probability vector infered for the site. If false it will set the 
	 * 				probability vector as is.
	 */
	private void fillInAlignmentMissingData( double minTreeLength, boolean sample) {
		alignment.encodeAsMatrix();
		
		Phylogeny workTree = ConservationUtils.pruneTree(ignoreSequences, model.getTree());
		//debugTree(workTree);
		int alnLength = alignment.length();
		for (int j = 0; j < alnLength; j++) {
			//System.out.println("site " + j);
			
			Map<String, Matrix> column = alignment.getColumnsAsVector(j, 1);
			List<String> gappedLeaves = ConservationUtils.getGappedSeqsInWindowMatrix(1, column, 0);
			if(gappedLeaves.size() == 0) {
				continue;
			}
			
			Phylogeny siteTree = ConservationUtils.pruneTree(gappedLeaves, workTree);//removeGappedSequences(1, workTree, column);
			double siteTreeLength = getTotalDistanceFromNode(siteTree.getRoot() );
			if(siteTreeLength >= minTreeLength) {
				//System.out.println("\tskipped ( treeLength to short): " + siteTreeLength);
				gibbsSample(column, workTree, new ArrayList<String>(gappedLeaves));
				Iterator<String> gapSeqsIt = gappedLeaves.iterator();
				while(gapSeqsIt.hasNext()) {
					String seq = gapSeqsIt.next();
					Matrix original = alignment.getAlignedSequence(seq).getVectorEncodedSequence();
					for(int i = 0; i < original.getRowDimension(); i++) {
						original.set(i, j, column.get(seq).get(i, 0));
					}
				}
			}
			
		}
		//model.fillInMissingData(alignment);
	}
	
	private void gibbsSample(Map<String, Matrix> column, Phylogeny tree, List<String> sequencesToSample ) {
		if(sequencesToSample.size() == 0) {
			return;
		}
		
		//Phylogeny siteTree = pruneTreeUsingTreeCache(sequencesToSample, ignoreSequences, tree);
		ConservationUtils.setUninformativeNodes(column, sequencesToSample);
		//System.out.println(siteTree.toNewHampshire(true));
		//debugTree(siteTree);
		model.clearCaches();
		double treelikelihood = model.pruneAndPeel(column, tree.getRoot(),  0);
		//System.out.println("Sufficient statistics calculation results for site i tree likelihood " + treelikelihood);
		Map<Integer, NodeLikelihoodParameters> nodeLikelihoodData = model.getNodeLikelihoodParameterMap();
		
		// This is to debug node statistics
		/*
		Iterator<Integer> nodeIdIt = nodeLikelihoodData.keySet().iterator();
		while(nodeIdIt.hasNext()) {
			int nodeId = nodeIdIt.next();
			NodeLikelihoodParameters nodeStats = nodeLikelihoodData.get(nodeId);
			Matrix nodeBaseProbs = nodeStats.computeProbabilityOfLetterAtNode();
			
			System.out.println("node " + nodeId );
			System.out.println("\talpha " );nodeStats.alpha.print(14, 12);
			System.out.println("\tbeta " ); nodeStats.beta.print(14,12);
			System.out.println("\tTransition Matrix "); nodeStats.transition.print(14, 12);
			System.out.println("\tBase probs at node "); nodeBaseProbs.print(25, 12);
			
			PhylogenyNode node = siteTree.getNode(nodeId);
			for(int i = 0; i < nodeBaseProbs.getRowDimension(); i++) {
				String base = "";
				switch(i) {
				case 0 : base = "A"; break;
				case 1 : base = "C"; break;
				case 2 : base = "G"; break;
				case 3 : base = "T"; break;
				}
				
				String tag  = "n"+nodeId + "_" + base;; 
				if(node.customTagExists(tag)) {
					node.removeCustomTagValue(tag);
				}
				TagValueUnit baseProbUnitVal = new TagValueUnit(tag , nodeBaseProbs.get(i, 0),"probs");
				node.addCustomTagValue(baseProbUnitVal, true);
			}
		}
		*/
		//System.out.println(siteTree.toPhyloXML(0));
		Random r = new Random();
		
		String seqToSample = sequencesToSample.get(r.nextInt(sequencesToSample.size()));
		
		Matrix inferredBaseProbs = getPosteriorProbabilityForLeaf(tree, nodeLikelihoodData, seqToSample); 
		
		double draw = r.nextDouble();
		
		// Sample for the sequence in question.
		boolean baseSet = false;
		double cummulativeProbability = 0;
		double previousCommulativeProbability = 0;
		Matrix newLeaf = new Matrix(inferredBaseProbs.getRowDimension(), 1);
		column.put(seqToSample, newLeaf);
		for(int i = 0; i < inferredBaseProbs.getRowDimension(); i++) {
			cummulativeProbability += inferredBaseProbs.get(i,0);
			
			if(!baseSet && draw >= previousCommulativeProbability && draw < cummulativeProbability ) {
				newLeaf.set(i,0,1);
				baseSet = true;
			} else {
				newLeaf.set(i,0,0);
			}
			
			previousCommulativeProbability = cummulativeProbability;
		}
		
		if(!baseSet) {
			System.out.println("Hummm base was not set I guess we drew a one, it was " + draw + ", the cummulative prob is " + cummulativeProbability);
			newLeaf.set(3,0,1);
		}
		
		sequencesToSample.remove(seqToSample);
		gibbsSample(column, tree, sequencesToSample);
		

			/*
			System.out.print("site parent base probs at ");
			parentBaseProbs.print(6, 10);

			System.out.print("P" + totalDistToParent + " ");
			Pt.print(4, 4);
			System.out.print("inferred base probs at site ");
			inferredBaseProbs.print(6, 3);
			*/
		//}
	}




	private Matrix getPosteriorProbabilityForLeaf(Phylogeny tree, Map<Integer, NodeLikelihoodParameters> nodeLikelihoodData, String seqToSample) {
		PhylogenyNode leaf = tree.getNode(seqToSample);
		PhylogenyNode parent = leaf.getParent();
		/*
		double totalDistToParent = 0;
		
		NodeLikelihoodParameters parentLikelihoods = null;
		while(parentLikelihoods == null && !parent.isRoot()) { //Life is hard: If the work tree has as root a child of the workTree root, we'll get no parent in this way.
			totalDistToParent += parent.getDistanceToParent();
			parent = parent.getParent();
			
			parentLikelihoods = nodeLikelihoodData.get(parent.getID());
			
		}
		
		if(parentLikelihoods == null) { // Case describbed in comment above... 
			parent = siteTree.getRoot();
			parentLikelihoods = nodeLikelihoodData.get(parent.getID());
			
			PhylogenyNode realTreeParentNode = tree.getNode(parent.getID());
			while(!realTreeParentNode.isRoot()) {
				totalDistToParent += realTreeParentNode.getDistanceToParent();
				realTreeParentNode = realTreeParentNode.getParent();
			}
		}
		*/
		NodeLikelihoodParameters parentLikelihoods = nodeLikelihoodData.get(parent.getID());
		Matrix parentBaseProbs = parentLikelihoods.computeProbabilityOfLetterAtNode();
		
		Matrix Pt = model.computeTransitions(leaf.getDistanceToParent());
		
		Matrix inferredBaseProbs = Pt.transpose().times(parentBaseProbs);
		return inferredBaseProbs;
	}


	private void generateNeutralAlignment(int colNums, String format, String out) throws IOException {
		MultipleAlignment ma = MultipleAlignmentFactory.create(format);
		for(int i = 0; i < colNums; i++) {
			Map<String, Short> col = model.sample();
			ma.addShortEncodedColumn(col);
		}
		BufferedWriter bw = new BufferedWriter(new FileWriter(out));
		ma.write(bw);
		bw.close();
	}

	private Map<Integer, Map<Double, Double>> scaleTreeWithRandomPrunnings(int window) {
		
		Map<Integer, Map<Double,Double>> scalings = new LinkedHashMap<Integer, Map<Double,Double>>();
		List<String> seqsNotInAlignment = new ArrayList<String>();
		Phylogeny alnTree = pruneToAlignment(model.getTree(), alignment, seqsNotInAlignment);
		double baseTreeLength =  getTotalDistanceFromNode(alnTree.getRoot());
		
		int siteNum = alignment.length();
		int refGaps = 0;
		Random random = new Random();
		for(int i = 0; i < siteNum - window + 1; i++) {
			//System.out.println("Column " + i + " ");
			Map<String, Matrix> column = alignment.getColumnsAsVector(i, window);
			//System.out.println("Referenceseqid " + alignment.getReferenceId());
			
			int referenceBaseColumnTotal = 0;
			Matrix refSeq = column.get(alignment.getReferenceId());
			for(int k = 0; k < refSeq.getRowDimension(); k++) {
				referenceBaseColumnTotal += refSeq.get(0,k);
			}
			if(referenceBaseColumnTotal == 0) {
				System.out.println("Reference gap at site " + i + " skipping...");
				refGaps++;
				continue;
			}
			Phylogeny siteBaseTree = alnTree;
			List<String> extricatedSeqs = new ArrayList<String>();
			Map<Double,Double> siteScalings = new LinkedHashMap<Double,Double>();
			double siteTreeLength = baseTreeLength;
			while(siteTreeLength >= minimumTreeLength) {
				
				OmegaFit fit = scaleColumn (window, column, siteBaseTree, extricatedSeqs);
				if(fit == null) {
					break;
				}
				siteTreeLength = fit.getTreeLength();
				siteScalings.put(fit.getOmega(), fit.getTotalTime());
				//System.out.println("\t site " + i + " window tree total dist " + scaling[0]);	
				String [] leaves = siteBaseTree.getAllExternalSeqNames();
				int leafToPrune = random.nextInt(leaves.length);	
				//System.out.println("Removing leaf " + leaves[leafToPrune]);
				extricatedSeqs.add(leaves[leafToPrune]);
				siteBaseTree = ConservationUtils.pruneTree(leaves[leafToPrune], siteBaseTree);
				column = alignment.getColumnsAsVector(i, window); // since it was modified we need to get it back complete.
				//System.out.println("prunned tree " + siteBaseTree.toNewHampshire(true));
			}
			if(siteScalings.size() > 0) {
				scalings.put(alignment.getReferenceStart() + i - refGaps, siteScalings);
			}
		}
		return scalings;
	}
	
	private OmegaFit scaleColumn(int window, 
			Map<String, Matrix> column,
			Phylogeny prunedAlnTree, 
			List<String> seqsInAlignmentNoInTree) {

		OmegaFit fit = null;
		
		//System.out.println("removed seqs from tree " + seqsInAlignmentNoInTree);
		Iterator<String> seqsInAlnNoInTreeIt = seqsInAlignmentNoInTree.iterator();
		while(seqsInAlnNoInTreeIt.hasNext()) {
			column.remove(seqsInAlnNoInTreeIt.next());
		}
		
		List<String> toPrune = ConservationUtils.getGappedSeqsInWindowMatrix(window, column, 0);
		//System.out.println("Gapped sequences " + toPrune);
		Iterator<String> pruneSeqsIt = toPrune.iterator();
		while(pruneSeqsIt.hasNext()) {
			column.remove(pruneSeqsIt.next());
		}
		
		
		Phylogeny siteTree = null;
		siteTree = ConservationUtils.pruneTree(toPrune,  prunedAlnTree);
		//System.out.println("Site base tree after prooning gapped seqs " + siteTree.toNewHampshire(true));
		double treeDist = getTotalDistanceFromNode(siteTree.getRoot());
		if(treeDist >= minimumTreeLength) {
			fit =  model.fitOmega(column, siteTree, window);
			fit.setTreeLength(treeDist);
		}
		return fit;
	}
	
	public OmegaFit scaleRegion(List<String> ignoreList,   GenomicAnnotation region) {
		Phylogeny alnTree = ConservationUtils.pruneTree(ignoreList, model.getTree());
		double alnTreeLength = getTotalDistanceFromNode(alnTree.getRoot());
		Map<String, Matrix> encodedAlignment = null;
		try {
		   encodedAlignment  = alignment.getColumnsAsVector(region.getStart(), region.length());
		} catch (ArrayIndexOutOfBoundsException aiobe) {
		    logger.error("Region "+ region.toString() + " is not in array");
		    return null;
		}
		double minTreeLength = alnTreeLength;
		/*
		try {
			MultipleAlignment aln = alignment.getSubAlignment(region.getStart(), region.getEnd() - 1, false);
			aln.setIOHelper(MultipleAlignmentIOFactory.create("PHYLIP"));
			StringWriter sb = new StringWriter();
			BufferedWriter bw = new BufferedWriter(sb);
			aln.write(bw);
			bw.flush();
			System.err.println("Alignment: ");
			System.err.println(sb.toString());
			aln.encodeAsMatrix();
			//encodedAlignment = aln.getColumnsAsVector(0, region.getLength() - 1);
		}catch (Exception e) {
			throw new RuntimeException(e);
		}
		*/
		for(int j = 0; j < region.length(); j++) {
			List<String> gappedSeqs = ConservationUtils.getGappedSeqsInWindowMatrix(1, encodedAlignment, j);
			//System.out.println("gapped seqs " + gappedSeqs);
			ConservationUtils.setUninformativeNodes(encodedAlignment, gappedSeqs, j);
			Phylogeny columnTree = ConservationUtils.pruneTree(gappedSeqs, alnTree);
			minTreeLength = Math.min(alnTreeLength, getTotalDistanceFromNode(columnTree.getRoot()));
		}
		OmegaFit fit =  model.fitOmega(encodedAlignment, alnTree, region.length()  );
		//System.err.print("Aln tree dist " + alnTreeLength);
		//System.err.println(" ... min tree dist " + minTreeLength );
		fit.setTreeLength(minTreeLength);
		return fit ;
	}
	
	
	public OmegaFit scaleRegion(Annotation region) {
		Phylogeny alnTree = model.getTree();
		double alnTreeLength = getTotalDistanceFromNode(alnTree.getRoot());
		Map<String, Matrix> encodedAlignment = null;
		try {
		   encodedAlignment  = alignment.getColumnsAsVector(region.getStart(), region.length());
		} catch (ArrayIndexOutOfBoundsException aiobe) {
			logger.error("Region "+ region.toString() + " is not in array");
		    return null;
		}
		double minTreeLength = alnTreeLength;
		
		for(int j = 0; j < region.length(); j++) {
			List<String> gappedSeqs = ConservationUtils.getGappedSeqsInWindowMatrix(1, encodedAlignment, j);
			//System.out.println("gapped seqs " + gappedSeqs);
			ConservationUtils.setUninformativeNodes(encodedAlignment, gappedSeqs, j);
			Phylogeny columnTree = ConservationUtils.pruneTree(gappedSeqs, alnTree);
			minTreeLength = Math.min(alnTreeLength, getTotalDistanceFromNode(columnTree.getRoot()));
		}
		OmegaFit fit =  model.fitOmega(encodedAlignment, alnTree, region.length()  );
		fit.setRegion(region);
		fit.setTreeLength(minTreeLength);
		return fit ;
	}


	public void scaleTree(int window, BufferedWriter bw, List<String> ignoreList, int overlap) throws IOException {
		//System.out.println("Alignment length : " + alignment.getAlignedSequenceIds().size() + 
		//		" Alignment start " + alignment.getReferenceStart() + " species aligned " + alignment.getAlignedSequenceIds() );

		Phylogeny alnTree = ConservationUtils.pruneTree(ignoreList, model.getTree());

		double alignmentTreeLength = getTotalDistanceFromNode(alnTree.getRoot());
		System.out.println("Using window " + window + " base tree total length " + alignmentTreeLength + " alignment length " + alignment.length());
		window = Math.min(window, alignment.length());
		System.out.println("TREE: " + alnTree.toNewHampshire(true));
		if(alignmentTreeLength < minimumTreeLength) {
			System.out.println("To few species aligned, alignment tree is too short " + alignmentTreeLength);
			return;
		}
		
		List<int[]> ungappedIslands = alignment.getUngappedReferenceIslands();
		Iterator<int []> ungappedRegionIt = ungappedIslands.iterator();
		//System.out.println("Alignment starts at " + alignment.getReferenceStart());
		while(ungappedRegionIt.hasNext()) {
			int [] region = ungappedRegionIt.next();
			if(region[1] - region[0] < window) {
				logger.info("\tjikes ungapped island  is small " + region[0] +"-"+region[1]);
			} else  { //If ungapped region is too small, just forget it.
				//System.out.println("\tgood ungapped island " + region[0]+"-"+region[1]);
				for(int i = region[0]; i < region[1] - window + 1; i = i + window - overlap) { //Go through the sites within ungapped ref
				//for(int i = region[0]; i < region[1] - window; i = i + window - overlap) {
					int refPosition = i + alignment.getReferenceStart();
					//System.out.println("\tref start  " + refPosition );
					//System.out.println("\tColumn " + i + " ");
					Map<String, Matrix> column = alignment.getColumnsAsVector(refPosition, window);

					double minTreeLength = alignmentTreeLength;
					for(int j = 0; j < window; j++) {
						List<String> gappedSeqs = ConservationUtils.getGappedSeqsInWindowMatrix(1, column, j);						
						ConservationUtils.setUninformativeNodes(column, gappedSeqs, j);
						Phylogeny columnTree = ConservationUtils.pruneTree(gappedSeqs, alnTree);
						minTreeLength = Math.min(alignmentTreeLength, getTotalDistanceFromNode(columnTree.getRoot()));
					}
					//TODO: the logic above is repeated for all but the new column in the window, we need to change things a bit to avoid redundanty compute the same things 
					//TODO: Cash this values for each tree, computation is not cheap.
					//System.out.println("\t site " + refPosition + " total dist " + minTreeLength);
					if(minTreeLength  > minimumTreeLength) {
						//TODO: update fitParameters to avoid 
						OmegaFit fit =  model.fitOmega(column, alnTree, window);
						fit.setTreeLength(minTreeLength);
						//System.out.println ("omega = " + fit.getOmega());
						writeSiteOmegaInfo(bw, refPosition, fit);
						//lastSiteWithOmega = i;
						if(i == region[1] - window && (region[1] - region[0] ) > (window - window/3)) { //if last base before to close to alignment end was omeagable set all remaining sites to this one.
							for(int j = 1; j + i< region[1]; j++) {
								writeSiteOmegaInfo(bw, refPosition + j, fit);
							}
						}
					} else {
						//System.out.println("total site tree branch length is too short " + treeDist + " skipping....");
					}
		
				}
			}
		}
			
	}
	
	public ArrayList<OmegaFit> scaleTree(int window, List<String> ignoreList, int overlap) throws IOException {
		//System.out.println("Alignment length : " + alignment.getAlignedSequenceIds().size() + 
		//		" Alignment start " + alignment.getReferenceStart() + " species aligned " + alignment.getAlignedSequenceIds() );

		Phylogeny alnTree = ConservationUtils.pruneTree(ignoreList, model.getTree());
		ArrayList<OmegaFit> rtrn=new ArrayList<OmegaFit>();
		
		double alignmentTreeLength = getTotalDistanceFromNode(alnTree.getRoot());
		//System.out.println("Using window " + window + " base tree total length " + alignmentTreeLength + " alignment length " + alignment.length());
		window = Math.min(window, alignment.length());
		//System.out.println("TREE: " + alnTree.toNewHampshire(true));
		if(alignmentTreeLength < minimumTreeLength) {
			//System.out.println("To few species aligned, alignment tree is too short " + alignmentTreeLength);
			return new ArrayList();
		}
		
		List<int[]> ungappedIslands = alignment.getUngappedReferenceIslands();
		Iterator<int []> ungappedRegionIt = ungappedIslands.iterator();
		//System.out.println("Alignment starts at " + alignment.getReferenceStart());
		while(ungappedRegionIt.hasNext()) {
			int [] region = ungappedRegionIt.next();
			if(region[1] - region[0] < window) {
				//System.err.println("\tjikes ungapped island  is small " + region[0] +"-"+region[1]);
			} else  { //If ungapped region is too small, just forget it.
				//System.out.println("\tgood ungapped island " + region[0]+"-"+region[1]);
				for(int i = region[0]; i < region[1] - window + 1; i = i + window - overlap) { //Go through the sites within ungapped ref
				//for(int i = region[0]; i < region[1] - window; i = i + window - overlap) {
					int refPosition = i + alignment.getReferenceStart();
					//System.out.println("\tref start  " + refPosition );
					//System.out.println("\tColumn " + i + " ");
					Map<String, Matrix> column = alignment.getColumnsAsVector(refPosition, window);
					double minTreeLength = alignmentTreeLength;
					for(int j = 0; j < window; j++) {
						List<String> gappedSeqs = ConservationUtils.getGappedSeqsInWindowMatrix(1, column, j);						
						ConservationUtils.setUninformativeNodes(column, gappedSeqs, j);
						Phylogeny columnTree = ConservationUtils.pruneTree(gappedSeqs, alnTree);
						minTreeLength = Math.min(alignmentTreeLength, getTotalDistanceFromNode(columnTree.getRoot()));
					}
					//TODO: the logic above is repeated for all but the new column in the window, we need to change things a bit to avoid redundanty compute the same things 
					//TODO: Cash this values for each tree, computation is not cheap.
					//System.out.println("\t site " + refPosition + " total dist " + minTreeLength);
					if(minTreeLength  > minimumTreeLength) {
						//TODO: update fitParameters to avoid 
						OmegaFit fit =  model.fitOmega(column, alnTree, window);
						fit.setTreeLength(minTreeLength);
						LightweightGenomicAnnotation reg=new BasicGenomicAnnotation("",chr, refPosition, refPosition+window);
						fit.setRegion(reg);
						rtrn.add(fit);
						//System.out.println ("omega = " + fit.getOmega());
						//writeSiteOmegaInfo(bw, refPosition, fit);
						//lastSiteWithOmega = i;
						
					} 
		
				}
			}
		}
			return rtrn;
	}
	
	public void scaleTreeWithSampling(int window, BufferedWriter bw, List<String> ignoreList,  int numberToAverage) throws IOException {
		//System.out.println("Alignment length : " + alignment.getAlignedSequenceIds().size() + 
		//		" Alignment start " + alignment.getReferenceStart() + " species aligned " + alignment.getAlignedSequenceIds() );
		System.out.println("Sampling scaling");
		Phylogeny alnTree = ConservationUtils.pruneTree(ignoreList, model.getTree());

		double alignmentTreeLength = getTotalDistanceFromNode(alnTree.getRoot());
		System.out.println("base tree total length " + alignmentTreeLength);
		if(alignmentTreeLength < minimumTreeLength) {
			//System.out.println("To few species aligned, alignment tree is too short " + alignmentTreeLength);
			return;
		}
		
		System.out.print("Getting ungapped reference islands ... ");
		List<int[]> ungappedIslands = alignment.getUngappedReferenceIslands();
		System.out.println(" Got " + ungappedIslands.size());
		Iterator<int []> ungappedRegionIt = ungappedIslands.iterator();
		//int refGaps = 0;
		//int priorUngappedRegEnd = 0;
		//int lastSiteWithOmega = 0;
		//int lastRegionEnd = 0;
		List<String> gappedLeaves = null;
		Map<String, Matrix> column = null;
		//System.out.println("Alignment starts at " + alignment.getReferenceStart());
		while(ungappedRegionIt.hasNext()) {
			int [] region = ungappedRegionIt.next();
			System.out.println("Doing region " + region[0] + "-" + region[1]);
			//refGaps += region[0] - lastRegionEnd;
			
			if(region[1] - region[0] < window) {
				//System.out.println("\tjikes ungapped island  is small " + region[0] +"-"+region[1]);
			}
			if(region[1] - region[0] >= window) { //If ungapped region is too small, just forget it.
				//System.out.println("\tgood ungapped island " + region[0]+"-"+region[1]);
				for(int i = region[0]; i < region[1] - window + 1; i++) { //Go through the sites within ungapped ref
					int refPosition = i + alignment.getReferenceStart();
					//System.out.println("\tref start  " + refPosition + " num gaps so far " + refGaps);
					//System.out.println("Column " + i + " ");
					//Map<String, short[]> column = alignment.getColumns(i, window);
					//TODO: Update to encode alignment as Matrix and avoid the getColumnsAsVector call.
					column = alignment.getColumnsAsVector(i, window);
					gappedLeaves = ConservationUtils.getGappedSeqsInWindowMatrix(window,column, 0);
					//Phylogeny siteTree = removeSequences( alnTree, column, gappedLeaves);
					Phylogeny siteTree = ConservationUtils.pruneTree(gappedLeaves, alnTree);//removeGappedSequences(window, alnTree, column);					
					//TODO: Cash this values for each tree, computation is not cheap.
					double treeDist = getTotalDistanceFromNode(siteTree.getRoot());
					
					//TODO: Cash this values for each tree, computation is not cheap.
					//double treeDist = getTotalDistanceFromNode(siteTree.getRoot());
					//System.out.println("\t site " + refPosition + " window tree " + siteTree.toNewHampshire(true) + " total dist " + treeDist);
					if(treeDist > minimumTreeLength) {
						//double omega =  model.fitParameters(column, siteTree, window);
						//TODO: update fitParameters to avoid
						//System.out.print("Total free memomry before fit " + Runtime.getRuntime().freeMemory());
						double omegaAvg = 0;
						double transitionAvg = 0;
						//System.out.println("Omega iterations ");
						int iterations = gappedLeaves.size() == 0 ? 1 : numberToAverage;
						int actualIterations = 0;
						OmegaFit averagedFit = new OmegaFit();
						for(int k = 0; k < iterations; k++) {
							gibbsSample(column, alnTree, new ArrayList<String>(gappedLeaves));
							OmegaFit fit =  model.fitOmega(column, alnTree, window);
							actualIterations++;
							averagedFit.setOmega(averagedFit.getOmega() + fit.getOmega());
							averagedFit.setFittedLogLikelihood(averagedFit.getFittedLogLikelihood() + fit.getFittedLogLikelihood());
							averagedFit.setInitialLogLikelihood(averagedFit.getInitialLogLikelihood() + fit.getInitialLogLikelihood());

							if(k ==0 && fit.getOmega() > MIN_INTERESTING_OMEGA) {
								break;
							}
							if(k < numberToAverage - 1) {
								column = alignment.getColumnsAsVector(i, window);
							}
							//System.out.print(" " + data[0]);
							//System.out.print(" it " + k + " free mem: " + Runtime.getRuntime().freeMemory());
						}
						averagedFit.setFittedLogLikelihood(averagedFit.getFittedLogLikelihood()/(double)actualIterations);
						averagedFit.setInitialLogLikelihood(averagedFit.getInitialLogLikelihood()/(double)actualIterations);
						averagedFit.setTreeLength(treeDist);
						averagedFit.setOmega(averagedFit.getOmega()/(double)(actualIterations));
						writeSiteOmegaInfo(bw, refPosition, averagedFit);
						//lastSiteWithOmega = i;
						if(i == region[1] - window ) { //if last base before to close to alignment end was omeagable set all remaining sites to this one.
							for(int j = 1; j + i< region[1]; j++) {
								writeSiteOmegaInfo(bw, refPosition + j, averagedFit);
							}
						}
						//System.out.println(" Finished iteration, free memory: " + Runtime.getRuntime().freeMemory());
					} else {
						//System.out.println("total site tree branch length is too short " + treeDist + " skipping....");
					}
		
				}
			}
		}
			
	}
	
	public void setAlignment(MultipleAlignment alignment) {
		this.alignment = alignment;
	}

	private Phylogeny removeGappedSequences(int window, Phylogeny alnTree, Map<String, Matrix> column) {
		List<String> toPrune = ConservationUtils.getGappedSeqsInWindowMatrix(window, column, 0);
		Phylogeny siteTree = removeSequences(alnTree, column, toPrune);
		return siteTree;
	}

	private Phylogeny removeSequences(Phylogeny alnTree, Map<String, Matrix> column, List<String> toPrune) {
		Iterator<String> pruneSeqsIt = toPrune.iterator();
		while(pruneSeqsIt.hasNext()) {
			column.remove(pruneSeqsIt.next());
		}
		Phylogeny siteTree = pruneTreeUsingTreeCache(toPrune, null, alnTree);
		return siteTree;
	}

	private void writeSiteOmegaInfo(BufferedWriter bw, int refPosition, OmegaFit fit) throws IOException {
		bw.write(String.valueOf((refPosition)));//bw.write(String.valueOf((alignment.getReferenceStart() + i - refGaps)));
		bw.write("\t");
		//bw.write(numberFormat.format(fit.getOmega()));
		bw.write(numberFormat.format(fit.getOmega()));
		bw.write("\t");
		bw.write(numberFormat.format(fit.getTreeLength()));
		//bw.write("\t");
		//bw.write(String.valueOf(transitions));
		/*
		if(model.binStatistics != null && model.binStatistics.size() > 0) {
			bw.write("\t");
			bw.write(String.valueOf(model.getPVal(fit, useLogOddsLikelihoodForPval)));
		}
		*/
		bw.write("\t");
		bw.write(numberFormat.format( fit.getLogOddsScore()));
		bw.write("\t");
		bw.write(tinnyNumberFormat.format(fit.getPVal()));
		bw.newLine();
	}
	
	private void writePIFitInfo(BufferedWriter bw, int refPosition, PiFit fit) throws IOException{
		bw.write(String.valueOf((refPosition)));//bw.write(String.valueOf((alignment.getReferenceStart() + i - refGaps)));
		if(fit != null) {
			for(int i = 0; i < fit.getPI().getRowDimension(); i++) {
				bw.write("\t");					
				bw.write(numberFormat.format(fit.getPI().get(i, 0)));
			}
			bw.write("\t");
			//double logLikelihoodRatio = fit.getLogLikelihoodRatio();
			//bw.write(numberFormat.format(logLikelihoodRatio < 0 ? 0 : logLikelihoodRatio));
			bw.write(numberFormat.format(fit.getLogLikelihoodRatio()));
			bw.write("\t");
			bw.write(numberFormat.format(fit.getTreeLength()));
			bw.write("\t");
			bw.write(String.valueOf(fit.getTreeBit()));
		} else {
			bw.write("\t0.25\t0.25\t0.25\t0.25\t0\t0");
		}
		bw.newLine();
	}
	
	private Phylogeny pruneToAlignment(Phylogeny tree, MultipleAlignment ma, List<String> speciesNotInAlignment) {
		String [] leaves = tree.getAllExternalSeqNames();
		for(int l = 0; l < leaves.length; l++) {
			if(alignment.getAlignedSequence(leaves[l]) == null) {
				speciesNotInAlignment.add(leaves[l]);
			}
		}
		//System.out.println("Species not in alignment " + speciesNotInAlignment);
		Phylogeny alnTree = pruneTreeUsingTreeCache(speciesNotInAlignment, null, model.getTree());
		return alnTree;
	}
	
	private Phylogeny pruneTreeUsingTreeCache(List<String> toPrune, List<String> otherCacheLookupLeaves, Phylogeny treeToPrune) {
		Phylogeny p = treeToPrune;
		if(toPrune.size() > 0) {
			Collections.sort(toPrune);
			StringBuilder key = new StringBuilder();
			List<String> allLeaves = toPrune;
			if(otherCacheLookupLeaves != null && otherCacheLookupLeaves.size() > 0) {
				allLeaves = new ArrayList<String>(toPrune.size() + otherCacheLookupLeaves.size());
				allLeaves.addAll(otherCacheLookupLeaves);
				allLeaves.addAll(toPrune);
			}
			Iterator<String> it = allLeaves.iterator();
			while(it.hasNext()) {
				String seqId = it.next();
				key.append(seqId);
			}
			//System.out.print("  key to pruned tree table " + key);
			p = prunnedTrees.get(key.toString());
			if(p == null) {
				p = ConservationUtils.pruneTree(toPrune, treeToPrune);
				//prunnedTrees.put(key.toString(), p);
			}
		}
		return p;
	}
	
	public double getTotalDistanceFromNode(PhylogenyNode n) {
		double dist = 0d;
		if(n != null && !n.isExternal()) {
			dist += n.getChildNode1().getDistanceToParent() + getTotalDistanceFromNode(n.getChildNode1());
			dist += n.getChildNode2().getDistanceToParent() + getTotalDistanceFromNode(n.getChildNode2());
			
		}
		
		return dist;
	}
	
	
	private double computeTotalDistance(Phylogeny t, List<String> seqs) {
		Iterator<String> seqIt = seqs.iterator();
		double dist = 0;
		while(seqIt.hasNext()) {
			String seq = seqIt.next();
			dist += t.getNode(seq).getDistanceToParent();
		}
		
		return dist;
	}

	public void setNeutralModel(File modelFile) throws IOException, ParseException {
		EvolutionaryModelParameters modelParams = new EvolutionaryModelParameters(modelFile);
		model = new EvolutionaryModel(modelParams);
	}

	public double getMinimumTreeLength() {
		return minimumTreeLength;
	}

	public void setMinimumTreeLength(double minimumTreeLength) {
		this.minimumTreeLength = minimumTreeLength;
	}
	
	private static List<String> processIgnoreListString(String ignoreListStr) {
		return ConservationUtils.commaSeparatedStringToList(ignoreListStr);
	}
	
	private void debugTree(Phylogeny tree) {
		debugSubtree(tree.getRoot(), 0);
	}
	
	private void debugSubtree(PhylogenyNode node, int tabs) {
		for(int i = 0; i < tabs; i++) {
			System.out.print("\t");
		}
		System.out.println(node.getID() + "[" + node.getSeqName() + "]");
		if(!node.isExternal()) {
			debugSubtree(node.getChildNode1(), tabs + 1);
			debugSubtree(node.getChildNode2(), tabs + 1);
		}
	}

	public EvolutionaryModel getModel() {
		return model;
	}

	public void setModel(EvolutionaryModel model) {
		this.model = model;
	}
}
