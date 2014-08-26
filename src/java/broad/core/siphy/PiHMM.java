package broad.core.siphy;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Date;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;
import java.util.Random;

import org.forester.phylogeny.Phylogeny;

import broad.core.hmm.BadModelException;
import broad.core.hmm.MarkovModel;
import broad.core.hmm.MarkovState;
import broad.core.multiplealignment.MultipleAlignment;
import broad.core.multiplealignment.MultipleAlignmentFactory;

import Jama.Matrix;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;

public class PiHMM extends MarkovModel<Map<String, Matrix>>{
	
	private EvolutionaryModel neutralModel;
	private double nu;
	private double mu;
	
	private static final double SMALL_PROB = 0.03;
	private static final double SMALL_NUM = 0.0000001;
	
	private static final double[][] conserved_pis = { 
		{1 - 3*SMALL_PROB, SMALL_PROB, SMALL_PROB, SMALL_PROB},
		{SMALL_PROB, 1 - 3*SMALL_PROB, SMALL_PROB, SMALL_PROB},
		{SMALL_PROB, SMALL_PROB, 1 - 3*SMALL_PROB, SMALL_PROB},
		{SMALL_PROB, SMALL_PROB, SMALL_PROB, 1 - 3*SMALL_PROB},
		
		{0.5 - SMALL_PROB, 0.5 - SMALL_PROB, SMALL_PROB, SMALL_PROB},
		{0.5 - SMALL_PROB, SMALL_PROB, 0.5 - SMALL_PROB, SMALL_PROB},
		{0.5 - SMALL_PROB, SMALL_PROB, SMALL_PROB, 0.5 - SMALL_PROB},
		{SMALL_PROB, 0.5 - SMALL_PROB, 0.5 - SMALL_PROB, SMALL_PROB},
		{SMALL_PROB, 0.5 - SMALL_PROB, SMALL_PROB, 0.5 - SMALL_PROB},
		{SMALL_PROB, SMALL_PROB, 0.5 - SMALL_PROB, 0.5 - SMALL_PROB}
	};
	
	private static final  String [] conserved_pi_names = {"A", "C", "G", "T", "AC", "AG", "AT", "CG", "CT", "GT"};
	
	
	public PiHMM(EvolutionaryModel neutral, int stateNumber) {
		super(stateNumber);
		this.neutralModel = neutral;
		addState(new PiStateModel(neutralModel, "neutral"));
	}
	
	public PiHMM(EvolutionaryModel neutral, double minExpectedElementLength, double expectedCoverage) {
		super(conserved_pis.length + 1);
		this.neutralModel = neutral;
		addState(new PiStateModel(neutralModel, "neutral"));
		setMuNuAndTransitions(minExpectedElementLength, expectedCoverage);
	}
	
	
	public void addConservedPi(Matrix pi, String name) {
		addState(new PiStateModel(neutralModel, pi, name));
	}
	
	public short[] viterbiMostLikelyEstimation(MultipleAlignment alignment) throws BadModelException {
		return super.viterbiMostLikelyEstimation(new AlignmentListWrapper(alignment) );
	}
	
	public double computePathLogLikelihood(MultipleAlignment alignment, short [] path, int initialStateIdx, int endStateIdx) {
		return super.computePathLogLikelihood(new AlignmentListWrapper(alignment), path, initialStateIdx, endStateIdx);
	}
	
	public ForwardResult<Map<String, Matrix>> runForwardAlgorithm(MultipleAlignment alignment) throws BadModelException {
		return super.runForwardAlgorithm(new AlignmentListWrapper(alignment) );
	}
	
	public BackwardResult<Map<String, Matrix>> runBackwardAlgorithm(MultipleAlignment alignment, ForwardResult<Map<String, Matrix>> forwardData) throws BadModelException {
		return super.runBackwardAlgorithm(new AlignmentListWrapper(alignment) , forwardData);
	}
	
	public static class PiStateModel implements MarkovState<Map<String, Matrix>> {
		EvolutionaryModel model;
		String name;
		double emissionProbOfUnalignRegion = 1;

		public PiStateModel(EvolutionaryModel neutralModel, String name) {
			model = neutralModel.copy();
			this.name = name;
		}
		
		public PiStateModel(EvolutionaryModel neutralModel, Matrix pi, String name) {
			model = neutralModel.copy();
			model.changeBackground(pi);
			this.name = name;
		}
		
		public double getEmissionLogProbability(Map<String, Matrix> alignmentColumn) {
			return Math.log(getEmissionProbability(alignmentColumn));
		}

		public double getEmissionProbability(Map<String, Matrix> alignmentColumn) {
			List<String> gappedLeaves = ConservationUtils.getGappedSeqsInWindowMatrix(1, alignmentColumn, 0);
			//if(gappedLeaves.size() <= alignmentColumn.keySet().size() - 1) {
			//	return emissionProbOfUnalignRegion;
			//}
			Phylogeny siteTree = ConservationUtils.pruneTree(gappedLeaves, model.getTree());
			//System.err.println("gapedLeaves: " + gappedLeaves);
			return model.computeLikelihood(alignmentColumn, siteTree.getRoot(), 0);
		}

		public String getName() {
			return name;
		}
		
		public void setEmissionProbOfUnalignRegion(double prob) {
			this.emissionProbOfUnalignRegion = prob;
		}

		public Map<String, Matrix> emitObservation() {
			Map<String, Short> column = model.sample();
			Map<String, Matrix> columnMatrix = new LinkedHashMap<String, Matrix>();
			for (String key : column.keySet()) {
				short val = column.get(key);
				Matrix valMat = new Matrix(model.alphabetSize, 1);
				valMat.set(val, 0, 1);
				columnMatrix.put(key, valMat);
			}
			return columnMatrix;
		}
		
		
		
	}
	
	public static class ConstrainedModel implements MarkovState<Map<String, Matrix>> {
		List<PiStateModel> conservedStates;
		String name;

		public ConstrainedModel( EvolutionaryModel neutralModel, String name) {
			this.name = name;
			conservedStates = new ArrayList<PiStateModel>();
			for(int i = 0; i < conserved_pis.length; i++) {
				Matrix pi = new Matrix(conserved_pis[i].length,1);
				for(int j = 0; j < conserved_pis[i].length; j++) {
					pi.set(j, 0, conserved_pis[i][j]);
				}
				PiStateModel piModel = new PiStateModel(neutralModel, pi, conserved_pi_names[i]);
				piModel.setEmissionProbOfUnalignRegion(0);
				conservedStates.add(piModel);;
			}
		}
		
		public double getEmissionLogProbability(Map<String, Matrix> alignmentColumn) {
			return Math.log(getEmissionProbability(alignmentColumn));
		}

		public double getEmissionProbability(Map<String, Matrix> alignmentColumn) {
			Iterator<PiStateModel> stateIt = conservedStates.iterator();
			double emissionProb = 0;
			while(stateIt.hasNext()) {
				PiStateModel state = stateIt.next();
				emissionProb = emissionProb + state.getEmissionProbability(alignmentColumn);
			}
			return emissionProb/(double)conservedStates.size();
		}

		public String getName() {
			return name;
		}
		
		public Map<String, Matrix> emitObservation() {
			//First pick a state, then emit from it
			Random r = new Random();
			int sampledStateIdx = r.nextInt(conservedStates.size());
			return conservedStates.get(sampledStateIdx).emitObservation();
		}
	}
	
	private static class AlignmentListWrapper implements List<Map<String, Matrix>> {
		MultipleAlignment alignment;
		
		public AlignmentListWrapper(MultipleAlignment alignment) {
			this.alignment = alignment;
		}
		
		public Map<String, Matrix> get(int columnNumber) {
			return alignment.getColumnsAsVector(columnNumber + alignment.getReferenceStart(), 1);
		}


		public boolean isEmpty() {
			return alignment.isEmpty();
		}
		
		public int size() {
			return alignment.length();
		}

		
		// Only implemented methods above, have not needed what follows.

		public boolean add(Map<String, Matrix> col) {
			return false;
		}

		public void add(int arg0, Map<String, Matrix> arg1) {

		}

		public boolean addAll(Collection<? extends Map<String, Matrix>> arg0) {
			return false;
		}

		public boolean addAll(int arg0, Collection<? extends Map<String, Matrix>> arg1) {
			return false;
		}

		public void clear() {
			// TODO Auto-generated method stub
			
		}

		public boolean contains(Object arg0) {
			return false;
		}

		public boolean containsAll(Collection<?> arg0) {
			// TODO Auto-generated method stub
			return false;
		}

		public int indexOf(Object arg0) {
			// TODO Auto-generated method stub
			return 0;
		}
		
		public Iterator<Map<String, Matrix>> iterator() {
			// TODO Auto-generated method stub
			return null;
		}

		public int lastIndexOf(Object arg0) {
			// TODO Auto-generated method stub
			return 0;
		}

		public ListIterator<Map<String, Matrix>> listIterator() {
			// TODO Auto-generated method stub
			return null;
		}

		public ListIterator<Map<String, Matrix>> listIterator(int arg0) {
			// TODO Auto-generated method stub
			return null;
		}

		public boolean remove(Object arg0) {
			// TODO Auto-generated method stub
			return false;
		}

		public Map<String, Matrix> remove(int arg0) {
			// TODO Auto-generated method stub
			return null;
		}

		public boolean removeAll(Collection<?> arg0) {
			// TODO Auto-generated method stub
			return false;
		}

		public boolean retainAll(Collection<?> arg0) {
			// TODO Auto-generated method stub
			return false;
		}

		public Map<String, Matrix> set(int arg0, Map<String, Matrix> arg1) {
			// TODO Auto-generated method stub
			return null;
		}

		public List<Map<String, Matrix>> subList(int arg0, int arg1) {
			// TODO Auto-generated method stub
			return null;
		}

		public Object[] toArray() {
			// TODO Auto-generated method stub
			return null;
		}

		public <T> T[] toArray(T[] arg0) {
			// TODO Auto-generated method stub
			return null;
		}
		
	}
	
	public static final String USAGE = "Usage: PiHMM TASK=<task_num> <task_args>\n" +
	"\tTasks:\n" +
	"\t\t1. Segment genome into pi-conserved and noncoserved.  -in <multiple alignment file>  -mod <Neutral Evolutionary model consisting of aminoacid background distribution, mutation matrix and neutral phylogenetic tree>" +
	"\n\t\t -l <The smoothness parameter, and should be the expected average conserved element length> and -gamma <The expected coverage parameter, how much of the genome is expected to be conserved>" +
	"\n\t\t -format <Alignment format default is FASTA is default> -ignore <comma separated species to ignore> -ref <reference sequence id, necessary if the alignment is not in MAF format>" +
	"\n\t\t  [-start <If MAF file, you may specify the reference start coordinate> -end <If MAF file, you may specify the reference end coordinate>]" +
	"\n\t\t2. Compute posterior probabilities of each site of being pi-conserved.  -in <multiple alignment file>  -mod <Neutral Evolutionary model consisting of aminoacid background distribution, mutation matrix and neutral phylogenetic tree>" +
	"\n\t\t -l <The smoothness parameter, and should be the expected minimum conserved element length> and -gamma <The expected coverage parameter, how much of the genome is expected to be conserved>" +
	"\n\t\t -format <Alignment format default is FASTA is default> -ignore <comma separated species to ignore> -ref <reference sequence id, necessary if the alignment is not in MAF format>" +
	"\n\t\t  [-start <If MAF file, you may specify the reference start coordinate> -end <If MAF file, you may specify the reference end coordinate>]" +
	"\n\t\t3. Compute log odds score for annotations in file. Scores reflect the log ratio of the probability of the path through each element being fully conserved or non conserved.  -in <Annotation file default format is assumed to be BED>  -mod <Neutral Evolutionary model consisting of aminoacid background distribution, mutation matrix and neutral phylogenetic tree>" +
	"\n\t\t -l <The smoothness parameter, and should be the expected minimum conserved element length> and -gamma <The expected coverage parameter, how much of the genome is expected to be conserved>" +
	"\n\t\t -format <Alignment format default is FASTA is default> -ignore <comma separated species to ignore> -ref <reference sequence id, necessary if the alignment is not in MAF format>" +
	"\n\t\t  [-start <If MAF file, you may specify the reference start coordinate> -end <If MAF file, you may specify the reference end coordinate>]" +
	"\n\t\t4. Simulate an alignment with for given HMM parameters. \n\t\tParameters:\n\t\t-mod <Base m,odel to use> \n\t\t-ignore <optional -- comma separated species to ignore in the given model tree> \n\t\t-colNum <Number of columns to sample> -\n\t\t-otherPI <Conserved or alternative PI distribution, as a comma separated list of the A,C,G,T frequencies, add as many -otherPI as different PI states are desired. > " +
	"\n\t\t -l <The smoothness parameter, how long in average should be the stretches of the generated alignment should be from the alternative PI model> and -gamma <The coverage parameter, how much of the sampled alignmnet should be made of the alternative PI model> \n\t\t-out <Generated alignment output file>" +
	"\n"; 
	
	private static final int CHUNK_SIZE = 200000;
	private static final int CHUNK_OVERLAP = 500;
	
	public static void main(String[] args) throws Exception {
		
		ArgumentMap argMap = CLUtil.getParameters(args, USAGE);
		
		if ("1".equals(argMap.getTask())) {	
			File modelFile = new File(argMap.getMandatory("mod"));
			String alnFile = argMap.getInput();
			String alnFileFormat = argMap.containsKey("format") ? argMap.get("format") : "FASTA";
			//String ref = argMap.getMandatory("ref");
			String ignoreListStr = argMap.get("ignore");
			double l = argMap.getDouble("l");
			double gamma = argMap.getDouble("gamma");
			
			
			EvolutionaryModelParameters modelParams = new EvolutionaryModelParameters(modelFile);
			List<String> ignoreList = ConservationUtils.commaSeparatedStringToList(ignoreListStr);
			if(!ignoreList.isEmpty()) {
				Phylogeny tree = modelParams.getTree();
				modelParams.setTree(ConservationUtils.pruneTree(ignoreList, tree));
			}
			EvolutionaryModel model = new EvolutionaryModel(modelParams);
			PiHMM hmm = createDefaultTwoStateChain(l, gamma, model);
			
			MultipleAlignment alignment = ConservationUtils.setUpAlignment(argMap, alnFile, alnFileFormat, ignoreList, model);
			//System.out.println("Alignment Size: " + alignment.length() + " human seq: " + alignment.getReference().getSequenceBases());
			
			List<int[]> ungappedIslands = alignment.getUngappedSequenceReferenceIslands();
			Iterator<int []> ungappedRegionIt = ungappedIslands.iterator();
			
			BufferedWriter bw = argMap. getOutputWriter();
			while(ungappedRegionIt.hasNext()) {
				int[] startEnd = ungappedRegionIt.next();
				int islandStart = startEnd[0] + alignment.getReferenceStart();
				int islandEnd   = startEnd[1] + alignment.getReferenceStart();
	 			int chunkStart = islandStart;
				//short [] path = new short [alignment.length()];
				while(chunkStart < islandEnd) {
					int chunkEnd = Math.min(chunkStart + CHUNK_SIZE, islandEnd);
					System.err.println("Chunk Start: " + chunkStart + " end " + chunkEnd);
					MultipleAlignment chunk = alignment.getSubAlignment(Math.max(islandStart, chunkStart - CHUNK_OVERLAP), chunkEnd, false);
					chunk.encodeAsMatrix();
					//System.out.println("Chunk Size: " + chunk.length() + " human seq: " + chunk.getReference().getSequenceBases());
					short [] chunkPath = hmm.viterbiMostLikelyEstimation(chunk);
					int shift = chunkStart == islandStart ? 0 : CHUNK_OVERLAP ;
					for(int i = shift; i < chunkPath.length; i++) {
						//path[chunkStart - alignment.getReferenceStart() + i] = chunkPath[i];
						bw.write((chunkStart + i - shift) + "\t" + chunkPath[i] );
						bw.newLine();
					}
					//System.out.println("chunk start " + chunkEnd + " alignment end " + alignment.getReferenceEnd());
					chunkStart = chunkEnd;
					
				}
			}
			bw.close();
		} else if ("2".equals(argMap.getTask())) {	
			File modelFile = new File(argMap.getMandatory("mod"));
			String alnFile = argMap.getInput();
			String alnFileFormat = argMap.containsKey("format") ? argMap.get("format") : "FASTA";
			//String ref = argMap.getMandatory("ref");
			String ignoreListStr = argMap.get("ignore");
			double l = argMap.getDouble("l");
			double gamma = argMap.getDouble("gamma");
			List<String> ignoreList = ConservationUtils.commaSeparatedStringToList(ignoreListStr);
			
			long initialTime = (new Date()).getTime();
			System.err.println("Start: ");
			EvolutionaryModelParameters modelParams = new EvolutionaryModelParameters(modelFile);
			if(!ignoreList.isEmpty()) {
				Phylogeny tree = modelParams.getTree();
				modelParams.setTree(ConservationUtils.pruneTree(ignoreList, tree));
			}
			EvolutionaryModel model = new EvolutionaryModel(modelParams);


			//PiHMM hmm = new PiHMM(model,l, gamma);
			PiHMM hmm = createDefaultTwoStateChain(l, gamma, model);
			
			MultipleAlignment alignment = ConservationUtils.setUpAlignment(argMap, alnFile, alnFileFormat, ignoreList, model);
			//alignment.encodeAsMatrix();
			System.err.println("Alignment loaded: " + ((new Date()).getTime() - initialTime));
			//System.out.println("Alignment Size: " + alignment.length() + " human seq: " + alignment.getReference().getSequenceBases());
			List<int[]> ungappedIslands = alignment.getUngappedSequenceReferenceIslands();
			Iterator<int []> ungappedRegionIt = ungappedIslands.iterator();
			
			BufferedWriter bw = argMap. getOutputWriter();
			while(ungappedRegionIt.hasNext()) {
				while(ungappedRegionIt.hasNext()) {
					int[] startEnd = ungappedRegionIt.next();
					int islandStart = startEnd[0] + alignment.getReferenceStart();
					int islandEnd   = startEnd[1] + alignment.getReferenceStart();
		 			int chunkStart = islandStart;
		 			while(chunkStart < islandEnd ) {
						int chunkEnd = Math.min(chunkStart + CHUNK_SIZE, islandEnd);
						System.err.println("Chunk Start: " + chunkStart + " end " + chunkEnd);
						MultipleAlignment chunk = alignment.getSubAlignment(Math.max(islandStart, chunkStart - CHUNK_OVERLAP), chunkEnd, false);
						chunk.encodeAsMatrix();
						System.err.println("Chunk encoded: " + ((new Date()).getTime() - initialTime));
						//System.out.println("Chunk Size: " + chunk.length() + " human seq: " + chunk.getReference().getSequenceBases());
						ForwardResult forward = hmm.runForwardAlgorithm(chunk);
						System.err.println("Chunk forward: " + ((new Date()).getTime() - initialTime));
						BackwardResult backward = forward.runBackwardAlgorithm();
						System.err.println("Chunk backward: " + ((new Date()).getTime() - initialTime));
						int shift = chunkStart == islandStart ? 0 : CHUNK_OVERLAP ;
						for(int i = shift; i < chunk.length(); i++) {
							bw.write((chunkStart + i - shift) + "\t" +(1 - backward.getPosteriorProbability(0,  i)));
							bw.newLine();
						}
						System.err.println("Chunk posterior written: " + ((new Date()).getTime() - initialTime));
						chunkStart = chunkEnd ;
		 			}
				}
			}
			bw.close();
		}else if("4".equals(argMap.getTask())) {
			File modelFile = new File(argMap.getMandatory("mod"));
			EvolutionaryModelParameters modelParams = new EvolutionaryModelParameters(modelFile);
			//String ref = argMap.getMandatory("ref");
			String ignoreListStr = argMap.get("ignore");
			double l = argMap.getDouble("l");
			double gamma = argMap.getDouble("gamma");
			String format = argMap.containsKey("format") ? argMap.get("format") : "FASTA";
			int colNum = argMap.getInteger("colNum");
			List<String> ignoreList = ConservationUtils.commaSeparatedStringToList(ignoreListStr);
			if(!ignoreList.isEmpty()) {
				Phylogeny tree = modelParams.getTree();
				modelParams.setTree(ConservationUtils.pruneTree(ignoreList, tree));
			}
			EvolutionaryModel model = new EvolutionaryModel(modelParams);
			//PiHMM hmm = createDefaultTwoStateChain(l, gamma, model);
			String otherPiStr = argMap.getMandatory("otherPi");

			String [] bgs = otherPiStr.split(",");
			double pA = Double.parseDouble(bgs[0]);
			double pC = Double.parseDouble(bgs[1]);
			double pG = Double.parseDouble(bgs[2]);
			double pT = Double.parseDouble(bgs[3]);
			if(Math.abs(pA + pC + pG + pT - 1) > SMALL_NUM ) {
				System.err.println("ERROR: pA + pC + pG + pT = " + (pA + pC + pG + pT) + " it is not a distribution");
				return;
			}
			Matrix piMat = new Matrix(model.alphabetSize, 1);
			piMat.set(0, 0, pA);
			piMat.set(1, 0, pC);
			piMat.set(2, 0, pG);
			piMat.set(3, 0, pT);

			PiHMM hmm = PiHMM.createTwoStateChain(l, gamma, model, piMat);
			List<Integer> emittedPath = hmm.emitPath(colNum);
			MultipleAlignment ma = hmm.generateAlignment(emittedPath, format);
			
			String outFile = argMap.getOutput();
			String pathOutFile = outFile +".path";
		
			BufferedWriter pBw = new BufferedWriter(new FileWriter(pathOutFile));
			pBw.write(emittedPath.toString());
			pBw.close();
		
			String [] externalSequenceNames = modelParams.getTree().getAllExternalSeqNames();
			List<String> seqOrder = new ArrayList<String>(externalSequenceNames.length);
			for(String s : externalSequenceNames) {
				seqOrder.add(s);
			}
			BufferedWriter bw = new BufferedWriter(new FileWriter(outFile));
			ma.write(bw, seqOrder);
			bw.close();

		}else {
			System.out.println(USAGE);
		}
	}

	protected MultipleAlignment generateAlignment(List<Integer> emittedPath, String format) {
		List<MarkovState<Map<String, Matrix>>> hmmStates = getStates();
		MultipleAlignment ma = MultipleAlignmentFactory.create(format);
		
		for(int stateIdx : emittedPath) {
			MarkovState<Map<String, Matrix>> state = hmmStates.get(stateIdx);
			Map<String, Matrix> col = state.emitObservation();
			ma.addShortEncodedColumnMatrix(col);
		}
		return ma;
	}

	public static PiHMM createDefaultTwoStateChain(double l, double gamma,EvolutionaryModel model) {
		PiHMM hmm = new PiHMM(model,2);
		ConstrainedModel cm = new ConstrainedModel(model,"constrained");
		hmm.addState(cm);
		hmm.setMuNuAndTransitions(l, gamma);
		hmm.setInitialStateTransitionProbability(0, 1-gamma);
		hmm.setInitialStateTransitionProbability(1, gamma);
		return hmm;
	}
	
	public static PiHMM createTwoStateChain(double l, double gamma, EvolutionaryModel model, Matrix otherPi) {
		PiHMM hmm = new PiHMM(model,2);
		PiStateModel newPiState = new PiStateModel(model, otherPi, "model2");
		hmm.addState(newPiState);
		hmm.setMuNuAndTransitions(l, gamma);
		hmm.setInitialStateTransitionProbability(0, 1-gamma);
		hmm.setInitialStateTransitionProbability(1, gamma);
		return hmm;
	}
	
	void useDefaultConservedStates() {
		for(int i = 0; i < conserved_pis.length; i++) {
			Matrix pi = new Matrix(conserved_pis[i].length,1);
			for(int j = 0; j < conserved_pis[i].length; j++) {
				pi.set(j, 0, conserved_pis[i][j]);
			}
			addConservedPi(pi, conserved_pi_names[i]);
		}
	}
	
	private void setMuNuAndTransitions(double l, double gamma) {
		mu = 1/l;
		nu = gamma*mu/((double)(getStates().size() - 1)*(1-gamma));
		
		setStateTransitionProbability(0, 0, 1 - (double)(getStates().size() - 1)*nu);
		for(int i = 1; i < getStates().size(); i++) {
			setStateTransitionProbability(0, i, nu);
			setStateTransitionProbability(i,0, mu);
			for(int j = 1; j < getStates().size(); j++) {
				setStateTransitionProbability(i,j,(1-mu)/(double)(getStates().size() - 1));
			}
		}
		

	}
	
}
