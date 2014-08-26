package broad.core.siphy;

import java.util.List;
import java.util.Map;

import broad.core.motif.PositionWeightColumn;
import broad.core.motif.PositionWeightMatrix;

import Jama.Matrix;

public class PositionWeightMatrixModel {
	EvolutionaryModel [] inducedModel;
	EvolutionaryModel baseModel;
	public PositionWeightMatrixModel(PositionWeightMatrix pwm, EvolutionaryModel baseModel) {
		this.baseModel = baseModel;
		inducedModel = new EvolutionaryModel[pwm.size()];
		for(int k = 0; k < pwm.size(); k++) {
			EvolutionaryModel mod = baseModel.copy();
			Matrix pi = new Matrix(mod.getAlphabetSize(),1);
			PositionWeightColumn col = pwm.get(k);
			for(int l = 0; l < col.getAlphabetSize(); l++) {
				pi.set(l, 0, col.getProbability(l));
			}
			mod.RPIDecomposition();
			mod.setPi(pi);
			inducedModel[k] = mod;
		}
	}
	public double score(Map<String, Matrix> alignmentWindow) {
		double sumLODS = 0;
		if(alignmentWindow == null || alignmentWindow.size() == 0) {
			return 0;
		}
		
		int alignmentWindowLength = alignmentWindow.values().iterator().next().getColumnDimension();
		if(alignmentWindowLength != inducedModel.length) {
			throw new IllegalArgumentException("Alignment window was " + alignmentWindowLength + " while PWM was of length " + inducedModel.length + " they must be of same length");
		}
		for(int i = 0; i < alignmentWindowLength; i++) {
			List<String> gappedSeqs = ConservationUtils.getGappedSeqsInWindowMatrix(1, alignmentWindow, i);						
			ConservationUtils.setUninformativeNodes(alignmentWindow, gappedSeqs, i);
			//Phylogeny columnTree = ConservationUtils.pruneTree(gappedSeqs, baseModel.);
			//baseModel.clearComputedTransitionsCache();
			//double neutralLikelihood = baseModel.pruneAndPeel(alignmentWindow, baseModel.getTree().getRoot(), i);
			double neutralLikelihood = baseModel.computeLikelihood(alignmentWindow, baseModel.getTree().getRoot(), i);
			//inducedModel[i].clearComputedTransitionsCache();
			//double colLikelihood     = inducedModel[i].pruneAndPeel(alignmentWindow, inducedModel[i].getTree().getRoot(), i);
			double colLikelihood     = inducedModel[i].computeLikelihood(alignmentWindow, inducedModel[i].getTree().getRoot(), i);
			//System.out.println("\tbase L " + neutralLikelihood + " pwm " + colLikelihood);			
			sumLODS += Math.log10(colLikelihood/neutralLikelihood);// - Math.log10(neutralLikelihood);
		}
		//System.out.println("Score: " + sumLODS);
		return sumLODS;
	}
}
