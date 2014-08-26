package broad.core.motif.meme;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;

import broad.core.motif.PositionWeightMatrix;
import broad.core.sequence.SequenceRegion;
import jaligner.Sequence;

public class MemeMotif {
	PositionWeightMatrix pwm;
	List<MatchedSequence> matches;
	int width;
	int sites;
	float llr;
	double eval;
	
	public MemeMotif(String string) {
		pwm = new PositionWeightMatrix(string);
		matches = new ArrayList<MatchedSequence>();
	}



	void setPwm(PositionWeightMatrix pwm) {
		this.pwm = pwm;
	}


	public String getName() { return pwm.getName();}

	void setMatches(List<MatchedSequence> matches) {
		this.matches = matches;
	}



	void setWidth(int width) {
		this.width = width;
	}



	void setSites(int sites) {
		this.sites = sites;
	}



	void setLlr(float llr) {
		this.llr = llr;
	}



	void setEval(double eval) {
		this.eval = eval;
	}



	public List<MatchedSequence> getMatches() {
		return matches;
	}



	public int getWidth() {
		return width;
	}



	public int getSites() {
		return sites;
	}



	public float getLlr() {
		return llr;
	}



	public double getEval() {
		return eval;
	}



	public void addMatch(Sequence match, int startInScannedSequence, double pval) {
		matches.add(new MatchedSequence(match, startInScannedSequence, pval));
	}
	
	public void addPWMColumn(double [] col) {
		pwm.addColumn(col);
	}
	
	public PositionWeightMatrix getPWM() { return pwm;}
	
	public static class MatchedSequence  {

		private int startOfMatch;
		private Sequence sequence;
		private double pvalue;
		MatchedSequence(Sequence seq,  int startOfMatch, double pval) {
			this.sequence = seq;
			this.startOfMatch = startOfMatch;
			this.pvalue = pval;
		}
		
		public int getStartOfMatch() {
			return startOfMatch;
		}
		public Sequence getSequence() {
			return sequence;
		}
		public double getPvalue() {
			return pvalue;
		}
	}

}
