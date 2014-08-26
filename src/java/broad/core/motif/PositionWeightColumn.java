package broad.core.motif;

import java.io.BufferedWriter;
import java.io.IOException;
import java.text.NumberFormat;

import broad.core.math.MathUtil;

public class PositionWeightColumn {
	public static final int DEFAULT_ALPHABET_SIZE = 4; //Nucleotides
	
	private static final double log1over4 = Math.log(0.25);
	private static final double LOG2 = Math.log(2);
	private static final double DEFAULT_PSEUDOCOUNT = 0.01;
	
	private int [] num; // in the order of ACGTN
	private double [] prob;
	private double [] logProb;
	private char baseIUB;
	private double negEntropy;


	public PositionWeightColumn() {
		init();
	}
	
	public PositionWeightColumn(String [] rawData) {
		double [] freqs = new double[rawData.length];
		int  []   counts = new int[rawData.length];
		boolean isFrequencies = true;
		for(int i = 0; i < rawData.length; i++) {
			freqs[i] = Double.parseDouble(rawData[i]);
			counts[i] = (int) freqs[i];
			isFrequencies = isFrequencies && (freqs[i] <= 1 || !(counts[i] == freqs[i]));
		}
		if(isFrequencies) {
			init(freqs);
		} else {
			setFromCounts(counts);
		}
		
	}

	public PositionWeightColumn(int [] n) {
		setFromCounts(n);
	}

	public PositionWeightColumn(int nA, int nC, int nG, int nT) {
		init();
		setNum(nA, nC, nG, nT);	

		toConsensus();
		negEntropy = getInformationContent();
	}

	public PositionWeightColumn(double [] freq) {
		init(freq);
	}
	
	private void init(double [] freq) {
		prob = new double[freq.length];
		num = new int[prob.length];

		double sum = 0;
		for (int i=0; i<prob.length; i++) {
			num[i] = (int)Math.round(freq[i]*100);
			prob[i] = freq[i];
			sum += prob[i];
		}

		logProb = new double[prob.length];
		// normalize and log prob
		for (int i=0; i<prob.length; i++) {
			prob[i] = prob[i]/sum;
			logProb[i] = Math.log(prob[i]+Double.MIN_VALUE);
		}

		toConsensus();
		negEntropy = getInformationContent();
	}

	public char getBaseIUB() { return baseIUB;}
	
	public void resetProbWithNullPriorProb() {
		double [] bgProb = new double[getAlphabetSize()];  // used for pseudo count
		bgProb[0]=bgProb[3]=0.1; // A and T
		bgProb[1]=bgProb[2]=0.1; // C and G
		double sumBg = 0; 
		for (int i=0; i<bgProb.length; i++) sumBg += bgProb[i];

		int totalNum = num[0]+num[1]+num[2]+num[3];

		for (int i=0; i<getAlphabetSize(); i++) {
			prob[i] = (num[i]+bgProb[i])/(double)(sumBg+totalNum);
			logProb[i] = Math.log(prob[i]+Double.MIN_VALUE);
		}
		toConsensus();
		negEntropy = getInformationContent();
	}



	public void setNum(int [] n) {
		for (int i=0; i<n.length; i++)
			num[i]=n[i];

		computeProb();
	}

	public void setNum(int nA, int nC, int nG, int nT) {
		num[0]=nA;
		num[1]=nC;
		num[2]=nG;
		num[3]=nT;

		computeProb();
	}


	public PositionWeightColumn getComplement() {
		PositionWeightColumn op = new PositionWeightColumn();

		op.num[0] = num[3]; op.num[3] = num[0];
		op.num[1] = num[2]; op.num[2] = num[1];
		op.prob[0] = prob[3]; op.prob[3] = prob[0];
		op.prob[1] = prob[2]; op.prob[2] = prob[1];
		op.logProb[0] = logProb[3]; op.logProb[3] = logProb[0];
		op.logProb[1] = logProb[2]; op.logProb[2] = logProb[1];

		op.toConsensus();
		op.negEntropy = op.getInformationContent();
		return op;
	}


	public void init(int alphabetSize) {
		num = new int[alphabetSize];
		prob = new double[alphabetSize];
		logProb = new double[alphabetSize];

		for (int i=0; i<alphabetSize; i++) {
			num[i]=0; prob[i]=0.25; logProb[i]=Math.log(prob[i]);
		}
	}

	public void init() {
		init(DEFAULT_ALPHABET_SIZE);
	}

	public double getLogRatio(char ch) {
		return getLogProb(ch)-log1over4;
	}

	public double getLogRatioWeightedByNegEntropy(char ch) {
		return (getLogProb(ch)-log1over4)*negEntropy/2;
	}

	public double getLogProb(char ch) {
		if (ch=='A' || ch=='a') 
			return logProb[0];
		else if (ch=='C' || ch=='c') 
			return logProb[1];
		else if (ch=='G' || ch=='g') 
			return logProb[2]; 
		else if (ch=='T' || ch=='t') 
			return logProb[3];
		else return log1over4;
	}
	
	public double getLogProb(int baseIdx) {
		return logProb[baseIdx];
	}


	public double getLogRatio(char i, boolean isNumSeq) {
		return getLogProb(i)-log1over4;
	}


	public double getLogRatioWeightedByNegEntropy(char i, boolean isNumSeq) {
		return (getLogProb(i)-log1over4)*negEntropy/2;
	}
	
	protected void setFromCounts(int[] n) {
		init();
		setNum(n);

		toConsensus();
		negEntropy = getInformationContent();
	}

	public void computeProb() {
		double [] bgProb = new double[getAlphabetSize()];  // used for pseudo count
		bgProb[0]=bgProb[3]=0.25; // A and T
		bgProb[1]=bgProb[2]=0.25; // C and G

		int totalNum = num[0]+num[1]+num[2]+num[3];

		for (int i=0; i<4; i++) {
			prob[i] = (num[i]+bgProb[i])/(double)(1+totalNum);
			logProb[i] = Math.log(prob[i]+Double.MIN_VALUE);
		}

		toConsensus();
	}



	public void printProb() {
		System.out.print(baseIUB + "\t");
		for (int i=0; i<getAlphabetSize(); i++)
			System.out.print(prob[i] + "\t");

		System.out.println();
	}



	public void print() {
		System.out.print(baseIUB + "\t");
		for (int i=0; i<getAlphabetSize(); i++)
			System.out.print(num[i] + "\t");
		for (int i=0; i<getAlphabetSize(); i++)
			System.out.print(prob[i] + "\t");

		System.out.println();
	}


	public char toConsensus() {
		if (prob[0]>0.7)
			baseIUB = 'A';
		else if (prob[1]>0.7)
			baseIUB = 'C';
		else if (prob[2]>0.7)
			baseIUB = 'G';
		else if (prob[3]>0.7)
			baseIUB = 'T';

		else if (prob[0]+prob[3]>0.8)
			baseIUB = 'W';// A or T
		else if (prob[0]+prob[1]>0.8) 
			baseIUB = 'M';// A or C
		else if (prob[0]+prob[2]>0.8) 
			baseIUB = 'R'; // G or A
		else if (prob[3]+prob[1]>0.8) 
			baseIUB = 'Y'; // C or T
		else if (prob[3]+prob[2]>0.8) 
			baseIUB = 'K'; // G or T
		else if (prob[1]+prob[2]>0.8) 
			baseIUB = 'S'; // C or G

		else if (prob[0]<0.2) 
			baseIUB = 'B'; // C or G or T
		else if (prob[3]<0.2) 
			baseIUB = 'V'; // A or C or G
		else if (prob[1]<0.2) 
			baseIUB = 'D'; // A or G or T
		else if (prob[2]<0.2) 
			baseIUB = 'H'; // A or C or T
		else
			baseIUB = 'N';

		return baseIUB;
	}

	// return score log(alphabetSize)-information
	public double getInformationContent() {
		double score = MathUtil.log2(getAlphabetSize());

		for (int i=0; i<getAlphabetSize(); i++)
			score += prob[i]*logProb[i]/LOG2;
		return score;
	}



	public void flatWeight(double infoTh) {
		if (getInformationContent() <= infoTh) 	    resetToUniformProb();
	}


	public void resetToUniformProb() {
		double uniformLetterProb = 1/(double)getAlphabetSize();
		double logUniformLetterProb = MathUtil.log2(uniformLetterProb);
		
		for (int i=0; i<prob.length; i++) {
			prob[i] = uniformLetterProb;
			logProb[i] = logUniformLetterProb;
		}
		negEntropy = 0;
		baseIUB = 'N';
	}

	public double getWeight(int letterIndex) {
		return prob[letterIndex];
	}

	public void write(BufferedWriter bw, NumberFormat formatter) throws IOException {
		for (int i = 0; i < getAlphabetSize(); i++) {
			bw.write(formatter.format(prob[i]));
			if(i < prob.length - 1) {
				bw.write("\t");
			}
		}
	}
	
	public String toString() {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < getAlphabetSize(); i++) {
			sb.append(prob[i]);
			if(i < getAlphabetSize() - 1) {
				sb.append("\t");
			}
		}
		return sb.toString();
	}
	
	public double getProbability(int baseIndex) { return prob[baseIndex];}
	
	public int getAlphabetSize() { return prob.length;}

	public void addPseudoCounts() {
		double [] newProbs = new double[prob.length];
		double total=0;
		for(int i = 0; i < prob.length; i++) {
			newProbs[i] =  (prob[i] + DEFAULT_PSEUDOCOUNT) * 1000;
			total += newProbs[i];
		}
		
		for(int i = 0; i < prob.length; i++) {
			newProbs[i] =  newProbs[i]/total;
		}

		init(newProbs);	
	}

	public double kullbackLeiber(PositionWeightColumn pwc) {
		double kl = 0;
		for (int i = 0; i < pwc.getAlphabetSize(); i++) {
			kl += getProbability(i)*(getLogProb(i) - pwc.getLogProb(i) + pwc.getProbability(i)*(pwc.getLogProb(i)- getLogProb(i)));
		}
		return kl;
	}
 
}


