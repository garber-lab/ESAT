package broad.core.motif;

import java.io.BufferedWriter;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Random;

import org.apache.commons.math3.stat.clustering.Clusterable;

import Jama.Matrix;
import broad.core.annotation.BED;
import broad.core.error.ParseException;
import broad.core.sequence.SequenceRegion;
import broad.core.sequence.WindowSlider;

public class PositionWeightMatrix extends ArrayList<PositionWeightColumn> implements Clusterable<PositionWeightMatrix>{

	private static final long serialVersionUID = 2484880500988774252L;
	//Matrix matrix;
	String name;
	
	static final int ALPHABET_SIZE = 4;
	static final String [] ALPHABET = {"A","C","G","T"};
	
	//Xiaohui's fields
    private int leftHighInfoStart;
    private int rightHighInfoStart;
    private boolean alignDir;
    
	
	public PositionWeightMatrix(String name ){
		super();
		this.name = name;
	}
	
	public void setMatrixFromRawData(List<String []> unparsedRows) throws ParseException{
		int alphabetSize = unparsedRows.size();
		
		if(unparsedRows.size() > 0) {
			int colNum = unparsedRows.get(0).length;
			for(int i = 0; i < colNum; i++) {
				Iterator<String []> rowIt = unparsedRows.iterator();
				String [] col = new String[alphabetSize];
				int j = 0;
				while(rowIt.hasNext()) {					
					col[j++] = rowIt.next()[i];
				}
				setUnparsedRow(i, col);
			}
		}
	}
	
	public void setMatrix(List<double []> sites) throws ParseException{
		if(sites.size() > 0) {
			double [] firstColumn = sites.get(0);
			int alphabetSize = firstColumn.length;
			Iterator<double []> siteIt = sites.iterator();
			while(siteIt.hasNext()) {
				double [] site = siteIt.next();
				if(site.length != alphabetSize) {
					throw new ParseException("Not all site column lenghts are the same, first was " + alphabetSize + " but position was " + site.length);
				}
				add(new PositionWeightColumn(site));
			}
		}
	}
	
	public void addColumn(double [] site) {
		add(new PositionWeightColumn(site));
	}
	
	public void setUnparsedRow(int colNum , String [] unparsedRow) throws ParseException {
		PositionWeightColumn column = new PositionWeightColumn(unparsedRow);
		add(colNum, column);
	}
	
	public void write(BufferedWriter bw, NumberFormat formatter)  throws IOException {
		bw.write(">" + name);
		bw.newLine();
		if(size() > 0) {
			int rows = get(0).getAlphabetSize();
			for(int i = 0; i < rows; i++) {
				Iterator<PositionWeightColumn> colIt = iterator();
				while(colIt.hasNext()) {
					PositionWeightColumn col = colIt.next();
					bw.write(formatter.format(col.getWeight(i)));
					if(colIt.hasNext()) {
						bw.write("\t");
					}
				}
				bw.newLine();
			}
		}
	}
	
	public void writeUNIPROBE(BufferedWriter bw, NumberFormat formatter)  throws IOException {
		if(size() > 0) {
			bw.write(name);
			bw.newLine();
			int rows = get(0).getAlphabetSize();
			for(int i = 0; i < rows; i++) {
				bw.write(ALPHABET[i]);bw.write(":\t");
				Iterator<PositionWeightColumn> colIt = iterator();
				while(colIt.hasNext()) {
					PositionWeightColumn col = colIt.next();
					bw.write(formatter.format(col.getWeight(i)));
					if(colIt.hasNext()) {
						bw.write("\t");
					}
				}
				bw.newLine();
			}
		}
	}
	
	public String getName() { return name;}
	
	public void setName(String name) {  this.name = name;}
	
	
    public int getLeftHighInfoStart() {
		return leftHighInfoStart;
	}

	public int getRightHighInfoStart() {
		return rightHighInfoStart;
	}

	// set left and right pos with minAlignNum
	public void setAlignPos(int minAlignNum) {
		int [] sortedIdx =new int[minAlignNum];
		double [] sortedIC =new double[minAlignNum];
		for(int i=0;i<minAlignNum;i++){
			sortedIdx[i]=-1;
			sortedIC[i]=0.0;
		}

		// bubble sort find minAlignNum largest scores, store in b1v, index in b1  
		for (int i=0; i< size(); i++) {
			double iIC = get(i).getInformationContent();
			for (int j=0;j<minAlignNum;j++) {
				if (iIC>sortedIC[j]) {
					for (int k=minAlignNum-1;k>j;k--) {
						sortedIC[k]=sortedIC[k-1];
						sortedIdx[k]=sortedIdx[k-1];
					}
					sortedIC[j]=iIC;
					sortedIdx[j]=i;
					break;
				}
			}
		}

		leftHighInfoStart=size();
		rightHighInfoStart=-1;
		for(int i=0;i<minAlignNum;i++){
			if(sortedIdx[i]>rightHighInfoStart){
				rightHighInfoStart=sortedIdx[i];
			}
			if(sortedIdx[i]<leftHighInfoStart) {
				leftHighInfoStart=sortedIdx[i];
			}
		}

		// consider ties to the worst position among minAlignNum
		for(int i=leftHighInfoStart-1;i>=0;i--) {
			double iIC = get(i).getInformationContent();
			if ( Math.abs(iIC-sortedIC[minAlignNum-1])<.0001 ) leftHighInfoStart=i;
		}
		for (int i=rightHighInfoStart+1;i< size();i++) {
			double i1 = get(i).getInformationContent();
			if (Math.abs(i1-sortedIC[minAlignNum-1])<.0001) rightHighInfoStart=i;
		}
	}
	
	/**
	 * Trims this PWM by removing start end ending columns with information content lesser than given.
	 * @param ic - Columns below this ic at PWM edges will be trimmed.
	 * @return The resulting PWM.
	 */
	public PositionWeightMatrix trimByInformationContent(double ic) {
		PositionWeightMatrix trimmed = new PositionWeightMatrix(getName() );
		
		int lastPos = size() - 1;
		while(lastPos >=0 && get(lastPos).getInformationContent() < ic) {
			lastPos--;
		}
		
		int startPos = 0;
		while(startPos < lastPos && get(startPos).getInformationContent() < ic) {
			startPos++;
		}		
		
		for(int i = startPos; i <= lastPos; i++) {
			trimmed.add(get(i));
		}
		
		return trimmed;
	}
	
	/**
	 * Compute Distribution of scores for motif.
	 *
	 * @param s
	 * @param isNumSeq
	 * @return
	 */
	public List<Double> computeScoreDistribution(float backgroundA, float backgroundC, float backgroundG, float backgroundT, int sampleSize) {
		double [] bg = {backgroundA,backgroundC, backgroundG, backgroundT};
		PositionWeightMatrix bgPWM = createIsoPWM(bg, "bg");
		//System.out.println("bg vector " + backgroundA + "," + backgroundC + "," + backgroundG + "," + backgroundT);
		List<Double> dist = new ArrayList<Double>( sampleSize);
		int [] kmer = new int[size()];
		for(int i = 0; i < sampleSize ; i++) {
			Random r = new Random();
			for(int j = 0; j < size(); j++) {
				kmer[j] = r.nextInt(ALPHABET_SIZE);
			}
			dist.add(getLogLikelihood(kmer) - bgPWM.getLogLikelihood(kmer));
			//System.out.println(printKmer(kmer)+ " -- " + (getLogLikelihood(kmer) - bgPWM.getLogLikelihood(kmer)));
		}
		return dist;
	}
	public String printKmer(int [] kmer) {
		StringBuilder sb = new StringBuilder("(");
		for(int i = 0; i < kmer.length; i++) {
			sb.append(kmer[i]).append(",");
		}
		sb.append(")");
		return sb.toString();
	}
	
	// log likelihood given by numerical seq: 0123(ACGT)
	public double getLogLikelihood(char [] s) {
		double lh = 0;
		for (int i=0; i<size(); i++) {
			lh += get(i).getLogProb(s[i]);
		}
		return lh;
	}
	
	public double getLogLikelihood(int [] s) {
		double lh = 0;
		for (int i=0; i<size(); i++) {
			lh += get(i).getLogProb(s[i]);
		}
		return lh;
	}
	
	public double getLogLikelihood(short[] encodedSequence, int i) {
		double lh = 0;
		for (int j=0; j<size(); j++) {
			//System.err.println("log prob of " + i + " val " + encodedSequence[j+i] + " is "  + get(j).getLogProb(encodedSequence[j+i]));
			lh += get(j).getLogProb(encodedSequence[j+i]); 
		}
		return lh;
	}

	/**
	 * Gets the next logical kmer. If k = 4 the first 5 4mers would be:
	 * (0,0,0,0), (0,0,0,1), (0,0,0,2), (0,0,0,3), (0,0,1,0).
	 * @param kmer - any valid kmer
	 * @return Next logical kmer.
	 */
	public static int []  getNextKmer(int [] kmer) {
		int [] next = kmer.clone();
		int alphabetSize = 4;
		int i = kmer.length - 1;
		while(i >= 0) {
			if(next[i] < alphabetSize - 1) {
				next[i]++;
				break;
			} else {
				next[i] = 0;
				i--;
			}
		}
		return next;
	}
	
	public int getNumCol() {
		return size();
	}
	
	public boolean isEmpty() { return size() == 0;}
	
	public PositionWeightMatrix permuteColumns(boolean preserveGCDinucleotides) {
		PositionWeightMatrix original = copy();
		List<Integer> dinucleotidesToPreserveFirstIdx = preserveGCDinucleotides ? original.gcDinucleotidesFistIdxs() : new ArrayList<Integer>();
		PositionWeightMatrix permutted = new PositionWeightMatrix(getName() + "_perm");
		Random r = new Random();
		List<Integer> idxList = new ArrayList<Integer>(original.size());
		for(int i = 0; i < original.size(); i++) {
			idxList.add(i);
		}
		while(idxList.size() > 0) {
			int idxOfIdx = r.nextInt(idxList.size());
			int idx = idxList.remove(idxOfIdx);
			if(dinucleotidesToPreserveFirstIdx.contains(idx-1)) {
				PositionWeightColumn col = original.get(idx-1);
				permutted.add(col);
				PositionWeightColumn nextCol = original.get(idx);
				idxList.remove(idxOfIdx - 1);
				permutted.add(nextCol);
			}else{ 
				PositionWeightColumn col = original.get(idx);
				permutted.add(col);			
				if(dinucleotidesToPreserveFirstIdx.contains(idx)) {
					PositionWeightColumn nextCol = original.get(idx+1); 
					idxList.remove(idxOfIdx);
					permutted.add(nextCol);
				}
			}

		}
		
		return permutted;
	}
	
	private List<Integer> gcDinucleotidesFistIdxs() {
		List<Integer> dinucleotidesToPreserveFirstIdx = new ArrayList<Integer>();
		for(int i = 0; i < size()-1; i++) {
			PositionWeightColumn thisCol = get(i);
			PositionWeightColumn nextCol = get(i+1);
			
			if(thisCol.getLogRatio('G')> 0 && nextCol.getLogRatio('C')>0) {
				dinucleotidesToPreserveFirstIdx.add(i);
			}
		}
		return dinucleotidesToPreserveFirstIdx;
	}

	public PositionWeightMatrix copy() {
		PositionWeightMatrix copy = new PositionWeightMatrix(getName());
		for(PositionWeightColumn c : this) {
			copy.add(c);
		}
		
		return copy;
		
	}



	/**
	 * Computes the information content of this position weith matris
	 * @return
	 */
	public double ic() {
		double ic = 0;
		for(PositionWeightColumn col : this) {
			ic += col.getInformationContent();
		}
		
		return ic;
	}

	public PositionWeightMatrix reverseComplement() {
		PositionWeightMatrix wm = new PositionWeightMatrix(getName());
		int numCol = size();
		for (int i=0; i<numCol; i++)
			wm.add(get(numCol-1-i).getComplement());

		wm.name = name;
		wm.rightHighInfoStart = numCol-1-leftHighInfoStart;
		wm.leftHighInfoStart = numCol-1-rightHighInfoStart;

		return wm;
	}
  
	public int getNumEndPosToErase(double infoThreshold) {
		double score;
		int num = 0;	
		int minNumCol = 6;
		for (int i=size()-1; i>=minNumCol; i--) {
			score = get(i).getInformationContent();
			if (score > infoThreshold) break;
			num++;
		}

		return num;
	}


	// get num of beginning position with low information content to erase
	public int getNumStartPosToErase(double infoThreshold) {
		double score;
		int num = 0;
		int minNumCol = 6;

		for (int i=0; i<size()-minNumCol; i++) {
			score = get(i).getInformationContent();
			if (score > infoThreshold) break;
			num++;
		}

		return num;
	}
    
    
	// only print refined weightmatrix,
	public void printRefinedWeightMatrix(double infoThreshold) {
		int startNumToDel = getNumStartPosToErase(infoThreshold);
		int endNumToDel = getNumEndPosToErase(infoThreshold);
		int newMotifLength = size() - startNumToDel - endNumToDel;
		if (newMotifLength < 6) newMotifLength = 6; // minimum number of positions
		newMotifLength = Math.min(newMotifLength, size()-startNumToDel);

		PositionWeightColumn [] tmp = new PositionWeightColumn[newMotifLength];
		
		for (int i=0; i<tmp.length; i++)
			tmp[i] = get(i+startNumToDel);

		for (int i=0; i<tmp.length; i++)
			tmp[i].print();
	}



	public void flatWeightMatrix(double infoTh) {
		for (int i=0; i<size(); i++)
			get(i).flatWeight(infoTh);
	}

	public String getConsensus() {
		char [] consensus = new char[size()];

		for (int i=0; i<size(); i++)
			consensus[i] = get(i).getBaseIUB();

		return new String(consensus);
	}


	public boolean getAlignDir() {
		return alignDir;
	}

	public List<BED> match(SequenceRegion region, double[] background, float minScore) {
		PositionWeightMatrix neutralPWM = createIsoPWM(background, "neutral");
		WindowSlider slider = WindowSlider.getSlider(region, size(), size() - 1);
		List<BED> scoredWindows = new ArrayList<BED>();
		
		while(slider.hasNext()) {
			SequenceRegion window = slider.next();
			//System.out.println(window + " seq: " + window.getSequenceBases());
			char[] windowChrs = window.getSequenceBases().toCharArray();
			double directScore = getLogLikelihood(windowChrs) - neutralPWM.getLogLikelihood(windowChrs);
			window.reverse();
			char [] reversedChrs = window.getSequenceBases().toCharArray();
			double reverseScore = getLogLikelihood(reversedChrs) - neutralPWM.getLogLikelihood(reversedChrs);
			double max = Math.max(directScore, reverseScore);
			if(max >= minScore) {
				BED scoredWindow = new BED(window);
				scoredWindow.setStart(scoredWindow.getStart() + region.getStart());
				scoredWindow.setEnd(scoredWindow.getEnd() + region.getStart());
				scoredWindow.setOrientation(directScore > reverseScore);
				scoredWindow.setScore(max);
				scoredWindow.setChromosome(region.getContainingSequenceId());
				scoredWindows.add(scoredWindow);
			}
		}
		return scoredWindows;
	}
	
	/**
	 * Creates a PWM of similar dimension to this PWM but with all columns set to the given vector,
	 * usually a neutral mutation vector.
	 * @param column
	 * @return
	 */
	public PositionWeightMatrix createIsoPWM(double [] column, String name) {
		PositionWeightMatrix pwm = new PositionWeightMatrix(name);
		 for(int i = 0; i < size(); i++) {
			 pwm.addColumn(column);
		 }
		 
		 return pwm;
	}

	public void addPseudoCounts() {
		for(PositionWeightColumn pwc : this) {
			pwc.addPseudoCounts();
		}
		
	}

	/**
	 * Computes the euclidean centroid, it may not be the formal centroid for different metrics but intuitively the
	 * average counts should provide a good cluster representative which is what this method intends to return
	 * @param pwmSet - Collection of pwms from which to compute the centroid0
	 * @return The euclidean centroid
	 * @throws IllegalArgumentException - When not all PWMs have the same dimension.
	 */
	public PositionWeightMatrix centroidOf(Collection<PositionWeightMatrix> pwmSet) throws IllegalArgumentException{
		Matrix centroidMatrix = null;
		Iterator<PositionWeightMatrix> pwmIt = pwmSet.iterator();
		if(pwmIt.hasNext()) {
			PositionWeightMatrix first = pwmIt.next();
			PositionWeightColumn firstCol = first.get(0);
			centroidMatrix = new Matrix(firstCol.getAlphabetSize(), first.getNumCol());
			addToCentroid(centroidMatrix, first);
		}
		while(pwmIt.hasNext()) {
			PositionWeightMatrix pwm = pwmIt.next();
			if (pwm.getNumCol() != centroidMatrix.getColumnDimension()) {
				throw new IllegalArgumentException("Error computing centroid. All PWMs in set should have the same dimension");
			}
			addToCentroid(centroidMatrix, pwm);
		}
		
		centroidMatrix.times(1/(double)pwmSet.size());
		
		PositionWeightMatrix centroid = new PositionWeightMatrix("centroid");
		for(int j = 0; j < centroidMatrix.getColumnDimension(); j++) {
			centroid.addColumn(centroidMatrix.getColumn(j));
		}
		return centroid;
	}

	protected void addToCentroid(Matrix m, PositionWeightMatrix pwm) {
		for(int i = 0; i < m.getColumnDimension(); i++) {
			PositionWeightColumn c = pwm.get(i);
			for (int j = 0; j < m.getRowDimension(); j++) {
				m.set(j, i, m.get(j, i) + c.getWeight(j));
			}
		}
	}

	/**
	 * Distance between two PWMs, uses the KL similarity
	 */
	public double distanceFrom(PositionWeightMatrix other) {

		return kullbackLeiber(other);
	}

	public double kullbackLeiber(PositionWeightMatrix other) {
		double kl = 0;
		for(int i = 0; i < other.getNumCol(); i++) {
			kl += get(i).kullbackLeiber(other.get(i));
		}
		return kl;
	}



}
