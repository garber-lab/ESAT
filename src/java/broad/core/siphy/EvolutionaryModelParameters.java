package broad.core.siphy;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.StringWriter;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.parsers.nhx.NHXParser;

import Jama.Matrix;
import broad.core.error.ParseException;
import broad.core.util.CLUtil;

public class EvolutionaryModelParameters {
	private Phylogeny tree;
	private double kappa;
	private double mu;
	private Matrix rateMatrix;
	private Matrix backgroundNucleotideFreqs;
	
	private static final String [] defaultAlphabet = {"A", "C", "G", "T" };
	
	public EvolutionaryModelParameters(Phylogeny tree, double kappa, double mu, Matrix Q, Matrix pi) {
		this.tree = tree;
		this.kappa = kappa;
		this.mu = mu;
		this.rateMatrix = Q;
		this.backgroundNucleotideFreqs = pi;
	}
	
	public EvolutionaryModelParameters(File sourceFile) throws IOException, ParseException {
		BufferedReader br = new BufferedReader(new FileReader(sourceFile));
		this.mu = 1;
		try {
			String line = null;
			boolean atMatrix = false;
			int rateMatrixRow = 0;
			Matrix alphabetTransform = Matrix.identity(4, 4);
			
			while((line = br.readLine()) != null) {
				line = line.replaceFirst("^\\s+", "").trim();
				if(line.startsWith("#") || line.trim().length() == 0) {
					continue;
				}
				if(line.startsWith("ALPHABET")) {
					atMatrix = false;
					line = line.replaceFirst("ALPHABET:\\s*", "").trim();
					String [] lineInfo = line.contains(",") ? line.split(",") : line.split("\\s+");
					
					if(lineInfo.length != defaultAlphabet.length) {
						throw new ParseException("Syntax error in alphabet provided, there should be " + defaultAlphabet.length + " comma separated letters but got " + line);
					}
					
					for (int j = 0; j < lineInfo.length; j++) {
						String letter = lineInfo[j];
						boolean found = false;
						for(int i = 0; i < defaultAlphabet.length; i++) {
							if(defaultAlphabet[i].equalsIgnoreCase(letter)) {
								found = true;
								alphabetTransform.set(i, j, 1);
							} else {
								alphabetTransform.set(i, j, 0);
							}
						}
						if(!found) {
							throw new ParseException("Alphabet provided has unknown letter " + letter);
						}
					}
					
					//System.out.print("alphabet transform: " );
					//alphabetTransform.print(4, 4);
				}else if(line.startsWith("BACKGROUND:")) {
					atMatrix = false;
					line = line.replaceFirst("BACKGROUND:\\s*", "").trim();
					String [] lineInfo = line.split("\\s+");
					backgroundNucleotideFreqs = new Matrix(lineInfo.length, 1);
					/*
					if(lineInfo.length != 4) {
						throw new ParseException("BACKGROUND record in model contains " + lineInfo.length + " rathern than 4 as it should");
					}
					*/
					for(int i = 0; i < backgroundNucleotideFreqs.getRowDimension(); i++) {
						backgroundNucleotideFreqs.set(i, 0, Double.parseDouble(lineInfo[i])) ;
					}
				} else if (line.startsWith("TREE:")) {
					line = line.replaceFirst("TREE:\\s*", "");
					NHXParser parser =  new NHXParser();
					parser.setSource(line);
					Phylogeny [] foundPhylogenies = parser.parse();
					tree = foundPhylogenies[0];
					atMatrix = false;
				} else if (line.startsWith("RATE_MAT:")) {
					atMatrix = true;
					rateMatrix = new Matrix(4,4);
				} else if (atMatrix) {
					String [] lineInfo = line.split("\\s+");
					if(lineInfo.length != 4) {
						throw new ParseException("Rate Matrix line " + line + " has " + lineInfo.length + " elements rather than 4 as it should");
					}
					for(int i = 0; i < 4; i++) {
						rateMatrix.set(rateMatrixRow, i, Double.parseDouble(lineInfo[i]));
					}
					rateMatrixRow++;
				} else if (line.startsWith("KAPPA:")) {
					atMatrix = false;
					line = line.replaceFirst("KAPPA:\\s*","").trim();
					this.kappa = Double.parseDouble(line);
				} else if (line.startsWith("MU:")) {
					line = line.replaceFirst("MU:\\s*","").trim();
					this.mu = Double.parseDouble(line);					
				}
				
			}
			
			if(backgroundNucleotideFreqs != null) {
				backgroundNucleotideFreqs = alphabetTransform.times(backgroundNucleotideFreqs);
			}
			
			if(rateMatrix != null) {
				rateMatrix = alphabetTransform.times(rateMatrix).times(alphabetTransform.transpose());
			}
		} finally {
			br.close();
		}
		
		if(rateMatrix == null  && kappa != 0 && backgroundNucleotideFreqs != null) {
			generateHKYQMatrix();
		}
		
	}
	
	public EvolutionaryModelParameters copy() {
		return  new EvolutionaryModelParameters(tree, kappa, mu, rateMatrix, backgroundNucleotideFreqs);
	}
	
	public EvolutionaryModelParameters(double kappa, double[] bg, double mu) {
		backgroundNucleotideFreqs = new Matrix(bg.length, 1);
		
		for(int i = 0; i < bg.length; i++) {
			backgroundNucleotideFreqs.set(i, 0, bg[i]);
		}
		generateHKYQMatrix(kappa, mu);
	}
	
	private void generateHKYQMatrix() {
		generateHKYQMatrix(kappa, mu);
	}

	public void generateHKYQMatrix(double kappa, double scalingFactor ) {
		rateMatrix = new Matrix(4,4);
		
		rateMatrix.set(0, 0, 0);
		rateMatrix.set(0,1,1);
		rateMatrix.set(0,2,kappa);
		rateMatrix.set(0,3,1);
		
		rateMatrix.set(1,0,1);
		rateMatrix.set(1,1,0);
		rateMatrix.set(1,2,1);
		rateMatrix.set(1,3,kappa);
		
		rateMatrix.set(2,0,kappa);
		rateMatrix.set(2,1,1);
		rateMatrix.set(2,2,0);
		rateMatrix.set(2,3,1);
		
		rateMatrix.set(3,0,1);
		rateMatrix.set(3,1,kappa);
		rateMatrix.set(3,2,1);
		rateMatrix.set(3,3,0);
		
		//rateMatrix.print(4, 5);
		
		rateMatrix.timesEquals(scalingFactor);
		//rateMatrix.print(4, 5);
		
		Matrix diagBG = Matrix.identity(4, 4);
		for(int i = 0; i < 4; i++) {
			diagBG.set(i, i, backgroundNucleotideFreqs.get(i, 0));
		}
		
		rateMatrix = rateMatrix.times(diagBG);
		
		for(int i = 0; i < rateMatrix.getColumnDimension(); i++) {
			double diagVal = 0;
			for(int j = 0; j < rateMatrix.getRowDimension(); j++) {
				if(i != j) {
					diagVal -= rateMatrix.get(i, j);
				}
			}
			
			rateMatrix.set(i, i, diagVal);
		}
		
	}
	
	public void write (BufferedWriter bw) throws IOException {
		if(backgroundNucleotideFreqs != null) {
			bw.write("BACKGROUND:");
			//backgroundNucleotideFreqs.print(30,20);
			for(int i = 0; i < backgroundNucleotideFreqs.getRowDimension(); i++) {
				bw.write(" " + backgroundNucleotideFreqs.get(i,0));
			}
			bw.newLine();
		}
		
		if(rateMatrix != null) {
			bw.write("RATE_MAT:\n");
			for(int i = 0; i < rateMatrix.getRowDimension(); i++) {
				for(int j = 0; j < rateMatrix.getColumnDimension(); j++) {
					CLUtil.writeRightJustifiedField(bw, String.valueOf(rateMatrix.get(i,j)), 20);
				}
				bw.newLine();
			}
		}
		
		if(tree != null) {
			bw.write("TREE: ");
			String treeStr = tree.toNewHampshire(true);
			//treeStr.replaceAll("\\Q,)\\E", ",");
			bw.write(treeStr);
			bw.newLine();
		}		
	}
	
	public String toString() {
		StringWriter sw = new StringWriter();
		try {
			BufferedWriter bw = new BufferedWriter(sw);
			write(bw);
			bw.flush();
		} catch (IOException ioe) {
			System.err.println("IOException while writing to String Writer... ");
			ioe.printStackTrace(System.err);
		}
		return sw.toString();
	}

	public double[] getBackgroundNucleotideFreqs() {		
		return backgroundNucleotideFreqs.getColumnPackedCopy();
	}

	public Matrix getRateMatrix() {
		if(rateMatrix == null && kappa > 0 && backgroundNucleotideFreqs != null) {
			generateHKYQMatrix();
		}
		return rateMatrix;
	}

	public Phylogeny getTree() {
		return tree;
	}
	
	public void setTree(Phylogeny tree) {
		this.tree = tree;
	}

	public void setRateMatrix(Matrix Q) {
		rateMatrix = Q;
	}
	
	public void setBackgroundNucleotideFreqs(Matrix dist) {
		backgroundNucleotideFreqs = dist;
	}

}
