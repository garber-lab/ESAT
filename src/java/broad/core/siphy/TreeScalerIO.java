package broad.core.siphy;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;
import java.util.Random;

import broad.core.annotation.BasicGenomicAnnotation;
import broad.core.math.EmpiricalDistribution;

public class TreeScalerIO {
	int estimationWindow;
	ArrayList<ScaledWindow> scalings;
	HashMap<Double, Integer> bins;
	private static final DecimalFormat DEFAULT_FORMATTER = new  DecimalFormat("##0.####"); 
	private DecimalFormat formatter = DEFAULT_FORMATTER;
	int totalInSources;
	private List<ScaledWindowParserListener>  listeners;
	
	
	public TreeScalerIO (int windowSize) {
		super();
		estimationWindow = windowSize;
		scalings = new ArrayList<ScaledWindow>();
		listeners = new ArrayList<ScaledWindowParserListener>();
	}
	
	public void loadIntoDistribution(File source, EmpiricalDistribution ed, boolean shortFormat) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(source));
		try {
			String line = null;
			while((line = br.readLine()) != null) {
				if(line.startsWith("#") || line.trim().length() == 0) {
					continue;
				}
				String [] lineInfo = line.trim().split("\t");
				try {
					if(shortFormat) {
						ed.add(Double.parseDouble(lineInfo[1]));
					} else {
						ed.add(Double.parseDouble(lineInfo[3]));
					}
				} catch(NumberFormatException nfe) {
					System.out.println("Uncomputed omega");
					continue;
				}
			}
		} finally {
			br.close();
		}
	}

	
	public TreeScalerIO() {
		this(1);
	}
	
	public void parse(File source, int shift, String chr) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(source));
		try {
			String line = null;
			while((line = br.readLine()) != null) {
				if(line.startsWith("#") || line.trim().length() == 0) {
					continue;
				}
				String [] lineInfo = line.trim().split("\t");
				ScaledWindow w = new ScaledWindow(lineInfo, chr, estimationWindow, shift);
				Iterator<ScaledWindowParserListener> listenerIt = listeners.iterator();
				while(listenerIt.hasNext()) {
					ScaledWindowParserListener listener = listenerIt.next();
					listener.newWindow(w);
				}
				
			}
		} finally {
			br.close();
		}
	}		
	
	public void load(InputStream is, int shift, String chr) throws IOException {
		BufferedReader br = new BufferedReader(new InputStreamReader(is));
		try {
			String line = null;
			while((line = br.readLine()) != null) {
				if(line.startsWith("#") || line.trim().length() == 0) {
					continue;
				}
				String [] lineInfo = line.trim().split("\t");
				try {
					scalings.add(new ScaledWindow(lineInfo, chr,  estimationWindow, shift));
				} catch (IllegalArgumentException iae) {
					System.err.println("ERROR: line " + line + " error: " + iae.getMessage() );
					
				}
				
			}
		} finally {
			br.close();
		}		
	}
	
	public void load(File source, int shift, String chr) throws IOException {
		load(new FileInputStream(source), shift, chr);
		
	}
	
	public void load(InputStream is) throws IOException {
		load(is, 0, "C");
	}
	
	public void loadIntegrated(File source, double minPVal) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(source));
		try {
			String line = null;
			while((line = br.readLine()) != null) {
				if(line.startsWith("#") || line.trim().length() == 0) {
					continue;
				}
				String [] lineInfo = line.trim().split("\t");
				ScaledWindow w = new ScaledWindow(lineInfo);
				totalInSources++;
				if(w.omegaPValue < minPVal) {
					scalings.add(w);
				}
								
				
			}
		} finally {
			br.close();
		}
		
		System.out.println("\tloaded so far " + scalings.size() + " seen " + totalInSources);
	}
	
	public void loadOldIntegrated(File source, int shift, String chr, double minPVal) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(source));
		try {
			String line = null;
			while((line = br.readLine()) != null) {
				if(line.startsWith("#") || line.trim().length() == 0) {
					continue;
				}
				String [] lineInfo = line.trim().split("\t");
				ScaledWindow w = new ScaledWindow(chr);
				w.setStart(Integer.parseInt(lineInfo[0]) + shift);
				w.setEnd(w.getStart() + estimationWindow);
				w.setChromosome(chr);
				w.setOmega(Double.parseDouble(lineInfo[1]));
				w.setOmegaPValue(Double.parseDouble(lineInfo[4]));
				totalInSources++;
				if(w.omegaPValue < minPVal) {
					scalings.add(w);
				}
								
				
			}
		} finally {
			br.close();
		}
		
		System.out.println("\tloaded so far " + scalings.size() + " seen " + totalInSources);
		
	}
	
	public void addListener(ScaledWindowParserListener listener) {
		listeners.add(listener);
	}
	
	public void loadAtMostRandom(File source, String chr, int maxNeutralOmegas) throws IOException{
		BufferedReader br = new BufferedReader(new FileReader(source));
		Random r = new Random();
		try {
			String line = null;
			while((line = br.readLine()) != null) {
				if(line.startsWith("#") || line.trim().length() == 0) {
					continue;
				}
				String [] lineInfo = line.trim().split("\t");
				ScaledWindow window = new ScaledWindow(lineInfo, chr,  estimationWindow, 0);
				if(scalings.size() == maxNeutralOmegas ) {
					int index = r.nextInt(maxNeutralOmegas + 10);
					if(index < maxNeutralOmegas) {
						scalings.set(index, window);
					}
				} else {
					scalings.add(window);
				}				
			}
			
		} finally {
			br.close();
		}		
	}
	
	
	public void stitch(int window, int skip) {
		ListIterator<ScaledWindow> it = scalings.listIterator();
		while(it.hasNext()) {
			ScaledWindow current = it.next();
			List<ScaledWindow> toStitch = new ArrayList<ScaledWindow>(window - 1);
			
			boolean jump = false;
			while(toStitch.size() < window - 1 && !jump && it.hasNext()) {
				ScaledWindow nextInWindow = it.next();
				if(nextInWindow.getStart() == current.getEnd() + toStitch.size()) {
					toStitch.add(nextInWindow);
				} else {
					jump = true;
				}
			}
			
			if(!jump && toStitch.size() == window - 1) {
				current.stitchToList(toStitch);
				for(int i = 0; i < window - 1; i++) {
					it.previous();
				}
				int deleted = 0;
				while(deleted < skip - 1 && it.hasNext()) {
					it.next();
					it.remove();
					deleted++;
				}
			} else {
				if(jump) {
					it.previous(); //If there was a jump we do not want to remove the window onto which we jumped.
				}
				for(int i = 0; i <= toStitch.size(); i++) {
					it.previous();
					it.remove();
				}
				//it.previous();
				//it.remove();
			}
 		}
	}
	
	public void writeFull(BufferedWriter bw, List<ScaledWindow> windows) throws IOException {
		Iterator<ScaledWindow> it = windows.iterator();
		while(it.hasNext()) {
			ScaledWindow sw = it.next();
			bw.write(sw.toFullString(formatter));
			bw.newLine();
		}
	}
	
	public void writeFull(BufferedWriter bw) throws IOException{
		writeFull(bw, scalings);
	}
	
	public void writeFullAsGenomeGraphUsingLogPval(BufferedWriter bw) throws IOException {
		Iterator<ScaledWindow> it = scalings.iterator();
		while(it.hasNext()) {
			ScaledWindow sw = it.next();
			bw.write(sw.toMidBaseLogPvalScoredGGraph(formatter));
			bw.newLine();
		}	
	}
 	
	
	public static class ScaledWindow extends BasicGenomicAnnotation implements Fit{
		double omega;
		double treeLength;
		double omegaPValue;
		double logOdds;
		double logOddsPValue;
		
		public ScaledWindow(String[] lineInfo, String chr, int windowSize, int shift) {
			super();
			if(lineInfo.length < 3 ) {
				StringBuilder buf = new StringBuilder("Line contained to few fields: ");
				for (int i = 0; i < lineInfo.length; i++) {
					buf.append(lineInfo[i]);
					if(i < lineInfo.length - 1) {
						buf.append(", ");
					}
				}
				throw new IllegalArgumentException(buf.toString());
			}
			
			
			if(lineInfo.length <= 5 ) {
				setChromosome(chr);
				setStart(Integer.parseInt(lineInfo[0]) + shift);
				setEnd(getStart() + windowSize);
				setOmega(Double.parseDouble(lineInfo[1]));
				setTreeLength(Double.parseDouble(lineInfo[2]));
				if(lineInfo.length >= 4) {
					setLogOdds(Double.parseDouble(lineInfo[3]));
				}
				if(lineInfo.length >=5) {
					if(!lineInfo[4].contains("?") && !lineInfo[4].contains("NaN")) {
						setLogOddsPValue(Double.parseDouble(lineInfo[4]));
					}
				}
			} else if (lineInfo.length <= 7) {
				if(chr != null) {
					setChromosome(chr);
				} else {
					setChromosome(lineInfo[0]);
				}
				setStart(Integer.parseInt(lineInfo[1]) + shift);
				setEnd(Integer.parseInt(lineInfo[2]) + shift);
				setOmega(Double.parseDouble(lineInfo[3]));
				//setOmegaPValue(Double.parseDouble(lineInfo[4]));
				if(lineInfo.length == 7) {
					setLogOdds(Double.parseDouble(lineInfo[5]));
					setLogOddsPValue(Double.parseDouble(lineInfo[6]));
				}
			}
			
		}
		public ScaledWindow(String[] lineInfo) {
			super("", lineInfo[0].replace("chr", ""), Integer.parseInt(lineInfo[1]), Integer.parseInt(lineInfo[2]));
			omega = Double.parseDouble(lineInfo[3]);
			try {
				omegaPValue = Double.parseDouble(lineInfo[4]);
			} catch(NumberFormatException nfe) {
				omegaPValue = 0.0001;
			}
			
			
		}
		
		public ScaledWindow(String chr) {
			super(chr);
		}
		
		public String toFullString() {
			return toFullString(DEFAULT_FORMATTER);
		}
		
		public String toFullString(DecimalFormat formatter) {
			StringBuilder sb = new StringBuilder(getChromosome().startsWith("chr") ? getChromosome() : "chr" + getChromosome());
			sb.append("\t").append(String.valueOf(getStart()))
			.append("\t").append(String.valueOf(getEnd()))
			.append("\t").append(formatter.format(getOmega()))
			.append("\t").append(formatter.format(getOmegaPValue()))
			.append("\t").append(formatter.format(getLogOdds()))
			.append("\t").append(formatter.format(getLogOddsPValue()));
			
			return sb.toString();
		}
		
		public String toString() {
			StringBuilder sb = new StringBuilder(String.valueOf(getStart()));
			sb.append("\t").append(DEFAULT_FORMATTER.format(getOmega()))
			.append("\t").append(DEFAULT_FORMATTER.format(getTreeLength()))
			.append("\t").append(DEFAULT_FORMATTER.format(getLogOdds()))
			.append("\t").append(DEFAULT_FORMATTER.format(getLogOddsPValue()));
			
			return sb.toString();			
			
		}
		
		public String toMidBaseLogPvalScoredGGraph(DecimalFormat formatter) {
			int midPoint = (getStart() + getEnd())/2;
			StringBuilder sb = new StringBuilder(getChromosome());
			sb.append("\t").append(String.valueOf(midPoint))
			.append("\t").append(formatter.format(-Math.log(getOmegaPValue())));
			
			return sb.toString();			
		}
		
		public void stitchToList(List<ScaledWindow> toStitch) {
			Iterator<ScaledWindow> it = toStitch.iterator();
			double avgOmega = omega;
			double avgpVal  = omegaPValue;
			double avgtreeLength = treeLength;
			
			while(it.hasNext()) {
				ScaledWindow sw = it.next();
				avgOmega += sw.omega;
				avgpVal  += sw.omegaPValue;
				avgtreeLength += sw.treeLength;
				super.stitchTo(sw);
			}
			omega      = avgOmega / (toStitch.size() + 1);
			treeLength = avgtreeLength / (toStitch.size() + 1);
			omegaPValue     = avgpVal / (toStitch.size() + 1);
		}
		
		public double getLogLikelihoodRatio(){return getLogOdds();}
		public int getPosition() {return getStart();}
		
		public double getLogOdds() {
			return logOdds;
		}
		public void setLogOdds(double logOdds) {
			this.logOdds = logOdds;
		}
		public double getLogOddsPValue() {
			return logOddsPValue;
		}
		public void setLogOddsPValue(double logOddsPValue) {
			this.logOddsPValue = logOddsPValue;
		}

		public double getOmegaPValue() {
			return omegaPValue;
		}
		public void setOmegaPValue(double value) {
			omegaPValue = value;
		}
		public double getOmega() {
			return omega;
		}
		public void setOmega(double omega) {
			this.omega = omega;
		}
		public double getTreeLength() {
			return treeLength;
		}
		public void setTreeLength(double treeLength) {
			this.treeLength = treeLength;
		}
		
		
	}


	public void addWindow(ScaledWindow sw) {
		scalings.add(sw);
	}

	public List<ScaledWindow> getScaledWindows() {
		return scalings;
	}

	public void clear() {
		scalings.clear();
		
	}
	
	public void setNumberFormatter (DecimalFormat formatter) {
		this.formatter = formatter;
	}

	public int getEstimationWindow() {
		return estimationWindow;
	}

	public void setEstimationWindow(int estimationWindow) {
		this.estimationWindow = estimationWindow;
	}

	public int getTotalInSources() {
		return totalInSources;
	}

	public void permute() {
		Random r = new Random();
		System.out.print("Permutting ... had " + scalings.size() + " scalings ");
		ArrayList<ScaledWindow> permuttedScalings = new ArrayList<ScaledWindow>(scalings.size());
		while(scalings.size() > 0) {
			int index = r.nextInt(scalings.size());
			ScaledWindow sw = scalings.remove(index);
			sw.setStart(permuttedScalings.size());
			sw.setEnd(sw.getStart() + 1);
			permuttedScalings.add(sw);
		}
		
		scalings = permuttedScalings;
		System.out.println(" ..  ended up with hopefully the same: " + scalings.size());
	}




}
