package broad.pda.feature.genome;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;

import org.apache.log4j.Logger;

import broad.core.annotation.BasicGenomicAnnotation;
import broad.core.annotation.GenomicAnnotation;
import broad.core.sequence.SequenceRegion;
import broad.pda.seq.protection.ReadSimulator2;

/**
 * Represents a Genome Assembly were each directory is in its own directory which includes the AGP file
 * It is a convenient way to access genomic data when it is stored in this standard way.
 * @author mgarber
 *
 */
public  class DirectoryInstalledGenomeAssembly {
	private static final Logger logger = Logger.getLogger(DirectoryInstalledGenomeAssembly.class.getName());
	protected File sequenceDirectory;
	private Chromosome X;
	private Chromosome Y;
	private Chromosome M;
	private List<Chromosome> Un;
	private ArrayList<Chromosome> autosomes;
	private ArrayList<Chromosome> random;
	private ArrayList<Chromosome> altHaplotypes;
	private File repeatDirectory;
	private File annotationDirectory;

	

	
	public DirectoryInstalledGenomeAssembly(File sequenceDirectory) throws Exception {
		if(!sequenceDirectory.isDirectory()) {
			throw new Exception("sequence directory " + sequenceDirectory +" is not a directory");
		}
		this.sequenceDirectory = sequenceDirectory;
		String seqDirPath = sequenceDirectory.getAbsolutePath();
		seqDirPath = seqDirPath.lastIndexOf("/") == seqDirPath.length() -1 
			? seqDirPath.substring(0, seqDirPath.length() -1)
			: seqDirPath;
			
		this.repeatDirectory   = new File(seqDirPath + "_repeatinfo/");
		this.annotationDirectory = new File(seqDirPath + "_annotations/");
		this.autosomes = new ArrayList<Chromosome>();
		this.random = new ArrayList<Chromosome>();
		this.Un = new ArrayList<Chromosome>();
		this.altHaplotypes = new ArrayList<Chromosome>();
		
		String [] seqSubDirs = sequenceDirectory.list();		
		for(int i = 0; i < seqSubDirs.length; i++) {
			File dir = new File(sequenceDirectory.getAbsoluteFile() + "/" + seqSubDirs[i]);
			if(! dir.isDirectory() || dir.getName().contains("tmp")) {
				continue;
			}
			
			if(dir.getName().matches("Un[0-9]+")) {
				continue;
			}
			if(useMergedUn() && dir.getName().startsWith("Un") ) {
				continue;
			}
			//System.out.println("Processing " + dir + " looking for agp " + seqDirPath + "/" + seqSubDirs[i] + "/chr" + seqSubDirs[i] + ".agp");
			String agpFile = seqDirPath + "/" + seqSubDirs[i] + "/chr" + seqSubDirs[i] + ".agp";
			String [] fileComponents = agpFile.split(File.separator);
			String name = fileComponents[fileComponents.length - 1].replace(".agp", "");
			
			Chromosome chr = new Chromosome(name, seqDirPath + "/" + seqSubDirs[i] + "/chr" + seqSubDirs[i] + ".agp");
			if(chr.getSymbol().endsWith("_random") ) {
				random.add((chr));
			} else if("X".equals(seqSubDirs[i])) {
				//System.out.println("Setting X" );
				X = chr;
			} else if ("Y".equals(seqSubDirs[i])){
				//System.out.println("Setting Y");
				Y = chr;
			} else if ("M".equals(seqSubDirs[i])){
				//System.out.println("Setting M");
				M = chr;
			} else if (seqSubDirs[i].contains("Un")) {
				Un.add(chr);
				random.add(chr);
			} else if(seqSubDirs[i].contains("hap")){
				altHaplotypes.add(chr);
			}else {
				//System.out.println("Adding chromosome " + chr.getSymbol() + " to the autosomes");
				autosomes.add(chr);
			}
		}

		Collections.sort(autosomes, new Comparator<Chromosome>() {

			public int compare(Chromosome arg0, Chromosome arg1) {
				return (int) (arg1.length() - arg0.length());
			}
			
		});
	}
	
	public List<Chromosome> getAllNonRandomChromosomes() {
		ArrayList<Chromosome> all = new ArrayList<Chromosome>(autosomes.size() + 3);		
		all.addAll(autosomes);
		if(X != null) {
			all.add(X);
		}
		
		if(Y != null) {
			all.add(Y);
		}
		
		if(M != null) {
			all.add(M);
		}
		
		return all;
	}
	
	public List<Chromosome> getAllChromosomes() {
		List<Chromosome> all = getAllNonRandomChromosomes();
		all.addAll(random);
		
		return all;
	}
	
	public Chromosome getX() { return X;}
	

	
	public void writeSizes(String outFile) throws IOException {
		BufferedWriter bw = new  BufferedWriter(new FileWriter(outFile));
		Iterator<Chromosome> chrIt = getAllChromosomes().iterator();
		while(chrIt.hasNext()) {
			Chromosome c = chrIt.next();
			bw.write("chr" + c.getName() + "\t" + c.length());
			bw.newLine();
		}
		bw.close();
	}
	
	public long getNonRandomTotalSize() {
		long size = 0;
		Iterator<Chromosome> chrIt = getAllNonRandomChromosomes().iterator();
		while(chrIt.hasNext()) {
			Chromosome c = chrIt.next();
			//System.out.println("Adding size of chr" + c.getSymbol() + ": " + c.getSize());
			size += c.length();
		}		
		return size;
	}
	
	public List<Chromosome> getAutosomes() {
		ArrayList<Chromosome> autosomesCopy = new ArrayList<Chromosome>(autosomes.size());
		autosomesCopy.addAll(autosomes);
		return autosomesCopy;
	}


	public int getAutosomeNumber() {
		return autosomes.size();
	}
	
	public Chromosome getChromosome(String chr) {
		String chrSymbol = chr;
		if(chrSymbol != null) {
			chrSymbol = chrSymbol.replace("chr", "");
		}
		Chromosome c = null;
		if(X != null && X.getSymbol().equals(chrSymbol)) {
			c = X;
		}else if((Y != null) && Y.getSymbol().equals(chrSymbol)) {
			c = Y;
		}else if((M != null) && M.getSymbol().equals(chrSymbol)) {
			c = M;
		} 
		
		if(c == null) {
			c = findChromosomeInList(autosomes, chrSymbol);
		} 
		
		if (c == null) {
			c = findChromosomeInList(random, chrSymbol);
		}
		return c;
	}
	
	public List<Chromosome> getUn() {
		return Un;
	}
	
	public long getGenomeSize() {
		long size = 0;
		if(X != null) {
			size += X.length();
		}
		
		if(Y != null) {
			size += Y.length();
		}
		
		Iterator<Chromosome> it = autosomes.iterator();
		while(it.hasNext()) {
			size += it.next().length();
		}
		return size;
	}
	
	public long getAutosomeGenomeSize() {
		long size = 0;

		Iterator<Chromosome> it = autosomes.iterator();
		while(it.hasNext()) {
			size += it.next().length();
		}
		return size;
	}
	
	/**
	 * Uses pseudo random number generation to draw a (non random) chromosome with
	 * chance proportional to the chromosome size
	 * @return Chromosome draw by random.
	 */
	public Chromosome  drawChromosome() {
		Chromosome chosen = null;
		double draw = Math.random();
		long genomeSize = getGenomeSize();
		long currentSize = 0;
		Iterator<Chromosome> it = getAllNonRandomChromosomes().iterator();
		while(it.hasNext()) {
			chosen = it.next();
			currentSize += chosen.length();
			if(currentSize/(double)genomeSize >= draw) {
				break;
			}
		}

		return chosen;
	}
	
	public SequenceRegion drawRandomRegion(int size) {
		Chromosome c = drawChromosome();
		return c.drawRandomRegion(size);
	}
	
	/**
	 * Uses pseudo random number generation to draw a (non random) autosome with
	 * chance proportional to the autosome size
	 * @return Chromosome draw by random.
	 */
	public Chromosome  drawAutosome() {
		Chromosome chosen = null;
		double draw = Math.random();
		long genomeSize = getAutosomeGenomeSize();
		long currentSize = 0;
		Iterator<Chromosome> it = getAutosomes().iterator();
		while(it.hasNext()) {
			chosen = it.next();
			currentSize += chosen.length();
			if(currentSize/(double)genomeSize >= draw) {
				break;
			}
		}

		return chosen;
	}
	
	public SequenceRegion drawRandomAutosomeRegion(int size) {
		Chromosome c = drawAutosome();
		return c.drawRandomRegion(size);
	}
	
	
	public File getSequenceDir() {
		return sequenceDirectory;
	}
	
	/**
	 * to be overrided if repeats are in directories
	 * @return
	 */
	protected boolean areRepeatsInDeepDirHierarchy() {
		return !(new File(repeatDirectory.getAbsolutePath() + "/chr1.fa.out").exists());
	}
	
	protected boolean useMergedUn() {
		return false;
	}
	
	/**
	 * subclasses should override if the repeats contain overlapping segments (i.e. are merged with trf)
	 * @return
	 */
	protected boolean doRepeatListNeedsCleanup() {return true;}
	
	public File getRepeatFile(String chrSymbol) {
		return new File(repeatDirectory.getAbsolutePath() + "/"  + 
				(areRepeatsInDeepDirHierarchy() ? chrSymbol + "/" : "") + 
				"chr"+chrSymbol + ".fa.out");
	}
	

	
	private Chromosome findChromosomeInList(List<Chromosome> chrs, String symbol) {
		Iterator<Chromosome> it = chrs.iterator();
		Chromosome c = null;
		while(c == null && it.hasNext()) {
			Chromosome chr = it.next();
			if(symbol.equals(chr.getSymbol())) {
				c = chr;
			}
		}
		return c;
	}
	
	public File getAnnotationDirectory() {
		return annotationDirectory;
	}
	
	public List<GenomicAnnotation> getReducedRepresentation(int totalSize, List<GenomicAnnotation> toAvoid) {
		long totalGenomicSize = getGenomeSize();
		List<Chromosome> nonRandomChrs = getAutosomes();
		Iterator<Chromosome>  nonRandomChrIt = nonRandomChrs.iterator();
		ArrayList<GenomicAnnotation> chunks = new ArrayList<GenomicAnnotation>(nonRandomChrs.size());
		Chromosome c = null;
		while(nonRandomChrIt.hasNext()) {
			c = nonRandomChrIt.next();
			double chrSizePrct = c.length()/(double)totalGenomicSize;
			int chunkSize = (int) Math.round(totalSize * chrSizePrct);
			BasicGenomicAnnotation chunk = new BasicGenomicAnnotation("chr" + c.getSymbol() + "_chunk");
			chunk.setChromosome(c.getSymbol());
			chunk.setStart(1);
			chunk.setEnd(1 + chunkSize);
			chunks.addAll(c.shuffle(chunk, toAvoid));
			System.out.println("Added chunk of size " + chunk.getLength() + " for chr" + c.getSymbol() + " out of intended total " + totalSize);
		}
		
		return chunks;
	}

	public void delete(SequenceRegion region) {
		Chromosome c = getChromosome(region.getChromosome());
		c.delete(region);
	}

	public void invert(SequenceRegion region) {
		Chromosome c = getChromosome(region.getChromosome());
		c.invert(region);
	}
	
	public void insert(SequenceRegion region, GenomicAnnotation insertionPoint) {
		Chromosome c = getChromosome(insertionPoint.getChromosome());
		c.insert(region, insertionPoint);
	}

	public void translocate(SequenceRegion region1, SequenceRegion region2) {
		Chromosome c1 = getChromosome(region1.getChromosome());
		Chromosome c2 = getChromosome(region2.getChromosome());
		
		c1.delete(region1);
		GenomicAnnotation insertPoint1 = new BasicGenomicAnnotation(region1.getName(), region1.getChromosome(), region1.getStart(), region1.getStart() + 1);
		c1.insert(region2, insertPoint1);
		
		c2.delete(region2);
		GenomicAnnotation insertPoint2 = new BasicGenomicAnnotation(region2.getName(), region2.getChromosome(), region2.getStart(), region2.getStart() + 1);
		c2.insert(region1, insertPoint2);
	}

}
