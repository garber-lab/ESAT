package broad.pda.feature.genome;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Random;

import org.apache.log4j.Logger;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.BasicAnnotation;

import broad.core.annotation.BasicGenomicAnnotation;
import broad.core.annotation.GenomicAnnotation;
import broad.core.annotation.LightweightGenomicAnnotation;
import broad.core.datastructures.IntervalTree;
import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.core.sequence.SequenceRegion;

import broad.pda.datastructures.Alignments;


/**
 * Represents a chromosome in a genome
 * @author engreitz
 *
 */
public class Chromosome {
	static Logger logger = Logger.getLogger(Chromosome.class.getName());
	
	final private String name;
	// private ShortBEDReader maskedRegions; // TODO
	private int length; //TODO
	Annotation centromere;
	List<AgpEntry> gaps;
	List<AgpEntry> clones;
	private AgpEntry shortArm;
	private File sequenceFile;
	
	Sequence sequence;
	
	/**
	 * @param The name of the chromosome    e.g. "chr12"
	 */
	public Chromosome(final String name) {
		this.name = name;
		gaps = new ArrayList<AgpEntry>();
		clones = new ArrayList<AgpEntry>();
	}
	
	/**
	 * @param The name of the chromosome    e.g. "chr12"
	 * @param An AGP file describing the chromosome assembly
	 */
	public Chromosome(String name, String agpFile) throws Exception {
		this(name);
		length = 0;
		File source = new File(agpFile);
		BufferedReader br = new BufferedReader(new FileReader(source));
		String line;
		String agpWithoutExt = agpFile.substring(0,agpFile.lastIndexOf("."));
		setSequenceFile(new File(agpWithoutExt + ".fa"));
		try {
			while((line = br.readLine()) != null) {
				if(line.startsWith("#") || line.trim().length() ==0){
					continue;
				}
				String[] lineSplit = line.split("\t");

				AgpEntry entry = AgpEntryFactory.createEntry(lineSplit);
				switch (entry.getType()) {
				case AgpEntry.CLONE_TYPE :
					clones.add(entry);
					break;
				case AgpEntry.GAP_TYPE :
					gaps.add(entry);
					break;
				case AgpEntry.CENTROMERE_TYPE :
					centromere = entry;
					break;
				case AgpEntry.SHORT_ARM_TYPE : 
					shortArm = entry;
					break;
				}
				length = Math.max(length, entry.getEnd());
			}
		}  finally {
			try {
				br.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}
	
	public boolean isSexChromosome() {
		return isSexChromosome(name);
	}
	
	public Annotation getShortArm() { return this.shortArm;}
	
	public static boolean isSexChromosome(String name) {
		return name.equals("chrX") || name.equals("chrY");
	}
	
	public boolean isAutosome() {
		return isAutosome(name);
	}
	
	public static boolean isAutosome(String name) {
		boolean result = false;
		try {
			Integer.parseInt(getChromosomeNumber(name));
			result = true;
		} catch (NumberFormatException n) {
			// do nothing
		}
		return result;
	}
	
	public String getName() {
		return name;
	}
	
	public String getSymbol() { return getName().replace("chr","");}
	
	public String getChromosomeNumber() {
		return getChromosomeNumber(name);
	}
	
	public static String getChromosomeNumber(String name) {
		String num = name.substring(3);
		return num;
	}
	
	public static int compareNames(String name1, String name2) {
		return new Chromosome(name1).compareTo(new Chromosome(name2));
	}
	
	public int compareTo(Object o) {
		if (!this.getClass().equals(o.getClass())) {
			throw new IllegalArgumentException("Cannot compare across classes.");
		}
		
		Chromosome other = (Chromosome) o;
		if (this.getName().equals(other.getName())) return 0;
	
		boolean isAutosome = this.isAutosome();
		boolean otherAutosome = other.isAutosome();
		String thisNum = this.getChromosomeNumber();
		String otherNum = other.getChromosomeNumber();
		
		int result;
		if (isAutosome && otherAutosome) {
			result = -1;
		} else if (!isAutosome() && otherAutosome) {
			result = 1;
		} else if (isAutosome && otherAutosome) {
			result = new Integer(thisNum).compareTo(new Integer(otherNum));
		} else {
			result = thisNum.compareTo(otherNum);
		}
		return result;
	}
	
	public long getUngappedSize() {
		long totalGaps = shortArm != null ? shortArm.getLength() : 0;
		
		Iterator<AgpEntry> gapIt = gaps.iterator();
		while(gapIt.hasNext()) {
			totalGaps += gapIt.next().getLength();
		}
		
		return length - totalGaps;
	}
	
	public void setSequenceFile(File file) {
		this.sequenceFile = file;
	}
	
	public int length() {return length;}
	
	public List<AgpEntry> gaps() { return gaps;}
	
	public static class AgpEntry extends SequenceRegion{

		static public final int CLONE_TYPE       = 0;
		static public final int CENTROMERE_TYPE  = 1;
		static public final int GAP_TYPE         = 2;
		static public final int CONTIG_TYPE      = 4;
		static public final int OTHER_TYPE		 = 3;
		public static final int SHORT_ARM_TYPE   = 5;
		private boolean reversedOrientation;
		private int type;
		private int number;
		
		public AgpEntry(String parentSequence) {
			super(parentSequence);
		}
		

		protected void setInReverseOrientation(boolean isInreverseOrientation) {
			this.reversedOrientation = isInreverseOrientation;
		}
		public boolean inReversedOrientation() {
			return reversedOrientation;
		}
			
		public void setType(int type) {
			this.type = type;
		}
		
		public int getType() {
			return type;
		}
		
		public int getLength() { return super.getLength() + 1;}
		
		public String toString() {
			StringBuffer buf = new StringBuffer("chr"+getContainingSequenceId());
			
			buf.append("\t")
				.append(getStart())
				.append("\t")
				.append(getEnd())
				.append("\t")
				.append(getNumber())
				.append("\t");
			
			switch (getType()) {
			case AgpEntry.CLONE_TYPE:
				buf.append("F\t");
				break;
			case AgpEntry.CENTROMERE_TYPE:
			case AgpEntry.SHORT_ARM_TYPE:
			case AgpEntry.GAP_TYPE:
				buf.append("N\t");
				break;
			case AgpEntry.CONTIG_TYPE:
				buf.append("W\t");
				break;
			default:
				buf.append("O\t");
				break;
			}
			
			buf.append(getName())
				.append("\t")
				.append("0\t")
				.append(getLength() + 1)
				.append("\t+");
			
			return buf.toString();
		}


		public int getNumber() {
			return number;
		}


		public void setNumber(int number) {
			this.number = number;
		}

	}
	
	public static class AgpEntryFactory {

		public static AgpEntry createEntry(String [] rawInfo) {
			AgpEntry entry = null;
			String parentName = rawInfo[0];
			
			entry = new AgpEntry(parentName);
			entry.setName(rawInfo[5]);
			entry.setInReverseOrientation(rawInfo.length == 9 && "-".equals(rawInfo[8]));	
			entry.setStart(Integer.parseInt(rawInfo[1]));
			entry.setEnd(Integer.parseInt(rawInfo[2]));
			entry.setNumber(Integer.parseInt(rawInfo[3]));
			if(rawInfo[0].startsWith("Un") || rawInfo[0].startsWith("un")) {
				entry.setChromosome("Un");
			} else {
				entry.setChromosome(rawInfo[0].length() < 4 ? rawInfo[0] : rawInfo[0].substring(3)); //try to handle both chrNN and NN notations for chromosomes
				//System.out.println("Entry's chromosome " + entry.getChromosome());
				if(entry.getChromosome().startsWith("0")) {
					entry.setChromosome(entry.getChromosome().substring(1));
				}
			}
			if ("F".equals(rawInfo[4])) {			
				entry.setName(rawInfo[5]);
				entry.setType(AgpEntry.CLONE_TYPE);
			} else if ("centromere".equalsIgnoreCase(rawInfo[6])) {
				entry.setName("centromere");
				entry.setType(AgpEntry.CENTROMERE_TYPE);	
			} else if ("short_arm".equalsIgnoreCase(rawInfo[6])) {
				entry.setName("short_arm");
				entry.setType(AgpEntry.SHORT_ARM_TYPE);	
			} else if ("N".equals(rawInfo[4])) {
				entry.setType(AgpEntry.GAP_TYPE);
			} else if ("W".equalsIgnoreCase(rawInfo[4])) {
				entry.setType(AgpEntry.CONTIG_TYPE);
			} else {
				entry.setName("other");
				entry.setType(AgpEntry.OTHER_TYPE);
				//System.out.print("Can't determine what this AGP entry is type<" + rawInfo[4] +"> name <  rawInfo<"+rawInfo[5]);
			}
			
			//System.out.println("created entry " + entry);
			return entry;
		}


	}
	
	public SequenceRegion drawRandomRegion(int size) {
		GenomicAnnotation region = new BasicGenomicAnnotation("initial",getSymbol(), 1, size + 1);
		List<? extends Annotation> randomized = shuffleWithinValid(region, getValidRegions(null).getBlocks());
		//TODO: Handle the case one more than one regions are returned!
		if(randomized.size() > 1) {
			logger.warn("WARNING: random region was actually randomized to 2 regions");
		}
		
		SequenceRegion toExtract = new SequenceRegion(getName(), randomized.get(0));
		if(sequence != null && sequence.getSequenceBases().length() > 0) {
			getRegion(toExtract);
		}
		return toExtract;
	}
	
	private Annotation getValidRegions(List<? extends Annotation> toAvoid) {
		if(toAvoid == null) {
			toAvoid = new ArrayList<Annotation>(); //initialize if null to avoid null pointer exceptions....
		}
		ArrayList<Annotation> toAvoidInChr = new ArrayList<Annotation>();
			
		if(centromere != null)
			toAvoidInChr.add(centromere);
		if(shortArm != null) 
			toAvoidInChr.add(shortArm);
		
		toAvoidInChr.addAll(gaps());
		
		Iterator<? extends Annotation> avoidIt = toAvoid.iterator();
		while(avoidIt.hasNext()) {
			Annotation avoid = avoidIt.next();
			if(getName().equals(avoid.getChr())) {
				toAvoidInChr.add(avoid);
				//System.out.println("\tAdded " + avoid + " to avoid in chr " + getSymbol());
			}
		}
		
		Annotation allOfMe = new BasicAnnotation(getName(), 1, length());
		
		Annotation valid = allOfMe.minus(toAvoidInChr);
		return valid;
	}
	
	public List<? extends GenomicAnnotation> shuffle(GenomicAnnotation toShuffle, List<? extends GenomicAnnotation> toAvoid) {
		ArrayList<GenomicAnnotation> oneMemberList = new ArrayList<GenomicAnnotation>(1);
		oneMemberList.add(toShuffle);
		return shuffle(oneMemberList, toAvoid);
	}
	
	public List<? extends Annotation> shuffleWithinValid(Annotation toShuffle, List<? extends Annotation> valid) {
		ArrayList<Annotation> oneMemberList = new ArrayList<Annotation>(1);
		oneMemberList.add(toShuffle);
		return shuffleWithinValid(oneMemberList, valid);
	}
	
	/**
	 * Randomizes (using a quasi uniform distribution) a list of genomic annotations
	 * @param toShuffle List of genomic annotations to shuffle
	 * @param toAvoid regions to avoid placing a shuffled annotation.
	 * @return
	 */
	public List<? extends GenomicAnnotation> shuffle(List<? extends GenomicAnnotation> toShuffle, List<? extends GenomicAnnotation> toAvoid) {
		Annotation valid = getValidRegions(toAvoid);
		//System.out.println("Regions to avoid: " + (toAvoid != null ? toAvoid.size() : 0) + " so we got valid regions " + (valid != null ? valid.size() : 0));
		return shuffleWithinValid(toShuffle,  valid.getBlocks());
	}

	public List<? extends GenomicAnnotation> shuffleWithinValid(List<? extends Annotation> toShuffle, List<? extends Annotation> valid) {
		ArrayList<GenomicAnnotation> randomized = new ArrayList<GenomicAnnotation>();
		Random randomizer = new Random(Math.round(Math.random()*1000000) );
		// The next hack attempts to simulate the trivial fact that large regions should be more likely
		// to be drawn by adding as many copies of it as its percent of the total length covered by
		// valid regions. The caveat is that we are rounding the percent thus we used 1000 rather than 100 
		// to account for small but not tiny regions... This may slow things down significantly when randomizing
		// al large list.
		//System.out.println("Total valid regions length " + totalValidLength);
		List<Annotation> normalizedValid = normalizeRegionsListBySize(valid);
		
		Iterator<? extends Annotation> toRandomizeIt = toShuffle.iterator();
		
		int normalizedValidNumber = normalizedValid.size();
		//System.out.println("Valid num: " + validNumber + " validNumber order of magnitude " + orderOfMagnitudeValid + " normlizedValidSize " + normalizingConstant + " total normalized valid " + normalizedValidNumber);
		while(toRandomizeIt.hasNext()) {
			Annotation original = toRandomizeIt.next();
			int randomValidRegion = randomizer.nextInt(normalizedValidNumber);
			Annotation validReg = normalizedValid.get(randomValidRegion);
			int randomizedStart = randomizer.nextInt(validReg.length()) + validReg.getStart();
			GenomicAnnotation randomizedOriginal = new BasicGenomicAnnotation(original);
			randomizedOriginal.setName(randomizedOriginal.getName() + "_randomized");
			randomized.add(randomizedOriginal);
			randomizedOriginal.setStart(randomizedStart);
			randomizedOriginal.setEnd(Math.min(randomizedStart + original.length(), validReg.getEnd()));
			GenomicAnnotation prior = randomizedOriginal;
			
			// Initially we got a random region from the NORMALIZED vector, from now on we work on the small valid region vector
			randomValidRegion = valid.indexOf(validReg);
			int i = 2;
			int left = original.length() - prior.length();
			while(left > 0) {
				//System.out.print("\tPrior valid region " + randomValidRegion);
				randomValidRegion = (randomValidRegion + 1) % valid.size();
				//System.out.print(" new valid region " + randomValidRegion);
				Annotation nextValidRegion = valid.get(randomValidRegion);
				//System.out.println(" region: " + nextValidRegion);
				GenomicAnnotation nextRandomizedChunk = new BasicGenomicAnnotation(original);
				nextRandomizedChunk.setName(original.getName() + "_randomized_" + i++);
				nextRandomizedChunk.setStart(nextValidRegion.getStart());
				nextRandomizedChunk.setEnd(Math.min(nextValidRegion.getStart() + left, nextValidRegion.getEnd()));
				randomized.add(nextRandomizedChunk);
				left = left - nextRandomizedChunk.length();
			}
		}
		return randomized;
	}
	
	private List<Annotation> normalizeRegionsListBySize(List<? extends Annotation> valid) {
		Iterator<? extends Annotation> validIt = valid.iterator();
		int totalValidLength = 0;
		while(validIt.hasNext()) {
			totalValidLength += validIt.next().length();
		}
		
		int valids = valid.size();	
		validIt = valid.iterator();
		double orderOfMagnitudeValid = Math.ceil(Math.log10(valids));
		int normalizingConstant =(int) Math.pow(10, (orderOfMagnitudeValid + 2));

		List<Annotation> normalizedValid = new ArrayList<Annotation>(normalizingConstant);
		while(validIt.hasNext()) { 
			Annotation reg = validIt.next();
			double roundedPctCoveratge = Math.round(reg.length()/(float)totalValidLength * normalizingConstant);
			//System.out.println("\tRegion " + reg + " is " + roundedPctCoveratge + " of total ");
			for(int i = 0; i < roundedPctCoveratge; i++) {
				//System.out.println("\t\tAdding region " + reg.getName() + ", total newValid size: " + newNormalizedValid.size());
				normalizedValid.add(reg);
			}
		}
		return normalizedValid;
	}

	
	public void getRegion(SequenceRegion region) {
		sequence.getRegion(region);
	}
	
	public void getRegion(SequenceRegion region, boolean softmask) {
		sequence.getRegion(region, softmask);
	}
	
	public void getRegion(SequenceRegion region, boolean softmask, Map<String, IntervalTree<Alignments>> okRepeats) {
		sequence.getRegion(region, softmask, okRepeats);
	}
	
	public List<SequenceRegion> getRegions(List<? extends GenomicAnnotation> annotations, int padding) throws IOException {
		FastaSequenceIO fsIO = new FastaSequenceIO(sequenceFile.getAbsolutePath());
		ArrayList<SequenceRegion> regs = new ArrayList<SequenceRegion>(annotations.size());
		Iterator<? extends GenomicAnnotation> it = annotations.iterator();
		while(it.hasNext()) {
			GenomicAnnotation annot = it.next();
			SequenceRegion reg = new SequenceRegion("chr"+getSymbol());
			reg.setRegionStart(annot.getStart() - padding);
			reg.setRegionEnd(annot.getEnd() + padding);
			reg.setId(annot.getName());
			regs.add(reg);
		}
		fsIO.extractRegions(regs);
		return regs;
	}
	
	public void extractRegions(List<? extends SequenceRegion> regions) throws IOException {
		if(sequence == null ){//|| sequence.getSequenceBases() == null || sequence.getSequenceBases().length() == 0) {
			loadSequence();
		}
		sequence.getRegions(regions);
	}
	
	public void extractRegion(SequenceRegion region) throws IOException {
		ArrayList<SequenceRegion> regionList = new ArrayList<SequenceRegion>(1);
		regionList.add(region);
		extractRegions(regionList);
	}
	
	public void loadSequence() throws IOException {
		FastaSequenceIO fsio = new FastaSequenceIO(sequenceFile);
		//System.err.println("  Extracting ... expted size " + size);
		fsio.extractRecordsWithIDLike(getSymbol(), false, length);
		sequence = fsio.loadAll().get(0);
		//size = sequence.getSequenceBases().length();
	}
	
	public void unloadSequence() { sequence.unloadAllSequences(); System.gc();} 
	
	public Sequence getSequence () { return sequence; }


	public void insert(SequenceRegion region, LightweightGenomicAnnotation insertPoint) {
		StringBuilder newSequence = new StringBuilder(sequence.getSequenceBases().substring(0,insertPoint.getStart()));		
		newSequence.append(region.getSequence().getSequenceBases());
		newSequence.append(sequence.getSequenceBases().substring(insertPoint.getStart()));
		sequence.setSequenceBases(newSequence.toString());
		
		Iterator<AgpEntry> gapIt = gaps.iterator();
		while(gapIt.hasNext()) {
			AgpEntry gap = gapIt.next();
			if(gap.getStart() >= insertPoint.getStart()) {
				gap.setStart(gap.getStart() + region.getLength());
				gap.setEnd(gap.getEnd() + region.getLength());
			}
		}
		
		if(centromere.getStart() >= insertPoint.getStart()) {
			centromere.setStart(centromere.getStart() + region.getLength());
			centromere.setEnd(centromere.getEnd() + region.getLength());
		}
	}
	
	public void delete(SequenceRegion region) {
		StringBuilder newSequence = new StringBuilder(sequence.getSequenceBases().substring(0, region.getStart()));
		newSequence.append(sequence.getSequenceBases().substring(region.getEnd() + 1));
		sequence.setSequenceBases(newSequence.toString());
		
		Iterator<AgpEntry> gapIt = gaps.iterator();
		while(gapIt.hasNext()) {
			AgpEntry gap = gapIt.next();
			if(gap.getStart() >= region.getStart()) {
				gap.setStart(gap.getStart() - region.getLength());
				gap.setEnd(gap.getEnd() - region.getLength());
			}
		}
		
		if(centromere.getStart() >= region.getStart()) {
			centromere.setStart(centromere.getStart() - region.getLength());
			centromere.setEnd(centromere.getEnd() - region.getLength());
		}
	}
	
	public void invert(SequenceRegion region) {
		StringBuilder newSequence = new StringBuilder(sequence.getSequenceBases().substring(0, region.getStart()));
		region.reverse();
		newSequence.append(region.getSequenceBases());
		newSequence.append(sequence.getSequenceBases().substring(region.getEnd() + 1));
		sequence.setSequenceBases(newSequence.toString());
	}
}
