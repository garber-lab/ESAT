package broad.pda.seq.protection;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Random;
import java.util.TreeSet;

import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.log4j.Logger;

import broad.core.datastructures.Pair;
import broad.core.sequence.Sequence;
import broad.pda.annotation.BEDFileParser;
import broad.pda.feature.genome.Chromosome;
import broad.pda.feature.genome.DirectoryInstalledGenomeAssembly;
import broad.pda.seq.alignment.AlignmentUtils;
import net.sf.samtools.BAMIndex;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMProgramRecord;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.SAMFileHeader.SortOrder;
import nextgen.core.alignment.Alignment;
import nextgen.core.alignment.SingleEndAlignment;
import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.BasicAnnotation;
import nextgen.core.annotation.Gene;
import nextgen.core.feature.GeneWindow;
import net.sf.picard.sam.BuildBamIndex;

public class ReadSimulator2 {

	private static final Logger logger = Logger.getLogger(ReadSimulator2.class.getName());
	private static final int MIN_DIST_BETWEEN_BOUND = 50; 
	private static final int BOUND_REGION_SIZE = 10; //TODO This should be variable per site
	private static final int readLength=26; //TODO Make this variable
	private double insertSizeMean= 30; //TODO Make this variable
	private double insertSizeSD=10; //TODO Make this variable
	private DirectoryInstalledGenomeAssembly genome;
	private SAMFileHeader genomeHeader;
	
	
	public ReadSimulator2(Gene gene, int numProteins, int coverage, double enrichment, double secondReadBias, String genomeDir, String save) throws Exception{
		genome = new DirectoryInstalledGenomeAssembly(new File(genomeDir));
		this.genomeHeader=this.getGenomeHeader(genome);
		
		//Step 1: Assign binding sites
		List<Annotation> boundRegions=generateBoundRegions(numProteins, gene);
		
		//Step 2: Generate reads from bound regions
		Pair<File> fastqFiles=generateReads(gene, boundRegions, coverage, enrichment, secondReadBias, save);
		
		//Write to file
		List<Annotation> boundRegionsGenome=write(save+".boundRegions.bed", boundRegions, gene);
		
		generateFragments(gene, boundRegionsGenome, coverage, enrichment, save);
		
		//Tophat Align the reads
		tophat(fastqFiles, save);
	}
	
	private void tophat(Pair<File> fastqFiles, String save) throws IOException, InterruptedException {
		String saveDir=fastqFiles.getValue1().getName().split("_")[0];
		String cmd="/seq/lincRNA/Pam/Software/tophat-2.0.8.Linux_x86_64/tophat2 --read-mismatches 2 --read-gap-length 2 --read-edit-dist 2 --mate-inner-dist 10 --mate-std-dev 20 --min-anchor-length 8 --splice-mismatches 0 --min-intron-length 70 --max-intron-length 500000 --max-insertion-length 3 --max-deletion-length 3 --solexa-quals --num-threads 1 --max-multihits 15 --no-discordant --no-mixed --library-type fr-firststrand --segment-mismatches 2 --segment-length 25 --min-segment-intron 50 --max-segment-intron 500000 --b2-sensitive --b2-N 0 --b2-L 20 --b2-i S,1,1.25 --b2-n-ceil L,0,0.15 --b2-gbar 4 --b2-mp 6,2 --b2-np 1 --b2-rdg 5,3 --b2-rfg 5,3 --b2-score-min L,-0.6,-0.6 --b2-D 15 --b2-R 2 --GTF /seq/lincRNA/data/mm9/scriptureIndex/RefSeqUCSCScripture.gff --no-novel-juncs --transcriptome-max-hits 1 --prefilter-multihits --transcriptome-index /seq/lincRNA/data/mm9/bowtie2/RefSeqUCSCScripture ";
		cmd+="-o "+saveDir+" /seq/lincRNA/data/mm9/bowtie2/dowloaded/mm9 "+fastqFiles.getValue1().getAbsolutePath()+" "+fastqFiles.getValue2().getAbsolutePath();
		
		logger.info(cmd);
		Runtime run=Runtime.getRuntime();
		Process p=run.exec(cmd);
		p.waitFor();
				
		//move
		cmd="mv "+saveDir+"/accepted_hits.bam "+ save+".tophat.bam";
		logger.info(cmd);
		p=run.exec(cmd);
		p.waitFor();
		
		//index
		cmd="java -jar /seq/mgarber/tools/picard-tools-1.73/BuildBamIndex.jar I="+save+".tophat.bam";
		logger.info(cmd);
		p=run.exec(cmd);
		p.waitFor();
	}

	private SAMFileHeader getTranscriptomeHeader(Gene gene) {
		// Setting header file for transcriptome alignment
		SAMFileHeader transcriptomeHeader = new SAMFileHeader();
		transcriptomeHeader.addProgramRecord(new SAMProgramRecord("Scripture-simulator"));
		transcriptomeHeader.addComment("Simulated reads");
		transcriptomeHeader.setSortOrder(SortOrder.coordinate);
		transcriptomeHeader.addSequence(new SAMSequenceRecord(gene.getName(), gene.length()));
		
		return transcriptomeHeader;
	}

	private void generateFragments(Gene gene, List<Annotation> boundRegions, int coverage, double enrichment, String save){
		//Take transcript and generate N random fragments across transcript
		Collection<GeneWindow> fragments=fragment(gene, new Double(this.insertSizeMean).intValue());
		//Break fragments into two bins: (i) fragments overlapping binding site, (ii) fragments not overlapping binding site
		Pair<List<GeneWindow>> bins=binFragments(fragments, boundRegions);
		//estimate frequency of bound to total fragments (p)
		double p=(double)bins.getValue1().size()/(double)(bins.getValue1().size()+bins.getValue2().size());
		logger.info("prob(bound)="+p+ " num bound="+bins.getValue1().size()+" num unbound="+bins.getValue2().size());
		//Sample fragments from bin (i) with frequency ep and from bin (ii) with frequency 1-ep
		double probBound=Math.min(1,  enrichment*p);
		logger.info("enrichment prob="+probBound);
		double probUnbound=1-probBound;
		
		int numBoundReads=new Double(probBound*coverage).intValue();
		int numUnboundReads=new Double(probUnbound*coverage).intValue();
		
		Collection<GeneWindow> boundReads=getReads(bins.getValue1(), numBoundReads);
		Collection<GeneWindow> unboundReads=getReads(bins.getValue2(), numUnboundReads);
		
		try {
			write(save+".fragments.bed", boundReads, unboundReads);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	private void write(String save, Collection<GeneWindow> boundReads, Collection<GeneWindow> unboundReads) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(GeneWindow read: boundReads){
			writer.write(read.toBED()+"\n");
		}
		
		for(GeneWindow read: unboundReads){
			writer.write(read.toBED()+"\n");
		}
		
		writer.close();
	}

	private Collection<GeneWindow> getReads(List<GeneWindow> fragments, int numReads) {
		Collection<GeneWindow> rtrn=new ArrayList<GeneWindow>();
		
		for(int i=0; i<numReads; i++){
			int index= new Double(Math.random()*fragments.size()).intValue();
			rtrn.add(fragments.get(index));
		}
		
		return rtrn;
	}

	private Pair<List<GeneWindow>> binFragments(Collection<GeneWindow> fragments, List<Annotation> boundRegions) {
		List<GeneWindow> bound=new ArrayList<GeneWindow>();
		List<GeneWindow> unbound=new ArrayList<GeneWindow>();
		
		for(GeneWindow fragment: fragments){
			if(fragment.overlaps(boundRegions)){bound.add(fragment);}
			else{unbound.add(fragment);}
		}
		
		Pair<List<GeneWindow>> rtrn=new Pair<List<GeneWindow>>(bound, unbound);
		return rtrn;
	}

	private boolean overlaps(GeneWindow fragment, Annotation proteinRegion) {
		if(!fragment.getChr().equalsIgnoreCase(proteinRegion.getChr())){return false;}
		
		for(Annotation block: fragment.getBlocks()){
			for(Annotation block2: proteinRegion.getBlocks()){
				if(overlap(block, block2)){return true;}
			}
		}
		return false;
	}

	private boolean overlap(Annotation block1, Annotation block2) {
		if(block1.getStart()>=block2.getStart() && block1.getStart()<=block2.getEnd()){return true;}
		if(block2.getEnd()>=block1.getStart() && block2.getEnd()<=block1.getEnd()){return true;}
		
		if(block2.getStart()>=block1.getStart() && block2.getStart()<=block1.getEnd()){return true;}
		if(block1.getEnd()>=block2.getStart() && block1.getEnd()<=block2.getEnd()){return true;}
		
		return false;
	}

	private Collection<GeneWindow> fragment(Gene gene, int size) {
		//Enumerate all kmers in gene
		Collection<GeneWindow> fragments=gene.getWindows(size);
		return fragments;
	}

	private Pair<File> generateReads(Gene gene, List<Annotation> boundRegions, int coverage,	double enrichment, double secondReadBias, String save) throws IOException, CloneNotSupportedException {
		File fastq1=new File(save+"_1.fq");
		File fastq2=new File(save+"_2.fq");
		Pair<File> rtrn=new Pair<File>(fastq1, fastq2);
		
		// Creating alignment writers
		FileWriter fastqWriter1=new FileWriter(fastq1);
		FileWriter fastqWriter2=new FileWriter(fastq2);
		
		SAMFileWriterFactory alnWriterFactory = new SAMFileWriterFactory();
		File genomeBamFile=new File(save+".bam");
		SAMFileWriter genomeAlignmentWriter = alnWriterFactory.makeSAMOrBAMWriter(genomeHeader, false, genomeBamFile);
		
		long totalReadsForTranscript = coverage; 
		double pctSecondReadAtCrosslink=secondReadBias;
		Chromosome chr = genome.getChromosome(gene.getChrNum());
		chr.loadSequence();
		gene.setSequenceFromChromosome(chr.getSequence());

		// We want to get a starting position from the 3 bimodal
		// distributions:
		// The probability that the inserts comes from a crosslinked
		// segment
		// The probability that when it does come its first read is from
		// the crosslinked site
		// If p = probability of background, and the enrichment is e,
		// then the probability of crosslink is e*p
		// therefore e+ep =1 and p = (1/(1+e))
		// double pctReadsFromCrosslink = 1-(1/(1+enrichmentAtBound));

		// MG: The calculations of probabilities were not what I think
		// the data actually looks like, rather:
		// If the background is uniform and bound regions contribute e
		// times as many reads then the total frequency is denom=
		// 1*(unbound positions) + e*(boundPositions)
		// Then probability of sampling from background is
		// p(unboundPosition)=(1/denom)
		// and the probability of sampling from the protein bound is
		// p(boundPosition)=(e/denom)
		
		
		int numBound = boundRegions.size();
		int numUnbound = gene.length() - numBound;
		double denom = numUnbound + (enrichment * numBound);

		// make quality string
		byte[] quality = new byte[readLength];
		for (int i = 0; i < readLength; i++) {
			quality[i] = 'I';
		}
		
		logger.info("numBound=" + numBound + " numUnbound="+ numUnbound);

		double pctReadsFromCrosslink = enrichment / denom;

		logger.info("p(insert from bound)=" + pctReadsFromCrosslink);

		// This is as close as the first read may get to the end of the gene
		int firstReadDistToGeneEnd = (int) Math.round(readLength + insertSizeMean + insertSizeSD);
		
		RandomDataGenerator rdg = new RandomDataGenerator();
		Random r = new Random();
		
		// OK lets begin generating reads
		for (int j = 0; j < totalReadsForTranscript; j++) {

			int pos = 0;
			double draw = rdg.nextUniform(0, 1);
			if (boundRegions.size() == 0 || draw > pctReadsFromCrosslink) { // If
				// no  crosslink OR draw from non crosslinked segments
				// while(true) {
				pos = r.nextInt(gene.length() - firstReadDistToGeneEnd);
				// if(distToBoundProt(pos, boundMidpoints) > insertSizeMean) {
				// break; //TODO: streamline this, it is very inneficient
				// }
				// }
			} else {

				int protein = boundRegions.size() == 1 ? 0 : rdg.nextInt(0,	boundRegions.size() - 1);
				boolean secondAtCrosslink = pctSecondReadAtCrosslink > rdg.nextUniform(0, 1);
				if (secondAtCrosslink) {
					pos = (int) boundRegions.get(protein).getMidpoint() + BOUND_REGION_SIZE / 2;
				} else {
					int insertSize = 0;
					while (insertSize <= 0) {
						insertSize = (int) Math.round(rdg.nextGaussian(insertSizeMean, insertSizeSD));
					}
					int lower = Math.max(0, (int) boundRegions.get(protein).getMidpoint() - insertSize);
					int upper = Math.min(gene.length() - readLength,	(int) boundRegions.get(protein).getMidpoint()	+ insertSize);
					logger.debug("insertSize " + insertSize + " lower "	+ lower + " upper " + upper);
					pos = rdg.nextInt(lower, upper);
				}

			}
			
			String name="@SCR_" + j + "_"	+ System.nanoTime();
			//Write fastq
			//Write SAM file
			Pair<SAMRecord> genomeAlignment=getSamRecord(gene, pos, name, quality, rdg);
			if(genomeAlignment!=null){
				genomeAlignmentWriter.addAlignment(genomeAlignment.getValue1());
				genomeAlignmentWriter.addAlignment(genomeAlignment.getValue2());
				fastqWriter1.write(toFastq(genomeAlignment.getValue1())+"\n");
				fastqWriter2.write(toFastq(genomeAlignment.getValue2())+"\n");
			}
			/*Pair<SAMRecord> genomeAlignment=getSamRecord(gene, pos, name, quality, rdg);
			if(genomeAlignment!=null){
				genomeAlignmentWriter.addAlignment(genomeAlignment.getValue1());
				genomeAlignmentWriter.addAlignment(genomeAlignment.getValue2());
			}*/
		}
		
		fastqWriter1.close();
		fastqWriter2.close();
		genomeAlignmentWriter.close();
		createIndex(genomeBamFile);
		
		return rtrn;
	}
	
	private String toFastq(SAMRecord value1) {
		String seq=value1.getReadString();
		if(value1.getReadNegativeStrandFlag()){
			seq=Sequence.reverseSequence(seq);
		}
		String name=value1.getReadName();
		String quality=value1.getBaseQualityString();
		return name+"\n"+seq+"\n"+"+\n"+quality;
	}

	private void createIndex(File file){
		File genomeBamIfxFile = new File(file.getAbsolutePath()+ BAMIndex.BAMIndexSuffix);
		if (genomeBamIfxFile.exists()) {
			genomeBamIfxFile.delete();
		}
		SAMFileReader reader = new SAMFileReader(file);
		BuildBamIndex.createIndex(reader, genomeBamIfxFile);
	}

	private SAMFileHeader getGenomeHeader(DirectoryInstalledGenomeAssembly org){
		// Setting up genome alignment header file
		SAMFileHeader genomeHeader = new SAMFileHeader();
		genomeHeader.addProgramRecord(new SAMProgramRecord("Scripture-simulator"));
		genomeHeader.addComment("Simulated reads");
		genomeHeader.setSortOrder(SortOrder.coordinate);
		List<Chromosome> chrs = org.getAllNonRandomChromosomes();
		for (Chromosome chr : chrs) {
			genomeHeader.addSequence(new SAMSequenceRecord(chr.getName(), chr.length()));
		}
		return genomeHeader;
	}

	private Pair<SAMRecord> getSamRecord(Gene gene, int pos, String readName, byte[] quality, RandomDataGenerator rdg) throws CloneNotSupportedException {
		SAMRecord genomeAln=new SAMRecord(genomeHeader);
		Gene read1=gene.trimGene(pos, pos+readLength);
		
		Chromosome chr = genome.getChromosome(read1.getChrNum());
		read1.setSequenceFromChromosome(chr.getSequence());
		//logger.info(readName+" "+read1.toUCSC()+" "+read1.getSequence());
		
		String seq=read1.getSequence();
		if(read1.isNegativeStrand()){
			seq=Sequence.reverseSequence(seq);
		}
		genomeAln.setReadName(readName);
		genomeAln.setReadString(seq);
		genomeAln.setBaseQualities(quality);
		genomeAln.setAlignmentStart(read1.getSAMStart());
		genomeAln.setReferenceName(read1.getChr());
		genomeAln.setCigarString(read1.toCigar());
		genomeAln.setMappingQuality(255);
		genomeAln.setReadNegativeStrandFlag(gene.isNegativeStrand());
		
		//logger.info(read1.toUCSC()+" "+read1.getSequence());
		
		//SAMRecord transcriptAln = new SAMRecord(transcriptomeHeader);
		//String cigar = readLength + "M";
		//transcriptAln.setReadName("@SCR_" + readNum + "_"	+ System.nanoTime());
		//transcriptAln.setReadString(gene.getSequence().substring(pos, pos + readLength));
		//transcriptAln.setBaseQualities(quality);
		//transcriptAln.setAlignmentStart(pos + 1);
		//transcriptAln.setReferenceName(gene.getName());
		//transcriptAln.setCigarString(cigar);
		//transcriptAln.setMappingQuality(255);

		int insertSize = 0;
		while (insertSize <= 0) {
			insertSize = (int) Math.round(rdg.nextGaussian(insertSizeMean, insertSizeSD));
		}

		int secondReadStart = pos + insertSize;
		int secondReadEnd=secondReadStart+readLength;
		
		
		if (secondReadEnd < gene.length()) {
			Pair<SAMRecord> rtrn;
			genomeAln.setReadPairedFlag(true);
			genomeAln.setProperPairFlag(true);
			genomeAln.setInferredInsertSize(insertSize);
			
			
			Gene read2=gene.trimGene(secondReadStart, secondReadEnd);
			read2.setSequenceFromChromosome(chr.getSequence());
			SAMRecord pair=(SAMRecord) genomeAln.clone();
			
			//Set pair info
			seq=read2.getSequence();
			if(gene.isNegativeStrand()){seq=Sequence.reverseSequence(seq);}
			pair.setReadString(seq);
			pair.setCigarString(read2.toCigar());
			pair.setAlignmentStart(read2.getSAMStart());
			pair.setProperPairFlag(true);
			pair.setInferredInsertSize(insertSize); //TODO Should this be -insert?
			
			//Set mate info
			genomeAln.setMateAlignmentStart(pair.getAlignmentStart());
			pair.setMateAlignmentStart(genomeAln.getAlignmentStart());
			genomeAln.setMateReferenceName(pair.getReferenceName());
			pair.setMateReferenceName(genomeAln.getReferenceName());
			
			//If gene is positive strand then:
			//R1 will be right of R2
			//R1 strand will be (-)
			//R2 strand will be (+)
			if(!gene.isNegativeStrand()){
				//Note: genomeAln is defined to be the left most read so this must be R2
				genomeAln.setFirstOfPairFlag(false);
				genomeAln.setSecondOfPairFlag(true);
				
				
				genomeAln.setReadNegativeStrandFlag(false);
				genomeAln.setMateNegativeStrandFlag(true);
				
				//Note: pair is R1
				pair.setFirstOfPairFlag(true);
				pair.setSecondOfPairFlag(false);
				
				pair.setReadNegativeStrandFlag(true);
				pair.setMateNegativeStrandFlag(false);
				rtrn=new Pair<SAMRecord>(pair, genomeAln);
			}
			
			//If gene is negative strand then:
			//R1 will be left of R2
			//R1 strand will be (+)
			//R2 strand will be (-)
			else{
				//Note: genomeAln is defined to be the left most read so this must be R1
				genomeAln.setFirstOfPairFlag(true);
				genomeAln.setSecondOfPairFlag(false);
				
				
				genomeAln.setReadNegativeStrandFlag(false);
				genomeAln.setMateNegativeStrandFlag(true);
				
				//Note: pair is R2
				pair.setFirstOfPairFlag(false);
				pair.setSecondOfPairFlag(true);
				
				pair.setReadNegativeStrandFlag(true);
				pair.setMateNegativeStrandFlag(false);
				rtrn=new Pair<SAMRecord>(genomeAln, pair);
			}
			
			return rtrn;
		}
		
		return null;
	}

	public static net.sf.samtools.SAMRecord mapSamRecord(net.sf.samtools.SAMRecord sam, Gene gene, SAMFileHeader header) {
		net.sf.samtools.SAMRecord  mappedSam = null;

		Alignment alignment = new SingleEndAlignment(sam);

		if(!alignment.getSpliceConnections().isEmpty() ) {
			throw new IllegalArgumentException("Transcript alignment " + sam.toString() + " contained an intron");
		}

		Gene samGene = AlignmentUtils.SAMFormatFullBED(sam);

		if(gene!=null && header.getSequence(gene.getChr())!= null){
			//System.err.print("Reference " + sam.getReferenceName() ); 
			//System.err.println(" length " + reference.getSequenceLength() + ", difference gene, ref "+(gene.getSize()-reference.getSequenceLength()));
			boolean flip=false;
			int relativeStart=samGene.getStart();
			int relativeEnd=samGene.getEnd();
			/*if(gene.isNegativeStrand()){
				//relativeStart=(gene.getSize())-(gene.getEnd());
				//relativeEnd=(gene.getSize())-(gene.getStart()); //+1
				flip=true;
			}*/
			//Gene geneSAM = gene.copy();

			
			//geneSAM.transcriptToGenomicPosition(relativeStart, relativeEnd); //put back end-1
			
			Gene geneSAM=gene.trimGene(relativeStart, relativeEnd);
			geneSAM.setName(samGene.getName());
			geneSAM.setSequence(samGene.getSequence());
			int seqIdx = header.getSequenceIndex(gene.getChr());

			logger.debug("gene chr: " + geneSAM.getChr() + " header idx: " + seqIdx);
			try {
				mappedSam = (net.sf.samtools.SAMRecord) sam.clone();
				mappedSam.setHeader(header);
				//geneSAM.updatePicardSAM(mappedSam, flip, seqIdx);
				String sequence=mappedSam.getReadString();
				//String quality=tokens[10];
				boolean reversedOrientation=mappedSam.getReadNegativeStrandFlag();
				if(flip){
					mappedSam.setReadString(Sequence.reverseSequence(sequence));
					mappedSam.setReadNegativeStrandFlag(!reversedOrientation);
				}

				mappedSam.setAlignmentStart(geneSAM.getStart()+1);
				//picardSAM.setAlignmentEnd(getEnd()+1);
				mappedSam.setCigarString(geneSAM.toCigar());
				mappedSam.setReferenceName(geneSAM.getChr());
				mappedSam.setReferenceIndex(seqIdx);
			} catch (CloneNotSupportedException e) {
				// TODO Handle this error!
				e.printStackTrace(System.err);
			}
		}

		else { logger.error("Could not map sam record to genome: gene " + gene + " header has chr? " + header.getSequence(gene.getChr()) + " record exons: " + samGene.getNumExons());}
		return mappedSam;
	}
		
		
	private List<Annotation> write(String save, List<Annotation> boundRegions, Gene gene) throws IOException {
		FileWriter writer=new FileWriter(save);
		List<Annotation> rtrn=new ArrayList<Annotation>();
		
		for(Annotation region: boundRegions){
			Gene toGenome = gene.trimGene(region.getStart(), region.getEnd());
			rtrn.add(toGenome);
			writer.write(toGenome+"\n");
		}
		
		writer.close();
		return rtrn;
	}

	/**
	 * Generate the locations of the binding sites
	 * @param numProteins the number of protein binding sites to simulate
	 * @param gene the gene to simulate on
	 * @return a collection of protein bound regions
	 */
	private List<Annotation> generateBoundRegions(int numProteins, Gene gene) {
		Random r = new Random();
		int numberOfBindings = numProteins;
		Gene annotation = gene;
		int transcriptLength = annotation.getSize();

		// Now we need to place the proteins on the transcript, we will drop
		// proteins if the do not fit.
		int maxProteinsFitting = (int) Math.floor((transcriptLength - 2 * MIN_DIST_BETWEEN_BOUND)	/ (double) MIN_DIST_BETWEEN_BOUND) - 1;
		if (numberOfBindings > maxProteinsFitting) {
			logger.info("Requested number of proteins for transcript "	+ annotation.getName() + " is too large, re-setting to " + maxProteinsFitting);
			numberOfBindings = maxProteinsFitting;
		}

		List<Annotation> bindings = new ArrayList<Annotation>(numberOfBindings);
		for (int k = 0; k < numberOfBindings; k++) {
			int iterations = 0;
			while (true) {
				int pos = r.nextInt(transcriptLength - MIN_DIST_BETWEEN_BOUND);
				if (pos > MIN_DIST_BETWEEN_BOUND && didNotLandNearExistingProt(bindings, pos)) {
					int start = pos - BOUND_REGION_SIZE / 2;
					int end = pos + BOUND_REGION_SIZE / 2;
					Annotation binding = new BasicAnnotation(annotation.getName(), start, end);
					binding.setName(annotation.getName() + "__" + k	+ "__binding");
					bindings.add(binding);
					break;
				}
				iterations++;
				if (iterations - 1 % 1000 == 0) {
					logger.info("getting hard to find bound midpoints, iterated "+ iterations+ " to find position for "+ k+ " protein in transcript "+ annotation.getName());
				}
			}
		}
		return bindings;
	}
	
	private boolean didNotLandNearExistingProt(List<Annotation> regions, int pos) {
		for (int i = 0; i < regions.size(); i++) {
			if (Math.abs(regions.get(i).getMidpoint() - pos) < MIN_DIST_BETWEEN_BOUND) {
				return false;
			}
		}
		return true;
	}
	
	public static void main(String[] args)throws Exception{
		if(args.length>6){
			Collection<Gene> genes=BEDFileParser.loadData(new File(args[0]));
			int numProteins=new Integer(args[1]);
			int numReads=new Integer(args[2]);
			double enrichment=new Double(args[3]);
			double biasFreq=new Double(args[4]);
			String genomeDir=args[5];
			String save=args[6];
			new ReadSimulator2(genes.iterator().next(), numProteins, numReads, enrichment, biasFreq, genomeDir, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=Genes \n args[1]=num proteins \n args[2]=number of reads \n args[3]=enrichment of bound \n args[4]=2nd read frequency \n args[5]=genome dir \n args[6]=save";
}
