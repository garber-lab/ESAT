package broad.pda.seq.protection;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Random;

import net.sf.picard.sam.BuildBamIndex;
import net.sf.samtools.BAMIndex;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileHeader.SortOrder;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMProgramRecord;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceRecord;
import nextgen.core.alignment.Alignment;
import nextgen.core.alignment.SingleEndAlignment;
import nextgen.core.annotation.Gene;

import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.log4j.Logger;
import org.broad.igv.Globals;

import broad.core.annotation.BED;
import broad.core.annotation.BEDReader;
import broad.core.annotation.BasicLightweightAnnotation;
import broad.core.annotation.GenomicAnnotationFilter;
import broad.core.annotation.LightweightGenomicAnnotation;
import broad.core.error.ParseException;
import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.annotation.BEDFileParser;
import broad.pda.annotation.GTFFileParser;
import broad.pda.feature.genome.Chromosome;
import broad.pda.feature.genome.DirectoryInstalledGenomeAssembly;
import broad.pda.seq.alignment.AlignmentUtils;

/**
 * Utility methods to get simulated reads of protect-seq like data
 * 
 * @author mgarber
 * 
 */

// TODO Add sort and index for visualization purposes
// TODO Make bound regions annotations a bed file so we can see them
// TODO Add RNase bias

public class ReadSimulator {
	private static final Logger logger = Logger.getLogger(ReadSimulator.class
			.getName());
	static final String USAGE = "\nTasks"
			+ "\n simulate \n\t\t-annotationFile <gene annotation file> \n\t\t-genomeDir <Genome directory set up in a directory structure fashion> "
			+ "\n\t\t-annotations <comma separated list of annotations to use> \n\t\t-annotationCoverages <comma separated list of coverages for each of the annotations provided (average read base coverage), list must be of the same size> "
			+ "\n\t\t-secondReadBias <fraction of second reads expected to land at the crosslinking site> "
			+ "\n\t\t-insertSizeMean <Mean of the insert size distribution \n\t\t-insertSizeSD <Standard Deviation of the insert size distribution"
			+ "\n\t\t-readLength <Read length> \n\t\t-outPrefix <output prefix for the genome and transcript centric alignment files generated"
			+ "\n\t\t-boundProteins <comma separated list of number of proteins bound to each of the annotations, list must of the same size as the annotation list>"
			+ "\n\t\t-enrichmentAtBound <Enrichment of protein bound fragments: [(#reads at bound sites)/(total bound region)]/[($reads on transcript)/(transcript lentgh)] "
			+ "\n";

	private static final int MIN_DIST_BETWEEN_BOUND = 50;
	private static final int BOUND_REGION_SIZE = 10;

	public static void main(String[] args) throws Exception {
		Globals.setHeadless(true);
		ArgumentMap argMap = CLUtil.getParameters(args, USAGE, "simulate");

		if ("simulate".equalsIgnoreCase(argMap.getTask())) {
			String annotationFile = argMap.getMandatory("annotationFile");
			String[] annotations = argMap.getMandatory("annotations")
					.split(",");
			String[] annotationCoveragesStr = argMap.getMandatory(
					"annotationCoverages").split(",");
			double[] annotationCoverages = new double[annotationCoveragesStr.length];
			HashMap<String, List<LightweightGenomicAnnotation>> annotationBoundProteins = new HashMap<String, List<LightweightGenomicAnnotation>>(
					annotations.length);
			if (annotations.length != annotationCoveragesStr.length) {
				System.err.println("Number of annotations "
						+ annotations.length
						+ " did not match number of coverages "
						+ annotationCoveragesStr.length);
				return;
			}

			File genomeDir = new File(argMap.getMandatory("genomeDir"));
			String prefix = argMap.getMandatory("outPrefix");
			int readLength = argMap.getInteger("readLength");
			double insertSizeMean = argMap.getDouble("insertSizeMean");
			double insertSizeSD = argMap.getDouble("insertSizeSD");

			double enrichmentAtBound = argMap.getDouble("enrichmentAtBound");
			double pctSecondReadAtCrosslink = argMap
					.getDouble("secondReadBias");

			// set up coverages
			for (int i = 0; i < annotationCoveragesStr.length; i++) {
				annotationCoverages[i] = Double
						.parseDouble(annotationCoveragesStr[i]);
			}

			BEDFileParser annotationParser = annotationFile.endsWith(".gff")
					|| annotationFile.endsWith(".gtf")
					|| annotationFile.endsWith(".GTF") ? new GTFFileParser(
							annotationFile) : new BEDFileParser(annotationFile);
							DirectoryInstalledGenomeAssembly org = new DirectoryInstalledGenomeAssembly(genomeDir);

							List<Gene> annotationRefSeqs = new ArrayList<Gene>();

							for (int i = 0; i < annotations.length; i++) {
								Gene annot = annotationParser.get(annotations[i]);
								if (annot == null) {
									throw new IllegalArgumentException(
											annotations[i]
													+ " was not found in the annotation file you provided "
													+ annotationFile);
								}
								annotationRefSeqs.add(annot);
							}

							if (argMap.containsKey("boundProteins")) {
								annotationBoundProteins = generateBoundRegions(argMap.getMandatory("boundProteins").split(","), annotationRefSeqs);
								writeBindings(annotationBoundProteins, prefix, annotationRefSeqs);
							} else {
								annotationBoundProteins = loadBoundRegions(prefix, annotations);
							}

							RandomDataGenerator rdg = new RandomDataGenerator();

							// make quality string
							byte[] quality = new byte[readLength];
							for (int i = 0; i < readLength; i++) {
								quality[i] = 'I';
							}
							String perfectMatchCigar = readLength + "M";
							// Setting header file for transcriptome alignment
							SAMFileHeader transcriptomeHeader = new SAMFileHeader();
							transcriptomeHeader.addProgramRecord(new SAMProgramRecord("Scripture-simulator"));
							transcriptomeHeader.addComment("Simulated reads");
							transcriptomeHeader.setSortOrder(SortOrder.coordinate);
							for (String annotName : annotations) {
								logger.debug("annotation " + annotName);
								Gene annot = annotationParser.get(annotName);
								transcriptomeHeader.addSequence(new SAMSequenceRecord(annotName, annot.getTranscriptLength()));
							}

							// Setting up genome alignment header file
							SAMFileHeader genomeHeader = new SAMFileHeader();
							genomeHeader.addProgramRecord(new SAMProgramRecord("Scripture-simulator"));
							genomeHeader.addComment("Simulated reads");
							genomeHeader.setSortOrder(SortOrder.coordinate);
							List<Chromosome> chrs = org.getAllNonRandomChromosomes();
							for (Chromosome chr : chrs) {
								genomeHeader.addSequence(new SAMSequenceRecord(chr.getName(), chr.length()));
							}

							// This is as close as the first read may get to the end of the gene
							int firstReadDistToGeneEnd = (int) Math.round(readLength
									+ insertSizeMean + insertSizeSD);

							// Creating alignment writers
							SAMFileWriterFactory alnWriterFactory = new SAMFileWriterFactory();

							File transcriptomeBamFile = new File(prefix + ".transcriptome.bam");
							File genomeBamFile = new File(prefix + ".genome.bam");
							SAMFileWriter transcriptomeAlignmentWriter = alnWriterFactory.makeSAMOrBAMWriter(transcriptomeHeader, false,	transcriptomeBamFile);
							SAMFileWriter genomeAlignmentWriter = alnWriterFactory.makeSAMOrBAMWriter(genomeHeader, false, genomeBamFile);

							Random r = new Random();
							Chromosome lastChr = null;

							String fastaOutput = prefix + ".transcripts.fa";
							File fastaOutFile = new File(fastaOutput);
							if (fastaOutFile.exists()) {
								fastaOutFile.delete();
							}
							FastaSequenceIO fsio = new FastaSequenceIO(fastaOutput);

							for (int i = 0; i < annotationRefSeqs.size(); i++) {
								logger.debug("annotation " + annotationRefSeqs.get(i).getName());
								Gene annot = annotationRefSeqs.get(i);
								List<LightweightGenomicAnnotation> bindings = annotationBoundProteins
										.get(annot.getName());

								int transcriptLength = annot.length();
								Chromosome chr = org.getChromosome(annot.getChr());
								if (lastChr == null
										|| !chr.getSymbol().equals(lastChr.getSymbol())) {
									if (lastChr != null) {
										lastChr.unloadSequence();
									}
									chr.loadSequence();
									lastChr = chr;
								}
								annot.setSequenceFromChromosome(chr.getSequence());

								fsio.append(annot.getSequenceObject(), fastaOutput);

								// coverage = readlength * #reads/transcript length. Thus we
								// need to sample coverage*transcriptLength/readLength reads.
								long totalReadsForTranscript = Math.round(transcriptLength
										* annotationCoverages[i] / readLength); // TODO Consider
								// making this
								// the input
								// param
								// directly

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

								int numBound = bindings.size();
								int numUnbound = transcriptLength - numBound;
								double denom = numUnbound + (enrichmentAtBound * numBound);

								logger.info("numBound=" + numBound + " numUnbound="
										+ numUnbound);

								double pctReadsFromCrosslink = enrichmentAtBound / denom;

								logger.info("p(insert from bound)=" + pctReadsFromCrosslink);

								// OK lets begin generating reads
								for (int j = 0; j < totalReadsForTranscript; j++) {

									int pos = 0;
									double draw = rdg.nextUniform(0, 1);
									if (bindings.size() == 0 || draw > pctReadsFromCrosslink) { // If
										// no  crosslink OR draw from non crosslinked segments
										// while(true) {
										pos = r.nextInt(transcriptLength - firstReadDistToGeneEnd);
										// if(distToBoundProt(pos, boundMidpoints) > insertSizeMean) {
										// break; //TODO: streamline this, it is very inneficient
										// }
										// }
									} else {

										int protein = bindings.size() == 1 ? 0 : rdg.nextInt(0,	bindings.size() - 1);
										boolean secondAtCrosslink = pctSecondReadAtCrosslink > rdg.nextUniform(0, 1);
										if (secondAtCrosslink) {
											pos = (int) bindings.get(protein).getMiddle() + BOUND_REGION_SIZE / 2;
										} else {
											int insertSize = 0;
											while (insertSize <= 0) {
												insertSize = (int) Math.round(rdg.nextGaussian(insertSizeMean, insertSizeSD));
											}
											int lower = Math.max(0, (int) bindings.get(protein).getMiddle() - insertSize);
											int upper = Math.min(transcriptLength - readLength,	(int) bindings.get(protein).getMiddle()	+ insertSize);
											logger.debug("insertSize " + insertSize + " lower "	+ lower + " upper " + upper);
											pos = rdg.nextInt(lower, upper);
										}

									}
									SAMRecord transcriptAln = new SAMRecord(transcriptomeHeader);
									transcriptAln.setReadName("SCR_" + j + "_"	+ System.nanoTime());
									transcriptAln.setReadString(annot.getSequence().substring(pos, pos + readLength));
									transcriptAln.setBaseQualities(quality);
									transcriptAln.setAlignmentStart(pos + 1);
									transcriptAln.setReferenceName(annotations[i]);
									transcriptAln.setCigarString(perfectMatchCigar);
									transcriptAln.setMappingQuality(255);

									int insertSize = 0;
									while (insertSize <= 0) {
										insertSize = (int) Math.round(rdg.nextGaussian(
												insertSizeMean, insertSizeSD));
									}

									int secondReadStart = transcriptAln.getAlignmentStart()
											+ insertSize + readLength;

									if (secondReadStart + readLength < transcriptLength) {
										logger.debug(transcriptAln.getReadName()
												+ " is paired ");
										transcriptAln.setMateNegativeStrandFlag(true);
										transcriptAln.setReadPairedFlag(true);
										transcriptAln.setProperPairFlag(true);
										transcriptAln.setInferredInsertSize(insertSize);
										transcriptAln.setMateReferenceName(transcriptAln
												.getReferenceName());
										transcriptAln.setSecondOfPairFlag(true);
										transcriptAln.setFirstOfPairFlag(false);

										SAMRecord transcriptAlnPair = (SAMRecord) transcriptAln
												.clone();
										transcriptAlnPair.setMateAlignmentStart(transcriptAln
												.getAlignmentStart());
										transcriptAlnPair.setSecondOfPairFlag(false);
										transcriptAlnPair.setFirstOfPairFlag(true);
										transcriptAlnPair.setReadString(/*
										 * Sequence.reverseSequence
										 * (
										 */annot.getSequence()
										 .substring(secondReadStart - 1,
												 secondReadStart + readLength - 1))/* ) */;
										transcriptAlnPair.setReadNegativeStrandFlag(true);
										transcriptAlnPair.setCigarString(perfectMatchCigar);
										transcriptAlnPair.setAlignmentStart(secondReadStart);
										transcriptAlnPair.setProperPairFlag(true);
										transcriptAlnPair.setInferredInsertSize(insertSize);
										transcriptAln.setMateAlignmentStart(transcriptAlnPair
												.getAlignmentStart());
										// transcriptAln.setMateNegativeStrandFlag(true);
										// transcriptAlnPair.setMateNegativeStrandFlag(false);

										transcriptomeAlignmentWriter
										.addAlignment(transcriptAln);
										transcriptomeAlignmentWriter
										.addAlignment(transcriptAlnPair);

										SAMRecord genomeAln = mapSamRecord(transcriptAln, annot, genomeHeader);
										SAMRecord genomeAlnPair = mapSamRecord(transcriptAlnPair, annot, genomeHeader);
										genomeAln.setMateReferenceName(genomeAln
												.getReferenceName());
										genomeAlnPair.setMateReferenceName(genomeAln
												.getReferenceName());
										// logger.debug("Going to write: " +
										// genomeAln.toString());
										genomeAlignmentWriter.addAlignment(genomeAln);
										genomeAlignmentWriter.addAlignment(genomeAlnPair);
									}
									// TODO MG: This is dangerous because it breaks all our
									// methods. The reads either have to all be paired or all
									// unpaired
									// Its OK to set the pair to unmapped
									// Better yet, we should set the pair to map exactly at the
									// end of the transcript
									else {
										logger.debug(transcriptAln.getReadName()
												+ " is NOT paired ");
										transcriptAln.setReadPairedFlag(true);
										transcriptAln.setSecondOfPairFlag(true);
										transcriptAln.setFirstOfPairFlag(false);
										transcriptAln.setMateUnmappedFlag(true);
										transcriptAln.setProperPairFlag(false);
										SAMRecord genomeAln = mapSamRecord(transcriptAln, annot, genomeHeader);
										transcriptomeAlignmentWriter.addAlignment(transcriptAln);
										genomeAlignmentWriter.addAlignment(genomeAln);
									}

								}
							}

							transcriptomeAlignmentWriter.close();
							logger.info("Generating index files");
							File transcriptomeBamIdxFile = new File(
									transcriptomeBamFile.getAbsolutePath()
									+ BAMIndex.BAMIndexSuffix);
							if (transcriptomeBamIdxFile.exists()) {
								transcriptomeBamIdxFile.delete();
							}
							SAMFileReader reader = new SAMFileReader(transcriptomeBamFile);
							BuildBamIndex.createIndex(reader, transcriptomeBamIdxFile);
							reader.close();

							genomeAlignmentWriter.close();
							File genomeBamIfxFile = new File(genomeBamFile.getAbsolutePath()
									+ BAMIndex.BAMIndexSuffix);
							if (genomeBamIfxFile.exists()) {
								genomeBamIfxFile.delete();
							}
							reader = new SAMFileReader(genomeBamFile);
							BuildBamIndex.createIndex(reader, genomeBamIfxFile);

		}
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
				relativeStart=(gene.getSize())-(gene.getEnd());
				relativeEnd=(gene.getSize())-(gene.getStart()); //+1
				flip=true;
			}*/
			Gene geneSAM = gene.copy();

			geneSAM.transcriptToGenomicPosition(relativeStart, relativeEnd); //put back end-1

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
		}/* else {
			SAMSequenceRecord reference = header.getSequence(sam.getReferenceName());
			logger.error("Could not map record gene was: " + gene.toBED() + " original record: " + sam.toString());
			//logger.info("Could not map back to genome read relative start-end: "+ relativeStart + "-" +relativeEnd + " ---- " + sam.getReadName() + " mapped to " + reference.getSequenceName() + " aligned reference length was " + reference.getSequenceLength() + " while annotation length was " + gene.getSize());
		}*/


		else { logger.error("Could not map sam record to genome: gene " + gene + " header has chr? " + header.getSequence(gene.getChr()) + " record exons: " + samGene.getNumExons());}
		return mappedSam;
	}

	private static void writeBindings(HashMap<String, List<LightweightGenomicAnnotation>> annotationBoundProteins,
			String prefix, List<Gene> annotationRefSeqs) throws IOException {
		BufferedWriter tBw = new BufferedWriter(new FileWriter(prefix + ".transcriptome.bound.bed"));
		BufferedWriter gBw = new BufferedWriter(new FileWriter(prefix + ".genome.bound.bed"));

		for (Gene annotation : annotationRefSeqs) {
			List<LightweightGenomicAnnotation> bindings = annotationBoundProteins.get(annotation.getName());
			for (LightweightGenomicAnnotation binding : bindings) {
				binding.setOrientation(".");
				tBw.write(new BED(binding).toShortString());
				tBw.newLine();
				
				Gene toGenome = annotation.transcriptToGenomicPosition(binding.getStart(), binding.getEnd());
				//toGenome.trim(binding.getStart(), binding.getEnd());
				toGenome.setName(binding.getName());
				gBw.write(toGenome.toBED());
				gBw.newLine();
			}
		}

		tBw.close();
		gBw.close();
	}

	private static HashMap<String, List<LightweightGenomicAnnotation>> loadBoundRegions(String prefix, String[] annotations) throws IOException,
			ParseException {
		HashMap<String, List<LightweightGenomicAnnotation>> rtrn = new HashMap<String, List<LightweightGenomicAnnotation>>();
		File boundRegionFile = new File(prefix + ".transcriptome.bound.bed");
		if (boundRegionFile.exists()) {

			BEDReader reader = new BEDReader();
			reader.load(boundRegionFile, new GenomicAnnotationFilter() {

				public boolean accept(LightweightGenomicAnnotation annotation) {
					return true;
				}

				public boolean isEnough(LightweightGenomicAnnotation annotation) {
					return false;
				}
			});

			List<BED> bindings = reader.getAnnotationList();
			for (BED binding : bindings) {
				String bindingName = binding.getName();
				String[] bindingNameInfo = bindingName.split("__");
				String annotationName = bindingNameInfo[0];

				boolean found = false;
				int i = 0;
				while (!found && i < annotations.length) {
					found = annotations[i].equals(annotationName);
					i++;
				}

				if (!found) {
					throw new IllegalArgumentException(
							boundRegionFile
							+ " exists, however it contains bindings for an annotation "
							+ annotationName
							+ " which is not in the list provided, the bindings file may be old. Regenerate bindings by specifying -boundProteins");
				}

				List<LightweightGenomicAnnotation> annotationBindings = null;
				if (!rtrn.containsKey(annotationName)) {
					annotationBindings = new ArrayList<LightweightGenomicAnnotation>();
					rtrn.put(annotationName, annotationBindings);
				}
				annotationBindings = rtrn.get(annotationName);
				annotationBindings.add(binding);

			}
		} else {
			throw new IllegalArgumentException(
					"Could not find previous run "
							+ prefix
							+ ".transcriptome.bound.bed, you must supply a bound proteins paramter -boundProteins so that random bindings may be generated");
		}
		// TODO Auto-generated method stub
		return rtrn;
	}

	private static HashMap<String, List<LightweightGenomicAnnotation>> generateBoundRegions(
			String[] annotationsBindEventsStr, List<Gene> annotations) {
		if (annotations.size() != annotationsBindEventsStr.length) {
			throw new IllegalArgumentException("Number of annotations "
					+ annotations.size()
					+ " did not match number of proteins bound to transcripts "
					+ annotationsBindEventsStr.length);
		}

		Random r = new Random();
		HashMap<String, List<LightweightGenomicAnnotation>> rtrn = new HashMap<String, List<LightweightGenomicAnnotation>>(
				annotations.size());

		for (int i = 0; i < annotationsBindEventsStr.length; i++) {
			int numberOfBindings = Integer.parseInt(annotationsBindEventsStr[i]);
			Gene annotation = annotations.get(i);
			int transcriptLength = annotation.getTranscriptLength();

			// Now we need to place the proteins on the transcript, we will drop
			// proteins if the do not fit.
			int maxProteinsFitting = (int) Math
					.floor((transcriptLength - 2 * MIN_DIST_BETWEEN_BOUND)	/ (double) MIN_DIST_BETWEEN_BOUND) - 1;
			if (numberOfBindings > maxProteinsFitting) {
				logger.info("Requested number of proteins for transcript "	+ annotation.getName() + " is too large, re-setting to " + maxProteinsFitting);
				numberOfBindings = maxProteinsFitting;
			}

			List<LightweightGenomicAnnotation> bindings = new ArrayList<LightweightGenomicAnnotation>(
					numberOfBindings);
			rtrn.put(annotation.getName(), bindings);
			for (int k = 0; k < numberOfBindings; k++) {
				int iterations = 0;
				while (true) {
					int pos = r.nextInt(transcriptLength - MIN_DIST_BETWEEN_BOUND);
					if (pos > MIN_DIST_BETWEEN_BOUND && didNotLandNearExistingProt(bindings, pos)) {
						int start = pos - BOUND_REGION_SIZE / 2;
						int end = pos + BOUND_REGION_SIZE / 2;
						LightweightGenomicAnnotation binding = new BasicLightweightAnnotation(annotation.getName(), start, end);
						binding.setName(annotation.getName() + "__" + k	+ "__binding");
						bindings.add(binding);
						break;
					}
					iterations++;
					if (iterations - 1 % 1000 == 0) {
						logger.info("getting hard to find bound midpoints, iterated "
								+ iterations
								+ " to find position for "
								+ k
								+ " protein in transcript "
								+ annotation.getName());
					}
				}
			}

		}
		return rtrn;
	}

	private static double distToBoundProt(int pos, int[] boundMidpoints) {
		double dist = Double.MAX_VALUE;
		for (int midpoint : boundMidpoints) {
			dist = Math.min(dist, Math.abs(pos - midpoint));
		}
		return dist;
	}

	private static boolean didNotLandNearExistingProt(
			List<LightweightGenomicAnnotation> regions, int pos) {
		for (int i = 0; i < regions.size(); i++) {
			if (Math.abs(regions.get(i).getMiddle() - pos) < MIN_DIST_BETWEEN_BOUND) {

				return false;
			}
		}
		return true;
	}

}