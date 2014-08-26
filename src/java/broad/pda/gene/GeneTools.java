package broad.pda.gene;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Gene;
import nextgen.core.annotation.Annotation.Strand;

import org.apache.log4j.Logger;

import broad.core.annotation.AnnotationReader;
import broad.core.annotation.AnnotationReaderFactory;
import broad.core.annotation.BED;
import broad.core.annotation.BEDReader;
import broad.core.annotation.BasicLightweightAnnotation;
import broad.core.annotation.GFF;
import broad.core.annotation.GenomicAnnotation;
import broad.core.annotation.LightweightGenomicAnnotation;
import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.annotation.BEDFileParser;
import broad.pda.annotation.GTFFileParser;
import broad.pda.annotation.Locus;
import broad.pda.datastructures.Alignments;
import broad.pda.seq.segmentation.AlignmentDataModel;
import broad.pda.seq.segmentation.AlignmentDataModelStats;
import broad.pda.seq.segmentation.ContinuousDataAlignmentModel;
import broad.pda.seq.segmentation.GenericAlignmentDataModel;
import broad.core.math.Statistics;


public class GeneTools {
	static String usage="Usage: GeneTools -task <task name> "+
	"\n\torient: \n\t\t  -regions <Specific regions to segment (BED format supported only)> \n\t\t-seqdir <Sequence directory (top level of a genome build)> \n\t\t-out <Output file>" +
		"\n\toverlap. Print out all elements in one set (set1) that overlap another (set2) while considering orientation. WORKS IN THE EXON LEVEL. : \n\t\t-set1 <File containing first set> \n\t\t-set2 <File containing second set in BED format, not assuming a full BED file here> \n\t\t -considerOrientation <optional>\n\t\t-out <Output file>" +
		"\n\tdifference. Take the difference of two annotation sets (all annotations in set1 that do not overlap set2). WORKS IN THE EXON LEVEL.: \n\t\t-set1 <File containing first set>  \n\t\t-set2 <File containing second set> \n\t\t -ignoreOrientation <optional> \n\t\t-out <Output file>" +
		"\n\tmerge. Takes all transcripts that overlap in the same orientation and produces a transcript that contains all exons of the overlappig transcripts: \n\t\t-in <Anntotaions to merge> \n\t\t-out <File to write merged annotations>"+
		"\n\tmergeByGtfGeneId. Takes all transcripts that share the same gene_id attribute in a GTF and merge them. \n\t\t -in <GTF> \n\t\t -out<file name> \n\t\t -outFormat <bed/gtf> -src <src if gtf output>\n"+
		"\n\texons. Extract exons from full transcripts: \n\t\t-in <Anntotaions> \n\t\t-out <File to writeextracted exons> \n\t\t -nonRedundant <optional flag if one wants a non-redundant set of exons>"+
		"\n\tintrons. Extract introns from full transcripts: \n\t\t-in <Anntotaions> \n\t\t-out <File to write extracted introns>"+
		"\n\tsampleIntrons. sample introns from full transcripts acording to the size distribution of exons in that sample: \n\t\t-in <Anntotaions> \n\t\t-out <File to write extracted introns>"+
		"\n\tpctOverlap. Reports the maximal percent of the transcripts in set2 that are overlapped by a transcript in set1: \n\t\t-set1 <First set file> \n\t\t-set2 <Second set file> \n\t\t-out <File to write extracted exons>"+
		"\n\tlastExon. Extract Last exons \n\t\t-in <Annotation set> [-three <If the 3' rather than 5' last exon is desired>]" +
		"\n\texonsToGene. Reverse the output of the exons command by reconstructing a full gene from a bed file containing its exons. Follow the conventation: a gene named 'geneX' will be reconstructed from the exons named geneX_0, geneX_1 etc': \n\t\t-in <Anntotaions> \n\t\t-out <File to write BED out file>\n"+
		"\n\ttoGFF. Converts a full BED into GFF format \n\t\t-in <Input file in full BED format> -source <Source name for the GFF> -out <Output file> -bundleOverlapingIsoforms <bundle overlapping isoforms with the same gene ID>\n" +
		"\n\tgtf2bed.\n\t\t -in <GTF> \n\t\t -out \n"+
		"\n\tisoformOverlap. Reports statistics on isoform overlaps in set : \n\t\t-in <First set file> \n\t\t-out <File to write extended merged bed file>\n"+
		"\n\toverlapByGenomicRegion. Print out all elements in one set (set1) that overlap another (set2). WORKS IN THE GENOMIC REGIONS LEVEL. : n\t\t-set1 <File containing first set> \n\t\t-set2 <File containing second set in BED format, not assuming a full BED file here> \n\t\t-out <Output file>\n\t\t-considerOrientation <true/false> \n\t\t -byChr <optional flag; defualt=false; for large files load each chr (1-22,X Y M) separately> \n\t\t -seqDir <name of genome directory that has a sizes file, mendatory if byChr is used> \n\t\t -set2ExtendUtr <length to extend both utr, optional>\n" +
		"\n\tdifferenceByGenomicRegion. Take the difference of two annotation sets (all annotations in set1 that do not overlap set2). WORKS IN THE GENOMIC REGIONS LEVEL.: \n\t\t-set1 <File containing first set>  \n\t\t-set2 <File containing second set> \n\t\t-out <Output file>\n\t\t-considerOrientation <true/false> -byChr <optional flag; defualt=false; for large files load each chr (1-22,X Y M) separately> \n\t\t  -seqDir <name of genome directory that has a sizes file, mendatory if byChr is used> \n" +
		"\n\tPickHighScoringIsoform. pick one high scoring isoforms out of all overlapping isoforms: -in <BED file> -out <BED file>\n" +
		"\n\tdedup. Filter out duplicate transcripts -in <Full BED file> -out <Dedupped transcript file>\n" +
		"\n\tcollapse. Collapse similar transcripts into a single one \n\t\t-in <Full BED or GTF file> \n\t\t-exonFudge <Maximum difference in the start/end of exons> " +
		"\n\textract. Extracts local mRNA coordinates in for all annotations in the input file -in <Full BED file> -mRNAStart <Start in the mRNA> -mRNAEnd <end in the mRNA> -fromThreePrime <If the coordinates are specified from the 3' end of the gene>\n" +
		"\n\tIsCompatible. For each reference transcript in set1 find a compatible isoform in set 2. Report the score of the compatible isoform: \n\t\t-set1 <reference bed file> \n\t\t-set2 <set2 file name> or -setLst <file with the names of sets to scan> \n\t\t -out<out file> \n"+
		"\n\tRmHighIsoformLoci. remove from the BED all isoforms in loci for which the number of isoforms is larger than a threshold: -in <BED file> -threshold <maximal number of isoforms in the remaining set> -out <BED file> \n" +
		"\n\tupdateScore. replace the bedScore of the input with that in the reference: -in <BED file> -ref <newScoreBedFile; must contain all annotation of -in> -out <BED file> \n" +
		"\n\tgetIntronTranscript. return a bed file with the transcript's introns as the exons: \n\t\t -in <BED> \n\t\t -out <BED>\n"+
		"\n\tfilterRepeats. report only transcripts that have no more than X precent ancestral repeats: \n\t\t -in <BED> \n\t\t -AR <ancestral repeat> \n\t\t -overlapThreshold <maximal overlap fraction> \n\t\t -out <BED>\n"+
		"\n\texactIntersection. report only transcripts that are in both sets: \n\t\t -set1 <BED> \n\t\t -set2 <BED> \n\t\t -out <BED>\n"+
		"\n\texactSubstraction. report only transcripts that are in set1 and not set2: \n\t\t -set1 <BED> \n\t\t -set2 <BED> \n\t\t -out <BED>\n"+
		"\n\tgetPromoters - writes gene promoters to a BED file \n\t\t-in <BED file> \n\t\t-out <output promoter file> \n\t\t-radius <Radius around TSS>  OR  -upstream <bp upstream to TSS> -downstream <bp downstream to TSS> \n" +
		"\n\ttoGenes. Takes a set of transcripts and identifies gene loci based on exon overla \n\t\t-in <Transcripts in BED or GTF format> \n\t\t-minOverlap <Minimum exon overlap for two transcripts to be associated with the same gene>"+
		"\n\textractSequence takes a full BED or GTF file and outputs fasta for each record \n\t\t-genes <Annotations in BED or GTF format> \n\t\t-sequenceDir <Directory in which each chromosome in its own directory> -\n\t\tout <output file>"+
		"\n\tsampleUniformaly. sample a random subset of the file. \n\t\t -in <BED> \n\t\t-size <number to sample> \n\t\t-minTranscriptSize<int> \n\t\t-out"+
		"\n\tmakeIntergenicRandomModel . sample uniformly transcripts from the unannotated intergenic space based while maintaining the transcript structures ditribution of a reference bed file. \n\t\t -refBed <BED> \n\t\t -chrSizes <Tab delimited, chr\tsize> \n\t\t -annotTofilter <BED> \n\t\t -centromers <BED> \n\t\t -numRand <int> \n\t\t -outprefix"+
		"\n\t BEDToRefFlat . \n\t\t -in <BED> \n\t\t-out \n"+
		"\n\tcountUniqLoci. Given a BED / GTF file, merge the transcripts in the file to unique loci and print to std out the number of unique loci \n\t\t-in <BED/GTF> \n\t\t-informat <BED/GTF (def is BED)> \n" +
		"\n\t toBED Writes an anontation file in BED format -in <Annotation file, the extension of the file dictates the format it is in> -out <output file or standard out if non is specified>"+
		"\n\tremoveRedundantRecords: Takes a BED/GTF file as input and outputs a BED file without redundant records \n\t\t-in <BED/GTF> \n\t\t-out <BED>"+
		"\n\t mapOrthologosIds. Given a input BED/GTF file, an ortholog id map and a reference bed file,   map transcripts ids in the input file to the id of their ortholog ; outputs .chip file. \n\t\t-in <BED/GTF> \n\t\t-informat <BED/GTF (def is BED)> \n\t\t-orthoMap <tab del thisSpeciesNId\tOtherSpeciesId> \n\t\t-out <outf> \n\t\t-refBed <transcripts for the curr species annotated by the ID1  > \n\t\t -useGeneName (true/false; use the gene symbol in ref GTF for mapping) \n\t\t-inSpMap <map  ID1 to the IDs in the orthoMap> \n"+
		"\n\t regionIntersect. Given a input BED/GTF file, an a region file (in BED format) report the transcripts intersected with the region set. \n\t\t-in <BED/GTF gene models> \n\t\t-regions <regions in BED format> -out <Intersected genes output file> \n"+
		"\n\t getLengths. Given a input BED/GTF file, an a region file (in BED format) report the lengths of all genes\n\t\t-in <BED/GTF gene models> \n\t\t-regions <regions in BED format> -out <Table of gene names and legnths> \n"+
		"\n";
	
	static final String DONOR_SITE="GT";
	static final String ACCEPTOR_SITE="AG";
		
	static final String RV_DONOR_SITE = "AC";
	static final String RV_ACCEPTOR_SITE = "CT";
	
	static final String NC_DONOR_SITE="GC";
	static final String NC_ACCEPTOR_SITE="AG";
	
	static final String RV_NC_DONOR_SITE="GC";
	static final String RV_NC_ACCEPTOR_SITE="CT";
	
	static final String NC2_DONOR_SITE="AT";
	static final String NC2_ACCEPTOR_SITE="AC";
	
	static final String RV_NC2_DONOR_SITE="AT";
	static final String RV_NC2_ACCEPTOR_SITE="GT";
	
	static Logger logger = Logger.getLogger(GeneTools.class.getName());
	
	public static void main(String [] args) throws Exception  {
		ArgumentMap argmap = CLUtil.getParameters(args, usage, "orient");
		
		if("orient".equalsIgnoreCase(argmap.getTask())) {
			String regionFile = argmap.getMandatory("regions");
			String seqdir = argmap.getMandatory("seqdir");
			String out = argmap.getOutput();
			
			Map<String, Collection<Gene>> regionMap =BEDFileParser.loadDataByChr(new File(regionFile));
			Map<String, Sequence> m = FastaSequenceIO.loadSequencesByNameFromDirectory(new File(seqdir));
			for(String cStr : regionMap.keySet()) {
				Sequence c = m.get(cStr);
				if(c != null) {
				} else {
					System.err.println("Skipping chromosome " + cStr + " we have no sequence for it");
					continue;
				}
				String cSeq = c.getSequenceBases();
				c.unloadSequence();
				
				Collection<Gene> chrRegions = regionMap.get(cStr);
				for(Gene g : chrRegions) {
					
					Alignments[] introns = g.getIntronsBlocks();
					if(introns.length == 0) {
						g.setOrientation(Strand.UNKNOWN);
					} else {
						String orientation = orientationFromSpliceSites(introns[0], c);
						//System.err.println("Orientation " + introns[0].toUCSC() + ": " + orientation);
						for(int i = 1; i < introns.length; i++) {
							String intronOrientation = orientationFromSpliceSites(introns[i], c);
							//System.err.println("Orientation " + introns[i].toUCSC() + ": " + intronOrientation);
							if("*".equals(orientation)) {
								orientation = intronOrientation;
							} else if (!"*".equals(intronOrientation) && !orientation.equals(intronOrientation)) {
								orientation = "*";
								break;
							}
						}
						if(orientation.equals("+")) g.setOrientation(Strand.POSITIVE);
						else if(orientation.equals("-")) g.setOrientation(Strand.NEGATIVE);
						else g.setOrientation(Strand.UNKNOWN);
					}
				}
			}
			
			writeFullBED(out, regionMap);
		} else if ("regionIntersect".equalsIgnoreCase(argmap.getTask())) {
			String geneFile = argmap.getInput();
			String regionFile = argmap.getMandatory("regions");
			boolean ignoreOrientation = argmap.containsKey("ignoreStrand");
			String out = argmap.getOutput();
			BEDFileParser geneParser = geneFile.endsWith(".gff") || geneFile.endsWith(".gtf") ? new GTFFileParser(geneFile) : new BEDFileParser(geneFile);
			BEDReader regionReader = new BEDReader(regionFile);
			
			BEDFileParser overlap = out.endsWith(".gtf") ? new GTFFileParser() : new BEDFileParser();
			Iterator<String> chrIt = geneParser.getChromosomeIterator();
			while(chrIt.hasNext()) {
				String chr = chrIt.next();
				logger.debug(chr);
				Iterator<GeneWithIsoforms> geneIt = geneParser.getChrTree(chr).valueIterator();
				while(geneIt.hasNext()){
					Gene g = geneIt.next();
					List<Annotation> overlapperRegions = regionReader.getOverlappers(g.getAlignment());
					logger.debug("gene: " + g.toBED() + " has " + overlapperRegions.size() + " overlappers");
					//System.err.println(g.toBED() + " has overlappers? " + !overlapperTree.isEmpty());	
					if(overlapperRegions.isEmpty()) {
						logger.debug("NO overlapper");
						overlap.addRefSeq(g);
					} else  {
						Gene geneFromRegions = new Gene(overlapperRegions);
						logger.debug("gene from regions: " + geneFromRegions.toBED());
						overlap.addRefSeq(g.getOverlap(geneFromRegions));
					}
				}
					
			}
			overlap.writeFullBed(out);
	
		}  else if ("overlap".equalsIgnoreCase(argmap.getTask())) {
			String set1In = argmap.getMandatory("set1");
			String set2In = argmap.getMandatory("set2");
			boolean ignoreOrientation = argmap.containsKey("ignoreStrand");
			String out = argmap.getOutput();
			
			BEDFileParser set1Parser = set1In.endsWith(".gff") || set1In.endsWith(".gtf") ? new GTFFileParser(set1In) : new BEDFileParser(set1In);
			BEDFileParser set2Parser = set2In.endsWith(".gff") || set2In.endsWith(".gtf") ? new GTFFileParser(set2In) : new BEDFileParser(set2In);
			
			BEDFileParser overlap = out.endsWith(".gtf") ? new GTFFileParser() : new BEDFileParser();
			Iterator<String> chrIt = set1Parser.getChromosomeIterator();
			while(chrIt.hasNext()) {
				String chr = chrIt.next();
				System.err.println(chr);
				Iterator<GeneWithIsoforms> geneIt = set1Parser.getChrTree(chr).valueIterator();
				while(geneIt.hasNext()){
					Gene g = geneIt.next();
					IntervalTree<GeneWithIsoforms> overlapperTree = set2Parser.getOverlappers(g);
					//System.err.println(g.toBED() + " has overlappers? " + !overlapperTree.isEmpty());	
					if(!overlapperTree.isEmpty()) {
						overlap.addRefSeq(g);
					} else if(ignoreOrientation) {
						Gene gCopy = g.copy(); //TODO: implement copy on RefSeqWithIsoforms
						g.setOrientation("+".equals(g.getOrientation()) ? "-" : "+");
						overlapperTree = set2Parser.getOverlappers(gCopy);
						if(!overlapperTree.isEmpty()) {
							overlap.addRefSeq(g);
						}
					}
				}
					
			}
			overlap.writeFullBed(out);
	
		}else if ("difference".equalsIgnoreCase(argmap.getTask())) {
			String set1In = argmap.getMandatory("set1");
			String set2In = argmap.getMandatory("set2");
			String out = argmap.getOutput();
			boolean ignoreOrientation = argmap.containsKey("ignoreOrientation");
			
			Map<String, Collection<Gene>> set1 =BEDFileParser.loadDataByChr(new File(set1In));
			BEDFileParser set2= new BEDFileParser(set2In);
			
			Map<String, Collection<Gene>> difference = new LinkedHashMap<String, Collection<Gene>>();
			for(String chr : set1.keySet()) {
				Set<Gene> chrOverlap = new TreeSet<Gene>();
				difference.put(chr, chrOverlap);
				Collection<Gene> chrList = set1.get(chr);
				for(Gene g : chrList) {
					IntervalTree<GeneWithIsoforms> overlapperTree = set2.getOverlappers(g);
					
					if(overlapperTree.isEmpty()) {
						boolean add = true;
						if(ignoreOrientation ) {
							g.setOrientation("+".equals(g.getOrientation()) ? "-" : "+");
							IntervalTree<GeneWithIsoforms> reverseOverlapperTree = set2.getOverlappers(g);
							add = reverseOverlapperTree.isEmpty();
						}
						if(add) {
							chrOverlap.add(g);
						}
					}
				}
				
			}
			writeFullBED(out, difference);			
		} else if ("merge".equalsIgnoreCase(argmap.getTask())) {
			String data = argmap.getInput();
			
			Map<String, Collection<Gene>> set = BEDFileParser.loadDataByChr(new File(data));
			
			Map<String, IntervalTree<Gene>> mergedSet = merge(set);
			
			BufferedWriter bw = argmap.getOutputWriter();
			for(String  chr : mergedSet.keySet() ) {
				writeFullBED(bw, mergedSet.get(chr));
			}
			bw.close();
		} else if ("exons".equalsIgnoreCase(argmap.getTask())) {
			String data = argmap.getInput();
			boolean nonRedundant=argmap.containsKey("nonRedundant")? true : false ;
			boolean oldVesrion=argmap.containsKey("oldVersion")? true : false ;
			BufferedWriter bw = argmap.getOutputWriter();
			if (oldVesrion==false){
				BEDFileParser set = new BEDFileParser(data);
				Set<Annotation> nonRedunSet= new HashSet<Annotation>();
				Collection<Annotation> RedunSet= new LinkedList<Annotation>();
				for(Gene g : set.GetGenes()) {
					if (nonRedundant)
						nonRedunSet.addAll(g.getExonSet());
					else
						RedunSet.addAll(g.getExonSet());
				}
			   if (nonRedundant){
				   for(Annotation exon : nonRedunSet) 
							bw.write(new BED(exon).toShortString()+"\n");
			   }
			   else{
				   for(Annotation exon : RedunSet) 
						bw.write(new BED(exon).toShortString()+"\n");
			   }
				bw.close();
			}
			else{ //old version of the code
				System.err.println("Run old code");
			Map<String, Collection<Gene>> set = BEDFileParser.loadDataByChr(new File(data));
			for(String  chr : set.keySet() ) {
				for(Gene g : set.get(chr)) {
					Set<? extends Annotation> exons = g.getExonSet();
					for(Annotation exon : exons) {
						bw.write(new BED(exon).toShortString());
						bw.newLine();
					}
				}
			}
			bw.close();
			}

		}
			else if ("sampleIntrons".equalsIgnoreCase(argmap.getTask())) {
				String data = argmap.getInput();
				BEDFileParser set = new BEDFileParser(data);
				BufferedWriter bw = argmap.getOutputWriter();
				sampleIntrons(set,bw);
				bw.close();
		}
			else if ("pctOverlap".equalsIgnoreCase(argmap.getTask())) {
			String set1In = argmap.getMandatory("set1");
			String set2In = argmap.getMandatory("set2");
			String out = argmap.getOutput();
			
			Map<String, Collection<Gene>> set1 = BEDFileParser.loadDataByChr(new File(set1In));
			Map<String, Collection<Gene>> set2 = BEDFileParser.loadDataByChr(new File(set2In));
			
			Map<String, IntervalTree<Gene>> set1merged = merge(set1);
			
			for(String chr : set2.keySet()) {
				Collection<Gene> chrList = set2.get(chr);
				IntervalTree<Gene> set1Tree = set1merged.get(chr);
				
				for(Gene g : chrList) {
					double pctOverlap = 0;
					if(set1Tree != null) {
						Iterator<Node<Gene>> overlapperIt = set1Tree.overlappers(g.getStart(), g.getEnd());
						while(overlapperIt.hasNext()) {
							Gene overlapper = overlapperIt.next().getValue();
							pctOverlap = Math.max(pctOverlap, g.percentOverlapping(overlapper));
						}
					}
					g.addExtraField(String.valueOf(pctOverlap));
				}
				
			}
			writeFullBED(out, set2);			
		}
		
		//Moran 17/3/2010 - TODO: make this work!

		else if ("togff".equalsIgnoreCase(argmap.getTask())) {
			String source = argmap.getMandatory("source");
			String file = argmap.getInput();
			Boolean bundle = argmap.containsKey("bundleOverlapingIsoforms") ;
			Boolean cufflinksForm = argmap.containsKey("cufflinksForm") ;
			BufferedWriter bw = argmap.getOutputWriter();
			BEDFileParser data = new BEDFileParser(file);
			if (bundle){
				BEDFileParser mergedData = new BEDFileParser(file);
				mergedData.merge();
				data.bed2gtf_bundleOverlappers(bw,source,mergedData);
			}
			else if (cufflinksForm){
				data.bed2CuffGtf(bw,source, false);
			}
			else
				data.bed2gtf(bw,source);
			bw.close();
		} else if ("tobed".equalsIgnoreCase(argmap.getTask())) {
			String in = argmap.getInput();
			BEDFileParser parser = in.endsWith(".gff") || in.endsWith(".gtf") ? new GTFFileParser(in) : new BEDFileParser(in);
			BufferedWriter bw = argmap.getOutputWriter();
			parser.writeFullBed(bw);
			bw.close();
		
		}else if ("gtf2bed".equalsIgnoreCase(argmap.getTask())){
			String file = argmap.getInput();
			BufferedWriter bw = argmap.getOutputWriter();
			BEDFileParser bed= new BEDFileParser();
			bed.loadGTF(new File(file),"");
			bed.writeFullBed(bw,false);
			//bed.bed2CuffGtf(bw, "", true); //test
			bw.close();
		}
		else if ("gtf2cuffGTF".equalsIgnoreCase(argmap.getTask())){
			String file = argmap.getInput();
			BufferedWriter bw = argmap.getOutputWriter();
			BEDFileParser bed= new BEDFileParser();
			bed.loadGTF(new File(file),"");
			bed.bed2CuffGtf(bw, "", true); //test
			bw.close();
		}
		else if ("mergeByGtfGeneId".equalsIgnoreCase(argmap.getTask())){
			String inFile = argmap.getInput();
			String outFormat=argmap.getMandatory("outformat");
			String outFile=argmap.getOutput();
			String src=argmap.containsKey("src")? argmap.get("src"):"unknown";
			mergeByGtfGeneId(inFile,outFile,outFormat,src,true);
			
	    	  
	    }
		//Merge set and calculate statistics on the overlap of isoform with the merged transcript. 
		//report by an extended BED file. (Moran, May 17th)
		else if ("isoformOverlap".equalsIgnoreCase(argmap.getTask())) {
			String data = argmap.getInput();
			//Map<String, Collection<RefSeqGene>> set = BEDFileParser.loadDataByChr(new File(data));
			//Map<String, IntervalTree<RefSeqGene>> mergedSet = merge(set);
			BEDFileParser set = new BEDFileParser(data);
			BEDFileParser mergedSet = new BEDFileParser(data);
			mergedSet.merge();
			String save = argmap.getOutput();
			Map<Gene, double[]> scores=isoformOverlapStatistics(mergedSet, set);
			
			writeExtendedFullBED(save, scores);
			
		
				
		
		} else if (argmap.getTask().toLowerCase().contains("dedup")) {
			String data = argmap.getInput();
			BEDFileParser set = new BEDFileParser(data);
			set.dedup();
			BufferedWriter bw = argmap.getOutputWriter();
			set.writeFullBed(bw);
			bw.close();
			
		} else if (argmap.getTask().equalsIgnoreCase("extract") ){ 
			String data = argmap.getInput();
			BEDFileParser set = new BEDFileParser(data);
			int mRNAStart = argmap.getInteger("mRNAStart");
			int mRNAEnd   = argmap.getInteger("mRNAEnd");
			boolean from3P = argmap.containsKey("fromThreePrime");
			BufferedWriter bw = argmap.getOutputWriter();
			Iterator<String> chrIt = set.getChromosomeIterator();
			while(chrIt.hasNext()) {
				String chr = chrIt.next();
				Iterator<GeneWithIsoforms> geneIt = set.getChrTree(chr).valueIterator();
				while(geneIt.hasNext()) {
					GeneWithIsoforms gene = geneIt.next();
					int start = from3P ? gene.getTranscriptLength() - mRNAEnd : mRNAStart;
					int end = from3P ? gene.getTranscriptLength() - mRNAStart : mRNAEnd;
					Gene mappedRegion = gene.trimGene(start, end);
					System.err.println("Gene " + gene.getName() + " - " + gene.toUCSC() + " region " + start + "-" + end+ " has trim: " + mappedRegion);
					mappedRegion.setOrientation(gene.getOrientation());
					mappedRegion.setName(gene.getName());
					bw.write(mappedRegion.toBED());
					bw.newLine();
				}
			}
			bw.close();
		} else if("IsCompatible".equalsIgnoreCase(argmap.getTask())){
			String set1In = argmap.getMandatory("set1");
			String set2In = argmap.getMandatory("set2");
			BufferedWriter bw = argmap.getOutputWriter();
			findCompatibleTranscripts (set1In,set2In,bw);
			bw.close();
			
		}
		else if ("updateScore".equalsIgnoreCase(argmap.getTask())){
			
			BEDFileParser set = new BEDFileParser(argmap.getMandatory("set"));
			BEDFileParser ref = new BEDFileParser(argmap.getMandatory("ref"));
			BufferedWriter bw = argmap.getOutputWriter();
			updateScore (set,ref,bw);
			bw.close();
			
		}
		
       else if ("getIntronTranscript".equalsIgnoreCase(argmap.getTask())){
			
    	   String data = argmap.getInput();
			BEDFileParser set = new BEDFileParser(data);
			BufferedWriter bw = argmap.getOutputWriter();
				for(Gene g : set.GetGenes()) {
					Gene g2=g.getIntronTranscript();
					if (g2!= null)
						bw.write(g2.toBED()+"\n");
				}
				
			bw.close();
			
		}
       else if ("filterRepeats".equalsIgnoreCase(argmap.getTask())){
    	   String infile = argmap.getInput();
    	   String Rfile= argmap.getMandatory("Repeat");
    	   Double maxT=argmap.getDouble("overlapThreshold");
    	   BufferedWriter bw = argmap.getOutputWriter();
    	   filterRepeats(infile,Rfile,maxT,bw);
    	   bw.close();
    	   
       }
       else if ("exactIntersection".equalsIgnoreCase(argmap.getTask())){
    	   String set1 = argmap.getMandatory("set1");
    	   String set2= argmap.getMandatory("set2");
    	   BufferedWriter bw = argmap.getOutputWriter();
    	   exactIntersection (set1,set2,bw,true);
    	   bw.close(); 
       }
       else if ("exactSubstraction".equalsIgnoreCase(argmap.getTask())){
    	   String set1 = argmap.getMandatory("set1");
    	   String set2= argmap.getMandatory("set2");
    	   BufferedWriter bw = argmap.getOutputWriter();
    	   exactIntersection (set1,set2,bw,false);
    	   bw.close(); 
       } 
       else if (argmap.getTask().toLowerCase().contains("getpromoters")) {
    	   int radius = argmap.containsKey("radius")?argmap.getInteger("radius"): 0 ;
    	   int upstream =argmap.containsKey("upstream")?argmap.getInteger("upstream"): 0 ;
    	   int downstream =  argmap.containsKey("downstream")?argmap.getInteger("downstream"): 0 ;
    	   if (upstream ==0  &  radius==0)
    	   		{ System.err.println("Provide either radius or upstream/downstream extension");}
    	   else {
    		   if (radius != 0){
    			   upstream = radius; downstream=radius;
    		   }
    		   
	    	   BEDFileParser geneSet = new BEDFileParser(argmap.getInput());
	    	   Iterator<String> chrIt = geneSet.getChromosomeIterator();
	    	   Set<Alignments> promoterSet = new TreeSet<Alignments>();
	    	   while(chrIt.hasNext()) {
	    		   String chr = chrIt.next();
	    		   Iterator<GeneWithIsoforms> geneIt = geneSet.getChrTree(chr).valueIterator();
	    		   while(geneIt.hasNext()) {
	    			   promoterSet.add(geneIt.next().getPromoter(upstream,downstream));
	    		   }
	    	   }
    	  
       
	    	   BufferedWriter bw = argmap.getOutputWriter();
	    	   for (Alignments p : promoterSet) {
	    		   bw.write(p.getChr());
	    		   bw.write("\t");
	    		   bw.write(String.valueOf(p.getStart()));
	    		   bw.write("\t");
	    		   bw.write(String.valueOf(p.getEnd()));
	    		   bw.write("\t");
	    		   bw.write(p.getName());
	    		   bw.newLine();
	    		   
	    	   }
	    	   bw.close();
    	   }
       }
       else if("togenes".equalsIgnoreCase(argmap.getTask())) {
    	   String annotationFile = argmap.getInput();
    	   double minOverlap = argmap.getDouble("minOverlap");
    	   BEDFileParser annotationParser =  annotationFile.endsWith(".gtf") || annotationFile.endsWith(".GTF")? new GTFFileParser(annotationFile) : new BEDFileParser(annotationFile);
    	   annotationParser.makeGenes(minOverlap);
    	   Iterator<String> chrIt = annotationParser.getChromosomeIterator();
    	   BufferedWriter bw = argmap.getOutputWriter();
    	   while(chrIt.hasNext()) {
    		   Iterator<GeneWithIsoforms> geneIt = annotationParser.getChrTree(chrIt.next()).valueIterator();
    		   while(geneIt.hasNext()) {
    			   bw.write(geneIt.next().toBED());
    			   bw.newLine();
    		   }
    	   }
    	   bw.close();
       } else if ("extractSequence".equalsIgnoreCase(argmap.getTask())) {
    	   int cache = 250;
    	   String annotationFile = argmap.getMandatory("genes");
    	   File sequenceDir = new File(argmap.getMandatory("sequenceDir"));
    	   String out = argmap.getOutput();
    	   BEDFileParser annotationParser =  annotationFile.endsWith(".gtf") || annotationFile.endsWith(".GTF")? new GTFFileParser(annotationFile) : new BEDFileParser(annotationFile);
    	   Map<String, Sequence> m=FastaSequenceIO.loadSequencesByNameFromDirectory(sequenceDir);
    	   Iterator<String> chrIt = annotationParser.getChromosomeIterator();
    	   FastaSequenceIO fsio = new FastaSequenceIO();
    	   
    	   while(chrIt.hasNext()) {
    		   String chr = chrIt.next();
    		   System.err.println("Processing "+chr);
    		   Sequence c = m.get(chr);
    		   IntervalTree<GeneWithIsoforms> chrGeneTree = annotationParser.getChrTree(chr);
    		   if(c != null && chrGeneTree != null) {
	    		   Iterator<GeneWithIsoforms> chrGeneIt = chrGeneTree.valueIterator();
	    		   int i = 0;
	    		   ArrayList<Sequence> tmpSequences = new ArrayList<Sequence>();
	    		   while(chrGeneIt.hasNext()) {
	    			   Gene gene = chrGeneIt.next();
	    			   System.err.println("Gene " + gene.getName() + " " + gene.toUCSC());
	    			   gene.setSequenceFromChromosome(c);
	    			   Sequence geneSeq = gene.getSequenceObject();
	    			   tmpSequences.add(geneSeq);
	    			   if(i>0 && i % cache == 0) {
	    				   System.err.println("Writing sequences " + (i-cache) + " to " + i);
	    				   fsio.append(tmpSequences, out);
	    				   tmpSequences.clear();
	    			   }
	    			   i++;
	    		   }
	    		   fsio.append(tmpSequences, out);
	    		   c.unloadSequence();
    		   }
    	   }
       }
       
       else if (argmap.getTask().equalsIgnoreCase("sampleUniformaly")){
    	   	String in = argmap.getInput();
    	   	Double size = argmap.getDouble("sampleSize");
    	   	Double minTranscriptSize = argmap.getDouble("minTranscriptSize");
			BufferedWriter bw = argmap.getOutputWriter();
			sampleUniformaly(in,size,minTranscriptSize,bw);
			bw.close();
       }
		
       else if (argmap.getTask().equalsIgnoreCase("makeIntergenicRandomModel")){
    	   
    	   BEDFileParser refBed = new BEDFileParser (argmap.getMandatory("refBed")); 
    	   String chrSizes = argmap.getMandatory("chrSizes");
    	   BEDFileParser  annotTofilter  = new BEDFileParser (argmap.getMandatory("annotTofilter")); 
    	   BEDFileParser  centromers = new BEDFileParser (argmap.getMandatory("centromers")); 
    	   int numRand =  argmap.getInteger("numRand");
    	   String outprefix = argmap.getMandatory("outprefix");
    	   String crossWithBed =argmap.containsKey("crossWithBed")? argmap.get("crossWithBed"):"";
    	   
    	   transcriptsNullModel nullmodel = new transcriptsNullModel (annotTofilter,centromers,chrSizes); 
    	   nullmodel.makeRandTranscriptList(numRand, refBed);
    	   if (crossWithBed.equalsIgnoreCase(""))
    	   		nullmodel.printRandModels(outprefix);
    	   else
    		   nullmodel.CrossRandModelsWithBed(crossWithBed,refBed,outprefix+"ransomOverlaps.txt");
    	   //nullmodel.writeIntergenicGenome(outprefix+"intergenicSpace.bed");
   		
    	   
       }
       else if (argmap.getTask().equalsIgnoreCase("BEDToRefFlat")){
    	   BEDFileParser bed = new BEDFileParser (argmap.getInput());
    	   BufferedWriter bw = argmap.getOutputWriter();
    	   
    	   for (Gene g: bed.GetGenes()){
    		   bw.write(g.toBED()+"\n");
    	   }
    	   bw.close();
       }
      else if ("removeRedundantRecords".equalsIgnoreCase(argmap.getTask())){
    	   String inputFile = argmap.getInput();
    	   String outputFile = argmap.getOutput();
    	   
    	   BEDFileParser annotations =  inputFile.endsWith(".gtf") || inputFile.endsWith(".GTF")? new GTFFileParser(inputFile) : new BEDFileParser(inputFile);
    	   List<Gene> allGenes = annotations.GetGenes();
    	   List<Gene> unduplicatedGenes = new ArrayList<Gene>();
    	   
    	   for(Gene gene:allGenes){
    		   boolean flag = false;
    		   for(Gene entry:unduplicatedGenes){
    			   if(gene.equals(entry)){
    				   flag = true;
    				   break;
    			   }
    		   }
    		   if(!flag){ 
    			   unduplicatedGenes.add(gene);
    		   }
    		   else{
    			   logger.info(gene.getName()+" is removed for the entry ");
    		   }
    		   
    	   }
    	   
    	   BufferedWriter bedBw = new BufferedWriter(new FileWriter(outputFile));
    	   for(Gene gene:unduplicatedGenes){
    		   bedBw.write(gene.toBED());
    		   bedBw.newLine();
    	   }
    	   bedBw.close();
       }
       else if("getlengths".equalsIgnoreCase(argmap.getTask())) {
    	   String annotationFile = argmap.getInput();
    	   BEDFileParser annotationParser =  annotationFile.endsWith(".gtf") || annotationFile.endsWith(".GTF")? new GTFFileParser(annotationFile) : new BEDFileParser(annotationFile);
    	   Iterator<String> chrIt = annotationParser.getChromosomeIterator();
    	   BufferedWriter bw = argmap.getOutputWriter();
    	   while(chrIt.hasNext()) {
    		   Iterator<GeneWithIsoforms> geneIt = annotationParser.getChrTree(chrIt.next()).valueIterator();
    		   while(geneIt.hasNext()) {
    			   GeneWithIsoforms gene = geneIt.next();
    			   bw.write(gene.getName()+"\t"+gene.getTranscriptLength());
    			   bw.newLine();
    		   }
    	   }
    	   bw.close();
       }
       
		else{System.err.println(usage);}

 }
	


	
	public static Map<String, IntervalTree<Gene>> merge(Map<String, Collection<Gene>> originalMap) {

		Map<String, IntervalTree<Gene>> mergedTreeMap = new LinkedHashMap<String, IntervalTree<Gene>>();

		Iterator<String> chrIt = originalMap.keySet().iterator();
		while(chrIt.hasNext()) {
			String chr = chrIt.next();
			Collection<Gene> original = originalMap.get(chr);
			mergedTreeMap.put(chr, merge(original));
		}

		return mergedTreeMap;
	}
	
	public static IntervalTree<Gene> merge(Collection<Gene> original){
		IntervalTree<Gene> newTree = new IntervalTree<Gene>();
		for(Gene currentElement : original) {
			Iterator<Gene> overlapperIt = new IntervalTree.ValuesIterator<Gene>(newTree.overlappers(currentElement.getStart(), currentElement.getEnd()));
			Gene overlapper = null;
			while(overlapperIt.hasNext() && overlapper == null) {
				Gene overlapperCandidate= overlapperIt.next();
				if(overlapperCandidate.getOrientation().equals(currentElement.getOrientation() )) {
					overlapper = overlapperCandidate;
				}
			}
			if(overlapper != null) {
				newTree.remove(overlapper.getStart(),overlapper.getEnd());
				Gene mergedElement= overlapper.takeUnion(currentElement);
				if (mergedElement.getStart() >=mergedElement.getEnd())
					System.err.println ("Start>end: " +mergedElement.toBED());
				else
					newTree.put(mergedElement.getStart(), mergedElement.getEnd(), mergedElement);
			} else {
				if (currentElement.getStart() >=currentElement.getEnd())
					System.err.println ("Start>end: " +currentElement.toBED());
				else
					newTree.put(currentElement.getStart(), currentElement.getEnd(), currentElement);
			}
		}
		return newTree;
	}

	
	
	public static String orientationForGene(Gene g, Sequence cSeq){
			Alignments[] introns = g.getIntronsBlocks();
			if(introns.length == 0) {
				g.setOrientation(Strand.UNKNOWN);
			} else {
				String orientation = orientationFromSpliceSites(introns[0], cSeq);
				for(int i = 1; i < introns.length; i++) {
					String intronOrientation = orientationFromSpliceSites(introns[i], cSeq);
					if("*".equals(orientation)) {
						orientation = intronOrientation;
					} else if (!"*".equals(intronOrientation) && !orientation.equals(intronOrientation)) {
						orientation = "*";
						break;
					}
				}
			
				if(orientation.equals("+")) g.setOrientation(Strand.POSITIVE);
				else if(orientation.equals("-")) g.setOrientation(Strand.NEGATIVE);
				else g.setOrientation(Strand.UNKNOWN);
				
			}
			return g.getOrientation().toString();
	}
	
	
	public static String orientationFromSpliceSites(Annotation introns, Sequence seq, boolean onlyCanonical) {
		String orientation = "*";
		
		if(seq != null) {
		
			String left = seq.getSubSequence("UP", introns.getStart(), introns.getStart() + 2).getSequenceBases().toUpperCase();
			String right = seq.getSubSequence("UP",introns.getEnd() - 2, introns.getEnd()).getSequenceBases().toUpperCase();
			
			//ystem.err.println("intron " + introns.toUCSC() + " left " + left + " right " + right);
			
			if(left.contains(DONOR_SITE) && right.contains(ACCEPTOR_SITE)){ //transcript in in direct orientation
				orientation = "+";
			} 
			
			else if(left.contains(RV_ACCEPTOR_SITE) && right.contains(RV_DONOR_SITE)){ //transcript in indirect orientation
				orientation = "-";
			} 	
			else if(onlyCanonical){orientation="*";}
			//MG: Added support for non-canonical splicing
			else if((left.contains(NC_DONOR_SITE) && right.contains(NC_ACCEPTOR_SITE)) || (left.contains(NC2_DONOR_SITE) && right.contains(NC2_ACCEPTOR_SITE))){
				orientation = "+";
			}
			else if((left.contains(RV_NC_ACCEPTOR_SITE) && right.contains(RV_NC_DONOR_SITE)) || (left.contains(RV_NC2_ACCEPTOR_SITE) && right.contains(RV_NC2_DONOR_SITE))){
				orientation="-";
			}
		} else {
			throw new IllegalArgumentException("ERROR: can't determine orientation from splice sites, sequence was null" );
		}
		
		return orientation;
	}
	
	public static String orientationFromSpliceSites(Annotation introns, Sequence chrSeq) {
		return orientationFromSpliceSites(introns, chrSeq, false);
	}
	
	private static void writeFullBED(String save, Map<String, Collection<Gene>>  genes)throws IOException{
		FileWriter writer=new FileWriter(save);
		for(String chr : genes.keySet()) {
			Collection<Gene> chrGenes = genes.get(chr);
			for(Gene gene: chrGenes){
				writer.write(gene.toBED()+"\n");
			}
		}
		writer.close();
	}
	
	
	
	private static void writeFullBED(BufferedWriter bw, IntervalTree<Gene>  genes)throws IOException{
		Iterator<Gene> geneIt = genes.valueIterator();
		while(geneIt.hasNext()) {
			bw.write(geneIt.next().toBED());
			bw.newLine();
		}
	}
	
	
	//Moran 18/3/2010
	/*TODO- debug exon addition : Boolean res=chrGenes.get(prefix).addExon(exon);

	private static Map<String, IntervalTree<RefSeqGene>> joinExonsToGene(
			Map<String, Collection<RefSeqGene>> exonSet) {
		
		Map<String, IntervalTree<RefSeqGene>> geneTreeMap = new LinkedHashMap<String, IntervalTree<RefSeqGene>>();
        
		for(String  chr : exonSet.keySet() ) {
			
			IntervalTree<RefSeqGene> newTree = new IntervalTree<RefSeqGene>();
			geneTreeMap.put(chr, newTree);
			Map <String,RefSeqGene> chrGenes= new LinkedHashMap<String,RefSeqGene> ();
			for(RefSeqGene exon : exonSet.get(chr)) {
				
				//take the prefix of the exon name and store it as part of the chr genes
				//if a second exon from the same gene was not added yet
				String exonName=exon.getName();
				String [] strs=exonName.split("_");
				String prefix= strs[0];
				exon.setName(prefix);
				if (! chrGenes.containsKey(prefix)){
					chrGenes.put(prefix, exon);
				}
				else {
					Boolean res=chrGenes.get(prefix).addExon(exon);
					if (!res){
						exon.setName(exonName);
						chrGenes.put(prefix, exon);
					}				
				}
			
			}// end for exons per chr
			
		//map to interval tree
			Iterator<RefSeqGene> geneIt = chrGenes.values().iterator();
			while(geneIt.hasNext()) {
				RefSeqGene g=geneIt.next();
				newTree.put(g.getStart(),g.getEnd(),g);
			}
		}//end chr
		return geneTreeMap;
		
	}

*/

 //return all elements in set 1 that overlap an element in set2
	public static Map<String, IntervalTree<Gene>> overlap
	(Map<String, IntervalTree<Gene>> set1,Map<String, IntervalTree<Gene>>set2, Boolean byExon) {
		
		Map<String, IntervalTree<Gene>> overlap = new LinkedHashMap<String, IntervalTree<Gene>>();
		for(String chr : set1.keySet()) {
			IntervalTree<Gene> chrOverlap = new IntervalTree<Gene>();
			overlap.put(chr, chrOverlap);
			IntervalTree<Gene> set1List = set1.get(chr);
			IntervalTree<Gene> set2List = set2.get(chr);
			if (set2List==null) continue;
			
			for(Gene g : set1List.toCollection()) {
				Iterator<Node<Gene>> overlapperIt = set2List.iterator(g.getStart(), g.getEnd());
				boolean overlaps = false;
				while(overlapperIt.hasNext() && !overlaps) {
					if (byExon){overlaps = g.overlapsExon(overlapperIt.next().getValue());}
					else {overlaps = g.overlapsGene(overlapperIt.next().getValue());}
				}
				if(overlaps) {
					chrOverlap.put(g.getStart(), g.getEnd(), g);
				}
			}
			
		}
		
		
	return overlap;
		
		
	}




	public static Map<String, IntervalTree<Gene>> mergeSets(
			Map<String, IntervalTree<Gene>> set1,
			Map<String, IntervalTree<Gene>> set2) {
		
		if (set1.isEmpty()) return set2;
		if (set2.isEmpty()) return set1;
		Map<String, IntervalTree<Gene>> mergedSet= new LinkedHashMap<String, IntervalTree<Gene>>();
		Set<String> allChrs=new TreeSet();
		allChrs.addAll(set1.keySet());
		boolean b=allChrs.addAll(set2.keySet());
		
		//for (String chr:set2.keySet()){if (allChrs.add(chr)) System.err.println(chr); /**/;} 
		for (String chr:allChrs ){
			IntervalTree<Gene> T= new IntervalTree<Gene>();
			if (set1.containsKey(chr))for (Gene e:set1.get(chr).toCollection()) {T.put(e.getStart(), e.getEnd(), e);}
			if (set2.containsKey(chr)) {for (Gene e:set2.get(chr).toCollection()) {T.put(e.getStart(), e.getEnd(), e);}}
			mergedSet.put(chr, T);
		}
		
		return mergedSet;
	}


	private static void writeExtendedFullBED(String save, Map<Gene, double[]> map )throws IOException{
		FileWriter writer=new FileWriter(save);
		
		for(Gene align: map.keySet()){
			double[] ps=map.get(align);
			
				writer.write(align.toBED());
				for(double p : ps) {
					writer.write("\t"+p);
				}
				writer.write("\n");
			
			
		}
		
		writer.close();
	}
	

	private static Map<Gene, double[]> isoformOverlapStatistics(BEDFileParser mergedSet,BEDFileParser set) {
		
		Map<Gene, double[]> rtrn=new TreeMap<Gene, double[]>();
		//go over every chr
		Iterator<String> chrIt = mergedSet.getChromosomeIterator();
		while(chrIt.hasNext()) {
			String chr = chrIt.next();
			IntervalTree<GeneWithIsoforms> mergedTree=mergedSet.getChrTree(chr);
			Iterator <GeneWithIsoforms> geneIt = mergedTree.valueIterator();
		//for each merged transcript get all of its overlapping isoforms and calc statistics
			while(geneIt.hasNext()) { 
				Gene mergedTranscript=geneIt.next();
				double[] scores= calcOverlapStatistics(mergedTranscript,set);
				rtrn.put(mergedTranscript, scores);
			}
		}
		return rtrn;
	}
	
	/*private static Map<RefSeqGene, double[]> isoformOverlapStatistics(
			Map<String, IntervalTree<RefSeqGene>> mergedSet,
			Map<String, Collection<RefSeqGene>> set) {
		
		Map<RefSeqGene, double[]> rtrn=new TreeMap<RefSeqGene, double[]>();
		//sort set into an interval tree
		Map<String, IntervalTree<RefSeqGene>> sortSet=sortBed(set);
		//go over every chr
		Iterator<String> chrIt = mergedSet.keySet().iterator();
		while(chrIt.hasNext()) {
			String chr = chrIt.next();
			IntervalTree<RefSeqGene> mergedTree=mergedSet.get(chr);
			Iterator <RefSeqGene> geneIt = mergedTree.valueIterator();
		//for each merged transcript get all of its overlapping isoforms and calc statistics
			while(geneIt.hasNext()) { 
				RefSeqGene mergedTranscript=geneIt.next();
				double[] scores= calcOverlapStatistics(mergedTranscript,sortSet);
				rtrn.put(mergedTranscript, scores);
						}
		}
		return rtrn;
	}*/


	private static double[] calcOverlapStatistics(Gene transcript,BEDFileParser set) {	
		double[] rtrn= new double[6];
		int cntr=0;
		int doubleExon=0;
		IntervalTree<GeneWithIsoforms> isoformsTree = set.getChrTree(transcript.getChr());
		Iterator<GeneWithIsoforms> overlapperIt = new IntervalTree.ValuesIterator<GeneWithIsoforms>(isoformsTree.overlappers(transcript.getStart(), transcript.getEnd()));
		List<Double> lst=new ArrayList<Double>();
		while(overlapperIt.hasNext() ) {
			Collection<Gene> overlapperCandidateSet= overlapperIt.next().getAllIsoforms();
			for (Gene overlapperCandidate : overlapperCandidateSet  ){
				if(overlapperCandidate.getOrientation().equals(transcript.getOrientation() )) {
					double pct=transcript.percentOverlapping(overlapperCandidate);
					lst.add(pct);
					cntr ++;
					if (overlapperCandidate.getNumExons()==2)
						doubleExon++;
				}
			}
		}
		
		rtrn[0]=cntr;
		rtrn[1]=Statistics.min(Statistics.toDoubleArray(lst));
		rtrn[2]=Statistics.max(Statistics.toDoubleArray(lst));
		rtrn[3]=Statistics.median(Statistics.toDoubleArray(lst));
		rtrn[4]=Statistics.mean(Statistics.toDoubleArray(lst));
		rtrn[5]=doubleExon;
		return rtrn;
	}





	/*private static BEDFileParser selectHighScoringIsoforms(BEDFileParser set,BEDFileParser mergedSet) {
		
		BEDFileParser rtrn=new BEDFileParser ();
		set.updateScrToBedScore();
		//go over every chr
		Iterator<String> chrIt = mergedSet.getChromosomeIterator();
		while(chrIt.hasNext()) {
			String chr = chrIt.next();
			IntervalTree<GeneWithIsoforms> mergedTree=mergedSet.getChrTree(chr);
			Iterator <GeneWithIsoforms> geneIt = mergedTree.valueIterator();
		//for each merged transcript get all of its overlapping isoforms and calc statistics
			while(geneIt.hasNext()) { 
				Gene mergedTranscript=geneIt.next();
				Gene res=selectHighScoringIsoform(mergedTranscript,set);
				if (res !=null)
					rtrn.addRefSeq(res);
			}
		}
		rtrn.updateScrToBedScore();
		return rtrn;
		
	}*/



	


	public static Gene selectHighScoringIsoform(Gene transcript, BEDFileParser set) {
		
		Gene currBest=null;
		IntervalTree<GeneWithIsoforms> isoformsTree = set.getChrTree(transcript.getChr());
		if (isoformsTree==null || ! isoformsTree.doesOverlap(transcript.getStart(), transcript.getEnd()) )
			return null;
		Iterator<GeneWithIsoforms> overlapperIt = new IntervalTree.ValuesIterator<GeneWithIsoforms>(isoformsTree.overlappers(transcript.getStart(), transcript.getEnd()));
		List<Double> lst=new ArrayList<Double>();
		double currScr=0;
		double currLength=0;
		Gene longestCandidate=null;
		while(overlapperIt.hasNext() ) {
			Collection<Gene> overlapperCandidateSet= overlapperIt.next().getAllIsoforms();
			
			for (Gene overlapperCandidate : overlapperCandidateSet  ){
				//keep track of the longest isoform, in case not isoform spans more than 1/2 the reference
				if (overlapperCandidate.getTranscriptLength() > currLength){
					longestCandidate=overlapperCandidate;
					currLength=longestCandidate.getTranscriptLength();
				}
				if(overlapperCandidate.getOrientation().equals(transcript.getOrientation() ) || overlapperCandidate.getOrientation().equals(Strand.UNKNOWN) ) {
					if (overlapperCandidate.getTranscriptLength() >= (0.5* transcript.getTranscriptLength())&& (overlapperCandidate.getBedScore() >currScr )){
						currBest=overlapperCandidate;
						currScr=overlapperCandidate.getBedScore();
					}
				}
			}
			//If there where no isoforms that where at least half of the 
			if (currScr==0){
				//System.err.println("Printed small segment overlapperCandidate: "+ longestCandidate.getName());
				currBest=longestCandidate;
			}
		}
		
		return currBest;
	}


	/*private static BEDFileParser RmHighIsoformLoci(BEDFileParser set,BEDFileParser mergedSet, Integer threshold) {
		
		BEDFileParser rtrn=new BEDFileParser ();
		set.updateScrToBedScore();
		//go over every chr
		Iterator<String> chrIt = mergedSet.getChromosomeIterator();
		while(chrIt.hasNext()) {
			String chr = chrIt.next();
			Iterator <GeneWithIsoforms> geneIt = mergedSet.getChrTree(chr).valueIterator();
		//for each merged transcript get all of its overlapping isoforms , if they are lower than threshold, add to rtrn
			while(geneIt.hasNext()) { 
				Gene mergedTranscript=geneIt.next();
				IntervalTree<GeneWithIsoforms> overlappers=set.getOverlappers(mergedTranscript);
				if (overlappers.size()<=threshold)
					rtrn.addRefSeqSet(overlappers.toCollection());
			}
		}
		
		return rtrn;
	}*/
	
	

	private static Map<String, IntervalTree<Gene>> sortBed(Map<String, Collection<Gene>> set) {
		
		Map<String, IntervalTree<Gene>> sortTreeMap = new LinkedHashMap<String, IntervalTree<Gene>>();
		Iterator<String> chrIt = set.keySet().iterator();
		while(chrIt.hasNext()) {
			String chr = chrIt.next();
			IntervalTree<Gene> newTree = new IntervalTree<Gene>();
			sortTreeMap.put(chr, newTree);
			Collection<Gene> oldTree = set.get(chr);
			for(Gene currentElement : oldTree) {
				newTree.put(currentElement.getStart(), currentElement.getEnd(), currentElement);
			}
		}
		
		return sortTreeMap;
	}
		

	
	private static void findCompatibleTranscripts(String set1In, String set2In,BufferedWriter bw) throws IOException {
		
		BEDFileParser set1=new BEDFileParser(set1In);
		BEDFileParser set2=new BEDFileParser(set2In);
		
		Iterator<String> chrIt=set1.getChromosomeIterator();
		Map <GeneWithIsoforms,Double> compatibleMap=new HashMap<GeneWithIsoforms,Double>();
		Map <GeneWithIsoforms,Integer> compatibleLengthMap=new HashMap<GeneWithIsoforms,Integer>();
		
		double compatible=0.0;
		double fullyCompatible=0.0;
		double refSize=0.0;
		
		while(chrIt.hasNext()){
			String chr=chrIt.next();
			IntervalTree<GeneWithIsoforms> chrTree= set1.getChrTree(chr);
			
			Iterator<GeneWithIsoforms> gIt=chrTree.valueIterator(); 
			while(gIt.hasNext()){
				GeneWithIsoforms gene=gIt.next();
				refSize++;
				IntervalTree<GeneWithIsoforms> overlapTree=set2.getOverlappers(gene);
				int[] numIntrons=new int[1];
				Gene compatibleGene= gene.findCompatibleGenes(overlapTree,numIntrons);
				Double scr=new Double(-1.0);
				Integer  scrLength=new Integer(0);
				if (compatibleGene!= null){
					scr=compatibleGene.getBedScore();
					scrLength=numIntrons[0];
					compatible++;
					if (scrLength==(gene.getNumExons()-1))
						fullyCompatible++;
					
				}
				compatibleMap.put(gene, scr);
				compatibleLengthMap.put(gene, scrLength);
			}
		}
		
		bw.write("1.2\n");
		bw.write(compatibleMap.size()+"\t"+2+"\n");
		bw.write("Name\tDescription\tRefScore\tCompatibleIsoScore\n");
		for (GeneWithIsoforms g:compatibleMap.keySet()){
			//bw.write(g.getName()+"\t"+g.getName()+"\t"+g.getBedScore()+"\t"+compatibleMap.get(g)+"\n");
			bw.write(g.getName()+"\t"+g.getName()+"\t"+(g.getNumExons()-1)+"\t"+compatibleLengthMap.get(g)+"\n");
		}
		
		System.err.println("Ref gene size: "+refSize +" Comaptible: "+compatible+" ("+(compatible/refSize)+") "+ "Fully comaptible: " +fullyCompatible+" ("+(fullyCompatible/refSize)+")");
	}

	private static void sampleIntrons(BEDFileParser set, BufferedWriter bw) throws IOException {
		
		for(Gene g : set.GetGenes()) {
			LinkedList<Annotation> introns = new LinkedList<Annotation>();
			LinkedList<Annotation> exons = new LinkedList<Annotation>();
			exons.addAll(g.getExonSet());
			introns.addAll( g.getIntronSet());
			for( int i=0; i<introns.size(); i++) {
				Annotation intron =introns.get(i);
				Annotation exon=exons.get(i);
				int delta = intron.getSize() - exon.getSize();
				if (delta <= 0)
					bw.write(new BED(intron).toShortString());
				else{
					int  randAdd= (int)Math.floor(Math.random()*(delta-1)); 
					Alignments nIntron= new Alignments (intron);
					nIntron.setStart(intron.getStart()+randAdd);
					nIntron.setEnd(intron.getStart()+randAdd+exon.getSize()-1);
					bw.write(new BED(nIntron).toShortString());
				}
					
				bw.newLine();
			}
		}
		
	}

	
	
	private static void updateScore(BEDFileParser set, BEDFileParser ref,BufferedWriter bw) throws IOException {
		
     for (Gene g: set.GetGenes()){
    	 Gene newG=ref.getExactIsoform(new GeneWithIsoforms(g));
    	 if (newG != null)
    		 g.setBedScore(newG.getBedScore());
    	 else
    		 g.setBedScore(Double.NaN);
    	 bw.write(g.toBED()+"\n");
     }
		
	}
	
	private static void filterRepeats(String infile, String Rfile,Double maxT, BufferedWriter bw) throws IOException {
		
		BEDFileParser bed=new BEDFileParser(infile);
		Iterator<String> chrIt=bed.getChromosomeIterator();
		while(chrIt.hasNext()){
			String chr=chrIt.next();
			BEDFileParser rBed= new BEDFileParser(Rfile,chr);
			rBed.merge(); //collapse to a uniq set of regions that is covered by AR
			for (GeneWithIsoforms gAll: bed.getChrTree(chr).toCollection()){
				for (Gene g : gAll.getAllIsoforms()){
					double pctOverlap= 0;
					Iterator<GeneWithIsoforms> it= rBed.getOverlappers(g).valueIterator();
					while(it.hasNext()){
						Gene r = it.next();
						pctOverlap+=g.percentOverlapping(r);
					}
					if (pctOverlap <= maxT)
						bw.write(g.toBED()+"\n");
				}
			}
		}
	}


	//reports the intersection if isIntersection==true, and the complement set otherwise.
	private static void exactIntersection(String set1, String set2,BufferedWriter bw, boolean isIntersection) throws IOException {
		
		BEDFileParser bed1=new BEDFileParser(set1);
		BEDFileParser bed2=new BEDFileParser(set2);
		BEDFileParser outBed=exactIntersection (bed1,bed2,isIntersection);
		outBed.writeFullBed(bw);
	}

	private static BEDFileParser exactIntersection(BEDFileParser bed1, BEDFileParser bed2,boolean isIntersection){
		
		BEDFileParser outBed=new BEDFileParser();
		Collection <Gene> allGenes=bed1.GetGenes();
		boolean wasFound=false;
		for (Gene g :allGenes){
			wasFound=false;
			Iterator<GeneWithIsoforms> overlappersIt=bed2.getOverlappers(g).valueIterator();
			while(overlappersIt.hasNext()){
				Collection <Gene> allI= overlappersIt.next().getAllIsoforms();
				for (Gene I:allI){
					if (g.equals(I)){
						wasFound=true;
						break;
					}
				}
				if (wasFound)
					break;
			}
			
			if (wasFound & isIntersection)
				outBed.addRefSeq(g);
			
			if (wasFound==false & isIntersection==false)
				outBed.addRefSeq(g);
			
				
		}
		return outBed;
	}
  
	public static BEDFileParser mergeByGtfGeneId(String inFile, String outFile,	String outFormat, String src,boolean writeMode) throws Exception {
		
		BEDFileParser inbed= new BEDFileParser();
		inbed.loadGTF(new File(inFile),"");
		BEDFileParser outBed =  inbed.mergeByGtfGeneId();
		if (writeMode){
			BufferedWriter  bw =new BufferedWriter(new FileWriter(outFile));
		if (outFormat.equalsIgnoreCase("bed"))
			
			outBed.writeFullBed(bw,false);
		else
			outBed.bed2gtf(bw, src);
		
		bw.close();
		}
		
		return outBed;
	}
	
	
	
	public static  BEDFileParser intersectGeneSets (Collection<Gene> set1 ,Collection<Gene> set2,boolean takeComplement ){
		
		BEDFileParser b1 = new BEDFileParser();
		BEDFileParser b2 = new BEDFileParser();
		b1.addRefSeqSet(set1);
		b2.addRefSeqSet(set2);
		BEDFileParser outBed=exactIntersection (b1,b2,!takeComplement);
		return outBed;
	}
	
	
	private static void sampleUniformaly(String in, Double size,
			Double minTranscriptSize, BufferedWriter bw) throws IOException {
		
		BEDFileParser bed = new BEDFileParser(in);
		List<Gene> genes = bed.GetGenes();
		if (minTranscriptSize >0 ){
			List<Gene> tmp= bed.GetGenes();
			genes = new LinkedList<Gene> ();
			for (Gene g : tmp ){
				if (g.getTranscriptLength()> minTranscriptSize)
					genes.add(g);
			}
		}
		double t = size/(new Double(genes.size()));
		System.err.println("size of cleane genes " + genes.size() +"threshlod "+ t );
		if (t<0.1)
			t=0.3;
		int i=0;
		double r =0;
		for (Gene g: genes ){
			 r = Math.random();
			 
			if (r <= t){
				//System.err.println("r is " + r +"threshlod "+ t );
				bw.write(g.toBED() +"\n");
				i++;
				
			}
			if (i == size)
				break;
		}
		
	}
	
	
	//Map the ids in the infile to their ortholog's id based on a map file and a reference bed file
		private static void mapOrthologosIds(String informat, String infile,
				String mapfile, String reffile, String outfile, String inSpMapf, boolean useGeneName) throws IOException {
			
			BEDFileParser bed ;
			if (informat.equalsIgnoreCase("BED")){
				   bed = new BEDFileParser (infile);
				   bed.merge();
		   }
		   else{
			   bed = new GTFFileParser(infile);
			   bed=bed.mergeByGtfGeneId();    	
		   }
		    
		   HashMap<String, String> orthoMap =  loadMap (mapfile,0,1);
		   
		   HashMap<String, String> inSpMap = new HashMap<String, String>() ;
		   BEDFileParser refbed= new BEDFileParser();
		   
		   if (! useGeneName) {
		    inSpMap =  loadMap (inSpMapf,1,2);
		    refbed = new BEDFileParser (reffile);
		   }
			
			BufferedWriter bw = new BufferedWriter (new FileWriter (outfile +".chip"));
			bw.write("probesetId\tgenSymbol\tGeneInfo\n");
				   
			List<Gene> set1genes = bed.GetGenes();
			int cntr =0;
			String orthoName, curr_name, ref_name,ref_name2;
			
			for (Gene g: set1genes){
				curr_name = g.getName();
				orthoName ="NA";
				ref_name2="";
				if (! useGeneName) {
					IntervalTree<GeneWithIsoforms> overlaps= refbed.getOverlappers(g);
					if (!overlaps.isEmpty()){
						Iterator <GeneWithIsoforms> it= overlaps.valueIterator();
						while (it.hasNext()){
							GeneWithIsoforms og = it.next();
							ref_name  = og.getName();
							if (inSpMap.containsKey(ref_name)) {
								ref_name2 =  inSpMap.get(ref_name);
								break;
							}
						}
					}
				}
				else //use the geneName
				{
					ref_name2 = g.getAttribute("gene_name");//.toUpperCase();
				}
				
				if ( orthoMap.containsKey(ref_name2))
				{
						orthoName = orthoMap.get(ref_name2);
						cntr++;
								
				}
				bw.write(curr_name+"\t"+orthoName+"\t"+orthoName+"\n");
					
			}
			
			bw.close();
			System.err.println("number of genes to map " + set1genes.size() );
			System.err.println("total mapped in species: " +  inSpMap.size() + " total mapped between species:" + orthoMap.size());
			System.err.println("number of input genes successfuly mapped to an ortholog" + cntr );
		}
		
		
		private static HashMap<String, String> loadMap(String mapfile, int i, int j) throws IOException {
			
			HashMap<String, String> map = new  HashMap<String, String>();
			BufferedReader br = new BufferedReader (new FileReader (mapfile));
			String nextLine;
			while ((nextLine = br.readLine()) != null) {
				String[] arr = nextLine.split("\t") ;
				if ((! map.containsKey(arr[i]) ||  !map.get(arr[i]).equalsIgnoreCase("")  ) & arr.length>j  )
				 { 	map.put(arr[i],arr[j]);}
				
			}
			br.close();
			return map;
		}
	


}

