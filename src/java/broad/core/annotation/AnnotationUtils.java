package broad.core.annotation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.Map;
import java.util.Random;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Gene;

import org.w3c.dom.traversal.NodeIterator;

import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.core.annotation.LightweightGenomicAnnotation;
import broad.core.annotation.BasicLightweightAnnotation;
import broad.core.error.ParseException;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.core.math.EmpiricalDistribution;
import broad.pda.annotation.BEDFileParser;
import broad.pda.annotation.GTFFileParser;
import broad.pda.gene.GeneWithIsoforms;

public class AnnotationUtils {
	public static final String USAGE = "Usage: AnnotationUtils TASK=<task_num> <task_args>\n" +
	"\tTasks:\n" +
	"\t\t1. Take the intersection two annotation sets: -set1 <File containing first set>  -set2 <File containing second set> -minScore <to filter by minimum annotation score>" +
	"\n\t\t\t [-set1format <[BED], GFF or generic> -set2format <[BED], GFF or generic> -outformat <desired output format BED, GFF or generic, default is se1format>" +
	"\n\t2. Take the difference of two annotation sets (all annotations in set1 that do not overlap set2): -set1 <File containing first set>  -set2 <File containing second set>  -minScore <to filter by minimum annotation score>" +
	"\n\t\t [-set1format <[BED], GFF or generic> -set2format <[BED], GFF or generic> -outformat <desired output format BED, GFF or generic, default is se1format>" +
	"\n\t3. Print out all elements in one set (set1) that overlap another (set2): -set1 <File containing first set>  -set2 <File containing second set> " +
	"\n\t\t [-set1format <[BED], GFF or generic> -set2format <[BED], GFF or generic> -outformat <desired output format BED, GFF or generic, default is se1format>" +
	"\n\t4. Merge annotation file -in <annotation file (dowa not support standard input input at this point)> -out <output file name or standard output by default> [-format <[BED], GFF or generic>]" +
	"\n\t5. Reformat annotation file  -in <annotation file (dowa not support standard input input at this point)> -out <output file name or standard output by default> [-informat <[BED], GFF or generic> -outformat <BED,[GFF],generic]" +
	"\n\t6. Shuffle scores -in <Annotation file> -out <output file name or standard output by default> [-informat <[BED], GFF or generic>]" +
	"\n\t7. Get overlap statistics of two sets (how much each of the elements in set1 is covered by set2). " +
	"\n\t\t [-set1format <[BED], GFF or generic> -set2format <[BED], GFF or generic> -outformat <desired output format BED, GFF or generic, default is se1format>" +
	"\n\t8. Take the complement of the the first annotation set in the second set (similar to 2 but elements incompletely overlapping will be broken into the regions that don't overlap): -set1 <File containing first set>  -set2 <File containing second set> " +
	"\n\t9. Filter any delimited file by extracting all records that overlap a given annotation file -in <Delimted file> -chrColumn <Column containing chrmosome information (first column is 0)> [-chr <Alternative if there is no chromosome information in the file, it is assumed that all positions are in a given chromosome>] -positionColumn <Column containing position information  (first column is 0)> -annotations <Annotation file> -annotationFormat <[BED], GFF, Generic> [-sepeartor <Character separting columns, tab is assumed by default>]"+
	"\n\t10. Print out all elements in one set (set1) that overlap another (set2) and are in the same orientation as the overlapping gene: -set1 <File containing first set>  -set2 <File containing second set> " +
	"\n\t\t [-set1format <[BED], GFF or generic> -set2format <[BED], GFF or generic> -outformat <desired output format BED, GFF or generic, default is se1format>" +
	"\n\t11. Compute the probability that an interval of given size overlaps a given annotation set by chance: -set <File containing first set>  -sizeFile <Size file with chromosome sizes> " +
	"\n\t12. Print distance metrics: for each element in first set, print the distance to the closest element in the second set." +
	"\n\t\t [-nonOverlapping <find the closest non-overlapping element> -set1format <[BED], GFF or generic> -set2format <[BED], GFF or generic> -chr <limit analysis to chr>" +
	"\n\tSlideAndCount\tSlide a window over the annotation map and count elements that overlap each window" +
	"\n\t\t -window <windowSize> -overlap <overlapSize> -chr <chr> -start <start> -end <end>" +
	"\n\tOverlapSignifiance. Calculate the significance of the overlap between two sets of annotations by permuting across the same chromosome." +
	"\n\t\t -set1 <set1> -set2 <set2> -nperm <number of permutations> -sizeFile <file of chromosome sizes>" +
	"\n";
	
	public static void main (String [] args) throws IOException, ParseException {
		ArgumentMap argMap = CLUtil.getParameters(args, USAGE);
		
		if ("1".equals(argMap.getTask())) {	
			String set1In = argMap.getMandatory("set1");
			String set2In = argMap.getMandatory("set2");
			String set1Format = argMap.containsKey("set1format") ? argMap.get("set1format") : "BED";
			String set2Format = argMap.containsKey("set2format") ? argMap.get("set2format") : "BED";
			double minScore = argMap.containsKey("minScore") ? argMap.getDouble("minScore") : -1;
			
			AnnotationReader<? extends GenomicAnnotation> set1 = AnnotationReaderFactory.create(set1In, set1Format, minScore);
			AnnotationReader<? extends GenomicAnnotation> set2 = AnnotationReaderFactory.create(set2In, set2Format, minScore);
			
			set1.intersect(set2.getAnnotationList());

			Iterator<String> chrIt = set1.getChromosomeAnnotationMap().keySet().iterator();
			BufferedWriter bw = argMap.getOutputWriter();
			while(chrIt.hasNext()) {
				String chr = chrIt.next();
				IntervalTree<? extends GenomicAnnotation> tree = set1.getChromosomeTree(chr);
				Iterator<? extends GenomicAnnotation> annotIt = tree.valueIterator();
				//System.out.println("chromsome: chr" + chr + " tree size " + tree.size() + " does iterator has a value " + annotIt.hasNext());
				while(annotIt.hasNext()) {
					LightweightGenomicAnnotation annot = annotIt.next();
					bw.write(annot.toString());
					bw.newLine();
				}
			}
			bw.close();
		} else if ("2".equals(argMap.getTask())) {	
			String set1In = argMap.getMandatory("set1");
			String set2In = argMap.getMandatory("set2");
			String set1Format = argMap.containsKey("set1format") ? argMap.get("set1format") : "BED";
			String set2Format = argMap.containsKey("set2format") ? argMap.get("set2format") : "BED";
			double minScore = argMap.containsKey("minScore") ? argMap.getDouble("minScore") : -1;
			
			AnnotationReader<? extends LightweightGenomicAnnotation> set1 = AnnotationReaderFactory.create(set1In, set1Format, minScore);
			AnnotationReader<? extends LightweightGenomicAnnotation> set2 = AnnotationReaderFactory.create(set2In, set2Format, minScore);
			
			BufferedWriter bw = argMap.getOutputWriter();
			Iterator<String> chrIt = set1.getChromosomeAnnotationMap().keySet().iterator();
			while(chrIt.hasNext()) {
				String chr = chrIt.next();
				IntervalTree<? extends LightweightGenomicAnnotation> tree = set1.getChromosomeTree(chr);
				IntervalTree<? extends LightweightGenomicAnnotation> set2Tree = set2.getChromosomeTree(chr);
				Iterator<? extends LightweightGenomicAnnotation> annotIt = tree.valueIterator();
				//System.out.println("chromsome: chr" + chr + " tree size " + tree.size() + " does iterator has a value " + annotIt.hasNext());
				while(annotIt.hasNext()) {
					LightweightGenomicAnnotation annot = annotIt.next();
					if(set2Tree == null || !set2Tree.overlappers(annot.getStart(), annot.getEnd()).hasNext()) {
						bw.write(annot.toString());
						bw.newLine();
					}
				}
			}
			
			
			
			bw.close();
		}else if ("3".equals(argMap.getTask())) {	
			String set1In = argMap.getMandatory("set1");
			String set2In = argMap.getMandatory("set2");
			String set1Format = argMap.containsKey("set1format") ? argMap.get("set1format") : "BED";
			String set2Format = argMap.containsKey("set2format") ? argMap.get("set2format") : "BED";
			double set1FilterScore = argMap.containsKey("set1FilterScore") ? argMap.getDouble("set1FilterScore") : Double.MIN_VALUE;
			double set2FilterScore = argMap.containsKey("set2FilterScore") ? argMap.getDouble("set2FilterScore") : Double.MIN_VALUE;
			
			AnnotationReader<? extends GenomicAnnotation> set1 = AnnotationReaderFactory.create(set1In, set1Format, set1FilterScore);
			AnnotationReader<? extends GenomicAnnotation> set2 = AnnotationReaderFactory.create(set2In, set2Format, set2FilterScore);

			set1.filterByOverlap(set2.getAnnotationList());
			BufferedWriter bw = argMap.getOutputWriter();
			Iterator<String> chrIt = set1.getChromosomeAnnotationMap().keySet().iterator();
			while(chrIt.hasNext()) {
				String chr = chrIt.next();
				IntervalTree<? extends GenomicAnnotation> tree = set1.getChromosomeTree(chr);
				Iterator<? extends GenomicAnnotation> annotIt = tree.valueIterator();
				//System.out.println("chromsome: chr" + chr + " tree size " + tree.size() + " does iterator has a value " + annotIt.hasNext());
				while(annotIt.hasNext()) {
					LightweightGenomicAnnotation annot = annotIt.next();
					bw.write(annot.toString());
					bw.newLine();
				}
			}
			bw.close();
		} else if ("4".equals(argMap.getTask())) {	
			String in = argMap.getInput();
			String format = argMap.containsKey("format") ? argMap.get("format") : "BED";
			
			AnnotationReader<? extends GenomicAnnotation> set = AnnotationReaderFactory.create(in, format);
			set.merge();

			BufferedWriter bw = argMap.getOutputWriter();
			Iterator<String> chrIt = set.getChromosomeAnnotationMap().keySet().iterator();
			while(chrIt.hasNext()) {
				String chr = chrIt.next();
				IntervalTree<? extends GenomicAnnotation> tree = set.getChromosomeTree(chr);
				Iterator<? extends GenomicAnnotation> annotIt = tree.valueIterator();
				while(annotIt.hasNext()) {
					LightweightGenomicAnnotation annot = annotIt.next();
					bw.write(annot.toString());
					bw.newLine();
				}
			}
			bw.close();
		}else if ("5".equals(argMap.getTask())) {	
			String in = argMap.getInput();
			String informat = argMap.containsKey("informat") ? argMap.get("informat") : "BED";
			String outformat = argMap.containsKey("outformat") ? argMap.get("outformat") : "GFF";
			
			AnnotationReader<? extends GenomicAnnotation> set = AnnotationReaderFactory.create(in, informat);

			AnnotationFactory<? extends GenomicAnnotation> factory  = AnnotationFactoryFactory.getFactory(outformat);
			BufferedWriter bw = argMap.getOutputWriter();
			Iterator<? extends GenomicAnnotation> it = set.getAnnotationList().iterator();
			while(it.hasNext()) {
				GenomicAnnotation annotation  = it.next();
				LightweightGenomicAnnotation reformattedAnnotation = factory.create(annotation);
				bw.write(reformattedAnnotation.toString());
				bw.newLine();
			}
			bw.close();
		} else if ("6".equals(argMap.getTask())) {	
			String in = argMap.getInput();
			String informat = argMap.containsKey("informat") ? argMap.get("informat") : "BED";
			
			AnnotationReader<? extends GenomicAnnotation> reader = AnnotationReaderFactory.create(in, informat);

			BufferedWriter bw = argMap.getOutputWriter();
			Iterator<String> chrIt = reader.getChromosomeIterator();
			while(chrIt.hasNext()) {
				String chr = chrIt.next();
				reader.shuffleScores(chr);
				Iterator<? extends GenomicAnnotation> it = reader.getChromosomeBEDs(chr).iterator();
				while(it.hasNext()) {
					LightweightGenomicAnnotation annotation  = it.next();
					bw.write(annotation.toString());
					bw.newLine();
				}
			}
			bw.close();
		} else if ("7".equals(argMap.getTask())) {	
			String set1In = argMap.getMandatory("set1");
			String set2In = argMap.getMandatory("set2");
			String set1Format = argMap.containsKey("set1format") ? argMap.get("set1format") : "BED";
			String set2Format = argMap.containsKey("set2format") ? argMap.get("set2format") : "BED";
			
			AnnotationReader<? extends GenomicAnnotation> set1 = AnnotationReaderFactory.create(set1In, set1Format);
			AnnotationReader<? extends GenomicAnnotation> set2 = AnnotationReaderFactory.create(set2In, set2Format);
			
			BufferedWriter bw = argMap.getOutputWriter();
			Iterator<String> chrIt = set1.getChromosomeAnnotationMap().keySet().iterator();
			while(chrIt.hasNext()) {
				String chr = chrIt.next();
				IntervalTree<? extends GenomicAnnotation> tree = set1.getChromosomeTree(chr);
				List<? extends GenomicAnnotation> chrSet1 = set1.getChromosomeBEDs(chr);
				for(GenomicAnnotation annotation : chrSet1) {
					List<? extends GenomicAnnotation> overlappers = set2.getOverlappers(annotation);
					List<Annotation> stitched = BasicGenomicAnnotation.stitchList(overlappers, 2);
					int overlap = 0;
					for (Annotation ga : stitched) {
						overlap  += annotation.getOverlap(ga);
					}
					annotation.setScore(overlap);
					bw.write(annotation.toString());
					bw.newLine();
					//System.out.println("chromsome: chr" + chr + " tree size " + tree.size() + " does iterator has a value " + annotIt.hasNext());
				}
			}
			bw.close();
		}  else if ("8".equals(argMap.getTask())) {	
			String set1In = argMap.getMandatory("set1");
			String set2In = argMap.getMandatory("set2");
			String set1Format = argMap.containsKey("set1format") ? argMap.get("set1format") : "BED";
			String set2Format = argMap.containsKey("set2format") ? argMap.get("set2format") : "BED";
			
			AnnotationReader<? extends GenomicAnnotation> set1 = AnnotationReaderFactory.create(set1In, set1Format);
			AnnotationReader<? extends GenomicAnnotation> set2 = AnnotationReaderFactory.create(set2In, set2Format);
			
			BasicAnnotationReader bar = set2.takeComplement(set1.getAnnotationList());
			BufferedWriter bw = argMap.getOutputWriter();
			Iterator<String> chrIt = bar.getChromosomeAnnotationMap().keySet().iterator();
			while(chrIt.hasNext()) {
				String chr = chrIt.next();
				IntervalTree<? extends GenomicAnnotation> tree = bar.getChromosomeTree(chr);
				Iterator<? extends GenomicAnnotation> annotIt = tree.valueIterator();
				//System.out.println("chromsome: chr" + chr + " tree size " + tree.size() + " does iterator has a value " + annotIt.hasNext());
				while(annotIt.hasNext()) {
					LightweightGenomicAnnotation annot = annotIt.next();
					bw.write("chr" +annot.getChromosome()+"\t"+ annot.getStart()+"\t" + annot.getEnd());
					bw.newLine();
				}
			}
			bw.close();
		} else if ("9".equals(argMap.getTask())) {
			String separator   = argMap.containsKey("separator") ? argMap.get("separator") : "\t";
			String annotFormat = argMap.containsKey("annotationFormat") ? argMap.get("annotationFormat") : "BED";
			String annotationFile = argMap.getMandatory("annotations");
			String chr = argMap.containsKey("chr") ? argMap.get("chr") : null;
			int chrCol = -1;
			if(chr == null) {
				chrCol = argMap.getInteger("chrColumn"); 
			} else {
				chr = chr.replace("chr","");
			}
			int posCol = argMap.getInteger("positionColumn");
			
			
			AnnotationReader<? extends GenomicAnnotation> reader = AnnotationReaderFactory.create(annotationFile, annotFormat);
			
			BufferedReader br = argMap.getInputReader();
			BufferedWriter bw = argMap.getOutputWriter();
			String line = null;
			while( (line = br.readLine())!= null) {
				String [] data = line.split(separator);
				int pos = Integer.parseInt(data[posCol]);
				IntervalTree<? extends GenomicAnnotation> tree =  reader.getChromosomeTree(chr == null ? data[chrCol] : chr);
				Iterator overlapperIt = tree.overlappers(pos, pos+1);
				if(overlapperIt.hasNext()) {
					bw.write(line);
					bw.newLine();
				}
				
			}
			br.close();
			bw.close();	
		}
			//Added by Moran
		else if ("10".equals(argMap.getTask())) {	
			String set1In = argMap.getMandatory("set1");
			String set2In = argMap.getMandatory("set2");
			String set1Format = argMap.containsKey("set1format") ? argMap.get("set1format") : "BED";
			String set2Format = argMap.containsKey("set2format") ? argMap.get("set2format") : "BED";
			
			AnnotationReader<? extends GenomicAnnotation> set1 = AnnotationReaderFactory.create(set1In, set1Format);
			AnnotationReader<? extends GenomicAnnotation> set2 = AnnotationReaderFactory.create(set2In, set2Format);
			
			set1.filterByOverlapAndOrientation(set2.getAnnotationList());
			BufferedWriter bw = argMap.getOutputWriter();
			Iterator<String> chrIt = set1.getChromosomeAnnotationMap().keySet().iterator();
			while(chrIt.hasNext()) {
				String chr = chrIt.next();
				IntervalTree<? extends GenomicAnnotation> tree = set1.getChromosomeTree(chr);
				Iterator<? extends GenomicAnnotation> annotIt = tree.valueIterator();
				//System.out.println("chromsome: chr" + chr + " tree size " + tree.size() + " does iterator has a value " + annotIt.hasNext());
				while(annotIt.hasNext()) {
					LightweightGenomicAnnotation annot = annotIt.next();
					bw.write(annot.toString());
					bw.newLine();
				}
			}
			bw.close();
		}
		
		else if ("11".equals(argMap.getTask())) {	
			String setFile = argMap.getMandatory("set");
			String sizeFile = argMap.getMandatory("sizeFile");
			int regionSize = argMap.getInteger("regionSize");
			int iterations = argMap.containsKey("iterations") ? argMap.getInteger("iterations") : 100000;
			String format = setFile.endsWith(".gff") || setFile.endsWith(".gtf") ? "GFF" : "BED";
			Random r = new Random(); 
			Map<String, Integer> chrSizes =  parseSizeFile(sizeFile, regionSize);
			
			AnnotationReader<? extends GenomicAnnotation> set = AnnotationReaderFactory.create(setFile, format);
			long setBaseCoverage = set.getBaseCoverage();
			long genomeSize = 0;
			for(int size : chrSizes.values()) {
				genomeSize += size;
			}
			
			LinkedHashMap<String, Double> chrProbs = new LinkedHashMap<String, Double>(chrSizes.size());
			double cummProb = 0;
			for(String chr : chrSizes.keySet()) {
				cummProb +=chrSizes.get(chr)/(double)genomeSize;
				chrProbs.put(chr, cummProb);
			}
			
			int overlapped = 0;
			for(int i = 0; i < iterations; i++) {
				
				String chr = getRandomChromosome(chrProbs, r.nextDouble());
				int position = r.nextInt(chrSizes.get(chr) - regionSize);
				
				LightweightGenomicAnnotation region = new BasicLightweightAnnotation(chr, position, position + regionSize); //YES, sometimes the end size will go past the chrSize, but it does not matter
				if (set.getOverlappers(region).size() > 0) {
					overlapped++;
				} 
			}
			System.out.println( "Set coverage: " + (setBaseCoverage)/(double)(genomeSize) + " Probability of overlap = "+(overlapped/(double)iterations));
		}
		
		/* 	
		 * Jesse - June 13, 2012
	   	 * Output the distance from each element in set1 to the closest element in set2
		 * Output format: set1 element + distance + set2 element 
	     */
		else if ("12".equals(argMap.getTask())) {
			String set1In = argMap.getMandatory("set1");
			String set2In = argMap.getMandatory("set2");
			String set1Format = argMap.containsKey("set1format") ? argMap.get("set1format") : "BED";
			String set2Format = argMap.containsKey("set2format") ? argMap.get("set2format") : "BED";
			String chromosome = argMap.containsKey("chr") ? argMap.get("chr") : "all";
			boolean nonOverlapping = argMap.containsKey("nonOverlapping");
			
			AnnotationReader<? extends GenomicAnnotation> set1 = AnnotationReaderFactory.create(set1In, set1Format, chromosome);
			System.out.println("Loaded " + set1.getAnnotationList().size() + " records for set1.");
			AnnotationReader<? extends GenomicAnnotation> set2 = AnnotationReaderFactory.create(set2In, set2Format, chromosome);
			System.out.println("Loaded " + set2.getAnnotationList().size() + " records for set2.");
			
			Map<? extends GenomicAnnotation,LightweightGenomicAnnotation> allClosest = set1.findClosestForAll(set2, nonOverlapping); 
			BufferedWriter bw = argMap.getOutputWriter();
			Iterator<? extends GenomicAnnotation> itr = allClosest.keySet().iterator();
			while (itr.hasNext()) {
				LightweightGenomicAnnotation annot = itr.next();
				LightweightGenomicAnnotation closest = allClosest.get(annot);
				if (closest != null) {
					bw.write(annot.toString() + "\t" + annot.getDistanceTo(closest) + "\t" + closest.toString() + "\n");
				}
			}
			bw.close();
		}	
		
		/**
		 *  Jesse June 19, 2012
		 */
		else if ("SlideAndCount".equals(argMap.getTask())) {
			String setFile = argMap.getInput();
			String setFormat = argMap.containsKey("informat") ? argMap.get("informat") : "BED";
			int window = argMap.getInteger("window");
			int overlap = argMap.getInteger("overlap");
			String chr = argMap.getMandatory("chr");
			int start = argMap.getInteger("start");
			int end = argMap.getInteger("end");
			
			BufferedWriter bw = argMap.getOutputWriter();
			
			// Load set and calculate general statistics
			AnnotationReader<? extends GenomicAnnotation> set = AnnotationReaderFactory.create(setFile, setFormat);
			int totalAnnotations = set.size();
			System.out.print(totalAnnotations);
			
			// Slide over region and count overlapping
			LightweightGenomicAnnotation region = new BasicLightweightAnnotation(chr, start, start + window);
			do {
				List<? extends GenomicAnnotation> overlappers = set.getOverlappers(region);
				int nOverlappers = overlappers.size();
				int nBases = getBaseCoverage(overlappers);
				double pctBases = (double) nBases / (double) window;
				bw.write(region.toUCSC() + "\t" + nOverlappers + "\t" + region.getChromosome() + "\t" + region.getStart() + "\t" + region.getEnd() +
						  "\t\t" + totalAnnotations + "\t" + nBases + "\t" + pctBases + "\n");
				region.setStart(region.getStart() + window - overlap);
				region.setEnd(region.getStart() + window);
			} while (region.getStart() < end);
			
			bw.close();
		}
		
		/**
		 * Jesse June 24, 2012
		 */
		else if ("OverlapSignificance".equals(argMap.getTask())) {
			String set1In = argMap.getMandatory("set1");
			String set2In = argMap.getMandatory("set2");
			String set1Format = argMap.containsKey("set1format") ? argMap.get("set1format") : "BED";
			String set2Format = argMap.containsKey("set2format") ? argMap.get("set2format") : "BED";
			String chromosome = argMap.containsKey("chr") ? argMap.get("chr") : "all";
			int n = argMap.getInteger("nperm");
			String sizeFile = argMap.getMandatory("sizeFile");
			Map<String, Integer> sizes = BEDFileParser.loadChrSizes(sizeFile);
			
			AnnotationReader<? extends GenomicAnnotation> set1 = AnnotationReaderFactory.create(set1In, set1Format, chromosome);
			System.out.println("Loaded " + set1.getAnnotationList().size() + " records for set1.");
			AnnotationReader<? extends GenomicAnnotation> set2 = AnnotationReaderFactory.create(set2In, set2Format, chromosome);
			System.out.println("Loaded " + set2.getAnnotationList().size() + " records for set2.");
			
			// Calculate distances and create an empirical distribution
			List<Integer> distances = new ArrayList<Integer>(set1.getDistancesToClosest(set2).values());
			int numOverlappers = countOverlappers(distances);
			
			// Permute many times and build empirical distribution of distances and overlappers
			Integer[] permutedOverlappers = new Integer[n];
			List<Integer> permutedDistances = new ArrayList<Integer>();
			for (int i = 0; i < n; i++) {
				set1.permuteAnnotationsOnChromosome(sizes);   // permute regions in place
				
				//BufferedWriter tmp = new BufferedWriter(new FileWriter(argMap.getOutput() + ".perm" + i));
				//tmp.write(set1.toStringForWriting());
				//tmp.close();
				
				Collection<Integer> currDistances = set1.getDistancesToClosest(set2).values();
				permutedOverlappers[i] = countOverlappers(currDistances);
				permutedDistances.addAll(currDistances);
			}
			
			// Calculate approximate permutation p-value
			int numMore = 0;
			for (int i = 0; i < permutedOverlappers.length; i++)
				if (permutedOverlappers[i].intValue() >= numOverlappers) numMore++;
			double pvalue = (double) numMore / (double) n;
			
			// Write results
			BufferedWriter bw = argMap.getOutputWriter();
			bw.write("# Size of set1: " + set1.size() + "\n");
			bw.write("# Size of set2: " + set2.size() + "\n");
			bw.write("# Number of elements in set1 that overlap elements in set2: " + numOverlappers + "\n");
			bw.write("# Permutations = " + n + "\n");
			bw.write("# Permutation p-value = " + pvalue + "\n");
			Collections.sort(Arrays.asList(permutedOverlappers));
			for (Integer x : permutedOverlappers)
				bw.write(x + "\n");
			bw.close();
			
			// Write distance files for observed and permuted versions of set1
			BufferedWriter ow = new BufferedWriter(new FileWriter(argMap.getOutput() + ".observed"));
			Collections.sort(distances);
			for (Integer i : distances) ow.write(i + "\n");
			ow.close();
			
			BufferedWriter ew = new BufferedWriter(new FileWriter(argMap.getOutput() + ".permuted"));
			Collections.sort(permutedDistances);
			for (Integer i : permutedDistances) ew.write(i + "\n");
			ew.close();
		}
		
		//Already implemented at GeneTools
		/*else if ("BED2GTF".equalsIgnoreCase(argMap.getTask())){
			String b = argMap.getMandatory("bed");
			String src = argMap.getMandatory("source");
			BufferedWriter bw = argMap.getOutputWriter();
			BEDFileParser bed=new BEDFileParser(b);
			bed.bed2gtf(bw,src);
			bw.close();
			
		}*/
			else {
			System.err.println("Task " + argMap.getTask() + " is invalid or no task was specified");
		}
	}

	private static String getRandomChromosome(LinkedHashMap<String, Double> chrProbs, double p) {
		String chr = null;
		Iterator<String> chrIt = chrProbs.keySet().iterator();
		double prevProb = 0;
		while(chr == null && chrIt.hasNext()) {
			String nextChr = chrIt.next();
			double nextProb = chrProbs.get(nextChr);
			if(prevProb <= p && nextProb > p) {
				chr = nextChr;
			}
		}
		
		return chr;
	}

	private static Map<String, Integer> parseSizeFile(String sizeFile, int regionSize) throws NumberFormatException, IOException {
		Map<String, Integer> map = new LinkedHashMap<String, Integer>();
		BufferedReader br = new BufferedReader(new FileReader(sizeFile));
		String line = null;
		while ((line = br.readLine())!= null) {
			String [] info = line.split("\t");
			int chrLength = Integer.parseInt(info[1]);
			if(chrLength > regionSize) {
				map.put(info[0], Integer.parseInt(info[1]));
			}
		}
		br.close();
		return map;
	}
	
	
	private static int countOverlappers(Collection<Integer> distances) {
		int overlappers = 0;
		for (Integer d : distances)
			if (d.intValue() == 0)
				overlappers++;
		return overlappers;
	}

	
	private static int getBaseCoverage(List<? extends GenomicAnnotation> annotations) {
		int cov = 0;
		Iterator<? extends GenomicAnnotation> it = annotations.iterator();
		while(it.hasNext()) {
			LightweightGenomicAnnotation annotation = it.next();
			cov = cov +  annotation.length();
		}
		return cov;
	}
	
	
	
	
}
