package umms.esat;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import org.apache.log4j.Logger;

import org.broad.igv.Globals;

import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.core.math.Statistics;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.annotation.BEDFileParser;

import net.sf.samtools.util.RuntimeIOException;
import umms.core.annotation.Gene;

public class CollapsedGeneLengths {

	private static String annotationFile;
	private static Map<Gene, Set<Gene>> collapsedGeneMap;
	private static Map<String, IntervalTree<Gene>> collapsedGenes;
	static Logger logger = Logger.getLogger(CollapsedGeneLengths.class.getName());
	static final double MIN_OVERLAP = 0.25;
	private static HashMap<Gene, String> duplicateNameMap;
	static List<String> rows = new ArrayList<String>();
	
	static final String usage = "Usage: CollapsedGeneLengths -task calculate "+
			"\n**************************************************************"+
			"\n\t\tMANDATORY arguments"+
			"\n\n\t\t-in <Annotation file for which to calculate gene lengths. [BED by default]> "+
			"\n\n\t\t-out <output name> ";

	
	public CollapsedGeneLengths(String[] args) throws IOException{
		
		/*
		 * Gives a log4j error. Check later.
		 */
		Globals.setHeadless(true);
		/*
		 * @param for ArgumentMap - size, usage, default task
		 * argMap maps the command line arguments to the respective parameters
		 */
		ArgumentMap argMap = CLUtil.getParameters(args,usage,"score3P");

		/*
		 * Read the annotation file
		 */
		annotationFile = argMap.getInput();		
		/*
		 * Check the format of the annotation file and call the GTF or BED parser accordingly
		 */
		if(annotationFile.endsWith(".gtf") || annotationFile.endsWith(".GTF")){
			logger.error("Please supply an annotation file in the BED format",new RuntimeIOException());
		}
		Map<String,Collection<Gene>> annotations =  BEDFileParser.loadDataByChr(new File(annotationFile));			
		//HashMap of gene name to the number of duplicates
		HashMap<String, Integer> duplicateMap = new HashMap<String, Integer>();
		//HashMap of gene to rowName
		duplicateNameMap = new HashMap<Gene, String>();
		int duplicates=0;
		for(String chr:annotations.keySet()){
			for(Gene gene : annotations.get(chr)) {
				if(!rows.contains(gene.getName())) {
					String name = gene.getName();
					if(duplicateNameMap.containsKey(gene)){
						logger.info("Entry for "+name+" already exists");
					}
					else{
						rows.add(name);
						duplicateMap.put(name, 1);
						duplicateNameMap.put(gene, name);
					}
				} 
				// If the gene name has another annotation
				else {
					
					if(duplicateNameMap.containsKey(gene)){
						logger.info("Entry for "+gene.getName()+" already exists in "+duplicateNameMap.get(gene));
					}
					else{
						if(!duplicateMap.containsKey(gene.getName()))
							duplicateMap.put(gene.getName(), 1);
						//Row name is now the geneName appended with the duplicate number
						duplicateMap.put(gene.getName(), (duplicateMap.get(gene.getName())+1));
						String name = (gene.getName()+"_"+duplicateMap.get(gene.getName()));
						rows.add(name);
						duplicateNameMap.put(gene, name);
						//logger.warn("Duplicated annotation : "+ gene.toBED());
						duplicates++;
					}
				}
			}
		}
		setCollapsedGeneMap(annotations);

		logger.info("Found " + duplicates + " duplicates, ignoring them, going to process " + rows.size() + " annotations");

		calculateGeneLengths(argMap.getOutput());
	}
	
	/**
	 * Calculates and output the minimum, maximum and mean transcript lengths of constituent isoforms used to collapse genes.
	 * @throws IOException 
	 */
	private void calculateGeneLengths(String outputFile) throws IOException {
		
		BufferedWriter bw = new BufferedWriter(new FileWriter(outputFile));

		bw.write("GeneName\tMinimum\tMaximum\tMean\tNumberOfIsoforms\n");
		for(String chr:collapsedGenes.keySet()){
			
			IntervalTree<Gene> tree=collapsedGenes.get(chr);
			Iterator<Gene> geneIt=tree.valueIterator();
			while(geneIt.hasNext()) {
				
				//For each collapsed gene				
				Gene gene = geneIt.next();
				bw.write(gene.getName()+"\t");
				double[] lengths = new double[collapsedGeneMap.get(gene).size()];
				//For all constituent genes
				int i=0;
				for(Gene isoform:collapsedGeneMap.get(gene)){
					lengths[i] = isoform.getSize();
					i++;
				}
				double min  = Statistics.min(lengths);
				double max = Statistics.max(lengths);
				double mean = Statistics.mean(lengths);
				
				bw.write(min+"\t"+max+"\t"+mean+"\t"+lengths.length+"\n");
			}
		}
		bw.close();
	}
	/**
	 * This function collapses the genes
	 * @param annotations
	 * @return
	 */
	private void setCollapsedGeneMap(Map<String,Collection<Gene>> annotations){
		
		collapsedGeneMap = new HashMap<Gene,Set<Gene>>();
		collapsedGenes = new HashMap<String, IntervalTree<Gene>>();
		String geneString = "gene_";
		
		for(String chr:annotations.keySet()){
			IntervalTree<Gene> oldTree = new IntervalTree<Gene>();
			for(Gene g:annotations.get(chr)){
				oldTree.put(g.getStart(), g.getEnd(), g);
			}
			
			IntervalTree<Gene> newTree = new IntervalTree<Gene>();
			collapsedGenes.put(chr, newTree);
			
			//For all genes in old tree
			for(Gene currentElement:oldTree.toCollection()) {
				//Find overlapping genes in new tree
				Iterator<Gene> overlapperIt = newTree.overlappingValueIterator(currentElement.getStart(), currentElement.getEnd());
				Gene overlapper = null;
				int overlapperNum = 0;
				//While there are more overlappers
				while(overlapperIt.hasNext()) {
					Gene overlapperCandidate= overlapperIt.next();
					//If overlapper and gene are compatible
					if(BEDFileParser.isOverlapCompatible(currentElement,overlapperCandidate, MIN_OVERLAP) ) {
						overlapper = overlapperCandidate;
					}				
					if(overlapper != null) {
						overlapperNum++;
						//Remove overlapper from new tree
						newTree.remove(overlapper.getStart(),overlapper.getEnd());
						if(!collapsedGeneMap.containsKey(overlapper)){
							logger.info("Not contains "+overlapper.getName());
						}
						Set<Gene> constituents = collapsedGeneMap.remove(overlapper);
						//Merge gene and overlapper
						Gene mergedElement=new Gene(overlapper.takeUnion(currentElement));
						//Set new name
						mergedElement.setName(overlapper.getName()+"_"+duplicateNameMap.get(currentElement));
						constituents.add(currentElement);
						constituents.add(currentElement);						
						collapsedGeneMap.put(mergedElement, constituents);
						
						//Put merged in new tree
						newTree.put(mergedElement.getStart(), mergedElement.getEnd(), mergedElement);
					} 
					overlapper = null;
				}
				//If no overlappers, add the gene itseld
				if(overlapperNum == 0) {
					//Bug fix : Dec 27 2010 - the function was not merging a refseq gene with isoforms that had several isoforms
					Gene mergedElement=new Gene(currentElement);
					//Set new name
					mergedElement.setName(geneString+duplicateNameMap.get(currentElement));
					Set<Gene> constituents = new TreeSet<Gene>();
					constituents.add(currentElement);
					collapsedGeneMap.put(mergedElement, constituents);
					newTree.put(mergedElement.getStart(), currentElement.getEnd(),  mergedElement);
				}
			}
		}
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		new CollapsedGeneLengths(args);

	}

}
