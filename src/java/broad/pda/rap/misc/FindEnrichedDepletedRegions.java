package broad.pda.rap.misc;

import broad.core.parser.CommandLineParser;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;

import org.apache.log4j.Logger;

import broad.pda.annotation.BEDFileParser;
import nextgen.core.annotation.Annotation;
import nextgen.core.coordinatesystem.GenomicSpace;
import nextgen.core.coordinatesystem.MaskedGenomicSpace;
import nextgen.core.feature.Window;
import nextgen.core.model.ScanStatisticDataAlignmentModel;
import nextgen.core.model.score.ScanStatisticScore;

public class FindEnrichedDepletedRegions {

	Logger logger=Logger.getLogger(FindEnrichedDepletedRegions.class.getName());
	int windowLength;
	
	/**
	 * Try and find large regions of enrichment and depletion that are unexpected given a  randomized model on a given chromosome
	 * @param model Data alignment model representing the underling read data
	 * @throws IOException 
	 */
	public FindEnrichedDepletedRegions(ScanStatisticDataAlignmentModel model, int windowLength, String chr, String save) throws IOException{
		this.windowLength=windowLength;
		
		//TODO We should get the local lambda
		
		//Idea 1: Look for enriched windows based on local lambda of chromosome (how to deal with depletion?)
		Iterator<ScanStatisticScore> scores=scanWindows(model, chr);
		write(save, scores, model);
		
		//for each point consider the left X bases and right X bases and determine the average difference
		//Collection<Annotation> changePoints=changePoints(model);
		//determine if difference is significant using parametric test
		
	}
	
	private void write(String save, Iterator<ScanStatisticScore> scores, ScanStatisticDataAlignmentModel model) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		while(scores.hasNext()){
			//logger.info("iterating");
			//TODO Use local lambda by chromosome
			ScanStatisticScore score=scores.next();
			Annotation feature=score.getAnnotation();
			feature.setName("e="+score.getEnrichment(model));
			logger.info(score.getAnnotation().getChr()+"\t"+score.getAnnotation().getStart()+"\t"+score.getAnnotation().getEnd()+"\t"+score.getScanPvalue());
			if(score.getScanPvalue()<0.05){writer.write(feature.toBED()+"\n");}
		}
		
		writer.close();
	}

	private Iterator<ScanStatisticScore> scanWindows(ScanStatisticDataAlignmentModel model, String chr) throws IOException {
		throw new UnsupportedOperationException("TODO");
		//Iterator<ScanStatisticScore> scores=model.scan(chr, this.windowLength, this.windowLength-1, true);
		//return scores;
	}


	public static void main(String[] args) throws IOException{
		CommandLineParser p = new CommandLineParser();
		p.addIntArg("-w", "window size", false, 500);
		p.addStringArg("-d1", "data bam file", true);
		p.addStringArg("-s", "chr size file", true);
		p.addStringArg("-c", "chromosome", false, "chrX");
		p.addStringArg("-o", "output file", true);
		p.addStringArg("-m", "mask file", false);
		p.parse(args);
		
		
		//MaskedGenomicSpace space=new MaskedGenomicSpace(BEDFileParser.loadChrSizes(p.getStringArg("-s")), BEDFileParser.loadAlignmentDataByChr(new File(p.getStringArg("-m"))), 500, false);
		
		GenomicSpace space=new GenomicSpace(BEDFileParser.loadChrSizes(p.getStringArg("-s")));
		
		
		System.err.println("Genome length= "+space.getLength());
		
		ScanStatisticDataAlignmentModel data=new ScanStatisticDataAlignmentModel(p.getStringArg("-d1"), space);
		int windowSize= p.getIntArg("-w");
		new FindEnrichedDepletedRegions(data, windowSize, p.getStringArg("-c"), p.getStringArg("-o"));
	}
	
}
