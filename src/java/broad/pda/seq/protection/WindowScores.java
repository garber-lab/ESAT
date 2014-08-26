package broad.pda.seq.protection;

import broad.core.parser.StringParser;

import java.io.IOException;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.apache.log4j.Logger;

import broad.core.annotation.MaximumContiguousSubsequence;
import broad.core.math.Statistics;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Gene;
import nextgen.core.annotation.Annotation.Strand;
import nextgen.core.coordinatesystem.TranscriptomeSpace;
import nextgen.core.model.TranscriptomeSpaceAlignmentModel;
import nextgen.core.model.score.CountScore;
import nextgen.core.model.score.CountScore;
import nextgen.core.model.score.ScanStatisticScore;
import nextgen.core.model.score.WindowProcessor;
import nextgen.core.model.score.WindowScoreIterator;
import nextgen.core.readFilters.FragmentLengthFilter;
import nextgen.core.readFilters.GenomicSpanFilter;
import nextgen.core.readFilters.MappingQualityFilter;
import nextgen.core.readFilters.NumHitsFilter;
import broad.pda.seq.protection.SampleData;

public class WindowScores {

	protected SampleData sample;
	protected SampleData ctrl;
	protected Map<Gene, CountScore> geneScores;
	protected Map<Gene, Double> geneAvgCoverage;
	protected Map<Gene, Map<Annotation, CountScore>> windowScores;
	protected int windowSize;
	protected int stepSize;
	protected static Logger logger = Logger.getLogger(SampleData.class.getName());
	private WindowProcessor<CountScore> processor;
	protected Map<String, Collection<Gene>> genesByChr;
	protected Map<String, Gene> genesByName;
	protected double expressionCutoffValue;
	private boolean gotWindowScoresFromFile;
	private static int DEFAULT_MAX_GENOMIC_SPAN = 100000;
	private static int DEFAULT_MAX_FRAGMENT_LENGTH = 150;
	protected boolean expressionByPval;
	private String originalBamFile;
	private boolean read1TranscriptionStrand;
	protected boolean fullyContainedReads;
	
	
	
}
