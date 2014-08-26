package broad.pda.seq.rap;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.io.File;
import java.io.IOException;
import java.util.List;

import broad.pda.annotation.BEDFileParser;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.Log;

import nextgen.core.coordinatesystem.CoordinateSpace;
import nextgen.core.coordinatesystem.GenomicSpace;
import nextgen.core.coordinatesystem.TranscriptomeSpace;
import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.AnnotationList;
import nextgen.core.annotation.BasicAnnotation;
import nextgen.core.model.AlignmentModel;
import nextgen.core.readFilters.*;

/**
 * @author engreitz
 * Class used to run various analyses on sequencing data alignments in genomic or transcriptome space.
 * Automatically filters improper pairs, chimeric reads and PCR duplicates
 * 
 */
public abstract class GenomeCommandLineProgram extends CommandLineProgram {
    private static final Log log = Log.getInstance(GenomeCommandLineProgram.class);
    
	@Option(doc="File specifying chromosome sizes.")
	public String SIZES = "/seq/lincRNA/data/mm9/sizes";
	
	@Option(doc="File containing masked regions.", optional=true)
	public String MASK_FILE = "/seq/mguttman/ChIPData/MaskFiles/MM9Segments/all.mask.mouse.n36.d2.bin";
	
	@Option(doc="Percent masked allowable per sliding window", optional=true)
	public double PCT_MASKED_ALLOWED = 50.0;
	
	@Option(doc="Region to process (e.g. chr1, chr1:5000-50230) Default (null) processes the entire genome", optional=true)
	public String REGION = null;
	
	@Option(doc="Specify region by gene name (from ANNOTATION file)", optional=true)
	public String GENE = null;
	
	@Option(doc="Maximum paired-end read fragment length to consider.", optional=true)
	public Integer MAX_FRAGMENT_LENGTH = 10000;
	
	@Option(doc="Minimum mapping quality for reads") 
	public Integer MIN_MAPPING_QUALITY = 30;
	
	@Option(doc="Specify genomic space or transcriptome space")
	public String COORD_SPACE = "genomic";
	
	@Option(doc="File containing gene annotations (BED)")
	public String ANNOTATION = "/seq/lincRNA/Shari/Annotations/RefSeq_LincV3_ChromatinWithNames_nonrandom_collapsed_overlappers.bed";
	
	protected CoordinateSpace coordinateSpace;

	
	
	@Override
	protected String[] customCommandLineValidation() {
		try {
			loadCoordinateSpace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return super.customCommandLineValidation();
	}
	
	protected void loadCoordinateSpace() throws IOException {
		
		if (COORD_SPACE.equalsIgnoreCase("genomic")){
			coordinateSpace = new GenomicSpace(SIZES, MASK_FILE, PCT_MASKED_ALLOWED);
			}
		else if (COORD_SPACE.equalsIgnoreCase("transcriptome")){
			coordinateSpace = new TranscriptomeSpace(BEDFileParser.loadDataByChr(new File(ANNOTATION)));
			}
		}
	
	public List<Annotation> getRegions() {
		List<Annotation> regions = new ArrayList<Annotation>();

		if (REGION != null) {
			
			if (coordinateSpace.hasChromosome(REGION)) {
				regions.add(coordinateSpace.getReferenceAnnotation(REGION));
			} else {
				try {
					regions.add(new BasicAnnotation(REGION));
				} catch (RuntimeException e) {
					throw new IllegalArgumentException("REGION is improperly formatted");
				}
			}
		} else if (GENE != null) {
			regions.add(getGenes());
		}
		else {
			regions.addAll(coordinateSpace.getReferenceAnnotations());
		}
		
		return regions;
	}
	
	public Map<String,Annotation> getRegionsMap() {
		Map<String,Annotation> regions = new HashMap<String,Annotation>();

		if (REGION != null) {
			if (coordinateSpace.hasChromosome(REGION)) {
				regions.put(REGION,coordinateSpace.getReferenceAnnotation(REGION));
			} else {
				try {
					regions.put(REGION,new BasicAnnotation(REGION));
				} catch (RuntimeException e) {
					throw new IllegalArgumentException("REGION is improperly formatted");
				}
			}
		} else {
			for(String chr:coordinateSpace.getReferenceNames()){
				regions.put(chr,coordinateSpace.getReferenceAnnotation(chr));
			}
		}
		return regions;
	}
	
	public AnnotationList<Annotation> getRegionSet() {
		AnnotationList<Annotation> annotations = new AnnotationList<Annotation>(coordinateSpace, getRegions());
		return annotations;
	}

	public Annotation getGenes() {
		
		if (GENE != null){
			try {
				return coordinateSpace.getReferenceAnnotation(GENE);
			} catch (RuntimeException e) {
				throw new IllegalArgumentException("Gene not in transcriptome");
			}
		} else {
			return null;
		}
	}
	
	protected CoordinateSpace getCoordinateSpace() { return coordinateSpace; }
	
	public Boolean isGenomicSpace() {
		if (COORD_SPACE.equalsIgnoreCase("genomic")) {
			return true;
		} else {
			return false;
		}
	}
	
	public AlignmentModel loadAlignmentModel(File bamFile) {
		return loadAlignmentModel(bamFile, true);
	}

	public AlignmentModel loadAlignmentModel(File bamFile, boolean pairedEnd) {
		// TODO: Extend to handle BED or other annotation files
		
		IoUtil.assertFileIsReadable(bamFile);
		AlignmentModel model = new AlignmentModel(bamFile.getAbsolutePath(), coordinateSpace, pairedEnd);
		if (COORD_SPACE.equalsIgnoreCase("genomic")) {
			model.addFilter(new GenomicSpanFilter(MAX_FRAGMENT_LENGTH));
		} else if (COORD_SPACE.equalsIgnoreCase("transcriptome")) {
			model.addFilter(new FragmentLengthFilter(coordinateSpace,MAX_FRAGMENT_LENGTH));
		}
		model.addFilter(new ChimeraFilter());
		model.addFilter(new DuplicateFilter());
		if (pairedEnd) model.addFilter(new ProperPairFilter());
		model.addFilter(new MappingQualityFilter(MIN_MAPPING_QUALITY));

		
		// TODO need to modify PairedEndWriter to save information about the other read
		return model;
	}
}
