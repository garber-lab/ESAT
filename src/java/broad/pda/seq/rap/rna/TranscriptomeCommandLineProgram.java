package broad.pda.seq.rap.rna;

import java.io.File;
import java.io.IOException;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.io.IoUtil;
import net.sf.picard.util.Log;

import nextgen.core.alignment.AbstractPairedEndAlignment.TranscriptionRead;
import nextgen.core.coordinatesystem.*;
import nextgen.core.model.AlignmentModel;
import nextgen.core.readFilters.ChimeraFilter;
import nextgen.core.readFilters.GenomicSpanFilter;
import nextgen.core.readFilters.MappingQualityFilter;
import nextgen.core.readFilters.ProperPairFilter;
import broad.pda.annotation.BEDFileParser;
import broad.pda.seq.rap.GenomeCommandLineProgram;


public abstract class TranscriptomeCommandLineProgram extends GenomeCommandLineProgram {
	private static final Log log = Log.getInstance(TranscriptomeCommandLineProgram.class);
	
	@Option(doc="BED file containing genes", shortName="GENES")
	public File GENE_BED;
	
	@Option(doc="Whether to consider strand, and which read is the stranded.  ")
	public TranscriptionRead TRANSCRIPTION_READ=TranscriptionRead.SECOND_OF_PAIR;
	
	protected TranscriptomeSpace ts;
	
	
	@Override
	protected String[] customCommandLineValidation() {
		loadCoordinateSpace();
		return super.customCommandLineValidation();
	}
	
	@Override
	protected void loadCoordinateSpace() {
		try {
			super.loadCoordinateSpace();
			ts = new TranscriptomeSpace(BEDFileParser.loadDataByChr(GENE_BED));
		} catch (IOException e) {
			throw new RuntimeException("Encountered IOException " + e);
		}
	}
	
	protected GenomicSpace getGenomicSpace() {
		return (GenomicSpace) getCoordinateSpace();
	}
	
	protected TranscriptomeSpace getTranscriptomeSpace() {
		return ts;
	}
	
	
	public AlignmentModel loadRNAAlignmentModel(File bamFile, boolean pairedEnd) {
		IoUtil.assertFileIsReadable(bamFile);
		AlignmentModel model = new AlignmentModel(bamFile.getAbsolutePath(), ts, pairedEnd, TRANSCRIPTION_READ);
		return model;
	}

	@Override
	public AlignmentModel loadAlignmentModel(File bamFile, boolean pairedEnd) {
		IoUtil.assertFileIsReadable(bamFile);
		AlignmentModel model = new AlignmentModel(bamFile.getAbsolutePath(), coordinateSpace, pairedEnd, TRANSCRIPTION_READ);
		model.addFilter(new GenomicSpanFilter(MAX_FRAGMENT_LENGTH));
		model.addFilter(new ChimeraFilter());
		model.addFilter(new ProperPairFilter());
		model.addFilter(new MappingQualityFilter(MIN_MAPPING_QUALITY));
		return model;
	}
}
