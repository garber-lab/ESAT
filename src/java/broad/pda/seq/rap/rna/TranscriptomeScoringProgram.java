package broad.pda.seq.rap.rna;

import java.io.File;

import net.sf.picard.cmdline.Option;
import net.sf.picard.util.Log;
import nextgen.core.model.AlignmentModel;
import nextgen.core.model.score.*;
import nextgen.core.annotation.*;
import nextgen.core.alignment.AbstractPairedEndAlignment.TranscriptionRead;

public abstract class TranscriptomeScoringProgram extends TranscriptomeCommandLineProgram {
	
	@Option(doc="Input SAM or BAM file", shortName="I")
	public File TARGET;

	@Option(doc="Control SAM or BAM file for normalization (use with SCORE=ratio).", optional=true)
	public File CONTROL = null;
	
	@Option(doc="Scoring function to use")
	public String SCORE = "count";
	
	@Option(doc="Whether to force paired end behavior")
	public boolean PAIRED_END=true;
	
	
	@Override
	protected String[] customCommandLineValidation() {
		if (SCORE.equalsIgnoreCase("ratio") && CONTROL == null) {
			return new String[] { "Must specify CONTROL when using SCORE=ratio." };
		}
		return super.customCommandLineValidation();
	}
	
	
	public WindowProcessor<? extends WindowScore> getWindowProcessor() {
		AnnotationCollection<? extends Annotation> target = loadRNAAlignmentModel(TARGET, PAIRED_END);
		AnnotationCollection<? extends Annotation> control = CONTROL == null ? null : loadRNAAlignmentModel(CONTROL, PAIRED_END);
		return getWindowProcessor(target, control);
	}

	
	public WindowProcessor<? extends WindowScore> getWindowProcessor(AnnotationCollection<? extends Annotation> target, AnnotationCollection<? extends Annotation> control) {
		WindowProcessor<? extends WindowScore> processor;
		
		if (SCORE.equalsIgnoreCase("count")) {
			processor = new CountScore.Processor(target);
		} else if (SCORE.equalsIgnoreCase("ratio")) {
			processor = new RatioScore.Processor(target, control);
		} else if (SCORE.equalsIgnoreCase("coverage")) {
			processor = new CoverageScore.Processor(target);
		} else if (SCORE.equalsIgnoreCase("length")) {
			processor = new LengthScore.Processor(target);
		} else if (SCORE.equalsIgnoreCase("sum")) {
			processor = new SumScore.Processor(target);
		} else if(SCORE.equalsIgnoreCase("all")){
			processor = new WindowAllScore.Processor((AlignmentModel)target, (AlignmentModel)control);
		} else {
			throw new IllegalArgumentException("Could not find scoring class for " + SCORE);
		}
		
		return processor;
	}
}
