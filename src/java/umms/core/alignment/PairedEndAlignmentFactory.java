package umms.core.alignment;

import umms.core.alignment.AbstractPairedEndAlignment.TranscriptionRead;

public class PairedEndAlignmentFactory {

	
	public AbstractPairedEndAlignment getAlignment(boolean fragment,SingleEndAlignment firstMate,
														SingleEndAlignment secondMate, TranscriptionRead strand){
		if(!fragment){
			return new PairedReadAlignment(firstMate, secondMate,strand);
		}
		else{
			return new FragmentAlignment(firstMate, secondMate,strand);
		}
	}
}
