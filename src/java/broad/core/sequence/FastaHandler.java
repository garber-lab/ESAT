package broad.core.sequence;

import java.io.IOException;


public interface FastaHandler {

	public void newSequence(AbstractFastaParser parser) throws IOException  ;
	public void newBase(AbstractFastaParser parser) throws IOException;
	public void eof(AbstractFastaParser parser) throws IOException ;

}
