package broad.core.sequence;

import java.io.IOException;


public interface FastaQualityHandler extends FastaHandler{

	public void newQuality(FastaQualityParser parser) throws IOException;

}
