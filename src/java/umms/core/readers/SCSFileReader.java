package umms.core.readers;

import java.io.File;

import net.sf.samtools.SAMFileReader;
//import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMRecordIterator;

public class SCSFileReader extends SAMFileReader {
	
	public SCSFileReader(final File file) {
        super(file, null, false);
    }

	public SCSRecordIterator iterator() {
        return super.iterator();
    }
}
