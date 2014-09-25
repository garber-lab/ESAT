package umms.core.readers;

import java.io.File;
import java.io.IOException;

import net.sf.samtools.*;
import umms.core.readers.SCSRecordIterator;
import umms.core.readers.SCSRecord;

import org.apache.log4j.Logger;

public class XamReaderTest {
	
		public static void main(String[] args) throws IOException {
		
			//SAMRecord xamRec = null;
			SCSRecord xamRec = null;
			
			/* sorted bam file with existing .bai file */
			String bamFile = "C:/Users/DerrA/Dropbox (UMASS MED - BIB)/islet_single_cell/DGE/shortened.20000000.sorted.wgu.bam";  

			/* sorted bam file with no index file */
//			String bamFile = "C:/Users/DerrA/Dropbox (UMASS MED - BIB)/islet_single_cell/DGE/shortened.20000000_noBai.sorted.wgu.bam";  

			/* sorted sam file */
//			String bamFile = "C:/Users/DerrA/Dropbox (UMASS MED - BIB)/islet_single_cell/DGE/shortened.20000000_noBai.sorted.wgu.sam";  

			File xamFile = new File(bamFile);
			if (xamFile.exists()) {
				System.out.println("File " + bamFile + " exists");
			}
			//SCSFileReader xamFileReader = new SCSFileReader(xamFile);
			
			//SCSRecordIterator xamIter = xamFileReader.iterator();
			
//			while (xamIter.hasNext()) {
//				xamRec = xamIter.next();
//				System.out.println("read name: " + xamRec.getReadName());
//				System.out.println(" read seq: " + xamRec.getReadString());
//			}
			
			//xamFileReader.close();

		}
}