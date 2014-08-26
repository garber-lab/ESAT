package broad.pda.seq.rap;

import java.io.*;
import java.util.*;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.util.Log;
import net.sf.samtools.*;

public class CalculateCoverageStats extends CommandLineProgram {
    private static final Log log = Log.getInstance(CommandLineProgram.class);

    @Option(doc="Input SAM or BAM file", shortName="I")
    public File INPUT;

	@Option(doc="Output file", shortName="O")
	public File OUTPUT;

	/**
	 * Stock main method.
	 *
	 * @param args main arguments
	 */
	public static void main(final String[] args) {
		System.exit(new CalculateCoverageStats().instanceMain(args));
	}
	
	

	@Override
	protected int doWork() {
		
		try {
			if (OUTPUT.exists()) {
				OUTPUT.delete();
			}
			
        	SAMFileReader reader = new SAMFileReader(INPUT);
        	final SAMFileHeader header = reader.getFileHeader();

			List<String> refs = new ArrayList<String>();
			Map<String, Integer> lengths = new HashMap<String, Integer>();
			Map<String, Integer> counts = new HashMap<String, Integer>();
			for (SAMSequenceRecord ref : header.getSequenceDictionary().getSequences()) {
				refs.add(ref.getSequenceName());
				lengths.put(ref.getSequenceName(), ref.getSequenceLength());
				counts.put(ref.getSequenceName(), 0);
			}
			
			Iterator<SAMRecord> itr = reader.iterator();
			while (itr.hasNext()) {
				SAMRecord rec = itr.next();
				if (counts.containsKey(rec.getReferenceName())) {
					counts.put(rec.getReferenceName(), counts.get(rec.getReferenceName()) + 1);
				}
			}
			
			BufferedWriter bw = new BufferedWriter(new FileWriter(OUTPUT));
			for (String ref : refs) {
				bw.write(ref + "\t" + lengths.get(ref) + "\t" + counts.get(ref) + "\n");
			}
			
			reader.close();
			bw.close();
        	
		} catch (Exception e) {
			log.error(e);
			return 1;
		}
		
		return 0;
	}
	
}
