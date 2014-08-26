package broad.pda.samtools;

import net.sf.samtools.*;
import net.sf.samtools.util.CloseableIterator;
import broad.core.datastructures.Pair;


/**
 * @author engreitz
 * This class iterators over a standard format SAM or BAM file, but gives an iterator over paired reads
 * with the same name.  Does not return reads that do not have a mate.
 */
public class PairedEndReadIterator implements CloseableIterator<Pair<SAMRecord>>{
	private SAMRecordIterator itr;
	private SAMFileHeader header;
	private Pair<SAMRecord> nextPair = null;
    private boolean isClosed = false;
    
    private int singleReadsSkipped = 0;
	
	public PairedEndReadIterator(SAMFileReader reader) {
		itr = reader.iterator();
		header = reader.getFileHeader();
		if (header.getSortOrder() != SAMFileHeader.SortOrder.queryname) {
			throw new IllegalArgumentException("At the moment, PairedEndReadIterator only handles SAM files sorted in queryname (read name) order.");
		}
		advance();
	}
	
	
	private void advance() {
		if (!itr.hasNext()) {
			nextPair = null;
			return;
		}
		
		SAMRecord first = itr.next();
		SAMRecord second;
		while (true) {
			if (!itr.hasNext()) {
				singleReadsSkipped++;
				nextPair = null;
				return;
			}
			
			second = itr.next();
			if (SAMRecordQueryNameComparator.compareReadNames(first.getReadName(), second.getReadName()) == 0) {
				break;
			} else {
				first = second;
				singleReadsSkipped++;
			}
		}

		nextPair = new Pair<SAMRecord>();
		nextPair.setValue1(first);
		nextPair.setValue2(second);
	}
	
	
	@Override
	public boolean hasNext() {
		if (isClosed) throw new IllegalStateException("Iterator has been closed");
		return nextPair != null;
	}
	
	@Override
	public Pair<SAMRecord> next() {
		if (isClosed) throw new IllegalStateException("Iterator has been closed");
		Pair<SAMRecord> result = nextPair;
		advance();
		return result;
	}
	
	@Override
	public void close() {
		if (!isClosed) {
			if (itr == null) {
				throw new IllegalStateException("Attempt to close non-current iterator");
			}
			itr.close();
			itr = null;
			isClosed = true;
		}
	}


	@Override
	public void remove() {
		throw new UnsupportedOperationException("Not supported: remove");
	}


	public int getSingleReadsSkipped() {
		return singleReadsSkipped;
	}
	
}
