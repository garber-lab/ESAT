/**
 * 
 */
package broad.pda.samtools;


import java.util.List;
import java.util.HashSet;
import java.util.Set;
import java.io.Closeable;
import java.io.File;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.util.CloseableIterator;

import org.apache.commons.lang3.StringUtils;
import org.broad.igv.sam.reader.AlignmentQueryReader;

import broad.pda.samtools.WrappedIterator;
import org.broad.igv.sam.Alignment;
import org.broad.igv.ui.util.MessageUtils;

/**
 * @author engreitz
 * Modeled after the Picard SAMFileReader
 * Currently this class does not just extend the Picard SAMFileReader because that class has a lot 
 * more functionality that would be a pain to implement in the paired-end format.
 * 
 * Also it's trying to accommodate iterators of both SAMRecords and Alignments to be compatible
 * with Scripture/IGV.  Annoying.
 * 
 * Note: With move to ScriptureV2, this class was replaced with nextgen.core.writers.PairedEndWriter 
 * and nextgen.core.readers.PairedEndReader
 */
public class SAMPairedEndFileReader implements AlignmentQueryReader {
	
	final static public int MAX_INSERT = 2000;			// perhaps make this a modifiable parameter
    private CloseableIterator<SAMRecord> mCurrentIterator = null;
	private SAMFileReader reader;
    SAMFileHeader header;
	
	public SAMPairedEndFileReader(final File file) {
		reader = new SAMFileReader(file);
	}
	
	public CloseableIterator<Alignment> iterator() {
		if (mCurrentIterator != null) {
            throw new IllegalStateException("Iteration in progress");
        }
		mCurrentIterator = new SAMPairedEndFileIterator();
		return new WrappedIterator(mCurrentIterator);
	}
	
	public CloseableIterator<SAMRecord> getSAMRecordIterator() {
		if (mCurrentIterator != null) {
            throw new IllegalStateException("Iteration in progress");
        }
		mCurrentIterator = new SAMPairedEndFileIterator();
		return mCurrentIterator;
	}
	
	public void close() {
		mCurrentIterator.close();
		reader.close();
	}
	
	public SAMFileReader getReader() {
		return reader;
	}
	
    public CloseableIterator<Alignment> query(String sequence, int start, int end, boolean contained) {
        SAMRecordIterator query = null;
        query = reader.query(sequence, start + 1, end, contained);
        return new WrappedIterator(new SAMPairedEndFileIterator(query));
    }
	
    public SAMFileHeader getHeader() {
        if (header == null) {
            loadHeader();
        }
        return header;
    }

    public Set<String> getSequenceNames() {
        SAMFileHeader header = getHeader();
        if (header == null) {
            return null;
        }
        Set<String> seqNames = new HashSet<String>();
        List<SAMSequenceRecord> records = header.getSequenceDictionary().getSequences();
        if (records.size() > 0) {
            for (SAMSequenceRecord rec : header.getSequenceDictionary().getSequences()) {
                String chr = rec.getSequenceName();
                seqNames.add(chr);
            }
        }
        return seqNames;
    }

    private void loadHeader() {
        header = reader.getFileHeader();
    }
    

    public boolean hasIndex() {
        return reader.hasIndex();
    }

	private class SAMPairedEndFileIterator implements CloseableIterator<SAMRecord> {
        private SAMRecord mNextRecord = null;
        private boolean isClosed = false;
        private CloseableIterator<SAMRecord> itr;
        
        public SAMPairedEndFileIterator() {
        	this(reader.iterator());
        }
        
        public SAMPairedEndFileIterator(CloseableIterator<SAMRecord> itr) {
        	this.itr = itr;
        	advance();
        }
        
		public void close() {
			if (!isClosed) {
				if (mCurrentIterator != null && this != mCurrentIterator) {
					throw new IllegalStateException("Attempt to close non-current iterator");
				}
				itr.close();
				mCurrentIterator = null;
				isClosed = true;
			}
		}
		
		public void remove() {
			throw new UnsupportedOperationException("Not supported: remove");
		}

		public boolean hasNext() {
			if (isClosed) throw new IllegalStateException("Iterator has been closed");
			return (mNextRecord != null);
		}

		public SAMRecord next() {
			if (isClosed) throw new IllegalStateException("Iterator has been closed");
			final SAMRecord result = mNextRecord;
			advance();
			return result;
		}
		
		private void advance() {
			if (!itr.hasNext()) mNextRecord = null;  // VERY IMPORTANT.  Otherwise infinite loop
			while (itr.hasNext()) {
				mNextRecord = convertToPairedEndFragment(itr.next());
				if (mNextRecord != null) break;
			}
        }
		
        /**
         * @return The record that will be return by the next call to next()
         */
        protected SAMRecord peek() {
            return mNextRecord;
        }
	}
	
	
	//Find the sequence of the other read and ame a real paired-end fragment
	public static SAMRecord convertToPairedEndFragment(SAMRecord rec) {
		int insertSize = Math.abs(rec.getInferredInsertSize());
		//TODO: possible to remove (insertSize > 0) clause? What does an insert size of 0 mean?
		if (rec.getReadPairedFlag() && rec.getProperPairFlag() && !rec.getReadUnmappedFlag() 
				&& (insertSize <= MAX_INSERT) && (rec.getAlignmentStart() < rec.getMateAlignmentStart())
				&& rec.getReferenceName() != "chrM") // current paired end representation doesn't do well with circular chromosomes
		{
			//We need to get the full fragment contained by read.getStart to pair.getEnd
			int readEnd=rec.getAlignmentEnd();
			int mateStart=rec.getMateAlignmentStart();
						
			//int extension = insertSize - rec.getReadLength();
			int extension=(mateStart-readEnd)+rec.getReadLength();
			if (extension <= MAX_INSERT) {
				String newRead = rec.getReadString() + StringUtils.repeat("N", extension);
				rec.setReadString(newRead);
				String newQual = StringUtils.repeat("A", newRead.length());
				if(!rec.getBaseQualityString().equals("*")) newQual = rec.getBaseQualityString() + StringUtils.repeat("A", extension);
				rec.setBaseQualityString(newQual);
				String newCigar=newRead.length()+"M";
				rec.setCigarString(newCigar);

				// Change attributes to represent single read
				rec.setMateReferenceName("*");
				rec.setMateAlignmentStart(0);
				rec.setFirstOfPairFlag(false);
				rec.setMateNegativeStrandFlag(false);
				rec.setMateUnmappedFlag(false);
				rec.setProperPairFlag(false);
				rec.setReadPairedFlag(false);
				rec.setSecondOfPairFlag(false);
			}
		} else {
			rec = null;
		}
		return rec;
	}
	


}
