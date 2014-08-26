package umms.core.readers;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecordIterator;

import net.sf.samtools.util.CloseableIterator;

import umms.core.readers.SCSRecord;

/**
 * A general interface that adds functionality to a CloseableIterator of
 * SAMRecords.  Currently, this interface is implemented by iterators that
 * want to validate as they are iterating that that the records in the
 * underlying SAM/BAM file are in a particular order.
 */
public interface SCSRecordIterator extends SAMRecordIterator {

    /**
     * Establishes that records returned by this iterator are expected to
     * be in the specified sort order.  If this method has been called,
     * then implementers must throw an IllegalStateException from next()
     * when a record is read that violates the sort order.  This method
     * may be called multiple times over the course of an iteration,
     * changing the expected sort, if desired -- from the time it is called,
     * it validates whatever sort is set, or stops validating if it
     * is set to null or SAMFileHeader.SortOrder.unsorted.  If this method
     * is not called, then no validation of the iterated records is done.
     *
     * @param sortOrder The order in which records are expected to be returned
     * @return  This SAMRecordIterator
     */
    public SCSRecordIterator assertSorted(SAMFileHeader.SortOrder sortOrder);
    
}
