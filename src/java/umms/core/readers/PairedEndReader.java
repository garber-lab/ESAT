package umms.core.readers;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import org.apache.log4j.Logger;
import org.broad.igv.Globals;

import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;


import net.sf.picard.util.Log;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.util.CloseableIterator;
import umms.core.alignment.Alignment;
import umms.core.alignment.ChromosomeInconsistencyException;
import umms.core.alignment.AbstractPairedEndAlignment.TranscriptionRead;
import umms.core.alignment.PairedEndAlignmentFactory;
import umms.core.alignment.SingleEndAlignment;
import umms.core.annotation.Annotation;
import umms.core.writers.PairedEndWriter;


//Reads the new PairedEnd BAM files into a queryable object

public class PairedEndReader {
    private static final Log log = Log.getInstance(PairedEndReader.class);
	
	public static enum AlignmentType { PAIRED_END, SINGLE_END };
    private boolean warned = false;
    private static Logger logger = Logger.getLogger(PairedEndIterator.class.getName());

	private CloseableIterator<Alignment> mCurrentIterator = null;
	private SAMFileReader reader;
	private SAMFileHeader header;
	
	private AlignmentType alignmentType;
	private TranscriptionRead strand;
	private boolean fragment=true;
	
	public PairedEndReader(File bam){
		//We are going to read the .pebam files with special attribute flags and create a new "Alignment" object
		this(bam,TranscriptionRead.UNSTRANDED);
	}

	public PairedEndReader(File bam,TranscriptionRead read){
		//We are going to read the .pebam files with special attribute flags and create a new "Alignment" object
		this(bam,read,true);
	}

	public PairedEndReader(File bam,TranscriptionRead read,boolean fra){
		//We are going to read the .pebam files with special attribute flags and create a new "Alignment" object
		this.reader=new SAMFileReader(bam);
		this.header=reader.getFileHeader();
		alignmentType = getAlignmentType(this.header);
		strand = read;
		if(alignmentType==AlignmentType.PAIRED_END){
			logger.info("Alignment is Paired");
		}
		else{
			logger.info("Alignment is Single");
		}
		
		setFragmentFlag(fra);
	}
	
	public static AlignmentType getAlignmentType(SAMFileHeader header) {
		//check if it is paired end with our modified record
		Object isPairedEndFormat=header.getAttribute(PairedEndWriter.mateLineFlag);
		return (isPairedEndFormat != null) ? AlignmentType.PAIRED_END : AlignmentType.SINGLE_END;
	}
	
	public AlignmentType getAlignmentType() { 
		return getAlignmentType(this.header);
	}
	
	public void setFragmentFlag(boolean fra){
		fragment = fra;
	}
	
	public void close() throws IOException {
		mCurrentIterator.close();
		reader.close();
	}
	
	public SAMFileHeader getHeader() {
		 if (header == null) {
	            loadHeader();
	        }
	        return header;
	}
	
	 private void loadHeader() {
	        header = reader.getFileHeader();
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


	public boolean hasIndex() {
        return reader.hasIndex();
    }

	public CloseableIterator<Alignment> iterator() {
		if (mCurrentIterator != null) {
            throw new IllegalStateException("Iteration in progress");
        }
		mCurrentIterator = new PairedEndIterator();
		return mCurrentIterator;
	}

	public CloseableIterator<Alignment> query(Annotation a, boolean contained) {
		SAMRecordIterator query = null;
	    query = reader.query(a.getReferenceName(), a.getSAMStart(), a.getSAMEnd(), contained);
	    mCurrentIterator = new PairedEndIterator(query);
	    return mCurrentIterator;	     
	}
	 	 
	
	private class PairedEndIterator implements CloseableIterator<Alignment> {
        private Alignment mNextRecord = null;
        private boolean isClosed = false;
        private CloseableIterator<SAMRecord> itr;
        
        public PairedEndIterator() {
        	this(reader.iterator());
        }
        
        public PairedEndIterator(CloseableIterator<SAMRecord> itr) {
        	this.itr = itr;
        	advance();
        }
        
		public void close() {
			if (!isClosed) {
				if (mCurrentIterator == null) {
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

		public Alignment next() {
			if (isClosed) throw new IllegalStateException("Iterator has been closed");
			final Alignment result = mNextRecord;
			advance();
			return result;
		}
		
		private void advance() {
			if (!itr.hasNext()) mNextRecord = null;  // VERY IMPORTANT.  Otherwise infinite loop
			while (itr.hasNext()) {
				SAMRecord r = itr.next();
				mNextRecord = samRecordToAlignment(r,strand,fragment);
				if (mNextRecord != null) {break;} 
				else {log.debug("samRecordToAlignment returned null for this record" + r.getSAMString()  );}
			}
        }
		
        /**
         * @return The record that will be return by the next call to next()
         */
        protected Alignment peek() {
            return mNextRecord;
        }
	}

	/**
	 * Parse the record
	 * Either it is our modified record with the mate information or a standard single end read
	 * @param record SamRecord
	 * @return Alignment
	 * @throws ChromosomeInconsistencyException 
	 */
	private Alignment samRecordToAlignment(SAMRecord record,TranscriptionRead transcriptionRead,boolean fragment) {
		Alignment rtrn;
		
		try {
			if (alignmentType == AlignmentType.PAIRED_END) {
				if (record.getReadPairedFlag() && !record.getMateUnmappedFlag()) {   //revised to read single end data @zhuxp 
					String name=record.getReadName();

					int mateStart=record.getMateAlignmentStart();
					String mateCigar=record.getAttribute(PairedEndWriter.mateCigarFlag).toString();
					String mateSequence=record.getAttribute(PairedEndWriter.mateSequenceFlag).toString();
					boolean mateNegativeStrand=record.getMateNegativeStrandFlag();
					String mateChr=record.getMateReferenceName();

					int readStart=new Integer(record.getAttribute(PairedEndWriter.readStartFlag).toString());
					String readCigar=record.getAttribute(PairedEndWriter.readCigarFlag).toString();
					String readSequence=record.getReadString();
					boolean readNegativeStrand=record.getReadNegativeStrandFlag();
					String readChr=record.getReferenceName();

					// Make the first mate
					record.setAlignmentStart(readStart);
					record.setReferenceName(readChr);
					record.setReadNegativeStrandFlag(readNegativeStrand);
					record.setCigarString(readCigar);
					record.setReadString(readSequence);
					record.setReadName(name);

					// Make the second mate
					SAMRecord record2=new SAMRecord(record.getHeader());
					record2.setAlignmentStart(mateStart);
					record2.setReferenceName(mateChr);
					record2.setReadNegativeStrandFlag(mateNegativeStrand);
					record2.setCigarString(mateCigar);
					record2.setReadString(mateSequence);
					record2.setReadName(name);

					SingleEndAlignment firstMate=new SingleEndAlignment(record);
					SingleEndAlignment secondMate=new SingleEndAlignment(record2);

					rtrn=new PairedEndAlignmentFactory().getAlignment(fragment, firstMate, secondMate, transcriptionRead);
					rtrn.setProperPairFlag(record.getProperPairFlag());
				
				} else {
					rtrn = new SingleEndAlignment(record,record.getFirstOfPairFlag());
				}
			} else {
				//check if paired end without our modified record --> throw warning
				if (record.getReadPairedFlag() && !record.getMateUnmappedFlag() && !warned) {
					log.warn("DEBUG: Paired end reads were found but file is not in our paired end format.  Processing as single-end reads ..."); 
					warned = true;
				}				
				//If first read is in direction of trasncription
				if(transcriptionRead.equals(TranscriptionRead.FIRST_OF_PAIR)){
					//This is the first read
					if(!record.getReadPairedFlag() || record.getFirstOfPairFlag()){
						//Orientation of fragment is same as that of read
					}
					//This is the other mate
					//Reverse its orientation
					else{
						record.setReadNegativeStrandFlag(!record.getReadNegativeStrandFlag());
					}
				}
				//Second read is the transcription read
				else if(transcriptionRead.equals(TranscriptionRead.SECOND_OF_PAIR)){
					//This is the first read
					//Reverse orientation
					if(record.getReadPairedFlag() && record.getFirstOfPairFlag()){
						record.setReadNegativeStrandFlag(!record.getReadNegativeStrandFlag());
					}
					//This is the other mate
					else{
						//NOTHING
					}
				}//UNSTRANDED
				else{
					//NOTHING
				}

				rtrn = new SingleEndAlignment(record, record.getReadPairedFlag() && record.getFirstOfPairFlag());
				//rtrn=new SingleEndAlignment(record);
			}
		} catch (RuntimeException e) {
			log.error("Failed on SAMRecord: " + record.toString());
			throw e;
		}

		return rtrn;
	}
	

	/**
	 * Get the lengths of the reference sequences
	 * @param chr the chromsoome to query size
	 * @return Map associating each reference name with sequence length
	 */
	public int getRefSequenceLengths(String chr) {
		Set<String> seqNames = getSequenceNames();
		for(String seq : seqNames) {
			SAMSequenceRecord rec = header.getSequence(seq);
			if(seq.equalsIgnoreCase(chr)){
				return Integer.valueOf(rec.getSequenceLength());
			}
		}
		return -99;
	}
	
	/**
	 * Get the lengths of the reference sequences
	 * @return Map associating each reference name with sequence length
	 */
	public Map<String, Integer> getRefSequenceLengths() {
		Map<String, Integer> rtrn=new TreeMap<String, Integer>();
		Set<String> seqNames = getSequenceNames();
		for(String seq : seqNames) {
			SAMSequenceRecord rec = header.getSequence(seq);
			rtrn.put(seq, Integer.valueOf(rec.getSequenceLength()));
			
		}
		return rtrn;
	}

	
	/**
	 * Returns {@code true} if the given SAM file is in paired end format (has the mateLine) attribute.
	 * Will also check the file with the default extension.
	 * @param fileToCheck
	 * @return
	 */
	public static boolean isPairedEndFormat(String fileToCheck) {
		return (getPairedEndFile(fileToCheck) != null);
	}
	
	
	/**
	 * Returns the file itself or the file plus the default {@code PAIRED_END_EXTENSION} if either
	 * is in paired end format.  If neither is in paired end format, returns null.
	 * @param fileToCheck
	 * @return
	 */
	public static String getPairedEndFile(String fileToCheck) {
		//Check the header to see if it is "Paired End Format"
		SAMFileReader reader = new SAMFileReader(new File(fileToCheck));
		if (getAlignmentType(reader.getFileHeader()) == AlignmentType.PAIRED_END) {
			return fileToCheck;
		}
				
		//else, lets make the logical extension and see if it is there
		String tmpBam = fileToCheck + PairedEndWriter.PAIRED_END_EXTENSION;
		boolean fileExists = new File(tmpBam).exists();
		
		if (fileExists) {
			//Lets test to make sure it has the correct line
			reader = new SAMFileReader(new File(tmpBam));
			if (getAlignmentType(reader.getFileHeader()) == AlignmentType.PAIRED_END) {
				return tmpBam;
			}
		}
		
		return null;
	}
	
	/**
	 * Finds a paired end file or creates one if it doesn't exist
	 * @param fileToCheck
	 * @return
	 */
	public static String getOrCreatePairedEndFile(String fileToCheck,TranscriptionRead txnRead) {
		String result = PairedEndReader.getPairedEndFile(fileToCheck);
		if (result == null) {
			result = PairedEndWriter.getDefaultFile(fileToCheck);
			PairedEndWriter writer = new PairedEndWriter(new File(fileToCheck), result);
			writer.convertInputToPairedEnd(txnRead);
		}
		return result;
	}
	
	public static void main(String[] args){
		Globals.setHeadless(true);
		/*
		 * @param for ArgumentMap - size, usage, default task
		 * argMap maps the command line arguments to the respective parameters
		 */
		ArgumentMap argMap = CLUtil.getParameters(args,usage,"makeFile");
		
		TranscriptionRead strand = TranscriptionRead.UNSTRANDED;
		if(argMap.get("strand").equalsIgnoreCase("first")){
			//System.out.println("First read");
			strand = TranscriptionRead.FIRST_OF_PAIR;
		}
		else if(argMap.get("strand").equalsIgnoreCase("second")){
			//System.out.println("Second read");
			strand = TranscriptionRead.SECOND_OF_PAIR;
		}
		else
			log.info("no strand");
		//String bamFile = argMap.getMandatory("alignment");
		String bamFile = getOrCreatePairedEndFile(argMap.getMandatory("alignment"),strand);
    	
		String file = PairedEndReader.getPairedEndFile(bamFile);
		if (file == null) {
			file = PairedEndWriter.getDefaultFile(bamFile);
			PairedEndWriter writer = new PairedEndWriter(new File(bamFile), file);
			writer.convertInputToPairedEnd(strand);
		}
	}
	
	static final String usage = "Usage: CreatePairedEndBamFile -task makeFile "+
			"\n**************************************************************"+
			"\n\t\tArguments"+
			"\n**************************************************************"+
			"\n\t\t-strand <VALUES: first, second, unstranded. Specifies the mate that is in the direction of transcription DEFAULT: Unstranded> "+
			"\n\n\t\t-alignment <Alignment file to be used for reconstruction. Required.> ";			
	
}
