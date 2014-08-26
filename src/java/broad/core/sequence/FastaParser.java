package broad.core.sequence;

import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;


public class FastaParser extends AbstractFastaParser {

	public void parse(InputStream in, FastaHandler handler) throws IOException {
		fileSize = 0;
		//RandomAccessFile raf = new RandomAccessFile(inputFastaFile, "r");
		//setRandomAccessFile(raf);
		//System.out.println("File size " + inputFastaFile.length());
		byte [] buf = new byte[BUFFER_SIZE];
		int i = 0;
		int read = 0;
		boolean readingSeqId = false;
		StringBuilder seqIdBuilder = new StringBuilder();
		boolean inComment = false;
		try {
			while(( read = in.read(buf)) > 0) {
				fileSize += read;
				for(int j = 0; j < read; j++) {
					char c = (char) buf[j];
					if('#' == c) {
						inComment = true;
					}
					//System.out.println("character " + c + " at file position " + raf.getFilePointer());
					if(inComment && (c != '\n' && c != '\r')) {
						continue;
					} else {
						inComment = false;
					}
					
					offset++;
					if(c == '>') {
						readingSeqId = true;
					} else if (readingSeqId) {
						if(c != '\n' & c != '\r') {
							seqIdBuilder.append(c);
						} else {
							setCurrentSequenceId(seqIdBuilder.toString());
							//System.out.println("new sequence " + getCurrentSequenceId());
							handler.newSequence(this);
							readingSeqId = false;
							seqIdBuilder = new StringBuilder(); //It used to be possible to reset a StringBuffer
						}
					} else {
						currentBase =  c;
						handler.newBase(this);
						i++;
					}
				}

			}
		} catch( EOFException eof) {
			handler.eof(this);
		} 
		handler.eof(this);
		
	}
	
}
