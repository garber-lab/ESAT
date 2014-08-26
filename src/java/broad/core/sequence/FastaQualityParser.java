package broad.core.sequence;

import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;


public class FastaQualityParser extends AbstractFastaParser {

	int currentQual = 0;
	public void parse(InputStream in, FastaQualityHandler handler) throws IOException {
		fileSize = 0;
		//RandomAccessFile raf = new RandomAccessFile(inputFastaFile, "r");
		//setRandomAccessFile(raf);
		//System.out.println("File size " + inputFastaFile.length());
		byte [] buf = new byte[BUFFER_SIZE];
		int i = 0;
		int read = 0;
		boolean readingSeqId = false;
		StringBuilder seqIdBuilder = new StringBuilder();
		try {
			while(( read = in.read(buf)) > 0) {
				fileSize += read;
				for(int j = 0; j < read; j++) {
					char c = (char) buf[j];
					//System.out.println("character " + c + " at file position " + raf.getFilePointer());
					offset++;
					if(c == '>') {
						readingSeqId = true;
					} else if (readingSeqId) {
						if(c != '\n' & c != '\r') {
							seqIdBuilder.append(c);
						} else {
							setCurrentSequenceId(seqIdBuilder.toString());
							System.out.println("new sequence " + getCurrentSequenceId());
							handler.newSequence(this);
							readingSeqId = false;
							seqIdBuilder = new StringBuilder(); //It used to be possible to reset a StringBuffer
						}
					} else if(c == ' ' || c == '\n' || c == '\r'){
						handler.newQuality(this);
						currentQual = 0;
						i = 0;
					} else {
						currentQual =  Character.getNumericValue(c)  + currentQual * (int)Math.pow(10, i);
						i++;
					}
				}

			}
		} catch( EOFException eof) {
			handler.eof(this);
		} 
		handler.eof(this);
		
	}

	/*
	public static void main(String [] args) throws Exception {
		String file = args[0];
		
		FastaQualityParser fqp = new FastaQualityParser();
		
		fqp.parse(new File(file), new FastaQualityHandler());
	}
	
	
			RandomAccessFile raf = new RandomAccessFile(inputFastaFile, "r");
		setRandomAccessFile(raf);
		System.out.println("File size " + inputFastaFile.length());
		int b = 0;
		int i = 0;
		try {
			while(b >= 0) {

				b  =  raf.read();
				char c = (char) b;
				//System.out.println("character " + c + " at file position " + raf.getFilePointer());

				if(c == '>') {
					setCurrentSequenceId(raf.readLine());
					handler.newSequence(this);
					System.out.println("new sequence " + getCurrentSequenceId());
				}else if(c == ' '){
					handler.newQuality(this);
					currentQual = 0;
					i = 0;
				} else if(c == '\n' || c == '\r'){
					continue;
				} else {
					currentQual =  Character.getNumericValue(c)  + currentQual * (int)Math.pow(10, i);
					i++;
				}


			}
		} catch( EOFException eof) {
			handler.eof(this);
		} finally {		
			raf.close();
		}
		handler.eof(this);
	*/
	
	public int getCurrentQuality() {
		return currentQual;
	}
	
}
