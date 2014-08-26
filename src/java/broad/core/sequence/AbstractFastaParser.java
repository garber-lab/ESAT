package broad.core.sequence;

import java.io.IOException;

public abstract class AbstractFastaParser {

	public static final int BUFFER_SIZE = 8192;
	private String currentSequenceId;
	//private RandomAccessFile raf;
	protected long fileSize;
	protected long offset;
	protected char currentBase;
	
	public String getCurrentSequenceId() {
		return currentSequenceId;
	}
	
	public long getOffset() throws IOException {
		return offset;
	}

	/*
	protected RandomAccessFile getRandomAccessFile() {
		return raf;
	}

	protected void setRandomAccessFile(RandomAccessFile raf) {
		this.raf = raf;
	}
	*/
	protected void setCurrentSequenceId(String currentSequenceId) {
		this.currentSequenceId = currentSequenceId;
	}
	

	public long getFileSize() {
		return fileSize;
	}
	
	public char getCurrentBase() {
		return currentBase;
	}

}