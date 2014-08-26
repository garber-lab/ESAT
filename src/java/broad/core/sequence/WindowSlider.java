package broad.core.sequence;

import java.util.Iterator;

public class WindowSlider implements Iterator<SequenceRegion> {
	
	private int windowSize;
	private int overlap;
	private int atPosition = 0; 
	private int seqStart;
	private Sequence seq;

	public static WindowSlider getSlider(Sequence seq, int windowSize, int overlap) {
		WindowSlider ws = new WindowSlider(windowSize, overlap);
		ws.seq = seq;
		ws.seqStart = 0;
		return ws;
	}
	
	public static WindowSlider getSlider(SequenceRegion seq, int windowSize, int overlap) {
		//System.out.println("Overloaded method used");
		WindowSlider ws = new WindowSlider(windowSize, overlap);
		ws.seq = seq;
		ws.atPosition = seq.getStart();
		ws.seqStart = seq.getStart();
		return ws;
	}
	
	WindowSlider(int windowSize, int overlap) {
		this.windowSize = windowSize;
		this.overlap    = overlap;
	}
	
	public boolean hasNext() {
		//System.out.println("atPosition " + atPosition  + " seq_end " + seq.getEnd());
		return atPosition < seq.getEnd() - windowSize;
	}

	public SequenceRegion next() {
		SequenceRegion window = seq.getRegion(atPosition - seqStart, atPosition - seqStart + windowSize);
		window.setStart(atPosition);
		window.setEnd(atPosition + windowSize);
		atPosition = atPosition + windowSize - overlap;
		return window;
	}

	public void remove() {
		// DO NOTHING

	}

}
