package broad.core.annotation;

public class BEDGraph extends ShortBED {
	
	public BEDGraph(String name) {
		super(name);
	}
	
	public BEDGraph(String chr, int start, int end) {
		super("", chr, start, end);
	}
	
	public BEDGraph(String chr, int start, int end, double score) {
		this(chr, start, end);
		setScore(score);
	}

	public BEDGraph(GenomicAnnotation annot) {
		super(annot);
	}
	
	public BEDGraph(String[] info) {
		this(info[0], Integer.parseInt(info[1]), Integer.parseInt(info[2]), Double.parseDouble(info[3]));
	}
}
