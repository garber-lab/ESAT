package broad.core.annotation;

import broad.core.error.ParseException;

public class RepeatMaskerAnnotation extends BasicGenomicAnnotation {
	private int bin, mismatches, deleted, inserted, basesLeftAfterMatch, repStart, repEnd, repLeft;
	private String repName, repClass, repFamily;

	public RepeatMaskerAnnotation() {
		super();
	}
	
	public RepeatMaskerAnnotation(String name) {
		super(name);
	}
	
	public RepeatMaskerAnnotation(GenomicAnnotation ga) {
		super(ga.getName(), ga.getChromosome(), ga.getStart(), ga.getEnd());
	}
	
	
	public RepeatMaskerAnnotation(String[] info) throws ParseException {
		this();
		
		this.bin = Integer.parseInt(info[0]);
		setScore(Integer.parseInt(info[1]));
		this.mismatches = Integer.parseInt(info[2]);
		this.deleted = Integer.parseInt(info[3]);
		this.inserted = Integer.parseInt(info[4]);
		setChromosome(info[5]);
		setStart(Integer.parseInt(info[6]));
		setEnd(Integer.parseInt(info[7]));
		this.basesLeftAfterMatch = Integer.parseInt(info[8]);
		setOrientation(info[9]);
		this.repName = info[10];
		this.repClass = info[11];
		this.repFamily = info[12];
		this.repStart = Integer.parseInt(info[13]);
		this.repEnd = Integer.parseInt(info[14]);
		this.repLeft = Integer.parseInt(info[15]);
	}
	
	public RepeatMaskerAnnotation(String chr, int start, int end, String strand, String repName, String repClass, String repFamily) {
		super(repName, chr, start, end);
		setOrientation(strand);
		this.repName = repName;
		this.repClass = repClass;
		this.repFamily = repFamily;
	}

	public RepeatMaskerAnnotation(int bin, int alignmentScore, int mismatches,
			int deleted, int inserted, String chr, int start, int end, int basesLeftAfterMatch, String strand, String repName, 
			String repClass, String repFamily, int repStart, int repEnd, int repLeft) {
		super(repName, chr, start, end);
		setOrientation(strand);
		this.bin = bin;
		setScore(alignmentScore);
		this.mismatches = mismatches;
		this.deleted = deleted;
		this.inserted = inserted;
		this.basesLeftAfterMatch = basesLeftAfterMatch;
		this.repName = repName;
		this.repClass = repClass;
		this.repFamily = repFamily;
		this.repStart = repStart;
		this.repEnd = repEnd;
		this.repLeft = repLeft;
	}
	
	public String getRepeatName() { return repName; }
	public String getRepeatClass() { return repClass; }
	public String getRepeatFamily() { return repFamily; }
	
	public String toString() {
		StringBuffer buf = new StringBuffer(getChromosome());
		buf.append("\t")
			.append(getStart())
			.append("\t")
			.append(getEnd())
			.append("\t")
			.append(getScore())
			.append("\t")
			.append(getRepeatName())
			.append("\t")
			.append(getRepeatClass())
			.append("\t")
			.append(getRepeatFamily());
		return buf.toString();
	}
	
}
