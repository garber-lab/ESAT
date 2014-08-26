package broad.core.annotation;


public class BasicTwoSubjectAnnotation implements TwoSubjectAnnotation {
	private LightweightGenomicAnnotation a;
	private LightweightGenomicAnnotation b;
	
	public BasicTwoSubjectAnnotation() {
		a = new BasicGenomicAnnotation("A");
		b = new BasicGenomicAnnotation("B");
	}
	public LightweightGenomicAnnotation getA() {
		return a;
	}

	public int getAEnd() {
		return a.getEnd();
	}

	public String getAName() {
		return a.getName();
	}

	public int getAStart() {
		return a.getStart();
	}

	public LightweightGenomicAnnotation getB() {
		return b;
	}

	public int getBEnd() {
		return b.getEnd();
	}

	public String getBName() {
		return b.getName();
	}

	public int getBStart() {
		return b.getStart();
	}
	
	public void setA(LightweightGenomicAnnotation a) { this.a = a;}
	public void setB(LightweightGenomicAnnotation b) { this.b = b;}

	public boolean isDirect() {
		return !a.inReversedOrientation() && ! b.inReversedOrientation();
	}
	
	public String toString() {
		return getA().getChromosome() + ":" + getA().getStart() + "-" + getA().getEnd() + "\t("+(isDirect() ? "+" : "-")+")\t" + getB().getChromosome() + ":" + getB().getStart() + "-" + getB().getEnd();
	}

}
