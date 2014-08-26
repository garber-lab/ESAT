package broad.core.annotation;

public interface TwoSubjectAnnotation {
	int getAStart();
	int getBStart();
	int getAEnd();
	int getBEnd();
	String getAName();
	String getBName();
	
	LightweightGenomicAnnotation getA();
	LightweightGenomicAnnotation getB();
	
	boolean isDirect();

}
