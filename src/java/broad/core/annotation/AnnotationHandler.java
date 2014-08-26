package broad.core.annotation;

public interface AnnotationHandler {

	void track(String line);

	void browserLine(String line);
	
	void annotation(GenomicAnnotation annotation);
	
	void eof();
	
	void begin();

}
