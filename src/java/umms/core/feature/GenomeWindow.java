package umms.core.feature;

import java.util.Collection;
import java.util.TreeSet;

import umms.core.annotation.Annotation;
import umms.core.annotation.BasicAnnotation;

public class GenomeWindow extends BasicAnnotation implements Window {

	private Collection<Annotation> containedRegions;
	
	public GenomeWindow(String chr, int fragmentStart, int fragmentEnd) {
		super(chr, fragmentStart, fragmentEnd);
		this.containedRegions=new TreeSet<Annotation>();
	}

	public GenomeWindow(BasicAnnotation annotation) {
		super(annotation);
		this.containedRegions=new TreeSet<Annotation>();
	}

	@Override
	public void addSourceAnnotation(Annotation annotation) {
		containedRegions.add(annotation);
	}

	@Override
	public Collection<? extends Annotation> getSourceAnnotations() {
		return this.containedRegions;
	}

	/**
	 * Get collection of windows spanning the annotation
	 * @param windowSize Window size
	 * @param stepSize Step size
	 * @return The collection of windows
	 */
	public Collection<Window> getWindows(int windowSize, int stepSize) {
		Collection<Window> subGenes=new TreeSet<Window>();
		for(int i=0; i<getSize(); i=i+stepSize){
			int start=i+this.getStart();
			int end=start+windowSize;
			Window w=new GenomeWindow(this.getChr(), start, end);
			subGenes.add(w);
		}
		return subGenes;
	}

}
