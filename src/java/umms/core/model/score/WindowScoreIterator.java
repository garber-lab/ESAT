package umms.core.model.score;

import java.util.Iterator;

import org.apache.log4j.Logger;

import umms.core.annotation.Annotation;
import umms.core.annotation.Gene;
import umms.core.feature.Window;
import net.sf.samtools.util.CloseableIterator;

public class WindowScoreIterator<T extends WindowScore> implements CloseableIterator<T> {

	Iterator<? extends Annotation> itr;
	WindowProcessor<T> processor;
	T previousScore = null;
	static Logger logger = Logger.getLogger(Gene.class.getName());

	public WindowScoreIterator(Iterator<? extends Annotation> windowIterator, WindowProcessor<T> processor, Annotation region){
		this.itr = windowIterator;
		this.processor = processor;
		
		// Don't bother initializing the region if the iterator is empty
		try {
			if (itr.hasNext()) processor.initRegion(region);
		} catch(NullPointerException e) {
			logger.info("Failing on " + region.toUCSC());
		}
	}

	@Override
	public boolean hasNext() {
		return this.itr.hasNext();
	}

	@Override
	public T next() {
		Annotation w = this.itr.next();
		T score= processor.processWindow(w, previousScore);
		//T score= processor.processWindow(w);
		previousScore=score;
		return score;
	}

	@Override
	public void remove() {
		throw new UnsupportedOperationException();
	}
	
	@Override
	public void close() {
		processor.finishedRegion();
	}
}
