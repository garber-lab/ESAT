package broad.pda.datastructures;


import java.util.Collection;
import java.util.Iterator;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;


import broad.core.annotation.*;
import broad.core.datastructures.IntervalTree;
import broad.core.error.ParseException;

/**
 * @author engreitz
 * Wrapper class to take advantage of all the code in AnnotationReader for Alignments
 */
@Deprecated

public class AlignmentsReader extends AnnotationReader<Alignments>{
	public AlignmentsReader() {
		super();
	}
	
	public AlignmentsReader(Collection<Alignments> annot) {
		super(annot);
	}
	
	public AlignmentsReader(BufferedReader br) throws ParseException, IOException {
		this();
		throw new UnsupportedOperationException();
	}
	
	public AlignmentsReader(String input) throws ParseException, IOException {
		this();
		throw new UnsupportedOperationException();
	}
	
	public void load(BufferedReader br) throws IOException, ParseException {
		throw new UnsupportedOperationException();
	}
	
	public int parse(BufferedReader br, GenomicAnnotationFilter<Alignments> filter, AnnotationHandler handler) throws ParseException {
		throw new UnsupportedOperationException();
	}

	public Alignments createAnnotation(GenomicAnnotation a) {
		return new Alignments(a);
	}

	public int parse(String file, GenomicAnnotationFilter<Alignments> filter,
			AnnotationHandler handler) throws ParseException, IOException {
		throw new UnsupportedOperationException();
	}
	
	/* (non-Javadoc)
	 * @see broad.core.annotation.AnnotationReader#write(java.io.BufferedWriter)
	 * Going to all this trouble to rewrite because I don't want to replace the toString() function for Alignments just yet.
	 * Otherwise could just call super.write(bw);  TODO
	 */
	@Override
	public void write(BufferedWriter bw) throws IOException {
		Iterator<String> chrIt = getChromosomeAnnotationMap().keySet().iterator();
		while (chrIt.hasNext()) {
			String chr = chrIt.next();
			IntervalTree<Alignments> tree = getChromosomeTree(chr);
			Iterator<Alignments> annotIt = tree.valueIterator();
			while (annotIt.hasNext()) {
				Alignments annot = annotIt.next();
				bw.write(annot.toStringWithScores());
				bw.newLine();
			}
		}
	}
}
