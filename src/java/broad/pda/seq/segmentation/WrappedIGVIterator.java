package broad.pda.seq.segmentation;

import net.sf.samtools.util.CloseableIterator;
import nextgen.core.alignment.SingleEndAlignment;

public class WrappedIGVIterator implements CloseableIterator<nextgen.core.alignment.Alignment>{

	
		CloseableIterator<org.broad.igv.sam.Alignment> iter;
		
		public WrappedIGVIterator(CloseableIterator<org.broad.igv.sam.Alignment> iter) {
			this.iter = iter;
		}
		
		 public void close() {
		        iter.close();
		    }

		    public boolean hasNext() {
		        return iter.hasNext();
		    }

		    public nextgen.core.alignment.Alignment next() {
		    	return convert(iter.next());
		    }

		    private nextgen.core.alignment.Alignment convert(org.broad.igv.sam.Alignment next) {
		    	//TODO Take the IGV Alignment and make one of our alignments
		    	nextgen.core.alignment.Alignment rtrn=new SingleEndAlignment(next);
		    	return rtrn;
		    }

			public void remove() {
		        iter.remove();
		    }
        
	
	
}
