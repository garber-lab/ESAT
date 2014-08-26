package umms.core.feature;

import java.util.Collection;
import java.util.TreeSet;

import umms.core.annotation.Annotation;
import umms.core.annotation.Gene;


public class GeneWindow extends Gene implements Window{

	private Collection<Gene> containedGenes;
	
	public GeneWindow(Gene gene) {
		super(gene);
		this.containedGenes = new TreeSet<Gene>();
	}
	
	
	public GeneWindow(GeneWindow window) {
		super(window);
		this.containedGenes=window.containedGenes;
	}
	
	public GeneWindow(Annotation window) {
		super(window);
		this.containedGenes=new TreeSet<Gene>();
	}

	public GeneWindow(Collection<? extends Annotation> exons) {
		super(exons);
		this.containedGenes=new TreeSet<Gene>();
	}


	/**
	 * Returns the size of the window
	 */
	public int getSize(){
		return super.getTranscriptLength();
	}

	@Override
	public void addSourceAnnotation(Annotation annotation) {
		Gene gene=new Gene(annotation);
		this.containedGenes.add(gene);
	}

	@Override
	public Collection<? extends Annotation> getSourceAnnotations() {
		return this.containedGenes; //This 
	}
	
	public boolean equals(Object o){
		if(! (o instanceof GeneWindow)) { return false;}
		if(super.equals(o)){
			GeneWindow other=(GeneWindow) o;
			//check the contained genes
			if(this.containedGenes==null || other.containedGenes==null){return false;}
			if(this.containedGenes.size()!= other.containedGenes.size()){return false;}
			//check each element against the other set
			for(Gene g: this.containedGenes){
				if(!other.containedGenes.contains(g)){return false;}
			}
			return true;
		}
		return false;
	}
	
	@Override
	public int compareTo(Annotation g){
		boolean done=false;
		if(this.equals(g)){return 0;}
		
		else{
			int result=super.compareTo(g);
			
			//if result is 0 then check the contained genes
			//MG: why does compare care about the underlying annotations?
			if(result==0  && (g instanceof GeneWindow) ){
				GeneWindow other = (GeneWindow) g;
				int size1=-1;
				int size2=-1;
				if(this.containedGenes!=null){size1=this.containedGenes.size();}
				if(other.containedGenes!=null){size2=other.containedGenes.size();}
				result = size2-size1;
				while(result==0 && !done){ //TODO: Is this the intended comparison? What if the set of contained ge
					for(Gene g1: this.containedGenes){
						for(Gene g2: other.containedGenes){
							result=g1.compareTo(g2);
						}
					}
					done=true;
				}
			}
			return result;
		}
		
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
			GeneWindow trimmed=this.trimGene(i, i+windowSize);
			subGenes.add(trimmed);
		}
		return subGenes;
	}

	
}
