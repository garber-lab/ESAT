package broad.pda.seq.graph;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import nextgen.core.alignment.Alignment;
import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Gene;

import org.apache.commons.collections15.Transformer;
import org.jgrapht.Graph;
import org.jgrapht.GraphPath;
import org.jgrapht.Graphs;

import broad.core.annotation.BasicLightweightAnnotation;
import broad.core.math.EmpiricalDistribution;
import broad.pda.datastructures.Alignments;
import broad.pda.seq.graph.ChromosomeWithBubblesJGraphT.BubbleEdge;
import broad.pda.seq.graph.ChromosomeWithBubblesJGraphT.EdgeSourceType;

public class Path implements Comparable<Path>,GraphPath<Annotation,BubbleEdge>{

	Collection<BubbleEdge> edges;
	Collection<Alignment> pairedEndEdges;
	Transformer<BubbleEdge, Number> transformer;
	String orientation;
	static ChromosomeWithBubblesJGraphT G;
	double spliceWeight=1;
	private double localLambda=0;
	//TODO: Lets store the Nodes and Edges
	int fillInDistance=20;
	private boolean isToComplex;	
	
	public Path(Collection<BubbleEdge> edges, double spliceWeight){
		this(edges, new StandardTransformer(), spliceWeight);
	}
	
	public Path(Collection<BubbleEdge> edges, Transformer<BubbleEdge, Number> transformer, double spliceWeight){
		this.spliceWeight=spliceWeight;
		this.transformer=transformer;
		this.edges=edges;
		this.pairedEndEdges=new ArrayList<Alignment>(); //TODO Consider making this a Set so that reads are unique
		setOrientation();
	}
	
	public Path(double spliceWeight){
		this(new TreeSet<BubbleEdge>(),spliceWeight);
	}
	public Path(){
		this(new TreeSet<BubbleEdge>(),1);
	}
	
	public Path(Gene gene, ChromosomeWithBubblesJGraphT graph){
		this(initialize(gene, graph),1);
	}
	
	private static Collection<BubbleEdge> initialize(Gene gene, ChromosomeWithBubblesJGraphT graph) {
		
		G = graph;
		//assume all are splice edges with count of 1
		Collection<BubbleEdge> edges=new TreeSet<BubbleEdge>();
		Collection<? extends Annotation> exonSet=gene.getSortedAndUniqueExons();
		Object[] exons=exonSet.toArray();
		for(int i=0; i<exons.length-1; i++){
			Alignments left=(Alignments)exons[i];
			Alignments right=(Alignments)exons[i+1];
			Alignments connection=new Alignments(left.getChr(), left.getEnd(), right.getStart());
			BubbleEdge be=new BubbleEdge(connection, left, right, 1, EdgeSourceType.SPLICED);
			be.setParent(graph);
			edges.add(be);
		}
		
		if(exons.length==1){
			Alignments exon=(Alignments)exons[0];
			BubbleEdge be=new BubbleEdge(exon, exon, exon, 1.0, EdgeSourceType.SELF);
			be.setParent(graph);
			edges.add(be);
		}
		return edges;
	}

	//??? DOESNT MAKE SENSE
	public boolean addPairedEnd(Alignment align){
		//check if the alignment overlaps one of the nodes
		for(BubbleEdge edge: this.edges){
			Annotation left=edge.getLeftNode();//??? DO WE NEED THIS?
			Annotation right=edge.getRightNode();//??? DO WE NEED THIS?
			//if so add and return true
			if(this.overlapsLeft(align) && this.overlapsRight(align)){
				return pairedEndEdges.add(align); 
			}
		}
		//else return false
		return false;
	}
	
	//REDUNDANT?????
	public double getScore(){
		double score=0;
		//go through bubble edge and score
		for(BubbleEdge edge: edges){
			score+=transformer.transform(edge).doubleValue();
		}
		return score;
	}
	
    /**
     * Returns the weight assigned to the path. Typically, this will be the sum
     * of the weights of the edge list entries (as defined by the containing
     * graph), but some path implementations may use other definitions.
     *
     * @return the weight of the path
     */
    public double getWeight(){
    	double score=0;
		//go through bubble edge and score
		for(BubbleEdge edge: edges){
			score+=transformer.transform(edge).doubleValue();
		}
		return score;
    }
	/**
	 *  TODO: This may be to slow, maybe compute start on the fly as Path is being constructed
	 * @return genomic start of the path.
	 */
	public int getStart() {
		if(edges.isEmpty()){return -1;}
		int start = Integer.MAX_VALUE;
		for(BubbleEdge edge: edges){
			start = Math.min(start, edge.getLeftNode().getStart());
			start = Math.min(start, edge.getRightNode().getStart());
		}
		return start;
	}
	
	/**
	 *  TODO: This may be to slow, maybe compute end on the fly as Path is being constructed
	 * @return genomic end of the path.
	 */
	public int getEnd() {
		if(edges.isEmpty()){return -1;}
		int end = -Integer.MAX_VALUE;
		for(BubbleEdge edge: edges){
			end = Math.max(end, edge.getLeftNode().getEnd());
			end = Math.max(end, edge.getRightNode().getEnd());
		}
		return end;
	}
	
	public String getChromosome() {
		return edges != null && edges.size() > 0 ? edges.iterator().next().getConnection().getChr() : null;
	}
	
	public boolean addEdge(BubbleEdge edge){
		boolean b= edges.add(edge);
		if(b){setOrientation();}
		return b;
	}
	
	public boolean addEdges(Collection<BubbleEdge> edges) {
		boolean ok = true;
		Iterator<BubbleEdge> edgeIt = edges.iterator();
		while(ok && edgeIt.hasNext()) {
			BubbleEdge edge = edgeIt.next();
			ok = addEdge(edge);
		}
		return ok;
	}
	
	//REDUNDANT???
	public Collection<BubbleEdge> getEdges() {
		return edges;
	}
	
	/**
     * Returns the edges making up the path. The first edge in this path is
     * incident to the start vertex. The last edge is incident to the end
     * vertex. The vertices along the path can be obtained by traversing from
     * the start vertex, finding its opposite across the first edge, and then
     * doing the same successively across subsequent edges; {@link
     * Graphs#getPathVertexList} provides a convenience method for this.
     *
     * <p>Whether or not the returned edge list is modifiable depends on the
     * path implementation.
     *
     * @return list of edges traversed by the path
     */
	public List<BubbleEdge> getEdgeList(){
    	return (List<BubbleEdge>) (new ArrayList<BubbleEdge>(edges));
    }

	
	public boolean removeEdge(BubbleEdge edge){
		boolean b= edges.remove(edge);
		if(b){setOrientation();}
		return b;
	}
	
	public int getNumberOfEdges(){return edges.size();}
	
	
	public Gene toGene(){
		Collection<Annotation> exons = getExons();
		if(exons.size()==0){System.err.println("EXON set is EMPTY");}
		Gene gene= new Gene(exons);
		gene.setOrientation(getOrientation());
		gene.setCountScore(getScore());
		return gene;
	}
	
	public Gene toGene(double[] extraFields){
		Collection<Annotation> exons = getExons();
		if(exons.size()==0){System.err.println("EXON set is EMPTY");}
		Gene gene= new Gene(exons);
		gene.setOrientation(getOrientation());
		gene.setCountScore(getScore());
		gene.setExtraFields(extraFields);
		return gene;
	}

	private Collection<Annotation> getExons() {
		Collection<Annotation> exons=new TreeSet<Annotation>();
		for(BubbleEdge edge: edges){
			//if(edge.getType().equals(EdgeSourceType.PAIRED)){System.err.println(edge.getType());}
			exons.add(edge.getLeftNode());
			exons.add(edge.getRightNode());
		}
		return exons;
	}
	
	private void setOrientation(){
		if(edges.isEmpty()){this.orientation="*".intern();}
		else{
			String orientation=edges.iterator().next().getOrientation();
			for(BubbleEdge edge: edges){
				if("*".equalsIgnoreCase(orientation)|| "*".equalsIgnoreCase(edge.getOrientation()) || orientation.equalsIgnoreCase(edge.getOrientation())){
					orientation=setOrientation(orientation, edge.getOrientation());
				}
				else{System.err.println("Throws "+orientation+" "+edge.getOrientation());throw new IllegalArgumentException("Path can not have conflicting orientations");}
			}
			this.orientation=orientation;
		}
	}
	
	private String setOrientation(String orientation1, String orientation2) {
		if(orientation1.equalsIgnoreCase("*")){return orientation2;}
		if(orientation2.equalsIgnoreCase("*")){return orientation1;}
		return orientation1;
	}
	
	private void setOrientation(String orientation) {
		this.orientation = orientation.intern();
	}

	public String getOrientation(){
		return this.orientation;
	}
	
	public static class StandardTransformer implements Transformer<BubbleEdge, Number>{
		public Number transform(BubbleEdge edge) {
			if(edge.getType().equals(EdgeSourceType.SPLICED)){return new Double(edge.getAllCounts());}
			else if(edge.getType().equals(EdgeSourceType.SELF)){}
			else if(edge.getType().equals(EdgeSourceType.PAIRED)){}
			else if(edge.getType().equals(EdgeSourceType.PAIRED_SUPPORT)){}
			return 0;
		}
	}

	public int compareTo(Path b) {
		/*int result = getAlignment().compareTo(b.getAlignment());
		if(result == 0) {
			List<LightweightGenomicAnnotation> aExons = new ArrayList<LightweightGenomicAnnotation>(getExons());
			List<LightweightGenomicAnnotation> bExons = new ArrayList<LightweightGenomicAnnotation>(b.getExons());
			
			int minLength = Math.min(aExons.size(), bExons.size());
			int idx = 0;
			while(idx < minLength && result == 0) {
				result = aExons.get(idx).compareTo(bExons.get(idx));
				idx++;
			}
			if(result == 0) {
				result = aExons.size() - bExons.size();
			}
		}

		
		return result;*/
		return toGene().compareTo(b.toGene());
	}
	
	public boolean equals(Object o) {
		Path a = (Path)o;
		if (!a.toGene().equals(this.toGene()))
			return false;
		if (a.spliceWeight != this.spliceWeight)
			return false;
		if (a.localLambda != this.localLambda)
			return false;
		if (a.fillInDistance != this.fillInDistance)
			return false;
		if (a.isToComplex != this.isToComplex)
			return false;
		
		return true;
	}
	
	public int hashCode() {
		return toGene().hashCode();
	}
	
	public Alignments getAlignment(){
		Alignments align=new Alignments(getChromosome(), getStart(), getEnd(), getOrientation());
		return align;
	}

	public Collection<Alignment> getPairedEndEdges() {
		return this.pairedEndEdges;
	}

	public void addPairedEndEdges(Collection<Alignment> pairedEndEdges2) {
		for(Alignment align: pairedEndEdges2){
			this.addPairedEnd(align);
		}
	}
	
	public String toString() {return toGene().toBED();}

	public Collection<Annotation> getSortedNodes() {
		Collection<Annotation> rtrn = getExons();
		
		return rtrn;
	}
	

	 /**
     * Returns the start vertex in the path.
     *
     * @return the start vertex
     */
    public Annotation getStartVertex(){
    	Collection<Annotation> sorted=getSortedNodes();
		return (Annotation) sorted.iterator().next();
    }
    
    /**
     * Returns the end vertex in the path.
     *
     * @return the end vertex
     */
    public Annotation getEndVertex(){
    	Collection<Annotation> sorted=getSortedNodes();
		return (Annotation)sorted.toArray()[sorted.size()-1];
    }
    
	//checks whether the path overlaps the left of the read
	//A node has to overlap the start of the read
	public boolean overlapsLeft(Alignment read) {
		for(BubbleEdge edge: this.edges){
			Annotation left=edge.getLeftNode();
			Annotation right=edge.getRightNode();
			//if so add and return true
			if(((left.getStart()<=read.getAlignmentStart()) && left.getEnd()>=read.getAlignmentStart()) || ((right.getStart()<=read.getAlignmentStart()) && right.getEnd()>=read.getAlignmentStart()) ){
				return true;
			}
		}
		return false;
	}
	
	public Collection<Annotation> getOverlappingNodes(Alignment read) {
		TreeSet<Annotation> result = new TreeSet<Annotation>();
		Annotation readAnnotation = new BasicLightweightAnnotation(read.getReferenceName(), read.getAlignmentStart(), read.getAlignmentEnd());
		for(BubbleEdge edge: this.edges){
			Annotation left=edge.getLeftNode();
			Annotation right=edge.getRightNode();
			if(left.overlaps(readAnnotation)) {
				result.add(left);
			}
			if(right.overlaps(readAnnotation)) {
				result.add(right);
			}
		}
		return result;
	}
	
	//checks whether the path overlaps the right of the read
	//A node has to overlap the end of the read
	public boolean overlapsRight(Alignment read) {
		for(BubbleEdge edge: this.edges){
			Annotation left=edge.getLeftNode();
			Annotation right=edge.getRightNode();
			//if so add and return true
			if(((left.getStart()<=read.getAlignmentEnd()) && left.getEnd()>=read.getAlignmentEnd()) || ((right.getStart()<=read.getAlignmentEnd()) && right.getEnd()>=read.getAlignmentEnd()) ){
				return true;
			}
		}
		return false;
	}

	public boolean overlaps(Path leftPath) {
		Alignments align1=new Alignments(this.getChromosome(), this.getStart(), this.getEnd());
		Alignments align2=new Alignments(leftPath.getChromosome(), leftPath.getStart(), leftPath.getEnd());
		return align1.overlaps(align2);
	}

	public double getPercentConsistent(EmpiricalDistribution pairedInsertDistribution) {
		Collection<Alignment> pairedEnds=this.getPairedEndEdges();
		double count=0;
		double total=0;
		
		for(Alignment pair: pairedEnds){
			int distance=0;
			Gene insert=this.toGene().trimAbsolute(pair.getAlignmentStart(), pair.getAlignmentEnd());
			if(insert==null){distance=-1;}
			else{
				distance=insert.getSize();
				double likelihood=pairedInsertDistribution.getCumulativeProbability(distance);
				System.err.println(distance+" "+likelihood);
				if(likelihood>.05 && likelihood<.95){count++;}
				total++;
			}
		}
		return count/total;
	}

	public Path copy() {
		Path copy = new Path(new ArrayList<BubbleEdge>(), transformer, this.spliceWeight);
		copy.orientation = orientation;
		for(BubbleEdge e : edges) {
			copy.addEdge(e);
		}
		
		for(Alignment a : pairedEndEdges) {
			copy.addPairedEnd(a);
		}
		copy.setOrientation(orientation);
		copy.setIsToComplex(isToComplex);
		return copy;
	}

	public double getCoverage() {
		//System.err.println("WEIGHT "+spliceWeight);
		double sum=0;
		Map<Annotation, Double> exons=new TreeMap<Annotation, Double>();
		for(BubbleEdge edge: edges){
			exons.put(edge.getLeftNode(), edge.getVertexCount(edge.getLeftNode()));
			exons.put(edge.getRightNode(), edge.getVertexCount(edge.getRightNode()));
			sum+=(edge.getSplicedCounts()*spliceWeight);
		}
		
		for(Annotation exon: exons.keySet()){
			sum+=exons.get(exon);
		}
		return sum;
	}

	public int getSize() {
		return this.toGene().getSize();
	}

	public void setLocalLambda(double localLambda) {
		this.localLambda=localLambda;
	}
	
	public double getLocalLambda(){return this.localLambda;}

	public Gene toFilledInGene() {
		Collection<Annotation> exons=new TreeSet<Annotation>();
		boolean addNodes=true;
		for(BubbleEdge edge: edges){
			if(edge.getType().equals(EdgeSourceType.PAIRED)){
				//System.err.println(edge.getType());
				//TODO Make a decision, is this likely to an exon? 
				int connectionLength=edge.getConnection().length();
				
				//TODO for now I'm going to hard code a distance but this should come from the insert dist
				if(connectionLength<this.fillInDistance){
					//make new alignments from left to right
					Annotation fillIn=mergeAlignments(edge.getLeftNode(), edge.getRightNode());
					exons.add(fillIn);
					addNodes=false;
				}
			}
			if(addNodes){
				exons.add(edge.getLeftNode());
				exons.add(edge.getRightNode());
			}
		}
		//exons=CollapseByIntersection.CollapseByIntersection(exons, false); //TODO get this to work
		Gene gene= new Gene(exons);
		gene.setOrientation(getOrientation());
		gene.setCountScore(getScore());
		return gene;
	}

	private Alignments mergeAlignments(Annotation leftNode, Annotation rightNode) {
		Alignments rtrn=new Alignments(leftNode.getChr(), Math.min(leftNode.getStart(), rightNode.getStart()), Math.max(leftNode.getEnd(), rightNode.getEnd()));
		return rtrn;
	}

	public double getNodeCount(Annotation v1) {
		for(BubbleEdge edge: this.edges){
			Annotation left=edge.getLeftNode();
			Annotation right=edge.getRightNode();
			if(left.equals(v1)){return edge.getVertexCount(left);}
			if(right.equals(v1)){return edge.getVertexCount(right);}
		}
		return 0.0;
	}

	public void setSpliceWeight(double weight){this.spliceWeight=weight;}

	public double getSpliceWeight() {
		return this.spliceWeight;
	}

	public void setIsToComplex(boolean isToComplex) {this.isToComplex = isToComplex;}
	
	public boolean isToComplex() {return isToComplex;}

	/**
     * Returns the graph over which this path is defined. The path may also be
     * valid with respect to other graphs.
     *
     * @return the containing graph
     */
    public Graph<Annotation, BubbleEdge> getGraph(){
    	return (Graph<Annotation, BubbleEdge>) G;
    }

}
