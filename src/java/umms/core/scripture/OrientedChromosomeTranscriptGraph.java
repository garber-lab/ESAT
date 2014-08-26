package umms.core.scripture;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.TreeSet;

import umms.core.annotation.Annotation;
import umms.core.annotation.Annotation.Strand;
import umms.core.annotation.AbstractAnnotation;
import umms.core.annotation.BasicAnnotation;
import umms.core.annotation.Gene;

import org.apache.log4j.Logger;
import org.jgrapht.EdgeFactory;
import org.jgrapht.GraphPath;
import org.jgrapht.Graphs;
import org.jgrapht.alg.ConnectivityInspector;
import org.jgrapht.alg.KShortestPaths;
import org.jgrapht.graph.DefaultDirectedWeightedGraph;
import org.jgrapht.graph.DefaultWeightedEdge;
import org.jgrapht.graph.GraphPathImpl;
import org.jgrapht.util.VertexPair;

import broad.core.annotation.LightweightGenomicAnnotation;
import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.pda.datastructures.Alignments;
import broad.pda.seq.graph.Path;

public class OrientedChromosomeTranscriptGraph extends DefaultDirectedWeightedGraph<Annotation, OrientedChromosomeTranscriptGraph.TranscriptGraphEdge > {
	private static final long serialVersionUID = 1302380140695950943L;
	private static Logger logger = Logger.getLogger(OrientedChromosomeTranscriptGraph.class.getName());
	private static final int MAX_PATHS = 20;

	private static final String quote="\"";
	private static final double scaleFactor = 0.0025;
	
	private String name;
	private Strand orientation;
	private IntervalTree<Annotation> vertices;
	private IntervalTree<OrientedChromosomeTranscriptGraph.TranscriptGraphEdge> edgeTree;


	protected OrientedChromosomeTranscriptGraph(EdgeFactory<Annotation, TranscriptGraphEdge> edgeFactory) {
		super(edgeFactory);
		vertices = new IntervalTree<Annotation>();
		edgeTree = new IntervalTree<OrientedChromosomeTranscriptGraph.TranscriptGraphEdge>();
		this.name = "noname";
	}

	public OrientedChromosomeTranscriptGraph(String name, String orientation) {
		this(new OrientedChromosomeTranscriptGraph.TranscriptGraphEdgeFactory()); 
		this.name = name;
		this.orientation = AbstractAnnotation.getStrand(orientation);
	}

	public TranscriptGraphEdge addEdge(Annotation v1, Annotation v2) {
		return addEdge(v1,v2,false);
	}

	public boolean addAnnotationToGraph(Annotation a) {
		List<? extends Annotation> exonsFromAnnotation = a.getBlocks(true);
		boolean result = a.isUnoriented() || a.getStrand() == this.orientation;
		
		if(!result) { return false;}
		
		for(int i = 0; i < exonsFromAnnotation.size() - 1; i++) {
			addEdge(exonsFromAnnotation.get(i), exonsFromAnnotation.get(i+1)); //TODO: update result
			
		}
		if(exonsFromAnnotation.size() == 1) {
			result = addVertex(a);
		}
		
		return result;
	}

	public boolean addVertex(Annotation v) { 
		//TODO: Why checking again?
		boolean result = canAdd(v);
		if(!result) { return false;}

		
		result = super.addVertex(v);

		if(result) {
			//logger.info(" Vertex is added.");
			vertices.put(v.getStart(), v.getEnd(), v);
		}


		return result;
	}

	/**
	 * Adds an exon to a graph
	 * @param v
	 * @return
	 */
	public int connectVertexToGraph(Annotation v) {
		int result = -1;
		boolean canAdd = canAdd(v);
		if(!canAdd) { return -1;}

		Collection<TranscriptGraphEdge> abuttingEdges = getAbuttingEdges(v);
		//logger.error("Size of the abutting edges = "+abuttingEdges.size());
		if(abuttingEdges.size() == 0) {
			boolean success = addVertex(v);
			if(success) { result = 0;}
		} else {
			//This is an ugly handling of unnoriented transcripts, but I could come with nothing better
			//If the annotation is oriented keep it as is
			Annotation vToAdd = v;
			if(v.isUnoriented()) {
				vToAdd = new Gene(v);
				vToAdd.setOrientation(this.orientation);
			} 
			result = 0;
			for (TranscriptGraphEdge e : abuttingEdges) {
				TranscriptGraphEdge ne =  addNewEdgeFromExistingByKeepingOneVertex(vToAdd, e);
				if(ne == null) {
					result = -1;
					break;
				}
				result ++;
			}
		}
		return result;

	}

	/**
	 * Returns true if the specified annotation is unoriented or has the same orientation as this graph
	 * @param v
	 * @return
	 */
	private boolean canAdd(Annotation v) {
		return v.isUnoriented() || v.getOrientation() == this.orientation;
	}

	/**
	 * ?????????????????????????
	 * @param v
	 * @return
	 */
	private Collection<TranscriptGraphEdge> getAbuttingEdges(Annotation v) {		
		Iterator<Node<TranscriptGraphEdge>> startOverlapIt = edgeTree.overlappers(v.getStart() - 1, v.getStart());
		Iterator<Node<TranscriptGraphEdge>> endOverlapIt = edgeTree.overlappers(v.getEnd(), v.getEnd()+1);
		
		Collection<TranscriptGraphEdge> abuttingEdges = new TreeSet<TranscriptGraphEdge> ();
		
		while (startOverlapIt.hasNext()) {
			Node<TranscriptGraphEdge> node = startOverlapIt.next();
			TranscriptGraphEdge overlappingEdge = node.getValue();
			if(overlappingEdge.getEnd() + 1 == v.getStart()) {
				abuttingEdges.addAll(node.getContainedValues());
			}
		}

		while (endOverlapIt.hasNext()) {
			Node<TranscriptGraphEdge> node = endOverlapIt.next();
			TranscriptGraphEdge overlappingEdge = node.getValue();
			if(overlappingEdge.getStart()  == v.getEnd()) {
				abuttingEdges.addAll(node.getContainedValues());
			}
		}
		return abuttingEdges;
	}

	public TranscriptGraphEdge addEdge(Annotation v1, Annotation v2, boolean augmentCount) {
		//The complexity is because if a vertex exists the same instance that is already in the graph must be used.
		
		Node<Annotation> v1InGraphNode = vertices.find(v1.getStart(), v1.getEnd());
		Node<Annotation> v2InGraphNode = vertices.find(v2.getStart(), v2.getEnd());
		Annotation v1InGraph = null;
		Annotation v2InGraph = null;

		if(v1InGraphNode == null) {
			addVertex(v1);
			v1InGraph = v1;
		} else {
			v1InGraph = v1InGraphNode.getValue();
		}
		if(v2InGraphNode == null) {
			addVertex(v2);
			v2InGraph = v2;
		} else {
			v2InGraph = v2InGraphNode.getValue();
		}

		TranscriptGraphEdge e = null;
		if(containsEdge(v1InGraph, v2InGraph)) {
			e = getEdge(v1InGraph, v2InGraph);
			if( augmentCount) {
				setEdgeWeight(e, getEdgeWeight(e)+1);
			}
			//logger.info("Adding new edge " + e.toUCSC() + " to graph between " + v1.toUCSC() + " and " + v2.toUCSC());
		}else {
			e = getEdgeFactory().createEdge(v1InGraph, v2InGraph);
			boolean success = super.addEdge(v1InGraph, v2InGraph, e);
			edgeTree.put(e.getStart(), e.getEnd(), e);
			//logger.info("Adding new edge " + e.toUCSC() + " to graph between " + v1.toUCSC() + " and " + v2.toUCSC());
		}	
		
		return e;
	}

	public List<GraphPath<Annotation, TranscriptGraphEdge>> getPaths() {
		List<GraphPath<Annotation, TranscriptGraphEdge>> paths = new ArrayList<GraphPath<Annotation, TranscriptGraphEdge>>();

		//ORPHAN PATHS
		//paths.addAll(getOrphanPaths());
		
		//NON-ORPHAN PATHS
		Collection<Annotation> sources = getSourceVertices();
		ConnectivityInspector<Annotation, TranscriptGraphEdge> ci = new ConnectivityInspector<Annotation, OrientedChromosomeTranscriptGraph.TranscriptGraphEdge>(this);
		
		//For each source
		for(Annotation s : sources) {
			KShortestPaths<Annotation, TranscriptGraphEdge> shortesPathAlg = new KShortestPaths<Annotation, OrientedChromosomeTranscriptGraph.TranscriptGraphEdge>(this, s, MAX_PATHS);
			java.util.Set<Annotation> connectedSetOfS =  ci.connectedSetOf(s);
			logger.debug("Connected set of " + s.toUCSC() + " is " + connectedSetOfS);
			for(Annotation t : connectedSetOfS) {
				if(isSink(t)) {
					List<GraphPath<Annotation, TranscriptGraphEdge>> sToTPaths = shortesPathAlg.getPaths(t);
					if(sToTPaths == null) {
						//logger.debug("Paths from " + s + " to " + t + " where null ");
					} else {
						//logger.debug("Paths from " + s + " to " + t + ": " + sToTPaths);
						paths.addAll(sToTPaths);
					}
				}
			}
		}

		return paths;
	}

	public static Gene pathToGene(GraphPath<Annotation, TranscriptGraphEdge> gp) {
		List<Annotation> pathVertices =  Graphs.getPathVertexList(gp);
		return new Gene(pathVertices);
	}

	/**
	 * Returns a collection of all the possible source/start vertices for the specified graph.
	 * 
	 * @return All vertices having allowed edge types.
	 */
	public Collection<Annotation> getSourceVertices() {
		//GET ALL VERTICES
		Collection<Annotation> vertices = vertexSet();
		Collection<Annotation> sources = new TreeSet<Annotation>();
		//FOR EACH VERTEX
		for(Annotation v: vertices) {
			/*
			 * IF THERE ARE NO INCOMING EDGES TO THIS VERTEX 
			 * AND
			 * THERE IS AT LEAST 1 OUTGOING EDGE FROM THIS VERTEX
			 */
			if(isSource(v) ) {
				sources.add(v);
			}
		}

		return sources;
	}


	/**
	 * Returns a collection of all the possible target/end vertices for this graph.
	 * 
	 * @return All vertices having allowed edge types.
	 */
	public Collection<Annotation> getSinkVertices() {

		//GET ALL VERTICES
		Collection<Annotation> vertices = vertexSet();
		Collection<Annotation> targets = new TreeSet<Annotation>();
		//FOR EACH VERTEX
		for(Annotation v: vertices) {
			//System.err.println(v.toUCSC()+" "+inOrientedDegree(v, validTypes)+" "+outOrientedDegree(v, validTypes) + " " + isSpanningGap(v.getStart(), v.getEnd()));
			/*
			 * IF THERE ARE NO OUTGOING EDGES FROM THIS VERTEX 
			 * AND
			 * THERE IS AT LEAST 1 INCOMING EDGE TO THIS VERTEX
			 */
			if(isSink(v) ) {
				//System.out.println("adding source");
				targets.add(v);
			}
		}

		return targets;
	}

	/**
	 * Returns a collection of all orphan vertices.
	 * 
	 * @return All vertices having allowed edge types.
	 */
	public Collection<Annotation> getOrphanVertices() {

		//GET ALL VERTICES
		Collection<Annotation> vertices = vertexSet();
		Collection<Annotation> targets = new TreeSet<Annotation>();
		//FOR EACH VERTEX
		for(Annotation v: vertices) {
			if(isOrphan(v) ) {
				//System.out.println("adding source");
				targets.add(v);
			}
		}

		return targets;
	}

	private boolean isSink(Annotation v) {
		return outgoingEdgesOf(v).size() == 0 && incomingEdgesOf(v).size()  >= 1;
	}

	private boolean isSource(Annotation v) {
		return incomingEdgesOf(v).size() == 0 && outgoingEdgesOf(v).size() >= 1;
	}

	private boolean isOrphan(Annotation v) {
		return incomingEdgesOf(v).size() == 0 && outgoingEdgesOf(v).size() == 0;
	}

	private TranscriptGraphEdge addNewEdgeFromExistingByKeepingOneVertex(Annotation newV, TranscriptGraphEdge e) {
		TranscriptGraphEdge addedEdge = null;
		Annotation target = getEdgeTarget(e);
		Annotation source = getEdgeSource(e);
		if(target.overlaps(newV)) {
			addedEdge = addEdge(source, newV);
		} else {
			addedEdge = addEdge(newV, target);
		}
		
		return addedEdge;
	}
	
	/**
	 * Writes the graph to a file
	 * @param save
	 * @param region
	 * @param setUpShortExonLabels
	 * @param graphName
	 * @throws IOException
	 */
	public void writeGraph(String save, Alignments region, boolean setUpShortExonLabels, String graphName) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		writer.write("digraph "+ graphName+" {\n rankdir=LR; \n node [shape = rectangle,color=red];\n");
		//go through each and node and write it and its edges
		//get all the graphs edges
		//start by writing an invisible connection between every node and every node in sorted order
		Iterator<Node<Annotation>> iter=vertices.iterator();
		int counter=0;
		while(iter.hasNext()){ //TODO: 1) Consider moving to the bottom 2) Consider not adding invisible edge to edges that already exist
			Annotation node=iter.next().getValue();
			if(node.overlaps(region)){
				if(counter>0){writer.write("->");}
				writer.write(quote+node.toUCSC()+quote);
				counter++;
			}
		}
		writer.write("[style=invis];\n"); //TODO: Scale the invisible link by the genomic size
		iter=vertices.iterator();
		//Then define the sizes of each node based on genomic sizes
		while(iter.hasNext()){
			Annotation exon=iter.next().getValue();
			//double exonCount=this.getCount(exon);
			//double localRate=this.getLocalRate(exon);
			if(exon.overlaps(region)){
				double normWidth=exon.length()*scaleFactor;
				writer.write(quote+exon.toUCSC()+quote+" [width="+normWidth+", fixedsize=true");//, label="+exonCount+", comment="+localRate); 
				writer.write("];\n");
			}
		}
		
		//This defines actual edges
		Collection<TranscriptGraphEdge> edges = edgeSet();
		for(TranscriptGraphEdge edge: edges){
			VertexPair<Annotation> nodes=getNodePair(edge);
			//System.err.println("WRITEGRAPH - edge" + edge.getConnection().toUCSC());
			Annotation first=nodes.getFirst();
			Annotation second=nodes.getSecond();
			if(first.overlaps(region) || second.overlaps(region)){
				writer.write(quote+first.toUCSC()+quote+"->"+quote+second.toUCSC()+quote);
				//writer.write(", minlen="+(edge.getConnection().length()*(scaleFactor/3))); //TODO: Add scaling for intron sizes (default no scaling)
				writer.write("];\n");
			}
		}
					
		writer.write("}");
		
		writer.close();
	}


	/**
	 * Returns the pair of vertices for a specified edge
	 * @param edge
	 * @return
	 */
	public VertexPair<Annotation> getNodePair(TranscriptGraphEdge edge) {
		return (new VertexPair<Annotation>(edge.getLeftNode(),edge.getRightNode()));
	}
	
	
	/**
	 * Returns all orphan paths in graph (self edges for orphan vertices)
	 * @return
	 */
	public Collection<GraphPath<Annotation, TranscriptGraphEdge>> getOrphanPaths(){
		List<GraphPath<Annotation, TranscriptGraphEdge>> paths = new ArrayList<GraphPath<Annotation, TranscriptGraphEdge>>();
		
		Iterator<Annotation> iter=getOrphanVertices().iterator();
		//Iterate through all vertices
		while(iter.hasNext()){
			Annotation align=iter.next();
			
			//Form new edge
			List<TranscriptGraphEdge> edge=new ArrayList<TranscriptGraphEdge>();
			edge.add(new TranscriptGraphEdge(new BasicAnnotation(align.getChr(), align.getEnd(), align.getEnd())));
			//edge.setParent(this);
			GraphPath<Annotation, TranscriptGraphEdge> path=new GraphPathImpl<Annotation, TranscriptGraphEdge>(this, align, align,edge, 0);
			//edge.setNodeCount(align, getCount(align));
			paths.add(path);	
		}
		
		return paths;
	}

	
	public static class TranscriptGraphEdgeFactory implements EdgeFactory<Annotation,  OrientedChromosomeTranscriptGraph.TranscriptGraphEdge> {

		public TranscriptGraphEdge createEdge(Annotation arg0, Annotation arg1) {
			TranscriptGraphEdge edge;
			if(!arg0.getChr().equals(arg1.getChr()) ) {
				throw new IllegalArgumentException("Cannot create an edge between annotations in different chromosomes: a1 " + arg0.toUCSC() + " a2 " + arg1.toUCSC());
			} 

			if(arg0.overlaps(arg1) ) {
				throw new IllegalArgumentException("Cannot create an edge between overlapping annotations " + arg0.toUCSC() + " a2 " + arg1.toUCSC());
			} 


			if (arg0.compareTo(arg1) <= 0) {
				edge = new TranscriptGraphEdge(new BasicAnnotation(arg0.getChr(), arg0.getEnd(), arg1.getStart()));
			} else {
				edge = new TranscriptGraphEdge(new BasicAnnotation(arg0.getChr(), arg1.getEnd(), arg0.getStart()));
			}
			return  edge;
		}



	}

	public static class TranscriptGraphEdge extends DefaultWeightedEdge implements Annotation{


		private static final long serialVersionUID = -2038056582331579372L;
		private Annotation annotation;
		private Annotation v1;
		private Annotation v2;

		public TranscriptGraphEdge (Annotation a) {
			this.annotation = a;
		}

		public TranscriptGraphEdge (Annotation a,Annotation v1,Annotation v2) {
			this.annotation = a;
			this.v1 = v1;
			this.v2 = v2;
		}
		
		public Annotation getLeftNode(){
			return v1;
		}
		
		public Annotation getRightNode(){
			return v2;
		}
		
		public int compareTo(Annotation o) {
			return annotation.compareTo(o);
		}

		public int getStart() {
			return annotation.getStart();
		}

		public int getSAMStart() {
			return annotation.getSAMStart();
		}

		public int getEnd() {
			return annotation.getEnd();
		}

		public int getSAMEnd() {
			return annotation.getSAMEnd();
		}

		public String getReferenceName() {
			return annotation.getReferenceName();
		}

		public String getChr() {
			return annotation.getChr();
		}

		public String getName() {
			return annotation.getName();
		}

		public Strand getOrientation() {
			return annotation.getOrientation();
		}

		public Strand getStrand() {
			return annotation.getStrand();
		}

		public boolean hasOrientation() {
			return annotation.hasOrientation();
		}

		public boolean isNegativeStrand() {
			return annotation.isNegativeStrand();
		}

		public boolean isUnoriented() {
			return annotation.isUnoriented();
		}

		public int numBlocks() {
			return annotation.numBlocks();
		}

		public List<? extends Annotation> getBlocks() {
			return annotation.getBlocks();
		}

		public List<? extends Annotation> getBlocks(boolean oriented) {
			return annotation.getBlocks(oriented);
		}

		public int length() {
			return annotation.length();
		}

		public int getSize() {
			return annotation.getSize();
		}

		public int size() {
			return annotation.size();
		}

		public int getLengthOnReference() {
			return annotation.getLengthOnReference();
		}

		public int getOrientedStart() {
			return annotation.getOrientedStart();
		}

		public int getOrientedEnd() {
			return annotation.getOrientedEnd();
		}


		public double getScore() {
			return annotation.getScore();
		}

		public int getReferenceCoordinateAtPosition(int positionInAnnotation) {
			return annotation.getReferenceCoordinateAtPosition( positionInAnnotation);
		}

		public int getPositionAtReferenceCoordinate(int referenceCoordinate) {
			return annotation.getPositionAtReferenceCoordinate(referenceCoordinate);
		}

		public int getReferenceCoordinateAtPosition(int positionInAnnotation,
				boolean ignoreOrientation) {
			return annotation.getReferenceCoordinateAtPosition(positionInAnnotation, ignoreOrientation);
		}

		public int getPositionAtReferenceCoordinate(int referenceCoordinate,boolean ignoreOrientation) {
			return annotation.getPositionAtReferenceCoordinate(referenceCoordinate, ignoreOrientation);
		}

		public void setStart(int start) {
			annotation.setStart(start);
		}

		public void setEnd(int end) {
			annotation.setEnd(end);
		}

		public void setOrientation(char orientation) {
			annotation.setOrientation(orientation);
		}

		public void setOrientation(Strand orientation) {
			annotation.setOrientation(orientation);
		}

		public void setOrientedStart(int orientedStart) {
			annotation.setOrientedStart(orientedStart);
		}

		public void setOrientedEnd(int orientedEnd) {
			annotation.setOrientedEnd(orientedEnd);
		}

		public void setReferenceName(String refName) {
			annotation.setReferenceName(refName);
		}

		public void setName(String name) {
			annotation.setName(name);
		}

		public void setScore(double score) {
			annotation.setScore(score);
		}

		public boolean equals(Annotation other) {
			return annotation.equals(other);
		}

		public void expand(int deltaStart, int deltaEnd) {
			annotation.expand(deltaStart, deltaEnd);
		}

		public Annotation trim(int deltaStart, int deltaEnd) {
			return annotation.trim(deltaStart, deltaEnd);
		}

		public void shift(int delta) {
			annotation.shift(delta);
		}

		public void moveToCoordinate(int coordinateInReference) {
			annotation.moveToCoordinate(coordinateInReference);
		}

		public Annotation copy() {
			return annotation.copy();
		}

		public List<Annotation> disect(Annotation a) {
			return annotation.disect(a);
		}

		public List<Annotation> disect(List<? extends Annotation> disectors) {
			return annotation.disect(disectors);
		}

		public Annotation minus(Annotation other) {
			return annotation.minus(other);
		}

		public Annotation minus(Collection<? extends Annotation> others) {
			return annotation.minus(others);
		}

		public int getDistanceTo(Annotation other) {
			return annotation.getDistanceTo(other);
		}

		public String toUCSC() {
			return annotation.toUCSC();
		}

		public String toBED() {
			return annotation.toBED();
		}

		public String toShortBED() {
			return annotation.toShortBED();
		}

		public String toBEDGraph() {
			return annotation.toBEDGraph();
		}

		public boolean overlaps(Annotation other, int buffer) {
			return annotation.overlaps(other,buffer);
		}

		public boolean overlaps(Collection<? extends Annotation> others,
				int buffer) {
			return annotation.overlaps(others, buffer);
		}

		public boolean overlaps(Annotation other) {
			return annotation.overlaps(other);
		}

		public boolean overlaps(Collection<? extends Annotation> others) {
			return overlaps(others);
		}


		public boolean overlaps(Annotation other, int buffer,boolean considerOrientation) {
			return annotation.overlaps(other, buffer, considerOrientation);
		}

		public boolean overlapsStranded(Annotation other) {
			return annotation.overlapsStranded(other);
		}

		public int getOverlap(Annotation other) {
			return annotation.getOverlap(other);
		}

		public boolean contains(Annotation other) {
			return annotation.contains(other);
		}

		public Annotation union(Annotation other) {
			return annotation.union(other);
		}

		public Annotation intersect(Annotation other) {
			return annotation.intersect(other);
		}

		public List<Annotation> intersect(List<? extends Annotation> annotations) {
			return intersect(annotations);
		}

		public int compareToAnnotation(Annotation b) {
			return annotation.compareToAnnotation(b);
		}

		public void stitchTo(Annotation next) {
			stitchTo(next);
		}

		public String toCigar() {
			return annotation.toCigar();
		}

		public boolean fullyContains(Annotation annotation) {
			return annotation.fullyContains(annotation);
		}

		public Annotation complement() {
			return annotation.complement();
		}

		public boolean equals(Annotation other, boolean useOrientation) {
			return annotation.equals(other, useOrientation);
		}

		public Collection<? extends Annotation> getSpliceConnections() {
			return annotation.getSpliceConnections();
		}

		@Override
		public String toBED(int r, int g, int b) {
			throw new UnsupportedOperationException("TODO");
		}

		@Override
		public boolean overlaps(Annotation other, boolean considerOrientation) {
			// TODO Auto-generated method stub
			return false;
		}

		@Override
		public int getMidpoint() {
			//get midpoint
			int mid=length()/2;
			//convert to reference space
			return getReferenceCoordinateAtPosition(mid);
		}

		@Override
		public String getFullInfoString() {
			// TODO Auto-generated method stub
			return null;
		}
		
	}
}
