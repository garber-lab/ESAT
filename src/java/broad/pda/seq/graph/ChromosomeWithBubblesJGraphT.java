/**
 * 
 */
package broad.pda.seq.graph;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.Iterator;
import java.util.TreeSet;

import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.Gene;
import nextgen.core.annotation.Annotation.Strand;

import org.apache.commons.collections15.Predicate;
import org.apache.log4j.Logger;
import org.broad.igv.sam.Alignment;
import broad.core.datastructures.IntervalTree.Node;
import broad.core.datastructures.ReversibleIterator;

import org.jgrapht.alg.ConnectivityInspector;
import org.jgrapht.alg.KShortestPaths;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.DirectedSubgraph;
import org.jgrapht.util.VertexPair;
import org.jgrapht.DirectedGraph;
import org.jgrapht.GraphPath;

import broad.core.annotation.BasicLightweightAnnotation;
import broad.core.datastructures.IntervalTree;
import broad.core.sequence.Sequence;
import broad.core.util.CollapseByIntersection;
import broad.pda.datastructures.Alignments;
import broad.pda.seq.graph.ChromosomeWithBubblesJGraphT.BubbleEdge;
import broad.pda.seq.segmentation.AlignmentDataModel;
import broad.pda.seq.segmentation.AlignmentDataModelStats;
import broad.pda.seq.segmentation.ReadFilter;

//THESE MUST BE REMOVED EVENTUALLY
import edu.uci.ics.jung.algorithms.shortestpath.DijkstraDistance;
//import edu.uci.ics.jung.algorithms.shortestpath.DijkstraShortestPath;
//import edu.uci.ics.jung.graph.Graph;

/**
 * Implements jung.Graph for Dijkstra (FOR NOW ONLY) 
 *
 */
public class ChromosomeWithBubblesJGraphT extends DefaultDirectedGraph<Annotation, BubbleEdge>{//, Graph<Annotation, BubbleEdge> {

	private static final long serialVersionUID = -6733460455807658588L;
	static Logger logger = Logger.getLogger(ChromosomeWithBubblesJGraphT.class.getName());
	private static final String quote="\"";
	private static final double scaleFactor = 0.0025;
	/*
	 * number of ranking paths between the start vertex and an end vertex
	 */
	private static final int MAX_GENE_PATHS = 3;
	private static final int MAX_ISOFORMS_REPORTED = 10;
	
	private String name;
	private int start;
	private int end;
	private double lambda;
	private double numberOfMarkers;
	private double numberOfReads;
	private IntervalTree<Annotation> vertices;
	private IntervalTree<BubbleEdge> edgeTree;
	private int chromosomeLength;
	private double rpkmConstant;
	private Map<Annotation, Double> localRate;
	
	private static final List<EdgeSourceType> defaultValidTypes = new ArrayList<EdgeSourceType>();
	
	private DijkstraDistance<Annotation, BubbleEdge> distanceAlgorithm; //TODO: This is not currently taking into account edge types nor edge counts Need a filter

	private Map<Annotation, Double> intronScores;
	private  Map<Annotation, Double> vertexCounts;
	
	private IntervalTree<Annotation> gapTree;

	private double spliceWeight=1;

	static {
		defaultValidTypes.add(EdgeSourceType.PAIRED);
		defaultValidTypes.add(EdgeSourceType.SPLICED);
		//defaultTypesForDegreeCounting.add(EdgeSourceType.SPLICED);
	}


	/**
	 * CONSTRUCTORS
	 * TODO: A constructor that instantiates the class with a graph. Application: subgraph
	 */
	public ChromosomeWithBubblesJGraphT ( String name) {
		this(name, 0, Integer.MAX_VALUE);
	}
	
	public ChromosomeWithBubblesJGraphT(String name, int start, int end) {
		super(BubbleEdge.class);
		this.name = name != null ? name.intern() : "noname".intern();
		this.vertexCounts=new TreeMap<Annotation, Double>();
		vertices = new IntervalTree<Annotation>();
		edgeTree = new IntervalTree<BubbleEdge>();
		//exonScores= new TreeMap<Annotation,Double>();
		
		//distanceAlgorithm = new DijkstraShortestPath<Annotation, BubbleEdge>(this, false);
		this.start = start;
		this.end   = end;
		//System.err.println("Initialized graph");
	}

	/**
	 * Constructs a chromosome with bubbles with a given set of nodes and connections between them. The graph will keep track
	 * of how many times a given edge has been added. Nodes in this case are though as exons of putative transcripts.
	 * @param name - Name  of the chromosome.
	 * @param vertices Collection of nodes, note than duplicate nodes will be ignored
	 * @param edges Collection of edges connecting nodes. We assume a close-open( [start, end) ) specification of the annotation, 
	 * 			      thus edge connects two nodes if the start of the edge is node1.getEnd() and the end of the edge is node2.getStart() - 1. 
	 * 				  <b>Edges must be sorted</b>. 
	 * @param edgeCounts Map of edges to their counts
	 * @param nodeCounts Map of vertices to their counts
	 * @param lambda 
	 * @param numMarkers
	 * @param numReads
	 * @param minSplices
	 * @throws IllegalArgumentException if a node or and edge are on a different chromsome than the chromosomeWithBubbles instance.
	 */
	public ChromosomeWithBubblesJGraphT (String name, Collection<? extends Annotation> vertices, Collection<? extends Annotation> edges, Map<? extends Annotation, Double> edgeCounts, Map<? extends Annotation, Double> nodeCounts, double lambda, double numMarkers, double numReads, int minSplices) throws IllegalArgumentException{
		
		/*
		 * Initialize the values
		 */
		this(name);
		this.lambda=lambda;
		this.numberOfMarkers=numMarkers;
		this.numberOfReads=numReads;
		this.vertexCounts=new TreeMap<Annotation, Double>();
		//this.exonScores=nodeCounts;
		
		/*
		 * Intron scores = edge weights. Map from intron (connection intron) to weight.
		 */
		//The hack below is to overcome the imposibility of adding a new item to a map<? extends Annotation
		intronScores = new TreeMap<Annotation, Double>(); 
		if(edgeCounts!=null){
			for(Annotation intron: edgeCounts.keySet()) {
				intronScores.put(intron, edgeCounts.get(intron));
			}
		}
		else{
			for(Annotation intron: edges) {
				intronScores.put(intron, 1.0);
			}
		}
		
		logger.debug("starting graph construction");
		
		// ADDING VERTICES TO GRAPH
		for(Annotation vertex : vertices) {
			/* 
			 * Check that the chromosome of the vertex is the same as the chromosome in question 
			 */
			if(!name.equals(vertex.getChr())) {
				throw new IllegalArgumentException("Trying to add node " + vertex.toUCSC() + " to a chromosome with bubbles " + name + " chromosome must be the same");
			}
			// Internal counter to check number of vertices
			//counter++;
			if(nodeCounts!=null){
				vertexCounts.put(vertex, nodeCounts.get(vertex));
			}
			//if(counter%1 ==0){System.err.println(node);}
			addVertex(vertex);
		}
		logger.debug("gone through all nodes");
		
		//ADDING EDGES TO GRAPH
		for(Annotation edge : edges) {
			/* 
			 * Check that the chromosome of the edge is the same as the chromosome in question 
			 */
			if(! name.equals(edge.getChr())) {
				throw new IllegalArgumentException("Trying to add edge " + edge.toUCSC() + " to a chromosome with bubbles " + name + " chromosome must be the same");
			}
			
			/*
			 * Get the splice count for this annotation/intron
			 * Will be null if no edge from the intron. 
			 */
			double edgeCount = getSpliceCount(edge);
			//System.err.println("Edge " + edge.toUCSC() + " intronScores is null? " +( intronScores==null )+ " if not null, does it contain edge? " + (intronScores != null && intronScores.containsKey(edge) ? intronScores.get(edge) : "edge not in map") );
	
			//Do not add edge, if edgeCount is less than the minimum threshold on splice counts
			if(intronScores != null && !intronScores.isEmpty() && edgeCount <= minSplices) {
				continue;
			}
			
			//FIND ALL VERTICES TOUCHING THE EDGE ON THE LEFT
			List<Annotation> lefts  = findAbuttingLeftVertices(edge);
			if(lefts.isEmpty()){
				//throw new IllegalStateException("Can't find left node for " + edge.toUCSC());
				//System.err.println("WARN: Can't find left node for " + edge.toUCSC());
				continue;
			}
			//FIND ALL VERTICES TOUCHING THE EDGE ON THE RIGHT
			List<Annotation>  rights = findAbuttingRightVertices(edge);
			if(rights.isEmpty()){
				//throw new IllegalStateException("Can't find rigth node for " + edge.toUCSC());
				//System.err.println("WARN: Can't find rigth node for " + edge.toUCSC());
				continue;
			}
			
			//FOR ALL THE "FROM" VERTICES FOR THIS EDGE
			for(Annotation left : lefts) {
				//FOR ALL THE "TO" VERTICES FOR THIS EDGE
				for(Annotation right : rights) {
					//System.err.println("Adding edge " + edge.toUCSC() + " between " + left.toUCSC() + " and " + right.toUCSC());
					BubbleEdge e = edge.getOrientation().equals(Strand.NEGATIVE) ? getEdge(right, left) : getEdge(left, right); 
					boolean couldAdd = false;
					//IF AN EDGE DOES NOT EXIST ALREADY
					if(e == null) {
						if(edge.getOrientation().equals(Strand.NEGATIVE)) {
							e = new BubbleEdge(edge, right, left, new Double(getSpliceCount(edge)).floatValue()); 
							e.setSpliceEdgeCount(edgeCount);
							couldAdd = addEdge(right, left,e);
							e.setVertexCount(right, getCount(right));
							e.setVertexCount(left, getCount(left));
							
						} else {
							e = new BubbleEdge(edge, left, right, new Double(edgeCount).floatValue());
							e.setSpliceEdgeCount(edgeCount);
							couldAdd = addEdge(left, right, e);
							e.setVertexCount(right, getCount(right));
							e.setVertexCount(left, getCount(left));
						}
					} else {
						logger.info("Something is not right, edge " + e.connection.toUCSC() + " was already in the graph. with count " + intronScores.get(e.connection));
						e.incrementCount();
						couldAdd= true;
					}
					//System.err.println("Added edge " + edge.toUCSC() + " between " + left.toUCSC() + " and " + right.toUCSC());
					if (!couldAdd) {
						logger.error("Adding of edge " + edge + " between " + left + " and " + right + " failed");
					}
				}
			}
		}
		System.err.println("gone through all edges");
	}

	/**
	 * Constructs a chromosome with bubbles with a given set of nodes and connections between them. The graph will keep track
	 * of how many times a given edge has been added. Nodes in this case are though as exons of putative transcripts.
	 * @param name - Name  of the chromosome.
	 * @param vertices Collection of nodes, note than duplicate nodes will be ignored
	 * @param edges Collection of edges connecting nodes. We assume a close-open( [start, end) ) specification of the annotation, 
	 * 			      thus edge connects two nodes if the start of the edge is node1.getEnd() and the end of the edge is node2.getStart() - 1. 
	 * 				  <b>Edges must be sorted</b>. 
	 * @throws IllegalArgumentException if a node or and edge are on a different chromsome than the chromosomeWithBubbles instance.
	 */
	public ChromosomeWithBubblesJGraphT (String name, Collection<? extends Annotation> vertices, Collection<? extends Annotation> edges) throws IllegalArgumentException{
		
		this(name, vertices, edges, 0, 0, 0, 0);
	}

	/**
	 * Constructs a chromosome with bubbles with a given set of nodes, edges.
	 * @param name Name  of the chromosome.
	 * @param nodes Collection of nodes, note than duplicate nodes will be ignored
	 * @param edges Collection of edges connecting nodes. We assume a close-open( [start, end) ) specification of the annotation, 
	 * 			      thus edge connects two nodes if the start of the edge is node1.getEnd() and the end of the edge is node2.getStart() - 1. 
	 * 				  <b>Edges must be sorted</b>. 
	 * @param lambda
	 * @param numMarkers
	 * @param numReads
	 * @param minEdgeSupport
	 */
	public ChromosomeWithBubblesJGraphT(String name, Map<? extends Annotation, Double> nodes, Map<BubbleEdge, Double> edges, double lambda, double numMarkers, double numReads, double minEdgeSupport){
		this(name);
		this.lambda=lambda;
		this.numberOfMarkers=numMarkers;
		this.numberOfReads=numReads;
		//this.exonScores=nodes;
		this.intronScores=makeIntronScoresFromBE(edges);
		
		for(Annotation node : nodes.keySet()) {
			
			if(! name.equals(node.getChr())) {
				throw new IllegalArgumentException("Trying to add node " + node.toUCSC() + " to a chromosome with bubbles " + name + " chromosome must be the same");
			}
			//node.setScore(nodes.get(node)); //TODO Orphans might not have scores
			addVertex(node);
			vertexCounts.put(node, nodes.get(node));
		}
		
		for(BubbleEdge edge: edges.keySet()){
			//edge.setCount(edges.get(edge));
			//System.err.println("Edge: " + edge.getConnection().toUCSC() + " type " + edge.getType() + " count " + edge.count + " map count " + edges.get(edge));
			if(edge.getType().ordinal() == EdgeSourceType.SPLICED.ordinal()){
				edge.setSpliceEdgeCount(this.getSpliceCount(edge.getConnection()));
			}
			if(edge.getType().ordinal() != EdgeSourceType.SPLICED.ordinal() || edge.getSplicedCounts() > minEdgeSupport) {
				edge.setParent(this);
				//double nodeCount=nodes.get(edge.getLeftNode())+nodes.get(edge.getRightNode());
				//edge.setNodeCount(nodeCount);
				//edge.setNodeCount(edge.getLeftNode(), this.getCount(edge.getLeftNode()));
				//edge.setNodeCount(edge.getRightNode(), this.getCount(edge.getRightNode()));
				addEdge(edge.getLeftNode(), edge.getRightNode(), edge);
			}
		}
	}
	
	/**
	 * //??? This constructor is not used. And does not return the correct paths
	 * @param name
	 * @param nodes
	 * @param edges
	 * @param fromBE		//??? //NOT USED
	 * @param lambda
	 * @param numMarkers
	 * @param numReads
	 */
	/*public ChromosomeWithBubblesJGraphT(String name, Collection<? extends Annotation> nodes, Collection<BubbleEdge> edges, boolean fromBE, double lambda, double numMarkers, double numReads){
		this(name);
		this.lambda=lambda;
		this.numberOfMarkers=numMarkers;
		
		this.numberOfReads=numReads;
		for(Annotation node : nodes) {
			if(! name.equals(node.getChromosome())) {
				throw new IllegalArgumentException("Trying to add node " + node.toUCSC() + " to a chromosome with bubbles " + name + " chromosome must be the same");
			}
			addVertex(node);
		}
		
		for(BubbleEdge edge: edges){
			addEdge(edge.getLeftNode(), edge.getRightNode(), edge);
		}
	}*/
	
	/**
	 * ???
	 * This will always return a Null Pointer Exception
	 * @param name
	 * @param nodes
	 * @param edges
	 * @param lambda
	 * @param numMarkers
	 * @param numReads
	 * @param minSplices
	 * @throws IllegalArgumentException
	 */
	public ChromosomeWithBubblesJGraphT (String name, Collection<? extends Annotation> nodes, Collection<? extends Annotation> edges, double lambda, double numMarkers, double numReads, int minSplices) throws IllegalArgumentException{
		this(name, nodes, edges, null, null, lambda, numMarkers, numReads, minSplices);
	}
	
	
	/**
	 * Adds an edge e between vertices v1 and v2
	 * @param e
	 * @param v1
	 * @param v2
	 * @return
	 */
	public boolean addEdge(Annotation v1, Annotation v2, BubbleEdge e) {
		boolean success = super.addEdge(v1, v2,e);
		//logger.info("Did it successfully add edge "+e.getConnection()+success);
		
		if(success) {
			edgeTree.put(e.getConnection().getStart(), e.getConnection().getEnd(),e);
			e.setParent(this);
		}
		
		return success;
	}
	
	/**
	 * Adds a new edge between v1 and v2 
	 * @param e
	 * @param v1
	 * @param v2
	 * @return
	 */
	public boolean addNewEdge(BubbleEdge e, Annotation v1, Annotation v2){
		//check if edge already exist
		boolean success = false;
		BubbleEdge previous=getEdge(v1, v2);
		if(previous==null){
			previous=getEdge(v1,v2);
		}
		
		if(previous==null){
			success = addEdge(v1, v2, e);
		}
	
		return previous != null || success;
	}

	/**
	 * Adds vertex v to the graph
	 * @param v
	 */
	public boolean addVertex(Annotation v) {
		boolean success = super.addVertex(v);
		if(success) {
			vertices.put(v.getStart(), v.getEnd(), v);
		}
		return success;
	}
	
	
	/**
	 * Returns a map of the intron to its score from a map of bubbleEdge to score. Only attributes scores if edge is spliced.
	 * @param edges2
	 * @return
	 */
	private Map<Annotation, Double> makeIntronScoresFromBE(Map<BubbleEdge, Double> edges2) {
		Map<Annotation, Double> rtrn=new TreeMap<Annotation, Double>();
		
		for(BubbleEdge edge: edges2.keySet()){
			if(edge.getType().equals(EdgeSourceType.SPLICED)){
				rtrn.put(edge.getConnection(), edges2.get(edge));
			}
		}
		
		return rtrn;
	}

	/**
	 * Returns the splice count for a given intron
	 * @param intron intron for which to return the splice counts
	 * @return splice counts
	 */
	private double getSpliceCount(Annotation intron){
		
		if(this.intronScores==null || !this.intronScores.containsKey(intron)){
			return 0.0;
		}
		else{
			return this.intronScores.get(intron);
		}
	}

	/**
	 * Returns all vertices on the left of the given edge.
	 * @param edge given edge
	 * @return All vertices touching edge on the left
	 */
	private List<Annotation> findAbuttingLeftVertices(Annotation edge) {
		//Node<Annotation> candidateNode = vertices.max(edge.getStart() , edge.getStart()+1);
		List<Annotation> abuttingvertices = new ArrayList<Annotation>();

		Iterator<Node<Annotation>> iter=vertices.overlappers(edge.getStart()-1, edge.getStart());
		
		while(iter.hasNext() ){
			//test if any match exactly, if so set it
			Annotation leftTemp=iter.next().getValue();
			if(leftTemp.getEnd() == edge.getStart() ){
				abuttingvertices.add(leftTemp);
			}
		}
		
		return abuttingvertices;
	}
	
	/**
	 * Returns all vertices on the right of the given edge.
	 * @param edge given edge
	 * @return All vertices touching edge on the right
	 */
	private List<Annotation> findAbuttingRightVertices(Annotation edge) {
		//Node<Annotation> candidateNode = vertices.min(edge.getEnd() , edge.getEnd()+1);
		List<Annotation> abuttingvertices = new ArrayList<Annotation>();

		Iterator<Node<Annotation>> iter=vertices.overlappers(edge.getEnd(), edge.getEnd()+1);
		while(iter.hasNext()){
			//test if any match exactly, if so set it
			Annotation rightTemp=iter.next().getValue();
			if(rightTemp.getStart() == edge.getEnd() ){
				abuttingvertices.add(rightTemp);
			}
		}		
		return abuttingvertices;
	}

	/**
	 * Returns a collection of paired end edges
	 * @return
	 */
	/*private Collection<BubbleEdge> getPairedEndEdges() {
		Collection<BubbleEdge> pairedEdges = new ArrayList<BubbleEdge>();
		Iterator<BubbleEdge> connectionIt = edgeTree.valueIterator();
		while(connectionIt.hasNext()) {
			BubbleEdge e = connectionIt.next();
			//System.err.println("Edge " + e.getConnection().toUCSC() + " type: " + e.getType());
			if(e.getType().equals(EdgeSourceType.PAIRED)) {
				//System.err.println("Added");
				pairedEdges.add(e);
			}	
		}
		return pairedEdges;
	}*/

	public Collection<Gene> getGenePaths(double alpha) {
		return getGenePaths(alpha, 0);
	}

	public Collection<Gene> getGenePaths(double alpha, double minSpliceFrequency) {
		Collection<Gene> genes=new TreeSet<Gene>();
		Collection<Path> paths=getPaths(alpha, minSpliceFrequency);
		for(Path path: paths){
			Path p = (Path) path;
			genes.add(p.toGene());
		}
		return genes;
	}
	
	/**
	 * Returns the collection of Paths through the graph. Will call getAllPaths with 0 minimum splicing frequency.
	 * @return collection of paths
	 */
	public Collection<Path> getAllPaths() {
		return getAllPaths(0);
	} 
	
	/**
	 * Returns the collection of Paths through the graph with the minimum splice frequency
	 * @param minSpliceFrequency minimum splice frequency for the edges
	 * @return
	 */
	public Collection<Path> getAllPaths(double minSpliceFrequency) {
		Collection<Path> rtrn=new TreeSet<Path>();
		
		rtrn.addAll(getOrphanPaths());//??
		rtrn.addAll(getPaths(1, minSpliceFrequency)); 
		
		return rtrn;
	}
	
	/**
	 * Returns all orphan nodes in graph (Vertices without any edges orphan vertices)
	 * @return
	 */
	public Collection<Gene> getOrphanNodes(){
		Collection<Gene> rtrn=new TreeSet<Gene>();
		
		Iterator<Annotation> iter=vertices.valueIterator();
		while(iter.hasNext()){
			Annotation align=iter.next();
			if(this.getIncidentEdges(align).size()==0){rtrn.add(new Gene(align));}
		}
		
		return rtrn;
	}
	
	/**
	 * Returns all orphan nodes in region (Vertices without any edges orphan vertices)
	 * @param region
	 * @return
	 */
	public Collection<Gene> getOrphanNodes(Alignments region){
		Collection<Gene> rtrn=new TreeSet<Gene>();
		
		Iterator<Node<Annotation>> iter=vertices.overlappers(region.getStart(), region.getEnd());
		while(iter.hasNext()){
			Annotation align=iter.next().getValue();
			if(this.getIncidentEdges(align).size()==0){rtrn.add(new Gene(align));}
		}
		
		return rtrn;
	}


	/**
	 * Returns all orphan paths in graph (self edges for orphan vertices)
	 * @return
	 */
	public Collection<Path> getOrphanPaths(){
		Collection<Path> rtrn=new TreeSet<Path>();
		
		Iterator<Annotation> iter=vertices.valueIterator();
		//Iterate through all vertices
		while(iter.hasNext()){
			Annotation align=iter.next();
			//if no edges to or from this vertex
			if(this.getIncidentEdges(align).size()==0){
				//Form new edge
				BubbleEdge edge=new BubbleEdge(align);
				edge.setParent(this);
				Path path=new Path(this.spliceWeight);
				//edge.setNodeCount(align, getCount(align));
				edge.setSpliceEdgeCount(0);
				path.addEdge(edge);
				rtrn.add(path);
			}
		}
		
		return rtrn;
	}

	/**
	 * Returns a collection of paths with alpha
	 * @param alpha
	 * @return
	 */
	public Collection<Path> getPaths(double alpha) {
		return getPaths(alpha, 0);
	}
	
	/**
	 * Returns a collection of paths with alpha and minimum spicing frequency
	 * @param alpha
	 * @param minSpliceFrequency	minimum splicing frequency
	 * @return
	 */
	public Collection<Path> getPaths(double alpha, double minSpliceFrequency) {
		Collection<Path> allPaths = new TreeSet<Path> ();
		//System.out.println("GETTING SOURCES");
		
		/*
		 * GET ALL SOURCES AND SINKS FOR THE GRAPH
		 */
//		Collection<Annotation> sources = getSourceVertices(getValidTypes()/*defaultTypesForDegreeCounting*/); // 06/01/2010 paired ends were ignored when looking for sources
//		System.out.println("# Sources: "+sources.size());
//		Collection<Annotation> allSinks = getTargetVertices(getValidTypes());
//		System.out.println("# Targets: "+allSinks.size());
		
		/*
		 * 1. GET ALL CONNECTED COMPONENTS
		 */
		ConnectivityInspector<Annotation, BubbleEdge> CI = new ConnectivityInspector<Annotation, BubbleEdge>(this);
		List<Set<Annotation>> comps = CI.connectedSets();
		//logger.info("# Connected Components: "+comps.size());
		/*
		 * FOR EACH COMPONENT
		 */
		for(Set<Annotation> component: comps){
			
			/*
			 * IF THE SIZE OF THE COMPONENT IS > 1
			 */
			if(component.size()>1){
				/*
				 * 2. 	MAKE 2 SUBGRAPHS: 1 for each orientation
				 * 		BASED ON THE ORIENTATIONS OF THE EDGES, REMOVE EDGES FROM THE SUBGRAPHS
				 */
				Set<BubbleEdge> subedges = this.getEdges(component);
				DirectedSubgraph<Annotation, BubbleEdge>[] strandSpecificSubgraphs = new DirectedSubgraph[2];
				
				// the [0] is the positive subgraph
				strandSpecificSubgraphs[0] = new DirectedSubgraph<Annotation, BubbleEdge>(this,component,subedges);
				// the [1] is the negative subgraph
				strandSpecificSubgraphs[1] = new DirectedSubgraph<Annotation, BubbleEdge>(this,component,subedges);				
				
				for(BubbleEdge e:subedges){
					if(e.inReversedOrientation()){
						strandSpecificSubgraphs[0].removeEdge(e);
					}
					else{
						strandSpecificSubgraphs[1].removeEdge(e);
					}
				}
				
				/*
				 * 3. FOR EACH POSITIVE AND NEGATIVE SUBGRAPHS
				 */
				for(DirectedSubgraph<Annotation, BubbleEdge> ds: strandSpecificSubgraphs){
					
					List<EdgeSourceType> validTypes = new ArrayList<EdgeSourceType>();
					validTypes.add(EdgeSourceType.SPLICED);
					ds = filterOutLessFreqEdges(ds, 0.1, validTypes);
					/*
					 * 		FIND THE SOURCES IN THE SUBGRAPH
					 * 		FIND THE TARGETS IN THE SUBGRAPH
					 */
					Collection<Annotation> subSources = getSourceVertices(getValidTypes(),ds);	
					Collection<Annotation> subSinks = getSinkVertices(getValidTypes(),ds);
					Collection<Path> subgraphPaths = new TreeSet<Path> ();
					/*
					 *  	FOR EACH SOURCE IN SUBGRAPH
					 */
					for(Annotation s : subSources ) {
						//logger.trace("Getting paths starting in " + s.toUCSC() + " ");
						Collection<Annotation> targets = findReachableTargetVertices(s, (new ArrayList<Annotation>()), subSinks, (new ArrayList<Annotation>()),0, ds);
						//logger.trace("to " + targets.size() + " targets ");
						long start=System.currentTimeMillis();
						for(Annotation e : targets ) {
							//logger.trace(" to target: " + e.toUCSC());
							if(s.equals(e)){
								//BubbleEdge ed = this.getEdge(s, e);
								//System.err.println(ed.getType().toString());
								logger.trace("Source and target vertex are the same. Do not process");//???
							}
							else{
								Collection<Path> sourcePaths = getPaths(s, e, alpha, minSpliceFrequency, ds);
								subgraphPaths.addAll(sourcePaths);
							}
						}
						long end = System.currentTimeMillis();
						//logger.debug("Time taken : "+(end-start));
					}
					
					ArrayList<Path> pathsToAdd = new ArrayList<Path>();
					for(Path p: subgraphPaths){
			            
						int i=0;
						boolean pathAdded = false;
			            //If there is no other path in the pathsToAddList
						for(;i<pathsToAdd.size();i++){
							
							// case when the new path is better than the path Py stored at index y
							if (p.getScore() > pathsToAdd.get(i).getScore()) {
								pathsToAdd.add(i, p);
			                    pathAdded = true;
			                    // ensures max size limit is not exceeded.
			                    if (pathsToAdd.size() > MAX_ISOFORMS_REPORTED) {
			                    	pathsToAdd.remove(MAX_ISOFORMS_REPORTED);
			                    }
			                    break;
			                }
							// case when the new path is of the same score as the path Py stored at index y
			                if (p.getScore() == pathsToAdd.get(i).getScore()) {
			                	pathsToAdd.add(i+1, p);
			                    pathAdded = true;
			                    // ensures max size limit is not exceeded.
			                    if (pathsToAdd.size() > MAX_ISOFORMS_REPORTED) {
			                    	pathsToAdd.remove(MAX_ISOFORMS_REPORTED);
			                    }
			                    break;
			                }
						}
						if (pathsToAdd.size() < MAX_ISOFORMS_REPORTED) {
			                    // the new path is inserted at the end of the list.
			                	pathsToAdd.add(p);
			                    pathAdded = true;
						} else {
			                    // max size limit is reached -> end of the loop over the
			                    // paths elements of the list at vertex v.
			                    break;
			            }
					}
					
					allPaths.addAll(pathsToAdd);
				}
			}
		}
		
		logger.info("There are "+allPaths.size()+" spliced paths");
		
		return allPaths;
	}
	

	/**
	 * Removes all outgoing edges from a node with less than X% of total number of reads emitted from node
	 * @param g
	 * @param alpha
	 * @return
	 */
	private DirectedSubgraph<Annotation, BubbleEdge> filterOutLessFreqEdges(DirectedSubgraph<Annotation,BubbleEdge> g, double alpha, Collection<EdgeSourceType> validTypes) {
		//For every node in the graph
		for(Annotation vertex : g.vertexSet()) {
			Collection<BubbleEdge> edges = g.outgoingEdgesOf(vertex);
			double total = 0.0;
			for(BubbleEdge e : edges) {
				total += e.getSplicedCounts();
			}
			for(BubbleEdge e : edges) {
				if(isEdgeValid(e,  alpha, validTypes)) {
					if(((double)e.getSplicedCounts()/(double)total)<alpha){
						g.removeEdge(e);
					}
				}
			}
			
		}
		return g;			
	}
	/**
	 * Returns all edges connecting the specified vertices 
	 * @param vertices
	 * @return
	 */
	public Set<BubbleEdge> getEdges(Collection<Annotation> vertices){
	
		Set<BubbleEdge> edges = new TreeSet<BubbleEdge>();
		
		for(Annotation v:vertices){
			edges.addAll(this.getIncidentEdges(v)); 
		}
		
		return edges;
	}
	
	public Collection<Annotation> findReachableTargetVertices(Annotation source, Collection<Annotation> alreadyTraversed,Collection<Annotation> sinkVertices, Collection<Annotation> reachableSinks, int counter, DirectedGraph<Annotation, BubbleEdge> g) {
		if(counter>150){
			System.err.println("Counter>150");
			return reachableSinks;
		}
		
		List<EdgeSourceType> validTypes = new ArrayList<EdgeSourceType>();
		validTypes.add(EdgeSourceType.SPLICED);
		//System.out.print(source.toUCSC()+" is connected to :");
		Collection<BubbleEdge> forwardEdges = getValidOutgoingEdges(source, 0.0, validTypes, g);
		if(sinkVertices.contains(source)){
		//if(outOrientedDegree(source, validTypes) == 0 && inOrientedDegree(source, validTypes) >= 1 && !isSpanningGap(source.getStart(), source.getEnd())){
			reachableSinks.add(source);
			//System.out.println(source.toUCSC()+" which is a sink.");
			return reachableSinks; 
		}
		else{
			if(forwardEdges.size() > 0) {
				for(BubbleEdge e : forwardEdges) {
					Annotation opposite = getOppositeVertex(source, e);
					if(!alreadyTraversed.contains(opposite) ){
						//System.out.print(opposite.toUCSC()+" , ");
						reachableSinks = findReachableTargetVertices(opposite, alreadyTraversed, sinkVertices,reachableSinks, counter++,g);
						alreadyTraversed.add(opposite);
					}
				}
			}
			return reachableSinks;
		}
	}

	/**
	 * Returns a collection of all paths starting from startNode with alpha and minimum splicing frequency
	 * @param startNode
	 * @param alpha
	 * @param minSpliceFreq	minimum splicing frequency
	 * @return
	 */
	public List<Path> getPaths(Annotation startNode, Annotation endNode, double alpha, double minSpliceFreq, DirectedGraph<Annotation, BubbleEdge> g){
		return getPaths(startNode,endNode, alpha, getValidTypes(), minSpliceFreq, g);
	}
	
	/**
	 * Returns a collection of all paths starting from startNode having the specified valid edge types
	 * @param startNode
	 * @param alpha
	 * @param validEdgeTypes
	 * @param minSpliceFreq	minimum splicing frequency
	 * @return
	 */
	public List<Path> getPaths(Annotation startNode, Annotation endNode,double alpha, Collection<EdgeSourceType> validEdgeTypes, double minSpliceFreq, DirectedGraph<Annotation, BubbleEdge> g) {
		return getPaths(startNode, endNode, alpha, validEdgeTypes, new HashMap<Annotation, Collection<Path>>(),1, 0, minSpliceFreq, g);
		//return getPaths(startNode, alpha, validEdgeTypes, new HashMap<Annotation, Collection<Path>>(),1);

	}
	
	/**
	 * Returns a collection of paths starting with given vertex, for the specified allowed edge types, given nodes visited so far, 
	 * the number of paths predicted so far, the direction and minimum splicing frequency
	 * @param startNode Starting vertex
	 * @param alpha
	 * @param validEdgeTypes
	 * @param nodesVisited	A map of vertices to a collection of paths for that node (this far)
	 * @param predictedPathsSoFar
	 * @param direction
	 * @param minFreq
	 * @return
	 */
	private List<Path> getPaths(Annotation startNode, Annotation endNode,double alpha, Collection<EdgeSourceType> validEdgeTypes, Map<Annotation, Collection<Path>> nodesVisited, int predictedPathsSoFar, int direction, double minFreq, DirectedGraph<Annotation, BubbleEdge> g) {
		List<Path> allPaths = new ArrayList<Path>();

		/*
		 * maximum number of edges of the calculated paths.
		 */
		int nMaxHops = Integer.MAX_VALUE;
		KShortestPaths<Annotation,BubbleEdge> KPaths = new KShortestPaths<Annotation,BubbleEdge>(g,startNode,MAX_GENE_PATHS,nMaxHops);
		allPaths = PWrapper(KPaths.getPaths(endNode));
		//System.err.println("\tAll paths: "+allPaths.size());
		return allPaths;
	}

	
	private List<Path> PWrapper(List<GraphPath<Annotation,BubbleEdge>> paths){
		
		List<Path> pw = new ArrayList<Path>();
		for(GraphPath<Annotation,BubbleEdge> p :paths){
			Path x = new Path();
			boolean sameOrientation = checkOrientation(p.getEdgeList());
			if(sameOrientation){
				x.addEdges(p.getEdgeList());
				pw.add(x);
			}
			else{
				//nothing
			}
			
		}
		if(pw.size()==0){
			return null;
		}
		return pw;
	}
	
	public boolean checkOrientation(Collection<BubbleEdge> edges) {
		boolean ok = true;
		Iterator<BubbleEdge> edgeIt = edges.iterator();
		if(edges.isEmpty()){
			ok=true;
		}
		else{
			String orientation=edgeIt.next().getOrientation();
			for(BubbleEdge edge: edges){
				if("*".equalsIgnoreCase(orientation)|| "*".equalsIgnoreCase(edge.getOrientation()) || orientation.equalsIgnoreCase(edge.getOrientation())){
					//value of ok doesnt change
				}
				else{
					ok = false;
					System.err.println("Throws "+orientation+" "+edge.getOrientation()+"Path can not have conflicting orientations. Not adding Path");
				}
			}
		}
		return ok;
	}
	/**
	 * Returns the opposite vertex of the edge e from the vertex v
	 * @param v
	 * @param e
	 * @return
	 */
	public Annotation getOppositeVertex(Annotation v, BubbleEdge e) {
		
		Annotation source = this.getEdgeSource(e);
		Annotation target = this.getEdgeTarget(e);
	        if (v.equals(source)) {
	            return target;
	        } else if (v.equals(target)) {
	            return source;
	        } else {
	            throw new IllegalArgumentException("no such vertex");
	        }
	}
	
	
	/**
	 * DO WE NEED THIS????
	 * Returns collection of genes between start and end
	 * @param start
	 * @param end
	 * @return
	 */
	public Collection<Gene> getPaths(int start, int end) {
		return getPaths(start, end, 1, true);
	}
	
	public Collection<Gene> getPaths(int start, int end, boolean forward) {
		return getPaths(start, end, 0, true);
	}
	
	public Collection<Gene> getPaths(int start, int end, double alpha, boolean forward) {		
		return forward? getPaths(start, end, alpha, getValidTypes()) : getBackwardPaths(start, end, alpha, getValidTypes());
	}
	
	public Collection<Gene> getPaths(int start, int end, double alpha, Collection<EdgeSourceType> validTypes) {
		List<Gene> paths = new ArrayList<Gene>();
		
		Iterator<Node<Annotation>> verticesOverlappingIt = vertices.overlappers(start, end);
		List<Annotation> verticesTraversedList = new ArrayList<Annotation>();
		
		while(verticesOverlappingIt.hasNext()) {
			Annotation v= verticesOverlappingIt.next().getValue();
			//System.err.println("New Vertex " + v.toUCSC() + " traversed list " +  verticesTraversedList + " can be reached? " + canBeReachedFrom(v, verticesTraversedList));
			if(v.getEnd() != end && !canBeReachedFrom(v, verticesTraversedList)) {
				// block one is from start to v 
				Alignments firstBlock = new  Alignments(name, start, v.getEnd()); // Change to support vertices not overlapping edges MG
				// find other blocks by following each edge to the next splice
				Collection<BubbleEdge> edges = getForwardEdges(v, alpha, validTypes);
				for(BubbleEdge edge : edges) {
					Annotation nextV = getOppositeVertex(v, edge);//getOpposite(v, edge);
					Alignments orientedFirstBlock  = new Alignments(name, start, v.getEnd()); // Change to support vertices not overlapping edges MG
					orientedFirstBlock.setOrientation(edge.getOrientation());
					Collection<Gene> nextPaths =  getPaths(nextV.getStart(), nextV.getStart() + (end - start) - firstBlock.length(), alpha, validTypes); //Recursion
					// append blocks found following each edge to the nascent transcript.
					if(nextPaths.size() == 0) {
						Alignments lastBlock = new Alignments(name, nextV.getStart(), nextV.getStart() + (end - start) - firstBlock.length());
						lastBlock.setOrientation(edge.getOrientation());
						//if(!lastBlock.getOrientation().equals(orientedFirstBlock)) {
							//System.err.println("WARN: Making an inconsistent gene, exon " + orientedFirstBlock.toUCSC() + "("+orientedFirstBlock.getOrientation()+") wheras " + lastBlock.toUCSC() + "("+lastBlock.getOrientation()+")");
						//}
						List<Alignments> exons = new ArrayList<Alignments>();
						exons.add(orientedFirstBlock);
						exons.add(lastBlock);
						Gene prependedPath = new Gene(exons);
						prependedPath.setOrientation(orientedFirstBlock.getOrientation());
						paths.add(prependedPath);
					} else {
						for(Gene newPath : nextPaths ) {
							Collection<? extends Annotation> newPathIntrons = newPath.getIntronSet()	;// Change to support vertices not overlapping edges MG
							if(newPathIntrons.size() == 0) {
								verticesTraversedList.add(new BasicLightweightAnnotation(name, newPath.getStart(), newPath.getStart()+1));
							} else  {
							
								Iterator<? extends Annotation> newPathIntronIt = newPathIntrons.iterator();
								while(newPathIntronIt.hasNext()) {
									Annotation intron = newPathIntronIt.next();
									verticesTraversedList.add(new BasicLightweightAnnotation(name, intron.getStart()-1, intron.getStart()));// Change to support vertices not overlapping edges MG
									verticesTraversedList.add(new BasicLightweightAnnotation(name, intron.getEnd(), intron.getEnd() + 1));// Change to support vertices not overlapping edges MG
								}
							}
							List<Annotation> exons = new ArrayList<Annotation>();
							//if(!newPath.getOrientation().equals(orientedFirstBlock)) {
							//	System.err.println("WARN: Making an inconsistent gene, exon " + orientedFirstBlock.toUCSC() + "("+orientedFirstBlock.getOrientation()+") wheras " + newPath.toString() );
							//}							
							exons.add(orientedFirstBlock); 
							exons.addAll(newPath.getExonSet()); //Changed to support vertices not overlapping edges
							Gene prependedPath = new Gene(exons);
							prependedPath.setOrientation(orientedFirstBlock.getOrientation());
							paths.add(prependedPath);
						}
					}
				}

			}
		}
		if(paths.isEmpty()){
			paths.add(new Gene(name, start, end));
		}
		
		//System.err.println(new Gene(name, start, end));
		
		return paths;
	}

	public Collection<Gene> getBackwardPaths(int start, int end, double alpha, Collection<EdgeSourceType> validTypes) {
		List<Gene> paths = new ArrayList<Gene>();
		
		Iterator<Node<Annotation>> verticesOverlappingIt = vertices.overlappers(start, end);
		List<Annotation> overlappingVertices = new ArrayList<Annotation>();
		while(verticesOverlappingIt.hasNext()) {
			overlappingVertices.add(verticesOverlappingIt.next().getValue());
		}
		Iterator<Annotation> vertexIt = new ReversibleIterator<Annotation>(overlappingVertices.listIterator(), false);
		
		List<Annotation> verticesTraversedList = new ArrayList<Annotation>();
		
		while(vertexIt.hasNext() ) {
			Annotation v= vertexIt.next();
			//System.err.println("New Vertex " + v.toUCSC() + " traversed list " +  verticesTraversedList + " can be reached? " + canBeReachedFrom(v, verticesTraversedList));
			if(v.getStart() != start && !canBeReachedFrom(v, verticesTraversedList)) {
				// block one is from start to v 
				Alignments lastBlock = new  Alignments(name, v.getStart(), end); // Change to support vertices not overlappin edges MG
				// find other blocks by following each edge to the next splice
				Collection<BubbleEdge> edges = getBackwardEdges(v, alpha, validTypes);
				//System.err.println("Vertex " + v.toUCSC() );
				for(BubbleEdge edge : edges) {
					Annotation prevV = getOppositeVertex(v, edge);
					//System.err.println("\tEdge, " + edge.getConnection().toUCSC() + edge.getConnection().getOrientation() + " Oposite v " + prevV.toUCSC());
					Alignments orientedLastBlock  = new Alignments(lastBlock); // Change to support vertices not overlappin edges MG
					if(edge.inReversedOrientation()) {
						orientedLastBlock.setOrientation(edge.inReversedOrientation() ? "-" : "+");
					}
					Collection<Gene> nextPaths =  getBackwardPaths(prevV.getStart()  - (end - start) - lastBlock.length(), prevV.getEnd(), alpha, validTypes); //Recursion
					// append blocks found following each edge to the nascent transcript.
					if(nextPaths.size() == 0) {
						Alignments firstBlock = new Alignments(name, prevV.getStart()  - (end - start) - lastBlock.length(), prevV.getEnd());
						firstBlock.setOrientation(edge.inReversedOrientation() ? "-" : "+");
						List<Alignments> exons = new ArrayList<Alignments>();
						exons.add(firstBlock);
						exons.add(orientedLastBlock);
						Gene prependedPath = new Gene(exons);
						paths.add(prependedPath);
					} else {
						for(Gene newPath : nextPaths ) {
							Collection<? extends Annotation> newPathIntrons = newPath.getIntronSet()	;// Change to support vertices not overlapping edges MG
							if(newPathIntrons.size() == 0) {
								verticesTraversedList.add(new BasicLightweightAnnotation(name, newPath.getStart(), newPath.getStart()+1));
							} else  {
							
								Iterator<? extends Annotation> newPathIntronIt = newPathIntrons.iterator();
								while(newPathIntronIt.hasNext()) {
									Annotation intron = newPathIntronIt.next();
									verticesTraversedList.add(new BasicLightweightAnnotation(name, intron.getStart()-1, intron.getStart()));// Change to support vertices not overlapping edges MG
									verticesTraversedList.add(new BasicLightweightAnnotation(name, intron.getEnd(), intron.getEnd() + 1));// Change to support vertices not overlapping edges MG
								}
							}
							List<Annotation> exons = new ArrayList<Annotation>();
							exons.addAll(newPath.getExonSet()); //Changed to support vertices not overlapping edges
							exons.add(lastBlock); 
							Gene prependedPath = new Gene(exons);
							paths.add(prependedPath);
						}
					}
					//System.err.println("\tpaths: " + paths);
				}

			}
		}
		if(paths.isEmpty()){
			paths.add(new Gene(name, start, end));
		}
		
		//System.err.println(new Gene(name, start, end));
		
		return paths;
	}
	
	public Collection<BubbleEdge> getBackwardEdges(Annotation vertex, double alpha, Collection<EdgeSourceType> validTypes) {
		Collection<BubbleEdge> backwardEdges = new ArrayList<BubbleEdge>();
		Collection<BubbleEdge> edges = getIncidentEdges(vertex);
		for(BubbleEdge e : edges) {
			if(getOppositeVertex(vertex, e).getEnd() < vertex.getStart() && isEdgeValid(e, alpha, validTypes)) {
				backwardEdges.add(e);
			}
		}

		return backwardEdges;			
	}

	private boolean canBeReachedFrom(Annotation v, List<Annotation> vertices) {
		boolean canBeReached = vertices.contains(v);
		
		Iterator<Annotation> vIt = vertices.iterator();
		while(!canBeReached && vIt.hasNext()) {
			Annotation traversedVertex = vIt.next();
			if(v.compareTo(traversedVertex) >= 0) {
				Number d = distanceAlgorithm.getDistance(v, traversedVertex);
				canBeReached = d != null;
			}
		}
		return canBeReached;
	}


	

	/**
	 * Returns a list of all types of allowed edge types: Paired & Spliced
	 * @return list of allowed edge types
	 */
	private List<EdgeSourceType> getValidTypes() {
		return defaultValidTypes;
	
	}

	/**
	 * Sets the splice weight for this graph
	 * @param spliceWeight Value for splice weight
	 */
	public void setSpliceWeight(double spliceWeight) {this.spliceWeight=spliceWeight;}
	
	/**
	 * Returns a collection of all incoming and outgoing edges of specified vertex
	 * @param vertex
	 * @return collection of all edges for specified vertex
	 */
	 public Collection<BubbleEdge> getIncidentEdges(Annotation vertex){
	        if (!containsVertex(vertex))
	            return null;
	        
	        Collection<BubbleEdge> incident_edges = new HashSet<BubbleEdge>();
	        incident_edges.addAll(incomingEdgesOf(vertex));
	        incident_edges.addAll(outgoingEdgesOf(vertex));
	        return Collections.unmodifiableCollection(incident_edges);
	 }
	 
	 /**
	  * Returns a collection of all the possible source/start vertices for the specified graph.
	  * @param validTypes List of allowed edge types.
	  * @return All vertices having allowed edge types.
	  */
	 public Collection<Annotation> getSourceVertices(Collection<EdgeSourceType> validTypes,DirectedGraph<Annotation, BubbleEdge> g) {
		 	//GET ALL VERTICES
			Collection<Annotation> vertices = g.vertexSet();
			Collection<Annotation> sources = new TreeSet<Annotation>();
			//FOR EACH VERTEX
			for(Annotation v: vertices) {
				//System.err.println(v.toUCSC()+" "+inOrientedDegree(v, validTypes)+" "+outOrientedDegree(v, validTypes) + " " + isSpanningGap(v.getStart(), v.getEnd()));
				//if(inOrientedDegree(v, validTypes) == 0 && outOrientedDegree(v, validTypes) > 0 && !isSpanningGap(v.getStart(), v.getEnd())) {
				/*
				 * IF THERE ARE NO INCOMING EDGES TO THIS VERTEX 
				 * AND
				 * THERE IS AT LEAST 1 OUTGOING EDGE FROM THIS VERTEX
				 * AND
				 * THE VERTEX IS NOT AN INTRON
				 */
				if(inOrientedDegree(v, validTypes,g) == 0 && outOrientedDegree(v, validTypes,g) >= 1 && !isSpanningGap(v.getStart(), v.getEnd())) {
					//System.out.println("adding source");
					sources.add(v);
					
				}
			}
			
			return sources;
		}
		
	 /**
	  * Returns a collection of all the possible target/end vertices for this graph.
	  * @param validTypes List of allowed edge types.
	  * @return All vertices having allowed edge types.
	  */
	 public Collection<Annotation> getSinkVertices(Collection<EdgeSourceType> validTypes,DirectedGraph<Annotation, BubbleEdge> g) {
		 	
		 	//GET ALL VERTICES
			Collection<Annotation> vertices = g.vertexSet();
			Collection<Annotation> targets = new TreeSet<Annotation>();
			//FOR EACH VERTEX
			for(Annotation v: vertices) {
				//System.err.println(v.toUCSC()+" "+inOrientedDegree(v, validTypes)+" "+outOrientedDegree(v, validTypes) + " " + isSpanningGap(v.getStart(), v.getEnd()));
				/*
				 * IF THERE ARE NO OUTGOING EDGES FROM THIS VERTEX 
				 * AND
				 * THERE IS AT LEAST 1 INCOMING EDGE TO THIS VERTEX
				 */
				if(outOrientedDegree(v, validTypes, g) == 0 && inOrientedDegree(v, validTypes,g) >= 1 && !isSpanningGap(v.getStart(), v.getEnd())) {
					//System.out.println("adding source");
					targets.add(v);
				}
			}
			
			return targets;
	}
	 
	public Collection<BubbleEdge> getForwardEdges(Annotation vertex, double alpha, Collection<EdgeSourceType> validTypes) {
			Collection<BubbleEdge> forwardEdges = new ArrayList<BubbleEdge>();
			Collection<BubbleEdge> edges = getIncidentEdges(vertex);
			for(BubbleEdge e : edges) {
				if(getOppositeVertex(vertex, e).getStart() > vertex.getEnd() && isEdgeValid(e,  alpha, validTypes)) {
					forwardEdges.add(e);
				}
			}

			return forwardEdges;			
	}
		 
	 public Collection<BubbleEdge> getValidOutgoingEdges(Annotation vertex, double alpha, Collection<EdgeSourceType> validTypes, DirectedGraph<Annotation, BubbleEdge> g) {
			Collection<BubbleEdge> outEdges = new ArrayList<BubbleEdge>();
			Collection<BubbleEdge> edges = g.outgoingEdgesOf(vertex);//this.getOutOrUnorientedEdges(vertex);
			for(BubbleEdge e : edges) {
				if(isEdgeValid(e,  alpha, validTypes)) {
					outEdges.add(e);
				}
			}

			return outEdges;			
	}
	 /**
	  * Returns the in degree of the given vertex for a list of allowed edge types
	  * @param v specific vertex
	  * @param validTypes Collections of allowed edge types
	  * @return in degree of given vertex
	  */
	 public int inOrientedDegree(Annotation v, Collection<EdgeSourceType> validTypes,DirectedGraph<Annotation, BubbleEdge> g) {
			int inDegree = 0;
			Collection<BubbleEdge> inEdges = g.incomingEdgesOf(v);
			for(BubbleEdge e : inEdges) {
				if(isEdgeValid(e, 1.1, validTypes)) {
					inDegree++;
				}
			}
			//System.err.println("#incident " + inDegree(v) + " #filtered indegree " + inDegree);
			return inDegree;
		}

	 /**
	  * Returns the out degree of the given vertex for a list of allowed edge types
	  * @param v specific vertex
	  * @param validTypes Collections of allowed edge types
	  * @return out degree of given vertex
	  */ 
	public int outOrientedDegree(Annotation v, Collection<EdgeSourceType> validTypes, DirectedGraph<Annotation, BubbleEdge> g) {
		int outDegree = 0;
		Collection<BubbleEdge> outEdges = g.outgoingEdgesOf(v);
		for(BubbleEdge e : outEdges) {
			if(isEdgeValid(e, 1.1, validTypes)) {
				outDegree++;
			}
		}	
		return outDegree;
	}

	/**
	 * Returns true if this region spans an intron. Also, returns true if the position starts in gap but goes past gap.
	 * @param start
	 * @param end
	 * @return
	 */
	public boolean isSpanningGap(int start, int end){
		boolean result = false;
		boolean isInAllIntrons = isInAllGaps(start, end);
		//System.err.print("Testing ("+start+","+end+")" + " is in all introns? " + isInAllIntrons + " ");
		if(isInAllIntrons) {
			Iterator<Node<BubbleEdge>> overlappers=edgeTree.overlappers(start, end);
			while(overlappers.hasNext() && !result){
				Node<BubbleEdge> next=overlappers.next();
				result = end > next.getEnd() || start < next.getStart();
				//System.err.print("["+next.getValue().getConnection().toUCSC()+"-"+result);
			}
		}
		//System.err.println(" ");
		return result;
	}

	/**
	 * Returns true if position starts in gap but goes past gap
	 * @param start
	 * @param end
	 * @return
	 */
	public boolean isGap(int start, int end){
		/*Iterator<Node<Annotation>> overlappers=edges.overlappers(start, end);
		boolean isGap = overlappers.hasNext();
		while(overlappers.hasNext() && isGap){
			Node<Annotation> next=overlappers.next();
			isGap = start >= next.getStart() && end<=next.getEnd();
		}
		return isGap;*/
		return false;
	}
	
	/**
	 * Returns true if the start and end positions are all found in gaps (introns??)
	 * @param start
	 * @param end
	 * @return
	 */
	public boolean isInAllGaps(int start, int end) {
		/*
		 * If gap tree has not been loaded yet
		 * 		Load it
		 */
		if( gapTree == null || gapTree.isEmpty() ) {
			loadChromosomeGapTree();
		}
		//return true if the positions overlap some gaps in tree.
		return gapTree.overlappers(start, end).hasNext();
	}
	
	/**
	 *	Loads the interval tree for all gaps for this chromosome.
	 */
	private void loadChromosomeGapTree() {
		
		/*
		 * GET ALL EDGES OF THE GRAPH
		 */
		//System.err.println("loading chromosome gap tree");
		Collection<Alignments> edgeAlignments = new ArrayList<Alignments>(edgeSet().size());
		Iterator<BubbleEdge> edgeIt = edgeTree.valueIterator();
		//System.err.println("\tEdges:");
		/*
		 * GET THE ALIGNMENTS BETWEEN THE START AND END OF EACH EDGE I.E. INTRONS
		 */
		while(edgeIt.hasNext()) {
			Annotation e = edgeIt.next().getConnection();
			/*BubbleEdge be = findEdge(e);
			if(be.count>5){
				System.err.println("\t\t"+e.toUCSC());
			*/	
				edgeAlignments.add(new Alignments(e.getChr(), e.getStart(), e.getEnd()));
			//}
			
		}
		/*
		 * COLLAPSE ALL THE INTRONS INTO A UNIQUE SET
		 */
		Set<Annotation> colapsedEdges = CollapseByIntersection.collapseByIntersection(edgeAlignments, true);
		//System.err.println("\tCollapsed edges:");
		gapTree = new IntervalTree<Annotation>();
		for(Annotation ce : colapsedEdges) {
			//System.err.println("\t\tAdding gap: " + ce.toUCSC());
			gapTree.put(ce.getStart(), ce.getEnd(), ce);
		}
		//System.err.println("Done loading chromosome gap tree");
	}
	
	/**
	 * Removes all nodes that span a gap
	 */
	public void removeSpanningGapNodes() {
		List<Annotation> vertexList = vertices.toList();
		for(Annotation v : vertexList) {;
			if(isSpanningGap(v.getStart(), v.getEnd())) {
				//logger.trace("Removing gap spanning vertex " + v.toUCSC());
				boolean couldRemove = removeVertex(v);
				if(!couldRemove) {logger.error("Could not remove gap spanning vertex " + v.toUCSC());} 
			}
		}
		
	}
	
	/**
	 * Removes vertex from graph
	 * @param vertex
	 */
	public boolean removeVertex(Annotation vertex) {
		Annotation removed = vertices.remove(vertex.getStart(), vertex.getEnd());
		boolean allGood = true;
		if(removed != null) {
			Collection<BubbleEdge> incidentEdges = new ArrayList<BubbleEdge>(incomingEdgesOf(vertex));
			for(BubbleEdge incidentEdge : incidentEdges) {
				//System.err.println("\tIn edge " + incidentEdge.getConnection().toUCSC());
				boolean couldRemove = removeEdge(incidentEdge);
				if(allGood) {allGood = couldRemove;}
				if(!allGood) {System.err.println("Could not remove in edge " + incidentEdge.getConnection().toUCSC());}
				
			}
			Collection<BubbleEdge> outEdges = new ArrayList<BubbleEdge>(outgoingEdgesOf(vertex));
			for(BubbleEdge outEdge : outEdges) {
				//System.err.println("\tOut edge " + outEdge.getConnection().toUCSC());
				boolean couldRemove = removeEdge(outEdge);
				if(allGood) {allGood = couldRemove;}
				if(!allGood) {System.err.println("Could not remove out edge " + outEdge.getConnection().toUCSC());}
			}
			boolean couldRemove = super.removeVertex(vertex);
			//System.err.println("Could remove " + vertex.toUCSC() + " from parent graph? " + couldRemove + " does parent graph contains verte? " + containsVertex(vertex));
			if(allGood) {allGood = couldRemove;}
		}
		
		return allGood;
	}
	
	/**
	 * Removes edge from graph
	 * @param e
	 * @return 
	 */
	public boolean removeEdge(BubbleEdge e) {
		BubbleEdge removed = edgeTree.remove(e.getConnection().getStart(), e.getConnection().getEnd());
		boolean ret = removed != null;
		if(ret || containsEdge(e)) { //The parent graph may contain duplicate edges according to us. So even if we do not find it, the parent may still have it.
			ret = super.removeEdge(e); 
		}
		
		return ret;
	}
	
	/**
	 * Returns edge
	 * @param e
	 * @return
	 */
	/*private BubbleEdge findEdge(Annotation e) {
		Annotation es = new BasicLightweightAnnotation(name, e.getStart()-1, e.getStart());
		Annotation ee = new BasicLightweightAnnotation(name, e.getEnd(), e.getEnd() + 1);
		
		return e.inReversedOrientation() ?  getEdge(ee, es) : getEdge(es, ee);
	}*/

	/**
	 * Returns true if the specified edge of a valid edge type.
	 * @param e specific edge
	 * @param alpha //???
	 * @param validTypes ollections of allowed edge types
	 * @return
	 * //NOT NEEDED : ALPHA
	 */
	private boolean isEdgeValid(BubbleEdge e, double alpha, Collection<EdgeSourceType> validTypes) {
		//System.err.print("\tEdge "+e.getCount()+" "+minReads + " type <" + e.getType()+">");
		
		boolean valid = validTypes == null || validTypes.isEmpty();
		
		if(!valid) {
			Iterator<EdgeSourceType> typeIt = validTypes.iterator();
			while(!valid && typeIt.hasNext()) {
				EdgeSourceType type = typeIt.next();
				valid = type.ordinal() == e.getType().ordinal();
			}
		}
		
		return valid && e.getPvalue() <= alpha;

	}

	/**
	 * Returns a collection of edges that are of the allowed edge types.
	 * @param edgeList
	 * @return filtered edges
	 */
	/*private Collection<BubbleEdge> filterValidEdges(Collection<BubbleEdge> edgeList) {
		Collection<BubbleEdge> filtered = new TreeSet<BubbleEdge>();
		for(BubbleEdge e : edgeList) {
			//System.out.println("Considering edge");
			if(isEdgeValid(e, 1.1, getValidTypes())) {
				//System.out.println("Edge is valid");
				filtered.add(e);
			}
		}
		
		return filtered;	
	}*/

	/**
	 * Returns a collection of edges that originate at the specified vertex (if oriented) 
	 * or have the specified vertex as an end point
	 * @param v Specified vertex
	 * @return collection of edges
	 */
	public Collection<BubbleEdge> getOutOrUnorientedEdges(Annotation v) {
		Collection<BubbleEdge> rtrn = new TreeSet<BubbleEdge>();
		/*
		 * ERROR MESSAGES IF NO OUTGOING EDGES FOR v
		 */
		if(outgoingEdgesOf(v) == null) {
			System.err.println("V " + v.toUCSC() + " has NULL out edges! ");
		} 
		/*
		 * ADD ALL OUTGOING EDGES TO COLLECTION
		 */
		else {
			//System.out.println("\tgetOutEdges(v): " + getOutEdges(v).size());
			rtrn.addAll(outgoingEdgesOf(v));
		}

		/*
		 * IF ORIENTATION OF ALL IN EDGES IS "*" (UNORIENTED), ADD THEM
		 */
		Collection<BubbleEdge> inEdges = incomingEdgesOf(v);
		for(BubbleEdge e : inEdges) {
			if("*".equals(e.getOrientation()) ){
					rtrn.add(e);
			}
		}
		return rtrn;
	}
	

	/**
	 * Returns a collection of edges that terminate at the specified vertex (if oriented) 
	 * or have the specified vertex as an end point
	 * @param v Specified vertex
	 * @return collection of edges
	 */
	public Collection<BubbleEdge> getInOrUnorientedEdges(Annotation v) {
		Collection<BubbleEdge> rtrn = new TreeSet<BubbleEdge>();
		/*
		 * ERROR MESSAGES IF NO INCOMING EDGES FOR v
		 */
		if(incomingEdgesOf(v) == null) {
			System.err.println("V " + v.toUCSC() + " has NULL in edges! ");
		} else {
			rtrn.addAll(incomingEdgesOf(v));
		}

		Collection<BubbleEdge> outEdges = outgoingEdgesOf(v);
		for(BubbleEdge e : outEdges) {
			if("*".equals(e.getOrientation()) ){
					rtrn.add(e);
			}
		}
		return rtrn;
	}

	
	/**
	 * Returns a collection of edges of the allowed edge types that originate at the specified vertex (if oriented) 
	 * or have the specified vertex as an end point
	 * @param v Specified vertex
	 * @return collection of edges
	 */
	/*private Collection<BubbleEdge> getOutOrUnorientedValidEdges(Annotation v) {
		return filterValidEdges(getOutOrUnorientedEdges(v));
		
	}*/
	
	
	
	/**
	 * Creates bubbles that connect the splice junctions specified in the data model. Nodes in the graphs (bubbles) are single base pairs 
	 * indicating the splice junction.
	 * @param model
	 * @param minAligmentGap
	 */
	

	
	/**
	 * Adds an edge to the graph, adding missing vertices to the graph.
	 * @param intron
	 * @param count
	 */
	public void addEdgeAddingMissingVertices(Annotation intron, float count) {
		Annotation start = new BasicLightweightAnnotation(intron.getChr(), intron.getStart() - 1, intron.getStart());//TODO: revise, should be [getStart()-1, getStart())
		Annotation end   = new BasicLightweightAnnotation(intron.getChr(), intron.getEnd(), intron.getEnd() + 1);
		
		
		BubbleEdge edge = intron.getOrientation().equals(Strand.NEGATIVE) ? getEdge(end, start) : getEdge(start, end);
		if(edge != null) {
			edge.setCount(edge.getAllCounts() + count); //update count if duplicated edge is added.
		} else {
			if(!containsVertex(start)) {
				addVertex(start);
			} 
			
			if(!containsVertex(end) ){
				addVertex(end);
			}
			
			if(intron.getOrientation().equals(Strand.NEGATIVE)) {
				edge = new BubbleEdge(intron, end, start, count);
				edge.setParent(this);
				edge.setVertexCount(end, this.getCount(end));
				edge.setVertexCount(start, this.getCount(start));
				//edge.setSpliceEdgeCount(this.getSpliceCount(intron));//MODIFIED BY MANUEL 1/4/2011 This looked like a bug. The splice count should be the one passed in as the paramter PLEASE REVIEW THIS
				intronScores.put(edge.connection, (double)count);
				addEdge(end, start, edge);
			} else {
				edge = new BubbleEdge(intron, start, end, count);
				edge.setParent(this);
				edge.setVertexCount(start, this.getCount(start));
				edge.setVertexCount(end, this.getCount(end));
				//edge.setSpliceEdgeCount(this.getSpliceCount(intron));//MODIFIED BY MANUEL 1/4/2011 This looked like a bug. The splice count should be the one passed in as the paramter PLEASE REVIEW THIS
				addEdge(start, end, edge);
				intronScores.put(edge.connection, (double)count);
			}		
		}
	}

	/**
	 * Returns edge between start and end
	 * @param start
	 * @param end
	 * @return
	 */
	public Annotation findEdge(int start, int end) {
		Node<BubbleEdge> edgeNode = edgeTree.find(start, end);
		return edgeNode != null ? edgeNode.getValue().getConnection() : null;
	} 
	
	/**
	 * Returns a list of all vertices overlapping region between start and end
	 * @param start
	 * @param end
	 * @return
	 */
	public List<Annotation> getOverlappingVertices(int start, int end) {
		ArrayList<Annotation> overlappers = new ArrayList<Annotation>();
		Iterator<Node<Annotation>> overlapperIt = vertices.overlappers(start, end);
		
		while(overlapperIt.hasNext()) {
			overlappers.add(overlapperIt.next().getValue());
		}
		return overlappers;
	}
	
	/**
	 * Creates an edge of type between v1 and v2 for connection alignment with weight
	 * @param connection
	 * @param v1
	 * @param v2
	 * @param count
	 * @param type
	 * @return
	 */
	public BubbleEdge createEdge(Alignments connection, Annotation v1, Annotation v2, Integer count, EdgeSourceType type) {
		BubbleEdge edge=new BubbleEdge(connection, v1, v2, count, type);
		edge.setParent(this);
		return edge;
	}
	
	/**
	 * Returns the pair of vertices for a specified edge
	 * @param edge
	 * @return
	 */
	public VertexPair<Annotation> getNodePair(BubbleEdge edge) {
		return (new VertexPair<Annotation>(edge.getLeftNode(),edge.getRightNode()));
	}

	public WindowIterator iterator(int windowSize, int start, int overlap) {
		return new WindowIterator(this, windowSize, start, overlap);
	}
	
	public WindowIterator iterator(int windowSize, int start, int overlap, boolean skipGapSpanningWindows) {
		return new WindowIterator(this, windowSize, start, overlap, skipGapSpanningWindows);
	}
	
	public WindowIterator iterator(int windowSize, int start) {
		return new WindowIterator(this, windowSize, start, windowSize - 1, true);
	}
	
	public double getLocalRate(Annotation exon){
		if(this.localRate==null || !this.localRate.containsKey(exon)){return this.lambda;}
		return Math.max(this.localRate.get(exon), lambda);
	}
	
	public double getLambda(){return this.lambda;}
	public int getRelativeGenomeLength(){return chromosomeLength;}	
	public double getNumberOfMarkers() {return this.numberOfMarkers;}
	public double getNumberOfReads() {return this.numberOfReads;}	
	public String getName() {return name;}	
	public double getRPKMConstant() {return this.rpkmConstant;}	
	public Map<Annotation, Double> getLocalRate(){return this.localRate;}	
	public int getEnd() {return end;}
	public int getStart() {return start;}
	
	/**
	 * Sets the start of this graph
	 * @param start
	 */
	public void setStart(int start){
		this.start = start;
	}
	
	/**
	 * Sets the end of this graph
	 * @param end
	 */
	public void setEnd(int end){
		this.end = end;
	}
	
	public double getSpliceWeight(){return this.spliceWeight;}
	public int getChromosomeLength() {return chromosomeLength;}

	public void setLocalRate(Collection<Path> paths){
		this.localRate=new TreeMap<Annotation,Double>();
		for(Path path: paths){
			Collection<? extends Annotation> exons=path.toGene().getExonSet();
			for(Annotation exon: exons){
				this.localRate.put(exon, path.getLocalLambda());
			}
		}
	}
	
	public void setLocalRate(Map<Annotation, Double> localRateMap) {
		this.localRate=localRateMap;		
	}
	
	public void setRPKMConstant(double num) {
		this.rpkmConstant=num;
	}
	
	public void setExonScores(Annotation node, Double score) {
		this.vertexCounts.put(node, score);
		
	}
	
	public void setChromosomeLength(int chromosomeLength) {
		this.chromosomeLength = chromosomeLength;
	}

	/**
	 * Write the graph as a DOT format
	 * @param save
	 * @param setUpShortExonLabels
	 * @throws IOException
	 */
	public void writeGraph(String save, boolean setUpShortExonLabels)throws IOException{
		Alignments chrRegion=new Alignments(this.getName(), 0, this.end);
		writeGraph(save, chrRegion, setUpShortExonLabels);
	}
	
	public void writeGraph(String save)throws IOException{
		writeGraph(save, true);
	}
	
	public void writeGraph(String save, Alignments region) throws IOException {
		writeGraph(save, region, true);
	}
	
	public void writeGraph(String save, Alignments region, boolean setUpShortExonLabels) throws IOException {

		writeGraph(save, region ,setUpShortExonLabels, getName()+"_"+lambda+"_"+this.numberOfMarkers+"_"+this.numberOfReads);
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
	//	Map<BubbleEdge, Pair<Annotation>> edgeMap= edges;
		
		//start by writing an invisible connection between every node and every node in sorted order
	//	IntervalTree<Annotation> nodeTree=this.vertices;
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
		//int exonCounter=0;
		while(iter.hasNext()){
			Annotation exon=iter.next().getValue();
			double exonCount=this.getCount(exon);
			double localRate=this.getLocalRate(exon);
			if(exon.overlaps(region)){
				double normWidth=exon.length()*scaleFactor;
				writer.write(quote+exon.toUCSC()+quote+" [width="+normWidth+", fixedsize=true, label="+exonCount+", comment="+localRate); 
				writer.write("];\n");
				//exonCounter++;
			}
		}
		
		//This defines actual edges
		Collection<BubbleEdge> edges = edgeSet();
		for(BubbleEdge edge: edges){
			VertexPair<Annotation> nodes=getNodePair(edge);
			//System.err.println("WRITEGRAPH - edge" + edge.getConnection().toUCSC());
			Annotation first=nodes.getFirst();
			Annotation second=nodes.getSecond();
			if(first.overlaps(region) || second.overlaps(region)){
				writer.write(quote+first.toUCSC()+quote+"->"+quote+second.toUCSC()+quote);
				writer.write("[label="+edge.getAllCounts()+", comment="+quote+edge.getType()+quote);
				if(edge.getType().equals(EdgeSourceType.PAIRED)){
					writer.write(", style=dotted");
				}
				//writer.write(", minlen="+(edge.getConnection().length()*(scaleFactor/3))); //TODO: Add scaling for intron sizes (default no scaling)
				writer.write("];\n");
			}
		}
					
		writer.write("}");
		
		writer.close();
	}

	private double getCount(Annotation align) {
		//TODO Consider putting back in
		/*if(this.exonScores==null){return 0.0;}
		if(!this.exonScores.containsKey(align)){return 0.0;}
		else{return this.exonScores.get(align);}*/
		
		if(this.vertexCounts==null){return 0.0;}
		else if(!this.vertexCounts.containsKey(align)){return 0.0;}
		else{return this.vertexCounts.get(align);}
	}

	/**
	 * This class represents the Edges of the graph
	 *
	 */
	public static class BubbleEdge implements Comparable<BubbleEdge>{
		
		private double count;
		private double pairedEndCounts;
		private double splicedCount;
		private double pvalue;
		private Annotation connection;
		private Annotation v1;
		private Annotation v2;
		private EdgeSourceType type;
		private Collection<Alignment> pairedEndSupport;
		private ChromosomeWithBubblesJGraphT parentGraph;
				
		/**
		 * This constructor creates an instance of a spliced edge (by default) from v1 to v2.
		 * Need to store source and destination vertices to be able to compare edges.
		 * @param connection
		 * @param v1
		 * @param v2
		 * @param connectionCount
		 */
		BubbleEdge(Annotation connection, Annotation v1, Annotation v2, double connectionCount) {
			//BY DEFAULT: The edge is spliced
			this(connection, v1, v2, connectionCount, EdgeSourceType.SPLICED);
		}


		/**
		 * This constructor creates an instance of an edge of "type" from v1 to v2.
		 * Need to store source and destination vertices to be able to compare edges.
		 * @param connection
		 * @param v1
		 * @param v2
		 * @param connectionCount
		 * @param type of connection.
		 */
		public BubbleEdge(Annotation connection, Annotation v1, Annotation v2, double connectionCount, EdgeSourceType type) {
			this.connection = connection;
			this.count = connectionCount;
			this.v1=(v1);
			this.v2=(v2);
			this.type = type;
			pairedEndSupport = new ArrayList<Alignment>(); //CANT BE A SET since paired ends may be duplicates
			if(this.connection.getOrientation()==null){
				this.connection.setOrientation(Strand.UNKNOWN);
			}
		}
		
		/**
		 * This function initiates the class instance as a copy of the input BubbleEdge
		 * @param toCopy The Edge to be copied
		 */
		BubbleEdge(BubbleEdge toCopy) {
			this.parentGraph = toCopy.parentGraph;
			this.connection = toCopy.connection;
			this.count = toCopy.count;
			this.v1=(toCopy.v1) ;
			this.v2=(toCopy.v2);
			this.type = toCopy.type;
			pairedEndSupport = toCopy.pairedEndSupport;
		}
		
		/**
		 * This constructor creates a self-edge from startVertex to startVertex with weight 0
		 * @param startVertex
		 */
		public BubbleEdge(Annotation startVertex) {
			this(startVertex, startVertex, startVertex, 0, EdgeSourceType.SELF);
		}
		
		/**
		 * Sets the parent graph of the edge
		 * @param parent Graph object to be set as the parent
		 */
		public void setParent(ChromosomeWithBubblesJGraphT parent) {
			this.parentGraph = parent;
		}
		
		/**
		 * Sets the spliced edge count for this edge.
		 * @param i Counts from spliced reads
		 */
		public void setSpliceEdgeCount(double i) {
			this.splicedCount=i;
		}
		
		/**
		 * Sets the counts for the vertex
		 * @param vertex Vertex for which counts must be set.
		 * @param i Counts for the vertex.
		 */
		public void setVertexCount(Annotation vertex, double i) {
			//System.err.println("node: " +node);
			//System.err.println("node " + node + " vertex counts " + (parentGraph == null ? null : parentGraph.vertexCounts) + " parent " + parentGraph);
			parentGraph.setExonScores(vertex, new Double(i));
			//graph.setNodeCounts(node, i);
		}

		/**
		 * Returns the counts for the vertex
		 * @param node Vertex for which counts must be returned.
		 * @return the counts for the vertex
		 */
		public double getVertexCount(Annotation vertex){
			if(EdgeSourceType.SELF.equals(getType())) {
				return parentGraph.getCount(connection);
			} else {
				return parentGraph.vertexCounts.get(vertex);		
			}
		}
			
		/**
		 * Returns the vertex opposite another vertex across an edge.
		 * @param v vertex at one end of edge
		 * @return vertex at opposite end of the edge
		 */
		public Annotation getOpposite(Annotation v) {
			Annotation opposite = null;
			if(v.equals(v1)) {
				opposite = v2;
			} else if (v.equals(v2)) {
				opposite = v1;
			}
			return opposite;
		}

		/**
		 * Get the type of edge.  
		 * @return type of edge. (SPLICED, PAIRED, PAIRED_SUPPORT, SELF)
		 */
		public EdgeSourceType getType() {
			return type;
		}

		/**
		 * Increments the count by 1
		 */
		public void incrementCount() {
			count++;
		}
		
		/**
		 * Increments the count by some amount
		 * @param amount amount by which to increment the count.
		 */
		public void incrementCount(float amount) {
			count += amount;
		}

		/**
		 * Returns the orientation of the edge.
		 * @return orientatio of edge.
		 */
		public String getOrientation(){return connection.getOrientation().toString();}

		//TODO: Fix the counts
		/**
		 * Returns the total paired end and spliced end counts for this edge.
		 * @return all counts for edge
		 */
		public double getAllCounts() {
			return getPairedEndCounts()+getSplicedCounts();
		}
		
		/**
		 * Returns the total paired end counts for this edge.
		 * @param count paired end counts for edge.
		 */
		public double getPairedEndCounts(){return this.pairedEndCounts;}
		
		/**
		 * Sets the total spliced end counts for this edge.
		 * @param  spliced end counts for edge.
		 */
		public void setPairedEndCounts(double count){this.pairedEndCounts=count;}
		
		/**
		 *  Returns the total spliced end counts for this edge.
		 * @return  spliced end counts for edge.
		 */
		public double getSplicedCounts(){return this.splicedCount;}
		
		/**
		 * Sets the counts for this edge.
		 * @param  counts for edge.
		 */
		public void setCount(double count) {
			this.count = count;
		}
		
		/**
		 * Returns the connection vertex for this edge
		 * @return Connection vertex for this edge
		 */
		public Annotation getConnection() {
			return connection;
		}
		
		/**
		 * Adds a paired end alignment to the collection of paired end alignments for this edge
		 * @param a paired end alignment to add to edge.
		 */
		public void addPairedEndSupport(Alignment a){
			pairedEndSupport.add(a);
		}
		
		/**
		 * Returns the collection of paired end alignments for this edge
		 * @return Collection of paired end alignments for edge
		 */
		public Collection<Alignment> getPairedEndSupport() { return pairedEndSupport;}
		
		/**
		 * Returns true if the connection node is on the negative strand
		 * @return true if connection node is on the negative strand
		 */
		public boolean inReversedOrientation() {return connection.getOrientation().equals(Strand.NEGATIVE);}

		/**
		 * Returns P-value for this edge
		 * @return p-value for edge
		 */
		public double getPvalue() {
			return pvalue;
		}
		
		/** 
		 * Sets the p-value for this edge
		 * @param pvalue p-value for edge
		 */
		public void setPvalue(double pvalue) {
			this.pvalue = pvalue;
		}

		/**
		 * Compares this edge to another edge.
		 */
		public int compareTo(BubbleEdge o) {
			int comparison = o.connection.compareTo(connection);
			if(comparison == 0) {
				if(connection.getOrientation().equals(Strand.NEGATIVE) != o.connection.getOrientation().equals(Strand.NEGATIVE)) {
					comparison = connection.getOrientation().equals(Strand.NEGATIVE) ? -1 : 1;
				}
				if(comparison == 0) {
					comparison = getLeftNode().compareTo(o.getLeftNode());
				}
				
				if(comparison == 0) {
					comparison = getRightNode().compareTo(o.getRightNode());
				}
			}
			return comparison;
		}
		
		/**
		 * Returns a hash code
		 */
		public int hashCode() {
			return (connection.hashCode() + getLeftNode().hashCode() + getRightNode().hashCode())* (connection.getOrientation().equals(Strand.NEGATIVE) ? -1 : 1);
		}
		
		/**
		 * Indicates whether some other object is "equal to" this one.
		 * @param Object to check equality with.
		 * @return true if the object is the same as the o argument
		 */
		public boolean equals(Object o) {
			boolean ret = false;
			if(o instanceof BubbleEdge) {
				ret = compareTo((BubbleEdge) o) == 0;
			}
			return ret;
		}

		/**
		 * Returns the from-vertex for this edge. 
		 * @return The left vertex (from vertex)
		 */
		public Annotation getLeftNode() {
			return v1;
		}

		/**
		 * Returns the to-vertex for this edge. 
		 * @return The right vertex (to vertex)
		 */
		public Annotation getRightNode() {
			return v2;
		}

	}

	public enum EdgeSourceType {
		SPLICED, PAIRED, PAIRED_SUPPORT, SELF; 

	}

	
	/**
	 * This class exercises the minimum length filter on reads
	 *
	 */
	public static class MinReadFilter implements ReadFilter<Alignment> {
		int minLength;
		public MinReadFilter(int minLength) {
			this.minLength = minLength;
		}

		public boolean passes(Alignments read) {
			return read.length() > minLength;
		}

		@Override
		public boolean passes(Alignment read) {
			// TODO Auto-generated method stub
			return false;
		}
		
	}

	/**
	 * Read filter overlapping region
	 * 
	 */
	public static class RegionReadFilter implements ReadFilter<Alignments> {
		int start;
		int end;
		
		public RegionReadFilter(int start, int end) {
			this.start = start;
			this.end  = end;
		}
		
		public boolean passes(Alignments read) {
			return read.getStart() < end && read.getEnd() > start; //If read overlaps region.
		}

	}

	/**
	 * An iterator over the collection of RefSeq genes
	 * 
	 */
	public static class WindowIterator implements Iterator<Collection<Gene>> {
		private ChromosomeWithBubblesJGraphT data;
		private int atPosition;
		private int windowSize;
		private int step;
		private int dataAccess;
		private int numWindowsDone;
		private int winWithNoJumps;
		
		private Collection<Gene> last;
		//private List<Annotation> lastOverlappingVertices;
		
		protected WindowIterator() {
			
		}
		
		public WindowIterator(ChromosomeWithBubblesJGraphT data, int windowSize, int start, int overlap){
			this(data,windowSize,start,overlap,false);
		}
		
		public WindowIterator(ChromosomeWithBubblesJGraphT data, int windowSize, int start, int overlap, boolean skipGapSpanningWindows) {
			this.data = data;
			this.windowSize = windowSize;
			this.atPosition = start;
			this.step = windowSize - overlap;
		}
		public boolean hasNext() {
			return atPosition < data.getEnd();
		}
		
		public int getDataAccess() { return dataAccess;}
		public int getNumWindowsDone() { return numWindowsDone;}
		public int getWinWithNoJumps() { return winWithNoJumps;}

		public Collection<Gene> next() {
			int end   = atPosition + windowSize;
			last = data.getPaths(atPosition, end, 0, true);
			for(Gene window : last) {
				Collection<? extends Annotation> windowIntrons = window.getIntronSet();
					for(Annotation intron: windowIntrons) {
						if(data.findEdge(intron.getStart(), intron.getEnd())==null) {throw new RuntimeException("BAD path intron: \n" + intron + "\npath:\n" + window.toBED());}
					}
			}
			atPosition += step;
			return last;
		}


		public void remove() {
			//
		}
		
		public Collection<Gene> next(int lowerBound, int higherBound) {
			atPosition = higherBound;
			int end = higherBound + windowSize;
			last = data.getPaths(atPosition, end, 0 , true);
			//lastOverlappingVertices = data.getOverlappingVertices(atPosition, end);
			return last;
		}
		
		public void jumpTo(int lowerBound, int higherBound) {
			atPosition = higherBound;
		}
		
	}
	
	/**
	 * FORWARD ITERATOR 
	 */
	public static class ForwardAndBackIterator extends WindowIterator {
		LinkedList<WindowIterator> iteratorQueue;
		
		public ForwardAndBackIterator(ChromosomeWithBubblesJGraphT data, int windowSize,  int overlap, boolean skipGapSpanningWindows) {
			iteratorQueue = new LinkedList<WindowIterator>();
			WindowIterator forwardIt = new WindowIterator(data, windowSize, data.getStart(), overlap, skipGapSpanningWindows);
			iteratorQueue.add(forwardIt);
			WindowIterator backIt = new ReverseWindowIterator(data, windowSize,data.getEnd(), overlap, skipGapSpanningWindows);
			iteratorQueue.add(backIt);
		}

		public boolean hasNext() {
			WindowIterator current = iteratorQueue.peek();
			boolean hasNext = false;
			if(current != null) {
				hasNext = current.hasNext();
				if(!hasNext) {
					iteratorQueue.pop();
					System.err.println("Switching iterator");
					hasNext = hasNext();
				}
			}
			return hasNext;
		}

		public Collection<Gene> next() {
			return hasNext() ? iteratorQueue.peek().next() : null;
		}

		public void remove() {
			// TODO Auto-generated method stub
			
		}
			
		public Collection<Gene> next(int lowerBound, int higherBound) {
			return hasNext() ? iteratorQueue.peek().next(lowerBound, higherBound) : null;
		}
	}

	/**
	 * REVERSE ITERATOR
	 */
	public static class ReverseWindowIterator extends WindowIterator {	
		public ReverseWindowIterator(ChromosomeWithBubblesJGraphT data, int windowSize, int start, int overlap, boolean skipGapSpanningWindows) {
			super(data, windowSize, start, overlap, skipGapSpanningWindows);
			//Move to next position that does not span a gap.
		}
		public boolean hasNext() {
			return super.atPosition > super.data.getStart();
		}
		
		public Collection<Gene> next() {
			int start   = super.atPosition - super.windowSize;
			super.last = super.data.getPaths(start, super.atPosition, 0, false);
			/*System.err.println("At position " + super.atPosition + " paths found: " + super.last.size());
			for(Gene window : super.last) {
				Collection<Alignments> windowIntrons = window.getIntronSet();
					for(Alignments intron: windowIntrons) {
						if(super.data.findEdge(intron.getStart(), intron.getEnd())==null) {throw new RuntimeException("BAD path intron: \n" + intron + "\npath:\n" + window.toBED());}
					}
			}*/
			super.atPosition -= super.step;
			return super.last;
		}


		public void remove() {
			//
		}
		
		public Collection<Gene> next(int lowerBound, int higherBound) {
			super.atPosition = lowerBound;
			int start = lowerBound - super.windowSize;
			super.last = super.data.getPaths(start, super.atPosition, 0, false);
			//super.lastOverlappingVertices = super.data.getOverlappingVertices(start, super.atPosition);
			return super.last;
		}
		
		public void jumpTo(int lowerBound, int higherBound) {
			super.atPosition = lowerBound;
		}
		
	}

}
