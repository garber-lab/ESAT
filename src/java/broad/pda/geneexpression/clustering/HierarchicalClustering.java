package broad.pda.geneexpression.clustering;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Stack;
import java.util.TreeSet;

import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;

import broad.core.datastructures.MatrixWithHeaders;
import broad.core.error.ParseException;
import broad.core.math.Statistics;
import broad.pda.geneexpression.agilent.AgilentUtils;


//Implements agglomerative clustering
//TODO: Should print out a cdt, gtr, and atr file

public class HierarchicalClustering {

	String metric="complete";
	MatrixWithHeaders data;
	ClusterDistanceFunction distanceFunction;
	ArrayList<Cluster>  initialColumnClusters;
	ArrayList<Cluster>  initialRowClusters;
	
	public HierarchicalClustering(MatrixWithHeaders data) {
		this(data, null, null);
	}
	
	public HierarchicalClustering(MatrixWithHeaders data1, Map<String, Collection<String>> columnGroups, Map<String, Collection<String>> rowGroups){
		
		this.data=data1;
		distanceFunction = new PearsonDistance();
		//Initial clusters will be each sample, unless a grouping is provided
		initialColumnClusters = makeInitialClusters(data, columnGroups, true);
		//Initialize clusters
		initialRowClusters = makeInitialClusters(data, rowGroups, false);
		
	}
	
	public void write(String save) throws IOException{
		//write re-ordered
		data.writeGCT(save);
	}
	
	public void setClusterDistanceFunction(ClusterDistanceFunction distanceFunction) throws IllegalArgumentException{
		this.distanceFunction = distanceFunction;
	}
	
	public void setClusterDistanceFunction(String distanceFunctionName) {
		if("euclidean".equalsIgnoreCase(distanceFunctionName)) {
			setClusterDistanceFunction(new EuclideanDistance());
		} else if ("abspearson".equalsIgnoreCase(distanceFunctionName)){
			setClusterDistanceFunction(new AbsolutePearsonFunction());
		}else if ("pearson".equalsIgnoreCase(distanceFunctionName)){
			setClusterDistanceFunction(new PearsonDistance());
		} else {
			throw new IllegalArgumentException ("Invalid distance function name: " + distanceFunction + " only pearson, agspearson and euclidean are supported");
		}
	}
	
	public void setLinkage(String linkage) throws IllegalArgumentException {
		if("complete".equalsIgnoreCase(linkage) || "single".equalsIgnoreCase(linkage) || "average".equalsIgnoreCase(linkage)) {
			this.metric = linkage;
		} else {
			throw new IllegalArgumentException("linkage was " + linkage +" but it must be one of complete, single or average");
		}
	}
	
	public void cluster (boolean clusterRow, boolean clusterColumn) {
		if(clusterColumn){
			System.err.println("Clustering columns");
			//clustering will go through each group pairwise and compute the distance until only 1 group exists
			Cluster columnCluster=createCluster(true);
			
			//Order them
			data=order(columnCluster, data, true);
		}
		
		if(clusterRow){
			System.err.println("Clustering rows");

			//cluster
			Cluster rowCluster=createCluster( false);
			
			//order
			data=order(rowCluster, data, false);
		}
	}
	
	private Cluster createCluster( boolean columns) {
		ArrayList<Cluster> initialClusters = columns? initialColumnClusters : initialRowClusters; 
		ArrayList<Cluster> clusters = initialClusters;

		System.err.print("Computing similarity matrix .. ");
		long pTime = System.nanoTime();
		//TODO: repeace similarityMat with a HashTable it is a much smaller structure when the groups are large 
		// and the number of pairwise distances compared is smaller.
		MatrixWithHeaders similarityMat = makeSimilarityMatrix(columns);
		System.err.println(" done  " +  (System.nanoTime() - pTime)/1000000);
		int level=0;
		while(clusters.size()>1){
			double minDistance=Double.MAX_VALUE;
			int iMin=-1;
			int jMin=-1;
			
			//Matrix similarityMatrix = new Matrix(initialClusters.size(), );
			
			System.err.print("level "+level);
			pTime = System.nanoTime();
			for(int i=0; i<initialClusters.size(); i++){
				for(int j=i; j<initialClusters.size(); j++){
					if(i!=j){
						Cluster cluster1=initialClusters.get(i);
						Cluster cluster2=initialClusters.get(j);
						double p=distance(cluster1, cluster2, similarityMat, columns);
						if(p<minDistance){
							minDistance=p;
							iMin=i;
							jMin=j;
						}
					}
				}
			}
			System.err.println(" done  " +  (System.nanoTime() - pTime)/1000000);
			Cluster cluster1=clusters.get(iMin);
			Cluster cluster2=clusters.get(jMin);
			//System.err.println(minDistance+" "+iMin+" "+jMin+" "+cluster1.getMembers()+" "+cluster2.getMembers());
			
			
			clusters.remove(cluster1);
			clusters.remove(cluster2);
			Cluster merged=merge(cluster1, cluster2);
			merged.setScore(minDistance);
			clusters.add(merged);
			level++;
		}
		
		return clusters.get(0);
	}
	
	private MatrixWithHeaders makeSimilarityMatrix(boolean columns)  throws IllegalStateException{
		List<String> rowColNames = columns? data.getColumnNames() : data.getRowNames();
		MatrixWithHeaders similarityMat = new MatrixWithHeaders(rowColNames, rowColNames);
		/*for(String r : rowColNames) {
			for(String c : rowColNames) {
				if(!c.equals(r)) {
					double[] vals1=getData(data, r, columns);
					double[] vals2=getData(data, c, columns);
					double d=distanceFunction.measure(vals1, vals2);
					similarityMat.set(r,c,d);
				}
			}
		}*/
		
		return similarityMat;
	}
	
	
	private double distance(Cluster cluster1, Cluster cluster2, MatrixWithHeaders similarityMat, boolean columns) throws IllegalStateException{
		double dist = 0;
		
		if(metric.equalsIgnoreCase("complete")){
			for(String member1: cluster1.getMembers()){
				for(String member2: cluster2.getMembers()){
					//double t = System.nanoTime();
					double d=getOrUpdateDistance(member1, member2, similarityMat, columns) ;
					//System.err.println("Too " + (System.nanoTime() - t)/1000000000 + " to compute dist ("+member1+","+member2+")");
					
					if(d>dist){dist=d;}
				}
			}
		}else if(metric.equalsIgnoreCase("single")){
			dist = Double.MAX_VALUE;
			for(String member1: cluster1.getMembers()){
				for(String member2: cluster2.getMembers()){
					double d=getOrUpdateDistance(member1, member2, similarityMat, columns);
					
					if(d<dist){dist=d;}
				}
			}
		}else if(metric.equalsIgnoreCase("average")){
			for(String member1: cluster1.getMembers()){
				for(String member2: cluster2.getMembers()){
					double d=getOrUpdateDistance(member1, member2, similarityMat, columns);
					dist += d;
				}
			}
			dist =  dist/(double)(cluster1.getMembers().size() * cluster2.getMembers().size());
		}
		else{
			throw new IllegalStateException("Linkage was " + metric + " it can only be one of complete, single or average, the HierarchicalCluster object is badly set up.");
		}
		
		return dist;	
		
	}
	
	
	private double getOrUpdateDistance(String member1, String member2, MatrixWithHeaders similarityMat, boolean columns) {
		double dist = similarityMat.get(member1, member2);
		if(dist == 0 && !member1.equals(member2)) {
			//System.err.println("Hummm... computing dist ("+member1+","+member2+")");
			double[] vals1=getData(data, member1, columns);
			double[] vals2=getData(data, member2, columns);
			dist=distanceFunction.measure(vals1, vals2);
			similarityMat.set(member1,member2,dist);
			similarityMat.set(member2, member1,dist);
		}
		
		return dist;
	}

	private double distanceOld(Cluster cluster1, Cluster cluster2, String metric, MatrixWithHeaders data, boolean columns) {
		if(metric.equalsIgnoreCase("complete")){
			double maxDist=-1;
			for(String sample1: cluster1.getMembers()){
				for(String sample2: cluster2.getMembers()){
					double[] vals1=getData(data, sample1, columns);
					double[] vals2=getData(data, sample2, columns);
					//double[] vals1=data.getColumn(sample1);
					//double[] vals2=data.getColumn(sample2);
					//double d=1-Statistics.pearsonDistance(vals1, vals2);
					double d=distanceFunction.measure(vals1, vals2);
					if(d>maxDist){maxDist=d;}
				}
				return maxDist;
			}
			
		}else if(metric.equalsIgnoreCase("single")){
			double minDist=Double.MAX_VALUE;
			for(String sample1: cluster1.getMembers()){
				for(String sample2: cluster2.getMembers()){
					double[] vals1=getData(data, sample1, columns);
					double[] vals2=getData(data, sample2, columns);
					//double[] vals1=data.getColumn(sample1);
					//double[] vals2=data.getColumn(sample2);
					//double d=1-Statistics.pearsonDistance(vals1, vals2);
					double d=distanceFunction.measure(vals1, vals2);
					if(d<minDist){minDist=d;}
				}
				return minDist;
			}
			
		}else if(metric.equalsIgnoreCase("average")){
			double avgDistance = 0;
			for(String sample1: cluster1.getMembers()){
				for(String sample2: cluster2.getMembers()){
					double[] vals1=getData(data, sample1, columns);
					double[] vals2=getData(data, sample2, columns);
					//double[] vals1=data.getColumn(sample1);
					//double[] vals2=data.getColumn(sample2);
					//double d=1-Statistics.pearsonDistance(vals1, vals2);
					double d=distanceFunction.measure(vals1, vals2);
					avgDistance += d;
				}
				return avgDistance/(double)(cluster1.getMembers().size() * cluster2.getMembers().size());
			}
			
		}
		else{
			throw new IllegalArgumentException("Linkage was " + metric + " it can only be one of complete, single or average");
		}
		return 0;
	}

	private double[] getData(MatrixWithHeaders data, String sample,boolean columns) {
		if(columns){return data.getColumn(sample);}
		else{return data.getRow(sample);}
	}

	private Cluster merge(Cluster cluster1, Cluster cluster2) {
		Collection<String> list=new TreeSet<String>();
		list.addAll(cluster1.getMembers());
		list.addAll(cluster2.getMembers());
		Cluster rtrn=new Cluster(list, cluster1.groupName+"_"+cluster2.groupName);
		rtrn.setSubcluster(cluster1, cluster2);
				
		return rtrn;
	}

	private ArrayList<Cluster> makeInitialClusters(MatrixWithHeaders data,	Map<String, Collection<String>> groups, boolean column) {
		ArrayList<Cluster> rtrn=new ArrayList<Cluster>();
		
		if(column){
			if(groups==null || groups.isEmpty()){
				//each sample
				for(String col: data.getColumnNames()){
					Cluster clust=new Cluster(col);
					clust.setScore(0);
					rtrn.add(clust);
				}
			}
			else{
				//each group is 1 cluster
				for(String group:groups.keySet()){
					Collection<String> vals=groups.get(group);
					Cluster clust=new Cluster(vals, group);
					clust.setScore(0);
					rtrn.add(clust);
				}
			}
		}
		else{
			
			if(groups==null || groups.isEmpty()){
				//each sample
				for(String row: data.getRowNames()){
					Cluster clust=new Cluster(row);
					clust.setScore(0);
					rtrn.add(clust);
				}
			}
			else{
				//each group is 1 cluster
				for(String group:groups.keySet()){
					Collection<String> vals=groups.get(group);
					Cluster clust=new Cluster(vals, group);
					clust.setScore(0);
					rtrn.add(clust);
				}
			}
		}
		
		return rtrn;
	}

	
	private MatrixWithHeaders order(Cluster topCluster, MatrixWithHeaders data, boolean columns){
		/*In hierarchical cluster displays, a decision is needed at each merge to specify which subtree
	    should go on the left and which on the right. Since, for n observations there are n-1 merges,
	    there are 2^{(n-1)} possible orderings for the leaves in a cluster tree, or dendrogram.
	    The algorithm used in hclust is to order the subtree so that the tighter cluster is on the left
	    (the last, i.e., most recent, merge of the left subtree is at a lower value than the last merge
	    of the right subtree). Single observations are the tightest clusters possible, and merges
	    involving two observations place them in order by their observation sequence number.*/
		
		
		
		List<String> order=topCluster.getOrdered();
		
		if(columns){
			MatrixWithHeaders rtrn=new MatrixWithHeaders(data.getRowNames(), order);
			
			for(String name: order){
				rtrn.setColumn(data.getColumn(name), name);
			}
			
			rtrn.setPIDToName(data.getPIDToName());
			return rtrn;
		}
		else{
			MatrixWithHeaders rtrn=new MatrixWithHeaders(order, data.getColumnNames());
			for(String name: order){
				rtrn.setRow(name, data.getRow(name));
			}
					
			rtrn.setPIDToName(data.getPIDToName());
			return rtrn;
			
		}
		
	}
	
	public static void main(String[] args)throws IOException, ParseException{
		if(args.length>4){
			MatrixWithHeaders data=new MatrixWithHeaders(args[0]);
			Map<String, Collection<String>> groups=AgilentUtils.parseExperimentInfoFileToGroups(new File(args[1]));
			String save=args[2];
			boolean clusterColumn=new Boolean(args[3]);
			boolean clusterRow=new Boolean(args[4]);
			HierarchicalClustering cluster=new HierarchicalClustering(data, groups, null);
			cluster.cluster(clusterRow, clusterColumn);
			cluster.write(save);
		}
		else{System.err.println(usage);}
	}
	static String usage=" args[0]=data \n args[1]=groups \n args[2]=save \n args[3]=cluster columns \n args[4]=cluster rows";

	public MatrixWithHeaders getMatrix() {
		return data;
	}
	
}
