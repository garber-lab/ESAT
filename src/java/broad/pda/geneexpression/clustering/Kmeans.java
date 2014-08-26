package broad.pda.geneexpression.clustering;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

import broad.core.datastructures.MatrixWithHeaders;
import broad.core.math.Statistics;

public class Kmeans {

	//String metric="complete";
	MatrixWithHeaders data;
	MatrixWithHeaders initialCentroids;
	KmeansInfo currentClusters;
	KmeansInfo newClusters;
	int K;
	String distMetric; //eclidean,pearson,JS
	SilhouetteInfo silhouette;
	
	
  public Kmeans(MatrixWithHeaders mat, MatrixWithHeaders initCentroids, int numOfClusters,String distanceMetric){
	 this.data=mat;
	 this.initialCentroids=initCentroids;
	 this.K=numOfClusters;
	 this.distMetric=distanceMetric;
	 this.silhouette= new SilhouetteInfo();
	 
	 if(distMetric.equalsIgnoreCase("JS"))
		 this.data=JSnormalize(mat);
	 	 
	 if (initCentroids==null)
		 this.initialCentroids=generateRandomCenters(mat,K);
	 //If the number of predefined centroids is smaller than k add random centroids
	 else if (initCentroids.rowDimension() < this.K)
		 addRandomCenters(initCentroids,mat);
	
	 MatrixWithHeaders distMat= new MatrixWithHeaders(mat.getRowNames(),this.initialCentroids.getRowNames());
	 this.currentClusters=new KmeansInfo(this.initialCentroids,distMat,this.K);
	 this.newClusters=new KmeansInfo(this.initialCentroids,distMat,this.K);

	 this.runKmeans();
  }

  
  

private MatrixWithHeaders JSnormalize(MatrixWithHeaders mat) {
		
	  MatrixWithHeaders res=new MatrixWithHeaders(mat,mat.getRowNames(),mat.getColumnNames());
	  for (int i=0; i<mat.rowDimension();i++)
		  res.setRow(mat.getRowName(i), normalizeToRelativeAbundance(mat.getRow(i)));
	  		
		return res;
	}

//assign random cluster to each data row, calc the random centroids
  private MatrixWithHeaders generateRandomCenters(MatrixWithHeaders mat, int k) {
	
	List<String> clusterNames=new LinkedList<String>();
	for (int i=1; i<=k ; i++)
		clusterNames.add(String.valueOf(i));
	MatrixWithHeaders distMat= new MatrixWithHeaders(mat.getRowNames(),clusterNames);
	MatrixWithHeaders centroids=new MatrixWithHeaders(clusterNames,mat.getColumnNames());	 
	KmeansInfo clusters =new KmeansInfo (centroids,distMat,k);
	clusters.assignRandomClusters(mat,k);
	MatrixWithHeaders res=calcCentroids(clusters);
	return res;
  }

  
  //This implementation runs one gene per iteration
  
public void runKmeans(){
  
	        // Asserts
	        if  (this.data.rowDimension() < this.K || this.K==0 ){
	        	System.err.println("USAGE: num of cluster shoud be postive and number of samples > num of clusters");
	        	return;
	        }
	        int count=1;
	        
	        //Initialize 
	        calcDistanceFromCentroids(this.currentClusters);
	        AssignSamplesToClosestCentroid(this.currentClusters);
	        //M
	        System.err.println("After initilaize Line 91");
	        this.currentClusters.printClusterSizes();
	        System.err.println();
	        
	        //Iterate
	        while(true){
	        	//kmeansIteration();
	        	kmeansByGeneIteration();
	        	System.err.println( count + "  kmeans iterations");
	        	 count+=1;
	        	//if (isStableClusters(this.currentClusters,this.newClusters,0.05))
	        		//break;
	        	if (this.isIdenticalClusters(this.currentClusters,this.newClusters))
	        		break;
	        	this.currentClusters = this.newClusters;
               
	        }
	        System.err.println("Run " + count + "  kmeans iterations");
	        this.currentClusters = this.newClusters;
  }
 
  
  private boolean isIdenticalClusters(KmeansInfo c1,KmeansInfo c2) {
	boolean res=true;
	for (String name:c1.sampleClusterAssignment.keySet()){
		if (! c1.sampleClusterAssignment.get(name).equalsIgnoreCase(c2.sampleClusterAssignment.get(name))){
			res=false; 
			break;
		}
	}
	return res;
 }

  //Check if all clusters didn't change in more than x precent
  private boolean isStableClusters(KmeansInfo c1,KmeansInfo c2,double pct) {
		boolean res=true;
		for (String name:c1.clusterSampleAssignment.keySet()){
			List<String> l1=c1.clusterSampleAssignment.get(name);
			if (! c2.clusterSampleAssignment.containsKey(name))
				continue;
			List<String> l2=c2.clusterSampleAssignment.get(name);
			double inter=0.0;
			for (String s:l1){
				if (l2.contains(s))
					inter++;
			}
			//if new cluster is larger than old in more than x% , or previous 
			//cluster got small in x% , cluster is not stable
			double len2=l2.size();
			double len1=l1.size();
			
			if (((len1-inter)/len1) > pct || ((len2-len1)/len1) > pct )
				res=false; break;
		}
		return res;
  }
  
 private void kmeansByGeneIteration(){
	  //Calculate centroids
	  MatrixWithHeaders new_centroids=calcCentroids(this.currentClusters);
	  MatrixWithHeaders distMat=new MatrixWithHeaders(this.data.getRowNames(),new_centroids.getRowNames());
	  this.newClusters=new KmeansInfo(new_centroids,distMat,this.K);
	  //Initialize cluster assigment information
	  this.newClusters.UpdateSampleClusterAssignment(this.currentClusters.sampleClusterAssignment);
	  
	  
	  //Iterate over genes- assign gene to the closest centroid,
	  //recalc the centroids of the clusters that changed
	  for (int g=0; g<this.data.rowDimension(); g++){
		 double[] dist= calcGeneDistanceFromCentroids(this.newClusters,g);
		 assignGeneToNewCluster(this.newClusters,g,dist);
	  }
	  
      
  }
  






private void kmeansIteration(){

	  //Calculate centroids
	  MatrixWithHeaders new_centroids=calcCentroids(this.currentClusters);
	  MatrixWithHeaders distMat=new MatrixWithHeaders(this.data.getRowNames(),new_centroids.getRowNames());
	  this.newClusters=new KmeansInfo(new_centroids,distMat,this.K);
      // Calc distance of all points from k centroids
      calcDistanceFromCentroids(this.newClusters);
      //Assign each sample to the closest centroid
      AssignSamplesToClosestCentroid(this.newClusters);
  }
  
  private void AssignSamplesToClosestCentroid(KmeansInfo clusters) { 
	  int j=0;
	  HashMap<String, List<String>> cs_map=new HashMap<String, List<String>>(); 
	  for (int i=0; i<this.data.rowDimension(); i++){
		  String name=this.data.getRowName(i);
		  double [] dists= clusters.centroidDist.getRow(i);
		  double minVal=Statistics.min(dists);
		  String cluster;
		  if (Double.isNaN(minVal) )
			  cluster=String.valueOf(this.K+1);
		  else{
			  for (j=0; j<dists.length; j++)
			  {if (dists[j]==minVal) break; }
			  cluster=clusters.centroids.getRowName(j);
		  }
		  clusters.sampleClusterAssignment.put(name, cluster );
		  if(! cs_map.containsKey(cluster)){
			  cs_map.put(cluster, new LinkedList<String>());
		  }
		  cs_map.get(cluster).add(name);
	  }
	  
	  //in case there are empty clusters - pick a random sample and make it a new cluster on 1 member
	  List<String> emptyList=clusters.areThereEmptyClusters(cs_map);
	  if( ! emptyList.isEmpty()){
		  
			  for(int i=0; i<emptyList.size(); i++){
				  String cluster=emptyList.get(i);
				  if (cluster.equalsIgnoreCase(String.valueOf(K+1))) // if the outliers cluster is empty- move on
						  continue;
				  String randname;
				  while(true){
					  int randSampleIx=new Double(Math.random()* (this.data.rowDimension()-1)).intValue();
					  randname=this.data.getRowName(randSampleIx);
					  if (cs_map.get(clusters.sampleClusterAssignment.get(randname)).size() >= 2)
						  break;
				  }
				  String prevCluster=clusters.sampleClusterAssignment.get(randname);
				  cs_map.get(prevCluster).remove(randname);
				  clusters.sampleClusterAssignment.put(randname, cluster );
				  if(! cs_map.containsKey(cluster)){
					  cs_map.put(cluster, new LinkedList<String>());
				  }
				  cs_map.get(cluster).add(randname);
			  }
	  }
	   
	  clusters.clusterSampleAssignment=cs_map;	  
  }

//calculate centroids - doesn't update clusters 
  private MatrixWithHeaders calcCentroids(KmeansInfo clusters){
	  
      MatrixWithHeaders centMat = new MatrixWithHeaders(clusters.centroids.getRowNames(),clusters.centroids.getColumnNames());
      for (int i=0; i<K; i++){
    	  String clusterName=centMat.getRowName(i);
    	  double[] avPattern= calcCentroids(clusters, clusterName);
    	  centMat.setRow(clusterName, avPattern);
      }
      return centMat;
          
  }
  
  private double[] calcCentroids(KmeansInfo clusters,String clusterName){
	  
	  double[] avPattern=new double[clusters.centroids.columnDimension()];
	  if (! clusters.clusterSampleAssignment.containsKey(clusterName))
		  System.err.println(clusterName + " was not found");
	  List<String> samples = clusters.clusterSampleAssignment.get(clusterName);
	  
	  
	  MatrixWithHeaders submat= this.data.submatrixByRowNames(samples);
	  if (submat != null)
		  avPattern=submat.getMeanOverAllRows();
	  else
	  {for (int j=0; j<avPattern.length; j++) avPattern[j]=Double.NaN; }
	  return  avPattern;
  }
  
  //update clusters
  private void calcDistanceFromCentroids(KmeansInfo clusters) {	
	  for (int i=0; i<data.rowDimension(); i++){
		 for (int k=0; k<clusters.centroids.rowDimension(); k++){
			 double val= distance(this.data.getRow(i),clusters.centroids.getRow(k));
			 clusters.centroidDist.set(i, k, val);
		 }
	  }
	return ;
  }
  
  
  private double[] calcGeneDistanceFromCentroids(KmeansInfo clusters, int geneNum) {
		double[] res=new double[this.K];
		for (int k=0; k<this.K; k++)
			 res[k]= distance(this.data.getRow(geneNum),clusters.centroids.getRow(k));	
		return res;
	}
  
  //assigns gene to a new cluster and updates the centroids of the leaving and entering clusters
  private void assignGeneToNewCluster(KmeansInfo clusters, int geneNum,double[] dists) {

	  double minVal=Statistics.min(dists);
	  int j=0;
	  String geneName=this.data.getRowName(geneNum);
	  String newCluster;
	  String oldCluster=clusters.sampleClusterAssignment.get(geneName);
	  if (Double.isNaN(minVal) )
		  newCluster=String.valueOf(this.K+1);
	  else{
		  for ( j=0; j<dists.length; j++)
		  {if (dists[j]==minVal) break; }
		  newCluster=clusters.centroids.getRowName(j);
	  }
	  //If the gene should be reassigned
	  if (! newCluster.equals(oldCluster)){
		  clusters.clusterSampleAssignment.get(oldCluster).remove(geneName);
		  clusters.clusterSampleAssignment.get(newCluster).add(geneName);
		  clusters.sampleClusterAssignment.put(geneName,newCluster);
		  //Recalc centroids
		  clusters.centroids.setRow(oldCluster,calcCentroids(clusters,oldCluster));
		  clusters.centroids.setRow(newCluster,calcCentroids(clusters,newCluster));
	  }
	 
		
	}


  private double distance(double[] a, double[] b) {

		return distance(a,b,this.distMetric);
}

  public static double distance(double[] a, double[] b, String metric) {

		double res= Double.NaN;
		if (metric.equalsIgnoreCase("euclidean"))
			res=Statistics.euclideanDistance(a,b);
		if (metric.equalsIgnoreCase("JS"))
			res=Statistics.JSDist(a,b);
		if (metric.equalsIgnoreCase("pearson"))
			res=1-Statistics.pearsonDistance(a, b);;
		
		return res;
  }
  private static double[] normalizeToRelativeAbundance(double[] vals) {
		int len=vals.length;
		double[] res= new double[len];
		double sum= 0;
		for (int i=0; i<len;i++) {sum+=vals[i];}
		for (int i=0; i<len;i++) {res[i]=(vals[i]/sum);}
		return res;
	}
  
  
  
  

  
  private void addRandomCenters(MatrixWithHeaders initCentroids, MatrixWithHeaders mat) {
	  
	    MatrixWithHeaders tmp= generateRandomCenters(mat,this.K);
		List<String> rowNames=tmp.getRowNames();
		List<String> fewNames=new LinkedList<String>();
		int d=this.K-initCentroids.rowDimension();
		for (int i=initCentroids.rowDimension(); i<this.K; i++)
			fewNames.add(rowNames.get(i));
		MatrixWithHeaders tmp2=tmp.submatrixByRowNames(fewNames);
				
		initCentroids.appendRows(tmp2);
	
}

  //The silhouette score of each gene is a measure of how close it is to genes within its cluster
  //in comparison to genes in the cluster that is second far from it.
  public void calcSilhouette() {
	  this.silhouette.initVals(this.currentClusters,this.data,this.distMetric,false);
	}

  public void calcLightSilhouette() {
	  this.silhouette.initVals(this.currentClusters,this.data,this.distMetric,true);
	} 
  public double getSilhouette() {
	  if (this.silhouette.isempty())
		  	return Double.NaN; 
		else
			return this.silhouette.getAvSi();
	}
	
  /*private MatrixWithHeaders addNaNcentroid(MatrixWithHeaders centroids) {
	  
	  List<String> RowList=centroids.getRowNames();
	  //List<String> DescList=centroids.getRowDescriptions();
	  RowList.add(String.valueOf(RowList.size()+1));
	  //DescList.add("NaNcluster");
	  String nanCluster=String.valueOf(RowList.size());
	  
	  MatrixWithHeaders res= new MatrixWithHeaders(RowList,centroids.getColumnNames());
	  
	  for (int i=0; i<centroids.rowDimension();i++)
		  res.setRow(centroids.getRowName(i), centroids.getRow(i));
	  double []vals =new double[centroids.columnDimension()];
	  for (int i=0;i<vals.length;i++) vals[i]=Double.NaN;
	  res.setRow(nanCluster, vals);
	  
	  this.K++;
	  return res;
}

*
*
*
*
*
* private MatrixWithHeaders calcCentroids(kmeansInfo clusters){
	  
      MatrixWithHeaders centMat = new MatrixWithHeaders(clusters.centroids.getRowNames(),clusters.centroids.getColumnNames());
      for (int i=0; i<K; i++){
    	  String clusterName=centMat.getRowName(i);
    	  double[] avPattern=new double[centMat.columnDimension()];
    	  List<String> samples = clusters.clusterSampleAssignment.get(clusterName);
    	  MatrixWithHeaders submat= this.data.submatrixByRowNames(samples);
    	  if (submat != null)
    		  avPattern=submat.getMeanOverAllRows();
    	  else
    	  {for (int j=0; j<avPattern.length; j++) avPattern[j]=Double.NaN; }
    	  centMat.setRow(clusterName, avPattern);
      }
      return centMat;
          
  }
*/
  
  
  
  public void writeClusters(BufferedWriter bw) throws IOException {
	for (int i=0; i<this.data.rowDimension();i++){
		String name=this.data.getRowName(i);
		String cluster=this.currentClusters.sampleClusterAssignment.get(name);
		bw.write(name+"\t"+cluster+"\n");
	}
	
  }

  public void writeClustersCentroids(BufferedWriter bw) throws IOException {
	this.currentClusters.centroids.writeGCT(bw);
	
  }

public void writeDistanceFromCentroids(BufferedWriter bw) throws IOException {
	calcDistanceFromCentroids(this.currentClusters);
	this.currentClusters.centroidDist.writeGCT(bw);
	
}


public void writeDistanceFromInputCentroids(BufferedWriter bw) throws IOException {

 KmeansInfo tmpClusters =new KmeansInfo(this.initialCentroids,this.currentClusters.centroidDist,this.K);
 tmpClusters.UpdateSampleClusterAssignment(this.currentClusters.sampleClusterAssignment);
 calcDistanceFromCentroids(tmpClusters);
 tmpClusters.centroidDist.writeGCT(bw);
	
}


public KmeansInfo getClusterInfo() { return currentClusters;}

public void writeDistanceFromInputCentroids(MatrixWithHeaders centroids,
		BufferedWriter bw) throws IOException {
	
	MatrixWithHeaders distMat= new MatrixWithHeaders(this.data.getRowNames(),centroids.getRowNames());
	KmeansInfo tmpClusters=new KmeansInfo(centroids,distMat,centroids.rowDimension());
    calcDistanceFromCentroids(tmpClusters);
	tmpClusters.centroidDist.writeGCT(bw);
}

public void writeSilhouetteRes(int currK, BufferedWriter bwGeneS,
		BufferedWriter bwAllS, BufferedWriter bwClusterS) throws IOException {

	if (this.silhouette.isempty())
		return;
	//For every gene write its si 
	for (int i=0; i<this.data.rowDimension();i++){
		String name=this.data.getRowName(i);
		double si=this.silhouette.getGeneSi(name);
		if(i<this.data.rowDimension()-1)
			bwGeneS.write(si+"\t");
		else
			bwGeneS.write(si+"\n");	
	}
	
	//Write total av K
	bwAllS.write(currK+"\t"+this.silhouette.getAvSi()+"\n");
	
	//write cluster K 
	bwClusterS.write(currK+"\t"+this.silhouette.getAvOfAvSi());
	HashMap<String,Double> siMap=this.silhouette.getClustersSi();
	for (Double d:siMap.values() )
		bwClusterS.write("\t"+d);
	bwClusterS.write("\n");
}


// Sub classes : 

	//*********
public static class KmeansInfo{


	MatrixWithHeaders centroids;
	MatrixWithHeaders centroidDist;
	HashMap <String,String> sampleClusterAssignment;
	HashMap <String,List<String>> clusterSampleAssignment;


	public KmeansInfo(MatrixWithHeaders initCentroids,MatrixWithHeaders initCentroidDist,int K) {
		this.centroids=initCentroids;
		this.centroidDist=initCentroidDist;
		sampleClusterAssignment=new HashMap <String,String>();
		clusterSampleAssignment=new HashMap <String,List<String>>();
		for (String cluster: centroids.getRowNames())
			clusterSampleAssignment.put(cluster, new LinkedList<String>());

		clusterSampleAssignment.put(String.valueOf(K+1), new LinkedList<String>());
	}


	public void printClusterSizes() {
		for (String cluster:this.clusterSampleAssignment.keySet())
			System.err.println(cluster +"\t"+ clusterSampleAssignment.get(cluster).size());

	}


	public boolean areThereEmptyClusters() {

		boolean res=false;
		for (String cluster:this.clusterSampleAssignment.keySet()){
			if (clusterSampleAssignment.get(cluster).isEmpty()){
				res=true;
				System.err.println("cluster " + cluster +" is empty");
			}
		}
		return res;
	}

	public List<String> areThereEmptyClusters(HashMap<String, List<String>> csMap) {

		List<String> res=new LinkedList<String>();
		for (String cluster:this.clusterSampleAssignment.keySet()){
			if ((!csMap.containsKey(cluster)) || csMap.get(cluster).isEmpty()){
				res.add(cluster);
			}
		}
		return res;
	}

	public void UpdateSampleClusterAssignment(HashMap<String, String> sampleClusterMap) {
		for (String key: sampleClusterMap.keySet()){
			String cluster=sampleClusterMap.get(key) ;
			this.sampleClusterAssignment.put(key,cluster);
			this.clusterSampleAssignment.get(cluster).add(key);

		}
	}

	public HashMap  <String,List<String>> getClusterSampleAssignment(){
		return this.clusterSampleAssignment;
	}

	public HashMap <String,String> getSampleClusterAssignment(){
		return this.sampleClusterAssignment;
	}

	public double[] getCentroid(String clustername) {
		if (this.centroids.containsRow(clustername))
			return this.centroids.getRow(clustername);
		return null;
	}


	public void assignRandomClusters(MatrixWithHeaders data,int K) {
		int j=0;
		HashMap<String, List<String>> cs_map=new HashMap<String, List<String>>(); 
		int clusterNum=0;
		for (int i=0; i<data.rowDimension(); i++){
			if (clusterNum>= K)
				clusterNum=0;
			String name=data.getRowName(i);
			if(containesNaN(data.getRow(i)) || allZero(data.getRow(i))){
				this.sampleClusterAssignment.put(name, String.valueOf(K+1) );
				continue;
			}
			String clusterName= this.centroids.getRowName(clusterNum);
			this.sampleClusterAssignment.put(name, clusterName );
			if(! cs_map.containsKey(clusterName)){
				cs_map.put(clusterName, new LinkedList<String>());
			}
			cs_map.get(clusterName).add(name);
			clusterNum++;
		}
		this.clusterSampleAssignment=cs_map;

	}


	private boolean allZero(double[] row) {
		int numzero=0;
		for (int i=0; i<row.length;i++){
			if ((row[i])==0)
				numzero++;
		}
		return (numzero==row.length);
	}


	private boolean containesNaN(double[] row) {
		for (int i=0; i<row.length;i++){
			if (Double.isNaN(row[i]))
				return true;
		}
		return false;
	}







}


//*********
class SilhouetteInfo{

	MatrixWithHeaders pairwiseDist;
	MatrixWithHeaders silhouetteScr;
	double avSilhouette=0; //mean of all genes (account for cluster size)
	double avAcrossClustersSilhouette=0; //mean of cluster means (don't account for cluster size)
	HashMap<String,Double> clusterAvSilhouette;

	public SilhouetteInfo (){
		this.pairwiseDist=null;
		this.silhouetteScr=null;
		this.clusterAvSilhouette=new HashMap<String,Double> ();
	}

	public double getAvOfAvSi() {return this.avAcrossClustersSilhouette;}


	HashMap<String,Double> getClustersSi() { return clusterAvSilhouette;}

	public double getAvSi() {return this.avSilhouette;}

	public double getGeneSi(String name) {return this.silhouetteScr.get(name,"si");	}

	public boolean isempty() {
		if ( this.silhouetteScr==null) 
			return true;
		return false;
	}

	//Calc the pairwise distance between all genes, and calc silhouette val for each gene
	public void initVals(KmeansInfo currentClusters,MatrixWithHeaders data, String distMetric,boolean lightVersion) {

		LinkedList<String> sampleNames=new LinkedList<String>();
		sampleNames.addAll(currentClusters.getSampleClusterAssignment().keySet());
		LinkedList<String> sScrs=new LinkedList<String>();
		sScrs.add("ai");  sScrs.add("bi");  sScrs.add("si"); 


		this.silhouetteScr=new MatrixWithHeaders(sampleNames,sScrs);	

		if(lightVersion){
			calcLightSilhouette(data,sampleNames,currentClusters,distMetric);
		}
		else{
			this.pairwiseDist=new MatrixWithHeaders(sampleNames,sampleNames);
			//calc pairwise dist
			calcPairwiseDist(sampleNames,data,  distMetric);
			//calc silhouette
			calcSilhouette(sampleNames,currentClusters);
		}
	}


	private void calcPairwiseDist(LinkedList<String> sampleNames,
			MatrixWithHeaders data, String distMetric) {

		for(int i=0; i<sampleNames.size(); i++){
			String name=sampleNames.get(i);
			double[]a=data.getRow(name);
			for(int j=0; j<sampleNames.size(); j++){
				String pair=sampleNames.get(j);
				//if (new Double(this.pairwiseDist.get(name,pair)).isNaN()){
				double[] b =data.getRow(pair);
				double res=Kmeans.distance(a,b,distMetric);
				this.pairwiseDist.set(name,pair,res);
				this.pairwiseDist.set(pair,name,res);
				//}


			}
		}

	}

	private void calcSilhouette(LinkedList<String> sampleNames, KmeansInfo currentClusters) {

		HashMap<String, List<String>> clusterSampleMap=currentClusters.getClusterSampleAssignment();
		for (int i=0; i<sampleNames.size(); i++){
			String name=sampleNames.get(i);
			String currCluster=currentClusters.getSampleClusterAssignment().get(name);
			ArrayList<Double> biArr=new ArrayList<Double>();
			double ai=calcAvDistanceFromCluster(name,clusterSampleMap.get(currCluster));
			for(String cluster: clusterSampleMap.keySet()){
				if (! cluster.equals(currCluster)){
					double bi= calcAvDistanceFromCluster(name,clusterSampleMap.get(cluster));
					if (! new Double(bi).isNaN())
						biArr.add(bi);
				}
			}

			double bi=Statistics.min(biArr);
			double si= calcSi(ai,bi);
			double[] arr={ai,bi,si};
			this.silhouetteScr.setRow(name,arr );
		}
		//update av silhouette
		this.avSilhouette=Statistics.mean(this.silhouetteScr.getColumn("si"));

		//update cluster av sillhouete
		for(String cluster: clusterSampleMap.keySet()){
			ArrayList<Double> siArr=new ArrayList<Double>();
			List<String> clusterSamples =clusterSampleMap.get(cluster);
			for (String sample:clusterSamples)
				siArr.add(this.silhouetteScr.get(sample,"si"));
			this.clusterAvSilhouette.put(cluster,Statistics.mean(siArr));
		}

	}

	private double calcAvDistanceFromCluster(String name, List<String> pairlist) {
		ArrayList<Double> dArr=new ArrayList<Double>();
		for (String pair:pairlist)
			dArr.add(this.pairwiseDist.get(name,pair));
		double res=Statistics.mean(dArr);
		return res;
	}

	//si= (bi-ai)/(max(ai,bi))
	private double calcSi(double ai, double bi) {

		double si= 0;
		if (ai>bi)
			si=bi/ai -1;
		if(bi>ai)
			si=1- ai/bi;
		return si;
	}


	public void calcLightSilhouette(MatrixWithHeaders data,LinkedList<String> sampleNames,
			KmeansInfo currentClusters,String distMetric) {


		HashMap<String, List<String>> clusterSampleMap=currentClusters.getClusterSampleAssignment();
		for (int i=0; i<sampleNames.size(); i++){
			String name=sampleNames.get(i);
			double[] sampleExp=data.getRow(name);
			String currCluster=currentClusters.getSampleClusterAssignment().get(name);
			ArrayList<Double> biArr=new ArrayList<Double>();

			double ai=Kmeans.distance(sampleExp,currentClusters.getCentroid(currCluster),distMetric);

			for(String cluster: clusterSampleMap.keySet()){
				if (! cluster.equals(currCluster)){
					double [] c=currentClusters.getCentroid(cluster);
					if (c == null)
						continue;
					double bi= Kmeans.distance(sampleExp,c,distMetric);
					if (! new Double(bi).isNaN())
						biArr.add(bi);
				}
			}

			double bi=Statistics.min(biArr);
			double si= calcSi(ai,bi);
			double[] arr={ai,bi,si};
			this.silhouetteScr.setRow(name,arr );
		}
		//update av silhouette
		this.avSilhouette=Statistics.mean(this.silhouetteScr.getColumn("si"));

		//update cluster av sillhouete
		ArrayList<Double> avOfAv =new ArrayList<Double>();
		for(String cluster: clusterSampleMap.keySet()){
			ArrayList<Double> siArr=new ArrayList<Double>();
			List<String> clusterSamples =clusterSampleMap.get(cluster);
			for (String sample:clusterSamples)
				siArr.add(this.silhouetteScr.get(sample,"si"));
			Double meanSi=Statistics.mean(siArr);
			this.clusterAvSilhouette.put(cluster,meanSi);
			if (! meanSi.isNaN())
				avOfAv.add(meanSi);
		}
		this.avAcrossClustersSilhouette=Statistics.mean(avOfAv);
	}



}


	

	


	

  




}
	        

