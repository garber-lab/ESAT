package broad.pda.geneexpression.clustering;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.TreeSet;

import broad.core.datastructures.Pair;

//A collection of samples or genes that represent a cluster
public class Cluster {

	Pair<Cluster> subclusters;
	
	Collection<String> members;
	String groupName;
	
	double score;
		
	public Cluster(String col) {
		this.members=new ArrayList<String>();
		this.members.add(col);
		this.subclusters=new Pair<Cluster>();
	}
	
	public Cluster(Collection<String> mem, String group) {
		this.members=mem;
		this.groupName=group;
		this.subclusters=new Pair<Cluster>();
	}
	
	public void setScore(double score){this.score=score;}
	public double getScore(){return this.score;}
	
	public boolean hasSubClusters(){
		if(this.subclusters==null){return false;}
		if(this.subclusters.isEmpty()){return false;}
		return true;
	}
	
	public void addToCluster(String member){this.members.add(member);}

	//Must be sorted by cluster score
	public Pair<Cluster> getSubclusters() {
		return this.subclusters;
	}


	public Collection<String> getMembers() {
		return this.members;
	}

	public void setSubcluster(Cluster cluster1, Cluster cluster2) {
		Pair<Cluster> pair=new Pair<Cluster>(cluster1, cluster2);
		this.subclusters=pair;
	}

	
	//Subcluster with smaller distance
	public Cluster getLeft() {
		double score1=this.getSubclusters().getValue1().getScore();
		double score2=this.getSubclusters().getValue2().getScore();
		if(score1<score2){return this.getSubclusters().getValue1();}
		else{return this.getSubclusters().getValue2();}
	}
	
	//Subcluster with larger distance
	public Cluster getRight() {
		double score1=this.getSubclusters().getValue1().getScore();
		double score2=this.getSubclusters().getValue2().getScore();
		if(score1>=score2){return this.getSubclusters().getValue1();}
		else{return this.getSubclusters().getValue2();}
	}

	//In order tree walk
	//TODO make sure it respect the initalization order when no other info exists
	public List<String> getOrdered() {
		List<String> rtrn=new ArrayList();
		if(this.hasSubClusters()){
			//System.err.println(this.members+" Subclusters "+this.getLeft()+" "+this.getRight());
			rtrn.addAll(this.getLeft().getOrdered());
			//System.err.println("traversed.. "+this.members);
			rtrn.addAll(this.getRight().getOrdered());
		}
		else{
			//System.err.println(this.members);
			return new ArrayList(this.members);
		}
		return rtrn;
	}
	
	/*public void getOrdered() {
		if(this.hasSubClusters()){
			//System.err.println(this.members+" Subclusters "+this.getLeft()+" "+this.getRight());
			this.getLeft().getOrdered();
			//System.err.println("traversed.. "+this.members);
			this.getRight().getOrdered();
		}
		else{System.err.println("END NODE "+this.members);}
		
	}*/
	
	public String toString(){
		return this.getMembers().toString();
	}
	
}
