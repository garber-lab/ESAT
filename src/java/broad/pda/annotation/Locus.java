package broad.pda.annotation;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import nextgen.core.annotation.Gene;

import org.broad.tribble.annotation.Strand;

import broad.core.datastructures.IntervalTree;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.NeighborAnalysis;
import broad.pda.gene.GeneWithIsoforms;

public class Locus {

	private static final String SetIsoformsMap = null;
	private GeneWithIsoforms refTranscript;
	//private IntervalTree<RefSeqGeneWithIsoforms> AllIsoforms;
	private BEDFileParser AllIsoformsBed; //This will store all the isoforms such that 2 isoforms that span the same genomic region will not override one another
	private HashMap<String,Collection<Gene>> SetIsoformMap ; 
	
	public Locus(GeneWithIsoforms g){
		AllIsoformsBed=new BEDFileParser();
		SetIsoformMap=new HashMap<String,Collection<Gene>> ();
		refTranscript=g;
	}

	public void addIsoform(GeneWithIsoforms iso,String set) {
		//AllIsoforms.put(iso.getStart(), iso.getEnd(), iso);
		Collection <Gene>allIso=iso.getAllIsoforms();
		AllIsoformsBed.addRefSeqSet(allIso);
		if (!SetIsoformMap.containsKey(set)){
			Collection<Gene> c=new LinkedList<Gene>();
			SetIsoformMap.put(set,c);
		}
		SetIsoformMap.get(set).addAll(allIso);
	}
	
	public GeneWithIsoforms getReferenceTranscript(){
		return this.refTranscript;
	}
	public HashMap<String,Integer> getNumIsoformsPerSet(){
		
		HashMap<String,Integer> res=new HashMap<String,Integer>();
		for (String s: SetIsoformMap.keySet()){
			res.put(s,this.SetIsoformMap.get(s).size() );
		}
		
		return res;
	}
	
	public HashMap<String,Double> getMaxScorePerSet(){
		HashMap<String,Double> res=new HashMap<String,Double>();
		for (String s: SetIsoformMap.keySet()){
			Collection<Gene> isoforms=SetIsoformMap.get(s);
			Double maxScr=0.0;
			boolean flag=true;
			for (Gene iso: isoforms){
				if (flag) {maxScr=iso.getBedScore(); flag=false;}
				else {maxScr=Math.max(maxScr, iso.getBedScore());}
			}
			res.put(s,maxScr);
		}
		return res;
	}
	
	public HashMap<String,Integer> getIsComaptibleWithRefPerSet(){
		HashMap<String,Integer> res=new HashMap<String,Integer>();
		for (String s: SetIsoformMap.keySet()){
			Collection<Gene> isoforms=SetIsoformMap.get(s);
			Integer maxScr=0;
			for (Gene iso: isoforms){
				if (refTranscript.numOfCompatibleIntrons(iso)== refTranscript.getNumExons()-1)
					maxScr++;
			}
			res.put(s,maxScr);
		}
		return res;	
	}
	
	public HashMap<String,Integer> getIsPartiallyComaptibleWithRefPerSet(){
		HashMap<String,Integer> res=new HashMap<String,Integer>();
		for (String s: SetIsoformMap.keySet()){
			Collection<Gene> isoforms=SetIsoformMap.get(s);
			Integer maxScr=0;
			for (Gene iso: isoforms){
				if (refTranscript.numOfCompatibleIntrons(iso)>0)
					maxScr++;
			}
			res.put(s,maxScr);
		}
		return res;	
	}
	
	public int getIsComaptibleBetweenSubsets(ArrayList<String> set1Lst,ArrayList<String> set2Lst,boolean isFully){
		
		Iterator<String> set1It=set1Lst.iterator();
		int compatible=0;
		while (set1It.hasNext()){
			String set1name=set1It.next();
			if (! SetIsoformMap.containsKey(set1name))
				{ //System.err.println(set1name);
				continue;}
			Collection<Gene> isoforms1=SetIsoformMap.get(set1name);
			for (Gene iso1:isoforms1){
					Iterator<String> set2It=set2Lst.iterator();
					while (set2It.hasNext()){
						String set2name=set2It.next();
						if (! SetIsoformMap.containsKey(set2name))
						{ //System.err.println(set2name);
						continue;}
						Collection<Gene> isoforms2=SetIsoformMap.get(set2name);
						for (Gene iso2: isoforms2){
								if (isFully && isFullyComaptible(iso1,iso2,false))
									compatible++;
								if ((isFully==false) && (iso1.numOfCompatibleIntrons(iso2)>0))
									compatible++;
							}
						}
				}
			}
	return compatible;			
	}
	
	public int getIsComaptibleBetweenAny2Sets(ArrayList<String> setLst,	boolean isFully) {
		int compatible=0;
		for (int i=0; i< setLst.size()-1;i++){
			String set1name=setLst.get(i);
			if (! SetIsoformMap.containsKey(set1name))
				continue;
			Collection<Gene> isoforms1=SetIsoformMap.get(set1name);
			for (Gene iso1:isoforms1){
					for (int j=i+1; j<setLst.size(); j++ ){
						String set2name=setLst.get(j);
						if (! SetIsoformMap.containsKey(set2name))
							continue;
						Collection<Gene> isoforms2=SetIsoformMap.get(set2name);
						for (Gene iso2: isoforms2){
								if (isFully && isFullyComaptible(iso1,iso2,false))
									compatible++;
								if ((isFully==false) && (iso1.numOfCompatibleIntrons(iso2)>0))
									compatible++;
						}
					}
			}
				
		}
		
		return compatible;		
	}
	
	
	public int getIsComaptibleBetweenSubsets(ArrayList<String> set1Lst,
			ArrayList<String> set2Lst, BEDFileParser uniqIsoBed) {
		
		Iterator<GeneWithIsoforms> unqRefsIt = uniqIsoBed.getOverlappers(this.refTranscript).valueIterator();
		int compatible=0;
		int i=0;
		ArrayList<Integer> compNum1=new ArrayList<Integer>();
		ArrayList<Integer> compNum2=new ArrayList<Integer>();
		//go over all unique reference isoforms
		while (unqRefsIt.hasNext()){
			Collection <Gene> refs = unqRefsIt.next().getAllIsoforms();
			for (Gene ref:refs){
				compNum1.add(i, numCompatibleInSet(set1Lst,ref));
				compNum2.add(i, numCompatibleInSet(set2Lst,ref));
				i++;
			}
		}
		
		for (int j=0; j<compNum1.size(); j++){
			int a=compNum1.get(j);
			int c= compNum2.get(j);
			if (a>0 & c >0){
				compatible=Math.max(compatible, a+c);
			}
		}
		return compatible;
		
	}
	
	public int getIsComaptibleBetweenAny2Sets(ArrayList<String> set3Lst,
			 BEDFileParser uniqIsoBed) {
		Iterator<GeneWithIsoforms> unqRefsIt = uniqIsoBed.getOverlappers(this.refTranscript).valueIterator();
		int compatible=0;
		int i=0;
		ArrayList<Integer> compNum1=new ArrayList<Integer>();
		while (unqRefsIt.hasNext()){
			Collection <Gene> refs = unqRefsIt.next().getAllIsoforms();
			for (Gene ref:refs){
				compNum1.add(i, numCompatibleInSet(set3Lst,ref));
				i++;
			}
		}
		for (int j=0; j<compNum1.size(); j++)
			compatible=Math.max(compatible, compNum1.get(j));
		
		if (compatible >1)
			return compatible;
		else
			return 0;
	}
	
    private int numCompatibleInSet(ArrayList<String> set1Lst, Gene ref) {
    	
    	int res=0;
    	Iterator<String> set1It=set1Lst.iterator();
		while (set1It.hasNext()){
			String set1name=set1It.next();
			if (! SetIsoformMap.containsKey(set1name))
				continue;
			Collection<Gene> isoforms1=SetIsoformMap.get(set1name);
			for (Gene iso1:isoforms1){
				if (isFullyComaptible(ref,iso1,false))
					res++;	
			}
		}
		return res;
	}

	
	
		
	public Gene getExonUnion(){
		
		boolean first=true;
		Gene mergedElement=refTranscript;
		
		if (! this.SetIsoformMap.isEmpty()){
			Iterator<GeneWithIsoforms> it=AllIsoformsBed.getChrTree(this.refTranscript.getChr()).valueIterator();
			
			Gene tmp;
			while(it.hasNext()){
				for (Gene iso: it.next().getAllIsoforms()){
					if (first) {mergedElement=iso; first=false;}
					else {mergedElement=  mergedElement.takeUnion(iso);}
				}
			}
		}
		return mergedElement;
	}
	
	//Output: First element specifies if genes are in the same orientation, and second element specifies the distance
	public double[] getDistanceTo3primeNeighbor(BEDFileParser geneSet)throws IOException{
		double[]res=new double[2];
		Gene consenzus=getExonUnion();
		BEDFileParser myBed=new BEDFileParser();
		myBed.addRefSeq(consenzus);
		NeighborAnalysis mNA=new NeighborAnalysis(myBed,geneSet);
		mNA.getNeighbors( new BEDFileParser(), null,false,-1);
		mNA.updateDistanceToNeighbors();
		int ix=1;
		if (consenzus.getOrientation().equals(Strand.NEGATIVE))
			ix=0;
		Gene leftNeighbor= mNA.getNeighbors(consenzus)[ix];
		if (leftNeighbor!= null) 
			res[0]= leftNeighbor.getOrientation().equals( consenzus.getOrientation())? 1:0;
		else
			res[0]=0;
		res[1]=mNA.getNeighborDistance(consenzus)[ix];
		
		return res;		
		
	}
	
	//Output: First element specifies if genes are in the same orientation, and second element specifies the distance
	public double[] getDistanceTo5primeNeighbor(BEDFileParser geneSet)throws IOException{
		double[]res=new double[2];
		Gene consenzus=getExonUnion();
		BEDFileParser myBed=new BEDFileParser();
		myBed.addRefSeq(consenzus);
		NeighborAnalysis mNA=new NeighborAnalysis(myBed,geneSet);
		mNA.getNeighbors( new BEDFileParser(), null,false,-1);
		mNA.updateDistanceToNeighbors();
		int ix=0;
		if (consenzus.getOrientation().equals(Strand.NEGATIVE))
			ix=1;
		Gene Neighbor= mNA.getNeighbors(consenzus)[ix];
		if (Neighbor!= null) 
			res[0]= Neighbor.getOrientation().equals( consenzus.getOrientation())? 1:0;
		else
			res[0]=0;
		res[1]=mNA.getNeighborDistance(consenzus)[ix];
		
		return res;		
		
	}
	
	public double[]distanceToOppositeStrand5primeNeighbor (BEDFileParser geneSet) throws IOException{
		
		double[]res=this.getDistanceTo5primeNeighbor(geneSet);
		if (res[0]!=0){ //Have the same orientation
			res[0]=0;
			res[1]=0;
		}
		return res;
	}

	public double getMaxScoreAcrossAllIso() {

		double maxScr=0;
		HashMap<String, Double> setMaxScoreMap=this.getMaxScorePerSet();
		Collection<Double> vals=setMaxScoreMap.values();
		for (Double d: vals){
			maxScr=Math.max(d,maxScr);
		}
		return maxScr;
	}

	//adds a suffix to the name of every isoform that specifies the loci position
	public void updateIsoformsNameWithLocusName() {
		
		String refName=";"+this.refTranscript.getChr()+":"+this.refTranscript.getStart()+"-"+this.refTranscript.getEnd();
		if (this.SetIsoformMap.isEmpty())
			return;
		
		/* This are supposed to be the same objects as in the sets themselves- with the exception of the transcript that initiates Super
		*/Iterator<GeneWithIsoforms> it =this.AllIsoformsBed.getChrTree(this.refTranscript.getChr()).valueIterator();
		while(it.hasNext()){
			GeneWithIsoforms g=it.next();
			g.addSuffixToName(refName);
			
		}
		
		
		for (String s:this.SetIsoformMap.keySet()){
			for (Gene g : this.SetIsoformMap.get(s))
				g.addSuffixToName(refName);
		}
		
	}

	public Collection<Gene> getSetIsoforms(String setName) {
		
		Collection<Gene> isoforms=new LinkedList<Gene> ();
		if ( SetIsoformMap.containsKey(setName))
			 isoforms=SetIsoformMap.get(setName);
		else
			System.err.println(setName);
	
		return isoforms;
	}

	public Collection<Gene> getAllIsoforms(String excludeSet) {
		
		LinkedList<Gene> lst= new LinkedList<Gene> ();
		if (! this.SetIsoformMap.isEmpty()){
			for (String s:this.SetIsoformMap.keySet()){
				if (!s.equalsIgnoreCase(excludeSet))
					lst.addAll(this.SetIsoformMap.get(s));
			}
		}
		return lst;
	}

	public Collection<Gene> getComaptibleBetweenAny2Sets(ArrayList<String> setLst, boolean isFully) {
		
		HashSet<Gene> res=new HashSet<Gene>();
		HashMap<String,Gene> map= new HashMap<String,Gene>();
		for (int i=0; i< setLst.size()-1;i++){
			String set1name=setLst.get(i);
			//System.err.println(set1name);
			if (! SetIsoformMap.containsKey(set1name))
				continue;
			Collection<Gene> isoforms1=SetIsoformMap.get(set1name);
			for (Gene iso1:isoforms1){
					for (int j=i+1; j<setLst.size(); j++ ){
						String set2name=setLst.get(j);
						if (! SetIsoformMap.containsKey(set2name))
							continue;
						Collection<Gene> isoforms2=SetIsoformMap.get(set2name);
						for (Gene iso2: isoforms2){
								if (isFully && isFullyComaptible(iso1,iso2,false)){
									map.put(iso1.getName(),iso1);
									//System.err.println("Found compatible " + iso1.getName()+ "From " +set1name +" and "+ set2name + "\n" + iso1.toBED());
								}
								if ((isFully==false) && (iso1.numOfCompatibleIntrons(iso2)>0)){
									map.put(iso1.getName(),iso1);
								}
						}
						
					}
			}
		}
		for(String s:map.keySet())
			res.add(map.get(s));
		return res;
	}

	public Collection<Gene> getUnqComaptibleBetweenAny2Sets(
			ArrayList<String> setLst, Collection<? extends Gene> unqIsoSet, boolean isFully) {
		
		LinkedList<Gene> candidates =new LinkedList <Gene>();
		for (Gene g:unqIsoSet){
			int scr=this.getNumberOfCompatibleSets(g,setLst,isFully);
			if (  scr>=2){
				g.setBedScore(scr);
				candidates.add(g);
			}
		}
				
		return candidates;
	}

	public Collection<Gene> getUnqComaptibleBetweenSubsets( ArrayList<String> set1Lst,ArrayList<String> set2Lst,
			Collection<Gene> unqIsoSet, boolean isFully) {
		
		LinkedList<Gene> candidates =new LinkedList <Gene>();
		for (Gene g:unqIsoSet){
			int scr1=this.getNumberOfCompatibleSets(g,set1Lst,isFully);
			int scr2=this.getNumberOfCompatibleSets(g,set2Lst,isFully);
			if ( scr1>0 && scr2>0){
				g.setBedScore(scr1+scr2);
				candidates.add(g);
			}
		}		
		return candidates;
	}
	
	public Collection<Gene> getComaptibleIsoforms(Gene g,boolean isFully) {
		
		Collection<Gene> lst=new LinkedList<Gene>();
		for (String setName:SetIsoformMap.keySet()){
			Collection<Gene> isoforms=SetIsoformMap.get(setName);
			for (Gene iso: isoforms){
				if (isFully && (g.numOfCompatibleIntrons(iso)==g.getNumExons()-1))
					lst.add(iso);
				if ((isFully==false) && (g.numOfCompatibleIntrons(iso)>0))
					lst.add(iso);
			}
		}
		return lst;
	}
	

	private int getNumberOfCompatibleSets(Gene g,ArrayList<String> setLst, boolean isFully) {
		
		int totalSet=0;
		for (String setName:setLst){
			if (! SetIsoformMap.containsKey(setName))
				continue;
			Collection<Gene> isoforms=SetIsoformMap.get(setName);
			boolean inSet=false;
			for (Gene iso: isoforms){
				if (isFully && (g.numOfCompatibleIntrons(iso)==g.getNumExons()-1))
					inSet=true;
				if ((isFully==false) && (g.numOfCompatibleIntrons(iso)>0))
					inSet=true;
			}
			if (inSet)
				totalSet++;
		}
		return totalSet;
		
		
	}
		
	
    public static Collection<Gene> SelectLongestIntronChainCandidate(Collection<Gene> betweenSubsetsIsoSet) {
			
    		int length=0;
    		Collection<Gene> lst= new LinkedList<Gene>();
			for (Gene g: betweenSubsetsIsoSet )
				length=Math.max(length,g.getNumExons());
			for (Gene g: betweenSubsetsIsoSet ){
				if (g.getNumExons()==length)
					lst.add(g);
			}
			
			return lst;
	}

	public static boolean isFullyComaptible(Gene iso1,Gene iso2,boolean stringent){
		int num1=iso1.numOfCompatibleIntrons(iso2);
		int num2=iso2.numOfCompatibleIntrons(iso1);
		boolean or =(num1 ==(iso1.getNumExons()-1) || num2 ==(iso2.getNumExons()-1) );
		boolean and=(num1 ==(iso1.getNumExons()-1) && num2 ==(iso2.getNumExons()-1) );
		if (stringent)
			return and;
		return or ; 
	}

	//RefSeq is a name of one set  stored in the locus that serves as 
	// the set of unique reference isoforms (like a cuff compare result)
	public Gene getLongestIntornChainIso(String refSet) {
		Gene res=null;
		/*if (! this.SetIsoformMap.containsKey(refSet)){
			System.err.println("No key , Ref set is: "+refSet + "keys are");
			Set<String> tmp=this.SetIsoformMap.keySet();
			for (String s:tmp)
				System.err.println(s);
		}*/
		if (refSet != null && this.SetIsoformMap.containsKey(refSet))
			res= longestIntornChainIso(this.SetIsoformMap.get(refSet));
		else{
			Collection <Gene> c=new LinkedList<Gene> ();
			for (String s:this.SetIsoformMap.keySet()){
				c.add(longestIntornChainIso(this.SetIsoformMap.get(s)));
			}
			res=longestIntornChainIso(c);
		}
		return res;
	}

	private Gene longestIntornChainIso(Collection<Gene> collection) {
		Gene res=null;
		HashMap<Gene,Double> map= new HashMap<Gene,Double> ();
		for (Gene g: collection)
			map.put(g,new Double(g.getNumExons()));
		Collection <Gene> longest =selectMaxScrIso (map,-1);
		HashMap<Gene,Double> map2= new HashMap<Gene,Double> ();
		for (Gene g: longest) //from the longest intron chain select the one that has the largest size
			map2.put(g,new Double(g.getTranscriptLength()));
		LinkedList <Gene> longest2 =(LinkedList<Gene>) selectMaxScrIso (map2,-1);
		if (longest2.size()>0)
			res=longest2.get(0);
		return res;
	}
		
	private static Collection <Gene> selectMaxScrIso(HashMap<Gene,Double> map,double initScr) {
		double scr=initScr;
		Collection <Gene>res=new LinkedList<Gene>(); 
		for (Gene g: map.keySet() )
			scr=Math.max(scr,map.get(g));
		for (Gene g: map.keySet() ){
			if (map.get(g)==scr)
				res.add(g);
		}
		return res;
	}

	public Gene getMerged() {
		
		List<Gene> genes= getAllMerged();
		return genes.get(0);
	}

	public List<Gene> getAllMerged() {
		
		BEDFileParser bed=new BEDFileParser();
		bed.addRefSeqSet(this.AllIsoformsBed.GetGenes());
		bed.merge();
		List<Gene> genes=bed.GetGenes();
		if (genes.size()>1){
			for (int i=0; i<genes.size();i++)
				System.err.println("merged in fragments: " + genes.get(i).toBED());
		}
		return genes;
	}


	public int numOfOverlapSets(ArrayList<String> setLst) {
		int res=0;
		for (int i=0; i<setLst.size() ;i++){
			boolean b=this.SetIsoformMap.containsKey(setLst.get(i));
			if (this.SetIsoformMap.containsKey(setLst.get(i)) &&  this.SetIsoformMap.get(setLst.get(i)).size()>0)
				res++;
		}
		return res;
	}

	public double getMaxPctCovered() {
		return getMaxCovered(0);
	}

	public double getMaxExonsCovered() {
		return getMaxCovered(1);
	}

	public double getMaxPctGenomeCovered() {
		return getMaxCovered(2);
	}
	
	private double getMaxCovered(int i){
		HashMap<Gene, ArrayList<Double>> overlapMap= pctOverlapWithRef();
		double maxOver=0;
		for (Gene g:overlapMap.keySet() )
			maxOver=Math.max(maxOver,overlapMap.get(g).get(i));
		return maxOver;
	}
	
	
	public HashMap<Gene, ArrayList<Double>> pctOverlapWithRef()
	{
		HashMap<Gene, ArrayList<Double>> res=new HashMap<Gene, ArrayList<Double>>();
	
		for (String s: SetIsoformMap.keySet()){
			Collection<Gene> isoforms=SetIsoformMap.get(s);
			ArrayList<Double> arr= new ArrayList<Double>();
			for (Gene iso: isoforms){
				arr.add(0,this.refTranscript.percentOverlapping(iso));
				arr.add(1,new Double(this.refTranscript.getMerged().numOverlappingExons(iso)));
				arr.add(2,this.refTranscript.getMerged().percentGenomeOverlapping(iso));
				res.put(iso,arr);
			}
			
		}
		return res;	
	
	}

	//TODO- decide how to implement this
	//Bias will be determined if only the first or last 50% of exons are covered 
	public int isTerminiCoverageBias() {
		int bias=0;
		Gene ref= this.getReferenceTranscript();
		int numEx=ref.getNumExons();
		Alignments[] aln = ref.getExons();
		for (String s: SetIsoformMap.keySet()){
			Collection<Gene> isoforms=SetIsoformMap.get(s);
			for (Gene iso: isoforms){
				int first=0;
				int last=0;
				for (int i=0;i<numEx; i++){
					
				}
				
			}
		}

		return bias;
	}

	public Collection<? extends Gene> getAllIsoforms() {
		//Don't exclude any set
		return getAllIsoforms("");
	}

	public Integer getNumberOfIsoforms() {
		
		return this.AllIsoformsBed.getNumberOfIsoforms();
	}

	public Collection<? extends Gene> selectRandIsoSubset(
			Integer maxIsoPerLoci) {
		
		
		Collection<Gene> rtrn = new TreeSet <Gene>();
		int i=0;
		for (Gene g: this.getAllIsoforms()){
			if (i< maxIsoPerLoci){
				rtrn.add(g); i++;
			}
			else
				break;
		}
			
		return rtrn;
	
	}

	
	
}
