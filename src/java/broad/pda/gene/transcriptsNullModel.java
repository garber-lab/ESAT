package broad.pda.gene;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.TreeSet;

import nextgen.core.annotation.AbstractAnnotation;
import nextgen.core.annotation.Annotation;
import nextgen.core.annotation.BasicAnnotation;
import nextgen.core.annotation.Gene;

import broad.core.datastructures.IntervalTree;
import broad.core.datastructures.IntervalTree.Node;
import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.NeighborsNullModel.chrSizeNode;


public class transcriptsNullModel {
	
	
	//This class generates a random model of intergenic transcripts.
	//Given a genome , a file of reference annotation and annotation of centromers
	//this class will generate a reference intergenic genome represented as an interval tree of segments
	//It can then sample uniformly sets of random transcript from this 
	//intergenic genome that maintain the gene-structures in a given input file
	static int SCALE_FACTOR=1000000;
	int GenomeSize;
	HashMap <String, IntervalTree<Alignments>> intergenicGenome;
	IntervalTree <Alignments> probFractionedChrTree; // chrs are held based on size in the genome
	HashMap <String, IntervalTree<Node<Alignments>>> probFractionedIntergenicGenome;
	HashMap <String , Integer> chrSizes;
	ArrayList <BEDFileParser> randTranscriptsList;
	
	//ctor
	public transcriptsNullModel(BEDFileParser refGenes, BEDFileParser centromers,String chrSizeF ) throws IOException {
		
		intergenicGenome = new HashMap <String, IntervalTree<Alignments>> ();
		probFractionedChrTree = new IntervalTree <Alignments> (); 
		probFractionedIntergenicGenome = new HashMap <String, IntervalTree<Node<Alignments>>> ();
		chrSizes = new HashMap <String , Integer> ();
		randTranscriptsList = new ArrayList <BEDFileParser> ();
		GenomeSize=0;
		
		initIntergenicGenome(chrSizeF); //make the original intergenicGenome
		removeRegionFromIntergenicGenome(refGenes); //remove known annotations
		removeRegionFromIntergenicGenome(centromers);
		makeProbFractionedIntergenicGenome(); //use the reduced intergenic genome to generate 
											// a probability for each  genome segment that is left
											//update chrSizes
		
	}
	
	/*
	private double initChrSize(String chrSizeF) throws IOException {

		
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(chrSizeF)));
		String nextLine;
		int IsoformMissCntr=0;
		double sum=0;
		HashMap<String,Double> sizemap=  new HashMap<String,Double>();
		
		while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
			String[] tokens=nextLine.split(" ");
			if (tokens[0].contains("rand") || tokens[0].contains("chrM"))
				continue;
			sizemap.put(tokens[0], Double.valueOf(tokens[1]));
			sum+=Double.valueOf(tokens[1]);
		}
		reader.close();
		HashMap<String,Integer> probmap= new HashMap<String,Integer>();
		for (String chr : sizemap.keySet()){
			int s=new Double((sizemap.get(chr)/sum)*100000).intValue(); ;
			probmap.put(chr, s);
		}
		int st=1;
		for (String chr : probmap.keySet()){
			int val=probmap.get(chr);
			int stop= (st +val)-1;
			chrSizeNode V=new chrSizeNode(chr,st,stop);
			this.chrSizes.put (chr,V );
			//Iterator<Node<chrSizeNode>> testIt=this.chrProb.overlappers(st+1, st+5);
			//boolean b=testIt.hasNext();
			st=stop+1;
		}
		return (st-1);
	}
	*/
	
	
	
	private void initIntergenicGenome(String chrSizeF) throws NumberFormatException, IOException {

		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(chrSizeF)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
			String[] tokens=nextLine.split(" ");
			if (tokens[0].contains("rand") || tokens[0].contains("chrM"))
				continue;
			String chr =tokens[0];
			Integer size= Integer.valueOf(tokens[1]);
			
			Alignments aln= new Alignments(chr,1,size);
			aln.setName(chr);
			IntervalTree<Alignments> alnTree= new IntervalTree<Alignments>();
			alnTree.put(1,size, aln);
			intergenicGenome.put(chr,alnTree);
		}
		
	}
	
	
	//takes the intergenic genome and remove input genic regions
	private void removeRegionFromIntergenicGenome(BEDFileParser refGenes) {

		Iterator <String> chrIt=refGenes.getChromosomeIterator();
		while(chrIt.hasNext()){
			String chr=chrIt.next();
			if (! intergenicGenome.containsKey(chr))
				continue;
			IntervalTree<GeneWithIsoforms> refTree=refGenes.getChrTree(chr);
			IntervalTree<Alignments> intergenicTree= this.intergenicGenome.get(chr);
			
			Iterator<GeneWithIsoforms> refIt=refTree.valueIterator();
			while (refIt.hasNext()){
				
				for (Gene currentElement : refIt.next().getAllIsoforms()){
					int c_st=currentElement.getStart();
					int c_end=currentElement.getEnd();
					Iterator<Alignments> overlapperIt = new IntervalTree.ValuesIterator<Alignments>(intergenicTree.overlappers(c_st,c_end));
					while(overlapperIt.hasNext()){
						Alignments overlapper = overlapperIt.next();
						//If all is with in the removed region, remove all
						if (overlapper.getStart()>= c_st && overlapper.getEnd()<=c_end )
							intergenicTree.remove(overlapper.getStart(),overlapper.getEnd());
						//If the gene is completly within the element, split the element
						if (overlapper.getStart()<= c_st ||  overlapper.getEnd()>=c_end ){
							intergenicTree.remove(overlapper.getStart(),overlapper.getEnd());
							if ( overlapper.getStart()< c_st ){
								Alignments stAln = new Alignments(chr,overlapper.getStart(),c_st);
								intergenicTree.put(overlapper.getStart(),c_st,stAln);
							}
							if ( overlapper.getEnd()> c_end ){
								Alignments endAln = new Alignments(chr,c_end,overlapper.getEnd());
								intergenicTree.put(c_end,overlapper.getEnd(),endAln);								
							}
							
						}
						intergenicGenome.put(chr,intergenicTree);
					}
					
					
				}
				
			}
		}
	}
	
				
	private void makeProbFractionedIntergenicGenome() {

		int totalSize=0;
		//go through the intergenic genome and calc sizes
		for (String chr: intergenicGenome.keySet()){
			int chr_size = calcChrTreeSize(intergenicGenome.get(chr));
			chrSizes.put(chr,chr_size);
			totalSize+= chr_size;
		}
		
		GenomeSize=totalSize;
		
		//make genome mapped by probability of each chr
		int st=0;
		int end=0;
		String fchr="";
		for (String chr: intergenicGenome.keySet()){
			int chr_size = chrSizes.get(chr);
			int prob =  (int)((new Double (chr_size) / new Double(totalSize)) * SCALE_FACTOR);
			st=end+1;
			end = Math.max( st + (prob-1), st+1);
			this.probFractionedChrTree.put(st, end, new Alignments(chr,st,end));
			fchr=chr;
		}
		this.probFractionedChrTree.remove(st, end);
		this.probFractionedChrTree.put(st, SCALE_FACTOR, new Alignments(fchr,st,SCALE_FACTOR));
		
		//make chrs mapped by probability of each interval
		for (String chr: intergenicGenome.keySet()){
			int chr_size = chrSizes.get(chr);
			Iterator <Alignments> chrIt = intergenicGenome.get(chr).valueIterator();
			 st=0;
			 end= 0;
			 IntervalTree<Node<Alignments>> tree = new IntervalTree<Node<Alignments>>();
			while (chrIt.hasNext()){
				Alignments aln=chrIt.next();
				int alnSize=aln.getEnd()-aln.getStart();
				st=end+1;
				end=st+alnSize-1;
				tree.put(st,end,new Node(st,end,aln));
			}
			System.err.println( "Chr is : " +chr + " lastIndex: " + end + " chrSize: "+ chr_size);
			this.probFractionedIntergenicGenome.put(chr,tree);
		}
		
	}
	
	


	private int calcChrTreeSize(IntervalTree<Alignments> alnChrTree) {
		
		int res=0;
		Iterator <Alignments> chrIt = alnChrTree.valueIterator();
		while (chrIt.hasNext()){
			Alignments aln=chrIt.next();
			res += aln.getEnd()-aln.getStart();
		}
		
		return res;
	}




    public void makeRandTranscriptList(int numRand, BEDFileParser ref){
    	
    	for (int i=1; i<=numRand; i++){
    		
    		if (i%100 == 0)
    			System.err.println("Making random set number " +i);
    		
    		BEDFileParser randSet= new BEDFileParser();
    		for (Gene g:  ref.GetGenes()){
    			Iterator<Node<Alignments>> chrIt =null;
    			boolean gotChr=false;
    			while(gotChr==false) {
    				Double randVal= Math.random();
    				int randChrIx = (int) (randVal* (SCALE_FACTOR-2));
    				//System.err.println(randChrIx);
    				 chrIt=probFractionedChrTree.overlappers(randChrIx,randChrIx+20);
    				 if (chrIt.hasNext()==true)
    					 gotChr=true;
    				 else
    					 System.err.println ("Did not find chr with probability range overlapping " +randVal + "  ("+randChrIx + ")");
    			}
    			
    			
    			String chr = chrIt.next().getValue().getChr();
    			int chr_size = chrSizes.get(chr);
    			IntervalTree<Node<Alignments>> chrTree =  this.probFractionedIntergenicGenome.get(chr); 
    			boolean addedGen=false;
    			Node<Alignments> segment = null; 
    			Iterator<Node<Node<Alignments>>> segmentIt=null;
    			int randIx=0;
    			while(addedGen==false){
    			
    				boolean gotSegment=false;
    				while (gotSegment==false) {
    					 randIx= (int) (Math.random()*(chr_size-2));
    					 segmentIt=chrTree.overlappers(randIx,randIx+20); 
    					 if (segmentIt.hasNext())
    						 gotSegment=true;
    					 else
    						 System.err.println ("Did not find segment with probability range overlapping " +randIx + " while chr size is: ("+ chr_size + ")" +  " chr Name:" + chr);
    				    	 
    				}
    				segment = segmentIt.next().getValue();
	    			int segentSize= segment.getEnd()-segment.getStart();
	    			int offset = randIx-segment.getStart();
	    			int delta = segment.getEnd() - randIx;
	    			if (delta >= g.getGenomicLength()){
	    				Alignments realSegment = segment.getValue();
	    				int startPos = realSegment.getStart() + offset; 
	    				Collection<Annotation> new_exons = new TreeSet<Annotation>();
	    				Collection<? extends Annotation> old_exons = g.getExonSet();
	    				for (Annotation old : old_exons) {
	    					
	    					int nst = startPos+ (old.getStart() - g.getStart());
	    					int nend =startPos+ (old.getEnd() - g.getStart());
	    					new_exons.add(new BasicAnnotation(chr,nst,nend));
	    				
	    				}
	    				Gene newGen = new Gene (new_exons);
	    				newGen.setName(g.getName());
	    				newGen.setOrientation (g.getOrientation());
	    				randSet.addRefSeq(newGen);
	    				addedGen=true;
	    			}
    			}
    			
    		}
    		randTranscriptsList.add(randSet);
    	}
    	
    }




	class chrSizeNode{
		String chr;
		int start;
		int stop;
		
		public chrSizeNode (String chr1,int st,int en){
			chr=chr1;
			start=st;
			stop=en;
		}
		
		public String getChr(){
			return chr;
		}
	}




	public void printRandModels(String outprefix) throws IOException {

		for (int i=0; i<randTranscriptsList.size(); i++)
		{
			BEDFileParser bed = randTranscriptsList.get(i);
			bed.writeFullBed(outprefix+"rand_"+i+".bed");
		}
		
	}

	public void writeIntergenicGenome(String outf) throws IOException {

		BufferedWriter bw = new BufferedWriter( new FileWriter(outf));
		for(String chr: intergenicGenome.keySet()){
			Iterator<Alignments> it =  intergenicGenome.get(chr).valueIterator();
			while(it.hasNext()){
				Alignments aln= it.next();
				Gene g = new Gene(aln);
				bw.write(g.toBED()+"\n");
			}
			
		}
		bw.close();
	}

	public ArrayList<BEDFileParser> getRandomSets() {
		
		return randTranscriptsList;
	}

	public void CrossRandModelsWithBed(String bedToCrossWith_fname, BEDFileParser refBed, String outfile) throws IOException {

		BEDFileParser bedToCrossWith = new BEDFileParser(bedToCrossWith_fname);
		//Count overlaps in ref bed
		BEDFileParser mergedBed =  refBed.getMergedCopy() ;
		int refBedCnt = countExonLevelOverlap (mergedBed,bedToCrossWith );
		BufferedWriter bw = new BufferedWriter(new FileWriter(outfile));
		
		double Cntr=0.0;
		double tot=0.0;
		int i=0;
		String [] val1 = new String[this.getRandomSets().size()];
		for (BEDFileParser b : this.getRandomSets()){
			
			BEDFileParser b_merged= b.getMergedCopy();
    		int cntRand=countExonLevelOverlap(b_merged,bedToCrossWith);
    		tot++;
    		if (cntRand >= refBedCnt & cntRand >0)
    			Cntr++;
    		bw.write(String.valueOf(cntRand)+"\n");
    		
    		
		}
		
		System.out.println("OverlapWith refbed: " + refBedCnt );
		System.out.println("Perm Pval: " + Cntr/tot + " numOfPerm: " +tot );
		bw.close();
		
	}
	
	
	private static int countExonLevelOverlap(BEDFileParser geneMergeSet,
			BEDFileParser elements) {
		int cnt=0;
		//count the number of genes that overlap the element across their exons
		for (Gene g: geneMergeSet.GetGenes()){
			IntervalTree<GeneWithIsoforms> overlap = elements.getOverlappers(g);
			if (! overlap.isEmpty())
				cnt++;
		}
		return cnt;
	}
	
}
