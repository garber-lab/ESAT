package broad.pda.gene;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Set;
import java.util.TreeSet;

import nextgen.core.annotation.Gene;

import org.apache.log4j.Logger;

import broad.core.math.Statistics;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.pda.annotation.BEDFileParser;
import broad.pda.gene.GeneTools;

public class functionalAnnotation {
	
	HashMap <String, HashMap<String,String>> goAnnotMap;
	HashMap <String,Integer> GoSizeMap;
	HashMap <String,Integer> annotatedGene;
	
	HashMap <String,String> geneSetNames;
	LinkedList <HashMap<String,String>> randomSetNames;
	HashMap <String,Double> termEnrichmentPermPval;
	HashMap <String,Double> termEnrichmentFischerPval;
	HashMap <String,Double> termEnrichmentFDRCorrectedPval;
	HashMap <String,Integer> termOverlap;
	
	static Logger logger = Logger.getLogger(GeneTools.class.getName());
	
	public functionalAnnotation () {
		this.goAnnotMap=new HashMap <String, HashMap<String,String>>();
		this.GoSizeMap= new HashMap <String,Integer>();
		this.annotatedGene=new HashMap <String,Integer>();
		this.geneSetNames =new HashMap<String,String> ();
		this.randomSetNames=new LinkedList <HashMap<String,String>> ();
		this.termEnrichmentPermPval=new HashMap <String,Double> ();
		this.termEnrichmentFDRCorrectedPval=new HashMap <String,Double> ();
		this.termEnrichmentFischerPval=new HashMap <String,Double> ();
		this.termOverlap=new HashMap <String,Integer> ();
		
	}
	
	public void loadGoMap(String gofile,String goformat, int goSizeThreshold) throws IOException{
		
		if(goformat.equals("davidOut"))
			loadGoMap_davidOut(gofile,goSizeThreshold);
		if(goformat.equals("go2nm"))
			loadGoMap_go2nm(gofile,goSizeThreshold);
	}
	
	private void loadGoMap_david(String gofile) throws IOException {

		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(gofile)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
			String[] tokens=nextLine.split("\t");
			String term=tokens[1];
			String id=tokens[0];
			
			if(!this.annotatedGene.containsKey(id))
				this.annotatedGene.put(id, 0);
			this.annotatedGene.put(id,this.annotatedGene.get(id)+1);
		
			if (! this.GoSizeMap.containsKey(term)){
				this.GoSizeMap.put(term,0);
				this.goAnnotMap.put(term,new HashMap<String,String>());
			}
				
			this.GoSizeMap.put(term,this.GoSizeMap.get(term)+1);
			this.goAnnotMap.get(term).put(id,"");
		}
		
	}

	private void loadGoMap_go2nm(String gofile, int goSizeThreshold) throws IOException {

		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(gofile)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
			String[] tokens=nextLine.split("\t");
			String term=tokens[0];
			String name=tokens[1];
			String finalTerm=term+"_"+name;
			String ids=tokens[3];
			String[] names= ids.split(",");
			if (names.length >= goSizeThreshold)
				loadTerm(finalTerm,names);
		}
		
		
	}
	
	private void loadGoMap_davidOut(String gofile, int goSizeThreshold) throws IOException {

		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(gofile)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
			String[] tokens=nextLine.split("\t");
			String term=tokens[1];
			String ids=tokens[5];
			String[] names= ids.split(",");
			if (names.length >= goSizeThreshold)
				loadTerm(term,names);
		}
		
		
	}
	private void loadTerm(String term, String[] names){
		HashMap<String,String> idHash=new HashMap<String,String>();
		for (int i=0; i<names.length; i++){
			String id=names[i];
			idHash.put(id,"");
			if(!this.annotatedGene.containsKey(id))
				this.annotatedGene.put(id, 0);
			this.annotatedGene.put(id,this.annotatedGene.get(id)+1);
		}
		this.GoSizeMap.put(term,names.length);
		this.goAnnotMap.put(term,idHash);
	}
	
	public void calcAllTermsEnrichment(Collection<Gene> geneSet, LinkedList<Collection<Gene>> randomGeneSets ) {
		
		this.geneSetNames=makeAnnotatedNameMap(geneSet);
		LinkedList<HashMap<String,String>> randSetNames=makeAnnotatedNameMap(randomGeneSets);
		for (String term: this.goAnnotMap.keySet()){
			Double[] permPval=new Double[1];
			Double[] fischerExact=new Double[1];
			Integer[] numOverlap=new Integer[1];
			LinkedList<Double> randomFischerExact=new LinkedList<Double>();
			LinkedList<Integer> randomOverlap=new LinkedList<Integer>();
		
			calcTermEnrichment(term,this.geneSetNames,randSetNames,permPval,fischerExact,numOverlap,randomFischerExact,randomOverlap);
		
			this.termOverlap.put(term,numOverlap[0]);
			this.termEnrichmentFischerPval.put(term,fischerExact[0]);
			this.termEnrichmentPermPval.put(term,permPval[0]);
			
		}
	}
	
	
 // Include only the elements f the gene set that are annotated by some term 
	private HashMap<String, String> makeAnnotatedNameMap(Collection<Gene> geneSet) {
		
		HashMap<String, String> res=new HashMap<String, String>();
		for (Gene g:geneSet){
			String name=g.getName();
			if(this.annotatedGene.containsKey(name))
				res.put(name, "");
		}
		return res;
	}
	
	private LinkedList<HashMap<String, String>> makeAnnotatedNameMap(LinkedList<Collection<Gene>> randomGeneSets) {

		LinkedList<HashMap<String, String>> res= new LinkedList<HashMap<String, String>>();
		for (int i=0; i<randomGeneSets.size(); i++){
			res.add(makeAnnotatedNameMap(randomGeneSets.get(i)));
		}

		return res;
	}

	//Input: a gene set and a list of random gene sets and a functional term
	//Description: the function calculates the overlap of the gene set with the term's annotated set
	//and than calculates a pval for the significance of this overlap , by comparing the number of overlaps
	//in the random sets. It also calculates a fischer exact 
	//Output: updates the the permPval,fischerExact and fischerExactRandom vars 
	public void calcTermEnrichment(String term, HashMap<String,String> geneSet, LinkedList<HashMap<String,String>> randomSets,
			Double[] permPval,Double[] fischerExact, Integer[] numOverlap,
			LinkedList<Double> randomFischerExact,LinkedList<Integer> randomOverlap ) {
		
		numOverlap[0]=new Integer(calcOverlapWithTerm(geneSet,term));
		fischerExact[0]=new Double(fischerExactTest(geneSet.size(),numOverlap[0], this.annotatedGene.size(),this.GoSizeMap.get(term)));
		double hits=0.0;
		for (int i=0; i<randomSets.size(); i++){
			HashMap<String,String> ref= randomSets.get(i);
			 double[] randRes= new double[2];
			 int randOver=calcOverlapWithTerm(ref,term);
			 //double randFisch=fischerExactTest(ref.size(),randOver, this.annotatedGene.size(),this.GoSizeMap.get(term));
			 if (randOver>=numOverlap[0])
				 hits++;
			 //randomFischerExact.add(randFisch); 
			 randomOverlap.add(randOver);
		}
		permPval[0]=new Double(hits/new Double(randomSets.size()));
		

	}

	//at the moment does a hypergeometric
	private double fischerExactTest(int setSize, int setHits, int populationSize,
			int termSize) {
		double p= Statistics.hypergeometric(populationSize,  setSize,termSize,  setHits);
		if (new Double(p).isNaN())
			System.err.println("P=NaN: "+ populationSize +", "+  setSize+", "+ termSize+", "+ setHits);
		return p;
	}

	private int calcOverlapWithTerm(HashMap<String, String> idSet,String term) {

		int hits=0;
		HashMap<String,String> termIds=this.goAnnotMap.get(term);
		for (String id:idSet.keySet()){
			if (termIds.containsKey(id))
				hits++;
		}
		return hits;
	}
	
	
	public void writeTermsEnrichments(String filename,int sizeThreshold, String str,ArrayList<String> termsSorted) throws IOException{
		BufferedWriter bw=new BufferedWriter(new FileWriter(filename));
		bw.write(str+"\n");	
		bw.write("Term\tPopulationSize\tTermSetSize\tGeneSetSize\tGeneTermOverlap\tFischerExact");
		if (! this.termEnrichmentPermPval.isEmpty())
			bw.write("\tPermPval");
		if (! this.termEnrichmentFDRCorrectedPval.isEmpty())
			bw.write("\tFDRCorrectedPval");
		bw.write("\n");
		
		//bw.write("AnnotatedGenes:\t"+this.annotatedGene.size()+"\tGeneSetSize:\t"+this.geneSetNames.size()+"\n");
		int popSize=this.annotatedGene.size();
		int geneSize=this.geneSetNames.size();
		
		if (termsSorted==null)
			termsSorted=new ArrayList<String>(); termsSorted.addAll(this.goAnnotMap.keySet());

		for (int i=0; i<termsSorted.size();i++){
			String term=termsSorted.get(i);
			if (this.goAnnotMap.get(term).size() >= sizeThreshold & this.termOverlap.get(term) >0 ){
				bw.write(term+"\t"+popSize+"\t"+geneSize+"\t"+this.goAnnotMap.get(term).size()+"\t");
				bw.write(this.termOverlap.get(term)+"\t");
				bw.write(String.valueOf(this.termEnrichmentFischerPval.get(term)));
				if (! this.termEnrichmentPermPval.isEmpty())
					bw.write("\t"+this.termEnrichmentPermPval.get(term));
				if (! this.termEnrichmentFDRCorrectedPval.isEmpty())
					bw.write("\t"+this.termEnrichmentFDRCorrectedPval.get(term));
				bw.write("\n");
			}
		}
		
		bw.close();
		
	} 
	
	
	//Given a:
	// 1) Ontology file .obo file for ontology description  2) Annotation file go-PID 
	// 3) PID->gene ID mapping  4) NM -> gene ID mapping  
	public static void parseGoRawFiles(String Ontology, String Annotation, String pid2gid, String nm2pidF,int nmCol, int pidCol, String outf) throws IOException{
		
		//HashMap<String, Set<String>> nm2pid= parseTable(nm2pidF,nmCol,pidCol);
		HashMap<String, Set<String>> pid2nm= parseTable(nm2pidF,pidCol,nmCol);
		HashMap <String,String> goCode2Name= parseOntologyFile(Ontology);
		
		HashMap <String,Set<String>> goCode2nm= parseGoAnntFile(Annotation,pid2nm);
		
		BufferedWriter bw= new BufferedWriter(new FileWriter(outf));
		for (String term: goCode2nm.keySet()){
			String annot="";
			if (goCode2Name.containsKey(term) )
				annot= goCode2Name.get(term);
			if (goCode2nm.containsKey(term)){
				Set<String> s=goCode2nm.get(term);
				int size= s.size();
				if (size>0){
					bw.write(term+"\t"+annot+"\t"+size+"\t");
					for (String name:s)
						bw.write(name+",");
					bw.write("\n");
				}
			}
		}
		bw.close();
		return ;
	}

	private static HashMap<String, Set<String>> parseGoAnntFile(
			String annotation, HashMap<String, Set<String>> pid2nm) throws IOException {
		
		HashMap<String, Set<String>> res= new HashMap<String, Set<String>>();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(annotation)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
			if (nextLine.contains("!"))
				continue;
			String[] tokens=nextLine.split("\t");
			String term=tokens[4];
			String pid=tokens[1];
			if (! res.containsKey(term))
				res.put(term, new TreeSet<String>());
			if (pid2nm.containsKey(pid))
				res.get(term).addAll(pid2nm.get(pid));
		}
		
		return res;
	}

	private static HashMap<String, String> parseOntologyFile(String ontology) throws IOException {


		HashMap<String, String> res= new HashMap<String, String>();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(ontology)));
		String nextLine;
		while ((nextLine = reader.readLine()) != null ) {
			if (nextLine.contains("[Term]")){
				String id= reader.readLine();//  id: GO:0000001
				String name= reader.readLine(); //name: mitochondrion inheritance
				String namespace= reader.readLine();  // namespace: biological_process
				id=id.replace("id: ", "");  
				name=name.replace("name: ","");
				namespace=namespace.replace("namespace: ", "");
				res.put(id,name);
			}
		}
			
		return res;
	}

	private static HashMap<String, Set<String>> parseTable(String filename,
			int c1, int c2) throws IOException {

		HashMap<String, Set<String>> res= new HashMap<String, Set<String>>();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(filename)));
		String nextLine;
		int top=Math.max(c1, c2);
		while ((nextLine = reader.readLine()) != null && (nextLine.trim().length() > 0)) {
			String[] tokens=nextLine.split("\t");
			if (tokens.length>=top){
				String n1=tokens[c1-1];
				String n2=tokens[c2-1];
				if (n1.equals("") || n2.equals("") )
					continue;
				if (! res.containsKey(n1))
					res.put(n1, new TreeSet<String>() );
				res.get(n1).add(n2);
			}
		}
		return res;
	}

	public Set<Gene> getTermGenes(String goTerm, BEDFileParser ref) {

		HashMap<String, Gene> nameGeneMap=ref.getNameGeneMap();
		Set<Gene> res= new TreeSet<Gene>();
		goTerm=goTerm.replace("\"" ,"");
		if (this.goAnnotMap.containsKey(goTerm)){
			for (String name : this.goAnnotMap.get(goTerm).keySet()){
				if (nameGeneMap.containsKey(name) & nameGeneMap.get(name) != null)
					res.add(nameGeneMap.get(name));
			}
		}
		return res;
	}
	
	private static void GOanalysis(String geneSetFile, String gofile,String goformat, String outDir,int goSizeThreshold) throws IOException {
		
		//File file= new File(geneSetFile);
		//String name= file.getName();
		//System.err.println(name);
		String[] fname=geneSetFile.split("\t");
		functionalAnnotation GO= new functionalAnnotation();
		
		GO.loadGoMap(gofile, goformat,goSizeThreshold);
		GO.calcTermEnrichments(fname[0]);
		ArrayList<String> termsSorted=GO.applyFDRCorrection(0.05);
		
		int covered= GO.numOfTotalCovered();
		int FDRCovered= GO.numOfFDRCovered();
		
		String str = fname[1]+" total covered :" +covered + " covered with FDR significant term : " + FDRCovered;
		
		
		String outf =  outDir  + fname[1] +"GoEnrich"  ;
		System.err.println(outf);
		
		GO.writeTermsEnrichments(outf,goSizeThreshold,str,termsSorted);
		
	}

	private int numOfFDRCovered() {
		HashMap<String, Integer> lst=new HashMap<String, Integer>();
		for (String term:this.termEnrichmentFDRCorrectedPval.keySet())
		{
			if (this.termEnrichmentFDRCorrectedPval.get(term)==1){
				for (String gene : this.goAnnotMap.get(term).keySet()){
					if (! lst.containsKey(gene))
						lst.put(gene,new Integer(0));
					lst.put(gene,(lst.get(gene)+1));
				}	
			}
		}
		return crossGeneSetWithList(lst);
	}

	private int numOfTotalCovered() { return crossGeneSetWithList(this.annotatedGene);}
	

	private int crossGeneSetWithList(HashMap<String, Integer> lst){
		int cnt=0;
		for (String gene: this.geneSetNames.keySet()){
			if (lst.containsKey(gene))
				cnt++;
		}
		return cnt;
	}

	private ArrayList<String> applyFDRCorrection(double alpha) {
		
		ArrayList<Double> pvals= new ArrayList<Double> ();
		System.err.println("size of pvals "+pvals.size());
		pvals.addAll(this.termEnrichmentFischerPval.values());
		ArrayList<String> g1=new ArrayList<String>();
		ArrayList<String> g2=new ArrayList<String>();
		double T=Statistics.FDRCorrect(pvals, alpha);
		for (String term: this.termEnrichmentFischerPval.keySet()){
			 double v= this.termEnrichmentFischerPval.get(term);
			Double pass=new Double(0.0);
			if (v<T){
				pass=1.0;
				g1.add(term);
			}
			else
				g2.add(term);
			this.termEnrichmentFDRCorrectedPval.put(term,pass);
			
		}
		
		g1.addAll(g2);
		return g1;
	}

	private static void GOanalysisBatch(String fileLst, String annot, String goFormat,
			String outDir,int goSizeThreshold) throws IOException {
		
		BufferedReader br = new BufferedReader(new FileReader(fileLst));
		String line;
		while((line = br.readLine()) != null) {
			line = line.trim();
			GOanalysis(line, annot,goFormat, outDir,goSizeThreshold) ;
		}		
	}

	
	private void calcTermEnrichments(String geneSetFile) throws IOException {

		loadGeneList(geneSetFile);
		int geneSetSize=this.geneSetNames.size();
		for (String term :this.goAnnotMap.keySet()) {
			Integer numOverlap=new Integer(calcOverlapWithTerm(this.geneSetNames,term));
			double fischerExact=new Double(fischerExactTest(geneSetSize,numOverlap, this.annotatedGene.size(),this.GoSizeMap.get(term)));
			this.termOverlap.put(term,numOverlap);
			this.termEnrichmentFischerPval.put(term,fischerExact);
		}
		
	}



	private void loadGeneList(String geneSetFile) throws IOException {
		
		BufferedReader br = new BufferedReader(new FileReader(geneSetFile));
		String line;
		while((line = br.readLine()) != null) {
			line = line.trim();
			this.geneSetNames.put(line,"" );
		}
	}

	private  void getTermGenes(String annot,String goformat, String interm,BufferedWriter bw) throws IOException {
		
		loadGoMap(annot, goformat,1);
		Set<String> genes= new HashSet<String>();
		for (String term:this.goAnnotMap.keySet() ){
			//String term=this.geneSetNames.get(termNum);
			if(term.contains(interm))
				genes.addAll(this.goAnnotMap.get(term).keySet());
		}
		for (String g:genes)
			bw.write(g+"\n");
	}


	static String usage="Usage: GeneTools -task <task name> "+
	"\n\t MakeGo2NMtable: process GO annot raw file to a table that associates GO terms with a comma seperated list of NM ids. \n\t\t -ontology <obo file>\n\t\t -annot <gene association file> \n\t\t -pid2nm <biomart table> \n\t\t -outf <out file name> "+
	"\n\t GOanalysis: given list of gene set file/ a gene set file , find enrichment of go terms: \n\t\t -fileLst <file with a list of gene set files> \n\t\t -annot <GO2NMtable; see tool above> \n\t\t -out <out dir name> "+
	"\n\t TermGeneList: write all genes tat are associated with GOterms that contain the word specified by term : \n\t\t -term <matching word to select terms> \n\t\t -annot <GO2NMtable; see tool above> \n\t\t -out <out dir name> "+
	"\n";
	
	public static void main(String [] args) throws Exception  {
		ArgumentMap argmap = CLUtil.getParameters(args, usage);
		
		logger.info(argmap.getTask() +"\n");
		
		if("MakeGo2NMtable".equalsIgnoreCase(argmap.getTask())) {
			String Ontology=argmap.getMandatory("ontology");
			String annot=argmap.getMandatory("annot");
			String pid2nm=argmap.getMandatory("pid2nm");
			String out=argmap.getOutput();
			if  (argmap.containsKey("Entrez")){
				logger.info("entrez");
				parseGoRawFiles(Ontology, annot ,pid2nm,  pid2nm,4, 5, out);}
			else{
				logger.info("nm");
				parseGoRawFiles(Ontology, annot ,pid2nm,  pid2nm,3, 5, out);}
				
		}
		if ("GOanalysis".equalsIgnoreCase(argmap.getTask())){
			String fileLst=argmap.getMandatory("fileLst");
			String annot=argmap.getMandatory("annot");
			String outDir=argmap.getOutput();
			int goSizeThreshold=argmap.containsKey("goSizeThreshold")? argmap.getInteger("goSizeThreshold"):5;
			GOanalysisBatch(fileLst,annot,"go2nm",outDir,goSizeThreshold);
		}
		if ("TermGeneList".equalsIgnoreCase(argmap.getTask())){
			String annot=argmap.getMandatory("annot");
			String term  =argmap.getMandatory("term");
			BufferedWriter bw  = argmap.getOutputWriter();
			functionalAnnotation GO= new functionalAnnotation();
			GO.getTermGenes(annot,"go2nm",term,bw);	
			bw.close();
		}
		else {
			System.err.println(usage);
		}
	}

	

	
	
	
}