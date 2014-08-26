package broad.pda.annotation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import nextgen.core.annotation.Gene;

import org.apache.log4j.Logger;

import broad.core.annotation.GFF;
import broad.pda.datastructures.Alignments;
import broad.pda.gene.GeneTools;
import broad.pda.gene.GeneWithIsoforms;

public class GTFFileParser extends BEDFileParser {
	
	static Logger logger = Logger.getLogger(GeneTools.class.getName());
	
	public GTFFileParser() {
		super();
	}
	
	public GTFFileParser(String fileName) throws IOException {
		super();
		loadIsoformDataByChrToTree(new File(fileName));
	}	


	public GTFFileParser(String fileName, String chr) throws IOException {
		super();
		loadIsoformDataByChrToTree(new File(fileName),chr);
	}
	
	
	
	@Override
	protected void loadIsoformDataByChrToTree(File file, String chr) throws IOException{
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		int i = 0;
		String nextLine;
		// currentExons maintains the set of exons for the latest transcript that has been read
		Collection<Alignments> currentExons = new ArrayList<Alignments>();
		// currentExtraFields maintains the extra fields for the latest transcript that has been read
		Map<String, String> currentExtraFields = new HashMap<String, String>();
		// maintain the start and stop codons so you can figure out CDS and UTR
		int translationStart = -1;
		int translationEnd = -1;
		while ((nextLine = reader.readLine()) != null ) {
			//System.err.println("Line " + nextLine + " looks like data? " + looksLikeData(nextLine, chr));
			String [] data = nextLine.split("\t");
			if(looksLikeData(nextLine, data, chr) && (nextLine.trim().length() > 0)){
				//System.err.println("Processing line");
				String featureChromosome = data[0];
				int start = 0;
				int end   = 0;
				String orientation = "*";
				boolean isOldFormat = false;
				// the fields read from this line of the file
				Map<String, String> fieldMap = new HashMap<String, String>();
				fieldMap.put("source", data[1]); // maintains the source of the annotation
				if(data.length == 9) {
					start = Integer.parseInt(data[3]);
					end = Integer.parseInt(data[4]);
					orientation = data[6];
					String [] extraFields = data[8].split("\\s*;\\s*");
					for(String extraField : extraFields) {
						String [] fieldInfo = extraField.split("\\s+");
						fieldMap.put(fieldInfo[0], fieldInfo[1]);
					}
				} else if (data.length == 6 || data.length == 7) {
					//System.err.println("Old format");
					start = Integer.parseInt(data[1]);
					end = Integer.parseInt(data[2]);
					orientation = data[3];
					fieldMap.put("transcript_id", data[5]);
					fieldMap.put("gene_id", data[5]);
					fieldMap.put("gene_name", data[4]);
					isOldFormat = true;
				} else {
					throw new IOException("BAD LIne " + nextLine + " has neither  9 (proper GTF) nor 6  (Old GTF) fields");
				}
				
				// a new transcript is being encountered
				if(currentExtraFields.isEmpty() || !fieldMap.get("transcript_id").equals(currentExtraFields.get("transcript_id"))) {
					if(!currentExtraFields.isEmpty()) {
						// create a new gene from the last set of fields
						//if (currentExons.isEmpty()) {   // Special case to treat a new type of transcript in Gencode v13 files - a "gene" type line with no exons to follow
							//do nothing - these lines will be ignored for now
						//}
						//else {
							Gene g = createRefSeqFromRawGTF(currentExons, currentExtraFields);
							if(translationStart > -1 && translationEnd > -1) {
								g.setCDSRegion(Math.min(translationStart,translationEnd), Math.max(translationStart,translationEnd));
							} 
							else {
								// if there is no start and stop codon, set CDS boundaries to exon boundaries
								int firstPos = Integer.MAX_VALUE;
								int lastPos = 0;
								for(Alignments exon : currentExons) {
									int s = exon.getStart();
									int e = exon.getEnd();
									if(s < firstPos) firstPos = s;
									if(e < firstPos) firstPos = e;
									if(s > lastPos) lastPos = s;
									if(e > lastPos) lastPos = e;
								}
								if(firstPos > lastPos) {
									throw new IllegalArgumentException("Cannot assign valid CDS start and end positions for gene " + g.getName());
								}
								g.setCDSRegion(firstPos, lastPos);
							}
							//System.err.println("new transcript_id " + fieldMap.get("transcript_id") + " prior transcript_id " + currentExtraFields.get("transcript_id")  +" Created gene " + g.toBED());
							addRefSeq(g);
						//}
					}
					// now set currentExtraFields to the new transcript
					currentExtraFields = fieldMap;
					// reset exons and start, stop codons
					currentExons.clear();
					translationStart = -1;
					translationEnd = -1;

				}
				// -1 on start point because gff is 1 shifted while BED is 0 shifted.
         	   //On the other hand, gff is inclusive while bed is not-> do not subtract 1 from the end coordinate
				
				// Now that currentExtraFields and currentExons are up to date, if this is an exon, add the exon to currentExons
				if ("exon".equals(data[2]) || isOldFormat) {
					currentExons.add(new Alignments(featureChromosome, start-1, end, orientation));
					//System.err.println("\tgotexon " + (new Alignments(data[0], Integer.parseInt(data[3]),Integer.parseInt(data[4]),data[6])));
				}
				
				
				// translation start is first position of start codon
				if(data[2].equals("start_codon")) {
					if(orientation.equals("+")) translationStart = Integer.parseInt(data[3]);
					else translationStart = Integer.parseInt(data[4]);
				}
				
				// translation end is before first position of stop codon
				if(data[2].equals("stop_codon")) {
					if(orientation.equals("+")) translationEnd = Integer.parseInt(data[3]) - 1;
					else translationEnd = Integer.parseInt(data[4]) + 1;
				}
				
				
			}
			i++;
			if(i%100000==0){System.err.println("Processed " + i + " lines.");}
		}
		reader.close();
		Gene g = createRefSeqFromRawGTF(currentExons, currentExtraFields);
		addRefSeq(g);
		
	}


	// Creates a RefSeqGene with name, chr, orientation, exon starts, exon ends, start = min exon start, end = max exon end
	private static Gene createRefSeqFromRawGTF(
			Collection<Alignments> currentExons,
			Map<String, String> currentExtraFields) {
		// Name is either geneName_transcriptName or transcriptName
		String name = currentExtraFields.containsKey("gene_name") ? currentExtraFields.get("gene_name") +"_" + currentExtraFields.get("transcript_id"): currentExtraFields.get("transcript_id");
		Gene g = new Gene(currentExons);
		g.setName(name);
		for(String key : currentExtraFields.keySet()) {
			g.addAttribute(key,currentExtraFields.get(key));
		}
		return g;
	}

	@Override
	protected void loadIsoformDataByChrToTree(File file) throws IOException {
		loadIsoformDataByChrToTree(file, null);
	}
	
	private static boolean looksLikeData( String line, String[] splitLine, String chr) {
		return splitLine.length <= 7 ||( line.contains("gene_id") && line.contains("transcript_id") && (chr == null || line.startsWith(chr)));
	}
	
	
	//Returns a hashmap that maps between shared gene ids to transcripts
	public HashMap<String, ArrayList<String>> getGTFGeneTranscriptMap() {
		
		HashMap<String, ArrayList<String>> res=  new HashMap<String, ArrayList<String>>();
		for (Gene g: this.GetGenes()){
			
			String gid=g.getAttribute("gene_id");
			String tid= g.getAttribute("transcript_id");
			if (! res.containsKey(gid))
				res.put(gid,new ArrayList<String>());
			res.get(gid).add(tid);
		}

		return res;
	}
	
	public void write(BufferedWriter bw, String source) throws IOException {
		Iterator<String> chrIt=this.getChromosomeIterator();
		while(chrIt.hasNext()){
			String chr=chrIt.next();
			Iterator<GeneWithIsoforms> geneIt = getChrTree(chr).valueIterator();
			while(geneIt.hasNext()) {
				GeneWithIsoforms gene = geneIt.next();
				Collection<Gene> isoforms = gene.getAllIsoforms();
				for(Gene isoform: isoforms){
					bw.write(isoform.toCufflinksGTF(source, isoform.getName().replaceAll("_[0-9]$", ""), isoform.getName(), ""));
					bw.newLine();
				}
			}
		}
	}

	
	/**
	 * Write a GTF using only the attributes in extraFields  (not modifying the gene_id, transcript_id as GTF2BED does) 
	 * @param bw  BufferedWriter
	 * @param Src string to annotate src field
	 * @param  UseExtraFields - if false call gtf2bed ; if true assumes has gene id and transcript id as maybe others
	 * @throws Exception 
	 **/
	public void GTFWrite(BufferedWriter bw, String Src, boolean UseExtraFields) throws Exception {
		
		if (!UseExtraFields)
			bed2CuffGtf(bw,Src,false);
		else {
			for (Gene gi : this.GetGenes()){
				String str;
				if (gi.getAttribute("gene_id") ==null  || gi.getAttribute("transcript_id") ==null) {
					logger.debug(gi.getName() + "Has no gene_id , transcript ID attrs ; make sure your input has all attrs!");
					
					return;
				}
				str=gi.toCufflinksGTF(Src,"", "","");
				bw.write(str);
			}
		}
				
				
		
	
		
	}
	
	
	public  HashMap <String,Locus> getGeneIDIsoMap (){
		HashMap <String,Locus> geneMap = new HashMap <String,Locus>();
			
		for (Gene iso: GetGenes()){
			String gene_id="";
			if (iso.getExtraFields() != null)
				gene_id=GFF.getGTFAttr(iso.getExtraFields()[0], "gene_id");
			else
				gene_id=iso.getAttribute("gene_id");
			//System.err.println(iso.getName()+"\n");
			if (!geneMap.containsKey(gene_id))
				geneMap.put(gene_id,new Locus(new GeneWithIsoforms(iso)) );
			geneMap.get(gene_id).addIsoform(new GeneWithIsoforms(iso), "def");
		}
		
		return geneMap;
	}
	
}
