package broad.pda.differentialExpression;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.core.datastructures.MatrixWithHeaders;
import broad.core.error.ParseException;
import broad.core.util.GMTParser;
import broad.core.util.ParseGCTFile;
import broad.pda.geneexpression.ExpressionExperimentInfo;
import broad.pda.geneexpression.agilent.AgilentUtils;
import broad.projection.gsa.GeneSetEnrichment;

public class RunGSEA {
	
	int minNum=5;
	private double alphaCutoff=.01;

	public RunGSEA(File gctFile, File experimentInfoFile, File chipFile, File gmtFile, String save) throws IOException, ParseException{
		MatrixWithHeaders data=new MatrixWithHeaders(gctFile.getAbsolutePath(), chipFile.getAbsolutePath());
		Map<String, ExpressionExperimentInfo> experimentInfo=AgilentUtils.parseExperimentInfoFile(experimentInfoFile);
		Map<String, Collection<String>> geneSets=GMTParser.ParseGeneSetFile(gmtFile, minNum);
		Map<String, Collection<String>> chipInfo=GMTParser.parseCHIPFileByName(chipFile);
		
		geneSets=convertGeneSets(geneSets, data, chipInfo, minNum);
		System.err.println("Finished converting gene lists");
		
		runGeneSetComparisons(data, experimentInfo, geneSets, save);
	}
	
	private Map<String, Collection<String>> convertGeneSets(Map<String, Collection<String>> geneSets, MatrixWithHeaders data, Map<String, Collection<String>> chipInfo, int minNum2) {
		Map<String, Collection<String>> rtrn=new TreeMap<String, Collection<String>>();
		
		for(String geneSet: geneSets.keySet()){
			Collection<String> genes=geneSets.get(geneSet);
			Collection<String> pids=convert(genes, data, chipInfo);
			if(pids.size()>minNum){rtrn.put(geneSet, pids);}
		}
		
		return rtrn;
	}

	private Collection<String> convert(Collection<String> genes, MatrixWithHeaders data, Map<String, Collection<String>> chipInfo) {
		//First check if the current gene names are in the row names
		//if so return itself
		//Else convert
		Collection<String> rtrn=new TreeSet<String>();
		
		int found=0;
		int notFound=0;
		
		for(String gene: genes){
			if(data.hasRow(gene)){rtrn.add(gene); found++;}
			else{
				Collection<Integer> indexes=data.getIndecesForRowDescription(gene);
				if(indexes!=null){
					for(Integer index: indexes){rtrn.add(data.getRowName(index));}
					found++;
				}
				else if(chipInfo.get(gene)!=null){
					Collection<String> pids=chipInfo.get(gene);
					for(String pid: pids){
						if(data.containsRow(pid)){rtrn.add(pid);} //Only if the gene is actually in the data matrix do we want to add it
						found++;
					}
				}
				else{notFound++;}
			}
			
		}
		//System.err.println("Found "+(((double)found/(found+notFound)))*100+"%");
		
		return rtrn;
	}

	private void runGeneSetComparisons(MatrixWithHeaders data, Map<String, ExpressionExperimentInfo> experimentInfo, Map<String, Collection<String>> geneSets, String save) throws IOException {
		Collection<String> negatives=new TreeSet<String>();
		Map<String, Collection<String>> samplesByName=new TreeMap<String, Collection<String>>();
				
		for(String barcode: data.getColumnNames()){
			if(experimentInfo.containsKey(barcode)){
				ExpressionExperimentInfo info=experimentInfo.get(barcode);
				String group=info.getSampleType();
				String name=info.getExperimentName();
				Collection<String> nameSamples=new ArrayList();
				if(samplesByName.containsKey(name)){nameSamples=samplesByName.get(name);}
				nameSamples.add(barcode);
				if(group.equalsIgnoreCase("Control")){negatives.add(barcode);}
				else{samplesByName.put(name, nameSamples);}
			}
		}
		
		Collection<String> genes=new TreeSet<String>();
		
		//TODO Write a GCT with all the enrichments for the significant gene sets
		MatrixWithHeaders maxMean=new MatrixWithHeaders(new ArrayList(geneSets.keySet()), new ArrayList(samplesByName.keySet()));
		MatrixWithHeaders KS=new MatrixWithHeaders(new ArrayList(geneSets.keySet()), new ArrayList(samplesByName.keySet()));
		
		for(String name: samplesByName.keySet()){
			System.err.println(name);
			GeneSetEnrichment diff=new GeneSetEnrichment(data, negatives, samplesByName.get(name), geneSets, true, true);
			diff.getNormalizedMaxMeanEnrichments().writeGCT(save+"."+name+".normalizedMaxMean.gct");
			diff.getNormalizedKSEnrichments().writeGCT(save+"."+name+".normalizedKS.gct");
			diff.getMaxMeanFDR().writeGCT(save+"."+name+".maxMeanFDR.gct");
			diff.getKSFDR().writeGCT(save+"."+name+".KSFDR.gct");
			
			for(String geneSet: geneSets.keySet()){
				double KSFDR=diff.getKSFDR().get(geneSet, 0);
				double MMFDR=diff.getMaxMeanFDR().get(geneSet, 0);
				double MMEnrich=diff.getNormalizedMaxMeanEnrichments().get(geneSet, 0);
				double KSEnrich=diff.getNormalizedKSEnrichments().get(geneSet, 0);
				if(MMFDR>alphaCutoff){MMEnrich=0;}
				if(KSFDR>alphaCutoff){KSEnrich=0;}
				maxMean.set(geneSet, name, MMEnrich);
				KS.set(geneSet, name, KSEnrich);
			}	
		}
		
		maxMean.writeGCT(save+".MaxMean.gct");
		KS.writeGCT(save+".KS.gct");
	}
	
	
	public static void main(String[] args)throws IOException, ParseException{
		if(args.length>5){
		File gctFile=new File(args[0]);
		File experimentInfo=new File(args[1]);
		File chipFile=new File(args[2]);
		String save=args[3];
		boolean preprocess=new Boolean(args[4]);
		File gmtFile=new File(args[5]);
		new RunGSEA(gctFile, experimentInfo, chipFile, gmtFile, save);
		}
		else{System.err.println(usage);}
	}

	static String usage=" args[0]=gct file \n args[1]=experimentInfo \n args[2]=chip file \n args[3]=save \n args[4]=preprocess \n args[5]=gmt file";
}
