package broad.projection.gsa;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import broad.core.datastructures.MatrixWithHeaders;
import broad.core.error.ParseException;
import broad.core.math.Statistics;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;
import broad.core.util.GMTParser;
import broad.core.util.ParseGCTFile;



public class GsaUtils {
	
	static String usage="Usage: GeneTools -task <task name> "+
	"\n\tcontinuousGSEA: Parameters \n\t\t -exp <gene expression file ( gct supported only)> \n\t\t -pattern <reference expression profiles, GCT format > \n\t\t -patternType <gct / cls> \n\t\t -gmt <Gene set file, .gmt format> \n\t\t -chip <mapping of probeID-GeneSymbol, .chip format> \n\t\t -setSize <minimal gene-set size> \n\t\t -out <out file name> \n\n"+
	"\n\tsubmitBatchContinuousGSEA: Parameters \n\t\t -exp <gene expression file (gct supported only)> \n\t\t -clsDir <name of cls> \n\t\t -gmt <Gene set file, .gmt format> \n\t\t -chip <mapping of probeID-GeneSymbol, .chip format> \n\t\t -setSize <minimal gene-set size> \n\t\t -outpath <out Dir> \n\t\t -q<LSF q> \n\t\t -scriptPath <path for executable> \n\t\t -submitFile <submitCommandsFile> \n\n"+
	"\n\tretrieveContinuousBatchGSEA: Parameters \n\t\t -resPath <path to submitBatchContinuousGSEA result dir>  -out <outfile> -suffix <suffix of result files> \n\n" +
	"\n\tsubmitBatchGSEA_Tool: Parameters \n\t\t -exp <gene expression file (gct supported only)> \n\t\t -clsDir <name of cls root directory; e.g cls/0/file.cls> \n\t\t -gmt <Gene set file, .gmt format> \n\t\t -chip <mapping of probeID-GeneSymbol, .chip format> \n\t\t -setSize <minimal gene-set size> \n\t\t -outpath <out Dir> \n\t\t -q <LSF q> \n\t\t -scriptPath <path for executable> \n\t\t -numOfPerm <number of permutations> \n\n"+
	"\n\tretrieveBatchGSEA_Tool: Parameters \n\t\t -resPath <path to submitBatchGSEA_tool result dir>  -out <outfile name with no extention> -gmt <Gene set file, .gmt format> \n\n" + 	
	"\n\tgetCorrMatrix: Parameters \n\t\t -exp <gene expression file ( gct supported only)>  -pattern <reference expression profiles, GCT format >  -out <outfile name >  \n\n" ;
	
	
	
	public static void main(String[] args)throws Exception{
		ArgumentMap argmap = CLUtil.getParameters(args, usage, "continuousGSEA");
		if("continuousGSEA".equalsIgnoreCase(argmap.getTask())) {
			String expressionFile = argmap.getMandatory("exp");
			String pattern = argmap.getMandatory("pattern");
			String patternType = argmap.getMandatory("patternType");
			String geneSets = argmap.getMandatory("gmt");
			String chipFile = argmap.getMandatory("chip");
			int minSetSize = Integer.parseInt(argmap.getMandatory("setSize"));
			String out = argmap.getOutput();
						
			MatrixWithHeaders gseaRes=continuousGSEA(expressionFile,pattern,geneSets,chipFile,minSetSize,patternType,out);
			
		}
		
		else if("submitBatchContinuousGSEA".equalsIgnoreCase(argmap.getTask())) {
			
			String expressionFile = argmap.getMandatory("exp");
			String geneSets = argmap.getMandatory("gmt");
			String chipFile = argmap.getMandatory("chip");
			String minSetSize = argmap.getMandatory("setSize");
			String queue = argmap.getMandatory("q");
			String script = argmap.getMandatory("scriptPath");
			String saveDirR = argmap.getMandatory("outpath");
			String jobs=argmap.getMandatory("submitFile");
			
			Runtime run=Runtime.getRuntime();
			
			FileWriter writer=new FileWriter(jobs);
			
			
			File[] clsDirs=new File(argmap.getMandatory("clsDir")).listFiles();
			for (int d=0; d<clsDirs.length; d++)
			{
				String subdir=clsDirs[d].getName();
				File[] clsFiles=new File(clsDirs[d].getAbsolutePath()).listFiles();
				
				File saveDir= new File (saveDirR+"/"+subdir);
				boolean exists = saveDir.exists();
			    if (!exists) {	run.exec("mkdir "+ saveDir.getAbsolutePath());}
				
				for(int i=0; i<clsFiles.length; i++){
					String save=saveDir.getAbsolutePath()+"/"+clsFiles[i].getName()+".gsea";
					String junk=saveDir.getAbsolutePath()+"/"+clsFiles[i].getName()+".junk";
					//-R \"rusage[mem=4]\"
					
					String command="bsub -q "+queue+ " -o "+ junk + " -e "+ junk+"2" + "  java -jar -Xmx3000m "+script+" -task continuousGSEA -exp " +expressionFile+ " -pattern " + clsFiles[i].getAbsolutePath() + " -patternType cls "+ " -gmt " + geneSets + " -chip " + chipFile + " -setSize "+ minSetSize + " -out " + save  ;
					//System.err.println(command);
					writer.write(command + "\n");
					run.exec(command);  
				}
			}
			writer.close();
		}
		else if("retrieveContinuousBatchGSEA".equalsIgnoreCase(argmap.getTask())) {
			
			String resPath=argmap.getMandatory("resPath");
			String suffix=argmap.getMandatory("suffix");
			String out = argmap.getOutput();
			
			File[] resDirs=new File(resPath).listFiles();
			boolean first=true;
			MatrixWithHeaders fdr=null;
			for (int d=0; d<resDirs.length; d++)
			{
				File[] resFiles=new File(resDirs[d].getAbsolutePath()).listFiles();
				for(int i=0; i<resFiles.length; i++){
					File file=resFiles[i];
					if (file.getName().contains(suffix)){
						MatrixWithHeaders fdrRes=new MatrixWithHeaders(file.getAbsolutePath());
						if (first){ fdr=fdrRes; first=false;}
						else {fdr.appendColumns(fdrRes);}
					}
				}
			}
			fdr.writeGCT(out);
		}
		
		else if("submitBatchGSEA_Tool".equalsIgnoreCase(argmap.getTask())) {
			
			String expressionFile = argmap.getMandatory("exp");
			String geneSets = argmap.getMandatory("gmt");
			String chipFile = argmap.getMandatory("chip");
			String minSetSize = argmap.getMandatory("setSize");
			String queue = argmap.getMandatory("q");
			String script = argmap.getMandatory("scriptPath");
			String saveDirR = argmap.getMandatory("outpath");
			String numPerm = argmap.getMandatory("numOfPerm");
			String numPlot = argmap.getMandatory("numPlot");
						
			Runtime run=Runtime.getRuntime();
						
			File[] clsDirs=new File(argmap.getMandatory("clsDir")).listFiles();
			for (int d=0; d<clsDirs.length; d++)
			{
				String subdir=clsDirs[d].getName();
				File[] clsFiles=new File(clsDirs[d].getAbsolutePath()).listFiles();
				
				File saveDir= new File (saveDirR+"/"+subdir);
				boolean exists = saveDir.exists();
			    if (!exists) {	run.exec("mkdir "+ saveDir.getAbsolutePath());}
				
				for(int i=0; i<clsFiles.length; i++){
					String save=saveDir.getAbsolutePath()+"/"+clsFiles[i].getName();
					String junk=saveDir.getAbsolutePath()+"/"+clsFiles[i].getName()+".junk";
					
					String command="bsub -q "+queue+ " -o "+ junk + "  java -cp "+
					" /ahg/regev/users/nmcabili/Src/scripts/GuiltByAssocaition/gsea2.jar  -Xmx3000m  xtools.gsea.Gsea -res " +
					expressionFile + " -cls " +  clsFiles[i].getAbsolutePath() + " -gmx  " + geneSets +
					" -collapse true -mode Max_probe -norm meandiv -nperm 100 -permute gene_set -rnd_type no_balance -scoring_scheme weighted -rpt_label " + 
					clsFiles[i].getName() +  " -metric Pearson -sort real -order descending -chip " + chipFile +
					" -include_only_symbols true -make_sets true -median false -num " + numPerm + 
					" -plot_top_x "+ numPlot + " -rnd_seed 149 -save_rnd_lists false -set_max 5000 -set_min " + minSetSize +
					" -zip_report false -out " + save + "  -gui false "	;
					
					run.exec(command);  
				}
			}
		}
		else if ("retrieveBatchGSEA_Tool".equalsIgnoreCase(argmap.getTask())){
			String resPath=argmap.getMandatory("resPath");
			String geneSetFile=argmap.getMandatory("gmt");
			String out = argmap.getOutput();
			retrieveBatchGSEA_Tool(resPath,geneSetFile, out);
		}
		else  if("getCorrMatrix".equalsIgnoreCase(argmap.getTask())) {
			String expressionFile = argmap.getMandatory("exp");
			String pattern = argmap.getMandatory("pattern");
			String out = argmap.getOutput();
			getCorrMatrix(expressionFile,pattern,out);
		}
		
	   else{System.err.println(usage);}
	}

	
	


	




	//Calculates the gene-set enrichment of each gene set with the a list ranked by its correlation 
	//to the pattern expression.  Calculate FDR to control for the testing multiple geneSets
	private static MatrixWithHeaders continuousGSEA(String expressionF,
			String patternF, String geneSetF, String chipF, int minSetSize, String patternType, String outfile) throws IOException, ParseException {
		
		File expressionFile=new File(expressionF);
		File patternFile=new File(patternF);
		File geneSetFile=new File(geneSetF);
		File chipFile=new File(chipF);
				
		MatrixWithHeaders expression=new MatrixWithHeaders(expressionFile.getAbsolutePath());
		MatrixWithHeaders patterns=null;
		if (patternType.equalsIgnoreCase("gct"))
			{patterns=new MatrixWithHeaders(patternFile.getAbsolutePath());}
		else
			{patterns=ParseGCTFile.parseNumericCls(patternFile,expression.getColumnNames()); }		
		Map<String, Collection<String>> preGeneSets=GMTParser.ParseGMTFile(geneSetFile, minSetSize);
		Map<String, String> probGenMap=ParseGCTFile.parseChipFile(chipFile);
		
		//Step1: extract sub expression matrix with probes that corresponds to a gene symbol
		MatrixWithHeaders expMat= extractGeneSymbSubmat (expression,probGenMap);
		//Reduce geneSets to effectiveGeneSet
		Map<String, Collection<String>> geneSets=preprocessGeneSets(preGeneSets,expMat);
		
		//Step2: compute GSEA
		GeneSetEnrichment[] gseaRes= new GeneSetEnrichment[patterns.getNumberRows()];
		for (int i=0; i<patterns.getNumberRows();i++){
		 MatrixWithHeaders tmp= patterns.submatrixByRowNames(patterns.getRowName(i));
		 gseaRes[i]= new GeneSetEnrichment(tmp,expMat ,geneSets, true, false);
		}	
		
		//Step3: write GSEA FDR result matrix lincRNAs*Gene-sets
		MatrixWithHeaders fdr=null;
		MatrixWithHeaders normalizedKS=null;
		MatrixWithHeaders minMax=null;
		MatrixWithHeaders minMaxFdr=null;
		for (int i=0; i<patterns.getNumberRows();i++){
			MatrixWithHeaders ksFdr=gseaRes[i].getKSFDR();
			MatrixWithHeaders nes=gseaRes[i].getNormalizedKSEnrichments();
			MatrixWithHeaders mnmx=gseaRes[i].getNormalizedMaxMeanEnrichments();
			MatrixWithHeaders mnmxfdr=gseaRes[i].getMaxMeanFDR();
			ArrayList<String> c= new ArrayList<String>();
			c.add(patterns.getRowName(i).concat(ksFdr.getColoumnName(0)));
			MatrixWithHeaders newksFDR=new MatrixWithHeaders(ksFdr.getData(),ksFdr.getRowNames(),c);
			MatrixWithHeaders newNES=new MatrixWithHeaders(nes.getData(),nes.getRowNames(),c);
			MatrixWithHeaders newMnmx=new MatrixWithHeaders(mnmx.getData(),mnmx.getRowNames(),c);
			MatrixWithHeaders newMnmxFdr=new MatrixWithHeaders(mnmxfdr.getData(),mnmxfdr.getRowNames(),c);
			
			
			if (i==0){ fdr=newksFDR; normalizedKS=newNES;  minMax=newMnmx; minMaxFdr=newMnmxFdr;   }
			else {fdr.appendColumns(newksFDR); normalizedKS.appendColumns(newNES);
					minMax.appendColumns(newMnmx); minMaxFdr.appendColumns(newMnmxFdr);  }		
		}
		
		//step 4: write files with FDR res and NES
		fdr.writeGCT(outfile+".FDR");
		normalizedKS.writeGCT(outfile+".NES");
		minMax.writeGCT(outfile+".minMax");
		minMaxFdr.writeGCT(outfile+".minMaxFdr");
		
		return fdr;
	}


	//remove genes that are not included in the expression matrix
	private static Map<String, Collection<String>> preprocessGeneSets(
			Map<String, Collection<String>> preGeneSets, MatrixWithHeaders expression) {
		
		Map<String, Collection<String>> geneSets=new TreeMap<String, Collection<String>>();
		for (String set: preGeneSets.keySet()){
			ArrayList<String> newset= new ArrayList<String>();
			for (String symb: preGeneSets.get(set)) {
				if (expression.containsRow(symb)) newset.add(symb);	}
			if (! newset.isEmpty()) geneSets.put(set,newset);
		}
		return geneSets;
	}


	private static MatrixWithHeaders extractGeneSymbSubmat(
		MatrixWithHeaders expression, Map<String, String> probGenMap) {
		
		List<String> colNames=expression.getColumnNames();
		ArrayList<String> rowNames=new ArrayList<String>();
		ArrayList<String> probeNames=new ArrayList<String>();
		for (String prob: probGenMap.keySet()){ 
			if (prob.equalsIgnoreCase("PROBE SET ID") || !expression.containsRow(prob)) continue;
			if (probGenMap.get(prob)!=null & !prob.equalsIgnoreCase(probGenMap.get(prob))) {
				probeNames.add(prob);
				rowNames.add(probGenMap.get(prob));
				}
			}
		if (probeNames==null) return null;
		MatrixWithHeaders subMat= new MatrixWithHeaders(rowNames,colNames);
		
		for (String prob: probeNames) {
			//System.err.println(prob);
			String symb=probGenMap.get(prob);
			//if(symb==null ) {System.err.println(prob+ " failed symb");}
			double [] vals=expression.getRow(prob);
			if(symb==null || vals==null) {System.err.println(prob+ " failed ");}
			subMat.setRow(symb,vals);
		}
		return subMat;
	}
	
	
	
	private static void retrieveBatchGSEA_Tool(String resPath,
			String geneSetFile, String out) throws IOException {
		
		Map<String, Collection<String>> preGeneSets=GMTParser.ParseGMTFile(new File(geneSetFile), 1);
		int numGeneSets= preGeneSets.keySet().size();
		ArrayList<String> geneSetsNames=new ArrayList<String>(preGeneSets.keySet());
		
		boolean first=true;
		MatrixWithHeaders ES=null;
		MatrixWithHeaders NES=null;
		MatrixWithHeaders PVAL=null;
		MatrixWithHeaders FDR=null;
		MatrixWithHeaders FWER=null;
		MatrixWithHeaders Sign=null;
		MatrixWithHeaders SetSize=new MatrixWithHeaders(geneSetsNames,new ArrayList<String>(Arrays.asList("size")) ,0);
		
		//0,1,2,...
		File[] resDirs1=new File(resPath).listFiles();
		for (int d1=0; d1<resDirs1.length; d1++)
		{   //probe.cls     
			File[] resDirs2=new File(resDirs1[d1].getAbsolutePath()).listFiles();
			for(int d2=0; d2<resDirs2.length; d2++){
				File dir=resDirs2[d2];
				String patternName=dir.getName();//name of the probe/continuous pattern
				//System.err.println(patternName);
				if (dir.isFile()) continue;
				File[] resDirs3=new File(dir.getAbsolutePath()).listFiles(); //inner dir				
				for(int d3=0; d3<resDirs3.length; d3++){
					if (resDirs3[d3].isFile()) continue;
					File[] resFiles=new File(resDirs3[d3].getAbsolutePath()).listFiles(); //res files
					
					ArrayList<String> colNames=new ArrayList<String>();
					colNames.add(patternName);
					//make default value matrices 
					MatrixWithHeaders iES=new MatrixWithHeaders(geneSetsNames,colNames,-1);
					MatrixWithHeaders iNES=new MatrixWithHeaders(geneSetsNames,colNames,-100);
					MatrixWithHeaders iPVAL=new MatrixWithHeaders(geneSetsNames,colNames,-1);
					MatrixWithHeaders iFDR=new MatrixWithHeaders(geneSetsNames,colNames,-1);
					MatrixWithHeaders iFWER=new MatrixWithHeaders(geneSetsNames,colNames,-1);
					MatrixWithHeaders iSign=new MatrixWithHeaders(geneSetsNames,colNames,0);
					
					for(int d4=0; d4<resFiles.length; d4++){
						File file=resFiles[d4];
						if (file.getName().contains("gsea_report_for") & file.getName().contains("xls")){
							int sign=-1;
							if (file.getName().contains("pos")) sign =1;
							//read text data and insert data
							String line;
							BufferedReader br= new BufferedReader(new FileReader(file));
							line = br.readLine(); //header
							while((line = br.readLine()) != null) {
								line = line.trim();
								String [] lineSplit = line.split("\t");
								String geneSet=lineSplit[0];
								if (! geneSetsNames.contains(geneSet)) continue;
								if (! lineSplit[4].isEmpty())iES.set(geneSet,patternName,new Double(lineSplit[4]).doubleValue());
								if (! lineSplit[5].isEmpty())iNES.set(geneSet,patternName,new Double(lineSplit[5]).doubleValue());
								if (! lineSplit[6].isEmpty())iPVAL.set(geneSet,patternName,new Double(lineSplit[6]).doubleValue());
								if (! lineSplit[7].isEmpty())iFDR.set(geneSet,patternName,new Double(lineSplit[7]).doubleValue());
								if (! lineSplit[8].isEmpty())iFWER.set(geneSet,patternName,new Double(lineSplit[8]).doubleValue());
								iSign.set(geneSet,patternName,sign);
								SetSize.set(geneSet,"size",new Double(lineSplit[3]).doubleValue());
							}	
						}	
					}
					
					//add to final matrix 
					if (first)
					{
						first=false;
						ES=iES; NES=iNES; PVAL=iPVAL; FDR=iFDR; FWER=iFWER; Sign=iSign;
					}
					else
					{
						ES.appendColumns(iES); NES.appendColumns(iNES); PVAL.appendColumns(iPVAL); 
						FDR.appendColumns(iFDR); FWER.appendColumns(iFWER); Sign.appendColumns(iSign);
					}
				}
			}
		}
		
		//write matrices
		ES.writeGCT(out+".ES.gct");
		NES.writeGCT(out+".NES.gct"); 
		PVAL.writeGCT(out+".PVAL.gct");
		FDR.writeGCT(out+".FDR.gct");
		FWER.writeGCT(out+".FWER.gct");
		Sign.writeGCT(out+".Sign.gct");
		SetSize.writeGCT(out+".SetSize.gct");
		
	}
						
		
	private static void getCorrMatrix(String expressionF, String patternF,
			String out) throws IOException, ParseException {

	    double Threshold=0.7; 
		FileWriter writer=new FileWriter(out + ".Sparse");
		
		File expressionFile=new File(expressionF);
		File patternFile=new File(patternF);
				
		MatrixWithHeaders expression=new MatrixWithHeaders(expressionFile.getAbsolutePath());
		MatrixWithHeaders patterns=new MatrixWithHeaders(patternFile.getAbsolutePath());
		
		if (expression.getNumberColumns() != patterns.getNumberColumns()) {System.err.println("getCorrMatrix: difference in the number of conditions between the pattern and expression files");} 
		
		MatrixWithHeaders corrMat=new MatrixWithHeaders (patterns.getRowNames(),expression.getRowNames());
		
		for (int i=0; i< patterns.getNumberRows(); i++ ){
			double[] p=patterns.getRow(i);
			for (int j=0; j< expression.getNumberRows(); j++ ){
				double c=Statistics.pearsonDistance(p,expression.getRow(j));
				
				corrMat.set(patterns.getRowName(i),expression.getRowName(j), c);
				if (c>Threshold) {
					String s=new Double(c).toString();
					//writer.write(patterns.getRowName(i)+"\t"+expression.getRowName(j)+"\t");
					writer.write(i+"\t"+j+"\t");
					if (s.length()<5) writer.write (s); 
					else writer.write (s,0,5);
					writer.write("\n");
					}
			}
		}
	
		//corrMat.writeGCT(out,corrMat.getColumnNames(),corrMat.getRowNames(),2);
		//corrMat.writeSparseMat(out,corrMat.getColumnNames(),corrMat.getRowNames(),2);
		writer.close();
		patterns.writeRowNamesToFile( out + ".RowNames.Sparse");
		expression.writeRowNamesToFile( out + ".ColNames.Sparse");
		return;
	}

	
	
}
