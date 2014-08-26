package broad.core.sequence;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

import org.broad.igv.Globals;

import broad.pda.annotation.BEDFileParser;
import broad.pda.seq.alignment.StoredAlignmentStats;
import broad.pda.seq.segmentation.AlignmentDataModelStats;
import broad.pda.seq.segmentation.ContinuousDataAlignmentModel;
import broad.pda.seq.segmentation.GenericAlignmentDataModel;

public class SequenceUtils {

	public static Sequence getChrSequence(String seqDir, String chr) throws IOException{
		String chrSeqFile=seqDir+"/"+chr.replaceAll("chr", "")+"/"+chr+".fa";
		if(new File(chrSeqFile).exists()){
			FastaSequenceIO fsio = new FastaSequenceIO(chrSeqFile);
			Sequence chrSeq=fsio.loadAll().get(0);
			return chrSeq;
		}
		return null;
	}
	
	public static String getChrSizesFile(String seqDir){
		return seqDir+"/sizes";
	}
	
	public static ContinuousDataAlignmentModel getDataModel(String bamFile, String sizeFile, String chr, boolean removeDuplicates) throws IOException{
		if(bamFile==null){return null;}
		Globals.setHeadless(true);
		
		AlignmentDataModelStats data;
		
		//check for .cda
		if(new File(bamFile+".cda").exists()){
			Map<String, Double> numberReads=parseCDA(bamFile+".cda");
			data=new AlignmentDataModelStats(new GenericAlignmentDataModel(bamFile, sizeFile, 0, removeDuplicates), null, numberReads);
		}
		else{
			data=new AlignmentDataModelStats(new GenericAlignmentDataModel(bamFile, sizeFile, 0, removeDuplicates), null, chr);
		}
		ContinuousDataAlignmentModel model=new ContinuousDataAlignmentModel(data);
		return model;
	}
	
	public static ContinuousDataAlignmentModel getDataModel(String bamFile, String sizeFile, boolean removeDuplicates) throws IOException{
		/*
		if(bamFile==null){return null;}
		Globals.setHeadless(true);
		
		AlignmentDataModelStats data;
		
		//check for .cda
		if(new File(bamFile+".cda").exists()){
			Map<String, Double> numberReads=parseCDA(bamFile+".cda");
			data=new AlignmentDataModelStats(new GenericAlignmentDataModel(bamFile, sizeFile, 0, removeDuplicates), null, numberReads);
		}
		else{
			data=new AlignmentDataModelStats(new GenericAlignmentDataModel(bamFile, sizeFile, 0, removeDuplicates));
		}
				
		ContinuousDataAlignmentModel model=new ContinuousDataAlignmentModel(data);
		return model;
		*/
		return getDataModel(bamFile, sizeFile, 0, removeDuplicates, false, null, false);
	}
	
	public  static ContinuousDataAlignmentModel getDataModel(String bamFile, String sizeFile, int minMappingQuality, 
			  boolean removeDuplicateFlags, boolean weighReadCounts, String strand, boolean loadPairsAsFragments) throws IOException {

		if (bamFile == null) return null;
		
		AlignmentDataModelStats dataModelStats = null;
		GenericAlignmentDataModel dataModel = new GenericAlignmentDataModel(bamFile, sizeFile, false, minMappingQuality, removeDuplicateFlags, weighReadCounts, strand, loadPairsAsFragments);
		
		//System.out.println(bamFile + ".cda");
		//System.out.println(new File (bamFile + ".cda").exists());
		
		if (new File(bamFile + ".cda").exists()) {
			Map<String, Double> numberReads = parseCDA(bamFile+".cda");
			dataModelStats = new AlignmentDataModelStats(dataModel, null, numberReads);
		} else {
			dataModelStats = new AlignmentDataModelStats(dataModel);
		}
		
		return new ContinuousDataAlignmentModel(dataModelStats);
	}
	

	private static Map<String, Double> parseCDA(String string) throws IOException {
		Map<String, Double> rtrn=new TreeMap<String, Double>();
		
		Collection<String> lines=BEDFileParser.loadList(string);
		
		for(String line: lines){
			String[] tokens=line.split("\t");
			rtrn.put(tokens[0], new Double(tokens[1]));
		}
		
		return rtrn;
	}
	
	
}
