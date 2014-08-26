package broad.pda.seq.segmentation;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map;

import broad.core.sequence.SequenceUtils;

public class PrecomputeDataAlignmentStats {

	public PrecomputeDataAlignmentStats(File in, String sizes, String save) throws IOException{
		ContinuousDataAlignmentModel cda=SequenceUtils.getDataModel(in.getAbsolutePath(), sizes, false);
		//just get the number of reads
		Map<String, Double> map=cda.getAlignmentDataModelStats().getNumberOfReads();
		FileWriter writer=new FileWriter(save);
		
		for(String chr: map.keySet()){
			writer.write(chr+"\t"+map.get(chr)+"\n");
		}
		writer.close();
	}
	
	public static void main(String[] args)throws IOException{
		if(args.length>1){
			File bam=new File(args[0]);
			String sizes=args[1];
			String save=bam+".cda";
			new PrecomputeDataAlignmentStats(bam, sizes, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=bam \n args[1]=sizes";
	
}
