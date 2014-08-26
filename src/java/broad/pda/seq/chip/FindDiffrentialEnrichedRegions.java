package broad.pda.seq.chip;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;

import nextgen.core.annotation.Annotation;

import org.apache.commons.math3.distribution.NormalDistribution;

import broad.pda.annotation.BEDFileParser;
import broad.pda.datastructures.Alignments;
import broad.pda.seq.segmentation.AlignmentDataModel;
import broad.pda.seq.segmentation.ContinuousDataAlignmentModel;
import broad.pda.seq.segmentation.GenericAlignmentDataModel;

public class FindDiffrentialEnrichedRegions {

	double alpha=.05;
	double minDiff=5.0;
	
	//TODO Do this with Skellam because normals are too liberal
	
	public FindDiffrentialEnrichedRegions(Collection<Annotation> regions, AlignmentDataModel data1, AlignmentDataModel data2, String save)throws IOException{
		ContinuousDataAlignmentModel model1=new ContinuousDataAlignmentModel(data1);
		ContinuousDataAlignmentModel model2=new ContinuousDataAlignmentModel(data2);
		
		Map<Annotation, double[]> differences=computeDiff(model1, model2, regions);
		
		write(save, differences);
	}

	private Map<Annotation, double[]> computeDiff(ContinuousDataAlignmentModel model1, ContinuousDataAlignmentModel model2, Collection<Annotation> regions) throws IOException {
		Map<Annotation, double[]> rtrn=new TreeMap<Annotation, double[]>();
		
		for(Annotation align: regions){
			if(model1.scoreSegment(align)[0]<alpha || model2.scoreSegment(align)[0]<alpha){
				double avg1=model1.getLambda(align.getChr());
				double avg2=model2.getLambda(align.getChr());
				
				double avg=avg1-avg2;
				double variance=avg1+avg2;
				
				NormalDistribution dist=new NormalDistribution(avg, 1);
				
				double val1=model1.count(align);
				double val2=model2.count(align);
				double val=val1-val2;
				
				double cdf=dist.cumulativeProbability(val);
				double p=Math.min(1-cdf, cdf)*2;
				//System.err.println(avg+" "+variance+" "+val+" "+cdf);
				double[] array={p, val};
				if(array[0]<alpha){rtrn.put(align, array);}
			}
		}
		return rtrn;
	}

	private void write(String save, Map<Annotation, double[]> differences) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(Annotation data: differences.keySet()){
			double[] vals=differences.get(data);
			writer.write(data+"\t"+vals[1]+"\n");
		}
		
		writer.close();
	}

	public static void main(String[] args)throws IOException {
		if(args.length>4){
			Collection<Annotation> regions=BEDFileParser.loadAlignmentData(new File(args[0]));
			AlignmentDataModel data1=new GenericAlignmentDataModel(args[1], args[3]);
			AlignmentDataModel data2=new GenericAlignmentDataModel(args[2], args[3]);
			String sizes=args[4];
			new FindDiffrentialEnrichedRegions(regions, data1, data2, sizes);
		}
		else{
			System.err.println(usage);
		}
	}
	
	static String usage=" args[0]=Regions \n args[1]=Alignment file 1 \n args[2]=Alignment file 2 \n args[3]=sizes \n args[4]=save";
}
