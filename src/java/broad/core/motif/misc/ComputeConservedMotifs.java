package broad.core.motif.misc;

import java.io.*;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import broad.core.annotation.BED;
import broad.core.motif.PositionWeightMatrix;
import broad.core.motif.PositionWeightMatrixIO;
import broad.core.siphy.TreeScaler;
import broad.pda.datastructures.Alignments;

public class ComputeConservedMotifs {

	
	public ComputeConservedMotifs(Alignments region, String pwmFile, String saveDir, String alignDir, String alignFormat, File modelFile)throws Exception{
		String alignFile=(alignDir+"/"+region.getChr()+".maf");
		TreeScaler scaler=new TreeScaler(region, modelFile, alignFile, alignFormat);
		List<BED>[] alignmentMatches=scoreMotifs(region, scaler, alignFormat, alignFile, pwmFile, true);
		String[] pwmNames=getMotifNames(pwmFile);
		writeWiggle(saveDir, alignmentMatches, pwmNames, region);
	}
	
	private void writeWiggle(String saveDir, List<BED>[] alignmentMatches, String[] pwmNames, Alignments align)throws IOException{
		for(int i = 0; i < alignmentMatches.length; i++) {
			String pwm = pwmNames[i];
			FileWriter writer=new FileWriter(saveDir+"/"+pwm.trim()+".wig");
			writer.write("browser position "+align.toUCSC()+"\n");
			writer.write("track type=wiggle_0 name="+pwm.trim()+" visibility=full autoScale=off \n");
			List<BED> list=alignmentMatches[i];
			for(BED bed: list) {
				writer.write(bed+"\n");
			}
			writer.close();
		}
		
	}
	
	public static String[] getMotifNames(String pwmFile)throws Exception{
		PositionWeightMatrixIO pwmIO = new PositionWeightMatrixIO();
		FileInputStream fis = new FileInputStream(pwmFile);
		pwmIO.load(fis);
		fis.close();
		List<PositionWeightMatrix> pwms = pwmIO.getMatrices();
		String[] rtrn=new String[pwms.size()];
		for(int i=0; i<rtrn.length; i++){
			rtrn[i]=pwms.get(i).getName().trim();
		}
		return rtrn;
	}
	
	public static List<BED>[] scoreMotifs(Alignments region, TreeScaler scaler, String alnFormat, String alignFile, String pwmFile, boolean reverseCompliment)throws Exception{
		try{
		PositionWeightMatrixIO pwmIO = new PositionWeightMatrixIO();
		FileInputStream fis = new FileInputStream(pwmFile);
		pwmIO.load(fis);
		fis.close();
		List<PositionWeightMatrix> pwms = pwmIO.getMatrices();
		List<PositionWeightMatrix> rPWMS = new ArrayList<PositionWeightMatrix>(pwms.size());
		Iterator<PositionWeightMatrix> it = pwms.iterator();
		
		List[] rtrn=new List[pwms.size()];
		int counter=0;
		while(it.hasNext()) {
			rtrn[counter++]=scaler.slidePWM(it.next(), -1000, 0); //Min seed score, num perm
		}
		return rtrn;
		}catch(NullPointerException ex){System.err.println("caught null pointer exception");}
		
		return null;
	}
	
	public static void main(String[] args)throws Exception{
		if(args.length>3){
			Alignments region=new Alignments(args[0]);
			String pwmFile=args[1];
			String saveDir=args[2];
			String alignDir=args[3];
			String alignFormat=args[4];
			File modelFile=new File(args[5]);
			new ComputeConservedMotifs(region, pwmFile, saveDir, alignDir, alignFormat, modelFile);
		}
		else{
			System.err.println("USAGE:\n args[0]=alignment Region\n args[1]=PWM file\n args[2]=save directory\n args[3]=alignment directory\n args[4]=align format \n args[5]=model file");
		}
	}
	
}
