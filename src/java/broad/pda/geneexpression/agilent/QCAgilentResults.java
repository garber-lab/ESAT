package broad.pda.geneexpression.agilent;

import java.util.*;
import java.io.*;

public class QCAgilentResults {

	public QCAgilentResults(File[] files, String save) throws IOException{
		FileWriter writer=new FileWriter(save);
		for(int i=0; i<files.length; i++){
			AgilentArrayStats stats=AgilentUtils.parseAgilentExperimentGlobalStats(files[i]);
			String passesQC=validate(stats);
			writer.write(files[i].getName()+"\t"+passesQC+"\n");
		}
		writer.close();
	}

	
	//This is hardcoded based on Julie's params
	//Return type is either Pass, Check, or Fail
	public static String validate(AgilentArrayStats stats) {
		
		boolean failed=false;
		boolean check=false;
		
		//Test 1: AnyColorPrcntFeatNonUnifOL
		if(stats.AnyColorPrcntFeatNonUnifOL>2){failed=true;}
		else if(stats.AnyColorPrcntFeatNonUnifOL>1){check=true;}
		
		//Test 2:
		if(stats.eQCOneColorSpikeDetectionLimit>2.5 || stats.eQCOneColorSpikeDetectionLimit<.001){failed=true;}
		else if(stats.eQCOneColorSpikeDetectionLimit>2 || stats.eQCOneColorSpikeDetectionLimit<.01){check=true;}
		
		//Test 3:
		if(stats.Metric_absGE1E1aSlope>1.5 || stats.Metric_absGE1E1aSlope<.75){failed=true;}
		else if(stats.Metric_absGE1E1aSlope>1.2 || stats.Metric_absGE1E1aSlope<.9){check=true;}
		
		//Test 4:
		if(stats.Metric_gE1aMedCVProcSignal>11){failed=true;}
		else if(stats.Metric_gE1aMedCVProcSignal>=8){check=true;}
		
		//Test 5:
		if(stats.gNegCtrlAveBGSubSig>7 || stats.gNegCtrlAveBGSubSig<-12){failed=true;}
		else if(stats.gNegCtrlAveBGSubSig>5 || stats.gNegCtrlAveBGSubSig<-10){check=true;}
		
		//Test 6:
		if(stats.Metric_gNegCtrlAveNetSig>45){failed=true;}
		else if(stats.Metric_gNegCtrlAveNetSig>40){check=true;}
		
		//Test 7
		if(stats.gNegCtrlSDevBGSubSig>12.5){failed=true;}
		else if(stats.gNegCtrlSDevBGSubSig>10){check=true;}
		
		//Test 8
		if(stats.Metric_gNonCntrlMedCVProcSignal>15){failed=true;}
		else if(stats.Metric_gNonCntrlMedCVProcSignal>10){check=true;}
		
		if(stats.Metric_gSpatialDetrendRMSFilteredMinusFit>18){failed=true;}
		else if(stats.Metric_gSpatialDetrendRMSFilteredMinusFit>15){check=true;}
		
		
		
		if(failed){return "failed";}
		if(check){return "check";}
		return "pass";
	}
	
	public static void main(String[] args)throws IOException{
		File[] files=new File(args[0]).listFiles();
		String save=args[1];
		new QCAgilentResults(files, save);
	}
	
}
