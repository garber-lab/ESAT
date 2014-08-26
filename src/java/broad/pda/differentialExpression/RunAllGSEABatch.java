package broad.pda.differentialExpression;

import java.io.File;
import java.io.IOException;

public class RunAllGSEABatch {

	public static void main(String[] args)throws IOException{
		if(args.length>7){
			File[] gctFiles=new File(args[0]).listFiles();
			File experimentInfo=new File(args[1]);
			File chipFile=new File(args[2]);
			String saveDir=args[3];
			boolean preprocess=new Boolean(args[4]);
			File gmtFile=new File(args[5]);
			
			String script=args[6];
			String queue=args[7];
			
			Runtime run=Runtime.getRuntime();
			for(int i=0; i<gctFiles.length; i++){
				String save=saveDir+"/"+gctFiles[i].getName()+".gsea";
				String command="bsub -q "+queue+" -o bsub.junk "+" java -jar -Xmx3000m "+script+" "+gctFiles[i].getAbsolutePath()+" "+experimentInfo.getAbsolutePath()+" "+chipFile.getAbsolutePath()+" "+save+" "+preprocess+" "+gmtFile.getAbsolutePath();
				System.err.println(command);
				run.exec(command);
			}
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=gct files \n args[1]=experiment info \n args[2]=chip file \n args[3]=save dir \n args[4]=preprocess \n args[5]=gmt file \n args[6]=script \n args[7]=queue";
	
}
