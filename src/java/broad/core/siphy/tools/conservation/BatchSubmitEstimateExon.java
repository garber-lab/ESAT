package broad.core.siphy.tools.conservation;

import java.io.*;
import java.util.*;

public class BatchSubmitEstimateExon {
	
	private static final String quote="\"";

	private static void writeError(InputStream errorStream) throws IOException {
		BufferedReader reader=	new BufferedReader(new InputStreamReader(errorStream));
		String nextLine;
		while ((nextLine = reader.readLine()) != null) {
			System.err.println(nextLine);
		}
		System.err.println();
	}
	
	public static void main(String[] args)throws IOException, InterruptedException{
		if(args.length>6){
		File[] files=new File(args[0]).listFiles();
		File modelFile=new File(args[1]);
		String alnDir=args[2];
		String alnFormat=args[3];
		String saveDir=args[4];
		String queue=args[5];
		String script=args[6];
		
		Runtime run=java.lang.Runtime.getRuntime();
		
		for(int i=0; i<files.length; i++){
			File file=files[i];
			String save=saveDir+"/"+files[i].getName();
			String command="bsub -q "+queue+" -o "+save+".bsub"+" java -jar -Xmx3000m "+script+" "+file.getAbsolutePath()+" "+modelFile.getAbsolutePath()+" "+alnDir+" "+alnFormat+" "+save;
			Process p=run.exec(command);
			p.waitFor();
			System.err.println(files[i]+" "+p.exitValue());
			writeError(p.getErrorStream());
			System.err.println(command);
		}
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=files \n args[1]=model file \n args[2]=alnDir \n args[3]=align format \n args[4]=save directory \n  args[5]=queue \n args[6]=script";
	
}
