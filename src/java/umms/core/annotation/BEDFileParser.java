package umms.core.annotation;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import org.apache.log4j.Logger;



public class BEDFileParser {
	static Logger logger = Logger.getLogger(BEDFileParser.class.getName());
	/**
	 * Read genes from the bed file and get genes by chromosome
	 * @param file Name of bed file
	 * @return Map of chromosome name to set of genes on chromosome
	 * @throws IOException
	 */
	public static Map<String, Collection<Gene>> loadDataByChr(String file) throws IOException{
		return loadDataByChr(new File(file));
	}
	
	public static Map<String, Collection<Gene>> loadDataByChr(File file) throws IOException{
		logger.info("Loading genes from file " + file.getName() + "...");
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(file)));
	
		Map<String, Collection<Gene>> rtrn=new TreeMap<String, Collection<Gene>>();
		String nextLine;
		int i=0;
		while ((nextLine = reader.readLine()) != null ) {
	
			if(looksLikeData(nextLine) ){
				//logger.info(nextLine);
				Gene gene = new Gene(nextLine, false);
				//System.err.println("Gene: " + gene.toBED());
				
				Collection<Gene> data=new TreeSet<Gene>();
				if(rtrn.containsKey(gene.getChr())){
					data.addAll(rtrn.get(gene.getChr()));
				}
				data.add(gene);
				rtrn.put(gene.getChr(), data);
			}
			i++;
			if(i%10000==0){logger.info("Loaded " + i + " genes.");}
		}

		reader.close();
		return rtrn;
	
	}
	
	private static boolean looksLikeData(String nextLine) {
		return nextLine.trim().length() > 0 && ! nextLine.startsWith("#") && !nextLine.startsWith("track") && !nextLine.startsWith("browser");
	}

    

}
