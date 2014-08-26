package broad.pda.differentialExpression;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import broad.core.datastructures.MatrixWithHeaders;
import broad.core.error.ParseException;

public class CollapseProjections {

	public CollapseProjections(MatrixWithHeaders data, String save) throws IOException{
		//go through the matrix and for shared gene set and UP/DOWN take the best score
		data=collapse(data);
		
		//convert signs for compatability
		//UP and + is + Down and - is plus UP and - is - and DOWN and - is plus
		data=convertSigns(data);
		
		data.writeGCT(save);
	}

	private MatrixWithHeaders convertSigns(MatrixWithHeaders data) {
		//if conflicting new signs will return 0
		
		Map<String, Collection<String>> groups=new TreeMap<String, Collection<String>>();
		for(String geneSet: data.getRowNames()){
			String name=geneSet.split(":")[0];
			Collection<String> set=new TreeSet<String>();
			if(groups.containsKey(name)){set=groups.get(name);}
			set.add(geneSet);
			groups.put(name, set);
		}
		
		MatrixWithHeaders rtrn=new MatrixWithHeaders(new ArrayList(groups.keySet()), data.getColumnNames());
		
		for(String geneSet: groups.keySet()){
			Collection<String> names=groups.get(geneSet);
			MatrixWithHeaders subset=data.submatrixByRowNames(names);
			for(String experiment: subset.getColumnNames()){
				double collapsed=collapse(subset, names, experiment);
				rtrn.set(geneSet, experiment, collapsed);
			}
		}
		
		return rtrn;
	}

	private double collapse(MatrixWithHeaders subset, Collection<String> names,	String experiment) {
		double[] vals=new double[names.size()];
		
		int i=0;
		for(String name: names){
			//System.err.println(name+" "+experiment);
			double rawVal=subset.get(name, experiment);
			double factor=1;
			if(name.split(":")[1].equalsIgnoreCase("DOWN")){
				if(rawVal<0){factor=1;}
				else{factor=-1;}
			}
			else if(name.split(":")[1].equalsIgnoreCase("UP")){
				if(rawVal>0){factor=1;}
				else{factor=-1;}
			}
			vals[i++]=(rawVal*factor);
		}
		
			
		return this.computeConsensusScore(vals);
	}

	private MatrixWithHeaders collapse(MatrixWithHeaders data) {
		Map<String, Collection<String>> compressedGenes=new TreeMap<String, Collection<String>>();
		
		for(String gene: data.getRowNames()){
			String coreName=getCore(gene);
			Collection<String> set=new TreeSet<String>();
			if(compressedGenes.containsKey(coreName)){set=compressedGenes.get(coreName);}
			set.add(gene);
			compressedGenes.put(coreName, set);
		}
		
		
		MatrixWithHeaders rtrn=new MatrixWithHeaders(new ArrayList<String>(compressedGenes.keySet()), data.getColumnNames());
		
		for(String gene: compressedGenes.keySet()){
			Collection<String> allGenes=compressedGenes.get(gene);
			MatrixWithHeaders temp=data.submatrixByRowNames(allGenes);
			for(String experiment: temp.getColumnNames()){
				double[] vals=temp.getColumn(experiment);
				double consensus=computeConsensusScore(vals);
				rtrn.set(gene, experiment, consensus);
			}
		}
		
		return rtrn;
	}
	
	private double computeConsensusScore(double[] vals) {
		//Return max if positive
		double max=max(vals);
		
		//Return min if negative
		double min=min(vals);
		
		return largestAbsValue(max, min);
		
	}

	private double max(double[] vals) {
		double max=-Double.MAX_VALUE;
		for(int i=0; i<vals.length; i++){
			max=Math.max(max, vals[i]);
		}
		return max;
	}
	
	private double min(double[] vals) {
		double max=Double.MAX_VALUE;
		for(int i=0; i<vals.length; i++){
			max=Math.min(max, vals[i]);
		}
		return max;
	}

	private double largestAbsValue(double max, double min) {
		double maxAbs=Math.abs(max);
		double minAbs=Math.abs(min);
		if(maxAbs>minAbs){return max;}
		return min;
	}

	private String getCore(String gene) {
		String[] tokens=gene.split("UP");
		if(tokens.length>1){return tokens[0]+":UP";}
		else{
			tokens=gene.split("DOWN");
			return tokens[0]+":DOWN";
		}
	}

	public static void main(String[] args)throws IOException, ParseException{
		if(args.length>1){
			MatrixWithHeaders data=new MatrixWithHeaders(args[0]);
			String save=args[1];
			new CollapseProjections(data, save);
		}
		else{System.err.println(usage);}
	}
	
	static String usage=" args[0]=data \n args[1]=save";
}
