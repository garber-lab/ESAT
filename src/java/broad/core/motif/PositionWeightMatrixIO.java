package broad.core.motif;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Stack;

import broad.core.annotation.BED;
import broad.core.error.ParseException;
import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.core.sequence.SequenceRegion;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;

public class PositionWeightMatrixIO {
	Stack<PositionWeightMatrix> matrices;

	
	public PositionWeightMatrixIO() {
		matrices = new Stack<PositionWeightMatrix>();
	}

	public void load(InputStream is) throws IOException, ParseException {
		BufferedReader br = new BufferedReader(new InputStreamReader(is));
		String line = null;
		List<String[]> rawMatrix = null;
		while ((line = br.readLine()) != null) {
			if(line.startsWith("#") || line.length() == 0) {
				System.err.println("Comment or empty line: " + line);
				continue;
			}
			//System.err.println("line: " + line);
			if(line.startsWith(">")) {
				
				if(rawMatrix != null) {
					PositionWeightMatrix pwm = matrices.peek();
					pwm.setMatrixFromRawData(rawMatrix);
				}
				
				PositionWeightMatrix pwm = new PositionWeightMatrix(line.substring(1));
				matrices.push(pwm);
				rawMatrix = new ArrayList<String[]>();
			} else {
				//PositionWeightColumn column = new PositionWeightColumn(line.split("\t"));
				//System.err.println("addeing " + column);
				//matrices.peek().add(column);
				//System.err.println("line:  " + line + " matrix " + rawMatrix);
				rawMatrix.add(line.split("\\s+"));
			}
		}
		
		if(rawMatrix != null) {
			matrices.peek().setMatrixFromRawData(rawMatrix);
		}
		
	}
	
	public void write(OutputStream os, NumberFormat formatter, List<PositionWeightMatrix> pwms) throws IOException{
		BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(os));
		Iterator<PositionWeightMatrix> it = pwms.iterator();
		while(it.hasNext()) {
			it.next().write(bw, formatter);
		}
		bw.flush();
	}
	
	public void writeUNIPROBE(OutputStream os, NumberFormat formatter, List<PositionWeightMatrix> pwms) throws IOException{
		BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(os));
		Iterator<PositionWeightMatrix> it = pwms.iterator();
		while(it.hasNext()) {
			it.next().writeUNIPROBE(bw, formatter);
		}
		bw.flush();
	}
	
	public void write(OutputStream os, NumberFormat formatter) throws IOException{
		write(os, formatter, matrices);
	}
	
	public List<PositionWeightMatrix> getMatrices() {
		return matrices;
	}
	
	public static final String USAGE = "Usage: PositionWeightMatrixIO TASK=<task_num> <task_args> -pwm <PWM file in FASTA fromat>\n" +
	"\tTasks:" +
	"\n\t\t1. Report matches of a PWMs to a sequence -in <Sequence file or standard input in FASTA format> -pwm <File containing PWMs> -minScore <Minimum Score to report>" +
	"\n\t\t\t-start <start matching at this point on the corresponding sequence include as many starts as sequenceIds> -end <Stop matching at this point on the corresponding sequence include as many ends as sequenceIds> " +
	"\n\t\t\t[-seqId <Sequence Ids in the fasta file to search on specify as many as targets you want to search on, if non is specified all start/ends will be assume to relate to first sequence in file>] "+
	"\n\t\tic Compute information content of all motifs on the file -pwm <File containing PWMs>" +
	"\n";
	
	public static void main(String [] args) throws Exception{
		ArgumentMap argMap = CLUtil.getParameters(args, USAGE);
		double [] neutral = {0.3,0.2,0.2,0.3};
		if("1".equals(argMap.getTask())) {
			FastaSequenceIO fsio = new FastaSequenceIO();
			//List<Sequence> seqs = fsio.loadAll(argMap.getInputStream());
			List<Integer> starts = argMap.containsKey("start") ? argMap.getIntegers("start") : null;
			List<Integer> ends = argMap.containsKey("end") ? argMap.getIntegers("end") : null;
			if(starts != null && starts.size() != ends.size()) { throw new IllegalArgumentException("Number of end and start positions specified should be the same");}
			List<String> seqIds = argMap.getAll("seqId");
			float minScore = argMap.containsKey("minScore") ? argMap.getFloat("minScore"): -1000;
			if(seqIds == null || seqIds.size() == 0) {
				BufferedReader br = new BufferedReader(new InputStreamReader(argMap.getInputStream()));
				try {
					String seqId = br.readLine();
					seqId = seqId.replaceFirst(">", "");
					seqIds = new ArrayList<String>(starts.size());
					for(int i = 0; i < ends.size(); i++) {
						seqIds.add(seqId);
					}
				} finally {
					if(br != null) { br.close();}
				}
			
			}
			InputStream is = argMap.getInputStream();
			List<SequenceRegion> regions = new ArrayList<SequenceRegion>();
			if(starts!= null && starts.size() > 0) { 
				for(int i = 0; i < ends.size(); i++) {
					SequenceRegion region = new SequenceRegion(seqIds.get(i));
					region.setStart(starts.get(i));
					region.setEnd(ends.get(i));
					regions.add(region);
				}

				fsio.extractRegions(regions, false, is);
			} else {
				List<Sequence> seqs = fsio.extractRecords(seqIds, is);
				Iterator<Sequence> seqIt = seqs.iterator();
				while(seqIt.hasNext()) {
					Sequence seq = seqIt.next();
					SequenceRegion region = new SequenceRegion(seq.getId());
					region.setStart(0);
					region.setEnd(seq.getLength() - 1);
					regions.add(region);
					region.setSequence(seq);
				}
			}
			is.close();
			String pwmFile = argMap.getMandatory("pwm");
			PositionWeightMatrixIO pwmIO = new PositionWeightMatrixIO();
			FileInputStream fis = new FileInputStream(pwmFile);
			try {
				pwmIO.load(fis);
			} finally {	
				if(fis != null) {fis.close();}
			}
			

			
			Map<PositionWeightMatrix, List<BED>> matches = new HashMap<PositionWeightMatrix, List<BED>>();
			Iterator<SequenceRegion> it = regions.iterator();
			while(it.hasNext()) {
				SequenceRegion region = it.next();
				Iterator<PositionWeightMatrix> pwmIt = pwmIO.matrices.iterator();
				while(pwmIt.hasNext()) {
					PositionWeightMatrix pwm = pwmIt.next();
					List<BED> pwmMatches = matches.get(pwm);
					if(pwmMatches == null) {
						pwmMatches = new ArrayList<BED>();
						matches.put(pwm, pwmMatches);
					}
					
					pwmMatches.addAll(pwm.match(region, neutral, minScore));
				}	
			}
			
			Iterator<PositionWeightMatrix> pwmIt = pwmIO.getMatrices().iterator();
			BufferedWriter bw = argMap.getOutputWriter();
			while(pwmIt.hasNext()) {
				PositionWeightMatrix pwm = pwmIt.next();
				bw.write("track name=\""+pwm.getName() + "\" visibility=2 \n");
				Iterator<BED> pwmMatchIt = matches.get(pwm).iterator();
				while(pwmMatchIt.hasNext()) {
					bw.write(pwmMatchIt.next().toString(false));
					bw.newLine();
				}
			}
			bw.close();
		} else if ("ic".equalsIgnoreCase(argMap.getTask())) {
			String pwmFile = argMap.getMandatory("pwm");
			PositionWeightMatrixIO pwmIO = new PositionWeightMatrixIO();
			FileInputStream fis = new FileInputStream(pwmFile);
			try {
				pwmIO.load(fis);
			} finally {	
				if(fis != null) {fis.close();}
			}
			List<PositionWeightMatrix> pwms = pwmIO.getMatrices();
			BufferedWriter bw = argMap.getOutputWriter();
			for(PositionWeightMatrix pwm : pwms) {
				bw.write(pwm.getName()+"\t"+pwm.ic());
				bw.newLine();
			}
		}
	}

	public void addPseudoCounts() {
		for(PositionWeightMatrix pwm : matrices) {
			pwm.addPseudoCounts();
		}
		
	}
}
