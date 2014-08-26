package broad.core.motif;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import broad.core.sequence.FastaSequenceIO;
import broad.core.sequence.Sequence;
import broad.core.sequence.SequenceRegion;
import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;

public class SequenceMotifCluster extends SequenceMotif {
	public static String USAGE = "Usage: SequenceMotifCluster TASK=<task_num> <task_args>\n" +
	"\tTasks:\n" +
	"\t\t1. Get list of motif matches: IN=<fasta sequence> MOTIF=<Core Motif e.g. YRNRY> SUSCCESSIVE=<number of consecutive core motifs> MINLENGTH=<minimum total length to report> MAXSPACING=<max spacing between motifs>\n";
	
	int maxDistance;
	Pattern degeneratedMotif;

	public SequenceMotifCluster(String corePattern, int minsuccessiveCores, int maxPatternDistance)
			throws SearchException {
		super(corePattern, minsuccessiveCores);
		maxDistance = maxPatternDistance;
		String spacing = ".{0,"+maxPatternDistance +"}";
		
		Pattern [] coreMotifs = super.getCoreMotifs();		
		
		StringBuffer coreMotif = new StringBuffer("(");
		for (int j = 0; j < coreMotifs.length; j++) {
			coreMotif.append(coreMotifs[j].pattern());
			if(j < coreMotifs.length - 1) {
				coreMotif.append("|");
			}
		}
		coreMotif.append(")+");
		
		StringBuffer buf = new StringBuffer("(");
		buf.append(coreMotif)
			.append("(")
			.append(spacing)
			.append(")")
			.append("){")
			.append((minsuccessiveCores - 1))
			.append(",}")
			.append(coreMotif);
		
		degeneratedMotif = Pattern.compile(buf.toString());
	}
	
	public List<SequenceRegion> match(Sequence seq) { 
		return match(seq, true);
	}

	public List<SequenceRegion> match(Sequence seq, boolean includeRegionSequence) { 
		//System.out.println("Degenerated Motif " + degeneratedMotif);
		ArrayList<SequenceRegion> matches = new ArrayList<SequenceRegion>();
		String seqBases = seq.getSequenceBases();
		if(includeRegionSequence) {
			seqBases = seqBases.toUpperCase();
		}
		//System.out.println(">"+seq.getId()+"\n"+seq.getSequenceBases());
		Matcher matcher = degeneratedMotif.matcher(seqBases);
		Matcher undegenerateMatcher = null;
		while(matcher.find()) {
			int start = matcher.start();
			int end   = matcher.end();
			SequenceRegion reg = seq.getRegion(start, end); //include last base.
			
			undegenerateMatcher = super.getMotif().matcher(reg.getSequenceBases());
			if(undegenerateMatcher.find()) {
				matches.add(reg);				
			} else {
				//System.out.println("Sequence " + start + "-" + end + "  --> " + reg.getSequenceBases() +" matched degenerate but did not contained undegenerated ");
			}
		}
		//System.out.println("Finished getting candidate zDNA");
		return matches;
	}
	/**
	 * @param args
	 * @throws SearchException 
	 * @throws NumberFormatException 
	 * @throws IOException 
	 */
	public static void main(String[] args) throws NumberFormatException, SearchException, IOException {
		ArgumentMap argMap = CLUtil.getParameters(args, USAGE);
		if ("1".equals(argMap.getTask())) {
			int minLength = argMap.getInteger("MINLENGTH");
			SequenceMotif motif = new SequenceMotifCluster(argMap.getMandatory("MOTIF"), argMap.getInteger("SUCCESSIVE"), argMap.getInteger("MAXSPACING"));
			FastaSequenceIO fsio = new FastaSequenceIO(argMap.getInput());
			Iterator<Sequence> seqIt = fsio.loadAll().iterator();
			Sequence seq = null;
			while(seqIt.hasNext()) {
				seq = seqIt.next();
				List<SequenceRegion> matches = motif.match(seq);
				//System.out.println("Sequence " + seq.getId());
				if(matches.size() == 0) {
					System.out.println("No matches");
				} else {
					Iterator<SequenceRegion> matchIt = matches.iterator();
					SequenceRegion seqReg = null;
					while(matchIt.hasNext()) {
						seqReg = matchIt.next();
						if(seqReg.getLength() >= minLength) {
							System.out.println((seqReg.getStart() + 1) + "\t" + (seqReg.getEnd())+ "\t" + seqReg.getSequenceBases());
						}
					}
				}
				
			}
		} else {
			System.err.println(USAGE);
		}

	}
}
