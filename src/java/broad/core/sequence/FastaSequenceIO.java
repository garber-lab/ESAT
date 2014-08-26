package broad.core.sequence;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Stack;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.log4j.Logger;

import broad.pda.rnai.ExtractSequence;


/**
 * Suppose to be a very generic and smart Fasta manipulation
 * Class. It is meant to handle both small and very large, multi
 * and single fasta files.
 * 
 * As ussual its initial version will be soooo much less ambitious.
 * 
 * @author MGarber
 *
 */
public class FastaSequenceIO {
	public static final int LINE_LENGTH = 60;
	private static Logger logger = Logger.getLogger(FastaSequenceIO.class.getName());
	
	public static void main(String [] args) throws IOException {
		String fileName = args[0];
		String seqId = args[1];
		String regions = args[2];
		String outFile = args[3];
		
		//System.out.println("looking for seq "+seqId+" in file "+fileName);
		
		String [] regStrings = regions.split(",");
		ArrayList<SequenceRegion> seqs = new ArrayList<SequenceRegion>(regStrings.length);
		for(int i = 0; i < regStrings.length; i++) {
			String [] startEnd = regStrings[i].split("\\.\\.");
			//System.out.println(regStrings[i]);
			SequenceRegion seq = new SequenceRegion(seqId);
			seq.setRegionStart(Integer.parseInt(startEnd[0]));
			seq.setRegionEnd(Integer.parseInt(startEnd[1]));
			seqs.add(seq);
		}
		
		FastaSequenceIO fsIOIn = new FastaSequenceIO(fileName);
		fsIOIn.extractRegions(seqs);
		FastaSequenceIO fsIOOut = new FastaSequenceIO(outFile);
		fsIOOut.write(seqs);
	}
	File file;
	public FastaSequenceIO(String fileName) {
		file = new File(fileName);
	}
	
	public FastaSequenceIO(File file) {
		this.file = file; 
	}
	
	public FastaSequenceIO() {
		super();
	}

	public List<Sequence> loadAll() throws IOException {
		logger.info("Loading fasta sequences from file " + file.getName() + "...");
		InputStream is = new FileInputStream(file);
		List<Sequence> seqs = loadAll(is);
		is.close();
		logger.info("Loaded " + seqs.size() + " sequences.");
		return seqs;
	}
	
	public static Collection<Sequence> loadSequences(File fasta) throws IOException {
		FastaSequenceIO fsio = new FastaSequenceIO(fasta);
		Collection<Sequence> rtrn = new ArrayList<Sequence>();
		rtrn.addAll(fsio.loadAll());
		return rtrn;
	}
	
	public static Map<String, Sequence> loadSequencesByNameFromDirectory(File sequenceDirectory) throws Exception {
		return loadSequencesByNameFromDirectory(sequenceDirectory, false);
	}
	
	public static Map<String, Sequence> loadSequencesByNameFromDirectory(File sequenceDirectory, boolean nonrandomOnly) throws Exception {
		if(!sequenceDirectory.isDirectory()) {
			throw new Exception("sequence directory " + sequenceDirectory +" is not a directory");
		}
		String seqDirPath = sequenceDirectory.getAbsolutePath();
		seqDirPath = seqDirPath.lastIndexOf("/") == seqDirPath.length() -1 
			? seqDirPath.substring(0, seqDirPath.length() -1)
			: seqDirPath;
			
		Map<String, Sequence> rtrn = new TreeMap<String, Sequence>();
		String [] seqFiles = sequenceDirectory.list();		
		for(int i = 0; i < seqFiles.length; i++) {
			File seqfile = new File(sequenceDirectory.getAbsoluteFile() + "/" + seqFiles[i]);
			if(! seqfile.getAbsolutePath().endsWith(".fa")) {
				continue;
			}
			if(nonrandomOnly && seqfile.getAbsolutePath().contains("random")) {
				continue;
			}
			
			Sequence chr = ExtractSequence.getFirstSequence(seqDirPath + "/" + seqFiles[i]);
			rtrn.put(chr.getId(), chr);
		}
		return rtrn;
	}
	
	
	public static Collection<String> getSequenceNames(String fasta) throws IOException {
		FileReader r = new FileReader(fasta);
		BufferedReader b = new BufferedReader(r);
		Collection<String> rtrn = new ArrayList<String>();
		while(b.ready()) {
			String line = b.readLine();
			if(line.contains(">")) {
				rtrn.add(line.replace(">", ""));
			}
		}
		r.close();
		b.close();
		return rtrn;
	}
	
	/**
	 * @param genomeFasta Fasta file of chromosomes
	 * @return Map of chromosome name to sequence
	 * @throws IOException
	 */
	public static Map<String, Sequence> getChrSequencesFromFasta(String genomeFasta) throws IOException {
		
		logger.info("Reading chromosome sequences from file " + genomeFasta + "...");
		FastaSequenceIO fsio = new FastaSequenceIO(genomeFasta);
		List<broad.core.sequence.Sequence> chrs = fsio.loadAll();
		Map<String, broad.core.sequence.Sequence> rtrn = new TreeMap<String, Sequence>();
		for(broad.core.sequence.Sequence chr : chrs) {
			rtrn.put(chr.getId(), chr);
		}
		logger.info("Loaded " + rtrn.size() + " chromosomes.");
		return rtrn;
				
	}

	public List<Sequence> loadAll(InputStream is) throws IOException {
		BufferedReader br = new BufferedReader(new InputStreamReader(is));
		String line = null;
		String currentSeqId = null;
		List<Sequence> seqs = new ArrayList<Sequence>();
		while((line = br.readLine()) != null) {
			if (line.startsWith(">")) {
				currentSeqId = line.substring(1);
				Sequence seq = new Sequence(currentSeqId);
				seqs.add(seq);
				continue;
			}
			
			seqs.get(seqs.size() - 1).append(line);
		}		
		
		br.close();		
		return seqs;		
	}
	public void extractRecordsIntoSequenceList(List<? extends Sequence> sequenceList) throws Exception {
		HashMap<String,Sequence> seqNameMap = new HashMap<String,Sequence>(sequenceList.size());
		Iterator<? extends Sequence> seqIt = sequenceList.iterator();
		while(seqIt.hasNext()) {
			Sequence seq = seqIt.next();
			seqNameMap.put(seq.getId(), seq);
		}
		Set<String> seqIds = seqNameMap.keySet();
		
		int seqsToFind = seqIds.size();
		BufferedReader br = new BufferedReader(new FileReader(file));
		String line = null;
		String currentSeqId = null;
		int found = 0;
		boolean inTargetSequence = false;
		while((line = br.readLine()) != null && (found < seqsToFind)) {
			if (line.startsWith(">")) {
				//If we just completed extracting a sequence do some cleanup.
				if(inTargetSequence) {
					inTargetSequence = false;
					found++;
				}
				String [] spaceSeparatedIds = line.substring(1).split("\\s");
				currentSeqId = spaceSeparatedIds[0];
				continue;
			}
			
			Sequence seq = seqNameMap.get(currentSeqId);
			if(seq != null) {
				seq.append(line);
			}

		}		
		//System.out.println("Found " + found + " sequences expected to find " + seqsToFind );
		br.close();		
	}
	
	public void mergeRecords(Sequence seq, String regEx) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(file));
		Pattern pattern = Pattern.compile(regEx);
		String line = null;
		String currentSeqId = null;
		int found = 0;
		boolean inTargetSequence = false;
		while((line = br.readLine()) != null) {
			if (line.startsWith(">")) {
				//If we just completed extracting a sequence do some cleanup.
				if(inTargetSequence) {
					inTargetSequence = false;
					found++;
				}
				currentSeqId = line.substring(1);
				//System.out.println("Seq Id <" + currentSeqId + "> looking for " + regEx);
				Matcher m = pattern.matcher(currentSeqId);
				if(m.find()) {
					//System.out.println("\tmatched!");
					inTargetSequence = true;
				} 
				continue;
			}
			
			if (inTargetSequence) {
				seq.append(line);
			}

		}		
		//System.out.println("Found " + found + " sequences ");
		br.close();	
	}
	
	public void writeRecordsWithMinLength(int minimumLength , InputStream inIs, BufferedWriter bw) throws IOException {
		BufferedReader br = new BufferedReader(new InputStreamReader(inIs));
		String line = null;
		String currentSeqId = null;
		int found = 0;
		boolean inTargetSequence = false;
		Sequence seq = null;
		while((line = br.readLine()) != null) {
			if (line.startsWith(">")) {
				//System.out.println("Old Sequence " + (seq == null ? "none" : seq.getId()) + " length ...  " + (seq == null ? "--" : seq.getLength()) );
				if(seq != null && seq.getLength() > minimumLength) {
					this.write(seq, bw);
				}
				//If we just completed extracting a sequence do some cleanup.
				if(inTargetSequence) {
					inTargetSequence = false;
					found++;
				}
				currentSeqId = line.substring(1);
				seq = new Sequence(currentSeqId);
			}else {
				seq.append(line);
			}

		}		
		br.close();	
	}
	
	public List<Sequence> extractRecordsWithIDsMatching (List<String> regExs, boolean encode) throws IOException {
		return extractRecordsWithIDsMatching(regExs, encode, 0);
	}
	
	public List<Sequence> extractRecordsWithIDsMatching (List<String> regExs, boolean encode, int expectedSize) throws IOException {
		Stack<Sequence> extracted = new Stack<Sequence>();
		BufferedReader br = new BufferedReader(new FileReader(file));
		List<Pattern> patterns = new ArrayList<Pattern>(regExs.size());
		Iterator<String> regExIt = regExs.iterator();
		while(regExIt.hasNext()) {
			Pattern pattern = Pattern.compile(regExIt.next());
			patterns.add(pattern);
		}
		String line = null;
		String currentSeqId = null;
		int found = 0;
		boolean inTargetSequence = false;
		while((line = br.readLine()) != null) {
			if (line.startsWith(">")) {
				//If we just completed extracting a sequence do some cleanup.
				if(inTargetSequence) {
					inTargetSequence = false;
					if(encode) {
						Sequence last = extracted.peek();
						if(last != null) {
							last.encodeSequence();
							last.unloadSequence();
						}
					}
				}
				currentSeqId = line.substring(1);
				Iterator<Pattern> patternIt = patterns.iterator();
				while(patternIt.hasNext()) {
					Pattern pattern = patternIt.next();
					//System.out.println("Seq Id <" + currentSeqId + "> looking for " + pattern.toString());
					Matcher m = pattern.matcher(currentSeqId);
					found++;
					if(m.find()) {
						//System.out.println("\tmatched!");
						Sequence seq = new Sequence(currentSeqId, expectedSize);
						inTargetSequence = true;
						extracted.add(seq);
						break;
					} 
				}
				continue;
			}
			
			if (inTargetSequence) {
				Sequence seq = extracted.peek();
				seq.append(line);
			}

		}		
		//System.out.println("Found " + found + " sequences were looking for " + regExs.size() + " extracted list size " + extracted.size());
		br.close();		
		return extracted;		
	}
	public List<Sequence> extractRecordsWithIDLike(String regEx, boolean encode) throws IOException {
		return extractRecordsWithIDLike(regEx, encode, 0);
	}
	
	public List<Sequence> extractRecordsWithIDLike(String regEx, boolean encode, int expectedSize) throws IOException {
		List<String> oneRegExList = new ArrayList<String>(1);
		oneRegExList.add(regEx);
		return extractRecordsWithIDsMatching(oneRegExList, encode, expectedSize);
	}
	
	public List<Sequence> extractRecords(Collection<String> recordIdList, InputStream is) throws IOException {
		return extractRecords(recordIdList, is, false);
	}
	
	public List<Sequence> extractRecords(Collection<String> recordIdList, InputStream is, boolean largeSequences) throws IOException {
		Stack<Sequence> extracted = new Stack<Sequence>();
		BufferedReader br = new BufferedReader(new InputStreamReader(is));
		String line = null;
		String currentSeqId = null;
		int found = 0;
		boolean inTargetSequence = false;
		System.out.println("to find: " + recordIdList);
		while((line = br.readLine()) != null && (found < recordIdList.size())) {
			if (line.startsWith(">")) {
				//If we just completed extracting a sequence do some cleanup.
				if(inTargetSequence) {
					inTargetSequence = false;
					found++;
				}
				//String [] spaceSeparatedIds = line.substring(1).split("\\s");
				//currentSeqId = spaceSeparatedIds[0];
				currentSeqId = line.substring(1);
				//System.out.println("Seq Id " + currentSeqId);
				if(recordIdList.contains(currentSeqId)) {
					//System.out.println("Found!");
					Sequence seq = new Sequence(currentSeqId);
					inTargetSequence = true;
					extracted.add(seq);
				} 
				continue;
			}
			
			if (inTargetSequence) {
				Sequence seq = extracted.pop();
				seq.append(line);
				extracted.push(seq);
			}

		}		
		//System.out.println("Found " + found + " sequences expected to find " + recordIdList.size() );
		br.close();		
		return extracted;		
	}
	public List<Sequence> extractRecords(Collection<String> recordIdList) throws IOException {
		return extractRecords(recordIdList, new FileInputStream(file));
	}
	
	
	public void extractRegion(SequenceRegion region) throws IOException {
		ArrayList<SequenceRegion> tmpList = new ArrayList<SequenceRegion>(1);
		tmpList.add(region);
		extractRegions(tmpList);
	}
	
	public void extractRegions(List<? extends SequenceRegion> regions) throws IOException{
		extractRegions(regions, false);
	}
	
	/**
	 * Extracts sequence regions from fasta file into 
	 * @param regions
	 * @param useRegionOrientation if true it reverses the sequence extracted.
	 * @throws IOException
	 */
	public void extractRegions(List<? extends SequenceRegion> regions, boolean useRegionOrientation) throws IOException{
		FileInputStream fis = new FileInputStream(file);
		extractRegions(regions, useRegionOrientation, fis);
		fis.close();
	}
	/**
	 * Extracts sequence regions from fasta file into 
	 * @param regions
	 * @param useRegionOrientation if true it reverses the sequence extracted.
	 * @throws IOException
	 */
	public void extractRegions(List<? extends SequenceRegion> regions, boolean useRegionOrientation, InputStream is) throws IOException{
		HashMap<String, List<SequenceRegion>> regionMap = new HashMap<String, List<SequenceRegion>>();
		Iterator<? extends SequenceRegion> it = regions.iterator();
		
		// Group things a bit
		while(it.hasNext()) {
			SequenceRegion seq = it.next();
			List<SequenceRegion> idSeqs = regionMap.get(seq.getContainingSequenceId());
			if(idSeqs == null) {
				idSeqs = new ArrayList<SequenceRegion>();
				regionMap.put(seq.getContainingSequenceId(), idSeqs);
			}
			idSeqs.add(seq);
		}

		//System.out.println(regionMap.keySet());
		// Finally proceed to open the FASTA file and extract sequences for each region
		InputStreamReader isr = new InputStreamReader(is);
		BufferedReader br = new BufferedReader(isr);
		//
		String line = null;
		String currentSeqId = null;
		long currentEndSeqPos = 0;
		List<SequenceRegion> seqRegions = null;
		Iterator<? extends SequenceRegion> rIt = null; 
		while((line = br.readLine()) != null) {
			if (line.startsWith(">")) {
				currentSeqId = line.substring(1).trim();
				//System.out.println("SeqID: <" + currentSeqId + ">");
				currentEndSeqPos = 0;
				seqRegions = regionMap.get(currentSeqId);
				continue;
			}
			
			if (seqRegions == null) {
				continue;
			}
			
			currentEndSeqPos += line.length();
			rIt = seqRegions.iterator();
			SequenceRegion reg = null;
			//System.out.println(" current end pos "+ currentEndSeqPos + " line length " + line.length());
			while(rIt.hasNext()) {
				reg = rIt.next();
				if (reg.getRegionStart() <= currentEndSeqPos  && reg.getRegionEnd() > currentEndSeqPos - line.length()) {
					//System.out.println(currentEndSeqPos + " ... appended ");
					reg.append(line.substring((int) Math.max(0,reg.getRegionStart() - (currentEndSeqPos - line.length() + 1)),
							(int) Math.min(line.length() - 1, reg.getRegionEnd() - (currentEndSeqPos - line.length() + 1))+1));
					//System.out.println(line.substring((int) Math.max(0,reg.getRegionStart() - (currentEndSeqPos - line.length() + 1)),
					//		(int) Math.min(line.length() - 1, reg.getRegionEnd() - (currentEndSeqPos - line.length() + 1))+1));
				}
			}
		}
		
		rIt = regions.iterator();
		while(rIt.hasNext()) {
			SequenceRegion reg = rIt.next();
			if((reg.getStart() + reg.getSequenceBases().length()) < reg.getEnd()) {
				System.out.println("Reset end was " + reg.getEnd() + " but sequence length is " + reg.getSequenceBases().length());
				reg.setEnd(reg.getStart() + reg.getSequenceBases().length() - 1);
			}
			
			if(useRegionOrientation && reg.inReversedOrientation()) {
				reg.reverse();
			}
		}
		
		//br.close();
	}
	
	public void write(List<? extends Sequence> seqs, int lineLength) throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(file));
		Iterator<? extends Sequence> it = seqs.iterator();
		while(it.hasNext()) {
			Sequence seq = it.next();
			write(seq, bw, lineLength);
		}
		bw.close();
	}
	
	public void write(List<? extends Sequence> seqs) throws IOException {
		write(seqs, LINE_LENGTH);
	}

	public void write(Sequence seq) throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(file));
		write(seq, bw);
		bw.close();	
	}
	
	public void write(Sequence seq, String fileName) throws IOException{
		BufferedWriter bw = new BufferedWriter(new FileWriter(fileName));
		write(seq, bw);
		bw.close();		
	}
	
	public void append(List<? extends Sequence> records, String fileName) throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(fileName, true));
		write(records, bw);
		bw.close();
	}
	
	public void append(Sequence seq, String output) throws IOException {
		ArrayList<Sequence> tmpList = new ArrayList<Sequence>(1);
		tmpList.add(seq);
		append(tmpList, output);
	}

	
	public void write(List<? extends Sequence> seqs, String fileName) throws IOException{
		write(seqs, fileName, LINE_LENGTH);
	}
	
	public void write(List<? extends Sequence> seqs, String fileName, int lineLength) throws IOException{
		BufferedWriter bw = new BufferedWriter(new FileWriter(fileName));
		write(seqs, bw, lineLength);
		bw.close();
	}

	public void write(List<? extends Sequence> seqs, BufferedWriter bw, int lineLength) throws IOException {
		Iterator<? extends Sequence> it = seqs.iterator();
		while(it.hasNext()) {
			write(it.next(), bw, lineLength);
		}
	}
	
	public void write(List<? extends Sequence> seqs, BufferedWriter bw) throws IOException {
		write(seqs, bw, LINE_LENGTH);
	}

	//FOR CSF
	public void write(List<? extends Sequence> seqs, BufferedWriter bw, Map<String, String> names, String refID, String pos) throws IOException {
		Iterator<? extends Sequence> it = seqs.iterator();
		while(it.hasNext()) {
			Sequence seq=it.next();
			String name=seq.getId();
			
			
			
			if(names.containsKey(name)){
				String newName=names.get(name);
				if(name.equalsIgnoreCase(refID)){
					newName=newName+"|"+pos+"|";
				}
				seq.setId(newName);
				write(seq, bw);
			}
			
			
		}
	}
	
	public void write(Sequence seq, BufferedWriter bw, int lineLength) throws IOException {
		StringBuilder sequenceBuilder = seq.getSequenceBuilder();
		if(seq == null || sequenceBuilder.length() == 0) {
			return;
		}
		
		if(sequenceBuilder.length() == 0) {
			return;
		}
		//logger.info("Writing sequence "+seq.getId() + " with sequence "+seq.getSequenceBases());
		int currentIndex = 0;
		writeSequenceId(seq, bw);		
		bw.newLine();
		while(currentIndex < sequenceBuilder.length()) {
			int toWrite = Math.min(lineLength, sequenceBuilder.length() - currentIndex) - 1;
			bw.write(sequenceBuilder.substring(currentIndex, currentIndex + toWrite + 1));
			bw.newLine();
			currentIndex = currentIndex + toWrite + 1;
		}
	}
	
	public void write(Sequence seq, BufferedWriter bw) throws IOException {
		write(seq, bw, LINE_LENGTH);
	}

	protected void writeSequenceId(Sequence seq, BufferedWriter bw) throws IOException {
		bw.write(">"+seq.getId());
	}
	/*
	// Make sure that for each sequence ID, the regions to extract are sorted
	// in ascending order based on their coordinates
	Iterator<String> regionIdIt = sequenceRegionIds.iterator();
	while(regionIdIt.hasNext()) {
		List<SequenceRegion> sequenceRegions = regionMap.get(it.next().getId());
		Collections.sort(sequenceRegions, new Comparator<SequenceRegion>() {
			public int compare(SequenceRegion a, SequenceRegion b) {
				long diff = a.getRegionStart() == b.getRegionStart() ? a.getRegionEnd() - b.getRegionEnd() : a.getRegionStart() - b.getRegionStart();
				int result = 0;
				if(diff > Integer.MAX_VALUE) {
					result = Integer.MAX_VALUE;
				} else if(diff < Integer.MIN_VALUE) {
					result = Integer.MIN_VALUE;
				} else {
					result = (int) diff;
				}
				return result;
			}
		});
	}
	*/

	public void unchunk(File[] inDirList, int chunkSize, String unchunkSeqName, String outputFile) throws IOException {
		Sequence unchunked = chunkSize * inDirList.length > 450000000 ? new Sequence(unchunkSeqName, true) : new Sequence(unchunkSeqName, false);

		for(int i = 0; i < inDirList.length; i++ ) {
			System.out.println("loading " + inDirList[i].getAbsolutePath());
			FastaSequenceIO chunkIO = new FastaSequenceIO(inDirList[i]);
			Sequence seq = chunkIO.loadAll().get(0);
			unchunked.append(seq.getSequenceBases());
			seq.unloadSequence();
		}
		
		write(unchunked, outputFile);
	}

	public void setSource(File newSource) {
		this.file = newSource;
	}

	public void breakUpMultifile(final int sizeToBreakup) throws IOException {
		
		FileInputStream fis = new FileInputStream(this.file);
		
		final String outPrefix = file.getAbsolutePath();
		
		
		FastaParser fp = new FastaParser();
		fp.parse(fis, new FastaHandler() {
			int seqNum = 0;
			int grpNum = 0;
			BufferedWriter out = null;
			int base = 0;
			
			public void eof(AbstractFastaParser parser) throws IOException {
				out.close();
			}

			public void newBase(AbstractFastaParser parser) throws IOException {
				out.write(parser.getCurrentBase());
				if(base > 0 && base % LINE_LENGTH == 0) {
					out.newLine();
				}
			}

			public void newSequence(AbstractFastaParser parser) throws IOException {
				if(seqNum % sizeToBreakup == 0) {
					if(out != null) {
						out.close();
					}
					
					out = new BufferedWriter(new FileWriter(outPrefix + "_" + grpNum++));
				}
				
				seqNum++;
				out.write(">" + parser.getCurrentSequenceId());
				out.newLine();
			}
			
		});
	}

	public static String createSizeFile(String refFasta) throws IOException {
		String outFile = refFasta + ".sizes";
		Map<String, Sequence> seqs = getChrSequencesFromFasta(refFasta);
		FileWriter w = new FileWriter(outFile);
		for(String name : seqs.keySet()) {
			w.write(name + "\t" + seqs.get(name).getLength() + "\n");
		}
		w.close();
		return outFile;
	}



}
