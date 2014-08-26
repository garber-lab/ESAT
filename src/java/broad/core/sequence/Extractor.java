package broad.core.sequence;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import broad.core.util.CLUtil;
import broad.core.util.CLUtil.ArgumentMap;

public class Extractor {
	public static String USAGE = "Usage: Extractor TASK=<task_num> <task_args>\n" +
	"\tTasks:\n" +
	"\t\t1. Batch Extract: IN=sequence id list file SEQ=sequence file OUT=output file\n" +
	"\t\t2. Extract Regions: IN=sequence file OUT=output file REGIONS=start1..end1,start2..end2.... SEQID=sequence id if multifasta(use -eq- for = and -s- for spaces REVERSE=Y/N\n UPPERCASE=<Y to uppercase all letters>\n" +
	"\t\t3. Extract records: -in <sequence file> -out <outout file> -seqids <seqId1,seqId2,seqId3,...> -reverse\n" +
	"\t\t4. Extract and Chunk: IN=<sequence file> PREFIX=<output file prefix> CHUNKSIZE=<chunk size> OVERLAP=<chunk overlap> REGION=start1..end1 SEQID=sequence id if multifasta(use -eq- for = and -s- for spaces and -lp- parenthesis and -rp- for right parenthesis> REVERSE=<Y/N> UPPERCASE=<Y to uppercase all letters>\n" +
	"\t\t5. Extract GC percent of all sequences in a multifasta file: IN=<sequence file>\n" +
	"\t\t6. Unchunk Fasta chunks: INDIR=<Directory with sequence chunks it assumes chunks files are named seqname_chunkStart_chunkEnd.extension (specify file extension with the EXT parameter) > \n"+
	"\t\t EXT=<chunk file extension> CHUNKSIZE=<chunk size> OVERLAP=<chunk overlap NOT HANDLING THIS YET>\n" +
	"\t\t7. Extract records matching motif in reverse or direct strand -in <multifasta file> -motif <sequence motif> -out <output file name>\n" +
	"\t\t8. Extract records with vClip preset parameters -in <multifasta>\n" +
	"\t\t9. Extract GC percent for the specified window size: -in <sequence file> -out <GC percent outout file> -windowSize <window size>\n" +
	"\t\t10. Reverse complement sequences in file: -in <sequence file or standard input> -out <file to output reverse-complemented sequences>" +
	"\n\t\t11. Extract sequences that are of at least the given length: -in <sequence file or standard input> -out <file to output or standard out> + -minLength <Minumum length of sequences to return>" +
	"\n\\t\t12. Extract sequences with names matching regular expression -in <sequence file or standard input> -out <file to output or standard out> -pattern <regex patterh>" +
    "\n";


    public Extractor() {
		super();
		// TODO Auto-generated constructor stub
	}

	public static void main(String [] args) throws Exception {
		ArgumentMap argMap = CLUtil.getParameters(args, USAGE);


		if("1".equals(argMap.getTask())) {
			getBatchRegions(argMap.getInput(), argMap.getMandatory("SEQ"), argMap.getOutput());
		} else if("2".equals(argMap.getTask())) {
			boolean reverse = argMap.get("REVERSE") == null ? false : "Y".equals(argMap.get("REVERSE"));
			String regions = argMap.getMandatory("REGIONS");
			String [] regStrings = regions.split(",");
			ArrayList<SequenceRegion> seqs = new ArrayList<SequenceRegion>(regStrings.length);
			for(int i = 0; i < regStrings.length; i++) {
				String [] startEnd = regStrings[i].split("\\.\\.");
				System.out.print("seqId" + argMap.getMandatory("SEQID"));
				String seqId = argMap.getMandatory("SEQID").replaceAll("-eq-","=");
				seqId = seqId.replaceAll("-s-"," ");
				seqId = seqId.replaceAll("-lp-", "(");
				seqId = seqId.replaceAll("-rp-", ")");
				System.out.println(" Transformed seqid <" + seqId +">");
				SequenceRegion seq = new SequenceRegion(seqId);
				seq.setRegionStart(Integer.parseInt(startEnd[0]));
				seq.setRegionEnd(Integer.parseInt(startEnd[1]));
				seqs.add(seq);
			}

			FastaSequenceIO fsIOIn = new FastaSequenceIO(argMap.getInput());
			fsIOIn.extractRegions(seqs);
			if(reverse) {
				Iterator<SequenceRegion> seqIt = seqs.iterator();
				while(seqIt.hasNext()) {
					seqIt.next().reverse();
				}
			}

			if("Y".equalsIgnoreCase(argMap.get("UPPERCASE"))) {
				Iterator<SequenceRegion> seqIt = seqs.iterator();
				while(seqIt.hasNext()) {
					seqIt.next().uppercase();
				}
			}
			FastaSequenceIO fsIOOut = new FastaSequenceIO(argMap.getOutput());
			fsIOOut.write(seqs);
		} else if("3".equals(argMap.getTask())){
			FastaSequenceIO fsio = new FastaSequenceIO(argMap.getInput());
			String [] ids = argMap.getMandatory("seqids").split(",");
			for(int i = 0; i < ids.length; i++) {
				ids[i] = ids[i].replaceAll("-eq-","=");
				ids[i] = ids[i].replaceAll("-s-"," ");
				ids[i] = ids[i].replaceAll("-lp-", "(");
				ids[i] = ids[i].replaceAll("-rp-", ")");
			}
			LinkedHashSet<String> seqIdSet = new LinkedHashSet<String>(ids.length);
			Collections.addAll(seqIdSet, ids);
			List<Sequence> seqs = fsio.extractRecords(seqIdSet);
			if(argMap.containsKey("reverse")) {
				System.out.println("Reversing......");
				Iterator<Sequence> seqIt = seqs.iterator();
				while(seqIt.hasNext()) {
					seqIt.next().reverse();
				}
			}
			fsio.write(seqs, argMap.getOutput());
		} else if("4".equals(argMap.getTask())) {
			boolean reverse = argMap.get("REVERSE") == null ? false : "Y".equals(argMap.get("REVERSE"));
			boolean uppercase = "Y".equalsIgnoreCase(argMap.get("UPPERCASE"));
			String region = argMap.getMandatory("REGION");
			String seqId = argMap.getMandatory("SEQID").replaceAll("-eq-","=");
			int chunkSize = argMap.getInteger("CHUNKSIZE");
			int chunkOverlap = argMap.getInteger("OVERLAP");
			String [] startEnd = region.split("\\.\\.");
			int regionStart = Integer.parseInt(startEnd[0]);
			int regionEnd   = Integer.parseInt(startEnd[1]);
			SequenceRegion seq = new SequenceRegion(seqId);
			seq.setRegionStart(regionStart);
			seq.setRegionEnd(regionEnd);

			FastaSequenceIO fsIOIn = new FastaSequenceIO(argMap.getInput());
			fsIOIn.extractRegion(seq);
			System.out.println("seq length " + seq.getLength());
			if(reverse) {
				seq.reverse();
			}
			if(uppercase) {
				seq.uppercase();
			}
			System.out.println ("Sequence extracted .... ");
			int numOfChunks = (int) Math.floor(seq.getLength()/((float)(chunkSize - chunkOverlap)));
			int chunkStart = 0;
			int chunkEnd   = 0;
			int genomicStart = 0;
			int genomicEnd   = 0;
			for(int i = 0; i< numOfChunks; i++) {
				chunkStart = i * (chunkSize - chunkOverlap);
				chunkEnd   = chunkStart + chunkSize;
				genomicStart = regionStart + chunkStart;
				genomicEnd   = regionStart + chunkEnd - 1;
				SequenceRegion chunk = seq.getRegion(chunkStart, chunkEnd);
				chunk.setId(argMap.getMandatory("PREFIX")+genomicStart + "-" + genomicEnd);
				fsIOIn.write(chunk, argMap.getMandatory("PREFIX") + "_" + genomicStart+ "-" + genomicEnd + ".fa");
				chunk.unloadSequence();
			}
			if(chunkEnd < seq.getLength()) {
				chunkStart = chunkEnd - chunkOverlap;
				genomicStart = regionStart + chunkStart;
				SequenceRegion lastChunk = seq.getRegion(chunkStart, seq.getLength());
				lastChunk.setId(argMap.getMandatory("PREFIX")+ genomicStart + "-" + regionEnd);
				fsIOIn.write(lastChunk, argMap.getMandatory("PREFIX") + "_" + genomicStart + "-" + regionEnd + ".fa");
			}
		} else if ("5".equals(argMap.getTask())) {
			FastaSequenceIO fsio = new FastaSequenceIO(argMap.getInput());
			List<Sequence> seqs = fsio.loadAll();
			Iterator<Sequence> it = seqs.iterator();
			Sequence seq = new Sequence("concatenated seq");

			while(it.hasNext()) {
				seq.append(it.next().getSequenceBases());
			}
			System.out.println("GC%: " + seq.gcContent());
		} else if("6".equals(argMap.getTask())) {
			final String chunkExt = argMap.getMandatory("EXT");
			final Pattern extPat = Pattern.compile("\\." + chunkExt.replace("\\.", "\\.") + "$");
			File inDir = new File(argMap.getInputDir());
			int chunkSize = argMap.getInteger("CHUNKSIZE");
			System.out.println("Ext matching pattern " + extPat);
			File [] inDirList = inDir.listFiles(new FilenameFilter() {

				public boolean accept(File dir, String fileName) {
					Matcher m = extPat.matcher(fileName);
					return m.find();
				}

			});

			Arrays.sort(inDirList, new Comparator<File>(){

				public int compare(File arg0, File arg1) {
					String fName0 = arg0.getName().replace("."+chunkExt, "");
					String fName1 = arg1.getName().replace("."+chunkExt, "");

					String [] arg0Info = fName0.split("_");
					int start0 = Integer.parseInt(arg0Info[arg0Info.length - 2]);
					String [] arg1Info = fName1.split("_");
					int start1 = Integer.parseInt(arg1Info[arg1Info.length - 2]);
					return start0 - start1;
				}

			});

			//Lets check the chunks to make sure we are not missing any.
			for(int i = 0; i < inDirList.length; i++) {
				String [] nameInfo = inDirList[i].getName().split("_");
				int chunkStart = Integer.parseInt(nameInfo[nameInfo.length - 2].replace("."+chunkExt, "").split("-")[0]);
				if(chunkStart != (chunkSize * i )) {
					throw new RuntimeException("Expected chunk starting at " + (chunkSize * i) + " but got " + inDirList[i].getName() + " which starts at " + chunkStart);
				}

			}

			// since we assume a naming convention we know how to get the orginal sequence name
			if(inDirList.length == 0 ) {
				System.err.println("No files in " + inDir.getAbsolutePath() + " matched the given extension " + chunkExt);
				return;
			}
			String seqName = inDirList[0].getName().split("_")[0];
			System.out.println("Seq name " + seqName);
			FastaSequenceIO fsio = new FastaSequenceIO();
			fsio.unchunk(inDirList, chunkSize, seqName, inDir.getAbsolutePath() + "/" + seqName + "." + chunkExt);
		} else if("7".equals(argMap.getTask())){
			String motif = argMap.getMandatory("motif");
			String out   = argMap.getOutput();
			FastaSequenceIO fsio = new FastaSequenceIO(argMap.getInput());

			List<Sequence> matches = new ArrayList<Sequence>();
			Iterator<Sequence> allIt = fsio.loadAll().iterator();
			while(allIt.hasNext()) {
				Sequence seq = allIt.next();
				if(seq.contains(motif)) {
					matches.add(seq);
				}
			}
			fsio.write(matches, out);
		}else if("8".equals(argMap.getTask())){
			String in = argMap.getInput();
			FastaSequenceIO fsio = new FastaSequenceIO(in);
			List<Sequence> all = fsio.loadAll();
			List<Sequence> pass = new ArrayList<Sequence>();

			Iterator<Sequence> allIt = all.iterator();
			int noVector = 0;
			int vectInLowQual = 0;
			while(allIt.hasNext()) {
				Sequence seq = allIt.next();
				String [] idInfo = seq.getId().split("\\s\\|\\s");
				if(idInfo.length < 2) {
					throw new Exception("seq " + seq.getId() + " has an id not containing vClip info ");
				}
				//System.out.println(idInfo[1]);
				String [] infoPairs = idInfo[1].split("\\s");
				HashMap<String, String> vClipInfo = new HashMap<String, String>(infoPairs.length);
				for(int i = 0; i < infoPairs.length; i++) {
					//System.out.println("seq " + seq.getId() + " pair " + infoPairs[i]);
					String[] pair = infoPairs[i].split("=");
					vClipInfo.put(pair[0], pair[1]);
				}

				if("None".equals(vClipInfo.get("LVClip")) || "None".equals(vClipInfo.get("RVClip")) ) {
					System.out.println("Skipping (No left vector or no right vector" + seq.getId());
					noVector++;
					continue;
				}

				int lqclip       = Integer.parseInt(vClipInfo.get("LQClip"));
				int lvclip		 = Integer.parseInt(vClipInfo.get("LVClip"));

				int rqclip       = Integer.parseInt(vClipInfo.get("RQClip"));
				int rvclip		 = Integer.parseInt(vClipInfo.get("RVClip"));

				if( (lqclip > lvclip) || (rqclip < rvclip) ) {
					System.out.println("Skipping (left/right quality clip starts after/before left/right vector clip" + seq.getId());
					vectInLowQual++;
					continue;
				}


				seq.setId(seq.getId().replaceAll(" \\|.*",""));
				pass.add(seq);

			}
			fsio.write(pass, in + ".filtered");

			System.out.println("Reads rejected: No vector " + noVector + " vector inside bad qual " + vectInLowQual);

		} else if("9".equals(argMap.getTask())) {

            String outFile   = argMap.getOutput();
            int windowSize = argMap.getInteger("windowSize");

            FastaSequenceIO fsio = new FastaSequenceIO(argMap.getInput());
            Iterator<Sequence> seqItr = fsio.loadAll().iterator();
            Sequence seq = new Sequence("Concatenated seq");

            while(seqItr.hasNext())
                seq.append(seqItr.next().getSequenceBases());

            BufferedWriter bw = new BufferedWriter(new FileWriter(outFile));

            int idx = 1;
            WindowSlider slider = WindowSlider.getSlider(seq, windowSize, 0);
            while (slider.hasNext())
                bw.append(String.valueOf(idx++)).append("\t").append(String.valueOf(slider.next().gcContent())).append("\n");

            bw.close();

        } else if("10".equals(argMap.getTask())) {
        	
        	InputStream inIs = argMap.getInputStream();
        	BufferedWriter bw = argMap.getOutputWriter();

            FastaSequenceIO fsio = new FastaSequenceIO();
            List<Sequence> sequences = fsio.loadAll(inIs);
            inIs.close();
            Iterator<Sequence> seqItr = sequences.iterator();

            while(seqItr.hasNext()) {
            	Sequence seq = seqItr.next();
            	seq.reverse();
            	seq.setId(seq.getId() + "_rev");
            }
                
            fsio.write(sequences, bw);
            bw.close();

        }else if("11".equals(argMap.getTask())) {
        	
        	InputStream inIs = argMap.getInputStream();
        	BufferedWriter bw = argMap.getOutputWriter();
        	int minLength = argMap.getInteger("minimumLength");

            FastaSequenceIO fsio = new FastaSequenceIO();
            fsio.writeRecordsWithMinLength(minLength,inIs, bw);
            inIs.close();
            bw.close();
        } else if("12".equals(argMap.getTask())){
			String regex = argMap.getMandatory("pattern");
			String out   = argMap.getOutput();
			FastaSequenceIO fsio = new FastaSequenceIO(argMap.getInput());
			Pattern p = Pattern.compile(regex);
		
			List<Sequence> matches = new ArrayList<Sequence>();
			Iterator<Sequence> allIt = fsio.loadAll().iterator();
			while(allIt.hasNext()) {
				Sequence seq = allIt.next();
				Matcher m = p.matcher(seq.getId());
				if (m.find()) {
					matches.add(seq);
				}
			}
			fsio.write(matches, out);
	
		}else {
			System.err.println(USAGE);
		}
	}

	/**
	 * @param seqIdListFile
	 * @param seqFile
	 * @param outFile
	 * @throws FileNotFoundException
	 * @throws IOException
	 */
	private static void getBatchRegions(String seqIdListFile, String seqFile, String outFile) throws FileNotFoundException, IOException {
		BufferedReader br = new BufferedReader(new FileReader(seqIdListFile));

		TreeSet<String> seqIdList = new TreeSet<String>();
		List<Sequence> extractedSeqs = null;
		String line = null;
		while((line = br.readLine()) != null) {
			if(line.startsWith("#")){
				continue;
			}
			String seqId = line.trim();
			seqIdList.add(seqId);
		}
		FastaSequenceIO fsio = new FastaSequenceIO(seqFile);
		extractedSeqs = fsio.extractRecords(seqIdList);
		fsio.write(extractedSeqs, outFile);
	}


}
