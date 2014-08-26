/**
 * 
 */
package broad.pda.countreads;

import broad.core.parser.CommandLineParser;
import broad.core.parser.StringParser;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.TreeSet;
import org.apache.log4j.Logger;


import net.sf.picard.sam.DuplicationMetrics;

/**
 * @author prussell
 *
 */
public class FastqLibraryStats {

	/**
	 * Constructor for single end reads
	 * @param singleReadFastq
	 * @throws IOException 
	 */
	public FastqLibraryStats(String singleReadFastq) throws IOException {
		this(singleReadFastq, null);
	}
	
	/**
	 * Constructor for paired end reads
	 * @param read1fastq
	 * @param read2fastq
	 * @throws IOException 
	 */
	public FastqLibraryStats(String read1fastq, String read2fastq) throws IOException {
		read1file = read1fastq;
		read2file = read2fastq;
		uniqueReads = new TreeSet<String>();
		duplicatedReads = new TreeSet<String>();
		estLibrarySize = -1;
		numUniqueReads = -1;
		totalReads = -1;
		pctDup = -1;
		readsPaired = false;
		if(read2fastq != null) readsPaired = true;
		collapseReadsAndCount();
	}
	
	private String read1file;
	private String read2file;
	private TreeSet<String> uniqueReads;
	private TreeSet<String> duplicatedReads;
	private int numUniqueReads;
	private double pctDup;
	private int totalReads;
	private long estLibrarySize;
	private boolean readsPaired;
	private static Logger logger = Logger.getLogger(FastqLibraryStats.class.getName());
	
	/**
	 * Get the number of unique reads
	 * @return The number of unique reads
	 */
	public int getNumUniqueReads() {
		return numUniqueReads;
	}
	
	/**
	 * Get the proportion of reads that are duplicates
	 * @return The proportion of reads that are duplicated
	 */
	public double getPercentDuplicated() {
		return pctDup;
	}
	
	/**
	 * Get the total number of reads
	 * @return The total number of reads
	 */
	public int getTotalReads() {
		return totalReads;
	}
	
	/**
	 * Get the estimated library size based on number of reads and percent duplicates
	 * @return The estimated library size
	 */
	public long getEstimatedLibrarySize() {
		return estLibrarySize;
	}
	
	/**
	 * Count unique reads by collapsing identical reads
	 * @param outFile Output file for statistics
	 * @param readsPaired Whether the reads are paired
	 * @throws IOException
	 */
	private void collapseReadsAndCount() throws IOException {
		
		// If paired reads, read read2 file. Otherwise just read read1 file.
		String read2 = read2file;
		if(!readsPaired) read2 = read1file;
		
		FileReader reader1 = new FileReader(read1file);
		BufferedReader buffered1 = new BufferedReader(reader1);
				
		FileReader reader2 = new FileReader(read2);
		BufferedReader buffered2 = new BufferedReader(reader2);
		
		int linesRead = 0;
		totalReads = 0;
		
		while(buffered1.ready()) {
			String read = buffered1.readLine();
			if(readsPaired) {
				read += "_" + buffered2.readLine();
			}
			linesRead++;
			if(linesRead % 4 == 2) {
				totalReads++;
				if(uniqueReads.contains(read)) duplicatedReads.add(read);
				uniqueReads.add(read);
			}
		}
		
		reader1.close();
		reader2.close();
		buffered1.close();
		buffered2.close();
		
		numUniqueReads = uniqueReads.size();
		try {
			estLibrarySize = DuplicationMetrics.estimateLibrarySize(totalReads, numUniqueReads).longValue();
		} catch (NullPointerException e) {
			logger.info("Warning: caught NullPointerException. Total reads = " + totalReads + ". Unique reads = " + numUniqueReads + ".");
		}
		pctDup = ((double)totalReads - (double)numUniqueReads)/totalReads;
	}
	
	/**
	 * Write fastq files of unique reads and duplicated reads
	 * @param outUniquePrefix Output fastq file of unique reads
	 * @param outDupPrefix Output fastq file of duplicated reads
	 * @param pairedReads Whether the reads are paired
	 * @throws IOException 
	 */
	public void writeSeparateFiles(String outUniquePrefix, String outDupPrefix, boolean pairedReads) throws IOException {
		
		if(duplicatedReads.isEmpty()) {
			throw new IllegalStateException("Duplicated read pair set is empty. Try calling collapseReads() first.");
		}
		
		// If paired reads, read read2 file. Otherwise just read read1 file.
		String read2;
		if(pairedReads) read2 = read2file;
		else read2 = read1file;
		
		FileReader reader1 = new FileReader(read1file);
		FileReader reader2 = new FileReader(read2);
		BufferedReader buffered1 = new BufferedReader(reader1);
		BufferedReader buffered2 = new BufferedReader(reader2);
		
		if(pairedReads) {
		
			FileWriter ou1 = new FileWriter(outUniquePrefix + "_1.fq");
			FileWriter ou2 = new FileWriter(outUniquePrefix + "_2.fq");
			FileWriter od1 = new FileWriter(outDupPrefix + "_1.fq");
			FileWriter od2 = new FileWriter(outDupPrefix + "_2.fq");
		
			logger.info("Writing unique reads to " + outUniquePrefix + "_1.fq" + " and " + outUniquePrefix + "_2.fq" + ".");
			logger.info("Writing duplicated reads to " + outDupPrefix + "_1.fq" + " and " + outDupPrefix + "_2.fq" + ".");
			
			while(buffered1.ready()) {
				String read1Line1 = buffered1.readLine();
				String read2Line1 = buffered2.readLine();
				String read1Line2 = buffered1.readLine();
				String read2Line2 = buffered2.readLine();
				String read1Line3 = buffered1.readLine();
				String read2Line3 = buffered2.readLine();
				String read1Line4 = buffered1.readLine();
				String read2Line4 = buffered2.readLine();
			
				String combined = read1Line2 + "_" + read2Line2;
				if(duplicatedReads.contains(combined)) {
					od1.write(read1Line1 + "\n");
					od1.write(read1Line2 + "\n");
					od1.write(read1Line3 + "\n");
					od1.write(read1Line4 + "\n");
					od2.write(read2Line1 + "\n");
					od2.write(read2Line2 + "\n");
					od2.write(read2Line3 + "\n");
					od2.write(read2Line4 + "\n");
				} else {
					ou1.write(read1Line1 + "\n");
					ou1.write(read1Line2 + "\n");
					ou1.write(read1Line3 + "\n");
					ou1.write(read1Line4 + "\n");
					ou2.write(read2Line1 + "\n");
					ou2.write(read2Line2 + "\n");
					ou2.write(read2Line3 + "\n");
					ou2.write(read2Line4 + "\n");
				}
			
			}
		
			buffered1.close();
			buffered2.close();
			ou1.close();
			ou2.close();
			od1.close();
			od2.close();
		
		}
		
		else {
			
			FileWriter ou1 = new FileWriter(outUniquePrefix + ".fq");
			FileWriter od1 = new FileWriter(outDupPrefix + ".fq");
		
			logger.info("Writing unique reads to " + outUniquePrefix + ".fq.");
			logger.info("Writing duplicated reads to " + outDupPrefix + ".fq.");
			
			while(buffered1.ready()) {
				String read1Line1 = buffered1.readLine();
				@SuppressWarnings("unused")
				String read2Line1 = buffered2.readLine();
				String read1Line2 = buffered1.readLine();
				String read2Line2 = buffered2.readLine();
				String read1Line3 = buffered1.readLine();
				@SuppressWarnings("unused")
				String read2Line3 = buffered2.readLine();
				String read1Line4 = buffered1.readLine();
				@SuppressWarnings("unused")
				String read2Line4 = buffered2.readLine();
			
				String combined = read1Line2 + "_" + read2Line2;
				if(duplicatedReads.contains(combined)) {
					od1.write(read1Line1 + "\n");
					od1.write(read1Line2 + "\n");
					od1.write(read1Line3 + "\n");
					od1.write(read1Line4 + "\n");
				} else {
					ou1.write(read1Line1 + "\n");
					ou1.write(read1Line2 + "\n");
					ou1.write(read1Line3 + "\n");
					ou1.write(read1Line4 + "\n");
				}
			
			}
		
			buffered1.close();
			buffered2.close();
			ou1.close();
			od1.close();

		}
		
	}
	
	/**
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-f", "Fastq list file. Each line: <sample_name> <fastq_1> <fastq_2>(optional)", true);
		p.addStringArg("-o", "Output table file", true);
		p.parse(args);
		String fastqList = p.getStringArg("-f");
		String outFile = p.getStringArg("-o");
		
		FileReader r = new FileReader(fastqList);
		BufferedReader b = new BufferedReader(r);
		FileWriter w = new FileWriter(outFile);
		StringParser s = new StringParser();
		
		String header = "sample_name\t";
		header += "total_reads\t";
		header += "unique_reads\t";
		header += "pct_duplicated\t";
		
		w.write(header + "\n");
		
		while(b.ready()) {
			String line = b.readLine();
			s.parse(line);
			int numFields = s.getFieldCount();
			if(numFields == 0) continue;
			if(numFields < 2 || numFields > 3) {
				throw new IllegalArgumentException("Each line in fastq list file must be of the form: <sample_name> <fastq_1> <fastq_2>(optional)");
			}
			String sampleName = s.asString(0);
			String lineToWrite = sampleName + "\t";
			logger.info("Processing sample " + sampleName);
			if(numFields == 2) {
				FastqLibraryStats d = new FastqLibraryStats(s.asString(1));
				lineToWrite += d.getTotalReads() + "\t";
				lineToWrite += d.getNumUniqueReads() + "\t";
				lineToWrite += d.getPercentDuplicated() + "\t";
			} else {
				FastqLibraryStats d = new FastqLibraryStats(s.asString(1), s.asString(2));
				lineToWrite += d.getTotalReads() + "\t";
				lineToWrite += d.getNumUniqueReads() + "\t";
				lineToWrite += d.getPercentDuplicated() + "\t";
			}
			w.write(lineToWrite + "\n");
		}
		
		r.close();
		b.close();
		w.close();
		
	}

}
