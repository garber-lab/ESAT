/**
 * 
 */
package broad.pda.samtools;

import broad.core.parser.CommandLineParser;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import net.sf.picard.sam.BuildBamIndex;
import net.sf.picard.sam.SortSam;
import net.sf.samtools.BAMFileWriter;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileHeader.SortOrder;
import net.sf.samtools.BAMIndex;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;



/**
 * @author prussell
 *
 */
public class SAMFragmentWriter {
	
	/**
	 * @param input Input sam or bam file
	 * @param output Output bam file
	 * @throws IOException
	 */
	public SAMFragmentWriter(String input, String output) throws IOException {
		this(input, output, null);
	}
	
	/**
	 * @param input Input sam or bam file
	 * @param output Output bam file
	 * @param samHeader Optional sam header if header is too large for SAMFileHeader object
	 * @throws IOException 
	 */
	public SAMFragmentWriter(String input, String output, String samHeader) throws IOException {
		// Establish the input reader
		SAMFileReader reader = new SAMFileReader(new File(input));
		SAMFileHeader header = reader.getFileHeader();
		SAMRecordIterator iter = reader.iterator();
		
		if(samHeader == null) {
			//We are going to write a BAM File directly
			BAMFileWriter writer=new BAMFileWriter(new File(output));
			writer.setSortOrder(SortOrder.coordinate, false);
			writer.setHeader(header);
			while(iter.hasNext()) {
				SAMRecord fragment = SAMPairedEndFileReader.convertToPairedEndFragment(iter.next());
				if(fragment != null) {
					writer.addAlignment(fragment);
					//writer.write(getSamString(fragment));
				}
			}
			writer.close();
			
			//Now build a BAM index
			File transcriptomeBamIdxFile = new File( output + BAMIndex.BAMIndexSuffix);
			if(transcriptomeBamIdxFile.exists()) { transcriptomeBamIdxFile.delete();}
			SAMFileReader reader2 = new SAMFileReader(new File(output));
			BuildBamIndex.createIndex(reader2,transcriptomeBamIdxFile);
			reader2.close();
		}
				
		else {
			//Else we are going to print a warning and write a SAMFile
			FileWriter writer=new FileWriter(output);
			// Get the header from user-defined file
			FileReader hr = new FileReader(samHeader);
			BufferedReader br = new BufferedReader(hr);
			while(br.ready()) {
				String line = br.readLine();
				if(line.substring(0, 1).equals("@")) writer.write(line + "\n");
			}
			hr.close();
			br.close();
			
			while(iter.hasNext()) {
				SAMRecord fragment = SAMPairedEndFileReader.convertToPairedEndFragment(iter.next());
				if(fragment != null) {
					writer.write(getSamString(fragment));
				}
			}
			writer.close();			
		}
		
		//Convert SAM To BAM
		//String cmd="java -jar /seq/mgarber/tools/picard-tools-1.73/SortSAM.jar I="+output+" O="+output+".sorted.bam SO=coordinate";
		
		reader.close();
	}

	/**
	 * Get sam line for record
	 * @param rec The record
	 * @return Sam formatted line ending in newline
	 */
	public static String getSamString(SAMRecord rec) {
		
		// Make string representation of record in sam format
		String rtrn = rec.getReadName() + "\t";
		
		// Get sam flags
		int flags = 0;
		// 0x1 template having multiple segments in sequencing
		if(rec.getReadPairedFlag()) flags += 1;
		// 0x4 segment unmapped
		if(rec.getReadUnmappedFlag()) flags += 4;
		// 0x10 SEQ being reverse complemented
		if(rec.getReadNegativeStrandFlag()) flags += 16;
		// 0x100 secondary alignment
		if(rec.getNotPrimaryAlignmentFlag()) flags += 256;
		// 0x200 not passing quality controls
		if(rec.getReadFailsVendorQualityCheckFlag()) flags += 512;
		// 0x400 PCR or optical duplicate
		if(rec.getDuplicateReadFlag()) flags += 1024;
		
		if(rec.getReadPairedFlag()) {
			// 0x2 each segment properly aligned according to the aligner
			if(rec.getProperPairFlag()) flags += 2;
			// 0x8 next segment in the template unmapped
			if(rec.getMateUnmappedFlag()) flags += 8;
			// 0x20 SEQ of the next segment in the template being reversed
			if(rec.getMateNegativeStrandFlag()) flags += 32;
			// 0x40 the first segment in the template
			if(rec.getFirstOfPairFlag()) flags += 64;
			// 0x80 the last segment in the template
			if(rec.getSecondOfPairFlag()) flags += 128;
		}
		
		rtrn += Integer.valueOf(flags).toString() + "\t";
		
		// The rest of the sam fields
		rtrn += rec.getReferenceName() + "\t";
		rtrn += rec.getAlignmentStart() + "\t";
		rtrn += rec.getMappingQuality() + "\t";
		rtrn += rec.getCigarString() + "\t";
		rtrn += rec.getMateReferenceName() + "\t";
		rtrn += rec.getMateAlignmentStart() + "\t";
		rtrn += rec.getInferredInsertSize() + "\t";
		rtrn += rec.getReadString() + "\t";
		rtrn += rec.getBaseQualityString() + "\n";
		
		return rtrn;
	}
	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		CommandLineParser p = new CommandLineParser();
		p.setProgramDescription("Read a standard paired end alignment file and write an aligmment file with full fragments instead of paired reads.");
		p.addStringArg("-i", "Input sam/bam file", true);
		p.addStringArg("-h", "Sam header file (e.g. if header is too large to be printed from a SAMFileHeader object", false);
		p.addStringArg("-o", "Output sam file", true);
		p.parse(args);
		String input = p.getStringArg("-i");
		String output = p.getStringArg("-o");
		String optionalHeader = p.getStringArg("-h");

		SAMFragmentWriter sfw = new SAMFragmentWriter(input, output, optionalHeader);	
	}

}
