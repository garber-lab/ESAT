package broad.core.sequence;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.RandomAccessFile;
import java.sql.Blob;
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;


public class FastaQualityBinaryIO {
	public static final int SEQ_ID_REC_LENGTH = 300;
	public static final int BUFFER_SIZE       = 8192;
	
	public void createBinaryFromAscii(InputStream in, String outFile) throws IOException {
		FastaQualityParser fqp = new FastaQualityParser();
		final RandomAccessFile  out = new RandomAccessFile(outFile, "rw");
		final HashMap<String, Long> sequenceOffsetMap = new HashMap<String, Long>();

		
		fqp.parse(in, new FastaQualityHandler() {
			int baseNum = 0;
			long sequenceStartOffset = 0;
			byte[] buffer = new byte [BUFFER_SIZE];
			int bufferFill = 0;
			
			public void eof(AbstractFastaParser parser) throws IOException {
				baseNum++;
				out.write(buffer, 0, bufferFill);
				out.seek(sequenceStartOffset + SEQ_ID_REC_LENGTH);
				out.writeInt(baseNum);
			}

			public void newQuality(FastaQualityParser parser) throws IOException {
				buffer[bufferFill++] = (byte) parser.getCurrentQuality();
				/*
				if(baseNum > 0 && baseNum % 26 == 0)  {
					System.out.print("\n");
				}
				*/
				//System.out.print(parser.getCurrentQuality() + " ");
				if(bufferFill == BUFFER_SIZE) {
					out.write(buffer);
					//System.out.println("Writing buffer: " + new String(buffer) + " ___end");
					bufferFill = 0;
				}
				baseNum++;
			}

			public void newSequence(AbstractFastaParser parser) throws IOException {
				
				long currentOffset = out.getFilePointer();
				if(out.getFilePointer() > 0) { 
					// Write current buffer first
					out.write(buffer, 0, bufferFill);
					out.writeChar('\n');
					bufferFill = 0;
					currentOffset = out.getFilePointer();
					
					//System.out.println("Writing prior sequence length " + baseNum + " at offset " + (sequenceStartOffset + SEQ_ID_REC_LENGTH));
					out.seek(sequenceStartOffset + SEQ_ID_REC_LENGTH);
					out.writeInt(baseNum);
					out.seek(currentOffset);
				}
				String seqId = parser.getCurrentSequenceId();
				char [] seqIdChars = seqId.toCharArray();
				int bytesToWrite = Math.min(SEQ_ID_REC_LENGTH, seqIdChars.length);
				sequenceStartOffset = currentOffset;
				
				sequenceOffsetMap.put(seqId, currentOffset);
				
				for(int i = 0; i < bytesToWrite; i++) {
					out.write((byte) seqIdChars[i]);
				} 
				
				for(int i = 0; i < SEQ_ID_REC_LENGTH - bytesToWrite; i++) {
					out.write((byte) ' ');
				}
				//System.out.println("Processed " + seqId + " current offset " + out.getFilePointer());
				out.skipBytes(4); // make room for sequence length!

				baseNum = 0;
			}

			public void newBase(AbstractFastaParser parser) throws IOException {
				// TODO Auto-generated method stub
				
			}			
		});
		
		out.close();
		
		BufferedWriter bw = new BufferedWriter(new FileWriter(outFile + ".idx"));
		List<String> sortedSeqIds = new ArrayList<String>(sequenceOffsetMap.keySet());
		Collections.sort(sortedSeqIds);
		Iterator<String> seqIdIt = sortedSeqIds.iterator();
		while(seqIdIt.hasNext()) {
			String seqId = seqIdIt.next();
			bw.write(seqId + "\t" + sequenceOffsetMap.get(seqId));
			bw.newLine();
		}
		bw.close();
	}
	
	
	public void toDBFromAscii(InputStream in, final Connection conn, String table) throws IOException, SQLException {
		String insSQL = "INSERT INTO " + table + "(scaffold, quality) VALUES (?, ?)";
		String selSQL = "SELECT quality FROM " + table + " WHERE scaffold = ?";
		String updSQL = "UPDATE " + table + 
			" SET quality = ?, gap = ?, q10 = ?, q15 = ?, q20 = ?, q25  =?, q30 = ? " + 
			" WHERE scaffold = ?";
		
		FastaQualityParser fqp = new FastaQualityParser();

		final PreparedStatement insertCS = conn.prepareStatement(insSQL);
		final PreparedStatement selectCS = conn.prepareStatement(selSQL);
		final PreparedStatement updateCS = conn.prepareStatement(updSQL);
		
		fqp.parse(in, new FastaQualityHandler() {
			int baseNum = 0;
			byte[] buffer = new byte [BUFFER_SIZE];
			int bufferFill = 0;
			int seqNum = 0;
			String workingSeqId = null;
			OutputStream blobOS = null;
			Blob qualityBlob = null;
			
			int gapNum = 0;
			int q10OrLessNum = 0;
			int q15OrLessNum = 0;
			int q20OrLessNum = 0;
			int q25OrLessNum = 0;
			int q30OrLessNum = 0;
			
			public void eof(AbstractFastaParser parser) throws IOException {
				baseNum++;
				blobOS.write(buffer, 0, bufferFill);
				blobOS.close();
		
				try {
					updateCS.setString(8, workingSeqId);
					updateCS.setBlob(1, qualityBlob);
					updateCS.setInt(2, gapNum);
					updateCS.setInt(3, q10OrLessNum);
					updateCS.setInt(4, q15OrLessNum + q10OrLessNum);
					updateCS.setInt(5, q20OrLessNum + q15OrLessNum + q10OrLessNum);
					updateCS.setInt(6, q25OrLessNum + q20OrLessNum + q15OrLessNum + q10OrLessNum);
					updateCS.setInt(7, q30OrLessNum + q25OrLessNum + q20OrLessNum + q15OrLessNum + q10OrLessNum);
					
					updateCS.execute();
					conn.commit();
				}catch (SQLException sqle) {
					System.out.println(sqle.getStackTrace());
					throw new IOException(sqle.getMessage());
				}
				System.out.println("Got EOF writing last buffer ");
			}

			public void newQuality(FastaQualityParser parser) throws IOException {
				int qual = parser.getCurrentQuality();
				buffer[bufferFill++] = (byte) qual;
				/*
				if(baseNum > 0 && baseNum % 26 == 0)  {
					System.out.print("\n");
				}
				*/
				//System.out.print(parser.getCurrentQuality() + " ");
				if(bufferFill == BUFFER_SIZE) {
					blobOS.write(buffer);
					//System.out.println("Writing buffer: " );//+ new String(buffer) + " ___end");
					bufferFill = 0;
				}
				
				if (qual == 0) {
					gapNum++;
				} else if(qual <= 10) {
					q10OrLessNum++;
				} else if(qual <= 15) {
					q15OrLessNum++;
				} else if (qual <= 20) {
					q20OrLessNum++;
				} else if (qual <= 25) {
					q25OrLessNum++;
				} else if (qual <= 30) {
					q30OrLessNum++;
				}
				
				baseNum++;
			}

			public void newSequence(AbstractFastaParser parser) throws IOException{
				try {
					if(seqNum > 0) { 
						System.out.println("Got new sequence, writing last " + workingSeqId + " buffer ");
						blobOS.write(buffer, 0, bufferFill);
						blobOS.flush();
						blobOS.close();
						
						updateCS.setString(8, workingSeqId);
						updateCS.setBlob(1, qualityBlob);
						updateCS.setInt(2, gapNum);
						updateCS.setInt(3, q10OrLessNum);
						updateCS.setInt(4, q15OrLessNum + q10OrLessNum);
						updateCS.setInt(5, q20OrLessNum + q15OrLessNum + q10OrLessNum);
						updateCS.setInt(6, q25OrLessNum + q20OrLessNum + q15OrLessNum + q10OrLessNum);
						updateCS.setInt(7, q30OrLessNum + q25OrLessNum + q20OrLessNum + q15OrLessNum + q10OrLessNum);
						
						updateCS.execute();
						conn.commit();
						bufferFill = 0;
					}
					
					String seqId = parser.getCurrentSequenceId();
					
					insertCS.setString(1, seqId);
					insertCS.setBytes(2, new byte[0]);
					//insertCS.setByte(2, (byte) 32);
					insertCS.execute();
					
					selectCS.setString(1, parser.getCurrentSequenceId());
					selectCS.executeQuery();
					ResultSet rs = selectCS.getResultSet();
					rs.next();
					qualityBlob = rs.getBlob(1);
					blobOS = qualityBlob.setBinaryStream(1);
					seqNum++;
					baseNum = 0;
					gapNum = 0;
					q10OrLessNum = 0;
					q15OrLessNum = 0;
					q20OrLessNum = 0;
					q25OrLessNum = 0;
					q30OrLessNum = 0;
					
					workingSeqId = seqId;
				} catch (SQLException sqle) {
					sqle.printStackTrace(System.out);
					throw new IOException(sqle.getMessage());
				}
			}

			public void newBase(AbstractFastaParser parser) throws IOException {
				// TODO Auto-generated method stub
				
			}			
		});
		

	}

	
	public void createAsciiFromBinary(String binaryQualityFile, int numQualBasesPerRow) throws IOException {
		RandomAccessFile raf = new RandomAccessFile(binaryQualityFile, "r");
		byte[] seqIdBytes = new byte[SEQ_ID_REC_LENGTH];
		System.out.println("offset " + raf.getFilePointer());
		raf.read(seqIdBytes);
		System.out.println("offset after reading seq name " + raf.getFilePointer());
		String seqId = new String(seqIdBytes).trim();
		int seqIdLength = raf.readInt();
		System.out.println("offset after reading seq length " + raf.getFilePointer());
		byte [] firstQual = new byte[seqIdLength];
		raf.read(firstQual);
		
		System.out.println("got seq id " + seqId + " of length " + seqIdLength);
		for (int i = 0; i < firstQual.length; i++) {
			if(i > 0 && i % 26 == 0)  {
				System.out.print("\n");
			}
			System.out.print((int) firstQual[i] + " ");
		}
		raf.close();
	}
	
	
	/**
	 * Retrieves quality base calls for a given sequence from a binary quality assembly file.
	 * @param Table - Database table containing the assembly quality info.
	 * @param seqId - Sequence id of interest.
	 * @param start - Zero based start of base (included)
	 * @param numOfBases -  number of bases to retrieve starting from base start.
	 * @return
	 * @throws SQLException 
	 */
	public int[] getQualityBases(Connection conn, String table , String seqId, int start, int numOfBases) throws SQLException {
		int [] qualities  = null;
		String selSQL = "SELECT quality FROM " + table + " WHERE scaffold = ?";
		final PreparedStatement selectCS = conn.prepareStatement(selSQL);
		
		selectCS.setString(1, seqId);
		ResultSet rs = selectCS.executeQuery();
		
		if(rs.next()) {
			Blob quality = rs.getBlob(1);
			byte [] qualBytes = quality.getBytes(start + 1, numOfBases);
			qualities = new int[qualBytes.length];
			for(int i = 0; i < qualBytes.length; i++) {
				qualities[i] = (int) qualBytes[i];
			}
		}
		return qualities;
	}
	
	/**
	 * Retrieves quality base calls for a given sequence from a binary quality assembly file.
	 * Use wisely, it returns an integer array with ALL quality scores.
	 * @param Table - Database table containing the assembly quality info.
	 * @param seqId - Sequence id of interest.
	 * @return
	 * @throws SQLException 
	 */
	public int[] getQualityBases(Connection conn, String table , String seqId) throws SQLException {
		int [] qualities  = null;
		String selSQL = "SELECT quality FROM " + table + " WHERE scaffold = ?";
		final PreparedStatement selectCS = conn.prepareStatement(selSQL);
		
		selectCS.setString(1, seqId);
		ResultSet rs = selectCS.executeQuery();
		
		if(rs.next()) {
			Blob quality = rs.getBlob(1);
			long size = quality.length();
			byte [] qualBytes = quality.getBytes(1, (int) size);
			qualities = new int[(int) size];
			for(int i = 0; i < qualBytes.length; i++) {
				qualities[i] = (int) qualBytes[i];
			}
		}
		return qualities;
	}
	
}
