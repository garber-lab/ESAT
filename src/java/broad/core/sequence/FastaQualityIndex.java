package broad.core.sequence;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;


public class FastaQualityIndex {
	private HashMap<String, SequenceFileLayout> indexMap;
	long fileSize = 0;
	public FastaQualityIndex() {
		indexMap = new HashMap<String, SequenceFileLayout>();
	}
	
	public void createIndexFromQualityFile(String fastaQualityFile, final long period) throws IOException {
		FastaQualityParser fqp = new FastaQualityParser();

		fqp.parse(new FileInputStream(fastaQualityFile), new FastaQualityHandler() {
			int baseNum = 0;
			
			public void eof(AbstractFastaParser parser) throws IOException {
				baseNum++;
				SequenceFileLayout sfl = indexMap.get(parser.getCurrentSequenceId());
				sfl.addOffset(baseNum, parser.getFileSize());
				sfl.setLength(baseNum);
				fileSize = parser.getFileSize();
			}

			public void newQuality(FastaQualityParser parser) throws IOException {
				if(baseNum % period == 0) {
					SequenceFileLayout sfl = indexMap.get(parser.getCurrentSequenceId());
					sfl.addOffset(baseNum, parser.getOffset());
					//System.out.println(parser.getCurrentSequenceId() + " base " + baseNum + " offset " + parser.getOffset());
				}
				baseNum++;
			}

			public void newSequence(AbstractFastaParser parser) {
				String seqId = parser.getCurrentSequenceId();
				SequenceFileLayout sfl = new SequenceFileLayout(seqId);
				indexMap.put(seqId, sfl);
				baseNum = 0;
			}

			public void newBase(AbstractFastaParser parser) throws IOException {
				// TODO Auto-generated method stub
				
			}			
		});
	}
	
	public void loadFromAscii(String indexFile) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(indexFile));
		
		try {
			String header = br.readLine();
			String [] headerInfo = header.split("\t");
			if("HEADER"  != headerInfo[0]) {
				throw new IOException("Index file format error, first line should  start with HEADER");
			}			
			fileSize = Long.parseLong(headerInfo[2]);
			
			String line;
			while((line = br.readLine()) != null) {
				String[] lineTopInfo = line.split("\\|");
				String seqId = lineTopInfo[0];
				SequenceFileLayout sfl = new SequenceFileLayout(seqId);
				indexMap.put(seqId, sfl);
				
				String [] offsetInfo = lineTopInfo[1].split(" ");
				for (int i = 0; i < offsetInfo.length; i++) {
					String [] basePositionInfo = offsetInfo[i].split("-");
					sfl.addOffset(Integer.parseInt(basePositionInfo[0]), Long.parseLong(basePositionInfo[1]));
				}
			}
			
		}finally {
			br.close();
		}
		
	}

	public void writeAsText(String outFile) throws IOException {
		List<String> seqList = new ArrayList<String>(indexMap.keySet());
		Collections.sort(seqList);
		Iterator<String> seqIt = seqList.iterator();
		BufferedWriter bw = new BufferedWriter(new FileWriter(outFile));
		
		try {
			bw.write("HEADER\t"+seqList.size() +"\t" + fileSize );
			bw.newLine();
			while(seqIt.hasNext()) {
				String seqId = seqIt.next();
				SequenceFileLayout sfl = indexMap.get(seqId);
				bw.write(sfl.toString());
				bw.newLine();
			}
		} finally {
			// should try/catch here but .... 
			bw.close();
		}
	}
	
	public static void main(String [] args) throws Exception {
		String file = args[0];
		int period  = Integer.parseInt(args[1]);
		
		FastaQualityIndex fqp = new FastaQualityIndex();
		fqp.createIndexFromQualityFile(args[0], period);
		fqp.writeAsText(file + ".idx");
	}
	public static class SequenceFileLayout {
		private String seqId;
		LinkedHashMap<Integer, Long> offsetMap;
		int length;
		
		public SequenceFileLayout(String seqId) {
			offsetMap = new LinkedHashMap<Integer, Long>();
			this.seqId = seqId;
		}

		public void setLength(int length) {
			this.length = length;
		}

		public void addOffset(int baseNum, long offset) {
			offsetMap.put(baseNum, offset);
		}
		
		public String toString() {
			StringBuilder sb = new StringBuilder(seqId);
			
			sb.append("|");
			Iterator<Integer> positionIt = offsetMap.keySet().iterator();
			while(positionIt.hasNext()) {
				int pos = positionIt.next();
				sb.append(pos).append("-").append(offsetMap.get(pos));
				sb.append(" ");
			}
			
			return sb.toString();
		}
		
	}

}
