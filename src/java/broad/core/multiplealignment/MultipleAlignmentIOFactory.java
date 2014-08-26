package broad.core.multiplealignment;

public class MultipleAlignmentIOFactory {
	public static MultipleAlignmentIO create(String format) {
		MultipleAlignmentIO maio = null;
		if("FASTA".equalsIgnoreCase(format)) {
			maio = new FastaMultipleAlignmentIO(); 
		} else if("PHYLIP".equalsIgnoreCase(format) ){
			maio = new PhylipInterleavedMultipleAlignmentIO();
		} else if("SEQPHYLIP".equalsIgnoreCase(format) ){
			maio = new PhylipSequencialMultipleAlignmentIO();
		} else if("MAF".equalsIgnoreCase(format) ){
			maio = new MAFIO();
		}else {
			throw new RuntimeException("Unsuported format " + format);
		}
		
		return maio;
	}
}
