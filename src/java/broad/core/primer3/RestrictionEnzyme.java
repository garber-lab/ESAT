package broad.core.primer3;

import java.util.List;
import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public  class RestrictionEnzyme {
	private String name;
	private String restrictionSite;
	private Pattern restrictionSitePattern;
	private int fivePrimeCutFromSite;
	private int threePrimerCutFromSize;

	public static String[][] WOBLE_PATTERN_MAP = {
		{"M",  "[AC]"},
		{"R",  "[AG]"},
		{"S",  "[CG]"},
		{"V",  "[AGC]"},
		{"W",  "[AT]"},
		{"Y",  "[CT]"},
		{"H",  "[ACT]"},
		{"K" , "[GT]"},
		{"D" , "[AGT]"},
		{"B" , "[TGC]"},
		{"N" , "[ATGC]"}
	};
	
	public static RestrictionEnzyme AcuI = new RestrictionEnzyme("AcuI", "CTGAAG", 17, 15);
	public static RestrictionEnzyme BpmI = new RestrictionEnzyme("BpmI", "CTGGAG", 16, 14);
	public static RestrictionEnzyme BpuEI = new RestrictionEnzyme("BpuEI", "CTTGAG", 16, 14);
	public static RestrictionEnzyme BsgI = new RestrictionEnzyme("BsgI", "GTGCAG", 16, 14);
	public static RestrictionEnzyme EcoP15I = new RestrictionEnzyme("EcoP15I", "CAGCAG", 27, 27);
	public static RestrictionEnzyme MmeI = new RestrictionEnzyme("MmeI", "TCCRAC", 20, 18);
	public static RestrictionEnzyme BbvI = new RestrictionEnzyme("BbvI", "GCAGC", 8, 12);
	public static RestrictionEnzyme BceAI = new RestrictionEnzyme("BceAI", "ACGGC", 12, 14);
	public static RestrictionEnzyme BsmFI = new RestrictionEnzyme("BsmFI", "GGGAC", 10, 14);
	public static RestrictionEnzyme BtgZI = new RestrictionEnzyme("BtgZI", "GCGATG", 10, 14);
	public static RestrictionEnzyme EciI = new RestrictionEnzyme("EciI", "GGCGGA", 11, 9);
	public static RestrictionEnzyme FokI = new RestrictionEnzyme("FokI", "GGATG", 9, 13);
	
	
	public RestrictionEnzyme(String name, String restrictionSite, int fivePrimeCutFromSite, int threePrimerCutFromSite) {
		super();
		this.name = name;
		this.restrictionSite = restrictionSite;		
		this.fivePrimeCutFromSite = fivePrimeCutFromSite;
		this.threePrimerCutFromSize = threePrimerCutFromSite;
		
		String pattern = restrictionSite;
		for(int i = 0; i < WOBLE_PATTERN_MAP.length; i++) {
			pattern = pattern.replaceAll(WOBLE_PATTERN_MAP[i][0], WOBLE_PATTERN_MAP[i][1]);
		}
		this.restrictionSitePattern = Pattern.compile(pattern);
		//System.out.println("Restrinction site: " + restrictionSite + " restrictionSitePattern " + restrictionSitePattern);
	}

	public List<String> digest(String seq) {
		ArrayList<String> digestion = new ArrayList<String>();
		String sequence = seq.toUpperCase();
		int lastCut = 0;
		Matcher m = restrictionSitePattern.matcher(sequence);
		System.out.println("Product " + sequence + " site " + restrictionSitePattern);
		while(m.find()) {			
			int cut = Math.min(m.end() + fivePrimeCutFromSite + 1, sequence.length());
			System.out.println("Found cutting site, will cut from " + lastCut + " to min ("+(m.end() + fivePrimeCutFromSite + 1)+","+( sequence.length() - 1));
			digestion.add(sequence.substring(lastCut, cut));
			lastCut = cut;
		}
		digestion.add(sequence.substring(lastCut));
		return digestion;
	}

	public String getName() {return name;}
	
	

}
