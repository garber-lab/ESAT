package broad.core.util;

import java.util.Comparator;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
//test
public class ChromosomeStringComparator implements Comparator<String> {
	private static final Pattern numberPattern = Pattern.compile("[0-9]{1,2}");
	public int compare(String a, String b) {
		
		int comparison = 0;
		String alean = a.replace("chr", "");
		String blean = b.replace("chr", "");
		
		Matcher am = numberPattern.matcher(alean);
		Matcher bm = numberPattern.matcher(blean);
		if(alean.length() < 3 && blean.length() < 3 && am.matches() && bm.matches()) { // This is dangerous as there could be also letters in the match. We assume for now that data is good.
			comparison =  Integer.valueOf(alean) - Integer.valueOf(blean);
		} else {
			comparison = alean.compareTo(blean); 
		}
		return comparison;
	}

}
