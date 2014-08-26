package broad.core.annotation;

import java.util.ArrayList;
import java.util.List;
import java.util.Iterator;

import broad.core.error.ParseException;

public class ShortBED extends BasicGenomicAnnotation {

	public ShortBED(String name) {
		super(name);
	}
	
	public ShortBED(LightweightGenomicAnnotation anot) {
		super (anot);
	}
	
	public ShortBED(String name, String chr, int start, int end) {
		super(name, chr, start, end);
	}
	
	
	public ShortBED(String [] info) throws ParseException {
		super("", info[0], Integer.parseInt(info[1]), Integer.parseInt(info[2]));
		
		if(info.length == 3) {
			setName(info[0] + ":" + info[1] + "-" + info[2]);
		} else {
			setName(info[3]);
		} 
		
		/*
		if(info.length > 4) {
			setScore(Double.parseDouble(info[4]));
		}
		if(info.length > 5) {
			setOrientation(info[5]);
		}
			
		
		if(info.length > 6) {
			for(int i = 6; i < info.length; i++) {
				addExtraScore(Double.valueOf(info[i]));
			}
		}*/
	}
	
	public String toString() { return toString(false);} 
	public String toString(boolean setNegativeScoresTo0) {
		StringBuffer buf = new StringBuffer(getChromosome());
		buf.append("\t")
			.append(getStart())
			.append("\t")
			.append(getEnd())
			.append("\t")
			.append(getName())
			.append("\t")
			.append(getScore() < 0 && setNegativeScoresTo0 ? "0" : getScore())
			.append("\t")
			.append(getOrientation());

		
		if(getExtraScores() != null ) {
			for(double score : getExtraScores()) {
				buf.append("\t").append(String.valueOf(score));
			}
		}
		return buf.toString();
 	}
	
	public String toShortString() {
		StringBuffer buf = new StringBuffer(getChromosome());
		buf.append("\t")
			.append(getStart())
			.append("\t")
			.append(getEnd())
			.append("\t")
			.append(getName())
			.append("\t")
			.append(getScore())
			.append("\t")
			.append(getOrientation());

		return buf.toString();
	}
	
	public String toWIGString() {
		StringBuffer buf = new StringBuffer("chr");
		buf.append(getChromosome())
			.append("\t")
			.append(getStart())
			.append("\t")
			.append(getEnd())
			.append("\t")
			.append(Math.round(getScore()) );
		
		return buf.toString();
	}



}
