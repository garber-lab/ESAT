package broad.core.annotation;

import java.util.ArrayList;
import java.util.List;
import java.util.Iterator;

import nextgen.core.annotation.Annotation;

import broad.core.error.ParseException;

public class BED extends BasicGenomicAnnotation {
	private String rgb;
	private int thickStart;
	private int thickEnd;
	private int [][] blockSizeStarts;
	private GenomicAnnotation region;
	
	public BED(String name) {
		super(name);
	}
	
	public BED(LightweightGenomicAnnotation anot) {
		super (anot);
	}
	
	public BED(String name, String chr, int start, int end) {
		super(name, chr, start, end);
	}
	
	public BED(BED anot) {
		super (anot);
		thickStart = anot.thickStart;
		thickEnd   = anot.thickEnd;
		blockSizeStarts = anot.blockSizeStarts;
	}
	
	public BED(Annotation anot) {
		super(anot);
	}
	
	public BED(String [] info) throws ParseException {
		super();
		setChromosome(info[0]); //try to handle both chrNN and NN notations for chromosomes
		setStart(Integer.parseInt(info[1]));
		setEnd(Integer.parseInt(info[2]));
		
		if(info.length == 3) {
			setName(info[0] + ":" + info[1] + "-" + info[2]);
		} else {
			setName(info[3]);
			
		} 
		
		if(info.length > 4) {
			setScore(Double.parseDouble(info[4]));
		}
		if(info.length > 5) {
			setOrientation(info[5]);
			
			if(info.length == 7) { // assume that color is the last of the fields if the line is short
				rgb = info[6];
			} else {
				thickStart = info.length > 7 && info[6].length() > 0 ? Integer.parseInt(info[6]) : 0;
				thickEnd   = info.length >= 8 && info[7].length() > 0 ? Integer.parseInt(info[7]) : 0;
				rgb        = info.length >= 9 && info[8].length() > 0 ? info[8] :  null;
				if(info.length > 9) {
					int blockCount = Integer.parseInt(info[9]);
					if(blockCount > 0) {
						blockSizeStarts = new int[blockCount][2];
						info[10] = info[10].replaceAll("\"", "");
						info[11] = info[11].replaceAll("\"", "");
						info[10] = info[10].replaceAll(" ", "");//Moran added
						info[11] = info[11].replaceAll(" ", "");//Moran added
					
						String [] starts = info[10].split(",");
						String [] ends   = info[11].split(",");
						if (starts.length < blockCount || ends.length < blockCount) {
							throw new ParseException("BAD BED FORMAT apparently the number of start ("+info[10] +") and end ("+info[11] +
								") items does not agree with te blockCount " + info[9]);
						}
						for(int i = 0; i < blockCount; i++) {
							blockSizeStarts[i][0] = Integer.parseInt(starts[i]);
							blockSizeStarts[i][1] = Integer.parseInt(ends[i]);
						}
					}
				}
			}
		} 
		if(info.length > 12) {
			for(int i = 12; i < info.length; i++) {
				addExtraScore(Double.valueOf(info[i]));
			}
		}
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
			.append(getOrientation())
			.append("\t")
			.append(thickStart > 0 ? thickStart : getStart())
			.append("\t")
			.append(thickEnd > 0 ? thickEnd : getEnd())
			.append("\t")
			.append(rgb != null ? rgb : "0,0,0");
		
		if(blockSizeStarts != null && blockSizeStarts.length > 0) {
			buf.append("\t").append(blockSizeStarts.length).append("\t");
			for(int i = 0; i < blockSizeStarts.length - 1; i++) {
				buf.append(blockSizeStarts[i][0]).append(",");
			}
			buf.append(blockSizeStarts[blockSizeStarts.length-1][0]).append("\t");
			
			for(int i = 0; i < blockSizeStarts.length - 1; i ++) {
				buf.append(blockSizeStarts[i][1]).append(",");
			}
			buf.append(blockSizeStarts[blockSizeStarts.length - 1][1]);
		}
		
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

	public int[][] getBlockSizeStarts() {
		return blockSizeStarts;
	}

	public void setBlockSizeStarts(int[][] blockSizeStarts) {
		this.blockSizeStarts = blockSizeStarts;
	}
	
	public List<GenomicAnnotation> getBlocks() {
		List<GenomicAnnotation> blocks = new ArrayList<GenomicAnnotation>();
		if(blockSizeStarts != null) {
			for(int i = 0; i < blockSizeStarts.length; i++) {
				BasicGenomicAnnotation block = new BasicGenomicAnnotation(getName() + "_" + i);
				block.setChromosome(getChromosome());
				block.setStart(getStart() + blockSizeStarts[i][1]);
				block.setEnd(block.getStart() + blockSizeStarts[i][0]);
				blocks.add(block);
			}
		}
		return blocks;
	}

	public String getRgb() {
		return rgb;
	}

	public void setRgb(String rgb) {
		this.rgb = rgb;
	}

	public int getThickEnd() {
		return thickEnd;
	}

	public void setThickEnd(int thickEnd) {
		this.thickEnd = thickEnd;
	}

	public int getThickStart() {
		return thickStart;
	}

	public void setThickStart(int thickStart) {
		this.thickStart = thickStart;
	}
	
	public void addBlock(String name, int start, int end) {
		int adjustedStart = start - getStart();
		int length = end -start;
		System.err.print("Mapped region: " + getLocationString() + "(" + getOrientation() + ") adding " + name + " " + start + "-" + end);
		if(blockSizeStarts == null || blockSizeStarts.length == 0) {
			blockSizeStarts = new int[1][2];
			blockSizeStarts[0][0] = length;
			blockSizeStarts[0][1] = adjustedStart;
		} else {
			int [][] newBlockSizes = new int[blockSizeStarts.length + 1][2];
			boolean newAdded = false;
			int i = 0;
			int j = 0;
			while(i < blockSizeStarts.length ){
				if(blockSizeStarts[i][1] <= adjustedStart || newAdded) {
					newBlockSizes[j][0] = blockSizeStarts[i][0];
					newBlockSizes[j][1] = blockSizeStarts[i][1];
					i++;
					j++;
				} else if (blockSizeStarts[i][1]  > adjustedStart && !newAdded) {
					newBlockSizes[j][0] = length;
					newBlockSizes[j][1] = adjustedStart;	
					newAdded = true;
					System.err.print(" added at idx " + j);
					j++;
				}
			}
			if(!newAdded) {
				newBlockSizes[blockSizeStarts.length][0] = length;
				newBlockSizes[blockSizeStarts.length][1] = adjustedStart;
			} 
			blockSizeStarts  = newBlockSizes;
		}
		System.err.println(" first block start " + blockSizeStarts[0][1]);
	}
	
	public String toBED() { 
		return toBED(0, 0, 0);
	}
	
	@Override
	public String toBED(int r, int g, int b){
		if(r < 0 || r > 255 || g < 0 || g > 255 || b < 0 || b > 255) {
			throw new IllegalArgumentException("RGB values must be between 0 and 255");
		}
		String rgb = r + "," + g + "," + b;
		List<? extends Annotation> exons = getBlocks();
		//String rtrn=this.getChr()+"\t"+this.getStart()+"\t"+this.getEnd()+"\t"+(name == null ? toUCSC() : this.name)+"\t"+getCountScore()+"\t"+orientation+"\t"+start+"\t"+stop+"\t"+rgb+"\t"+exons.size();
		String rtrn=getReferenceName()+"\t"+getStart()+"\t"+getEnd()+"\t"+(getName() == null ? toUCSC() : getName())+"\t"+getScore()+"\t"+getOrientation()+"\t"+getThickStart()+"\t"+getThickEnd()+"\t"+rgb+"\t"+exons.size();
		String sizes="";
		String starts="";
		for(Annotation exon : exons){
			sizes=sizes+(exon.length())+",";
			starts=starts+(exon.getStart()-getStart())+",";
		}
		rtrn=rtrn+"\t"+sizes+"\t"+starts;
		return rtrn;
	}
	

	public boolean mayHaveBlocks() {
		return true;
	}

	//Moran : Feb 19th 2010
	//checks if the blocks of this object intersect a GenomicAnnotation
	public boolean IntersectBlocks( GenomicAnnotation annot ) {
		
		List<GenomicAnnotation> blocks=getBlocks();
		Iterator <GenomicAnnotation> it = blocks.iterator();
		while(it.hasNext()) {
			GenomicAnnotation currBlock=it.next();
			if (currBlock.overlaps(annot))
			{return true;}
			
		}
		
		return false;
	}

}
