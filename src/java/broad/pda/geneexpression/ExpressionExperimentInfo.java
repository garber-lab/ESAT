package broad.pda.geneexpression;

import java.io.File;

public class ExpressionExperimentInfo {
	
	String barcodeName;
	String experimentName;
	String plateName;
	String platePosition;
	String hybDate;
	boolean passedQC;
	String sampleType;
	int timePoint;
	String batch;
	
	public ExpressionExperimentInfo(String line){
		String[] tokens=line.split("\t");
		this.barcodeName=tokens[0];
		this.experimentName=tokens[1];
		this.plateName=tokens[2];
		this.platePosition=tokens[3];
		this.hybDate=tokens[4];
		this.passedQC=(tokens[5].equalsIgnoreCase("no"));
		this.sampleType=tokens[6];
		this.timePoint=new Integer(tokens[7]);
		this.batch=tokens[9];
	}
	
	public String getBatch(){return this.batch;}
	public String getExperimentBarcode(){return this.barcodeName;}
	public String getPlateName(){return this.plateName;}
	public String getExperimentPlatePosition(){return this.platePosition;}
	public String getExperimentName(){return this.experimentName;}
	public String getExperimentHybDate(){return this.hybDate;}
	public boolean passedQC(){return this.passedQC;}
	public String getSampleType(){return this.sampleType;}
	public int getTimePoint(){return this.timePoint;}

	public boolean isControl() {
		if(this.sampleType.equalsIgnoreCase("Control")){return true;}
		return false;
	}
	
}
