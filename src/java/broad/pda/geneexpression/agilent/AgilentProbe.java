package broad.pda.geneexpression.agilent;


public class AgilentProbe{

	int featureNumber;
	int row;
	int column;
	String accessions;
	String chrCoordinates;
	int subtypeMask;
	String subtypeName;
	int start;
	String sequence;
	int probeUID;
	int controlUID;
	String probeName;
	String geneName;
	String systematicName;
	String description;
	double positionX;
	double positionY;
	boolean geneIsFound;
	double geneProcessedSignal;
	double geneProcessedError;
	double geneMeanSignal;
	double geneMedianSignal;
	double geneBackgroundMean;
	double geneBackgroundMedian;
	boolean geneIsWellAboveBackground;
	int controlType;
	boolean probeNonUniform;
	private boolean probePopulationOutlier;
	private boolean isSaturated;
	
	
	public AgilentProbe(String line, String headerLine){
		String[] tokens=line.split("\t");
		String[] header=headerLine.split("\t");
		
		for(int i=0; i<tokens.length; i++){
			setParam(tokens[i], header[i]);
		}
	}


	private void setParam(String val, String key) {
		if(key.equalsIgnoreCase("FeatureNum")){this.featureNumber=new Integer(val);}
		else if(key.equalsIgnoreCase("Row")){this.row=new Integer(val);}
		else if(key.equalsIgnoreCase("Col")){this.column=new Integer(val);}
		else if(key.equalsIgnoreCase("accessions")){this.accessions=val;}
		else if(key.equalsIgnoreCase("chr_coord")){this.chrCoordinates=val;}
		else if(key.equalsIgnoreCase("SubTypeMask")){this.subtypeMask=new Integer(val);}
		else if(key.equalsIgnoreCase("SubTypeName")){this.subtypeName=val;}
		else if(key.equalsIgnoreCase("Start")){this.start=new Integer(val);}
		else if(key.equalsIgnoreCase("Sequence")){this.sequence=val;}
		else if(key.equalsIgnoreCase("ProbeUID")){this.probeUID=new Integer(val);}
		else if(key.equalsIgnoreCase("ControlType")){this.controlType=new Integer(val);}
		else if(key.equalsIgnoreCase("ProbeName")){this.probeName=val;}
		else if(key.equalsIgnoreCase("GeneName")){this.geneName=val;}
		else if(key.equalsIgnoreCase("SystematicName")){this.systematicName=val;}
		else if(key.equalsIgnoreCase("Description")){this.description=val;}
		else if(key.equalsIgnoreCase("PositionX")){this.positionX=new Double(val);}
		else if(key.equalsIgnoreCase("PositionY")){this.positionY=new Double(val);}
		else if(key.equalsIgnoreCase("gSurrogateUsed")){}
		else if(key.equalsIgnoreCase("gIsFound")){
			if(val.equalsIgnoreCase("1")){this.geneIsFound=true;}
			else{this.geneIsFound=false;}
		}
		else if(key.equalsIgnoreCase("gProcessedSignal")){this.geneProcessedSignal=new Double(val);}
		else if(key.equalsIgnoreCase("gProcessedSigError")){this.geneProcessedError=new Double(val);}
		else if(key.equalsIgnoreCase("gNumPixOLHi")){}
		else if(key.equalsIgnoreCase("gNumPixOLLo")){}
		else if(key.equalsIgnoreCase("gNumPix")){}
		else if(key.equalsIgnoreCase("gMeanSignal")){this.geneMeanSignal=new Double(val);}
		else if(key.equalsIgnoreCase("gMedianSignal")){this.geneMedianSignal=new Double(val);}
		else if(key.equalsIgnoreCase("gPixSDev")){}
		else if(key.equalsIgnoreCase("gPixNormIQR")){}
		else if(key.equalsIgnoreCase("gBGNumPix")){}
		else if(key.equalsIgnoreCase("gBGMeanSignal")){this.geneBackgroundMean=new Double(val);}
		else if(key.equalsIgnoreCase("gBGMedianSignal")){this.geneBackgroundMedian=new Double(val);}
		else if(key.equalsIgnoreCase("gBGPixSDev")){}
		else if(key.equalsIgnoreCase("gBGPixNormIQR")){}
		else if(key.equalsIgnoreCase("gNumSatPix")){}
		else if(key.equalsIgnoreCase("gIsSaturated")){
			if(val.equalsIgnoreCase("1")){this.isSaturated=true;}
			else{this.isSaturated=false;}
		}
		else if(key.equalsIgnoreCase("gIsFeatNonUnifOL")){
			if(val.equalsIgnoreCase("1")){this.probeNonUniform=true;}
			else{this.probeNonUniform=false;}
		}
		else if(key.equalsIgnoreCase("gIsBGNonUnifOL")){}
		else if(key.equalsIgnoreCase("gIsFeatPopnOL")){
			if(val.equalsIgnoreCase("1")){this.probePopulationOutlier=true;}
			else{this.probePopulationOutlier=false;}
		}
		else if(key.equalsIgnoreCase("gIsBGPopnOL")){}
		else if(key.equalsIgnoreCase("IsManualFlag")){}
		else if(key.equalsIgnoreCase("gBGSubSignal")){}
		else if(key.equalsIgnoreCase("gBGSubSigError")){}
		else if(key.equalsIgnoreCase("gIsPosAndSignif")){}
		else if(key.equalsIgnoreCase("gPValFeatEqBG")){}
		else if(key.equalsIgnoreCase("gNumBGUsed")){}
		else if(key.equalsIgnoreCase("gIsWellAboveBG")){
			if(val.equalsIgnoreCase("1")){this.geneIsWellAboveBackground=true;}
			else{this.geneIsWellAboveBackground=false;}
		}
		else if(key.equalsIgnoreCase("gBGUsed")){}
		else if(key.equalsIgnoreCase("gBGSDUsed")){}
		else if(key.equalsIgnoreCase("ErrorModel")){}
		else if(key.equalsIgnoreCase("gSpatialDetrendIsInFilteredSet")){}
		else if(key.equalsIgnoreCase("gSpatialDetrendSurfaceValue")){}
		else if(key.equalsIgnoreCase("SpotExtentX")){}
		else if(key.equalsIgnoreCase("SpotExtentY")){}
		else if(key.equalsIgnoreCase("gNetSignal")){}
		else if(key.equalsIgnoreCase("gMultDetrendSignal")){}
		else if(key.equalsIgnoreCase("gProcessedBackground")){}
		else if(key.equalsIgnoreCase("gProcessedBkngError")){}
		else if(key.equalsIgnoreCase("IsUsedBGAdjust")){}
		else if(key.equalsIgnoreCase("gInterpolatedNegCtrlSub")){}
		else if(key.equalsIgnoreCase("gIsInNegCtrlRange")){}
		else if(key.equalsIgnoreCase("gIsUsedInMD")){}		
	}

	public String getGeneName() {
		return this.geneName;
	}

	public String getProbename() {
		return this.probeName;
	}

	public double getProcessedSignal() {
		return this.geneProcessedSignal;
	}

	public int getColumnNum() {
		return this.column;
	}


	public String getGeneDescription() {
		return this.description;
	}

	public int getRownum() {
		return this.row;
	}

	public boolean isControlProbe() {
		boolean control=(this.controlType!=0);
		return control;
	}

	public boolean isProbeFound() {
		return this.geneIsFound;
	}

	public boolean isProbeNonUniform() {
		return this.probeNonUniform;
	}

	public boolean isProbePopulationOutlier() {
		return this.probePopulationOutlier;
	}

	public boolean isProbeWellAboveBackground() {
		return this.geneIsWellAboveBackground;
	}

	public boolean isSaturated() {
		return this.isSaturated;
	}


	public double getNormalizedSignal(double normFactor, double spikeInLowerLimit, double saturationLevel) {
		double floor=Math.pow(10, spikeInLowerLimit);
		double val=geneProcessedSignal;
		if(!this.geneIsWellAboveBackground || val<floor){
			return -999;
			//val=floor;
			//return Double.NaN;
		}
		double rtrn=Math.log(val/normFactor)/Math.log(2);
		return rtrn;
	}
	
	public double getNormalizedSignal(double spikeInLowerLimit, double saturationLevel) {
		double floor=Math.pow(10, spikeInLowerLimit);
		double val=geneProcessedSignal;
		if(!this.geneIsWellAboveBackground || val<floor){
			return -999;
			//val=floor;
			//return Double.NaN;
		}
		double rtrn=Math.log(val)/Math.log(2);
		return rtrn;
	}


	public boolean isFlagged(double spikeInLowerLimit) {
		double floor=Math.pow(10, spikeInLowerLimit);
		double val=geneProcessedSignal;
		if(!this.geneIsFound || !this.geneIsWellAboveBackground || this.isProbeNonUniform() || this.isProbePopulationOutlier() || val<floor){return true;}
		return false;
	}


	public double getLogProcessedSignal() {
		return Math.log(geneProcessedSignal);
	}
	
}
