package broad.pda.geneexpression.agilent;


public class AgilentArrayStats{

	double normalizationFactor;	
	double saturationValue;
	int numberOfSaturatedFeatures;
	int numberOfStatuaredNonControls;
	double percentile99;
	double percentile50;
	double percentile1;
	int numberOfNonUniformFeatures;
	int numberOfNonUniformGenes;
	int totalNumberOfReplicatedGenes;
	int numberOfNonUniformPopulation;
	int numberOfNonUniformFeaturesBackground;
	int numberOfNonUniformPopulationBackground;
	int totalNumFeatures;
	int numberOfFoundFeatures;
	double spikeInDetectionLimit;
	double medianPercentileCVProcessedSignal;
	int numberFeaturesWellAboveBackground;
	
	//Names match header
	double gPercentileIntensityProcessedSignal;
	double gSaturationValue;
	int gNumSatFeat;
	int gNumFeatureNonUnifOL;
	int gNumPopnOL;
	int gNumNonUnifBGOL;
	int gNumPopnBGOL;
	int gNonCtrlNumSatFeat;
	double gNonCtrl99PrcntNetSig;
	double gNonCtrl50PrcntNetSig;
	double gNonCtrl1PrcntNetSig;
	int TotalNumFeatures;
	int NumFoundFeat;
	double eQCOneColorSpikeDetectionLimit;
	double gMedPrcntCVProcSignal;
	int TotalNumberOfReplicatedGenes;
	int gNonCtrlNumWellAboveBG;
	
	
	float	gDarkOffsetAverage;
	float	gDarkOffsetMedian;
	float	gDarkOffsetStdDev;
	int	gDarkOffsetNumPts;
	float	gAvgSig2BkgeQC;
	float	gAvgSig2BkgNegCtrl;
	float	gRatioSig2BkgeQC_NegCtrl;
	float	gLocalBGInlierNetAve;
	float	gLocalBGInlierAve;
	float	gLocalBGInlierSDev;
	int	gLocalBGInlierNum;
	float	gGlobalBGInlierAve;
	float	gGlobalBGInlierSDev;
	int	gGlobalBGInlierNum;
	float	gOffsetUsed;
	float	gGlobalFeatInlierAve;
	float	gGlobalFeatInlierSDev;
	int	gGlobalFeatInlierNum;
	float	AnyColorPrcntFeatNonUnifOL;
	float	AnyColorPrcntBGNonUnifOL;
	float	AnyColorPrcntFeatPopnOL;
	float	AnyColorPrcntBGPopnOL;
	float	TotalPrcntFeatOL;
	int	gNumNegBGSubFeat;
	int	gNonCtrlNumNegFeatBGSubSig;
	float	gSpatialDetrendRMSFit;
	float	gSpatialDetrendRMSFilteredMinusFit;
	float	gSpatialDetrendSurfaceArea;
	float	gSpatialDetrendVolume;
	float	gSpatialDetrendAveFit;
	int	gCtrleQCNumSatFeat;
	float	gCtrleQC99PrcntNetSig;
	float	gCtrleQC50PrcntNetSig;
	float	gCtrleQC1PrcntNetSig;
	float	geQCMedPrcntCVBGSubSig;
	float	geQCSig2BkgLow2;
	int	gNegCtrlNumInliers;
	float	gNegCtrlAveNetSig;
	float	gNegCtrlSDevNetSig;
	float	gNegCtrlAveBGSubSig;
	float	gNegCtrlSDevBGSubSig;
	float	gAveNumPixOLLo;
	float	gAveNumPixOLHi;
	float	gPixCVofHighSignalFeat;
	int	gNumHighSignalFeat;
	float	AddErrorEstimateGreen;
	float	ROIHeight;
	float	ROIWidth;
	float	CentroidDiffX;
	float	CentroidDiffY;
	float	MaxNonUnifEdges;
	float	MaxSpotNotFoundEdges;
	float	gMultDetrendRMSFit;
	float	gMultDetrendSurfaceAverage;
	String	eQCLowSigName2;
	float	eQCOneColorLogLowSignal;
	float	eQCOneColorLogLowSignalError;
	float	eQCOneColorLogHighSignal;
	float	eQCOneColorLinFitLogLowConc;
	float	eQCOneColorLinFitLogLowSignal;
	float	eQCOneColorLinFitLogHighConc;
	float	eQCOneColorLinFitLogHighSignal;
	float	eQCOneColorLinFitSlope;
	float	eQCOneColorLinFitIntercept;
	float	eQCOneColorLinFitRSQ;
	float	gNonCtrl50PrcntBGSubSig;
	float	gCtrleQC50PrcntBGSubSig;
	float	geQCMedPrcntCVProcSignal;
	float	gOutlierFlagger_Auto_FeatB_Term;
	float	gOutlierFlagger_Auto_FeatC_Term;
	float	gOutlierFlagger_Auto_BgndB_Term;
	float	gOutlierFlagger_Auto_BgndC_Term;
	float	OutlierFlagger_FeatChiSq;
	float	OutlierFlagger_BgndChiSq;
	int	GriddingStatus;
	int	NumGeneNonUnifOL;
	int	ExtractionStatus;
	String	QCMetricResults;
	boolean	GridHasBeenOptimized;
	float	gNegCtrlSpread;
	float	Metric_AnyColorPrcntFeatNonUnifOL;
	boolean	Metric_AnyColorPrcntFeatNonUnifOL_IsInRange;
	float	Metric_DetectionLimit;
	boolean	Metric_DetectionLimit_IsInRange;
	float	Metric_absGE1E1aSlope;
	boolean	Metric_absGE1E1aSlope_IsInRange;
	float	Metric_gE1aMedCVProcSignal;
	boolean	Metric_gE1aMedCVProcSignal_IsInRange;
	float	Metric_gNegCtrlAveBGSubSig;
	boolean	Metric_gNegCtrlAveBGSubSig_IsInRange;
	float	Metric_gNegCtrlAveNetSig;
	boolean	Metric_gNegCtrlAveNetSig_IsInRange;
	float	Metric_gNegCtrlSDevBGSubSig;
	boolean	Metric_gNegCtrlSDevBGSubSig_IsInRange;
	float	Metric_gNonCntrlMedCVProcSignal;
	boolean	Metric_gNonCntrlMedCVProcSignal_IsInRange;
	float	Metric_gSpatialDetrendRMSFilteredMinusFit;
	boolean	Metric_gSpatialDetrendRMSFilteredMinusFit_IsInRange;

	
	
	public AgilentArrayStats(String line, String headerLine){
		String[] tokens=line.split("\t");
		String[] header=headerLine.split("\t");
		
		for(int i=0; i<tokens.length; i++){
			setParam(tokens[i], header[i]);
		}
	}


	private void setParam(String val, String key) {
		if(key.equalsIgnoreCase("gPercentileIntensityProcessedSignal")){this.normalizationFactor=new Float(val); this.gPercentileIntensityProcessedSignal=new Double(val);}
		else if(key.equalsIgnoreCase("gSaturationValue")){this.saturationValue=new Double(val); this.gSaturationValue=new Double(val);}
		else if(key.equalsIgnoreCase("gNumSatFeat")){this.numberOfSaturatedFeatures=new Integer(val); this.gNumSatFeat=new Integer(val);}
		else if(key.equalsIgnoreCase("gNumFeatureNonUnifOL")){this.numberOfNonUniformFeatures=new Integer(val); this.gNumFeatureNonUnifOL=new Integer(val);}
		else if(key.equalsIgnoreCase("gNumPopnOL")){this.numberOfNonUniformPopulation=new Integer(val);this.gNumPopnOL=new Integer(val);}
		else if(key.equalsIgnoreCase("gNumNonUnifBGOL")){this.numberOfNonUniformFeaturesBackground=new Integer(val); this.gNumNonUnifBGOL=new Integer(val);}
		else if(key.equalsIgnoreCase("gNumPopnBGOL")){this.numberOfNonUniformPopulationBackground=new Integer(val); this.gNumPopnBGOL=new Integer(val);}
		else if(key.equalsIgnoreCase("gNonCtrlNumSatFeat")){this.numberOfStatuaredNonControls=new Integer(val); this.gNonCtrlNumSatFeat=new Integer(val);}
		else if(key.equalsIgnoreCase("gNonCtrl99PrcntNetSig")){this.percentile99=new Double(val); this.gNonCtrl99PrcntNetSig=new Double(val);}
		else if(key.equalsIgnoreCase("gNonCtrl50PrcntNetSig")){this.percentile50=new Double(val);this.gNonCtrl50PrcntNetSig=new Double(val);}
		else if(key.equalsIgnoreCase("gNonCtrl1PrcntNetSig")){this.percentile1=new Double(val); this.gNonCtrl1PrcntNetSig=new Double(val);}
		else if(key.equalsIgnoreCase("TotalNumFeatures")){this.totalNumFeatures=new Integer(val); this.TotalNumFeatures=new Integer(val);}
		else if(key.equalsIgnoreCase("NumFoundFeat")){this.numberOfFoundFeatures=new Integer(val); this.NumFoundFeat=new Integer(val);}
		else if(key.equalsIgnoreCase("eQCOneColorSpikeDetectionLimit")){this.spikeInDetectionLimit=new Double(val); this.eQCOneColorSpikeDetectionLimit=new Double(val);}
		else if(key.equalsIgnoreCase("gMedPrcntCVProcSignal")){this.medianPercentileCVProcessedSignal=new Double(val); this.gMedPrcntCVProcSignal=new Double(val);}
		else if(key.equalsIgnoreCase("TotalNumberOfReplicatedGenes")){this.totalNumberOfReplicatedGenes=new Integer(val); this.TotalNumberOfReplicatedGenes=new Integer(val);}
		else if(key.equalsIgnoreCase("gNonCtrlNumWellAboveBG")){this.numberFeaturesWellAboveBackground=new Integer(val); this.gNonCtrlNumWellAboveBG=new Integer(val);}
		else if(key.equalsIgnoreCase("AnyColorPrcntFeatNonUnifOL")){this.AnyColorPrcntFeatNonUnifOL=new Float(val);}
		else if(key.equalsIgnoreCase("eQCOneColorSpikeDetectionLimit")){this.eQCOneColorSpikeDetectionLimit=new Double(val);}
		else if(key.equalsIgnoreCase("Metric_absGE1E1aSlope")){this.Metric_absGE1E1aSlope=new Float(val);}
		else if(key.equalsIgnoreCase("Metric_gE1aMedCVProcSignal")){this.Metric_gE1aMedCVProcSignal=new Float(val);}
		else if(key.equalsIgnoreCase("gNegCtrlAveBGSubSig")){this.gNegCtrlAveBGSubSig=new Float(val);}
		else if(key.equalsIgnoreCase("Metric_gNegCtrlAveNetSig")){this.Metric_gNegCtrlAveNetSig=new Float(val);}
		else if(key.equalsIgnoreCase("gNegCtrlSDevBGSubSig")){this.gNegCtrlSDevBGSubSig=new Float(val);}
		else if(key.equalsIgnoreCase("Metric_gNonCntrlMedCVProcSignal")){this.Metric_gNonCntrlMedCVProcSignal=new Float(val);}
		else if(key.equalsIgnoreCase("Metric_gSpatialDetrendRMSFilteredMinusFit")){this.Metric_gSpatialDetrendRMSFilteredMinusFit=new Float(val);}
	}


	
	public double getNonControl1PercentileSignal() {
		return this.percentile1;
	}


	public double getNonControl50PercentileSignal() {
		return this.percentile50;
	}

	public double getNonControl99PercentileSignal() {
		return this.percentile99;
	}

	public int getNumberOfBackgroundNonUniformFeaturesOverlimit() {
		return this.numberOfNonUniformFeaturesBackground;
	}

	public int getNumberOfBackgroundPopulationOverlimit() {
		return this.numberOfNonUniformPopulationBackground;
	}

	public int getNumberOfFoundFeatures() {
		return this.numberOfFoundFeatures;
	}

	public int getNumberOfNonUniformFeaturesOverlimit() {
		return this.numberOfNonUniformFeatures;
	}

	public int getNumberOfPopulationOverlimit() {
		return this.numberOfNonUniformPopulation;
	}

	public int getNumberOfSaturatedProbes() {
		return this.numberOfSaturatedFeatures;
	}

	public double getSaturationSignalValue() {
		return this.saturationValue;
	}

	
	public double getSpikeInDetectionLimit() {
		return this.spikeInDetectionLimit;
	}

	public int getTotalNumberOffeatures() {
		return this.totalNumFeatures;
	}


	public double getPercentFeaturesWellAboveBackground() {
		return this.numberFeaturesWellAboveBackground/(double)this.totalNumFeatures;
	}


	public double getNormalizationFactor() {
		return this.normalizationFactor;
	}

}
