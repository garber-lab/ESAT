package broad.core.math;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.StringTokenizer;
import java.util.Vector;

public class Loess {
	 private int arraySize=25000;
	 
	 private HashMap hash1, hash2, hashM, hashA, reverseAMap;
	 private double smoothingParameter=.33333333;
	 private float band=.2f;
	 private double maxIter=3;
	 private String outFile;
	 private  ArrayList M;
	 private ArrayList A;
	 private File file;
	 private int numberInNeighborhood;
	 private double pointOfEstimation;
	 private double delta;
	 private double scaledDistance;
	 private Map MapM;
	 private Map finalMValues;
	 
	 public Loess(String name)throws IOException{
	  this.file=new File(name);
	  finalMValues=new HashMap();
	  hash1= new HashMap();
	  hash2= new HashMap();
	  hashM= new HashMap();
	  hashA= new HashMap();
	  reverseAMap=new HashMap();
	  float band=.2f;
	  Vector temp=Loess.readDstFile(file);
	  StringTokenizer token;
	  String key;
	  for(int i=0; i<temp.size(); i++){
	   token= new StringTokenizer((String)temp.get(i), "\t");
	   key=token.nextToken();
	   hash1.put(key,new Double(token.nextToken()));
	   hash2.put(key,new Double(token.nextToken()));
	  }
	  this.arraySize=hash1.size()*2;
	  calculateMA();
	  MapM=Lowess();
	  
	 }
	 
	 public Loess(HashMap hash1, HashMap hash2){
	 	finalMValues=new HashMap();
	 	hashM= new HashMap();
	 	hashA= new HashMap();
	 	reverseAMap=new HashMap();
	 	float band=.2f;
	 	this.hash1=hash1;
	 	this.hash2=hash2;
	 	this.arraySize=hash1.size()*2;
	 	calculateMA();
	 	MapM=Lowess();
	 }
	 
	 public Loess(HashMap hash1, HashMap hash2, double bandwidth, double smoothing){
	 	finalMValues=new HashMap();
	 	hashM= new HashMap();
	 	hashA= new HashMap();
	 	reverseAMap=new HashMap();
	 	band=new Double(bandwidth).floatValue();
	 	this.smoothingParameter=smoothing;
	 	this.hash1=hash1;
	 	this.hash2=hash2;
	 	this.arraySize=hash1.size()*2;
	 	calculateMA();
	 	MapM=Lowess();
	 }
	 
	 public Loess(HashMap hashA, HashMap hashM, boolean flank){
	 	finalMValues=new HashMap();
	 	this.hashA=hashA;
	 	this.hashM=hashM;
	 	reverseAMap=new HashMap();
	 	float band=.2f;
	 	this.arraySize=hashM.size()*2;
	 	Iterator iter=hashA.keySet().iterator();
	 	while(iter.hasNext()){
	 		Object key=iter.next();
	 		Object val=this.hashA.get(key);
	 		ArrayList list=new ArrayList();
	 		list.add(val);
	 		reverseAMap.put(list, key);
	 	}
	 	MapM=Lowess();
	 }

	 
	 private void calculateMA(){
	  Set keySet=hash1.keySet();
	  Iterator keyIter=keySet.iterator();
	  while (keyIter.hasNext()){
	   Object key=keyIter.next();
	   Double aa=(Double)hash1.get(key);
	   Double bb=(Double)hash2.get(key);
	   //double MValue=Math.log(aa.doubleValue()/bb.doubleValue());
	   double MValue=log2(aa.doubleValue()/bb.doubleValue());
	   hashM.put(key, new Double(MValue));
	   double AValue=(log2(aa.doubleValue())+log2(bb.doubleValue()))/2;
	   //double AValue=(Math.log(aa.doubleValue())+Math.log(bb.doubleValue()))/2;
	   hashA.put(key,new Double(AValue));
	   if(reverseAMap.containsKey(new Double(AValue))){
	    //System.out.println("has duplicate A values");
	    ArrayList temp=(ArrayList)reverseAMap.get(new Double(AValue));
	    temp.add(key);
	    reverseAMap.put(new Double(AValue), temp);
	   }else{
	    ArrayList temp= new ArrayList();
	    temp.add(key);
	    reverseAMap.put(new Double(AValue), temp);}
	  }
	 }
	 
	     
	 public double log2(double val){return Math.log(val)/Math.log(2);}    
	 
	     private Map Lowess(){
	       int counter=0;
	       int start=0;
	       double pointOfEstimation;
	       double regressionMedian=0;
	       double[] intensityArray;
	       Object[] keyArray;
	       int maxDistance;
	       //double[] neighborhood=new double[this.numberInNeighborhood];
	       double[] neighborhood=new double[arraySize];
	       Object[] neighborhood_M= new Object[arraySize];
	       double[] neighborhood_regression= new double[arraySize];
	       double[] distance= new double[arraySize];
	       double[] regression_value_storage= new double[arraySize];
	       double[] w=new double[arraySize];
	       Map normalizedM=new HashMap();
	       double scaled_distance;
	       double slope;
	       double intercept;
	       double[] regression_value;
	       double[] regression_abs_value;
	       
	       
	       
	       //sort by A Values
	       Collection collec=hashA.values();
	       Object[] AValues=collec.toArray();
	       Arrays.sort(AValues);
	       intensityArray=new double[arraySize];
	       keyArray= new Object[arraySize];
	       regression_abs_value= new double[arraySize];
	       regression_value= new double[arraySize];
	       
	       for(int i=0; i<AValues.length; i++){
	        Double temp=(Double)AValues[i];
	        intensityArray[i]=temp.doubleValue();
	       }       
	       
	       for(int i=0; i<AValues.length; ){
	        ArrayList list=(ArrayList)reverseAMap.get(AValues[i]);
	        for(int j=0; j<list.size(); j++){
	         keyArray[i]=list.get(j);
	         i++;
	        }
	       }
	       
	       
	       numberInNeighborhood=new Double(this.smoothingParameter*AValues.length).intValue();
	       //System.out.println(AValues.length+" "+this.numberInNeighborhood);
	       
	       //double[] neighborhood=new double[this.numberInNeighborhood+ AValues.length];
	       //Object[] neighborhood_M= new Object[this.numberInNeighborhood+ AValues.length];
	       //double[] neighborhood_regression= new double[this.numberInNeighborhood+ AValues.length];
	       //double[] distance= new double[this.numberInNeighborhood];
	       //double[] regression_value_storage= new double[this.numberInNeighborhood];
	       //double[] w=new double[this.numberInNeighborhood];
	       regression_value_storage= new double[arraySize];
	       
	       for(int iter=0; iter<this.maxIter; iter++){
	        start=0;
	        for(int j=0; j<AValues.length; j++){
	         //if(j%100==0){System.out.println("j="+j);}
	         pointOfEstimation=intensityArray[j];
	         
	         if((j>=this.numberInNeighborhood/2) && (j<AValues.length-(this.numberInNeighborhood/2))){start++;}
	         maxDistance=0;
	         for(int i=0; i<this.numberInNeighborhood; i++){
	          
	          //neighborhood[i]=intensityArray[i];
	          neighborhood[i]=intensityArray[i+start];
	          //neighborhood_M[i]=this.hashM.get(keyArray[i]);
	          neighborhood_M[i]=this.hashM.get(keyArray[i+start]);
	          //neighborhood_regression[i]=regression_value_storage[i];
	          neighborhood_regression[i]=regression_value_storage[i+start];
	          
	          if(pointOfEstimation>neighborhood[i]){distance[i]=pointOfEstimation-neighborhood[i];}
	          else{distance[i]=neighborhood[i]-pointOfEstimation;}
	          
	          if(distance[i]>maxDistance){maxDistance=new Double(distance[i]+1).intValue();}
	          //System.out.println("maxDist "+maxDistance+" "+distance[i]);
	         }
	         
	         for(int i=0; i<this.numberInNeighborhood; i++){
	          scaledDistance=distance[i]/maxDistance;
	          if(iter==1){delta=1;}
	          else{
	           if(regressionMedian>0){delta=K(neighborhood_regression[i]/(6*regressionMedian));}
	           else{delta=0;}
	          }
	          w[i]=delta*Math.pow((1-(Math.pow(scaledDistance,3))),3);
	          //if(w[i]!=0){System.out.println(w[i]+" "+delta+" "+Math.pow((1-(Math.pow(scaledDistance,3))),3)+" "+scaledDistance+" "+distance[i]+" "+maxDistance);}
	         }
	         double a=0;
	         double b=0;
	         double c=0;
	         double d=0;
	         double x=0;
	         double y=0;
	         
	         
	         for(int i=0; i<this.numberInNeighborhood; i++){
	          a=a+w[i]*Math.pow(neighborhood[i],2);
	          b=b+w[i]*neighborhood[i];
	          c=c+w[i]*neighborhood[i];
	          d=d+w[i];
	          if(neighborhood_M[i]!=null){
	          x=x+w[i]*neighborhood[i]*((Double)neighborhood_M[i]).doubleValue();
	          y=y+w[i]*((Double)neighborhood_M[i]).doubleValue();}
	          else{x=x+w[i]*neighborhood[i]; y=y+w[i];}
	          //System.out.println(neighborhood[i]+" "+w[i]+" a "+a+" b "+b+" c "+c+" d "+d+" x "+x+" y "+y);
	          if(w[i]!=0){
	            //System.out.println("not 0");
	            //System.out.println(neighborhood[i]+" "+w[i]+" a "+a+" b "+b+" c "+c+" d "+d+" x "+x+" y "+y);
	          }
	         }
	         if((b*c-a*d)!=0){
	          slope=(b*y-d*x)/(b*c-a*d);
	          intercept=(c*x-a*y)/(b*c-a*d);
	          regression_value[j]=slope*pointOfEstimation+intercept;
	          //System.out.println("slope: "+slope+" pointEstimation "+pointOfEstimation+" intercept "+intercept+" regressionVal "+regression_value[j]);
	          //System.out.println("b*c-a*d "+(b*c-a*d));
	         }
	         else{regression_value[j]=regression_value[j];}
	        }
	        for(int j=0; j<AValues.length; j++){
	         regression_value_storage[j]=regression_value[j];
	         regression_abs_value[j]=Math.abs(regression_value[j]);
	        }
	        regressionMedian=computeMedian(regression_abs_value);
	        //System.out.println("regression median "+regressionMedian);
	       }
	       
	       
	       for(int j=0; j<AValues.length; j++){
	        Object key=keyArray[j];
	        Double MVal=(Double)this.hashM.get(key);
	        Double AVal=(Double)this.hashA.get(key);
	        Double newVal;
	        if(regression_value[j]!=Double.NaN){
	         newVal=new Double(MVal.doubleValue()-regression_value[j]);}
	        else{newVal=MVal;}
	        //System.out.println("final regression val "+regression_value[j]+" new MVal "+newVal.toString()+ "old MVal "+MVal.toString());
	        //if(newVal.isNaN()){normalizedM.put(key, new MA(AVal,MVal));
	       // finalMValues.put(key, MVal);}
	       //else{
	        normalizedM.put(key, new MA(AVal,newVal));
	        finalMValues.put(key, newVal);
	       }
	       return normalizedM;
	      }
	     
	     private double K(double med){
	       double u=med;
	       double rtrn;
	       
	       if(Math.abs(u)<=1){rtrn=(15/16)*Math.pow((1-u*u),2);}
	       else {rtrn=0;}
	       
	       return rtrn;
	     }
	     
	     private double computeMedian(double[] arr){
	       double[] vector=arr;
	       int len=vector.length;
	       double[] temp1=new double[len];
	       double[] temp2=new double[len];
	       int ind=0;
	       double m;
	       
	       for(int ijk=0; ijk<vector.length; ijk++){
	        //try{temp1[ind]=vector[ijk];ind++;}catch(ArrayIndexOutOfBoundsException e){}
	         temp1[ind]=vector[ijk];
	       }
	       int aaa=temp1.length;
	       Arrays.sort(temp1);
	       aaa=temp1.length;
	       
	       if(aaa%2==1){m=temp1[(aaa-1)/2];}
	       else{m=(temp1[(aaa-2)/2]+temp1[(aaa/2)])/2;}
	       return m;
	      }
	     
	     public static Vector readDstFile(File dstFile) throws IOException{
	   FileInputStream fileInput;
	   BufferedReader buf = null;
	   Vector temp = new Vector(2000, 500);
	   String aLine = new String("");
	      try{
	        fileInput = new FileInputStream(dstFile);
	        buf = new BufferedReader(new InputStreamReader (fileInput));
	      } catch (IOException ex){
	        return temp;
	      }

	      while(true) {
	        try {
	          aLine = buf.readLine();
	          if(aLine == null) break;
	          temp.add(aLine);
	        } catch (IOException e) {
	          temp.removeAllElements();
	          return temp;
	        }
	      }
	     return temp;
	    }
	 
	     public static void main(String[] args) throws IOException{
	       FileWriter writer= new FileWriter("/Users/mguttman/Desktop/35260-JAVALOESS.txt");
	       Loess loess= new Loess("/Users/mguttman/Desktop/35260.txt");
	       Map map=loess.getMap();
	       Set keys=map.keySet();
	       Iterator iter=keys.iterator();
	       while(iter.hasNext()){
	        Object key=iter.next();
	        Object val=map.get(key);
	        Double logVal=new Double(loess.log2(((Double)val).doubleValue()));
	        //if(new Double(val.toString()).doubleValue()>4|| new Double(val.toString()).doubleValue()<-4){val=new Double(0);}
	        Object AVal=loess.hashA.get(key);
	        Object MVal=loess.hashM.get(key);
	        writer.write(key.toString()+"\t"+AVal.toString()+"\t"+val.toString()+"\t"+MVal.toString()+"\n");
	       }
	       
	       writer.flush();
	       //for(int i=0; i<arrayX.length; i++){//System.out.println("x value "+arrayX[i]);}
	       //for(int i=0; i<arrayX.length; i++){//System.out.println("y value "+arrayYEst[i]);}
	     }
	     public Map getMap(){return MapM;}
	     
	     public Map getFinalMValues(){return finalMValues;}
	     
	     public Map getOriginalMValues(){return this.hashM;}
	     
	     public Map getAMap(){return this.hashA;}
	     
	     
	     class MA{
	     	Double AVal;
	     	Double MVal;
	     	
	     	MA(Double A, Double M){
	     		AVal=A;
	     		MVal=M;
	     	}
	     	
	     	public String toString(){return AVal+"\t"+MVal;}
	     }
}
