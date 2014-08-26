package broad.core.annotation;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import umms.core.annotation.Annotation;

public class GFF extends BasicGenomicAnnotation {
	private String source;
	private String feature;
	private Integer    frame;
	private String group;
	
	private LinkedHashMap<String,List<String>> attributes;

	public GFF() {
		super();
		attributes = new LinkedHashMap<String, List<String>>();
	}
		
	public GFF(String name) {
		super(name);
		attributes = new LinkedHashMap<String, List<String>>();
	}
	
	public GFF(LightweightGenomicAnnotation annotation) {
		super(annotation);
		attributes = new LinkedHashMap<String, List<String>>();
		setSource("Unknown");
		setFeature("GenomicAnnotation");
	}
	public GFF(Annotation annotation) {
		super(annotation);
		attributes = new LinkedHashMap<String, List<String>>();
		setSource("Unknown");
		setFeature("GenomicAnnotation");
	}
	
	public GFF(GFF gff) {
		super(gff);
		attributes = new LinkedHashMap<String, List<String>>(gff.attributes.size());
		Iterator<String> gffAttrKeyIt = gff.attributes.keySet().iterator();
		while(gffAttrKeyIt.hasNext()) {
			String key = gffAttrKeyIt.next();
			attributes.put(key, gff.getValues(key));
		}
		
		setSource(gff.getSource());
		frame  = gff.getFrame();
		setFeature(gff.getFeature());
	}
	
	public GFF(String [] rawFields) {
		attributes = new LinkedHashMap<String, List<String>>();
		int i = 0;
		super.setChromosome(rawFields[i++]);
		setSource(rawFields[i++]);
		setFeature(rawFields[i++]);
		setStart(rawFields[i++]);
		setEnd(Integer.parseInt(rawFields[i++]));
		String scoreStr = rawFields[i++];
		if(!".".equals(scoreStr)) {
			setScore(new Float(Float.parseFloat(scoreStr)));
		}
		
		setStrand(rawFields[i++]);
		String frameStr = rawFields[i++];
		if(!".".equals(frameStr)) {
			setFrame(Integer.parseInt(frameStr));
		}
		
		setAttributes(rawFields[i]);

	}
	
	public void setAttributes(String rawAttributeData) {
		String [] attributeArray = rawAttributeData.split(" *; *");
		String firstAttributeVal = null;
		for(int i = 0; i < attributeArray.length; i++) {
			String [] attributeData = attributeArray[i].split(" +|=");
			if(i==0 && attributeData.length > 1) {
				firstAttributeVal = attributeData[1].replace("\"", "");
			}
			for(int j = 1; j < attributeData.length; j++) {
				addAttribute(attributeData[0], attributeData[j].replace("\"", ""));
			}
		}
		
		if(getFirstValue("name") != null) {
			setName(getFirstValue("name"));
		}else if (firstAttributeVal != null ) { 
			setName(firstAttributeVal);
		}else if (attributeArray.length == 1) {
			setName(attributeArray[0]);
		} else {
			setName("no_name");
		}
		
	}

	public Strand getStrand() {return super.getOrientation();}
	public void setStrand(String strand) {setReversedOrientation("-".equals(strand) || "-1".equals(strand));} 
	
	public String getFeature() {
		return feature;
	}

	public void setFeature(String feature) {
		this.feature = feature.intern();
	}

	public int getFrame() {
		return frame;
	}

	public void setFrame(int frame) {
		this.frame = frame;
	}

	public String getSource() {
		return source;
	}
	
	public void addAttribute(String attributeName, String attributeValue) {
		List<String> values = attributes.get(attributeName);
		if(values == null) {
			values = new ArrayList<String>();
			attributes.put(attributeName.intern(), values);
		}
		values.add(attributeValue);
	}
	
	public void clearAttribute(String attributeName) {
		attributes.remove(attributeName);
	}
	
	public void addAttributeValues(String attributeName, List<String> values) {
		List<String> attValues = attributes.get(attributeName);
		if(values == null) {
			attValues = new ArrayList<String>();
			attributes.put(attributeName.intern(), attValues);
		}
		values.addAll(values);		
	}
	
	public Map<String, List<String>> getAttributes() { return attributes;}
	
	public List<String> getValues(String attributeName) { return attributes.get(attributeName); }
	
	public String getFirstValue(String attributeName) {
		String val = null;
		if(attributes.containsKey(attributeName) && (attributes.get(attributeName).size() > 0)) {
			val = attributes.get(attributeName).get(0);
		}
		return val;
	}
	
	public String getSequenceName() {
		return ("chr" + super.getChromosome()).intern();
	}
	
	public void setSource(String source) {this.source = source.intern();}
	
	public String toString() {
		return toString(false, false);
	}
	
	public String toString(boolean quoteValues) {
		return toString(quoteValues, false);
	}
	
	public String toString(boolean quoteValues, boolean useEqualSign) {
		StringBuffer sb = new StringBuffer(getGFFFieldString(getChromosome()));
		
		sb.append("\t")
			.append(getGFFFieldString(getSource()))
			.append("\t")
			.append(getGFFFieldString(getFeature()))
			.append("\t")
			.append(getStart()+1) //BUG FIX : MORAN AUG 17TH, ADD +1
			.append("\t")
			.append(getEnd()) //GFFs is inclusive end //BUG FIX : MORAN AUG 17TH, REMOVE -1
			.append("\t")
			.append(getGFFFieldString(getScore()))
			.append("\t")
			.append(getStrand())
			.append("\t")
			.append(getGFFFieldString(frame));
			//.append("\t");
			//.append(getName());
		
		//sb.append(getAttributeValString("name", getName(), quoteValues));
		
		Iterator<String> attributeNameIt = attributes.keySet().iterator();
		String attributeName = null;
		if(attributeNameIt.hasNext()) {
			sb.append("\t");
			while(attributeNameIt.hasNext()) {
				attributeName = attributeNameIt.next();
				if("name".equals(attributeName)) {
					continue;
				}
				//sb.append(" ; ");
				if(!useEqualSign) {
					sb.append(getAttributeValString(attributeName, attributes.get(attributeName), quoteValues));
				} else {
					String val = attributes.get(attributeName).size() > 0 ? attributes.get(attributeName).get(0) : "0";
					sb.append(getAttributeValString(attributeName, val, quoteValues, useEqualSign));
				}
				if(attributeNameIt.hasNext()) {
					sb.append(" ; ");
				}
				
			}
		} 

		
		return sb.toString();
	}
	
	
	
	public String toCufflinksString(boolean quoteValues) {
		StringBuffer sb = new StringBuffer(getGFFFieldString(getChromosome()));
		
		sb.append("\t")
			.append(getGFFFieldString(getSource()))
			.append("\t")
			.append(getGFFFieldString(getFeature()))
			.append("\t")
			.append(getStart()+1) //BUG FIX : MORAN AUG 17TH, ADD +1
			.append("\t")
			.append(getEnd()) //GFFs is inclusive end //BUG FIX : MORAN AUG 17TH, REMOVE -1
			.append("\t")
			.append(getGFFFieldString(getScore()))
			.append("\t")
			.append(getStrand())
			.append("\t")
			.append(getGFFFieldString(frame));
			//.append("\t");
			//.append(getName());
		
		//sb.append(getAttributeValString("name", getName(), quoteValues));
		
		//gene_id and transcript_id should appear first based on UCSC genome browser
		
		sb.append("\t");
		sb.append(getAttributeValString("gene_id",attributes.get("gene_id"), quoteValues)); 
		sb.append(" ; ");
		sb.append(getAttributeValString("transcript_id", attributes.get("transcript_id"), quoteValues));
		sb.append(" ; ");
		
		Iterator<String> attributeNameIt = attributes.keySet().iterator();
		String attributeName = null;
		if(attributeNameIt.hasNext()) {
			while(attributeNameIt.hasNext()) {
				attributeName = attributeNameIt.next();
				if("gene_id".equals(attributeName) || "transcript_id".equals(attributeName) ) {
					continue;
				}
				//sb.append(" ; ");
				sb.append(getAttributeValString(attributeName, attributes.get(attributeName), quoteValues));
				if(attributeNameIt.hasNext()) {
					sb.append(" ; ");
				}
				
			}
		} 

		
		return sb.toString();
	}
	
	private String getAttributeValString(String key, String val, boolean quoteValues, boolean useEqualSign) {
		StringBuilder sb = new StringBuilder(key);
		sb.append(useEqualSign ? "=" : " ");
		if(quoteValues) {sb.append("\"");};
		sb.append(val);
		if(quoteValues) {sb.append("\"");};
		return sb.toString();
	}
	
	private String getAttributeValString(String key, Collection<String> vals, boolean quoteValues) {
		StringBuilder sb = new StringBuilder(key);
		
		for(String val : vals) {
			sb.append(" ");
			boolean alreadyQuoted = val.contains("\"");  //remove double quote
			if(quoteValues && !alreadyQuoted) {sb.append("\"");};
			sb.append(val);
			if(quoteValues && !alreadyQuoted ) {sb.append("\"");};
		}
		return sb.toString();
	}

	private String getGFFFieldString(Object obj) {
		return obj == null ? "." : obj.toString();
	}

	public String getGroup() {
		return group;
	}

	public void setGroup(String group) {
		this.group = group.intern();
	}

	public static String getGTFAttr(String rawFields, String attrName){
	    String [] attributeArray =rawFields.split(";");
	    return getGTFAttr(attributeArray, attrName);
	}
	
	public static String getGTFAttr(String [] attributeArray, String attrName){
	    String res="";
	    for (int i=0;  i < attributeArray.length; i++) {
		String [] attributeData = attributeArray[i].split(" +|=");
		if (attributeData[0].equalsIgnoreCase(attrName)){
			res=attributeData[1];
			break;
		}
		else if (attributeData[1].equalsIgnoreCase(attrName)){
			res=attributeData[2];
			break;
		}
	    }
	return res;
	}	
}
