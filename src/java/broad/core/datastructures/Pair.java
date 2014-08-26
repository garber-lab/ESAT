package broad.core.datastructures;

public class Pair<T1>{

	T1 value1;
	T1 value2;
	
	public Pair(T1 v1, T1 v2){
		this.value1=v1;
		this.value2=v2;
	}
	
	public Pair(){}

	
	public void setValue1(T1 v1){this.value1=v1;}
	public void setValue2(T1 v2){this.value2=v2;}
	
	public T1 getValue1(){return value1;}
	public T1 getValue2(){return value2;}
	
	public boolean hasValue2(){
		if(value2==null){return false;}
		return true;
	}

	public boolean hasValue1(){
		if(value1==null){return false;}
		return true;
	}
	
	public boolean isEmpty() {
		if(value1==null && value2==null){return true;}
		return false;
	}
	
	public boolean equals(Pair<T1> other){
		if(other.value1.equals(value1) && other.value2.equals(value2)){return true;}
		return false;
	}


	public boolean equals(Object other){
		Pair<T1> t=(Pair<T1>)other;
		return equals(t);
	}
	
	public int hashCode() {
		String h = Integer.valueOf(value1.hashCode()).toString() + "_" + Integer.valueOf(value2.hashCode()).toString();
		return h.hashCode();
	}
	

	/*public int compareTo(Object o) {
		return compareTo((Pair<T1>)o);
	}
	
	public int compareTo(Pair<T1> other) {
		//System.err.println("Comparing for uniqueness");
		//first check pair 1, if equal then return pair2
		int order=((Comparable)value1).compareTo(other.getValue1());
		//else return pair1
		if(order==0){return ((Comparable)value2).compareTo(other.getValue2());}
		else{return order;}
	}*/
	
}
