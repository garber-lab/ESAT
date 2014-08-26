package broad.core.math;


import java.io.IOException;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collection;


//--------------------------------------
// Systematically generate combinations.
//--------------------------------------

// From http://www.merriampark.com/comb.htm#Source
//
// This version modified by JTE to use int instead of BigInteger
//   merely for readability purposes.
//   Refer to the URL above to obtain the original source which the
//   author offers "free for you to use in whatever way you wish"
//

public class CombinationGenerator {

  private int[] a;
  private int n;
  private int r;
  private BigInteger  numLeft;
  private BigInteger  total;
 // private int[][] allPossible;
  
  //------------
  // Constructor
  //------------

  public CombinationGenerator (int n, int r) {
    if (r > n) {
      throw new IllegalArgumentException ();
    }
    if (n < 1) {
      throw new IllegalArgumentException ();
    }
    this.n = n;
    this.r = r;
    a = new int[r];
    BigInteger nFact = getFactorial(n);
    BigInteger rFact = getFactorial(r);
    BigInteger nminusrFact = getFactorial(n - r);
    total = nFact.divide(rFact.multiply(nminusrFact));
    reset();
    //this.allPossible=this.getAllPossible();
  }

  //------
 // Reset
 //------

 public void reset () {
   for (int i = 0; i < a.length; i++) {
     a[i] = i;
   }
   numLeft = new BigInteger (total.toString ());
 }

 //------------------------------------------------
 // Return number of combinations not yet generated
 //------------------------------------------------

 public BigInteger getNumLeft () {
   return numLeft;
 }

 //-----------------------------
 // Are there more combinations?
 //-----------------------------

 public boolean hasMore () {
   return numLeft.compareTo (BigInteger.ZERO) == 1;
 }

 //------------------------------------
 // Return total number of combinations
 //------------------------------------

 public int getTotal () {
   return total.intValue();
 }

 //------------------
 // Compute factorial
 //------------------

 private static BigInteger getFactorial (int n) {
   BigInteger fact = BigInteger.ONE;
   for (int i = n; i > 1; i--) {
     fact = fact.multiply (new BigInteger (Integer.toString (i)));
   }
   return fact;
 }

 //--------------------------------------------------------
 // Generate next combination (algorithm from Rosen p. 286)
 //--------------------------------------------------------

 public int[] getNext () {

	 //TODO: Fix this, why are they all the same???
	 
   if (numLeft.equals (total)) {
     numLeft = numLeft.subtract (BigInteger.ONE);
     //System.err.println("Am I always here??");
     return a;
   }

   //System.err.println("Or here?");
   int i = r - 1;
   while (a[i] == n - r + i) {
     i--;
   }
   a[i] = a[i] + 1;
   for (int j = i + 1; j < r; j++) {
     a[j] = a[i] + j - i;
   }

   numLeft = numLeft.subtract (BigInteger.ONE);
   return a;

 }
 
 private int[][] getAllPossible(){
	 int[][] rtrn=new int[this.getTotal()][];
	 
	 int counter=0;
	 while(hasMore()){
		 int[] next=getNext();
		 rtrn[counter++]=(next);
	 }
	 
	 reset();
	 
	 return rtrn;
 }
 
 public static void main(String[] args)throws IOException{
	 CombinationGenerator comb=new CombinationGenerator(12, 2);
	 while(comb.hasMore()){
		 int[] next=comb.getNext();
		 for(int i=0; i<next.length; i++){System.err.print(next[i]+"\t");}
		 System.err.println();
	 }
 }

 	public int[] getNextRandom() {
 		ArrayList<Integer> all=new ArrayList<Integer>();
 		for(int i=0; i<n; i++){all.add(i);}
 		
 		int[] rtrn=new int[r];
 		for(int i=0; i<r; i++){
 			rtrn[i]=all.remove(new Double(Math.random()*all.size()).intValue());
 		}
 		return rtrn;
 	}
 
}
