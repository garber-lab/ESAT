package broad.core.annotation;

public class MaximumContiguousSubsequence {

	
	//static private int seqStart = 0;
    //static private int seqEnd = -1;

	
	/**
     * Linear-time maximum contiguous subsequence sum algorithm.
     * seqStart and seqEnd represent the actual best sequence.
     */
    public static int[] maxSubSum3(int[] a)
    {
        int maxSum = 0;
        int thisSum = 0;
        
        int seqStart=0;
        int seqEnd=-1;

        for( int i = 0, j = 0; j < a.length; j++ )
        {
            thisSum += a[ j ];

            if( thisSum > maxSum )
            {
                maxSum = thisSum;
                seqStart = i;
                seqEnd   = j;
            }
            else if( thisSum < 0 )
            {
                i = j + 1;
                thisSum = 0;
            }
        }

        int[] rtrn={maxSum, seqStart, seqEnd};
        
        return rtrn;
    }

    
    public static double[] maxSubSum3(double[] a)
    {
        double maxSum = 0;
        double thisSum = 0;
        
        double seqStart=0;
        double seqEnd=-1;

        for( int i = 0, j = 0; j < a.length; j++ )
        {
            thisSum += a[ j ];

            if( thisSum > maxSum )
            {
                maxSum = thisSum;
                seqStart = i;
                seqEnd   = j;
            }
            else if( thisSum < 0 )
            {
                i = j + 1;
                thisSum = 0;
            }
        }

        double[] rtrn={maxSum, seqStart, seqEnd};
        
        return rtrn;
    }
    	 
    
    public static int contiguousStartSubSequenceOverMin(double [] a, double min) {
    	int i = 0;
    	while (i < a.length && a[i++] < min) { 
    		;
    	}
    	
    	return i;
    }
    
    public static int contiguousEndSubSequenceOverMin(double [] a, double min) {
    	int i = a.length;
    	while (--i >=0 && a[i] < min) { 
    		;
    	}
    	
    	return i;
    }
    
    
    public static void main( String [ ] args )
    {
        double a[ ] = {-.05,-.05,-.05,1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,1};
        double[] maxSum = maxSubSum3( a );
        System.out.println( "Max sum is " + maxSum[0] + "; it goes"
                       + " from " + maxSum[1] + " to " + maxSum[2] );
    }


	
}
