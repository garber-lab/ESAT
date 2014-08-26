package umms.core.alignment;

/**
 * This is an exception class for the PairedEndAlignment classes. 
 * The ChromosomeInconsistencyException is thrown whenever the chromosomes of the 
 * individual pair reads do not match
 * @author skadri
 *
 */
public class ChromosomeInconsistencyException extends Exception{

	 
    public ChromosomeInconsistencyException() {
    	
    }

    /**
     * This constructor prints the specified message while throwing the exception
     * @param message
     */
    public ChromosomeInconsistencyException(String message)
    {
       super(message);
    }

}
