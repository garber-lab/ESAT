package umms.core.exception;

/**
 * This is an exception class thrown when a given method is declared but not completely implemented
 * @author skadri
 *
 */
public class IncompleteMethodImplementationException extends Exception{

	 
    public IncompleteMethodImplementationException() {
    	
    	System.err.println("Method is not implemented");
    }

    /**
     * This constructor prints the specified message while throwing the exception
     * @param message
     */
    public IncompleteMethodImplementationException(String message)
    {
       super(message);
    }

}
