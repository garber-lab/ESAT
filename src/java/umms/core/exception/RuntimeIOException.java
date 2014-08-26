package umms.core.exception;

import java.io.IOException;

public class RuntimeIOException extends RuntimeException {
	 
    public RuntimeIOException() {
    	super();
    }
    
    public RuntimeIOException(IOException e) {
    	this(e.getMessage());
    }

    /**
     * This constructor prints the specified message while throwing the exception
     * @param message
     */
    public RuntimeIOException(String message)
    {
       super(message);
    }
}
