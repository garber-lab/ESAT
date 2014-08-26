package broad.core.primer3;

/*
* $
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2003 by the
* Broad Institute/Massachusetts Institute of Technology.
* All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support
* whatsoever. Neither the Broad Institute nor MIT can be responsible for its
* use, misuse, or functionality.
* 
* $
*/

/**
 * Exception thrown when there is a problem generating
 * primers using primer3.
 * <br><br>
 * @author Andrew R. Zimmer
 */
public class Primer3Exception extends Exception {
	private static final long serialVersionUID = -379477790476217676L;

	/**
     * Creates a new Primer3Exception
     * with the the given message.
     *
     * @param message An error message describing
     *                the nature of the problem.
     */
    public Primer3Exception(String message) {
        super(message);
    }

    /**
     * Creates a new Primer3Exception
     * with a Throwable that was the
     * underlying cause of the problem.
     *
     * @param cause the Throwable that was the
     *              underlying cause of the problem.
     */
    public Primer3Exception(Throwable cause) {
        super(cause);
    }

}