package broad.core.primer3;

import java.util.Collection;

/*
* $Id: Primer3PrimerPickingResult.java 24844 2006-03-29 23:23:36Z mgarber $
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
* $Header$
*/

/**
 * TODO: document me!
 */
public interface Primer3PrimerPickingResult {

    public Collection getPrimers();

    public Primer3PickingSummary getPrimerPickingSummary();
}
