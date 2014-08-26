package broad.core.primer3;

/*
* $Id: Primer3PickingSummary.java 24844 2006-03-29 23:23:36Z mgarber $
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
 * A Primer3PickingSummary is a summary of how many
 * primers primer3 picked in a primer picking session.
 * It includes detailed summary information about
 * why primers were rejected
 * <br><br>
 * It's really a java wrapper around the primer3
 * oligo_struts typedef struct.
 * @author Andrew R. Zimmer
 */
public interface Primer3PickingSummary {

    /* Total number of tested oligos of given type   */
    public int getNumOligosTested();

    /* Number of oligos rejected because of Ns       */
    public int getNumNs();

    public int getNumOverlappingTargets();

    public int getNumOverlappingExcludedRegions();

    /* number of primers that failed due to gc content */
    public int getNumGC();

    /* number of primers that failed because of the incorrect number of GC's at 3' end*/
    public int getNumGCClampFailures();

    /* number of primer failures due to a low melting temp*/
    public int getNumColdMeltingTemp();

    /* number of primer failures due to a high melting temp*/
    public int getNumHotMeltingTemp();

    /* number of primer failures due to self complentarity being too high anywhere in the primer */
    public int getNumSelfComplementarityAnywhere();

    /* number of primer failures due to self complentarity being too high at the 3' end */
    public int getNum3PrimeComplementarity();

    /* number of primer failures due to similarity with mispriming library */
    public int getNumMispriming();

    /* number of primer failures due to long mononucleotide runs*/
    public int getNumPolyX();

    /* number of primer failures due to low sequence quality */
    public int getNumSeqQuality();

    /* number of primer failures due to stability of 3' end */
    public int getNum3PrimerStability();

    /* number of acceptable primers */
    public int getNumOk();
}
