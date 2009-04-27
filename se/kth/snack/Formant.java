/*
 * This software has been licensed to the Centre of Speech Technology, KTH
 * by AT&T Corp. and Microsoft Corp. with the terms in the accompanying
 * file BSD.txt, which is a BSD style license.
 *
 *    "Copyright (c) 1987-1990  AT&T, Inc.
 *    "Copyright (c) 1986-1990  Entropic Speech, Inc. 
 *    "Copyright (c) 1990-1994  Entropic Research Laboratory, Inc. 
 *                   All rights reserved"
 *
 * Written by:  David Talkin
 * Revised by: John Shore
 */

package se.kth.snack;


class Formant implements Command {

    // debug levels
    static final int DEB_PAUSE = 8;
    static final int DEB_LPC_PARS = 4;
    static final int DEB_PARAMS = 2;
    static final int DEB_ENTRY = 1;
    static final int MAXFORMANTS = 7;

    /** for storing the crosscorrelation information */
    class CROSS {
        /** rms energy in the reference window */
        double rms;
        /** 1st-order autoregressive flattening constant. */
        double k1;
        /** max in the crosscorr. fun. q15 */
        double maxval;
        /** lag # at which max occured */
        short maxloc;
        /** the number of correlation lags computed */
        short nlags;
        /** the first non-zero lag computed */
        short firstlag;
        /** the normalized corsscor. fun. q15 */
        short[] correl;
    }

    /** for storing the DP information */
    class DPREC {
        /** # of candidate pitch intervals in the frame */
        short ncands;
        /** locations of the candidates */
        short[] locs;
        /** peak values of the candidates */
        short[] pvals;
        /** modified peak values of the candidates */
        double[] mpvals;
        /** pointers to best previous cands. */
        short[] prept;
        /** cumulative error for each candidate */
        double[] dpvals;
    }

    /* end of structure definitions for the pitch tracker */

    /* */

    /* Structure definitions for the formant tracker.. */

    /** structure of a DP lattice node for formant tracking */
    class FORM {
        /** # of candidate mappings for this frame */
        short ncand;
        /** pole-to-formant map-candidate array */
        short[][] cand;
        /** backpointer array for each frame */
        short[] prept;
        /** cum. errors associated with each cand. */
        double[] cumerr;
    }

    /** structure to hold raw LPC analysis data */
    class POLE {
        /** rms for current LPC analysis frame */
        double rms;
        /** rms for current F0 analysis frame */
        double rms2;
        /** fundamental frequency estimate for this frame */
        double f0;
        /** probability that frame is voiced */
        double pv;
        /** spec. distance between current and prev. frames */
        double change;
        /** # of complex poles from roots of LPC polynomial */
        short npoles;
        /** array of complex pole frequencies (Hz) */
        double[] freq;
        /** array of complex pole bandwidths (Hz) */
        double[] band;
    }

    int debug = 0;

    int w_verbose = 0;

    /* dpform.c */

    /* a formant tracker based on LPC polynomial roots and dynamic programming */

    /*
     * At each frame, the LPC poles are ordered by increasing frequency. All
     * "reasonable" mappings of the poles to F1, F2, ... are performed. The cost
     * of "connecting" each of these mappings with each of the mappings in the
     * previous frame is computed. The lowest cost connection is then chosen as
     * the optimum one. At each frame, each mapping has associated with it a
     * cost based on the formant bandwidths and frequencies. This "local" cost
     * is finally added to the cost of the best "connection." At end of
     * utterance (or after a reasonable delay like .5sec) the best mappings for
     * the entire utterance may be found by retracing back through best
     * candidate mappings, starting at end of utterance (or current frame).
     */

    // Here are the major fudge factors for tweaking the formant tracker.

    /** maximum number of candidate mappings allowed */
    static final int MAXCAN = 300;

    /** equivalent delta-Hz cost for missing formant */
    static final double MISSING = 1;
    /** equivalent bandwidth cost of a missing formant */
    static final int NOBAND = 1000;
    /** cost for proportional frequency changes */
    static final double DF_FACT = 20.0;
    /** with good "stationarity" function: */
    /** DF_FACT = 80.0, */
    /* cost for proportional frequency changes */
    /** cost for proportional dev. from nominal freqs. */
    static final double DFN_FACT = 0.3;
    /** cost per Hz of bandwidth in the poles */
    static final double BAND_FACT = .002;
    /** F_BIAS = 0.0004, bias toward selecting low-freq. poles */
    /** bias toward selecting low-freq. poles */
    static final double F_BIAS = 0.000;
    /** cost of mapping f1 and f2 to same frequency */
    static final double F_MERGE = 2000.0;

    /** "nominal" freqs. */
    private double[] fre, fnom = {
        500, 1500, 2500, 3500, 4500, 5500, 6500
    },
    /** frequency bounds */
    fmins = {
        50, 400, 1000, 2000, 2000, 3000, 3000
    },
    /** for 1st 5 formants */
    fmaxs = {
        1500, 3500, 4500, 5000, 6000, 6000, 8000
    };

    /** number of poles to consider */
    private int maxp;
    /** number of formants to find */
    int maxf;
    int ncan;
    int domerge = 1;

    private short[][] pc;

    /** can this pole be this freq.? */
    private boolean canbe(int pnumb, int fnumb) {
        return fre[pnumb] >= fmins[fnumb] && fre[pnumb] <= fmaxs[fnumb];
    }

    /**
     * This does the real work of mapping frequencies to formants.
     * 
     * @param cand candidate number being considered
     * @param pnumb pole number under consideration
     * @param fnumb formant number under consideration
     */
    private void candy(int cand, int pnumb, int fnumb) {
        int i, j;

        if (fnumb < maxf) {
            pc[cand][fnumb] = -1;
        }
        if ((pnumb < maxp) && (fnumb < maxf)) {
            // Sysem.err.printf("\ncan:%3d  pnumb:%3d  fnumb:%3d", cand, pnumb, fnumb);
            if (canbe(pnumb, fnumb)) {
                pc[cand][fnumb] = (short) pnumb;
                // allow for f1 , f2 merger
                if (domerge != 0 && (fnumb == 0) && (canbe(pnumb, fnumb + 1))) {
                    ncan++;
                    pc[ncan][0] = pc[cand][0];
                    candy(ncan, pnumb, fnumb + 1); // same pole, next formant
                }
                candy(cand, pnumb + 1, fnumb + 1); // next formant; next pole
                if (((pnumb + 1) < maxp) && canbe(pnumb + 1, fnumb)) {
                    // try other frequencies for this formant
                    ncan++; // add one to the candidate index/tally
                    // Sysem.err.printf("\n%4d  %4d  %4d", ncan, pnumb + 1, fnumb);
                    for (i = 0; i < fnumb; i++) {
                        // clone the lower formants
                        pc[ncan][i] = pc[cand][i];
                    }
                    candy(ncan, pnumb + 1, fnumb);
                }
            } else {
                candy(cand, pnumb + 1, fnumb);
            }
        }
        // If all pole frequencies have been examined without finding one which
        // will map onto the current formant, go on to the next formant leaving
        // the current formant null.
        if ((pnumb >= maxp) && (fnumb < maxf - 1) && (pc[cand][fnumb] < 0)) {
            if (fnumb != 0) {
                j = fnumb - 1;
                while ((j > 0) && pc[cand][j] < 0) {
                    j--;
                }
                i = ((j = pc[cand][j]) >= 0) ? j : 0;
            } else {
                i = 0;
            }
            candy(cand, i, fnumb + 1);
        }
    }

    /**
     * Given a set of pole frequencies and allowable formant frequencies for
     * nform formants, calculate all possible mappings of pole frequencies to
     * formants, including, possibly, mappings with missing formants.
     * 
     * poles ordered by increasing FREQUENCY
     */
    private void get_fcand(int npole, int nform, short[][] pcan, double[] freq, double[] band) {
        ncan = 0;
        pc = pcan;
        fre = freq;
        maxp = npole;
        maxf = nform;
        candy(ncan, 0, 0);
        ncan++; // (converts ncan as an index to ncan as a candidate count)
    }

    /** */
    private void set_nominal_freqs(double f1) {
        for (int i = 0; i < MAXFORMANTS; i++) {
            fnom[i] = ((i * 2) + 1) * f1;
            fmins[i] = fnom[i] - ((i + 1) * f1) + 50.0;
            fmaxs[i] = fnom[i] + (i * f1) + 1000.0;
        }
    }

    /**
     * find the maximum in the "stationarity" function (stored in rms)
     */
    private double get_stat_max(POLE[] pole, int nframes) {
        int i, poleP = 0;
        double amax, t;

        for (i = 1, amax = (pole[poleP++]).rms; i++ < nframes;) {
            if ((t = (pole[poleP++]).rms) > amax) {
                amax = t;
            }
        }

        return amax;
    }

    /** */
    private Sound dpform(Sound ps, int nform, double nom_f1) {
        double pferr, conerr, minerr, dffact, ftemp, berr, ferr, bfact, ffact, rmsmax, fbias, rmsdffact, merger = 0.0, merge_cost, FBIAS;
        double[][] fr, ba;
        int i, j, k, l, ic, ip, mincan = 0;
        short[][] pcan;
        FORM[] fl;
        POLE[] pole; // raw LPC pole data structure array
        Sound fbs;
        int dmaxc, dminc, dcountc, dcountf;

        if (ps != null) {
            if (nom_f1 > 0.0) {
                set_nominal_freqs(nom_f1);
            }
            pole = (POLE[]) ps.extHead;
            rmsmax = get_stat_max(pole, ps.length);
            FBIAS = F_BIAS / (.01 * ps.samprate);
            // Setup working values of the cost weights.
            // keep dffact scaled to frame rate
            dffact = (DF_FACT * .01) * ps.samprate; 
            bfact = BAND_FACT / (.01 * ps.samprate);
            ffact = DFN_FACT / (.01 * ps.samprate);
            merge_cost = F_MERGE;
            if (merge_cost > 1000.0) {
                domerge = 0;
            }

            // Allocate space for the formant and bandwidth arrays to be passed
            // back.
            if ((debug & DEB_ENTRY) != 0) {
                System.err.printf("Allocating formant and bandwidth arrays in dpform()\n");
            }
            fr = new double[nform * 2][];
            ba = fr + nform;
            for (i = 0; i < nform * 2; i++) {
                fr[i] = new double[ps.length];
            }
//          cp = new_ext(ps.name,"fb");
//          if (fbs = new_signal(cp,SIG_UNKNOWN, dup_header(ps.header), fr,ps.length, ps.samprate, nform * 2)) {
            if (true) {
                // Allocate space for the raw candidate array.
                if ((debug & DEB_ENTRY) != 0) {
                    System.err.printf("Allocating raw candidate array in dpform()\n");
                }
                pcan = new short[MAXCAN][];
                for (i = 0; i < MAXCAN; i++) {
                    pcan[i] = new short[nform];
                }

                // Allocate space for the dp lattice
                if ((debug & DEB_ENTRY) != 0) {
                    System.err.printf("Allocating DP lattice structure in dpform()\n");
                }
                fl = new FORM[ps.length];
                for (i = 0; i < ps.length; i++) {
                    fl[i] = new FORM();
                }

                //
                // main formant tracking loop
                //
                if ((debug & DEB_ENTRY) != 0) {
                    System.err.printf("Entering main computation loop in dpform()\n");
                }
                for (i = 0; i < ps.length; i++) { // for all analysis frames...

                    // initialize candidate mapping count to 0
                    ncan = 0;

                    // moderate the cost of frequency jumps by the relative
                    // amplitude
                    rmsdffact = pole[i].rms;
                    rmsdffact = rmsdffact / rmsmax;
                    rmsdffact = rmsdffact * dffact;

                    // Get all likely mappings of the poles onto formants for
                    // this frame.
                    // if there ARE pole frequencies available...
                    if (pole[i].npoles != 0) {
                        get_fcand(pole[i].npoles, pole[i].freq, pole[i].band, nform, pcan);

                        // Allocate space for this frame's candidates in the dp
                        // lattice.
                        fl[i].prept = new short[ncan];
                        fl[i].cumerr = new double[ncan];
                        fl[i].cand = new short[ncan][];
                        // allocate cand. slots and install candidates
                        for (j = 0; j < ncan; j++) {
                            fl[i].cand[j] = new short[nform];
                            for (k = 0; k < nform; k++) {
                                fl[i].cand[j][k] = pcan[j][k];
                            }
                        }
                    }
                    fl[i].ncand = (short) ncan;
                    // compute the distance between the current and previous
                    // mappings
                    // for each CURRENT mapping...
                    for (j = 0; j < ncan; j++) {
                        if (i != 0) { // past the first frame?
                            minerr = 0;
                            if (fl[i - 1].ncand != 0) {
                                minerr = 2.0e30;
                            }
                            mincan = -1;
                            // for each PREVIOUS map...
                            for (k = 0; k < fl[i - 1].ncand; k++) {
                                for (pferr = 0.0, l = 0; l < nform; l++) {
                                    ic = fl[i].cand[j][l];
                                    ip = fl[i - 1].cand[k][l];
                                    if ((ic >= 0) && (ip >= 0)) {
                                        ftemp = 2.0 * Math.abs(pole[i].freq[ic] - pole[i - 1].freq[ip]) / (pole[i].freq[ic] + pole[i - 1].freq[ip]);
//                                        ftemp = pole[i].freq[ic] - pole[i - 1].freq[ip];
//                                        if (ftemp >= 0.0) {
//                                            ftemp = ftemp / pole[i - 1].freq[ip];
//                                        } else {
//                                            ftemp = ftemp / pole[i].freq[ic];
//                                        }
                                        // cost prop. to SQUARE of deviation to
                                        // discourage large jumps
                                        pferr += ftemp * ftemp;
                                    } else {
                                        pferr += MISSING;
                                    }
                                }
                                // scale delta-frequency cost and add in prev.
                                // cum. cost
                                conerr = (rmsdffact * pferr) + fl[i - 1].cumerr[k];
                                if (conerr < minerr) {
                                    minerr = conerr;
                                    mincan = k;
                                }
                            } // end for each PREVIOUS mapping...
                        } else { // (i.e. if this is the first frame... )
                            minerr = 0;
                        }

                        // point to best previous mapping
                        fl[i].prept[j] = (short) mincan;
                        // (Note that mincan=-1 if there were no candidates in
                        // prev. fr.)
                        // Compute the local costs for this current mapping.
                        for (k = 0, berr = 0, ferr = 0, fbias = 0; k < nform; k++) {
                            ic = fl[i].cand[j][k];
                            if (ic >= 0) {
                                if (k == 0) { // F1 candidate?
                                    ftemp = pole[i].freq[ic];
                                    merger = (domerge != 0 && (ftemp == pole[i].freq[fl[i].cand[j][1]])) ? merge_cost : 0.0;
                                }
                                berr += pole[i].band[ic];
                                ferr += (Math.abs(pole[i].freq[ic] - fnom[k]) / fnom[k]);
                                fbias += pole[i].freq[ic];
                            } else { // if there was no freq. for this formant
                                fbias += fnom[k];
                                berr += NOBAND;
                                ferr += MISSING;
                            }
                        }

                        // Compute the total cost of this mapping and best
                        // previous.
                        fl[i].cumerr[j] = (FBIAS * fbias) + (bfact * berr) + merger + (ffact * ferr) + minerr;
                    } // end for each CURRENT mapping...

                    if ((debug & DEB_LPC_PARS) != 0) {
                        System.err.printf("\nFrame %4d  # candidates:%3d stat:%f prms:%f", i, ncan, rmsdffact, pole[i].rms);
                        for (j = 0; j < ncan; j++) {
                            System.err.printf("\n	");
                            for (k = 0; k < nform; k++) {
                                if (pcan[j][k] >= 0) {
                                    System.err.printf("%6.0f ", pole[i].freq[fl[i].cand[j][k]]);
                                } else {
                                    System.err.printf("  NA   ");
                                }
                            }
                            System.err.printf("  cum:%7.2f pp:%d", fl[i].cumerr[j], fl[i].prept[j]);
                        }
                    }
                } // end for all analysis frames...
                // Pick the candidate in the final frame with the lowest cost.
                // Starting with that min.-cost cand., work back thru the
                // lattice.
                if ((debug & DEB_ENTRY) != 0) {
                    System.err.printf("Entering backtrack loop in dpform()\n");
                }
                dmaxc = 0;
                dminc = 100;
                dcountc = dcountf = 0;
                for (mincan = -1, i = ps.length - 1; i >= 0; i--) {
                    if ((debug & DEB_LPC_PARS) != 0) {
                        System.err.printf("\nFrame:%4d mincan:%2d ncand:%2d ", i, mincan, fl[i].ncand);
                    }
                    // need to find best starting candidate?
                    if (mincan < 0) {
                        if (fl[i].ncand != 0) {
                            // have candidates at this frame?
                            minerr = fl[i].cumerr[0];
                            mincan = 0;
                            for (j = 1; j < fl[i].ncand; j++) {
                                if (fl[i].cumerr[j] < minerr) {
                                    minerr = fl[i].cumerr[j];
                                    mincan = j;
                                }
                            }
                        }
                    }
                    if (mincan >= 0) {
                        // if there is a "best" candidate at this frame
                        if ((j = fl[i].ncand) > dmaxc) {
                            dmaxc = j;
                        } else if (j < dminc) {
                            dminc = j;
                        }
                        dcountc += j;
                        dcountf++;
                        for (j = 0; j < nform; j++) {
                            k = fl[i].cand[mincan][j];
                            if (k >= 0) {
                                fr[j][i] = pole[i].freq[k];
                                if ((debug & DEB_LPC_PARS) != 0) {
                                    System.err.printf("%6.0f", fr[j][i]);
                                }
                                ba[j][i] = pole[i].band[k];
                            } else { // IF FORMANT IS MISSING...
                                if (i < ps.length - 1) {
                                    // replicate backwards
                                    fr[j][i] = fr[j][i + 1];
                                    ba[j][i] = ba[j][i + 1];
                                } else {
                                    // or insert neutral values
                                    fr[j][i] = fnom[j];
                                    ba[j][i] = NOBAND;
                                }
                                if ((debug & DEB_LPC_PARS) != 0) {
                                    System.err.printf("%6.0f", fr[j][i]);
                                }
                            }
                        }
                        mincan = fl[i].prept[mincan];
                    } else {
                        // if no candidates, fake with "nominal" frequencies.
                        for (j = 0; j < nform; j++) {
                            fr[j][i] = fnom[j];
                            ba[j][i] = NOBAND;
                            if ((debug & DEB_LPC_PARS) != 0) {
                                System.err.printf("%6.0f", fr[j][i]);
                            }
                        }
                    } // note that mincan will remain =-1 if no candidates
                } // end unpacking formant tracks from the dp lattice
                // Deallocate all the DP lattice work space.
//                if (debug & DEB_ENTRY) {
//                    System.err.printf("%s complete; max. cands:%d  min. cands.:%d average cands.:%f\n", fbs.name, dmaxc, dminc, ((double) dcountc) / dcountf);
//                    System.err.printf("Entering memory deallocation in dpform()\n");
//                }
                fl = null;

                fbs = Snack_NewSound(ps.samprate, Sound.SNACK_FLOAT, nform * 2);
                Snack_ResizeSoundStorage(fbs, ps.length);
                for (i = 0; i < ps.length; i++) {
                    for (j = 0; j < nform * 2; j++) {
                        Snack_SetSample(fbs, j, i, (float) fr[j][i]);
                    }
                }
                fbs.length = ps.length;

                return fbs;
            } else {
                System.err.printf("Can't create a new Signal in dpform()\n");
            }
        } else {
            System.err.printf("Bad data pointers passed into dpform()\n");
        }
        return null;
    }

    //
    // lpc_poles.c
    //
    // computation and I/O routines for dealing with LPC poles
    //

    private static final int MAXORDER = 30;

    /** */
    private double integerize(double time, double freq) {
        int i = (int) (.5 + (freq * time));
        return i / freq;
    }

    /** Round the argument to the nearest integer. */
//    private int eround(double flnum) {
//        return flnum >= 0.0 ? (int) (flnum + 0.5) : (int) (flnum - 0.5);
//    }

    /** */
    private Sound lpc_poles(Sound sp, int lpc_ord, int lpc_type, int w_type, double wdur, double frame_int, double preemp) {
        int i, j, size, step, nform, nfrm;
        boolean init;
        POLE[] pole;
        double lpc_stabl = 70.0, energy, normerr;
        double[] lpca = new double[MAXORDER];
        double[] bap = null, frp = null, rhp = null;
        short[] datap, dporg;
        Sound lp;

        if (lpc_type == 1) { // force "standard" stabilized covariance (ala bsa)
            wdur = 0.025;
            preemp = Math.exp(-62.831853 * 90. / sp.samprate); // exp(-1800piT)
        }
        if ((lpc_ord > MAXORDER) || (lpc_ord < 2)/* || (! ((short)sp.data)[0]) */) {
            return null;
        }
        // np = (char)new_ext(sp.name,"pole");
        wdur = integerize(wdur, sp.samprate);
        frame_int = integerize(frame_int, sp.samprate);
        nfrm = 1 + (int) (((((double) sp.length) / sp.samprate) - wdur) / (frame_int));
        if (nfrm >= 1/* lp.buff_size >= 1 */) {
            size = (int) (.5 + (wdur * sp.samprate));
            step = (int) (.5 + (frame_int * sp.samprate));
            pole = new POLE[nfrm/* lp.buff_size */];
            datap = dporg = new short[sp.length];
            for (i = 0; i < Snack_GetLength(sp); i++) {
                datap[i] = (short) Snack_GetSample(sp, 0, i);
            }
            for (j = 0, init = true/*, datap = ((short)sp.data)[0] */;
                 j < nfrm /* lp.buff_size */;
                 j++, datap += step) {
                pole[j] = new POLE();
                pole[j].freq = frp = new double[lpc_ord];
                pole[j].band = bap = new double[lpc_ord];

                switch (lpc_type) {
                case 0:
                    if (!lpc(lpc_ord, lpc_stabl, size, datap, lpca, rhp, null, normerr, energy, preemp, w_type)) {
                        System.err.printf("Problems with lpc in lpc_poles()");
                        break;
                    }
                    break;
                case 1:
                    if (!lpcbsa(lpc_ord, lpc_stabl, size, datap, lpca, rhp, null, normerr, energy, preemp)) {
                        System.err.printf("Problems with lpc in lpc_poles()");
                        break;
                    }
                    break;
                case 2: {
                    int Ord = lpc_ord;
                    double alpha, r0;

                    w_covar(datap, Ord, size, 0, lpca, alpha, r0, preemp, 0);
                    if ((Ord != lpc_ord) || (alpha <= 0.0))
                        System.err.printf("Problems with covar(); alpha:%f  Ord:%d\n", alpha, Ord);
                    energy = Math.sqrt(r0 / (size - Ord));
                }
                    break;
                }
                pole[j].change = 0.0;
                /* don't waste time on low energy frames */
                if ((pole[j].rms = energy) > 1.0) {
                    formant(lpc_ord, (double) sp.samprate, lpca, nform, frp, bap, init);
                    pole[j].npoles = (short) nform;
                    init = false; /* use old poles to start next search */
                } else { /* write out no pole frequencies */
                    pole[j].npoles = 0;
                    init = true; /* restart root search in a neutral zone */
                }

//                 if (debug & 4) {
//                    System.err.printf("\nfr:%4d np:%4d rms:%7.0f  ", j, pole[j].npoles, pole[j].rms);
//                    for (k = 0; k < pole[j].npoles; k++)
//                        System.err.printf(" %7.1f", pole[j].freq[k]);
//                    System.err.printf("\n                   ");
//                    for (k = 0; k < pole[j].npoles; k++)
//                        System.err.printf(" %7.1f", pole[j].band[k]);
//                    System.err.printf("\n");
//                }

            } /* end LPC pole computation for all lp.buff_size frames */
            /* lp.data = (caddr_t)pole; */
            lp = Snack_NewSound((int) (1.0 / frame_int), Sound.LIN16, lpc_ord);
            Snack_ResizeSoundStorage(lp, nfrm);
            for (i = 0; i < nfrm; i++) {
                for (j = 0; j < lpc_ord; j++) {
                    Snack_SetSample(lp, j, i, (float) pole[i].freq[j]);
                }
            }
            lp.length = nfrm;
            lp.extHead = (byte[]) pole;
            return (lp);
        } else {
            System.err.printf("Bad buffer size in lpc_poles()\n");
        }
        return null;
    }

    /** */
    private double frand() {
        return Math.random() / Integer.MAX_VALUE;
    }

    //
    // a quick and dirty interface to bsa's stabilized covariance LPC
    //

    /** max lpc order */
    static final int NPM = 30;

    private boolean lpcbsa(int np, int wind, short[] data, int dataP, double[] lpc, double[] rho, double[] nul1, double[] nul2, double[] energy, double lpc_stabl, double preemp) {
        /* static */int i, mm, owind = 0, wind1;
        /* static */double[] w = new double[1000];
        double[] rc = new double[NPM], phi = new double[NPM * NPM], shi = new double[NPM], sig = new double[1000];
        double xl = .09, fham, amax;
        int psp1, psp3, pspl;

        if (owind != wind) { // need to compute a new window?
            fham = 6.28318506 / wind;
            for (psp1 = 0, i = 0; i < wind; i++, psp1++) {
                w[psp1] = .54 - .46 * Math.cos(i * fham);
            }
            owind = wind;
        }
        wind += np + 1;
        wind1 = wind - 1;

        for (psp3 = 0; psp3 < wind; psp3++) {
            sig[psp3] = data[dataP++] + .016 * frand() - .008;
        }
        for (psp3 = +1; psp3 < wind; psp3++) {
            sig[psp3 - 1] = sig[psp3] - preemp * sig[psp3 - 1];
        }
        for (amax = 0., psp3 = np; psp3 < wind1; psp3++) {
            amax += sig[psp3] * sig[psp3];
        }
        energy[0] = Math.sqrt(amax / owind);
        amax = 1.0 / (energy[0]);

        for (psp3 = 0; psp3 < wind1; psp3++) {
            sig[psp3] *= amax;
        }
        if ((mm = dlpcwtd(sig, wind1, lpc, np, rc, phi, shi, xl, w)) != np) {
            System.err.printf("LPCWTD error mm<np %d %d\n", mm, np);
            return false;
        }
        return true;
    }

    //
    // Copyright (c) 1987, 1988, 1989 AT&T
    // All Rights Reserved
    //
    // THIS IS UNPUBLISHED PROPRIETARY SOURCE CODE OF AT&T
    // The copyright notice above does not evidence any
    // actual or intended publication of such source code.
    //
    // downsample.c
    // a quick and dirty downsampler
    //

    /**
     * create the coefficients for a symmetric FIR lowpass filter using the
     * window technique with a Hanning window.
     */
    private boolean lc_lin_fir(double fc, double coef[], double[] nf) {
        int i;
        double n;
        double twopi, fn, c;

        if (((nf[0] % 2) != 1) || (nf[0] > 127)) {
            if (nf[0] <= 126) {
                nf[0] = nf[0] + 1;
            } else {
                nf[0] = 127;
            }
        }
        n = (nf[0] + 1) / 2;

        // compute part of the ideal impulse response
        twopi = Math.PI * 2.0;
        coef[0] = 2.0 * fc;
        c = Math.PI;
        fn = twopi * fc;
        for (i = 1; i < n; i++) {
            coef[i] = Math.sin(i * fn) / (c * i);
        }

        // Now apply a Hanning window to the (infinite) impulse response.
        fn = twopi / (nf[0] - 1);
        for (i = 0; i < n; i++) {
            coef[i] *= (.5 + (.5 * Math.cos(fn * (i))));
        }

        return true;
    }

    /**
     * ic contains 1/2 the coefficients of a symmetric FIR filter with unity
     * passband gain. This filter is convolved with the signal in buf. The
     * output is placed in buf2. If invert != 0, the filter magnitude response
     * will be inverted.
     */
    private void do_fir(short[] buf, int bufP, short[] bufo, int bufoP, short ic[], int in_samps, int ncoef, int invert) {
        int buft, bufp, bufp2;
        short stem;
        short[] co = new short[256], mem = new short[256];
        int i, j, k, l, m, sum, integral;

        for (i = ncoef - 1, bufp = +ncoef - 1, bufp2 = 0, buft = ((ncoef - 1) * 2), integral = 0; i-- > 0;) {
            if (invert == 0) {
                co[buft--] = co[bufp2++] = ic[bufp--];
            } else {
                integral += (stem = ic[bufp--]);
                co[buft--] = co[bufp2++] = (short) -stem;
            }
        }
        if (invert == 0) {
            co[buft--] = co[bufp2++] = ic[bufp--]; // point of symmetry
        } else {
            integral *= 2;
            integral += ic[bufp];
            co[buft--] = (short) (integral - ic[bufp]);
        }
//      for (i = (ncoef * 2) - 2; i >= 0; i--) {
//          System.err.printf("\n%4d%7d", i, co[i]);
//      }
        for (i = ncoef - 1, buft = 0; i-- > 0;) {
            mem[buft++] = 0;
        }
        for (i = ncoef; i-- > 0;) {
            mem[buft++] = buf[bufP++];
        }
        l = 16384;
        m = 15;
        k = (ncoef << 1) - 1;
        for (i = in_samps - ncoef; i-- > 0;) {
            for (j = k, buft = 0, bufp = 0, bufp2 = 1, sum = 0; j-- > 0; mem[buft++] = mem[bufp2++]) {
                sum += (((co[bufp++] * mem[buft]) + l) >> m);
            }
            mem[--buft] = buf[bufP++]; // new data to memory
            bufo[bufoP++] = (short) sum;
        }
        for (i = ncoef; i-- > 0;) { // pad data end with zeros
            for (j = k, buft = 0, bufp = 0, bufp2 = 1, sum = 0; j-- > 0; mem[buft++] = mem[bufp2++]) {
                sum += (((co[bufp++] * mem[buft]) + l) >> m);
            }
            mem[--buft] = 0;
            bufo[bufoP++] = (short) sum;
        }
    }

    /** */
    private int get_abs_maximum(short[] d, int n) {
        int i;
        short amax, t;

        if ((t = d[0]++) >= 0) {
            amax = t;
        } else {
            amax = (short) -t;
        }

        for (i = n - 1; i-- > 0;) {
            if ((t = d[0]++) > amax) {
                amax = t;
            } else {
                if (-t > amax) {
                    amax = (short) -t;
                }
            }
        }
        return amax;
    }

    /** */
    private boolean dwnsamp(short[] buf, int bufP, short[][] buf2, int in_samps, int[] out_samps, int insert, int decimate, int ncoef, int[] smin, int[] smax, short ic[]) {
        int bufp, bufp2;
        short[] buft;
        int i, j, k, l, m;
        int imax, imin;

        buf2[0] = buft = new short[insert * in_samps];

        k = imax = get_abs_maximum(buf, in_samps);
        if (k == 0) {
            k = 1;
        }
        if (insert > 1) {
            k = (32767 * 32767) / k; // prepare to scale data
        } else {
            k = (16384 * 32767) / k;
        }
        l = 16384;
        m = 15;

        // Insert zero samples to boost the sampling frequency and scale the
        // signal to maintain maximum precision.
        for (i = 0, bufp = 0, bufp2 = bufP; i < in_samps; i++) {
            buft[bufp++] = (short) (((k * (buf[bufp2++])) + l) >> m);
            for (j = 1; j < insert; j++) {
                buft[bufp++] = 0;
            }
        }

        do_fir(buft, in_samps * insert, buft, ncoef, ic, 0);

        /* Finally, decimate and return the downsampled signal. */
        out_samps[0] = j = (in_samps * insert) / decimate;
        k = decimate;
        for (i = 0, bufp = 0, bufp2 = 0, imax = imin = buft[bufp]; i < j; bufp += k, i++) {
            buft[bufp2++] = buft[bufp];
            if (imax < buft[bufp]) {
                imax = buft[bufp];
            } else if (imin > buft[bufp]) {
                imin = buft[bufp];
            }
        }
        smin[0] = imin;
        smax[0] = imax;
        buf2[0] = new short[out_samps[0]];
        return true;
    }

    /** */
    private boolean ratprx(double a, int[] l, int[] k, int qlim) {
        double aa, af, q, em, qq = 0, pp = 0, ps, e;
        int ai, ip, i;

        aa = Math.abs(a);
        ai = (int) aa;
//        af = Math.mod(aa,1.0);
        i = (int) aa;
        af = aa - i;
        q = 0;
        em = 1.0;
        while (++q <= qlim) {
            ps = q * af;
            ip = (int) (ps + 0.5);
            e = Math.abs((ps - ip) / q);
            if (e < em) {
                em = e;
                pp = ip;
                qq = q;
            }
        }
        k[0] = (int) ((ai * qq) + pp);
        k[0] = a > 0 ? k[0] : -k[0];
        l[0] = (int) qq;
        return true;
    }

    /** */
    private Sound Fdownsample(double freq2, Sound s, int start, int end) {
        short[] bufin;
        short[][] bufout;
        /* static */double beta = 0.0;
        /* static */double[] b = new double[256];
        double ratio_t, maxi, ratio, beta_new, freq1;
        /* static */double[] ncoeff = new double[] { 127 };
        /* static */int ncoefft = 0, nbits = 15;
        /* static */short[] ic = new short[256];
        int[] insert = new int[1], decimate = new int[1];
        int out_samps, smin, smax;
        Sound so;

        int i;
        double j;

        freq1 = s.samprate;

        bufout = new short[1][];
        bufin = new short[end - start + 1];
        for (i = start; i <= end; i++) {
            bufin[i - start] = (short) Snack_GetSample(s, 0, i);
        }

        ratio = freq2 / freq1;
        ratprx(ratio, insert, decimate, 10);
        ratio_t = ((double) insert[0]) / ((double) decimate[0]);

        if (ratio_t > .99) {
            return (s);
        }

        freq2 = ratio_t * freq1;
        beta_new = (.5 * freq2) / (insert[0] * freq1);

        if (beta != beta_new) {
            beta = beta_new;
            if (!lc_lin_fir(beta, ncoeff, b)) {
                System.err.printf("\nProblems computing interpolation filter\n");
                return null;
            }
            maxi = (1 << nbits) - 1;
            j = (ncoeff[0] / 2) + 1;
            for (ncoefft = 0, i = 0; i < j; i++) {
                ic[i] = (short) (0.5 + (maxi * b[i]));
                if (ic[i] != 0) {
                    ncoefft = i + 1;
                }
            }
        } // endif new coefficients need to be computed

        if (dwnsamp(bufin, end - start + 1, bufout, out_samps, insert, decimate[0], ncoefft, ic, smin, smax)) {
            /* so.buff_size = so.file_size = out_samps; */
            so = Snack_NewSound(0, Sound.LIN16, s.nchannels);
            Snack_ResizeSoundStorage(so, out_samps);
            for (i = 0; i < out_samps; i++) {
                Snack_SetSample(so, 0, i, (float) bufout[0][i]);
            }
            so.length = out_samps;
            so.samprate = (int) freq2;
            return so;
        } else {
            System.err.printf("Problems in dwnsamp() in downsample()\n");
        }

        return null;
    }

    /** */
    private Sound highpass(Sound s) {

        short[] datain, dataout;
        /* static */short[] lcf;
        /* static */int len = 0;
        double scale, fn;
        int i;
        Sound so;

        // Headerh,dup_header();

        final int LCSIZ = 101;
        // This assumes the sampling frequency is 10kHz and that the FIR is a
        // Hanning function of (LCSIZ/10)ms duration.

        datain = new short[s.length];
        dataout = new short[s.length];
        for (i = 0; i < Snack_GetLength(s); i++) {
            datain[i] = (short) Snack_GetSample(s, 0, i);
        }

        if (len == 0) { // need to create a Hanning FIR?
            lcf = new short[LCSIZ];
            len = 1 + (LCSIZ / 2);
            fn = Math.PI * 2.0 / (LCSIZ - 1);
            scale = 32767.0 / (.5 * LCSIZ);
            for (i = 0; i < len; i++) {
                lcf[i] = (short) (scale * (.5 + (.4 * Math.cos(fn * i))));
            }
        }
        do_fir(datain, s.length, dataout, len, lcf, 1); // in downsample.c
        so = Snack_NewSound(s.samprate, Sound.LIN16, s.nchannels);
        if (so == null) {
            return null;
        }
        Snack_ResizeSoundStorage(so, s.length);
        for (i = 0; i < s.length; i++) {
            Snack_SetSample(so, 0, i, (float) dataout[i]);
        }
        so.length = s.length;
        return so;
    }

    private enum subOptions {
        START,
        END,
        PROGRESS,
        FRAME,
        PREEMP,
        NUMFORM,
        ORDER,
        WINLEN,
        WINTYPE,
        LPCTYPE,
        DSFREQ,
        NOMFREQ
    }

    public int exec(Sound s, Tcl.Interp interp, int objc, Tcl.Obj objv[]) {
        int nform, i, j, lpc_ord, lpc_type, w_type;
        String w_type_str = null;
        double frame_int, wdur, ds_freq, nom_f1 = -10.0, preemp;
        double cor_wdur;
        Sound dssnd = null, hpsnd = null, polesnd = null;
        Sound formantsnd = null, hpsrcsnd, polesrcsnd;
        Tcl.Obj list;
        int arg, startpos = 0, endpos = -1;
        final String[] subOptionStrings = {
            "-start", "-end", "-progress", "-framelength", "-preemphasisfactor", "-numformants", "-lpcorder", "-windowlength", "-windowtype", "-lpctype", "-ds_freq", "-nom_f1_freq", null
        };

        lpc_ord = 12;
        lpc_type = 0; // use bsa's stabilized covariance if != 0
        w_type = 2; // window type: 0=rectangular; 1=Hamming; 2=cos4
        ds_freq = 10000.0;
        wdur = .049; // for LPC analysis
        cor_wdur = .01; // for crosscorrelation F0 estimator
        frame_int = .01;
        preemp = .7;
        nform = 4;

        for (arg = 2; arg < objc; arg += 2) {
            int index;

            if (Tcl.GetIndexFromObj(interp, objv[arg], subOptionStrings, "option", 0, index) != 0) {
                return -1;
            }

            if (arg + 1 == objc) {
                Tcl.AppendResult(interp, "No argument given for ", subOptionStrings[index], " option", (String) null);
                return -1;
            }

            switch (subOptions.values()[index]) {
            case START: {
                if (Tcl.GetIntFromObj(interp, objv[arg + 1], startpos) != 0)
                    return -1;
                break;
            }
            case END: {
                if (Tcl.GetIntFromObj(interp, objv[arg + 1], endpos) != 0)
                    return -1;
                break;
            }
            case PROGRESS: {
                String str = Tcl.GetStringFromObj(objv[arg + 1], null);

                if (str.length() > 0) {
                    Tcl.IncrRefCount(objv[arg + 1]);
                    s.cmdPtr = objv[arg + 1];
                }
                break;
            }
            case FRAME: {
                if (Tcl.GetDoubleFromObj(interp, objv[arg + 1], frame_int) != 0)
                    return -1;
                break;
            }
            case PREEMP: {
                if (Tcl.GetDoubleFromObj(interp, objv[arg + 1], preemp) != 0)
                    return -1;
                break;
            }
            case NUMFORM: {
                if (Tcl.GetIntFromObj(interp, objv[arg + 1], nform) != 0)
                    return -1;
                break;
            }
            case ORDER: {
                if (Tcl.GetIntFromObj(interp, objv[arg + 1], lpc_ord) != 0)
                    return -1;
                break;
            }
            case WINLEN: {
                if (Tcl.GetDoubleFromObj(interp, objv[arg + 1], wdur) != 0)
                    return -1;
                break;
            }
            case WINTYPE: {
                w_type_str = Tcl.GetStringFromObj(objv[arg + 1], null);
                break;
            }
            case LPCTYPE: {
                if (Tcl.GetIntFromObj(interp, objv[arg + 1], lpc_type) != 0)
                    return -1;
                break;
            }
            case DSFREQ: {
                if (Tcl.GetDoubleFromObj(interp, objv[arg + 1], ds_freq) != 0)
                    return -1;
                break;
            }
            case NOMFREQ: {
                if (Tcl.GetDoubleFromObj(interp, objv[arg + 1], nom_f1) != 0)
                    return -1;
                break;
            }
            }
        }
        if (startpos < 0) {
            startpos = 0;
        }
        if (endpos >= (s.length - 1) || endpos == -1) {
            endpos = s.length - 1;
        }
        if (startpos > endpos) {
            return 0;
        }

        // Check for errors in specifying parameters

        if (nform > (lpc_ord - 4) / 2) {
            Tcl.AppendResult(interp, "Number of formants must be <= (lpc order - 4)/2", null);
            return -1;
        }

        if (nform > MAXFORMANTS) {
            Tcl.AppendResult(interp, "A maximum of 7 formants are supported at this time", null);
            return -1;
        }

        if (s.storeType != Sound.SOUND_IN_MEMORY) {
            Tcl.AppendResult(interp, "formant only works with in-memory sounds", (String) null);
            return -1;
        }

        if (w_type_str != null) {
            int len = w_type_str.length();
            if (w_type_str.equalsIgnoreCase("rectangular") || w_type_str.equalsIgnoreCase("0")) {
                w_type = 0;
            } else if (w_type_str.equalsIgnoreCase("hamming") || w_type_str.equalsIgnoreCase("1")) {
                w_type = 1;
            } else if (w_type_str.equalsIgnoreCase("cos^4") || w_type_str.equalsIgnoreCase("2")) {
                w_type = 2;
            } else if (w_type_str.equalsIgnoreCase("hanning") || w_type_str.equalsIgnoreCase("3")) {
                w_type = 3;
            } else {
                Tcl.AppendResult(interp, "unknown window type: ", w_type_str, (String) null);
                return -1;
            }
        }

        Utils.Snack_ProgressCallback(s.cmdPtr, interp, "Computing formants", 0.05);

        if (ds_freq < s.samprate) {
            dssnd = Fdownsample(ds_freq, s, startpos, endpos);
        }

        Utils.Snack_ProgressCallback(s.cmdPtr, interp, "Computing formants", 0.5);

        hpsrcsnd = (dssnd != null ? dssnd : s);

        if (preemp < 1.0) { /* be sure DC and rumble are gone! */
            hpsnd = highpass(hpsrcsnd);
        }

        Utils.Snack_ProgressCallback(s.cmdPtr, interp, "Computing formants", 0.6);

        polesrcsnd = (hpsnd != null ? hpsnd : s);

        if (!(polesnd = lpc_poles(polesrcsnd, wdur, frame_int, lpc_ord, preemp, lpc_type, w_type))) {
            Tcl.AppendResult(interp, "Problems in lpc_poles()", null);
            return -1;
        }

        Utils.Snack_ProgressCallback(s.cmdPtr, interp, "Computing formants", 0.7);

        // LPC poles are now available for the formant estimator.
        if ((formantsnd = dpform(polesnd, nform, nom_f1)) == null) {
            Tcl.AppendResult(interp, "Problems in dpform()", null);
            return -1;
        }

        Utils.Snack_ProgressCallback(s.cmdPtr, interp, "Computing formants", 0.95);

        // SaveSound(formantsnd, interp, "outt.wav", null, 0, null, 0,
        // formantsnd.length, Sound.WAV_STRING);

        if (dssnd != null) {
            Utils.Snack_DeleteSound(dssnd);
        }
        if (hpsnd != null) {
            Utils.Snack_DeleteSound(hpsnd);
        }
        Utils.Snack_DeleteSound(polesnd);

        list = Tcl.NewListObj(0, null);

        for (i = 0; i < formantsnd.length; i++) {
            Tcl.Obj frameList;
            frameList = Tcl.NewListObj(0, null);
            Tcl.ListObjAppendElement(interp, list, frameList);
            for (j = 0; j < nform * 2; j++) {
                Tcl.ListObjAppendElement(interp, frameList, Tcl.NewDoubleObj((double) Snack_GetSample(formantsnd, j, i)));
            }
        }

        Utils.Snack_DeleteSound(formantsnd);

        Utils.Snack_ProgressCallback(s.cmdPtr, interp, "Computing formants", 1.0);

        Tcl.SetObjResult(interp, list);

        return -1;
    }
}

/* */
