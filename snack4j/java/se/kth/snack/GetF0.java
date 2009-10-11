/*
 * This software has been licensed to the Centre of Speech Technology, KTH
 * by Microsoft Corp. with the terms in the accompanying file BSD.txt,
 * which is a BSD style license.
 *
 *    "Copyright (c) 1990-1996 Entropic Research Laboratory, Inc. 
 *                   All rights reserved"
 *
 * Written by:  Derek Lin
 * Checked by:
 * Revised by:  David Talkin
 *
 * Brief description:  Estimates F0 using normalized cross correlation and
 *   dynamic programming.
 *
 */

package se.kth.snack;


class GetF0 {

    class F0_params {
        /** only correlation peaks above this are considered */
        float cand_thresh,
        /** degree to which shorter lags are weighted */
        lag_weight,
        /** weighting given to F0 trajectory smoothness */
        freq_weight,
        /** fixed cost for a voicing-state transition */
        trans_cost,
        /** amplitude-change-modulated VUV trans. cost */
        trans_amp,
        /** spectral-change-modulated VUV trans. cost */
        trans_spec,
        /** fixed bias towards the voiced hypothesis */
        voice_bias,
        /** cost for octave F0 jumps */
        double_cost,
        /** talker-specific mean F0 (Hz) */
        mean_f0,
        /** weight to be given to deviations from mean F0 */
        mean_f0_weight,
        /** min. F0 to search for (Hz) */
        min_f0,
        /** max. F0 to search for (Hz) */
        max_f0,
        /** inter-frame-interval (sec) */
        frame_step,
        /** duration of correlation window (sec) */
        wind_dur;

        /** max. # of F0 cands. to consider at each frame */
        int n_cands,
        /** Specify optional signal pre-conditioning. */
        conditioning;
    }

    static final int BIGSORD = 100;

    /** for storing the crosscorrelation information */
    class Cross {
        /** rms energy in the reference window */
        float rms;
        /** max in the crosscorr. fun. q15 */
        float maxval;
        /**lag # at which max occured */
        short maxloc;
        /** the first non-zero lag computed */
        short firstlag;
        /** the normalized corsscor. fun. q15 */
        float[] correl;
    }

    /** for storing the DP information */
    class Dprec {
        /** # of candidate pitch intervals in the frame */
        short ncands;
        /** locations of the candidates */
        @SuppressWarnings("hiding")
        short[] locs;
        /** peak values of the candidates */
        float[] pvals;
        /** modified peak values of the candidates */
        float[] mpvals;
        /** pointers to best previous cands. */
        short[] prept;
        /** cumulative error for each candidate */
        float[] dpvals;
    }

    /** for lpc stat measure in a window */
    class Windstat {
        float[] rho = new float[BIGSORD + 1];
        float err;
        float rms;
    }

    /** for stationarity measure */
    class Stat {
        @SuppressWarnings("hiding")
        float[] stat;
        float[] rms;
        float[] rms_ratio;
    }

    class Frame {
        Cross cp;
        Dprec dp;
        float rms;
        Frame next;
        Frame prev;
    }

    int debug_level = 0;

    enum subOptions {
        START,
        END,
        F0MAX,
        F0MIN,
        PROGRESS,
        FRAME,
        METHOD,
        WINLEN
    };

    int Get_f0(Sound sound, Tcl.Interp interp, int objc, Tcl.Obj objv[]) {
        float[] fdata;
        boolean done;
        long[] buff_size = new long[1];
        long actsize;
        double sf, start_time;
        F0_params par/* ,read_f0_params() */;
        float[] f0p, vuvp, rms_speech, acpkp;
        int i, vecsize;
        /* int init_dp_f0(), dp_f0(); */
        /* static */int framestep = -1;
        long[] sdstep = new long[] {
            0
        };
        long total_samps;
        int ndone = 0;
        Tcl.Obj list;
        double framestep2 = 0.0, wind_dur;

        int arg, startpos = 0, endpos = -1, fmax, fmin;
        /* static */final String[] subOptionStrings = {
            "-start", "-end", "-maxpitch", "-minpitch", "-progress", "-framelength", "-method", "-windowlength", null
        };

        if (sound.cmdPtr != null) {
            Tcl.DecrRefCount(sound.cmdPtr);
            sound.cmdPtr = null;
        }

        par = new F0_params();
        par.cand_thresh = 0.3f;
        par.lag_weight = 0.3f;
        par.freq_weight = 0.02f;
        par.trans_cost = 0.005f;
        par.trans_amp = 0.5f;
        par.trans_spec = 0.5f;
        par.voice_bias = 0.0f;
        par.double_cost = 0.35f;
        par.min_f0 = 50;
        par.max_f0 = 550;
        par.frame_step = 0.01f;
        par.wind_dur = 0.0075f;
        par.n_cands = 20;
        par.mean_f0 = 200; // unused
        par.mean_f0_weight = 0.0f; // unused
        par.conditioning = 0; // unused

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
            case F0MAX: {
                if (Tcl.GetIntFromObj(interp, objv[arg + 1], fmax) != 0)
                    return -1;
                par.max_f0 = fmax;
                break;
            }
            case F0MIN: {
                if (Tcl.GetIntFromObj(interp, objv[arg + 1], fmin) != 0)
                    return -1;
                par.min_f0 = fmin;
                break;
            }
            case PROGRESS: {
                String str = Tcl.GetStringFromObj(objv[arg + 1], null);

                if (str.length() > 0) {
                    Tcl.IncrRefCount(objv[arg + 1]);
                    sound.cmdPtr = objv[arg + 1];
                }
                break;
            }
            case FRAME: {
                if (Tcl.GetDoubleFromObj(interp, objv[arg + 1], framestep2) != 0)
                    return -1;
                par.frame_step = (float) framestep2;
                break;
            }
            case METHOD: {
                break;
            }
            case WINLEN: {
                if (Tcl.GetDoubleFromObj(interp, objv[arg + 1], wind_dur) != 0)
                    return -1;
                par.wind_dur = (float) wind_dur;
                break;
            }
            }
        }
        if (startpos < 0) {
            startpos = 0;
        }
        if (endpos >= (sound.length - 1) || endpos == -1) {
            endpos = sound.length - 1;
        }
        if (startpos > endpos) {
            return 0;
        }

        sf = sound.samprate;

        if (framestep > 0) {
            /* If a value was specified with -S, use it. */
            par.frame_step = (float) (framestep / sf);
        }
        start_time = 0.0f;
        if (check_f0_params(interp, par, sf) != 0) {
            throw new IllegalArgumentException("invalid/inconsistent parameters -- exiting.", null);
        }

        total_samps = endpos - startpos + 1;
        if (total_samps < ((par.frame_step * 2.0) + par.wind_dur) * sf) {
            throw new IllegalArgumentException("input range too small for analysis by get_f0.", null);
        }
        // Initialize variables in get_f0.c; allocate data structures; determine
        // length and overlap of input frames to read.
        if (init_dp_f0(sf, par, buff_size, sdstep) != 0 || buff_size[0] > Integer.MAX_VALUE || sdstep[0] > Integer.MAX_VALUE) {
            throw new IllegalArgumentException("problem in init_dp_f0().");
        }

        if (debug_level != 0) {
            System.err.printf("init_dp_f0 returned buff_size %ld, sdstep %ld.\n", buff_size, sdstep);
        }

        if (buff_size[0] > total_samps) {
            buff_size[0] = total_samps;
        }

        actsize = Math.min(buff_size[0], sound.length);
        fdata = new float[(int) Math.max(buff_size[0], sdstep[0])];
        list = Tcl.NewListObj(0, null);
        Utils.Snack_ProgressCallback(sound.cmdPtr, interp, "Computing pitch", 0.0);
        ndone = startpos;

        while (true) {
            done = (actsize < buff_size[0]) || (total_samps == buff_size[0]);
            Snack_GetSoundData(sound, ndone, fdata, actsize);
            /* if (sound.debug > 0) Snack_WriteLog("dp_f0...\n"); */
            if (dp_f0(fdata, (int) actsize, (int) sdstep[0], sf, par, f0p, vuvp, rms_speech, acpkp, vecsize, done)) {
                Tcl.AppendResult(interp, "problem in dp_f0().", null);
                return -1;
            }
            /* if (sound.debug > 0) Snack_WriteLogInt("done dp_f0",vecsize); */
            for (i = vecsize - 1; i >= 0; i--) {
                Tcl.Obj frameList;
                frameList = Tcl.NewListObj(0, null);
                Tcl.ListObjAppendElement(interp, list, frameList);
                Tcl.ListObjAppendElement(interp, frameList, Tcl.NewDoubleObj(f0p[i]));
                Tcl.ListObjAppendElement(interp, frameList, Tcl.NewDoubleObj(vuvp[i]));
                Tcl.ListObjAppendElement(interp, frameList, Tcl.NewDoubleObj(rms_speech[i]));
                Tcl.ListObjAppendElement(interp, frameList, Tcl.NewDoubleObj(acpkp[i]));
            }

            if (done) {
                break;
            }

            ndone += sdstep[0];
            actsize = Math.min(buff_size[0], sound.length - ndone);
            total_samps -= sdstep[0];

            if (actsize > total_samps) {
                actsize = total_samps;
            }

            if (true) {
                int res = Utils.Snack_ProgressCallback(sound.cmdPtr, interp, "Computing pitch", (double) ndone / sound.length);
                if (res != 0) {
                    return -1;
                }
            }
        }

        Utils.Snack_ProgressCallback(sound.cmdPtr, interp, "Computing pitch", 1.0);

        free_dp_f0();

        Tcl.SetObjResult(interp, list);

        return 0;
    }

    /**
     * Some consistency checks on parameter values. Return a positive integer if
     * any errors detected, 0 if none.
     */
    private int check_f0_params(Tcl.Interp interp, F0_params par, double sample_freq) {
        int error = 0;
        double dstep;

        if ((par.cand_thresh < 0.01) || (par.cand_thresh > 0.99)) {
            Tcl.AppendResult(interp, "ERROR: cand_thresh parameter must be between [0.01, 0.99].", null);
            error++;
        }
        if ((par.wind_dur > .1) || (par.wind_dur < .0001)) {
            Tcl.AppendResult(interp, "ERROR: wind_dur parameter must be between [0.0001, 0.1].", null);
            error++;
        }
        if ((par.n_cands > 100) || (par.n_cands < 3)) {
            Tcl.AppendResult(interp, "ERROR: n_cands parameter must be between [3,100].", null);
            error++;
        }
        if ((par.max_f0 <= par.min_f0) || (par.max_f0 >= (sample_freq / 2.0)) || (par.min_f0 < (sample_freq / 10000.0))) {
            Tcl.AppendResult(interp, "ERROR: min(max)_f0 parameter inconsistent with sampling frequency.", null);
            error++;
        }
        dstep = ((int) (0.5 + (sample_freq * par.frame_step))) / sample_freq;
        if (dstep != par.frame_step) {
            if (debug_level != 0) {
                Tcl.AppendResult(interp, "Frame step set to exactly match signal sample rate.", null);
            }
            par.frame_step = (float) dstep;
        }
        if ((par.frame_step > 0.1) || (par.frame_step < (1.0 / sample_freq))) {
            Tcl.AppendResult(interp, "ERROR: frame_step parameter must be between [1/sampling rate, 0.1].", null);
            error++;
        }

        return error;
    }

    /** */
    void get_fast_cands(float[] fdata, float[] fdsdata, float[] engref, float[] maxval, float[] peaks, int size, int start, int[] nlags, int[] maxloc, int[] locs, int[] ncand, int ind, int step, int dec, Cross cp, F0_params par) {
        int decind, decstart, decnlags, decsize, i, j;
        int lp;
        int corp;
        float xp, yp, lag_wt;
        int pe;

        lag_wt = par.lag_weight / nlags[0];
        decnlags = 1 + (nlags[0] / dec);
        if ((decstart = start / dec) < 1) {
            decstart = 1;
        }
        decind = (ind * step) / dec;
        decsize = 1 + (size / dec);
        corp = 0;

        crossf(fdsdata, decind, decsize, decstart, decnlags, engref, maxloc, maxval, corp);
        cp.maxloc = (short) maxloc[0]; // location of maximum in correlation
        cp.maxval = maxval[0]; // max. correlation value (found at maxloc)
        cp.rms = (float) Math.sqrt(engref[0] / size); // rms in reference window
        cp.firstlag = (short) decstart;

        // return high peaks in xcorr
        get_cand(cp, peaks, locs, decnlags, ncand, par.cand_thresh);

        // Interpolate to estimate peak locations and values at high sample
        // rate.
        for (i = ncand[0], lp = 0, pe = 0; i-- > 0; pe++, lp++) {
            j = locs[lp] - decstart - 1;
            peak(cp.correl[corp + j], xp, yp);
            // refined lag
            locs[lp] = (locs[lp] * dec) + (int) (0.5 + (xp * dec));
            peaks[pe] = yp * (1.0f - (lag_wt * locs[lp])); // refined amplitude
        }

        if (ncand[0] >= par.n_cands) { // need to prune candidates?
            int loc, locm;
            int lt;
            float smaxval;
            int pem;
            int outer, inner, lim;
            for (outer = 0, lim = par.n_cands - 1; outer < lim; outer++) {
                for (inner = ncand[0] - 1 - outer, pe = (ncand[0]) - 1, pem = pe - 1, loc = (ncand[0]) - 1, locm = loc - 1; inner-- > 0; pe--, pem--, loc--, locm--) {
                    if ((smaxval = peaks[pe]) > peaks[pem]) {
                        peaks[pe] = peaks[pem];
                        peaks[pem] = smaxval;
                        lt = locs[loc];
                        locs[loc] = locs[locm];
                        locs[locm] = lt;
                    }
                }
            }
            // leave room for the unvoiced hypothesis
            ncand[0] = par.n_cands - 1;
        }
        crossfi(fdata, ind * step, size, start, nlags, 7, engref, maxloc, maxval, corp, locs, ncand[0]);

        cp.maxloc = (short) maxloc[0]; // location of maximum in correlation
        cp.maxval = maxval[0]; // max. correlation value (found at maxloc)
        cp.rms = (float) Math.sqrt(engref[0] / size); // rms in reference window
        cp.firstlag = (short) start;
        // return high peaks in xcorr
        get_cand(cp, peaks, locs[0], nlags, ncand, (int) par.cand_thresh);
        if (ncand[0] >= par.n_cands) { // need to prune candidates again?
            int loc, locm;
            int lt;
            float smaxval;
            int /*pe,*/ pem;
            int outer, inner, lim;
            for (outer = 0, lim = par.n_cands - 1; outer < lim; outer++) {
                for (inner = ncand[0] - 1 - outer, pe = (ncand[0]) - 1, pem = pe - 1, loc = (ncand[0]) - 1, locm = loc - 1; inner-- > 0; pe--, pem--, loc--, locm--) {
                    if ((smaxval = peaks[pe]) > peaks[pem]) {
                        peaks[pe] = peaks[pem];
                        peaks[pem] = smaxval;
                        lt = locs[loc];
                        locs[loc] = locs[locm];
                        locs[locm] = lt;
                    }
                }
            }
            // leave room for the unvoiced hypothesis
            ncand[0] = par.n_cands - 1;
        }
    }

    /** */
    float[] downsample(double freq, float[] input, int samsin, int[] samsout, int decimate, int state_idx, int first_time, int last_time) {
        /* static */float[] b = new float[2048];
        /* static */float[] foutput = null;
        float beta = 0.0f;
        /* static */float[] ncoeff = new float[] { 127 };
        /* static */int ncoefft = 0;
        int init;

        if (input != null && (samsin > 0) && (decimate > 0) && samsout[0] != 0) {
            if (decimate == 1) {
                return input;
            }

            if (first_time != 0) {
                int nbuff = (int) ((samsin / decimate) + (2 * ncoeff[0]));

                ncoeff[0] = ((int) (freq * .005)) | 1;
                beta = .5f / decimate;
                foutput = new float[nbuff];
                // assert foutput : "Can't allocate foutput in downsample";
                for (; nbuff > 0;) {
                    foutput[--nbuff] = 0.0f;
                }

                if (!lc_lin_fir(beta, ncoeff, b)) {
                    System.err.printf("\nProblems computing interpolation filter\n");
                    return null;
                }
                ncoefft = (int) ((ncoeff[0] / 2) + 1);
            } // endif new coefficients need to be computed

            if (first_time != 0) {
                init = 1;
            } else if (last_time != 0) {
                init = 2;
            } else {
                init = 0;
            }

            if (downsamp(input, foutput, samsin, samsout, state_idx, decimate, ncoefft, b, init)) {
                return (foutput);
            } else {
                System.err.printf("Problems in downsamp() in downsample()\n");
            }
        } else {
            System.err.printf("Bad parameters passed to downsample()\n");
        }

        return null;
    }

    /** Get likely candidates for F0 peaks. */
    private void get_cand(Cross cross, float[] peak, float cand_thresh, int[] loc, int[] ncand, int nlags) {
        int i, lastl;
        int t;
        float o, p, q, clip;
        int r, s;
        int start, ncan, maxl;

        clip = cand_thresh * cross.maxval;
        maxl = cross.maxloc;
        lastl = nlags - 2;
        start = cross.firstlag;

        r = 0;
        o = cross.correl[r++]; // first point
        q = cross.correl[r++]; // middle point
        p = cross.correl[r++];
        s = 0;
        t = 0;
        ncan = 0;
        for (i = 1; i < lastl; i++, o = q, q = p, p = cross.correl[r++]) {
            if ((q > clip) && // is this a high enough value?
                (q >= p) && (q >= o)) { // NOTE: this finds SHOLDERS and
                                        // PLATEAUS
                // as well as peaks (is this a good idea?)
                peak[s++] = q; // record the peak value
                loc[t++] = i + start; // and its location
                ncan++; // count number of peaks found
            }
        }
//        o = q;
//        q = p;
//        if ((q > clip) && (q >= 0)) {
//            s++ = q;
//            t++ = i + start;
//            ncan++;
//        }
         
        ncand[0] = ncan;
    }

    /**
     * buffer-to-buffer downsample operation
     * This is STRICTLY a decimator! (no upsample)
     */
    private boolean downsamp(float[] in, float[] out, int samples, int[] outsamps, int decimate, int ncoef, int state_idx, float fc[], int init) {
        if (in != null && out != null) {
            do_ffir(in, samples, out, outsamps, state_idx, ncoef, fc, 0, decimate, init);
            return true;
        } else {
            System.err.printf("Bad signal(s) passed to downsamp()\n");
        }
        return false;
    }

    /**
     * fc contains 1/2 the coefficients of a symmetric FIR filter with unity
     * passband gain. This filter is convolved with the signal in buf. The
     * output is placed in buf2. If(invert), the filter magnitude response will
     * be inverted. If(init&1), beginning of signal is in buf; if(init&2), end
     * of signal is in buf. out_samps is set to the number of output points
     * placed in bufo.
     */
    private void do_ffir(float[] buf, float[] bufo, float[] fc, int in_samps, int ncoef, int invert, int skip, int init, int[] out_samps, int idx) {
        int dp1, dp2, dp3;
        float sum, integral;
        /* static */float[] co = null, mem = null;
        /* static */float[] state = new float[1000];
        /* static */int fsize = 0, resid = 0;
        int i, j, k, l;
        int sp;
        float[] buf1;
        int bufP = 0, bufoP = 0; // TODO

        buf1 = buf;
        if (ncoef > fsize) {
            // allocate memory for full coeff. array and filter memory
            fsize = 0;
            i = (ncoef + 1) * 2;
            co = new float[i];
            mem = new float[i];
            fsize = ncoef;
        }

        // fill 2nd half with data
        for (i = ncoef, dp1 = ncoef - 1; i-- > 0;) {
            mem[dp1++] = buf[bufP++];
        }

        if ((init & 1) != 0) { // Is the beginning of the signal in buf?
            // Copy the half-filter and its mirror image into the coefficient
            // array.
            for (i = ncoef - 1, dp3 = ncoef - 1, dp2 = 0, dp1 = ((ncoef - 1) * 2), integral = 0.0f; i-- > 0;)
                if (invert == 0) {
                    co[dp1--] = co[dp2++] = fc[dp3--];
                } else {
                    integral += (sum = fc[dp3--]);
                    co[dp1--] = co[dp2++] = -sum;
                }
            if (invert == 0) {
                co[dp1] = fc[dp3]; // point of symmetry
            } else {
                integral *= 2;
                integral += fc[dp3];
                co[dp1] = integral - fc[dp3];
            }

            for (i = ncoef - 1, dp1 = 0; i-- > 0;) {
                mem[dp1++] = 0;
            }
        } else {
            for (i = ncoef - 1, dp1 = 0, sp = 0; i-- > 0;) {
                mem[dp1++] = state[sp++];
            }
        }

        i = in_samps;
        resid = 0;

        k = (ncoef << 1) - 1; // inner-product loop limit

        if (skip <= 1) { // never used
            // *out_samps = i;
            // for( ; i-- > 0; ) {
            // for(j=k, dp1=mem, dp2=co, dp3=mem+1, sum = 0.0; j-- > 0;
            // *dp1++ = *dp3++ )
            // sum += *dp2++ * *dp1;
            //
            // *--dp1 = *buf++;
            // *bufo++ = (sum < 0.0)? sum -0.5 : sum +0.5;
            // }
            // if(init & 2) {
            // for(i=ncoef; i-- > 0; ) {
            // for(j=k, dp1=mem, dp2=co, dp3=mem+1, sum = 0.0; j-- > 0;
            // *dp1++ = *dp3++ )
            // sum += *dp2++ * *dp1;
            // *--dp1 = 0.0;
            // *bufo++ = (sum < 0)? sum -0.5 : sum +0.5;
            // }
            // *out_samps += ncoef;
            // }
            // return;
        } else { // skip points (e.g. for downsampling)
            // the buffer end is padded with (ncoef-1) data points
            for (l = 0; l < out_samps[0]; l++) {
                for (j = k - skip, dp1 = 0, dp2 = 0, dp3 = skip, sum = 0.0f; j-- > 0; mem[dp1++] = mem[dp3++]) {
                    sum += co[dp2++] * mem[dp1];
                }
                // new data to memory
                for (j = skip; j-- > 0; mem[dp1++] = buf[bufP++]) {
                    sum += co[dp2]++ * mem[dp1];
                }
                bufo[bufoP++] = (sum < 0.0) ? sum - 0.5f : sum + 0.5f;
            }
            if ((init & 2) != 0) {
                resid = in_samps - out_samps[0] * skip;
                for (l = resid / skip; l-- > 0;) {
                    for (j = k - skip, dp1 = 0, dp2 = 0, dp3 = skip, sum = 0.0f; j-- > 0; mem[dp1++] = mem[dp3++]) {
                        sum += co[dp2++] * mem[dp1];
                    }
                    for (j = skip; j-- > 0; mem[dp1++] = 0.0f) {
                        sum += co[dp2++] * mem[dp1];
                    }
                    bufo[bufoP++] = (sum < 0.0) ? sum - 0.5f : sum + 0.5f;
                    (out_samps[0])++;
                }
            } else {
                for (dp3 = idx - ncoef + 1, l = ncoef - 1, sp = 0; l-- > 0;) {
                    state[sp++] = buf1[dp3++];
                }
            }
        }
    }

    /**
     * create the coefficients for a symmetric FIR lowpass filter using the
     * window technique with a Hanning window.
     */
    private boolean lc_lin_fir(
    float fc, float[] coef, float[] nf) {
        int i, n;
        double twopi, fn, c;

        if (((nf[0] % 2) != 1)) {
            nf[0] = nf[0] + 1;
        }
        n = (int) ((nf[0] + 1) / 2);

        // Compute part of the ideal impulse response (the sin(x)/x kernel).
        twopi = Math.PI * 2.0;
        coef[0] = (float) (2.0 * fc);
        c = Math.PI;
        fn = twopi * fc;
        for (i = 1; i < n; i++) {
            coef[i] = (float) (Math.sin(i * fn) / (c * i));
        }

        // Now apply a Hanning window to the (infinite) impulse response.
        // (Probably should use a better window, like Kaiser...)
        fn = twopi / nf[0];
        for (i = 0; i < n; i++) {
            coef[n - i - 1] *= (float) ((.5 - (.5 * Math.cos(fn * (i + 0.5)))));
        }

        return true;
    }

    /**
     * Use parabolic interpolation over the three points defining the peak
     * vicinity to estimate the "true" peak.
     * 
     * @param y vector of length 3 defining peak
     * @param xp x,y values of parabolic peak fitting the input points.
     */
    private void peak(float[] y, float[] xp, float[] yp) {
        float a, c;

        a = (float) ((y[2] - y[1]) + (.5 * (y[0] - y[2])));
        if (Math.abs(a) > .000001) {
            xp[0] = c = (float) ((y[0] - y[2]) / (4.0 * a));
            yp[0] = y[1] - (a * c * c);
        } else {
            xp[0] = 0.0f;
            yp[0] = y[1];
        }
    }

    /*
     * A fundamental frequency estimation algorithm using the normalized cross
     * correlation function and dynamic programming. The algorithm implemented
     * here is similar to that presented by B. Secrest and G. Doddington, "An
     * integrated pitch tracking algorithm for speech systems", Proc. ICASSP-83,
     * pp.1352-1355. It is fully described by D. Talkin,
     * "A robust algorithm for ptich tracking (RAPT)", in W. B. Kleijn & K. K.
     * Paliwal (eds.) Speech Coding and Synthesis, (New York: Elsevier, 1995).
     */

    /*
     * For each frame, up to par.n_cands cross correlation peaks are considered
     * as F0 intervals. Each is scored according to its within- frame properties
     * (relative amplitude, relative location), and according to its
     * connectivity with each of the candidates in the previous frame. An
     * unvoiced hypothesis is also generated at each frame and is considered in
     * the light of voicing state change cost, the quality of the cross
     * correlation peak, and frequency continuity.
     */

    /*
     * At each frame, each candidate has associated with it the following items:
     * its peak value its peak value modified by its within-frame properties its
     * location the candidate # in the previous frame yielding the min. err.
     * (this is the optimum path pointer!) its cumulative cost: (local cost +
     * connectivity cost + cumulative cost of its best-previous-frame-match).
     */

    /*
     * Dynamic programming is then used to pick the best F0 trajectory and
     * voicing state given the local and transition costs for the entire
     * utterance.
     */

    /*
     * To avoid the necessity of computing the full crosscorrelation at the
     * input sample rate, the signal is downsampled; a full ccf is computed at
     * the lower frequency; interpolation is used to estimate the location of
     * the peaks at the higher sample rate; and the fine-grained ccf is computed
     * only in the vicinity of these estimated peak locations.
     */

    /*
     * READ_SIZE: length of input data frame in sec to read DP_CIRCULAR:
     * determines the initial size of DP circular buffer in sec DP_HIST: stored
     * frame history in second before checking for common path DP_CIRCULAR >
     * READ_SIZE, DP_CIRCULAR at least 2 times of DP_HIST DP_LIMIT: in case no
     * convergence is found, DP frames of DP_LIMIT secs are kept before output
     * is forced by simply picking the lowest cost path
     */

    static final double READ_SIZE = 0.2;
    static final double DP_CIRCULAR = 1.5;
    static final double DP_HIST = 0.5;
    static final double DP_LIMIT = 1.0;

    /*
     * stationarity parameters - STAT_WSIZE: window size in sec used in
     * measuring frame energy/stationarity STAT_AINT: analysis interval in sec
     * in measuring frame energy/stationarity
     */
    static final double STAT_WSIZE = 0.030;

    static final double STAT_AINT = 0.020;

    /*
     * headF points to current frame in the circular buffer, tailF points to the
     * frame where tracks start cmpthF points to starting frame of converged
     * path to backtrack
     */

    private Frame headF = null, tailF = null, cmpthF = null;

    /* array for backtracking in convergence check */
    private int[] pcands = null;

    private int cir_buff_growth_count = 0;

    /* # of frames in circular DP buffer */
    private int size_cir_buffer,
    /* # of frames required before convergence test */
    size_frame_hist,
    /* # of frames before forcing output */
    size_frame_out,
    /* # of frames from tailF to headF */
    num_active_frames,
    /* # of frames allocated to output buffers */
    output_buf_size;

    /*
     * DP parameters
     */
    private float tcost, tfact_a, tfact_s, frame_int, vbias, fdouble, wdur, ln2, freqwt, lagwt;

    private int step, size, nlags, start, stop, ncomp;

    int[] locs = null;

    private short maxpeaks;

    /* number of windows seen before resued */
    private int wReuse = 0;

    private Windstat[] windstat = null;

    private float[] f0p = null, vuvp = null, rms_speech = null, acpkp = null, peaks = null;

    private int first_time = 1, pad;

    /** */
    int get_Nframes(long buffsize, int pad, int step) {
        if (buffsize < pad) {
            return 0;
        } else {
            return (int) (buffsize - pad) / step;
        }
    }

    /** */
    int init_dp_f0(double freq, F0_params par, long[] buffsize, long[] sdstep) {
        int nframes;
        int i;
        int stat_wsize, agap, ind, downpatch;

        // reassigning some constants

        tcost = par.trans_cost;
        tfact_a = par.trans_amp;
        tfact_s = par.trans_spec;
        vbias = par.voice_bias;
        fdouble = par.double_cost;
        frame_int = par.frame_step;

        step = (int) Math.round(frame_int * freq);
        size = (int) Math.round(par.wind_dur * freq);
        frame_int = (float) ((step) / freq);
        wdur = (float) ((size) / freq);
        start = (int) Math.round(freq / par.max_f0);
        stop = (int) Math.round(freq / par.min_f0);
        nlags = stop - start + 1;
        // # of samples required by xcorr comp. per fr.
        ncomp = size + stop + 1;
        // maximum number of "peaks" findable in ccf
        maxpeaks = (short) (2 + (nlags / 2));
        ln2 = (float) Math.log(2.0);
        size_frame_hist = (int) (DP_HIST / frame_int);
        size_frame_out = (int) (DP_LIMIT / frame_int);

        // SET UP THE D.P. WEIGHTING FACTORS: The intent is to make the
        // effectiveness of the various fudge factors independent of frame rate
        // or sampling frequency.

        // Lag-dependent weighting factor to emphasize early peaks (higher
        // freqs)
        lagwt = par.lag_weight / stop;

        // Penalty for a frequency skip in F0 per frame
        freqwt = par.freq_weight / frame_int;

        i = (int) (READ_SIZE * freq);
        if (ncomp >= step) {
            nframes = ((i - ncomp) / step) + 1;
        } else {
            nframes = i / step;
        }

        /*
         * buffsize is the number of samples needed to make F0 computation of
         * nframes DP frames possible. The last DP frame is patched with enough
         * points so that F0 computation on it can be carried. F0 computaion on
         * each frame needs enough points to do
         * 
         * 1) xcross or cross correlation measure: enough points to do xcross -
         * ncomp
         * 
         * 2) stationarity measure: enough to make 30 msec windowing possible -
         * ind
         * 
         * 3) downsampling: enough to make filtering possible -- downpatch
         * 
         * So there are nframes whole DP frames, padded with pad points to make
         * the last frame F0 computation ok.
         */

        // last point in data frame needs points of 1/2 downsampler filter
        // length long, 0.005 is the filter length used in downsampler
        downpatch = (((int) (freq * 0.005)) + 1) / 2;

        stat_wsize = (int) (STAT_WSIZE * freq);
        agap = (int) (STAT_AINT * freq);
        ind = (agap - stat_wsize) / 2;
        i = stat_wsize + ind;
        pad = downpatch + ((i > ncomp) ? i : ncomp);
        buffsize[0] = nframes * step + pad;
        sdstep[0] = nframes * step;

        // Allocate space for the DP storage circularly linked data structure

        size_cir_buffer = (int) (DP_CIRCULAR / frame_int);

        // creating circularly linked data structures
        tailF = alloc_frame(nlags, par.n_cands);
        headF = tailF;

        // link them up
        for (i = 1; i < size_cir_buffer; i++) {
            headF.next = alloc_frame(nlags, par.n_cands);
            headF.next.prev = headF;
            headF = headF.next;
        }
        headF.next = tailF;
        tailF.prev = headF;

        headF = tailF;

        // Allocate sscratch array to use during backtrack convergence test.
        if (pcands == null) {
            pcands = new int[par.n_cands];
//          assert pcands : "can't allocate pathcands";
        }

        // Allocate arrays to return F0 and related signals.

        // Note: remember to comparevecsize with size_frame_out, because
        // size_cir_buffer is not constant
        output_buf_size = size_cir_buffer;
        rms_speech = new float[output_buf_size];
//      assert rms_speech : "rms_speech ckalloc failed");
        f0p = new float[output_buf_size];
//      assert f0p : "f0p ckalloc failed");
        vuvp = new float[output_buf_size];
//      assert vuvp : "vuvp ckalloc failed");
        acpkp = new float[output_buf_size];
//      assert acpkp : "acpkp ckalloc failed");

        // Allocate space for peak location and amplitude scratch arrays.
        peaks = new float[maxpeaks];
//      assert peaks,"peaks ckalloc failed");
        locs = new int[maxpeaks];
//      assert locs, "locs ckalloc failed");

        // Initialise the retrieval/saving scheme of window statistic measures
        wReuse = agap / step;
        if (wReuse != 0) {
            windstat = new Windstat[wReuse];
//          assert windstat : "windstat ckalloc failed";
            for (i = 0; i < wReuse; i++) {
                windstat[i].err = 0;
                windstat[i].rms = 0;
            }
        }

        if (debug_level != 0) {
            System.err.printf("done with initialization:\n");
            System.err.printf(" size_cir_buffer:%d  xcorr frame size:%d start lag:%d nlags:%d\n", size_cir_buffer, size, start, nlags);
        }

        num_active_frames = 0;
        first_time = 1;

        return 0;
    }

    /**
     * analysis control parameters
     */
    int dp_f0(float[] fdata, int buff_size, int sdstep,
              double freq, F0_params par,
              float[][] f0p_pt, float[][] vuvp_pt,
              float[][] rms_speech_pt, float[][] acpkp_pt,
              int[] vecsize, int last_time) {
        float maxval, engrefl;
        float[] sta, rms_ratio, dsdata/* ,downsample() */;
        float ttemp, ftemp, ft1, ferr, err, errmin;
        int i, j, k, loc1, loc2;
        int nframes, maxloc, ncand, ncandp, minloc, decimate, samsds;

        Stat stat = null;

        nframes = get_Nframes(buff_size, pad, step); /* # of whole frames */

        if (debug_level != 0) {
            System.err.printf("******* Computing %d dp frames ******** from %d points\n", nframes, buff_size);
        }

        // Now downsample the signal for coarse peak estimates.

        decimate = (int) (freq / 2000.0); // downsample to about 2kHz
        if (decimate <= 1) {
            dsdata = fdata;
        } else {
            samsds = ((nframes - 1) * step + ncomp) / decimate;
            dsdata = downsample(fdata, buff_size, sdstep, freq, samsds, decimate, first_time, last_time);
            if (dsdata == null) {
                throw new IllegalStateException("can't get downsampled data.\n");
            }
        }

        // Get a function of the "stationarity" of the speech signal.

        stat = get_stationarity(fdata, freq, buff_size, nframes, step, first_time);
        if (stat == null) {
            throw new IllegalStateException("can't get stationarity\n");
        }
        sta = stat.stat;
        rms_ratio = stat.rms_ratio;

        // MAIN FUNDAMENTAL FREQUENCY ESTIMATION LOOP
        if (first_time == 0 & nframes > 0) {
            headF = headF.next;
        }

        for (i = 0; i < nframes; i++) {

            // NOTE: This buffer growth provision is probably not necessary. It
            // was put in (with errors) by Derek Lin and apparently never
            // tested. My tests and analysis suggest it is completely
            // superfluous. DT 9/5/96

            // Dynamically allocating more space for the circular buffer
            if (headF == tailF.prev) {
                Frame frm;

                if (cir_buff_growth_count > 5) {
                    throw new IllegalStateException(String.format("too many requests (%d) for dynamically allocating space.\n   There may be a problem in finding converged path.\n", cir_buff_growth_count));
                }
                if (debug_level != 0) {
                    System.err.printf("allocating %d more frames for DP circ. buffer.\n", size_cir_buffer);
                }
                frm = alloc_frame(nlags, par.n_cands);
                headF.next = frm;
                frm.prev = headF;
                for (k = 1; k < size_cir_buffer; k++) {
                    frm.next = alloc_frame(nlags, par.n_cands);
                    frm.next.prev = frm;
                    frm = frm.next;
                }
                frm.next = tailF;
                tailF.prev = frm;
                cir_buff_growth_count++;
            }

            headF.rms = stat.rms[i];
            get_fast_cands(fdata, dsdata, i, step, size, decimate, start, nlags, engref, maxloc, maxval, headF.cp, peaks, locs, ncand, par);

            // Move the peak value and location arrays into the dp structure
            {
                int ftp1, ftp2;
                int sp1;
                int sp2;

                for (ftp1 = 0, ftp2 = 0, sp1 = 0, sp2 = 0, j = ncand; j-- > 0;) {
                    headF.dp.pvals[ftp1++] = peaks[ftp2++];
                    headF.dp.locs[sp1++] = (short) locs[sp2++];
                }
                headF.dp.locs[sp1] = -1; // distinguish the UNVOICED candidate
                headF.dp.pvals[ftp1] = maxval;
                // (high cost if cor. is high)
                headF.dp.mpvals[ncand] = vbias + maxval;
            }

            // Apply a lag-dependent weight to the peaks to encourage the
            // selection of the first major peak. Translate the modified peak
            // values into costs (high peak ==> low cost).
            for (j = 0; j < ncand; j++) {
                ftemp = 1.0f - (locs[j] * lagwt);
                headF.dp.mpvals[j] = 1.0f - (peaks[j] * ftemp);
            }
            ncand++; // include the unvoiced candidate
            headF.dp.ncands = (short) ncand;

            // COMPUTE THE DISTANCE MEASURES AND ACCUMULATE THE COSTS.

            ncandp = headF.prev.dp.ncands;
            // for each of the current candidates...
            for (k = 0; k < ncand; k++) {
                minloc = 0;
                errmin = Float.MAX_VALUE;
                if ((loc2 = headF.dp.locs[k]) > 0) { // current cand. is voiced
                    // for each PREVIOUS candidate...
                    for (j = 0; j < ncandp; j++) {
                        // Get cost due to inter-frame period change.
                        loc1 = headF.prev.dp.locs[j];
                        if (loc1 > 0) { /* prev. was voiced */
                            ftemp = (float) Math.log(((double) loc2) / loc1);
                            ttemp = Math.abs(ftemp);
                            ft1 = (fdouble + Math.abs(ftemp + ln2));
                            if (ttemp > ft1) {
                                ttemp = ft1;
                            }
                            ft1 = (fdouble + Math.abs(ftemp - ln2));
                            if (ttemp > ft1) {
                                ttemp = ft1;
                            }
                            ferr = ttemp * freqwt;
                        } else { // prev. was unvoiced
                            ferr = tcost + (tfact_s * sta[i]) + (tfact_a / rms_ratio[i]);
                        }
                        // Add in cumulative cost associated with previous peak.
                        err = ferr + headF.prev.dp.dpvals[j];
                        if (err < errmin) { // find min. cost
                            errmin = err;
                            minloc = j;
                        }
                    }
                } else { // this is the unvoiced candidate
                    // for each PREVIOUS candidate...
                    for (j = 0; j < ncandp; j++) {

                        // Get voicing transition cost.
                        if (headF.prev.dp.locs[j] > 0) { /* previous was voiced */
                            ferr = tcost + (tfact_s * sta[i]) + (tfact_a * rms_ratio[i]);
                        } else {
                            ferr = 0.0f;
                        }
                        // Add in cumulative cost associated with previous peak.
                        err = ferr + headF.prev.dp.dpvals[j];
                        if (err < errmin) { // find min. cost
                            errmin = err;
                            minloc = j;
                        }
                    }
                }
                // Now have found the best path from this cand. to prev. frame
                if (first_time != 0 && i == 0) { // this is the first frame
                    headF.dp.dpvals[k] = headF.dp.mpvals[k];
                    headF.dp.prept[k] = 0;
                } else {
                    headF.dp.dpvals[k] = errmin + headF.dp.mpvals[k];
                    headF.dp.prept[k] = (short) minloc;
                }
            } // END OF THIS DP FRAME

            if (i < nframes - 1) {
                headF = headF.next;
            }

            if (debug_level >= 2) {
                System.err.printf("%d engref:%10.0f max:%7.5f loc:%4d\n", i, engref, maxval, maxloc);
            }

        } // end for (i ...)

        // DONE WITH FILLING DP STRUCTURES FOR THE SET OF SAMPLED DATA NOW FIND
        // A CONVERGED DP PATH

        vecsize[0] = 0; // # of output frames returned

        num_active_frames += nframes;

        if (num_active_frames >= size_frame_hist || last_time != 0) {
            Frame frm;
            int num_paths, best_cand, frmcnt, checkpath_done = 1;
            float patherrmin;

            if (debug_level != 0) {
                System.err.printf("available frames for backtracking: %d\n", num_active_frames);
            }

            patherrmin = Float.MAX_VALUE;
            best_cand = 0;
            num_paths = headF.dp.ncands;

            // Get the best candidate for the final frame and initialize the
            // paths' backpointers.
            frm = headF;
            for (k = 0; k < num_paths; k++) {
                if (patherrmin > headF.dp.dpvals[k]) {
                    patherrmin = headF.dp.dpvals[k];
                    // index indicating the best candidate at a path
                    best_cand = k;
                }
                pcands[k] = frm.dp.prept[k];
            }

            // Input data was exhausted. force final outputs.
            if (last_time != 0) {
                cmpthF = headF; // Use the current frame as starting point.
            } else {
                // Starting from the most recent frame, trace back each
                // candidate's best path until reaching a common candidate at
                // some past frame.
                frmcnt = 0;
                while (true) {
                    frm = frm.prev;
                    frmcnt++;
                    checkpath_done = 1;
                    for (k = 1; k < num_paths; k++) { // Check for convergence.
                        if (pcands[0] != pcands[k])
                            checkpath_done = 0;
                    }
                    // Prepare for checking at prev. frame.
                    if (checkpath_done == 0) {
                        for (k = 0; k < num_paths; k++) {
                            pcands[k] = frm.dp.prept[pcands[k]];
                        }
                    } else { // All paths have converged.
                        cmpthF = frm;
                        best_cand = pcands[0];
                        if (debug_level != 0) {
                            System.err.printf("paths went back %d frames before converging\n", frmcnt);
                        }
                        break;
                    }
                    if (frm == tailF) { /* Used all available data? */
                        // Delay some more?
                        if (num_active_frames < size_frame_out) {
                            // Yes, don't backtrack at this time.
                            checkpath_done = 0;
                            cmpthF = null;
                        } else {
                            // No more delay! Force best-guess output.
                            checkpath_done = 1;
                            cmpthF = headF;
//                          System.err.printf("WARNING: no converging path found after going back %d frames, will use the lowest cost path\n", num_active_frames);
                        }
                        break;
                    } // end if (frm ...)
                } // end while (1)
            } // end if (last_time) ... else

            // BACKTRACKING FROM cmpthF (best_cand) ALL THE WAY TO tailF
            i = 0;
            frm = cmpthF; // Start where convergence was found (or faked).
            while (frm != tailF.prev && checkpath_done != 0) {
                if (i == output_buf_size) { // Need more room for outputs?
                    output_buf_size *= 2;
                    if (debug_level != 0) {
                        System.err.printf("reallocating space for output frames: %d\n", output_buf_size);
                    }
                    rms_speech = new float[output_buf_size];
//                  assert rms_speech : "rms_speech realloc failed in dp_f0()");
                    f0p = new float[output_buf_size];
//                  assert f0p : "f0p realloc failed in dp_f0()");
                    vuvp = new float[output_buf_size];
//                  assert vuvp : "vuvp realloc failed in dp_f0()");
                    acpkp = new float[output_buf_size];
//                  assert acpkp : "acpkp realloc failed in dp_f0()");
                }
                rms_speech[i] = frm.rms;
                acpkp[i] = frm.dp.pvals[best_cand];
                loc1 = frm.dp.locs[best_cand];
                vuvp[i] = 1.0f;
                best_cand = frm.dp.prept[best_cand];
                ftemp = loc1;
                if (loc1 > 0) { // Was f0 actually estimated for this frame?
                    // loc1 must be a local maximum.
                    if (loc1 > start && loc1 < stop) {
                        float cormax, cprev, cnext, den;

                        j = loc1 - start;
                        cormax = frm.cp.correl[j];
                        cprev = frm.cp.correl[j + 1];
                        cnext = frm.cp.correl[j - 1];
                        den = (float) (2.0 * (cprev + cnext - (2.0 * cormax)));
                        
                        // Only parabolic interpolate if cormax is indeed a
                        // local turning point. Find peak of curve that goes
                        // though the 3 points

                        if (Math.abs(den) > 0.000001) {
                            ftemp += 2.0f - ((((5.0f * cprev) + (3.0f * cnext) - (8.0f * cormax)) / den));
                        }
                    }
                    f0p[i] = (float) (freq / ftemp);
                } else {
                    // No valid estimate; just fake some arbitrary F0.
                    f0p[i] = 0;
                    vuvp[i] = 0.0f;
                }
                frm = frm.prev;

                if (debug_level >= 2) {
                    System.err.printf(" i:%4d%8.1f%8.1f\n", i, f0p[i], vuvp[i]);
                }
                // f0p[i] starts from the most recent one
                // Need to reverse the order in the calling function
                i++;
            } // end while()
            if (checkpath_done != 0) {
                vecsize[0] = i;
                tailF = cmpthF.next;
                num_active_frames -= vecsize[0];
            }
        } // end if()

        if (debug_level != 0) {
            System.err.printf("writing out %d frames.\n", vecsize[0]);
        }

        f0p_pt[0] = f0p;
        vuvp_pt[0] = vuvp;
        acpkp_pt[0] = acpkp;
        rms_speech_pt[0] = rms_speech;
//      acpkp_pt = acpkp;

        if (first_time != 0) {
            first_time = 0;
        }
        return 0;
    }

    /** */
    Frame alloc_frame(int nlags, int ncands) {
        Frame frm;
        int j;

        frm = new Frame();
        frm.dp = new Dprec();
//      assert (frm.dp,"frm.dp ckalloc failed in alloc_frame";
        frm.dp.ncands = 0;
        frm.cp = new Cross();
//      assert frm.cp : "frm.cp ckalloc failed in alloc_frame";
        frm.cp.correl = new float[nlags];
//      assert frm.cp.correl, "frm.cp.correl ckalloc failed";
        // Allocate space for candidates and working arrays.
        frm.dp.locs = new short[ncands];
//      assert frm.dp.locs : "frm.dp.locs ckalloc failed in alloc_frame()";
        frm.dp.pvals = new float[ncands];
//      assert frm.dp.pvals : "frm.dp.pvals ckalloc failed in alloc_frame()";
        frm.dp.mpvals = new float[ncands];
//      assert frm.dp.mpvals : "frm.dp.mpvals ckalloc failed in alloc_frame();
        frm.dp.prept = new short[ncands];
//      assert frm.dp.prept : "frm.dp.prept ckalloc failed in alloc_frame()";
        frm.dp.dpvals = new float[ncands];
//      assert frm.dp.dpvals : "frm.dp.dpvals ckalloc failed in alloc_frame()";

        // Initialize the cumulative DP costs to zero
        for (j = ncands - 1; j >= 0; j--) {
            frm.dp.dpvals[j] = 0.0f;
        }

        return frm;
    }

    /** push window stat to stack, and pop the oldest one */
    private int save_windstat(float[] rho, int order, float err, float rms) {
        int i, j;

        if (wReuse > 1) { // push down the stack
            for (j = 1; j < wReuse; j++) {
                for (i = 0; i <= order; i++) {
                    windstat[j - 1].rho[i] = windstat[j].rho[i];
                }
                windstat[j - 1].err = windstat[j].err;
                windstat[j - 1].rms = windstat[j].rms;
            }
            for (i = 0; i <= order; i++) {
                windstat[wReuse - 1].rho[i] = rho[i]; // save
            }
            windstat[wReuse - 1].err = err;
            windstat[wReuse - 1].rms = rms;
            return 1;
        } else if (wReuse == 1) {
            for (i = 0; i <= order; i++) {
                windstat[0].rho[i] = rho[i]; // save
            }
            windstat[0].err = err;
            windstat[0].rms = rms;
            return 1;
        } else {
            return 0;
        }
    }

    /** */
    private int retrieve_windstat(float[] rho, int order, float[] err, float[] rms) {
        Windstat wstat;
        int i;

        if (wReuse != 0) {
            wstat = windstat[0];
            for (i = 0; i <= order; i++)
                rho[i] = wstat.rho[i];
            err[0] = wstat.err;
            rms[0] = wstat.rms;
            return 1;
        } else {
            return 0;
        }
    }

    /** */
    private float get_similarity(int order, int size, float[] pdata, float[] cdata, float[] rmsa, float[] rms_ratio, float pre, float stab, int w_type, int init) {
        float[] rho3 = new float[BIGSORD + 1], a2 = new float[BIGSORD + 1], rho1 = new float[BIGSORD + 1], a1 = new float[BIGSORD + 1], b = new float[BIGSORD + 1];
        float err3, rms3, rmsd3, b0, t, err1, rms1, rmsd1;

        // (In the lpc() calls below, size-1 is used, since the windowing and
        // preemphasis function assumes an extra point is available in the input
        // data array. This condition is apparently no longer met after Derek's
        // modifications.)

        // get current window stat
        xlpc(order, stab, size - 1, cdata, a2, rho3, (float[]) null, err3, rmsd3, pre, w_type);
        rms3 = wind_energy(cdata, size, w_type);

        if (init == 0) {
            /* get previous window stat */
            if (!retrieve_windstat(rho1, order, err1, rms1)) {
                xlpc(order, stab, size - 1, pdata, a1, rho1, (float[]) null, err1, rmsd1, pre, w_type);
                rms1 = wind_energy(pdata, size, w_type);
            }
            xa_to_aca(a2 + 1, b, b0, order);
            t = xitakura(order, b, b0, rho1 + 1, err1) - .8f;
            if (rms1 > 0.0) {
                rms_ratio[0] = (0.001f + rms3) / rms1;
            } else if (rms3 > 0.0) {
                rms_ratio[0] = 2.0f; // indicate some energy increase
            } else {
                rms_ratio[0] = 1.0f; // no change
            }
        } else {
            rms_ratio[0] = 1.0f;
            t = 10.0f;
        }
        rmsa[0] = rms3;
        save_windstat(rho3, order, err3, rms3);
        return ((float) (0.2 / t));
    }

    /*
     * This is an ad hoc signal stationarity function based on Itakura distance
     * and relative amplitudes.
     */

    private Stat stat = null;

    private float[] mem = null;

    /**
     * This illustrates the window locations when the very first frame is read.
     * It shows an example where each frame step | . | is 10 msec. The frame
     * step size is variable. The window size is always 30 msec. The window
     * centers '*' is always 20 msec apart. The windows cross each other right
     * at the center of the DP frame, or where the '.' is.
     * <pre>
     * ------------------ current window
     * 
     * ------------------ previous window
     * 
     * | . | . | . | . | . | . | . | . | . | ^ ^ ^ ^ ^ ^ ^ ^ fdata ^ ^ ^ q p
     * 
     * --- ind
     * </pre>
     * fdata, q, p, ind, are variables used below.
     */
    private Stat get_stationarity(float[] fdata, double freq, int buff_size, int nframes, int frame_step, int first_time) {
        /* static */int nframes_old = 0, memsize;
        float preemp = 0.4f, stab = 30.0f;
        int p, q, r, datend;
        int ind, i, j, m, size, order, agap, w_type = 3;

        agap = (int) (STAT_AINT * freq);
        size = (int) (STAT_WSIZE * freq);
        ind = (agap - size) / 2;

        if (nframes_old < nframes || stat == null || first_time != 0) {
            // move this to init_dp_f0() later
            nframes_old = nframes;
            stat = new Stat();
//          assert stat : "stat ckalloc failed in get_stationarity";
            stat.stat = new float[nframes];
//          assert stat.stat : "stat.stat ckalloc failed in get_stationarity";
            stat.rms = new float[nframes];
//          assert stat.rms : "stat.rms ckalloc failed in get_stationarity";
            stat.rms_ratio = new float[nframes];
//          assert stat.rms_ratio : "stat.ratio ckalloc failed in get_stationarity";
            memsize = (int) (STAT_WSIZE * freq) + (int) (STAT_AINT * freq);
            mem = new float[memsize];
//          assert mem : "mem ckalloc failed in get_stationarity()";
            for (j = 0; j < memsize; j++) {
                mem[j] = 0;
            }
        }

        if (nframes == 0) {
            return stat;
        }

        q = ind;
        datend = buff_size;

        if ((order = (int) (2.0 + (freq / 1000.0))) > BIGSORD) {
            System.err.printf("Optimim order (%d) exceeds that allowable (%d); reduce Fs\n", order, BIGSORD);
            order = BIGSORD;
        }

        // prepare for the first frame
        for (j = memsize / 2, i = 0; j < memsize; j++, i++) {
            mem[j] = fdata[i];
        }

        // never run over end of frame, should already taken care of when read

        for (j = 0, p = q - agap; j < nframes; j++, p += frame_step, q += frame_step) {
            if ((p >= 0) && (q >= 0) && (q + size <= datend)) {
                stat.stat[j] = get_similarity(order, size, p, q, (stat.rms[j]), (stat.rms_ratio[j]), preemp, stab, w_type, 0);
            } else {
                if (first_time != 0) {
                    if ((p < 0) && (q >= 0) && (q + size <= datend)) {
                        stat.stat[j] = get_similarity(order, size, null, q, (stat.rms[j]), (stat.rms_ratio[j]), preemp, stab, w_type, 1);
                    } else {
                        stat.rms[j] = 0.0f;
                        stat.stat[j] = 0.01f * 0.2f; /* a big transition */
                        stat.rms_ratio[j] = 1.0f; /* no amplitude change */
                    }
                } else {
                    if ((p < 0) && (q + size <= datend)) {
                        stat.stat[j] = get_similarity(order, size, mem, mem + (memsize / 2) + ind, (stat.rms[j]), (stat.rms_ratio[j]), preemp, stab, w_type, 0);
                        /* prepare for the next frame_step if needed */
                        if (p + frame_step < 0) {
                            for (m = 0; m < (memsize - frame_step); m++) {
                                mem[m] = mem[m + frame_step];
                            }
                            r = q + size;
                            for (m = 0; m < frame_step; m++) {
                                mem[memsize - frame_step + m] = fdata[r++];
                            }
                        }
                    }
                }
            }
        }

        // last frame, prepare for next call
        for (j = (memsize / 2) - 1, p = (nframes * frame_step) - 1; j >= 0 && p >= 0; j--) {
            mem[j] = fdata[p--];
        }
        return stat;
    }

    /** Round the argument to the nearest integer. */
//    int eround(double flnum) {
//        return ((flnum >= 0.0) ? (int) (flnum + 0.5) : (int) (flnum - 0.5));
//    }
     
    /** */
    void free_dp_f0() {
        int i;
        Frame frm, next;
        pcands = null;
        rms_speech = null;
        f0p = null;
        vuvp = null;
        acpkp = null;
        peaks = null;
        locs = null;
        if (wReuse != 0) {
            windstat = null;
        }
        frm = headF;
        for (i = 0; i < size_cir_buffer; i++) {
            next = frm.next;
            frm = next;
        }
        headF = null;
        tailF = null;
        stat = null;
        mem = null;
    }

    /** */
    int cGet_f0(Sound sound, Tcl.Interp interp, float[][] outlist, int[] length) {
        float[] fdata;
        boolean done;
        long[] buff_size = new long[1];
        long actsize;
        double sf, start_time;
        F0_params par/* ,read_f0_params() */;
        float[] f0p, vuvp, rms_speech, acpkp;
        int i, vecsize;
        /* static */int framestep = -1;
        long[] sdstep = new long[] {
            0
        };
        long total_samps;
        int ndone = 0;
        Tcl.Obj list;
        float[] tmp = new float[5 + sound.length / 80];
        int count = 0;
        int startpos = 0, endpos = -1;

        if (sound.cmdPtr != null) {
            Tcl.DecrRefCount(sound.cmdPtr);
            sound.cmdPtr = null;
        }

        par = new F0_params();
        par.cand_thresh = 0.3f;
        par.lag_weight = 0.3f;
        par.freq_weight = 0.02f;
        par.trans_cost = 0.005f;
        par.trans_amp = 0.5f;
        par.trans_spec = 0.5f;
        par.voice_bias = 0.0f;
        par.double_cost = 0.35f;
        par.min_f0 = 50;
        par.max_f0 = 550;
        par.frame_step = 0.01f;
        par.wind_dur = 0.0075f;
        par.n_cands = 20;
        par.mean_f0 = 200; // unused
        par.mean_f0_weight = 0.0f; // unused
        par.conditioning = 0; // unused

        if (startpos < 0) {
            startpos = 0;
        }
        if (endpos >= (sound.length - 1) || endpos == -1) {
            endpos = sound.length - 1;
        }
        if (startpos > endpos) {
            return 0;
        }

        sf = sound.samprate;

        if (framestep > 0) { /* If a value was specified with -S, use it. */
            par.frame_step = (float) (framestep / sf);
        }
        start_time = 0.0f;
        if (check_f0_params(interp, par, sf) != 0) {
            Tcl.AppendResult(interp, "invalid/inconsistent parameters -- exiting.", null);
            return -1;
        }

        total_samps = endpos - startpos + 1;
        if (total_samps < ((par.frame_step * 2.0) + par.wind_dur) * sf) {
            Tcl.AppendResult(interp, "input range too small for analysis by get_f0.", null);
            return -1;
        }

        // Initialize variables in get_f0.c; allocate data structures; determine
        // length and overlap of input frames to read.

        if (init_dp_f0(sf, par, buff_size, sdstep) != 0 || buff_size[0] > Integer.MAX_VALUE || sdstep[0] > Integer.MAX_VALUE) {
            Tcl.AppendResult(interp, "problem in init_dp_f0().", null);
            return -1;
        }

        if (debug_level != 0) {
            System.err.printf("init_dp_f0 returned buff_size %ld, sdstep %ld.\n", buff_size, sdstep);
        }

        if (buff_size[0] > total_samps) {
            buff_size[0] = total_samps;
        }

        actsize = Math.min(buff_size[0], sound.length);
        fdata = new float[(int) Math.max(buff_size[0], sdstep[0])];
        list = Tcl.NewListObj(0, null);
        /* Snack_ProgressCallback(sound.cmdPtr, interp, "Computing pitch", 0.0); */
        ndone = startpos;

        while (true) {
            done = (actsize < buff_size[0]) || (total_samps == buff_size[0]);
            Snack_GetSoundData(sound, ndone, fdata, actsize);
            // if (sound.debug > 0) Snack_WriteLog("dp_f0...\n");
            if (dp_f0(fdata, (int) actsize, (int) sdstep[0], sf, par, f0p, vuvp, rms_speech, acpkp, vecsize, done)) {
                Tcl.AppendResult(interp, "problem in dp_f0().", null);
                return -1;
            }
            // if (sound.debug > 0) Snack_WriteLogInt("done dp_f0",vecsize);
            for (i = vecsize - 1; i >= 0; i--) {
                tmp[count] = f0p[i];
                count++;
            }

            if (done) {
                break;
            }

            ndone += sdstep[0];
            actsize = Math.min(buff_size[0], sound.length - ndone);
            total_samps -= sdstep[0];

            if (actsize > total_samps) {
                actsize = total_samps;
            }

//            if (true) {
//                int res = Snack_ProgressCallback(sound.cmdPtr, interp, "Computing pitch", (double) ndone / sound.length);
//                if (res != 0) {
//                    return -1;
//                }
//            }
        }

//      Snack_ProgressCallback(sound.cmdPtr, interp, "Computing pitch", 1.0);

        free_dp_f0();

        outlist[0] = tmp;
        length[0] = count;
//      Tcl.SetObjResult(interp, list);

        return 0;
    }
}

/* */
