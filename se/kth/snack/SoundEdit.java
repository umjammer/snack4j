/* 
 * Copyright (C) 1997-2005 Kare Sjolander <kare@speech.kth.se>
 *
 * This file is part of the Snack Sound Toolkit.
 * The latest version can be found at http://www.speech.kth.se/snack/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

package se.kth.snack;

import vavi.util.win32.WAVE;

class SoundEdit {

    void SnackCopySamples(Sound dest, int to, Sound src, int from, int len) {
        if (dest.storeType == Sound.SOUND_IN_MEMORY) {
            int sn, si, dn, di, tot = 0, blklen;

            to *= src.nchannels[0];
            from *= src.nchannels[0];
            len *= src.nchannels[0];

            if (dest == src && from < to) {
                tot = len;
                if (dest.precision == Sound.SNACK_SINGLE_PREC) {
                    while (tot > 0) {
                        sn = (from + tot) >> Sound.FEXP;
                        si = (from + tot) - (sn << Sound.FEXP);
                        dn = (to + tot) >> Sound.FEXP;
                        di = (to + tot) - (dn << Sound.FEXP);

                        if (di == 0) {
                            blklen = si;
                        } else if (si == 0) {
                            blklen = di;
                        } else {
                            blklen = Math.min(si, di);
                        }

                        blklen = Math.min(blklen, tot);
                        si -= blklen;
                        di -= blklen;

                        if (si < 0) {
                            si = Sound.FBLKSIZE + si;
                            sn--;
                        }
                        if (di < 0) {
                            di = Sound.FBLKSIZE + di;
                            dn--;
                        }
                        if (sn >= src.nblks) {
                            /* Reached end of allocated blocks for src */
                            return;
                        }
                        if (dn >= dest.nblks) {
                            /* Reached end of allocated blocks for dest */
                            return;
                        }
                        System.arraycopy(src.blocks[sn], si, dest.blocks[dn], di, blklen);
                        tot -= blklen;
                    }
                } else {
                    while (tot > 0) {
                        sn = (from + tot) >> Sound.DEXP;
                        si = (from + tot) - (sn << Sound.DEXP);
                        dn = (to + tot) >> Sound.DEXP;
                        di = (to + tot) - (dn << Sound.DEXP);

                        if (di == 0) {
                            blklen = si;
                        } else if (si == 0) {
                            blklen = di;
                        } else {
                            blklen = Math.min(si, di);
                        }

                        blklen = Math.min(blklen, tot);
                        si -= blklen;
                        di -= blklen;

                        if (si < 0) {
                            si = Sound.DBLKSIZE + si;
                            sn--;
                        }
                        if (di < 0) {
                            di = Sound.DBLKSIZE + di;
                            dn--;
                        }
                        if (sn >= src.nblks) {
                            // Reached end of allocated blocks for src
                            return;
                        }
                        if (dn >= dest.nblks) {
                            // Reached end of allocated blocks for dest
                            return;
                        }
                        System.arraycopy(src.blocks[sn], si, dest.blocks[dn], di, blklen);
                        tot -= blklen;
                    }
                }
            } else {
                if (dest.precision == Sound.SNACK_SINGLE_PREC) {
                    while (tot < len) {
                        sn = (from + tot) >> Sound.FEXP;
                        si = (from + tot) - (sn << Sound.FEXP);
                        dn = (to + tot) >> Sound.FEXP;
                        di = (to + tot) - (dn << Sound.FEXP);
                        blklen = Math.min(Sound.FBLKSIZE - si, Sound.FBLKSIZE - di);
                        blklen = Math.min(blklen, len - tot);
                        if (sn >= src.nblks) {
                            /* Reached end of allocated blocks for src */
                            return;
                        }
                        if (dn >= dest.nblks) {
                            /* Reached end of allocated blocks for dest */
                            return;
                        }
                        System.arraycopy(src.blocks[sn], si, dest.blocks[dn], di, blklen);
                        tot += blklen;
                    }
                } else {
                    while (tot < len) {
                        sn = (from + tot) >> Sound.DEXP;
                        si = (from + tot) - (sn << Sound.DEXP);
                        dn = (to + tot) >> Sound.DEXP;
                        di = (to + tot) - (dn << Sound.DEXP);
                        blklen = Math.min(Sound.DBLKSIZE - si, Sound.DBLKSIZE - di);
                        blklen = Math.min(blklen, len - tot);
                        if (sn >= src.nblks) {
                            // Reached end of allocated blocks for src
                            return;
                        }
                        if (dn >= dest.nblks) {
                            // Reached end of allocated blocks for dest
                            return;
                        }
                        System.arraycopy(src.blocks[sn], si, dest.blocks[dn], di, blklen);
                        tot += blklen;
                    }
                }
            }
        } else if (dest.storeType == Sound.SOUND_IN_FILE) {
        }
    }

    void SnackSwapSoundBuffers(Sound s1, Sound s2) {
        int tmpInt;
        Number[][] tmpBlocks;

        tmpBlocks = s1.blocks;
        s1.blocks = s2.blocks;
        s2.blocks = tmpBlocks;

        tmpInt = s1.nblks;
        s1.nblks = s2.nblks;
        s2.nblks = tmpInt;

        tmpInt = s1.exact;
        s1.exact = s2.exact;
        s2.exact = tmpInt;

        tmpInt = s1.maxlength;
        s1.maxlength = s2.maxlength;
        s2.maxlength = tmpInt;
    }

    static void Snack_PutSoundData(Sound s, int pos, byte[] buf, int nSamples) {
        int dn, di, tot = 0, blklen;

        if (s.storeType != Sound.SOUND_IN_MEMORY) {
            return;
        }

        if (s.precision == Sound.SNACK_SINGLE_PREC) {
            while (tot < nSamples) {
                dn = (pos + tot) >> Sound.FEXP;
                di = (pos + tot) - (dn << Sound.FEXP);
                blklen = Math.min(Sound.FBLKSIZE - di, nSamples - tot);
                if (dn >= s.nblks) {
                    // Reached end of allocated blocks for s
                    return;
                }
                System.arraycopy(buf, tot, s.blocks[dn], di, blklen);
                tot += blklen;
            }
        } else {
            while (tot < nSamples) {
                dn = (pos + tot) >> Sound.DEXP;
                di = (pos + tot) - (dn << Sound.DEXP);
                blklen = Math.min(Sound.DBLKSIZE - di, nSamples - tot);
                if (dn >= s.nblks) {
                    // Reached end of allocated blocks for s
                    return;
                }
                System.arraycopy(buf, tot, s.blocks[dn], di, blklen);
                tot += blklen;
            }
        }
    }

    static void Snack_GetSoundData(Sound s, int pos, Object[] buf, int nSamples) {
        int sn, si, tot = 0, blklen;

        if (s.storeType == Sound.SOUND_IN_MEMORY) {
            if (s.precision == Sound.SNACK_SINGLE_PREC) {
                while (tot < nSamples) {
                    sn = (pos + tot) >> Sound.FEXP;
                    si = (pos + tot) - (sn << Sound.FEXP);
                    blklen = Math.min(Sound.FBLKSIZE - si, nSamples - tot);
                    if (sn >= s.nblks) {
                        // Reached end of allocated blocks for s
                        return;
                    }
                    System.arraycopy(s.blocks[sn], si, buf, tot, blklen);
                    tot += blklen;
                }
            } else {
                while (tot < nSamples) {
                    sn = (pos + tot) >> Sound.DEXP;
                    si = (pos + tot) - (sn << Sound.DEXP);
                    blklen = Math.min(Sound.DBLKSIZE - si, nSamples - tot);
                    if (sn >= s.nblks) {
                        // Reached end of allocated blocks for s
                        return;
                    }
                    System.arraycopy(s.blocks[sn], si, buf, tot, blklen);
                    tot += blklen;
                }
            }
        } else if (s.storeType == Sound.SOUND_IN_FILE) {
            int i;

            if (s.linkInfo.linkCh == null) {
                OpenLinkedFile(s, s.linkInfo);
            }

            for (i = 0; i < nSamples; i++) {
                if (s.precision == Sound.SNACK_SINGLE_PREC) {
                    buf[i] = (float) GetSample(s.linkInfo, pos + i);
                } else {
                    buf[i] = (double) GetSample(s.linkInfo, pos + i);
                }
            }
        }
    }

    int lengthCmd(Sound s, Tcl.Interp interp, int objc, Tcl.Obj objv[]) {
        int arg, len, type = 0, newlen = -1, i;
        String string = null;

        if (s.debug > 0) {
            System.err.println("Enter lengthCmd\n");
        }

        if (objc >= 3) {
            for (arg = 2; arg < objc; arg++) {
                string = Tcl.GetStringFromObj(objv[arg], len);

                if (string.equals("-units")) {
                    string = Tcl.GetStringFromObj(objv[arg + 1], len);
                    if (string.equalsIgnoreCase("seconds"))
                        type = 1;
                    if (string.equalsIgnoreCase("samples"))
                        type = 0;
                    arg++;
                } else if (Tcl.GetIntFromObj(interp, objv[2], newlen) != 0) {
                    return -1;
                }
            }
        }

        if (newlen < 0) {
            if (type == 0) {
                Tcl.SetObjResult(interp, Tcl.NewIntObj(s.length));
            } else {
                Tcl.SetObjResult(interp, Tcl.NewDoubleObj((float) s.length / s.samprate));
            }
        } else {
            if (s.storeType != Sound.SOUND_IN_MEMORY) {
                Tcl.AppendResult(interp, "setting sound length only works with", " in-memory sounds", (String) null);
                return -1;
            }
            if (type == 1) {
                newlen *= s.samprate;
            }
            if (newlen > s.length) {
                if (Utils.Snack_ResizeSoundStorage(s, newlen) != 0) {
                    return -1;
                }
                for (i = s.length * s.nchannels[0]; i < newlen * s.nchannels[0]; i++) {
                    switch (s.encoding[0]) {
                    case Sound.LIN16:
                    case Sound.LIN24:
                    case Sound.LIN32:
                    case Sound.ALAW:
                    case Sound.MULAW:
                    case Sound.LIN8:
                    case Sound.SNACK_FLOAT:
                        if (s.precision == Sound.SNACK_SINGLE_PREC) {
                            Sound.FSAMPLE(s, i, 0.0f);
                        } else {
                            Sound.DSAMPLE(s, i, 0.0);
                        }
                        break;
                    case Sound.LIN8OFFSET:
                        if (s.precision == Sound.SNACK_SINGLE_PREC) {
                            Sound.FSAMPLE(s, i, 128.0f);
                        } else {
                            Sound.DSAMPLE(s, i, 128.0);
                        }
                        break;
                    }
                }
            }
            s.length = newlen;
            Utils.Snack_UpdateExtremes(s, 0, s.length, Sound.SNACK_NEW_SOUND);
            Utils.Snack_ExecCallbacks(s, Sound.SNACK_NEW_SOUND);
        }

        if (s.debug > 0) {
            System.err.println("Exit lengthCmd\n");
        }

        return 0;
    }

    int lastIndexCmd(Sound s, Tcl.Interp interp, int objc, Tcl.Obj objv[]) {
        if (s.debug > 0) {
            System.err.println("Enter lastIndexCmd\n");
        }

        if (objc != 2) {
            Tcl.WrongNumArgs(interp, 1, objv, "lastIndex");
            return -1;
        }

        Tcl.SetObjResult(interp, Tcl.NewIntObj(s.length - 1));

        if (s.debug > 0) {
            System.err.println("Exit lastIndexCmd\n");
        }

        return 0;
    }

    enum subOptions {
        START,
        END
    };

    int insertCmd(Sound s, Tcl.Interp interp, int objc, Tcl.Obj objv[]) {
        Sound ins;
        int inspoint, arg, startpos = 0, endpos = -1;
        String string;
        final String[] subOptionStrings = {
            "-start", "-end", null
        };

        if (s.storeType != Sound.SOUND_IN_MEMORY) {
            Tcl.AppendResult(interp, "insert only works with in-memory sounds", (String) null);
            return -1;
        }

        if (objc < 4) {
            Tcl.WrongNumArgs(interp, 1, objv, "insert sound sample");
            return -1;
        }

        string = Tcl.GetStringFromObj(objv[2], null);

        if ((ins = Snack_GetSound(interp, string)) == null) {
            return -1;
        }

        if (Tcl.GetIntFromObj(interp, objv[3], inspoint) != 0) {
            return -1;
        }

        if (inspoint < 0 || inspoint > s.length) {
            Tcl.AppendResult(interp, "Insertion point out of bounds", null);
            return -1;
        }

        if (s.encoding != ins.encoding || s.nchannels != ins.nchannels) {
            Tcl.AppendResult(interp, "Sound format differs: ", string, null);
            return -1;
        }

        for (arg = 4; arg < objc; arg += 2) {
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
            }
        }
        if (startpos < 0)
            startpos = 0;
        if (endpos >= (ins.length - 1) || endpos == -1)
            endpos = ins.length - 1;
        if (startpos > endpos)
            return 0;

        if (Utils.Snack_ResizeSoundStorage(s, s.length + ins.length) != 0) {
            return -1;
        }
        SnackCopySamples(s, inspoint + endpos - startpos + 1, s, inspoint, s.length - inspoint);
        SnackCopySamples(s, inspoint, ins, startpos, endpos - startpos + 1);
        s.length += (endpos - startpos + 1);
        Utils.Snack_UpdateExtremes(s, 0, s.length, Sound.SNACK_NEW_SOUND);
        Utils.Snack_ExecCallbacks(s, Sound.SNACK_NEW_SOUND);

        return 0;
    }

    int cropCmd(Sound s, Tcl.Interp interp, int objc, Tcl.Obj objv[]) {
        int startpos = 0, endpos = 0, totlen;

        if (s.storeType != Sound.SOUND_IN_MEMORY) {
            Tcl.AppendResult(interp, "crop only works with in-memory sounds", (String) null);
            return -1;
        }

        if (objc != 4) {
            Tcl.WrongNumArgs(interp, 1, objv, "crop start end");
            return -1;
        }
        if (Tcl.GetIntFromObj(interp, objv[2], startpos) != 0) {
            return -1;
        }
        if (Tcl.GetIntFromObj(interp, objv[3], endpos) != 0) {
            return -1;
        }

        if ((endpos >= s.length - 1) || endpos < 0) {
            endpos = s.length - 1;
        }
        if (startpos >= endpos) {
            return 0;
        }
        if (startpos < 0) {
            startpos = 0;
        }
        totlen = endpos - startpos + 1;

        SnackCopySamples(s, 0, s, startpos, totlen);
        s.length = totlen;
        Utils.Snack_ExecCallbacks(s, Sound.SNACK_NEW_SOUND);

        return 0;
    }

    enum subOptions2 {
        START,
        END
    };

    int copyCmd(Sound s, Tcl.Interp interp, int objc, Tcl.Obj objv[]) {
        int arg, startpos = 0, endpos = -1;
        Sound master;
        String string;
        final String[] subOptionStrings = {
            "-start", "-end", null
        };

        if (s.storeType != Sound.SOUND_IN_MEMORY) {
            Tcl.AppendResult(interp, "copy only works with in-memory sounds", (String) null);
            return -1;
        }

        if (objc < 3) {
            Tcl.WrongNumArgs(interp, 1, objv, "copy sound");
            return -1;
        }

        string = Tcl.GetStringFromObj(objv[2], null);

        if ((master = Snack_GetSound(interp, string)) == null) {
            return -1;
        }

        if (master.storeType != Sound.SOUND_IN_MEMORY) {
            Tcl.AppendResult(interp, "can only copy from in-memory sounds", (String) null);
            return -1;
        }

        for (arg = 3; arg < objc; arg += 2) {
            int index;

            if (Tcl.GetIndexFromObj(interp, objv[arg], subOptionStrings, "option", 0, index) != 0) {
                return -1;
            }

            if (arg + 1 == objc) {
                Tcl.AppendResult(interp, "No argument given for ", subOptionStrings[index], " option", (String) null);
                return -1;
            }

            switch (subOptions2.values()[index]) {
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
            }
        }
        if (startpos < 0)
            startpos = 0;
        if (endpos >= (master.length - 1) || endpos == -1)
            endpos = master.length - 1;
        if (startpos > endpos)
            return 0;

        s.samprate = master.samprate;
        s.encoding = master.encoding;
        s.sampsize = master.sampsize;
        s.nchannels = master.nchannels;
        s.length = endpos - startpos + 1;
        if (Utils.Snack_ResizeSoundStorage(s, s.length) != 0) {
            return -1;
        }
        SnackCopySamples(s, 0, master, startpos, s.length);
        s.maxsamp = master.maxsamp;
        s.minsamp = master.minsamp;
        s.abmax = master.abmax;
        Utils.Snack_ExecCallbacks(s, Sound.SNACK_NEW_SOUND);

        return 0;
    }

    enum subOptions3 {
        START,
        END,
        MIXSCALE,
        PRESCALE,
        PROGRESS
    };

    int mixCmd(Sound s, Tcl.Interp interp, int objc, Tcl.Obj objv[]) {
        int arg, startpos = 0, endpos = -1, i, c, j;
        double mixscale = 1.0, prescale = 1.0;
        Sound mixsnd;
        String string;
        final String[] subOptionStrings = {
            "-start", "-end", "-mixscaling", "-prescaling", "-progress", null
        };

        if (s.storeType != Sound.SOUND_IN_MEMORY) {
            Tcl.AppendResult(interp, "mix only works with in-memory sounds", (String) null);
            return -1;
        }

        if (objc < 3) {
            Tcl.WrongNumArgs(interp, 1, objv, "mix sound");
            return -1;
        }

        string = Tcl.GetStringFromObj(objv[2], null);

        if ((mixsnd = Snack_GetSound(interp, string)) == null) {
            return -1;
        }

        if (mixsnd.storeType != Sound.SOUND_IN_MEMORY) {
            Tcl.AppendResult(interp, "can only mix from in-memory sounds", (String) null);
            return -1;
        }

        if (s.encoding != mixsnd.encoding || s.nchannels != mixsnd.nchannels) {
            Tcl.AppendResult(interp, "Sound format differs: ", string, null);
            return -1;
        }

        if (s.cmdPtr != null) {
            Tcl.DecrRefCount(s.cmdPtr);
            s.cmdPtr = null;
        }

        for (arg = 3; arg < objc; arg += 2) {
            int index;

            if (Tcl.GetIndexFromObj(interp, objv[arg], subOptionStrings, "option", 0, index) != 0) {
                return -1;
            }

            if (arg + 1 == objc) {
                Tcl.AppendResult(interp, "No argument given for ", subOptionStrings[index], " option", (String) null);
                return -1;
            }

            switch (subOptions3.values()[index]) {
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
            case MIXSCALE: {
                if (Tcl.GetDoubleFromObj(interp, objv[arg + 1], mixscale) != 0)
                    return -1;
                break;
            }
            case PRESCALE: {
                if (Tcl.GetDoubleFromObj(interp, objv[arg + 1], prescale) != 0)
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
        if (endpos - startpos + 1 > mixsnd.length) {
            endpos = startpos + mixsnd.length - 1;
        }

        Utils.Snack_ProgressCallback(s.cmdPtr, interp, "Mixing sound", 0.0);

        for (i = startpos, j = 0; i <= endpos; i++, j++) {
            for (c = 0; c < s.nchannels[0]; c++) {
                float outsmp = (float) (prescale * Sound.FSAMPLE(s, (i * s.nchannels[0] + c)) + mixscale * Sound.FSAMPLE(mixsnd, (j * s.nchannels[0] + c)));
                if (outsmp > 32767.0f) {
                    outsmp = 32767.0f;
                }
                if (outsmp < -32768.0f) {
                    outsmp = -32768.0f;
                }
                Sound.FSAMPLE(s, (i * s.nchannels[0] + c), outsmp);
            }
            if ((i % 100000) == 99999) {
                int res = Utils.Snack_ProgressCallback(s.cmdPtr, interp, "Mixing sound", (double) i / (endpos - startpos));
                if (res != 0) {
                    return -1;
                }
            }
        }

        Utils.Snack_ProgressCallback(s.cmdPtr, interp, "Mixing sound", 1.0);

        Utils.Snack_UpdateExtremes(s, startpos, endpos, Sound.SNACK_NEW_SOUND);
        Utils.Snack_ExecCallbacks(s, Sound.SNACK_NEW_SOUND);

        return 0;
    }

    enum subOptions4 {
        RATE,
        FREQUENCY,
        SKIPHEAD,
        BYTEORDER,
        CHANNELS,
        ENCODING,
        FORMAT,
        START,
        END,
        FILEFORMAT,
        GUESSPROPS
    };

    int appendCmd(Sound s, Tcl.Interp interp, int objc, Tcl.Obj objv[]) {
        Sound t, dummy;
        int arg, startpos = 0, endpos = -1, length = 0;
        String filetype, str;
        final String[] subOptionStrings = {
            "-rate", "-frequency", "-skiphead", "-byteorder", "-channels", "-encoding", "-format", "-start", "-end", "-fileformat", "-guessproperties", null
        };

        if (objc < 3) {
            Tcl.WrongNumArgs(interp, 1, objv, "append variable");
            return -1;
        }
        if (s.storeType != Sound.SOUND_IN_MEMORY) {
            Tcl.AppendResult(interp, "append only works with in-memory sounds", (String) null);
            return -1;
        }

        if ((t = Snack_NewSound(s.samprate, s.encoding, s.nchannels)) == null) {
            Tcl.AppendResult(interp, "Couldn't allocate new sound!", null);
            return -1;
        }
        t.debug = s.debug;
        t.guessEncoding = -1;
        t.guessRate = -1;
        t.swap = 0;

        for (arg = 3; arg < objc; arg += 2) {
            int index;

            if (Tcl.GetIndexFromObj(interp, objv[arg], subOptionStrings, "option", 0, index) != 0) {
                return -1;
            }

            if (arg + 1 == objc) {
                Tcl.AppendResult(interp, "No argument given for ", subOptionStrings[index], " option", (String) null);
                return -1;
            }

            switch (subOptions4.values()[index]) {
            case RATE:
            case FREQUENCY: {
                if (Tcl.GetIntFromObj(interp, objv[arg + 1], t.samprate) != 0) {
                    return -1;
                }
                t.guessRate = 0;
                break;
            }
            case SKIPHEAD: {
                if (Tcl.GetIntFromObj(interp, objv[arg + 1], t.skipBytes) != 0) {
                    return -1;
                }
                break;
            }
            case BYTEORDER: {
                str = Tcl.GetStringFromObj(objv[arg + 1], length);

                if (str.equalsIgnoreCase("littleEndian")) {
                    SwapIfBE(t);
                } else if (str.equalsIgnoreCase("bigEndian")) {
                    SwapIfLE(t);
                } else {
                    Tcl.AppendResult(interp, "-byteorder option should be bigEndian or littleEndian", null);
                    return -1;
                }
                t.guessEncoding = 0;
                break;
            }
            case CHANNELS: {
                if (GetChannels(interp, objv[arg + 1], t.nchannels) != 0)
                    return -1;
                break;
            }
            case ENCODING:
            case FORMAT: {
                if (GetEncoding(interp, objv[arg + 1], t.encoding, t.sampsize) != 0)
                    return -1;
                t.guessEncoding = 0;
                break;
            }
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
            case FILEFORMAT: {
                if (Tcl.GetStringFromObj(objv[arg + 1], null).length() > 0) {
                    if (GetFileFormat(interp, objv[arg + 1], t.fileType) != 0)
                        return -1;
                    t.forceFormat = 1;
                }
                break;
            }
            case GUESSPROPS: {
                int guessProps;
                if (Tcl.GetBooleanFromObj(interp, objv[arg + 1], guessProps) != 0)
                    return -1;
                if (guessProps != 0) {
                    if (t.guessEncoding == -1)
                        t.guessEncoding = 1;
                    if (t.guessRate == -1)
                        t.guessRate = 1;
                }
                break;
            }
            }
        }
        if (t.guessEncoding == -1)
            t.guessEncoding = 0;
        if (t.guessRate == -1)
            t.guessRate = 0;
        if (startpos < 0)
            startpos = 0;
        if (startpos > endpos && endpos != -1)
            return 0;

        str = Tcl.GetStringFromObj(objv[2], length);
        if (length < 10 && (dummy = Snack_GetSound(interp, str)) != null) {
            Tcl.AppendResult(interp, "You must use the concatenate command instead", null);
            return -1;
        }

        filetype = LoadSound(t, interp, objv[2], startpos, endpos);
        if (filetype == null) {
            Utils.Snack_DeleteSound(t);
            return -1;
        }
        if (s.encoding != t.encoding || s.nchannels != t.nchannels) {
            Utils.Snack_DeleteSound(t);
            Tcl.AppendResult(interp, "Sound format differs: ", null);
            return -1;
        }

        if (Utils.Snack_ResizeSoundStorage(s, s.length + t.length) != 0) {
            return -1;
        }
        SnackCopySamples(s, s.length, t, 0, t.length);
        s.length += t.length;
        Utils.Snack_UpdateExtremes(s, s.length - t.length, s.length, Sound.SNACK_MORE_SOUND);
        Utils.Snack_ExecCallbacks(s, Sound.SNACK_MORE_SOUND);
        Utils.Snack_DeleteSound(t);

        return 0;
    }

    enum subOptions5 {
        SMOOTH
    };

    int concatenateCmd(Sound s, Tcl.Interp interp, int objc, Tcl.Obj objv[]) {
        Sound app;
        String string;
        int i, arg, offset = 0, smoothjoin = 0;
        final String[] subOptionStrings = {
            "-smoothjoin", null
        };
        if (s.storeType != Sound.SOUND_IN_MEMORY) {
            Tcl.AppendResult(interp, "concatenate only works with in-memory sounds", (String) null);
            return -1;
        }

        if (objc < 3) {
            Tcl.WrongNumArgs(interp, 1, objv, "concatenate sound");
            return -1;
        }

        string = Tcl.GetStringFromObj(objv[2], null);

        if ((app = Snack_GetSound(interp, string)) == null) {
            return -1;
        }

        if (s.encoding != app.encoding || s.nchannels != app.nchannels) {
            Tcl.AppendResult(interp, "Sound format differs: ", string, null);
            return -1;
        }

        for (arg = 3; arg < objc; arg += 2) {
            int index;

            if (Tcl.GetIndexFromObj(interp, objv[arg], subOptionStrings, "option", 0, index) != 0) {
                return -1;
            }

            if (arg + 1 == objc) {
                Tcl.AppendResult(interp, "No argument given for ", subOptionStrings[index], " option", (String) null);
                return -1;
            }

            switch (subOptions5.values()[index]) {
            case SMOOTH: {
                if (Tcl.GetIntFromObj(interp, objv[arg + 1], smoothjoin) != 0)
                    return -1;
                break;
            }
            }
        }

        if (s.length < smoothjoin) {
            Tcl.AppendResult(interp, "First sound is too short", null);
            return -1;
        }
        if (app.length < 2 * smoothjoin) {
            Tcl.AppendResult(interp, "Second sound is too short", null);
            return -1;
        }
        if (smoothjoin > 0) {
            offset = 80;
            if (s.length < offset)
                offset = s.length - 1;
            for (i = 0; i < offset; i++) {
                float z = (float) ((0.5 + 160.0 / 2 - 1 - i) * 3.141592653589793 / 160);
                float win = (float) Math.exp(-3.0 * z * z);

                Sound.FSAMPLE(s, s.length - offset + i, (float) ((1.0 - win) * Sound.FSAMPLE(s, s.length - offset + i) + win * Sound.FSAMPLE(app, i)));
            }
        } else {
            offset = 0;
        }

        if (Utils.Snack_ResizeSoundStorage(s, s.length + app.length - offset) != 0) {
            return -1;
        }
        SnackCopySamples(s, s.length, app, offset, app.length - offset);
        Utils.Snack_UpdateExtremes(s, s.length, s.length + app.length - offset, Sound.SNACK_MORE_SOUND);
        s.length += (app.length - offset);
        Utils.Snack_ExecCallbacks(s, Sound.SNACK_MORE_SOUND);

        return 0;
    }

    int cutCmd(Sound s, Tcl.Interp interp, int objc, Tcl.Obj objv[]) {
        int start = 0, end = 0;

        if (s.storeType != Sound.SOUND_IN_MEMORY) {
            Tcl.AppendResult(interp, "cut only works with in-memory sounds", (String) null);
            return -1;
        }

        if (objc != 4) {
            Tcl.WrongNumArgs(interp, 1, objv, "cut start end");
            return -1;
        }

        if (Tcl.GetIntFromObj(interp, objv[2], start) != 0) {
            return -1;
        }

        if (Tcl.GetIntFromObj(interp, objv[3], end) != 0) {
            return -1;
        }

        if (start < 0 || start > s.length - 1) {
            Tcl.AppendResult(interp, "Start point out of bounds", null);
            return -1;
        }

        if (end < start || end > s.length - 1) {
            Tcl.AppendResult(interp, "End point out of bounds", null);
            return -1;
        }

        SnackCopySamples(s, start, s, end + 1, s.length - end - 1);
        s.length = s.length - (end - start + 1);
        Utils.Snack_UpdateExtremes(s, 0, s.length, Sound.SNACK_NEW_SOUND);
        Utils.Snack_ExecCallbacks(s, Sound.SNACK_NEW_SOUND);

        return 0;
    }

    int Lowpass(Sound s, Tcl.Interp interp, int rate, int hirate) {
        int c, i;
        float outsmp;
        double insmp = 0.0, last;
        double a = 6.28318530718 * rate / hirate;
        double b = Math.exp(-a / hirate);
        double out;

        for (c = 0; c < s.nchannels[0]; c++) {
            last = 0.0;
            for (i = 0; i < s.length; i++) {
                insmp = Sound.FSAMPLE(s, (i * s.nchannels[0] + c));

                out = insmp * a + last * b;
                last = insmp;
                outsmp = (float) (0.4 * out);

                if (outsmp > 32767.0f) {
                    outsmp = 32767.0f;
                }
                if (outsmp < -32768.0f) {
                    outsmp = -32768.0f;
                }

                Sound.FSAMPLE(s, (i * s.nchannels[0] + c), outsmp);

                if ((i % 100000) == 99999) {
                    int res = Utils.Snack_ProgressCallback(s.cmdPtr, interp, "Converting rate", 0.5 + 0.5 * ((double) (i + c * s.length) / (s.length * s.nchannels[0])));
                    if (res != 0) {
                        return -1;
                    }
                }
            }
        }

        return 0;
    }

    private int Resample(Sound s, Sound t, Tcl.Interp interp) {
        int i, j, c, res, pos;
        float leftsmp = 0.0f, rightsmp;
        double f, frac, dj;

        frac = (double) s.samprate / (double) t.samprate;

        Utils.Snack_ProgressCallback(s.cmdPtr, interp, "Converting rate", 0.0);
        for (c = 0; c < s.nchannels[0]; c++) {

            for (i = 0; i < t.length; i++) {

                dj = frac * i;
                j = (int) dj;
                f = dj - j;

                pos = j * s.nchannels[0] + c;
                leftsmp = Sound.FSAMPLE(s, pos);
                rightsmp = Sound.FSAMPLE(s, pos + s.nchannels[0]);

                Sound.FSAMPLE(t, (i * s.nchannels[0] + c), (float) (leftsmp * (1.0 - f) + rightsmp * f));

                if ((i % 100000) == 99999) {
                    res = Utils.Snack_ProgressCallback(s.cmdPtr, interp, "Converting rate", (0.5 * (i + c * t.length)) / (t.length * s.nchannels[0]));
                    if (res != 0) {
                        Utils.Snack_DeleteSound(t);
                        return -1;
                    }
                }
            }
        }
        res = Lowpass(t, interp, (int) (0.425 * Math.min(t.samprate, s.samprate)), s.samprate);
        if (res != 0) {
            return -1;
        }
        Utils.Snack_ProgressCallback(s.cmdPtr, interp, "Converting rate", 1.0);

        return 0;
    }

    enum subOptions6 {
        RATE,
        FREQUENCY,
        CHANNELS,
        ENCODING,
        FORMAT,
        PROGRESS
    }

    int convertCmd(Sound s, Tcl.Interp interp, int objc, Tcl.Obj objv[]) {
        int arg, i, j;
        int samprate = s.samprate, nchannels = s.nchannels[0];
        int encoding = s.encoding[0], sampsize = s.sampsize[0];
        int snchan = s.nchannels[0];
        Sound t = null;
        final String[] subOptionStrings = {
            "-rate", "-frequency", "-channels", "-encoding", "-format", "-progress", null
        };

        if (s.debug > 0) {
            System.err.println("Enter convertCmd\n");
        }

        if (s.storeType != Sound.SOUND_IN_MEMORY) {
            Tcl.AppendResult(interp, "convert only works with in-memory sounds", (String) null);
            return -1;
        }

        if (objc < 4) {
            Tcl.WrongNumArgs(interp, 1, objv, "convert option value");
            return -1;
        }

        if (s.cmdPtr != null) {
            Tcl.DecrRefCount(s.cmdPtr);
            s.cmdPtr = null;
        }

        for (arg = 2; arg < objc; arg += 2) {
            int index;

            if (Tcl.GetIndexFromObj(interp, objv[arg], subOptionStrings, "option", 0, index) != 0) {
                return -1;
            }

            switch (subOptions6.values()[index]) {
            case RATE:
            case FREQUENCY: {
                if (Tcl.GetIntFromObj(interp, objv[arg + 1], samprate) != 0) {
                    return -1;
                }
                if (samprate < 1) {
                    Tcl.AppendResult(interp, "Rate must be > 1", null);
                    return -1;
                }
                break;
            }
            case CHANNELS: {
//                if (GetChannels(interp, objv[arg + 1], nchannels) != 0) {
//                    return -1;
//                }
                break;
            }
            case ENCODING:
            case FORMAT: {
//                if (GetEncoding(interp, objv[arg + 1], encoding, sampsize) != 0) {
//                    return -1;
//                }
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
            }
        }

        if ((t = Snack_NewSound(samprate, encoding, s.nchannels)) == null) {
            Tcl.AppendResult(interp, "Couldn't allocate temporary sound!", null);
            return -1;
        }
        t.debug = s.debug;
        t.length = (int) (s.length * (float) samprate / s.samprate);
        if (Utils.Snack_ResizeSoundStorage(t, t.length) != 0) {
            Tcl.AppendResult(interp, "Couldn't allocate temporary sound!", null);
            return -1;
        }
        if (samprate != s.samprate) {
            if (s.length > 0) {
                if (Resample(s, t, interp) != 0) {
                    Utils.Snack_DeleteSound(t);
                    return -1;
                }
                SnackSwapSoundBuffers(s, t);
            }
            s.length = t.length;
            s.samprate = t.samprate;
        }
        if (Utils.Snack_ResizeSoundStorage(t, t.length * nchannels) != 0) {
            Tcl.AppendResult(interp, "Couldn't allocate temporary sound!", null);
            return -1;
        }
        t.nchannels[0] = nchannels;

        if (encoding != s.encoding[0]) {
            Utils.Snack_ProgressCallback(s.cmdPtr, interp, "Converting encoding", 0.0);
            for (i = 0; i < s.length * snchan; i++) {
                float value = 0.0f;

                switch (s.encoding[0]) {
                case Sound.LIN16:
                case Sound.SNACK_FLOAT:
                case Sound.ALAW:
                case Sound.MULAW:
                    value = Sound.FSAMPLE(s, i);
                    break;
                case Sound.LIN32:
                    value = Sound.FSAMPLE(s, i) / 65536.0f;
                    break;
                case Sound.LIN24:
                    value = Sound.FSAMPLE(s, i) / 256.0f;
                    break;
                case Sound.LIN8OFFSET:
                    value = (Sound.FSAMPLE(s, i) - 128.0f) * 256.0f;
                    break;
                case Sound.LIN8:
                    value = Sound.FSAMPLE(s, i) * 256.0f;
                    break;
                }

                switch (encoding) {
                case Sound.LIN16:
                case Sound.SNACK_FLOAT:
                    Sound.FSAMPLE(t, i, value);
                    break;
                case Sound.LIN32:
                    Sound.FSAMPLE(t, i, value * 65536.0f);
                    break;
                case Sound.LIN24:
                    Sound.FSAMPLE(t, i, value * 256.0f);
                    break;
                case Sound.ALAW:
                    Sound.FSAMPLE(t, i, (float) Snack_Alaw2Lin(Snack_Lin2Alaw((short) value)));
                    break;
                case Sound.MULAW:
                    Sound.FSAMPLE(t, i, (float) Snack_Mulaw2Lin(Snack_Lin2Mulaw((short) value)));
                    break;
                case Sound.LIN8OFFSET:
                    Sound.FSAMPLE(t, i, (value / 256.0f) + 128.0f);
                    break;
                case Sound.LIN8:
                    Sound.FSAMPLE(t, i, value / 256.0f);
                    break;
                }

                if ((i % 100000) == 99999) {
                    int res = Utils.Snack_ProgressCallback(s.cmdPtr, interp, "Converting encoding", (double) i / (s.length * snchan));
                    if (res != 0) {
                        Utils.Snack_DeleteSound(t);
                        return -1;
                    }
                }
            }
            Utils.Snack_ProgressCallback(s.cmdPtr, interp, "Converting encoding", 1.0);
            SnackSwapSoundBuffers(s, t);
            s.encoding = t.encoding;
            s.sampsize = t.sampsize;
        }

        if (nchannels != snchan) {
            if (nchannels > 1 && snchan > 1) {
                Tcl.AppendResult(interp, "Can only convert n . 1 or 1 . n channels", (String) null);
                Utils.Snack_DeleteSound(t);
                return -1;
            }
            Utils.Snack_ProgressCallback(s.cmdPtr, interp, "Converting channels", 0.0);
            if (nchannels == 1) {
                for (i = 0; i < s.length; i++) {
                    float value = 0.0f;

                    for (j = 0; j < snchan; j++) {
                        value += Sound.FSAMPLE(s, i * snchan + j);
                    }
                    value = value / snchan;

                    Sound.FSAMPLE(t, i, value);

                    if ((i % 100000) == 99999) {
                        int res = Utils.Snack_ProgressCallback(s.cmdPtr, interp, "Converting channels", (double) i / s.length);
                        if (res != 0) {
                            Utils.Snack_DeleteSound(t);
                            return -1;
                        }
                    }
                }
            }
            if (snchan == 1) {
                for (i = s.length - 1; i >= 0; i--) {
                    for (j = 0; j < nchannels; j++) {
                        Sound.FSAMPLE(t, i * nchannels + j, Sound.FSAMPLE(s, i));
                    }
                    if ((i % 100000) == 99999) {
                        int res = Utils.Snack_ProgressCallback(s.cmdPtr, interp, "Converting channels", (double) (s.length - i) / s.length);
                        if (res != 0) {
                            Utils.Snack_DeleteSound(t);
                            return -1;
                        }
                    }
                }
            }
            Utils.Snack_ProgressCallback(s.cmdPtr, interp, "Converting channels", 1.0);
            SnackSwapSoundBuffers(s, t);
            s.nchannels = t.nchannels;
        }
        Utils.Snack_DeleteSound(t);
        Utils.Snack_UpdateExtremes(s, 0, s.length, Sound.SNACK_NEW_SOUND);
        Utils.Snack_ExecCallbacks(s, Sound.SNACK_NEW_SOUND);

        if (s.debug > 0) {
            System.err.println("Exit convertCmd\n");
        }

        return 0;
    }

    enum subOptions7 {
        START,
        END,
        PROGRESS
    };

    int reverseCmd(Sound s, Tcl.Interp interp, int objc, Tcl.Obj objv[]) {
        int arg, startpos = 0, endpos = -1, i, j, c;
        final String[] subOptionStrings = {
            "-start", "-end", "-progress", null
        };

        if (s.storeType != Sound.SOUND_IN_MEMORY) {
            Tcl.AppendResult(interp, "reverse only works with in-memory sounds", (String) null);
            return -1;
        }

        if (s.cmdPtr != null) {
            Tcl.DecrRefCount(s.cmdPtr);
            s.cmdPtr = null;
        }

        if (objc < 2) {
            Tcl.WrongNumArgs(interp, 1, objv, "reverse");
            return -1;
        }

        for (arg = 2; arg < objc; arg += 2) {
            int index;

            if (Tcl.GetIndexFromObj(interp, objv[arg], subOptionStrings, "option", 0, index) != 0) {
                return -1;
            }

            if (arg + 1 == objc) {
                Tcl.AppendResult(interp, "No argument given for ", subOptionStrings[index], " option", (Sound) null);
                return -1;
            }

            switch (subOptions7.values()[index]) {
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
            }
        }
        if (startpos < 0)
            startpos = 0;
        if (endpos >= (s.length - 1) || endpos == -1)
            endpos = s.length - 1;
        if (startpos > endpos)
            return 0;

        if (s.writeStatus == Sound.WRITE) {
            Snack_StopSound(s, interp);
        }

        Utils.Snack_ProgressCallback(s.cmdPtr, interp, "Reversing sound", 0.0);

        for (i = startpos, j = endpos; i <= startpos + (endpos - startpos) / 2; i++, j--) {
            for (c = 0; c < s.nchannels[0]; c++) {
                float swap = Sound.FSAMPLE(s, i * s.nchannels[0] + c);
                Sound.FSAMPLE(s, i * s.nchannels[0] + c, Sound.FSAMPLE(s, j * s.nchannels[0] + c));
                Sound.FSAMPLE(s, j * s.nchannels[0] + c, swap);
                if ((i % 100000) == 99999) {
                    int res = Utils.Snack_ProgressCallback(s.cmdPtr, interp, "Reversing sound", (double) i / (startpos + (endpos - startpos) / 2));
                    if (res != 0) {
                        return -1;
                    }
                }
            }
        }

        Utils.Snack_ProgressCallback(s.cmdPtr, interp, "Reversing sound", 1.0);

        Utils.Snack_ExecCallbacks(s, Sound.SNACK_NEW_SOUND);

        return 0;
    }

    int sampleCmd(Sound s, Tcl.Interp interp, int objc, Tcl.Obj objv[]) {
        int i, n, val;
        double fval;
        String buf;

        if (objc < 3) {
            Tcl.WrongNumArgs(interp, 1, objv, "sample index ?val? ...");
            return -1;
        }

        if (Tcl.GetIntFromObj(interp, objv[2], i) != 0)
            return -1;
        if (i < 0 || i >= s.length) {
            Tcl.AppendResult(interp, "Index out of bounds", null);
            return -1;
        }
        if (objc > 3 && objc > s.nchannels[0] + 3) {
            Tcl.AppendResult(interp, "Too many samples given", null);
            return -1;
        }

        i *= s.nchannels[0];

        if (objc < 4) {
            if (s.storeType != Sound.SOUND_IN_MEMORY && s.linkInfo.linkCh == null) {
                OpenLinkedFile(s, s.linkInfo);
            }

            for (n = 0; n < s.nchannels[0]; n++, i++) {
                switch (s.encoding[0]) {
                case Sound.LIN16:
                case Sound.LIN32:
                case Sound.LIN24:
                case Sound.ALAW:
                case Sound.MULAW:
                case Sound.LIN8OFFSET:
                case Sound.LIN8:
                    if (s.storeType == Sound.SOUND_IN_MEMORY) {
                        if (s.precision == Sound.SNACK_SINGLE_PREC) {
                            buf = String.format("%d", (int) Sound.FSAMPLE(s, i));
                        } else {
                            buf = String.format("%d", (int) Sound.DSAMPLE(s, i));
                        }
                    } else {
                        buf = String.format("%d", (int) GetSample(s.linkInfo, i));
                    }
                    break;
                case Sound.SNACK_FLOAT:
                case Sound.SNACK_DOUBLE:
                    if (s.storeType == Sound.SOUND_IN_MEMORY) {
                        if (s.precision == Sound.SNACK_SINGLE_PREC) {
                            buf = String.format("%f", Sound.FSAMPLE(s, i));
                        } else {
                            buf = String.format("%.12f", Sound.DSAMPLE(s, i));
                        }
                    } else {
                        buf = String.format("%f", GetSample(s.linkInfo, i));
                    }
                    break;
                }
                if (n < s.nchannels[0] - 1) {
                    Tcl.AppendResult(interp, buf, " ", null);
                } else {
                    Tcl.AppendResult(interp, buf, null);
                }
            }
        } else {
            if (s.storeType != Sound.SOUND_IN_MEMORY) {
                Tcl.AppendResult(interp, "setting sample values only works with in-memory sounds", (String) null);
                return -1;
            }
            for (n = 3; n < 3 + s.nchannels[0]; n++, i++) {
                String str;
                int len;

                if (n >= objc)
                    break;
                str = Tcl.GetStringFromObj(objv[n], len);
                if (str.equals("?"))
                    continue;
                if (s.encoding[0] == Sound.SNACK_FLOAT || s.encoding[0] == Sound.SNACK_DOUBLE) {
                    if (Tcl.GetDoubleFromObj(interp, objv[n], fval) != 0)
                        return -1;
                    /*
                     * if (fval < -32768.0 || fval > 32767.0) {
                     * Tcl.AppendResult(interp,
                     * "Sample value not in range -32768, 32767", null); return
                     * -1; }
                     */
                } else {
                    if (Tcl.GetIntFromObj(interp, objv[n], val) != 0)
                        return -1;
                }
                switch (s.encoding[0]) {
                case Sound.LIN16:
                case Sound.ALAW:
                case Sound.MULAW:
                    if (val < -32768 || val > 32767) {
                        Tcl.AppendResult(interp, "Sample value not in range -32768, 32767", null);
                        return -1;
                    }
                case Sound.LIN32:
                case Sound.LIN24:
                    if (val < -8388608 || val > 8388607) {
                        Tcl.AppendResult(interp, "Sample value not in range -8388608, 8388607", null);
                        return -1;
                    }
                    if (s.precision == Sound.SNACK_SINGLE_PREC) {
                        Sound.FSAMPLE(s, i, val);
                    } else {
                        Sound.DSAMPLE(s, i, val);
                    }
                    break;
                case Sound.SNACK_FLOAT:
                case Sound.SNACK_DOUBLE:
                    if (s.precision == Sound.SNACK_SINGLE_PREC) {
                        Sound.FSAMPLE(s, i, (float) fval);
                    } else {
                        Sound.DSAMPLE(s, i, fval);
                    }
                    break;
                case Sound.LIN8OFFSET:
                    if (val < 0 || val > 255) {
                        Tcl.AppendResult(interp, "Sample value not in range 0, 255", null);
                        return -1;
                    }
                    if (s.precision == Sound.SNACK_SINGLE_PREC) {
                        Sound.FSAMPLE(s, i, val);
                    } else {
                        Sound.DSAMPLE(s, i, val);
                    }
                    break;
                case Sound.LIN8:
                    if (val < -128 || val > 127) {
                        Tcl.AppendResult(interp, "Sample value not in range -128, 127", null);
                        return -1;
                    }
                    if (s.precision == Sound.SNACK_SINGLE_PREC) {
                        Sound.FSAMPLE(s, i, val);
                    } else {
                        Sound.DSAMPLE(s, i, val);
                    }
                    break;
                }
            }
        }

        return 0;
    }

    int swapCmd(Sound s, Tcl.Interp interp, int objc, Tcl.Obj objv[]) {
        Sound t;
        String string;
        int tmpInt;
        float tmpFloat;

        if (s.storeType != Sound.SOUND_IN_MEMORY) {
            Tcl.AppendResult(interp, "reverse only works with in-memory sounds", (String) null);
            return -1;
        }

        if (objc < 3) {
            Tcl.WrongNumArgs(interp, 1, objv, "swap sound");
            return -1;
        }

        string = Tcl.GetStringFromObj(objv[2], null);

        if ((t = Snack_GetSound(interp, string)) == null) {
            return -1;
        }

        if (s.encoding != t.encoding || s.nchannels != t.nchannels || s.samprate != t.samprate) {
            Tcl.AppendResult(interp, "Sound format differs: ", string, null);
            return -1;
        }

        SnackSwapSoundBuffers(s, t);

        tmpFloat = s.maxsamp;
        s.maxsamp = t.maxsamp;
        t.maxsamp = tmpFloat;

        tmpFloat = s.minsamp;
        s.minsamp = t.minsamp;
        t.minsamp = tmpFloat;

        tmpFloat = s.abmax;
        s.abmax = t.abmax;
        t.abmax = tmpFloat;

        tmpInt = s.length;
        s.length = t.length;
        t.length = tmpInt;

        Utils.Snack_ExecCallbacks(s, Sound.SNACK_NEW_SOUND);
        Utils.Snack_ExecCallbacks(t, Sound.SNACK_NEW_SOUND);

        return 0;
    }

    int byteswapCmd(Sound s, Tcl.Interp interp, int objc, Tcl.Obj objv[]) {

        // float *block = (float *) ckalloc(1000000);
        // float *block2 = (float *) ckalloc(1000000);
        // int i, j, jmax = 100;
        // double time;
        //
        // time = SnackCurrentTime();
        // for (j = 0; j < jmax ; j++) {
        // memcpy(block, block2, 1000000);
        // }
        // Snack_WriteLogInt("memcpy",
        // (int)(1000000*(SnackCurrentTime()-time)));
        //
        // time = SnackCurrentTime();
        // for (j = 0; j < jmax ; j++) {
        // for (i = 0; i < 250000; i++) {
        // block[i] = block2[i];
        // }
        // }
        // Snack_WriteLogInt(" =[] ", (int)(1000000*(SnackCurrentTime()-time)));
        //
        // time = SnackCurrentTime();
        // for (j = 0; j < jmax ; j++) {
        // float *p = block, *q = block2;
        // for (i = 0; i < 250000; i++) {
        // *p++ = *q++;
        // }
        // }
        // Snack_WriteLogInt(" =++ ", (int)(1000000*(SnackCurrentTime()-time)));

        return 0;
    }

    /* byte reverse command for Snack qzhou@lucent.com 2-3-2000 */

    /* byte reverse (bit swap) macro */
    static byte FLIP_BITS(byte byte_) {
        return (byte) ((((byte_ >> 7) & 0x01) | ((byte_ >> 5) & 0x02) | ((byte_ >> 3) & 0x04) | ((byte_ >> 1) & 0x08) | ((byte_ << 1) & 0x10) | ((byte_ << 3) & 0x20) | ((byte_ << 5) & 0x40) | ((byte_ << 7) & 0x80)) & 0xFF);
    }

    /* Function to add flipBits command, not updated for 2.0 ! */
    int flipBitsCmd(Sound s, Tcl.Interp interp, int objc, Tcl.Obj objv[]) {
        /*
         * int i; unsigned charsampCh;
         */

        if (s.storeType != Sound.SOUND_IN_MEMORY) {
            Tcl.AppendResult(interp, "flipBits only works with in-memory sounds", (String) null);
            return -1;
        }

        if (objc != 2) {
            Tcl.WrongNumArgs(interp, 1, objv, "flipBits");
            return -1;
        }

        if (s.encoding[0] == Sound.MULAW) {
            Tcl.AppendResult(interp, "flipBits only works with Mulaw sounds", (String) null);
            return -1;
        }

        /* stop writing */
        if (s.writeStatus == Sound.WRITE) {
            Snack_StopSound(s, interp);
        }

        /* flip bits for every byte */
        /*
         * for ( i = 0; i < s.length; i++ ) { sampCh = &(UCSAMPLE(s, i));sampCh
         * = FLIP_BITS(sampCh ); }
         */
        Utils.Snack_UpdateExtremes(s, 0, s.length, Sound.SNACK_NEW_SOUND);
        Utils.Snack_ExecCallbacks(s, Sound.SNACK_NEW_SOUND);

        return 0;
    }

    /*
     * 
     * The following functions add interoperability between Snack and the CSLU
     * Speech Toolkit. Two functions are provided to convert Snack sound objects
     * into CSLUsh wave objects and vice versa.
     */

// #ifdef SNACK_CSLU_TOOLKIT
    int fromCSLUshWaveCmd(Sound s, Tcl.Interp interp, int objc, Tcl.Obj objv[]) {
        WAVE w;
        String handle;

        if (s.storeType != Sound.SOUND_IN_MEMORY) {
            Tcl.AppendResult(interp, "fromCSLUshWave only works with in-memory sounds", (String) null);
            return -1;
        }

        if (objc != 3) {
            Tcl.WrongNumArgs(interp, 1, objv, "fromCSLUshWave waveObj");
            return -1;
        }

        handle = Tcl.GetStringFromObj(objv[2], null);
        if (!(w = Obj_GetData(interp, WAVE, handle))) {
            Tcl.AppendResult(interp, "Failed getting data from waveObj: ", handle, null);
            return -1;
        }
        if (w.attr[WAVE_TYPE] != WAVE_TYPE_LINEAR) {
            Tcl.AppendResult(interp, "waveObj must be WAVE_TYPE_LINEAR", null);
            return -1;
        }
        if (s.writeStatus == Sound.WRITE) {
            Snack_StopSound(s, interp);
        }
        s.samprate = (int) w.attr[WAVE_RATE];
        s.encoding[0] = Sound.LIN16;
        s.sampsize[0] = 2; // TODO sizeof(short)
        s.nchannels[0] = 1;
        s.length = w.len;
        if (Utils.Snack_ResizeSoundStorage(s, s.length) != 0) {
            return -1;
        }
        Snack_PutSoundData(s, 0, w.samples, w.len * s.sampsize);
        Utils.Snack_UpdateExtremes(s, 0, s.length, Sound.SNACK_NEW_SOUND);
        Utils.Snack_ExecCallbacks(s, Sound.SNACK_NEW_SOUND);

        return 0;
    }

    int toCSLUshWaveCmd(Sound s, Tcl.Interp interp, int objc, Tcl.Obj objv[]) {
        WAVE w;
        float[] attr = new float[WAVE_ATTRIBUTES];

        if (s.storeType != Sound.SOUND_IN_MEMORY) {
            Tcl.AppendResult(interp, "toCSLUshWave only works with in-memory sounds", (String) null);
            return -1;
        }

        if (objc != 2) {
            Tcl.WrongNumArgs(interp, 1, objv, "toCSLUshWave");
            return -1;
        }
        if (s.encoding[0] != Sound.LIN16 || s.nchannels[0] != 1) {
            Tcl.AppendResult(interp, "Sorry, only implemented for lin16, mono sounds", null);
            return -1;
        }

        attr[WAVE_RATE] = s.samprate;
        attr[WAVE_TYPE] = 0;
        /* Doesn't handle large sounds yet > 512kB */
        if (createWave(interp, s.length, WAVE_ATTRIBUTES, attr, s.blocks[0], w) != 0)
            return -1;

        return 0;
    }
// #endif
}
