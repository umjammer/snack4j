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

import java.nio.ByteBuffer;
import java.util.HashMap;
import java.util.Map;


class Utils {

    int Snack_AddCallback(Sound s, Callback proc, Object cd) {
        Callback cb = new Callback();

        proc.clientData = cd;
        proc.id = s.firstCB.size();

        if (s.debug > 1) {
            System.err.printf("  Snack_AddCallback", cb.id);
        }

        return cb.id;
    }

    static void Snack_RemoveCallback(Sound s, int id) {
        Callback cbGoner = null;

        if (s.debug > 1) {
            System.err.printf("  Snack_RemoveCallback", id);
        }

        if (id == -1) {
            return;
        }

        for (Callback cb : s.firstCB) {
            if (cb.id == id) {
                cbGoner = cb;
                return;
            }
        }
    }

    static void Snack_ExecCallbacks(Sound s, int flag) {

        if (s.debug > 1) {
            System.err.printf("  Enter Snack_ExecCallbacks\n");
        }

        for (Callback cb : s.firstCB) {
            if (s.debug > 2) {
                System.err.printf("    Executing callback", cb.id);
            }
            cb.exec(cb.clientData, flag);
            if (s.debug > 2) {
                System.err.printf("    callback done\n");
            }
        }

        if (s.changeCmdPtr != null) {
            Tcl.Obj cmd = null;

            cmd = Tcl.NewListObj(0, null);
            Tcl.ListObjAppendElement(s.interp, cmd, s.changeCmdPtr);

            if (flag == Sound.SNACK_NEW_SOUND) {
                Tcl.ListObjAppendElement(s.interp, cmd, Tcl.NewStringObj("New", -1));
            } else if (flag == Sound.SNACK_DESTROY_SOUND) {
                Tcl.ListObjAppendElement(s.interp, cmd, Tcl.NewStringObj("Destroyed", -1));
            } else {
                Tcl.ListObjAppendElement(s.interp, cmd, Tcl.NewStringObj("More", -1));
            }
            if (Tcl.GlobalEvalObj(s.interp, cmd) != 0) {
                Tcl.AddErrorInfo(s.interp, "\n    (\"command\" script)");
                Tcl.BackgroundError(s.interp);
            }
        }
    }

    static void Snack_GetExtremes(Sound s, SnackLinkedFileInfo info, int start, int end, int chan, float[] pmax, float[] pmin) {
        int i, inc;
        float maxs, mins;

        if (s.length == 0) {
            if (s.encoding[0] == Sound.LIN8OFFSET) {
                pmax[0] = 128.0f;
                pmin[0] = 128.0f;
            } else {
                pmax[0] = 0.0f;
                pmin[0] = 0.0f;
            }
            return;
        }

        if (chan == -1) {
            inc = 1;
            chan = 0;
        } else {
            inc = s.nchannels[0];
        }

        start = start * s.nchannels[0] + chan;
        end = end * s.nchannels[0] + chan;

        switch (s.encoding[0]) {
        case Sound.LIN8OFFSET:
            maxs = 0.0f;
            mins = 255.0f;
            break;
        case Sound.LIN8:
            maxs = -128.0f;
            mins = 127.0f;
            break;
        case Sound.LIN24:
        case Sound.LIN24PACKED:
            maxs = -8388608.0f;
            mins = 8388607.0f;
            break;
        case Sound.LIN32:
            maxs = -2147483648.0f;
            mins = 2147483647.0f;
            break;
        default:
            maxs = -32768.0f;
            mins = 32767.0f;
        }

        if (s.precision == Sound.SNACK_SINGLE_PREC) {
            if (s.storeType == Sound.SOUND_IN_MEMORY) {
                for (i = start; i <= end; i += inc) {
                    float tmp = Sound.FSAMPLE(s, i);
                    if (tmp > maxs) {
                        maxs = tmp;
                    }
                    if (tmp < mins) {
                        mins = tmp;
                    }
                }
            } else {
                for (i = start; i <= end; i += inc) {
                    float tmp = GetSample(info, i);
                    if (tmp > maxs) {
                        maxs = tmp;
                    }
                    if (tmp < mins) {
                        mins = tmp;
                    }
                }
            }
        } else {
            if (s.storeType == Sound.SOUND_IN_MEMORY) {
                for (i = start; i <= end; i += inc) {
                    float tmp = (float) Sound.DSAMPLE(s, i);
                    if (tmp > maxs) {
                        maxs = tmp;
                    }
                    if (tmp < mins) {
                        mins = tmp;
                    }
                }
            } else {
                for (i = start; i <= end; i += inc) {
                    float tmp = GetSample(info, i);
                    if (tmp > maxs) {
                        maxs = tmp;
                    }
                    if (tmp < mins) {
                        mins = tmp;
                    }
                }
            }
        }
        if (maxs < mins) {
            maxs = mins;
        }
        if (mins > maxs) {
            mins = maxs;
        }

        pmax[0] = maxs;
        pmin[0] = mins;
    }

    static void Snack_UpdateExtremes(Sound s, int start, int end, int flag) {
        float maxs, mins;
        float[] newmax = new float[1], newmin = new float[1];

        if (flag == Sound.SNACK_NEW_SOUND) {
            s.maxsamp = -32768.0f;
            s.minsamp = 32767.0f;
        }

        maxs = s.maxsamp;
        mins = s.minsamp;

        Snack_GetExtremes(s, null, start, end - 1, -1, newmax, newmin);

        if (newmax[0] > maxs) {
            s.maxsamp = newmax[0];
        } else {
            s.maxsamp = maxs;
        }
        if (newmin[0] < mins) {
            s.minsamp = newmin[0];
        } else {
            s.minsamp = mins;
        }
        if (s.maxsamp > -s.minsamp) {
            s.abmax = s.maxsamp;
        } else {
            s.abmax = -s.minsamp;
        }
    }

    static void Snack_DeleteSound(Sound s) {
        Snack_FileFormat ff;

        if (s.debug > 1) {
            System.err.printf("  Enter Snack_DeleteSound\n");
        }

        Snack_ResizeSoundStorage(s, 0);
        if (s.storeType == Sound.SOUND_IN_FILE && s.linkInfo.linkCh != null) {
            CloseLinkedFile(s.linkInfo);
        }

        for (ff = snackFileFormats; ff != null; ff = ff.nextPtr) {
            if (strcmp(s.fileType, ff.name) == 0) {
                if (ff.freeHeaderProc != null) {
                    ff.freeHeaderProc(s);
                }
            }
        }

        Snack_ExecCallbacks(s, Sound.SNACK_DESTROY_SOUND);
        for (Callback currCB : s.firstCB) {
            if (s.debug > 1) {
                System.err.printf("  Freed callback", currCB.id);
            }
        }

        if (s.changeCmdPtr != null) {
            Tcl.DecrRefCount(s.changeCmdPtr);
        }

        if (s.cmdPtr != null) {
            Tcl.DecrRefCount(s.cmdPtr);
        }

        if (s.debug > 1) {
            System.err.printf("  Sound object freed\n");
        }
    }

    static int Snack_ResizeSoundStorage(Sound s, int len) {
        int neededblks, i, blockSize, sampSize;

        if (s.debug > 1) {
            System.err.printf("  Enter ResizeSoundStorage", len);
        }

        if (s.precision == Sound.SNACK_SINGLE_PREC) {
            blockSize = Sound.FBLKSIZE;
            sampSize = 4;
        } else {
            blockSize = Sound.DBLKSIZE;
            sampSize = 8;
        }
        neededblks = 1 + (len * s.nchannels[0] - 1) / blockSize;

        if (len == 0) {
            neededblks = 0;
            s.exact = 0;
        }

        if (neededblks > s.maxblks) {
            ByteBuffer tmp = ByteBuffer.allocate(neededblks);
            s.maxblks = neededblks;
            s.blocks[0] = tmp.asFloatBuffer().array();
        }

        if (s.maxlength == 0 && len * s.nchannels[0] < blockSize) {

            // Allocate exactly as much as needed.

            if (s.debug > 2) {
                System.err.printf("    Allocating minimal block", len * s.nchannels[0] * 4);
            }

            s.exact = len * s.nchannels[0] * sampSize;
            s.blocks[0] = new Number[s.exact];
            i = 1;
            s.maxlength = len;
        } else if (neededblks > s.nblks) {
            Number[] tmp = s.blocks[0];

            if (s.debug > 2) {
                System.err.printf("    Allocating full block(s)", neededblks - s.nblks);
            }

            // Do not count exact block, needs to be re-allocated
            if (s.exact > 0) {
                s.nblks = 0;
            }

            for (i = s.nblks; i < neededblks; i++) {
                s.blocks[i] = new Number[Sound.CBLKSIZE];
            }
            if (i < neededblks) {
                if (s.debug > 2) {
                    System.err.printf("    block alloc failed", i);
                }
                return -1;
            }

            // Copy and de-allocate any exact block
            if (s.exact > 0) {
                System.arraycopy(tmp, 0, s.blocks, 0, s.exact);
                s.exact = 0;
            }

            s.maxlength = neededblks * blockSize / s.nchannels[0];
        } else if (neededblks == 1 && s.exact > 0) {

            // Reallocate to one full block

            Number[] tmp = new Number[Sound.CBLKSIZE];

            if (s.debug > 2) {
                System.err.printf("    Reallocating full block", Sound.CBLKSIZE);
            }

            if (tmp != null) {
                System.arraycopy(s.blocks, 0, tmp, 0, s.exact);
                s.blocks[0] = tmp;
                s.maxlength = blockSize / s.nchannels[0];
            }
            s.exact = 0;
        }

        if (neededblks < s.nblks) {
            s.maxlength = neededblks * blockSize / s.nchannels[0];
        }

        s.nblks = neededblks;

        if (s.debug > 1) {
            System.err.printf("  Exit ResizeSoundStorage", neededblks);
        }

        return 0;
    }

    static final String[] encs = {
        "", "Lin16", "Alaw", "Mulaw", "Lin8offset", "Lin8", "Lin24", "Lin32", "Float", "Double", "Lin24packed"
    };

    int GetChannels(Tcl.Interp interp, Tcl.Obj obj, int[] nchannels) {
        int length, val;
        String str = Tcl.GetStringFromObj(obj, length);

        if (str.equalsIgnoreCase("MONO")) {
            nchannels[0] = Sound.SNACK_MONO;
            return 0;
        }
        if (str.equalsIgnoreCase("STEREO")) {
            nchannels[0] = Sound.SNACK_STEREO;
            return 0;
        }
        if (str.equalsIgnoreCase("QUAD")) {
            nchannels[0] = Sound.SNACK_QUAD;
            return 0;
        }
        if (Tcl.GetIntFromObj(interp, obj, val) != 0)
            return -1;
        if (val < 1) {
            Tcl.AppendResult(interp, "Number of channels must be >= 1", null);
            return -1;
        }
        nchannels[0] = val;
        return 0;
    }

    int GetEncoding(Tcl.Interp interp, Tcl.Obj obj, int[] encoding, int[] sampsize) {
        int length;
        String str = Tcl.GetStringFromObj(obj, length);

        if (str.equalsIgnoreCase("LIN16")) {
            encoding[0] = Sound.LIN16;
            sampsize[0] = 2;
        } else if (str.equalsIgnoreCase("LIN24")) {
            encoding[0] = Sound.LIN24;
            sampsize[0] = 4;
        } else if (str.equalsIgnoreCase("LIN24PACKED")) {
            encoding[0] = Sound.LIN24PACKED;
            sampsize[0] = 3;
        } else if (str.equalsIgnoreCase("LIN32")) {
            encoding[0] = Sound.LIN32;
            sampsize[0] = 4;
        } else if (str.equalsIgnoreCase("FLOAT")) {
            encoding[0] = Sound.SNACK_FLOAT;
            sampsize[0] = 4;
        } else if (str.equalsIgnoreCase("DOUBLE")) {
            encoding[0] = Sound.SNACK_DOUBLE;
            sampsize[0] = 8;
        } else if (str.equalsIgnoreCase("ALAW")) {
            encoding[0] = Sound.ALAW;
            sampsize[0] = 1;
        } else if (str.equalsIgnoreCase("MULAW")) {
            encoding[0] = Sound.MULAW;
            sampsize[0] = 1;
        } else if (str.equalsIgnoreCase("LIN8")) {
            encoding[0] = Sound.LIN8;
            sampsize[0] = 1;
        } else if (str.equalsIgnoreCase("LIN8OFFSET")) {
            encoding[0] = Sound.LIN8OFFSET;
            sampsize[0] = 1;
        } else {
            Tcl.AppendResult(interp, "Unknown encoding", null);
            return -1;
        }
        return 0;
    }

    void SwapIfBE(Sound s) {
        if (Utils.littleEndian != 0) {
            s.swap = 0;
        } else {
            s.swap = 1;
        }
    }

    void SwapIfLE(Sound s) {
        if (Utils.littleEndian != 0) {
            s.swap = 1;
        } else {
            s.swap = 0;
        }
    }

    Command info = new Command() {
        public int exec(Sound s, Tcl.Interp interp, int objc, Tcl.Obj objv[]) {
            Tcl.Obj[] objs = new Tcl.Obj[8];

            objs[0] = Tcl.NewIntObj(s.length);
            objs[1] = Tcl.NewIntObj(s.samprate);
            if (s.encoding[0] == Sound.SNACK_FLOAT) {
                objs[2] = Tcl.NewDoubleObj(s.maxsamp);
                objs[3] = Tcl.NewDoubleObj(s.minsamp);
            } else {
                objs[2] = Tcl.NewIntObj((int) s.maxsamp);
                objs[3] = Tcl.NewIntObj((int) s.minsamp);
            }
            objs[4] = Tcl.NewStringObj(encs[s.encoding[0]], -1);
            objs[5] = Tcl.NewIntObj(s.nchannels[0]);
            objs[6] = Tcl.NewStringObj(s.fileType, -1);
            objs[7] = Tcl.NewIntObj(s.headSize);

            Tcl.SetObjResult(interp, Tcl.NewListObj(8, objs));
            return 0;
        }
    };

    enum subOptions {
        START,
        END,
        CHANNEL
    }

    Command max = new Command() {

        public int exec(Sound s, Tcl.Interp interp, int objc, Tcl.Obj objv[]) {
            int startpos = 0, endpos = s.length - 1, arg, channel = -1;
            float[] maxsamp = new float[1], minsamp = new float[1];
            SnackLinkedFileInfo info;
            final String[] subOptionStrings = {
                "-start", "-end", "-channel", null
            };

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
                case CHANNEL: {
                    String str = Tcl.GetStringFromObj(objv[arg + 1], null);
                    if (GetChannel(interp, str, s.nchannels[0], channel) != 0) {
                        return -1;
                        break;
                    }
                }
                }
            }
            if (endpos < 0) {
                endpos = s.length - 1;
            }

            if (startpos < 0 || (startpos >= s.length && startpos > 0)) {
                Tcl.AppendResult(interp, "Start value out of bounds", null);
                return -1;
            }
            if (endpos >= s.length) {
                Tcl.AppendResult(interp, "End value out of bounds", null);
                return -1;
            }

            if (objc == 2) {
                if (s.encoding[0] == Sound.SNACK_FLOAT) {
                    Tcl.SetObjResult(interp, Tcl.NewDoubleObj(s.maxsamp));
                } else {
                    Tcl.SetObjResult(interp, Tcl.NewIntObj((int) s.maxsamp));
                }
            } else {
                if (s.storeType != Sound.SOUND_IN_MEMORY) {
                    OpenLinkedFile(s, info);
                }
                Snack_GetExtremes(s, info, startpos, endpos, channel, maxsamp, minsamp);
                if (s.storeType != Sound.SOUND_IN_MEMORY) {
                    CloseLinkedFile(info);
                }
                if (s.encoding[0] == Sound.SNACK_FLOAT) {
                    Tcl.SetObjResult(interp, Tcl.NewDoubleObj(maxsamp[0]));
                } else {
                    Tcl.SetObjResult(interp, Tcl.NewIntObj((int) maxsamp[0]));
                }
            }

            return 0;
        }
    };

    enum subOptions2 {
        START,
        END,
        CHANNEL
    }

    Command min = new Command() {
        public int exec(Sound s, Tcl.Interp interp, int objc, Tcl.Obj objv[]) {
            int startpos = 0, endpos = s.length - 1, arg, channel = -1;
            float[] maxsamp = new float[1], minsamp = new float[1];
            SnackLinkedFileInfo info;
            final String[] subOptionStrings = {
                "-start", "-end", "-channel", null
            };

            for (arg = 2; arg < objc; arg += 2) {
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
                case CHANNEL: {
                    String str = Tcl.GetStringFromObj(objv[arg + 1], null);
                    if (GetChannel(interp, str, s.nchannels, channel) != 0) {
                        return -1;
                    }
                    break;
                }
                }
            }
            if (endpos < 0)
                endpos = s.length - 1;

            if (startpos < 0 || (startpos >= s.length && startpos > 0)) {
                Tcl.AppendResult(interp, "Start value out of bounds", null);
                return -1;
            }
            if (endpos >= s.length) {
                Tcl.AppendResult(interp, "End value out of bounds", null);
                return -1;
            }

            if (objc == 2) {
                if (s.encoding[0] == Sound.SNACK_FLOAT) {
                    Tcl.SetObjResult(interp, Tcl.NewDoubleObj(s.minsamp));
                } else {
                    Tcl.SetObjResult(interp, Tcl.NewIntObj((int) s.minsamp));
                }
            } else {
                if (s.storeType != Sound.SOUND_IN_MEMORY) {
                    OpenLinkedFile(s, info);
                }
                Snack_GetExtremes(s, info, startpos, endpos, channel, maxsamp, minsamp);
                if (s.storeType != Sound.SOUND_IN_MEMORY) {
                    CloseLinkedFile(info);
                }
                if (s.encoding[0] == Sound.SNACK_FLOAT) {
                    Tcl.SetObjResult(interp, Tcl.NewDoubleObj(minsamp[0]));
                } else {
                    Tcl.SetObjResult(interp, Tcl.NewIntObj((int) minsamp[0]));
                }
            }

            return 0;
        }
    };

    Command changed = new Command() {
        public int exec(Sound s, Tcl.Interp interp, int objc, Tcl.Obj objv[]) {
            if (objc != 3) {
                Tcl.WrongNumArgs(interp, 1, objv, "changed new|more");
                return -1;
            }
            if (s.storeType == Sound.SOUND_IN_MEMORY) {
                Snack_UpdateExtremes(s, 0, s.length, Sound.SNACK_NEW_SOUND);
            }
            if (objc > 2) {
                String string = Tcl.GetStringFromObj(objv[2], null);

                if (string.equalsIgnoreCase("new")) {
                    Snack_ExecCallbacks(s, Sound.SNACK_NEW_SOUND);
                    return 0;
                }
                if (string.equalsIgnoreCase("more")) {
                    Snack_ExecCallbacks(s, Sound.SNACK_MORE_SOUND);
                    return 0;
                }
                Tcl.AppendResult(interp, "unknow option, must be new or more", (String) null);
                return -1;
            }

            return 0;
        }
    };

    Command destroy = new Command() {
        public int exec(Sound s, Tcl.Interp interp, int objc, Tcl.Obj objv[]) {
            String string = Tcl.GetStringFromObj(objv[0], null);
            int debug = s.debug;

            if (debug > 0) {
                System.err.printf("Enter destroyCmd\n");
            }

            if (s.writeStatus == Sound.WRITE) {
                s.destroy = 1;
            }
            s.length = 0;
            if (wop == Sound.IDLE) {
                Snack_StopSound(s, interp);
            }
            s.soundTable.remove(string);

            Tcl.DeleteCommand(interp, string);

            // The sound command and associated Sound struct are now deallocated
            // because SoundDeleteCmd has been called as a result of
            // Tcl.DeleteCommand().

            if (debug > 0) {
                System.err.printf("Exit destroyCmd\n");
            }

            return 0;
        }
    };

    Command flush = new Command() {
        public int exec(Sound s, Tcl.Interp interp, int objc, Tcl.Obj objv[]) {
            if (s.storeType != Sound.SOUND_IN_MEMORY) {
                Tcl.AppendResult(interp, "flush only works with in-memory sounds", (String) null);
                return -1;
            }

            Snack_StopSound(s, interp);
            Snack_ResizeSoundStorage(s, 0);
            s.length = 0;
            s.maxsamp = 0.0f;
            s.minsamp = 0.0f;
            s.abmax = 0.0f;
            Snack_ExecCallbacks(s, Sound.SNACK_NEW_SOUND);

            return 0;
        }
    };

    enum options {
        OPTLOAD,
        OPTFILE,
        CHANNEL,
        RATE,
        FREQUENCY,
        CHANNELS,
        ENCODING,
        FORMAT,
        BYTEORDER,
        BUFFERSIZE,
        SKIPHEAD,
        GUESSPROPS,
        PRECISION,
        CHGCMD,
        FILEFORMAT,
        OPTDEBUG
    };

    Command configure = new Command() {
        public int exec(Sound s, Tcl.Interp interp, int objc, Tcl.Obj objv[]) {
            int arg, filearg = 0, newobjc;
            Tcl.Obj[] newobjv = null;
            final String[] optionStrings = {
                "-load", "-file", "-channel", "-rate", "-frequency", "-channels", "-encoding", "-format", "-byteorder", "-buffersize", "-skiphead", "-guessproperties", "-precision", "-changecommand", "-fileformat", "-debug", null
            };
            Snack_FileFormat ff;

            if (s.debug > 0) {
                System.err.printf("Enter configureCmd\n");
            }

            Snack_RemoveOptions(objc - 2, objv[2], optionStrings, newobjc, (Tcl.Obj[]) newobjv);
            if (newobjc > 0) {
                for (ff = snackFileFormats; ff != null; ff = ff.nextPtr) {
                    if (strcmp(s.fileType, ff.name) == 0) {
                        if (ff.configureProc != null) {
                            if (ff.configureProc(s, interp, objc, objv)) {
                                return 0;
                            }
                        }
                    }
                }
            }
            for (arg = 0; arg < newobjc; arg++) {
                Tcl.DecrRefCount(newobjv[arg]);
            }

            if (objc == 2) { // get all options
                Tcl.Obj[] objs = new Tcl.Obj[6];

                objs[0] = Tcl.NewIntObj(s.length);
                objs[1] = Tcl.NewIntObj(s.samprate);
                if (s.encoding[0] == Sound.SNACK_FLOAT) {
                    objs[2] = Tcl.NewDoubleObj(s.maxsamp);
                    objs[3] = Tcl.NewDoubleObj(s.minsamp);
                } else {
                    objs[2] = Tcl.NewIntObj((int) s.maxsamp);
                    objs[3] = Tcl.NewIntObj((int) s.minsamp);
                }
                objs[4] = Tcl.NewStringObj(encs[s.encoding[0]], -1);
                objs[5] = Tcl.NewIntObj(s.nchannels[0]);

                Tcl.SetObjResult(interp, Tcl.NewListObj(6, objs));

                return 0;
            } else if (objc == 3) { // get option
                int index;

                if (Tcl.GetIndexFromObj(interp, objv[2], optionStrings, "option", 0, index) != 0) {
                    return -1;
                }

                switch (options.values()[index]) {
                case OPTLOAD: {
                    if (s.storeType == Sound.SOUND_IN_MEMORY) {
                        Tcl.SetObjResult(interp, Tcl.NewStringObj(s.fcname, -1));
                    } else {
                        Tcl.SetObjResult(interp, Tcl.NewStringObj("", -1));
                    }
                    break;
                }
                case OPTFILE: {
                    if (s.storeType == Sound.SOUND_IN_FILE) {
                        Tcl.SetObjResult(interp, Tcl.NewStringObj(s.fcname, -1));
                    } else {
                        Tcl.SetObjResult(interp, Tcl.NewStringObj("", -1));
                    }
                    break;
                }
                case CHANNEL: {
                    if (s.storeType == Sound.SOUND_IN_CHANNEL) {
                        Tcl.SetObjResult(interp, Tcl.NewStringObj(s.fcname, -1));
                    } else {
                        Tcl.SetObjResult(interp, Tcl.NewStringObj("", -1));
                    }
                    break;
                }
                case RATE:
                case FREQUENCY: {
                    Tcl.SetObjResult(interp, Tcl.NewIntObj(s.samprate));
                    break;
                }
                case CHANNELS: {
                    Tcl.SetObjResult(interp, Tcl.NewIntObj(s.nchannels[0]));
                    break;
                }
                case ENCODING:
                case FORMAT: {
                    Tcl.SetObjResult(interp, Tcl.NewStringObj(encs[s.encoding[0]], -1));
                    break;
                }
                case BYTEORDER:
                    if (s.sampsize[0] > 1) {
                        if (littleEndian!=0) {
                            if (s.swap != 0) {
                                Tcl.SetObjResult(interp, Tcl.NewStringObj("bigEndian", -1));
                            } else {
                                Tcl.SetObjResult(interp, Tcl.NewStringObj("littleEndian", -1));
                            }
                        } else {
                            if (s.swap != 0) {
                                Tcl.SetObjResult(interp, Tcl.NewStringObj("littleEndian", -1));
                            } else {
                                Tcl.SetObjResult(interp, Tcl.NewStringObj("bigEndian", -1));
                            }
                        }
                    } else {
                        Tcl.SetObjResult(interp, Tcl.NewStringObj("", -1));
                    }
                    break;
                case BUFFERSIZE: {
                    Tcl.SetObjResult(interp, Tcl.NewIntObj(s.buffersize));
                    break;
                }
                case SKIPHEAD: {
                    Tcl.SetObjResult(interp, Tcl.NewIntObj(s.skipBytes));
                    break;
                }
                case GUESSPROPS:
                    break;
                case PRECISION: {
                    if (s.precision == Sound.SNACK_DOUBLE_PREC) {
                        Tcl.SetObjResult(interp, Tcl.NewStringObj("double", -1));
                    } else {
                        Tcl.SetObjResult(interp, Tcl.NewStringObj("single", -1));
                    }
                    break;
                }
                case CHGCMD: {
                    Tcl.SetObjResult(interp, s.changeCmdPtr);
                    break;
                }
                case FILEFORMAT: {
                    Tcl.SetObjResult(interp, Tcl.NewStringObj(s.fileType, -1));
                    break;
                }
                case OPTDEBUG: {
                    Tcl.SetObjResult(interp, Tcl.NewIntObj(s.debug));
                    break;
                }
                }
            } else { // set option

                s.guessEncoding = -1;
                s.guessRate = -1;

                for (arg = 2; arg < objc; arg += 2) {
                    int index;

                    if (Tcl.GetIndexFromObj(interp, objv[arg], optionStrings, "option", 0, index) != 0) {
                        return -1;
                    }

                    if (arg + 1 == objc) {
                        Tcl.AppendResult(interp, "No argument given for ", optionStrings[index], " option", (String) null);
                        return -1;
                    }

                    switch (options.values()[index]) {
                    case OPTLOAD: {
                        filearg = arg + 1;
                        s.storeType = Sound.SOUND_IN_MEMORY;
                        break;
                    }
                    case OPTFILE: {
                        filearg = arg + 1;
                        s.storeType = Sound.SOUND_IN_FILE;
                        break;
                    }
                    case CHANNEL: {
                        filearg = arg + 1;
                        s.storeType = Sound.SOUND_IN_CHANNEL;
                        break;
                    }
                    case RATE:
                    case FREQUENCY: {
                        if (Tcl.GetIntFromObj(interp, objv[arg + 1], s.samprate) != 0)
                            return -1;
                        s.guessRate = 0;
                        break;
                    }
                    case CHANNELS: {
                        int oldn = s.nchannels[0];

                        if (GetChannels(interp, objv[arg + 1], s.nchannels) != 0)
                            return -1;
                        if (oldn != s.nchannels[0]) {
                            s.length = s.length * oldn / s.nchannels[0];
                        }
                        break;
                    }
                    case ENCODING:
                    case FORMAT: {
                        if (GetEncoding(interp, objv[arg + 1], s.encoding, s.sampsize) != 0) {
                            return -1;
                        }
                        s.guessEncoding = 0;
                        break;
                    }
                    case BYTEORDER: {
                        int length;
                        String str = Tcl.GetStringFromObj(objv[arg + 1], length);
                        if (str.equalsIgnoreCase("littleEndian")) {
                            SwapIfBE(s);
                        } else if (str.equalsIgnoreCase("bigEndian")) {
                            SwapIfLE(s);
                        } else {
                            Tcl.AppendResult(interp, "-byteorder option should be bigEndian", " or littleEndian", null);
                            return -1;
                        }
                        s.guessEncoding = 0;
                        break;
                    }
                    case BUFFERSIZE: {
                        if (Tcl.GetIntFromObj(interp, objv[arg + 1], s.buffersize) != 0)
                            return -1;
                        break;
                    }
                    case SKIPHEAD: {
                        if (Tcl.GetIntFromObj(interp, objv[arg + 1], s.skipBytes) != 0)
                            return -1;
                        break;
                    }
                    case GUESSPROPS: {
                        int guessProps;
                        if (Tcl.GetBooleanFromObj(interp, objv[arg + 1], guessProps) != 0)
                            return -1;
                        if (guessProps != 0) {
                            if (s.guessEncoding == -1)
                                s.guessEncoding = 1;
                            if (s.guessRate == -1)
                                s.guessRate = 1;
                        }
                        break;
                    }
                    case PRECISION: {
                        int length;
                        String str = Tcl.GetStringFromObj(objv[arg + 1], length);
                        if (str.equalsIgnoreCase("double")) {
                            s.precision = Sound.SNACK_DOUBLE_PREC;
                        } else if (str.equalsIgnoreCase("single")) {
                            s.precision = Sound.SNACK_SINGLE_PREC;
                        } else {
                            Tcl.AppendResult(interp, "-precision option should be single", " or double", null);
                            return -1;
                        }
                        break;
                    }
                    case CHGCMD: {
                        if (s.changeCmdPtr != null) {
                            Tcl.DecrRefCount(s.changeCmdPtr);
                        }
                        s.changeCmdPtr = objv[arg + 1];
                        Tcl.IncrRefCount(s.changeCmdPtr);
                        break;
                    }
                    case FILEFORMAT: {
                        if (Tcl.GetStringFromObj(objv[arg + 1], null).length() > 0) {
                            if (GetFileFormat(interp, objv[arg + 1], s.fileType) != 0) {
                                return -1;
                            }
                            s.forceFormat = 1;
                        }
                        break;
                    }
                    case OPTDEBUG: {
                        if (arg + 1 == objc) {
                            Tcl.AppendResult(interp, "No debug flag given", null);
                            return -1;
                        }
                        if (Tcl.GetIntFromObj(interp, objv[arg + 1], s.debug) != 0) {
                            return -1;
                        }
                        break;
                    }
                    }
                }
                if (s.guessEncoding == -1) {
                    s.guessEncoding = 0;
                }
                if (s.guessRate == -1) {
                    s.guessRate = 0;
                }

                if (filearg > 0) {
                    if (SetFcname(s, interp, objv[filearg]) != 0) {
                        return -1;
                    }
                }

                if (filearg > 0 && s.fcname.length() > 0) {
                    if (s.storeType == Sound.SOUND_IN_MEMORY) {
                        String type = LoadSound(s, interp, null, 0, -1);

                        if (type == null) {
                            return -1;
                        }
                        Snack_UpdateExtremes(s, 0, s.length, Sound.SNACK_NEW_SOUND);
                    } else if (s.storeType == Sound.SOUND_IN_FILE) {
//                      Snack_FileFormat ff;

                        if (s.linkInfo.linkCh != null) {
                            CloseLinkedFile(s.linkInfo);
                            s.linkInfo.linkCh = null;
                        }
                        for (ff = snackFileFormats; ff != null; ff = ff.nextPtr) {
                            if (strcmp(s.fileType, ff.name) == 0) {
                                if (ff.freeHeaderProc != null) {
                                    ff.freeHeaderProc(s);
                                }
                            }
                        }
                        if (GetHeader(s, interp, null) != 0) {
                            s.fileType = NameGuessFileType(s.fcname);
                        }
                        Snack_ResizeSoundStorage(s, 0);
                        if (s.encoding[0] == Sound.LIN8OFFSET) {
                            s.maxsamp = 128.0f;
                            s.minsamp = 128.0f;
                        } else {
                            s.maxsamp = 0.0f;
                            s.minsamp = 0.0f;
                        }
                    } else if (s.storeType == Sound.SOUND_IN_CHANNEL) {
                        int mode = 0;

                        Snack_ResizeSoundStorage(s, 0);
                        s.rwchan = Tcl.GetChannel(interp, s.fcname, mode);
                        if (!(mode & Tcl.READABLE)) {
                            s.rwchan = null;
                        }
                        if (s.rwchan != null) {
                            Tcl.SetChannelOption(interp, s.rwchan, "-translation", "binary");
                            // #ifdef Tcl.81_API
                            // Tcl.SetChannelOption(interp, s.rwchan,
                            // "-encoding",
                            // "binary");
                            // #endif
                        }
                    }
                }
                if (filearg > 0 && s.fcname.length() == 0) {
                    if (s.storeType == Sound.SOUND_IN_FILE) {
                        s.length = 0;
                    }
                }
                Snack_ExecCallbacks(s, Sound.SNACK_NEW_SOUND);
            }
            if (s.debug > 0) {
                System.err.printf("Exit configureCmd\n");
            }

            return 0;
        }
    };

    enum options2 {
        OPTLOAD,
        OPTFILE,
        CHANNEL,
        RATE,
        FREQUENCY,
        CHANNELS,
        ENCODING,
        FORMAT,
        BYTEORDER,
        BUFFERSIZE,
        SKIPHEAD,
        GUESSPROPS,
        PRECISION,
        CHGCMD,
        FILEFORMAT,
        OPTDEBUG
    }

    Command cget = new Command() {
        public int exec(Sound s, Tcl.Interp interp, int objc, Tcl.Obj objv[]) {
            /* static */String[] optionStrings = {
                "-load", "-file", "-channel", "-rate", "-frequency", "-channels", "-encoding", "-format", "-byteorder", "-buffersize", "-skiphead", "-guessproperties", "-precision", "-changecommand", "-fileformat", "-debug", null
            };

            if (objc == 2) {
                Tcl.WrongNumArgs(interp, 1, objv, "cget option");
                return -1;
            } else if (objc == 3) { // get option
                int index;

                if (Tcl.GetIndexFromObj(interp, objv[2], optionStrings, "option", 0, index) != 0) {
                    return -1;
                }

                switch (options2.values()[index]) {
                case OPTLOAD: {
                    if (s.storeType == Sound.SOUND_IN_MEMORY) {
                        Tcl.SetObjResult(interp, Tcl.NewStringObj(s.fcname, -1));
                    } else {
                        Tcl.SetObjResult(interp, Tcl.NewStringObj("", -1));
                    }
                    break;
                }
                case OPTFILE: {
                    if (s.storeType == Sound.SOUND_IN_FILE) {
                        Tcl.SetObjResult(interp, Tcl.NewStringObj(s.fcname, -1));
                    } else {
                        Tcl.SetObjResult(interp, Tcl.NewStringObj("", -1));
                    }
                    break;
                }
                case CHANNEL: {
                    if (s.storeType == Sound.SOUND_IN_CHANNEL) {
                        Tcl.SetObjResult(interp, Tcl.NewStringObj(s.fcname, -1));
                    } else {
                        Tcl.SetObjResult(interp, Tcl.NewStringObj("", -1));
                    }
                    break;
                }
                case RATE:
                case FREQUENCY: {
                    Tcl.SetObjResult(interp, Tcl.NewIntObj(s.samprate));
                    break;
                }
                case CHANNELS: {
                    Tcl.SetObjResult(interp, Tcl.NewIntObj(s.nchannels[0]));
                    break;
                }
                case ENCODING:
                case FORMAT: {
                    Tcl.SetObjResult(interp, Tcl.NewStringObj(encs[s.encoding[0]], -1));
                    break;
                }
                case BYTEORDER:
                    if (s.sampsize[0] > 1) {
                        if (Utils.littleEndian != 0) {
                            if (s.swap != 0) {
                                Tcl.SetObjResult(interp, Tcl.NewStringObj("bigEndian", -1));
                            } else {
                                Tcl.SetObjResult(interp, Tcl.NewStringObj("littleEndian", -1));
                            }
                        } else {
                            if (s.swap != 0) {
                                Tcl.SetObjResult(interp, Tcl.NewStringObj("littleEndian", -1));
                            } else {
                                Tcl.SetObjResult(interp, Tcl.NewStringObj("bigEndian", -1));
                            }
                        }
                    } else {
                        Tcl.SetObjResult(interp, Tcl.NewStringObj("", -1));
                    }
                    break;
                case BUFFERSIZE: {
                    Tcl.SetObjResult(interp, Tcl.NewIntObj(s.buffersize));
                    break;
                }
                case SKIPHEAD: {
                    Tcl.SetObjResult(interp, Tcl.NewIntObj(s.skipBytes));
                    break;
                }
                case GUESSPROPS:
                    break;
                case CHGCMD: {
                    Tcl.SetObjResult(interp, s.changeCmdPtr);
                    break;
                }
                case PRECISION: {
                    if (s.precision == Sound.SNACK_DOUBLE_PREC) {
                        Tcl.SetObjResult(interp, Tcl.NewStringObj("double", -1));
                    } else {
                        Tcl.SetObjResult(interp, Tcl.NewStringObj("single", -1));
                    }
                    break;
                }
                case FILEFORMAT: {
                    Tcl.SetObjResult(interp, Tcl.NewStringObj(s.fileType, -1));
                    break;
                }
                case OPTDEBUG: {
                    Tcl.SetObjResult(interp, Tcl.NewIntObj(s.debug));
                    break;
                }
                }
            }

            return 0;
        }
    };

    /* NOTE: NSOUNDCOMMANDS needs updating when new commands are added. */

    /** */
    Map<String, Command> sndCmdProcs = new HashMap<String, Command>();

    /** */
    Map<String, Command> audioCmdProcs = new HashMap<String, Command>();

    /** */
    Map<String, Command> mixerCmdProcs = new HashMap<String, Command>();

    /** */
    int SoundCmd(Object clientData, Tcl.Interp interp, int objc, Tcl.Obj objv[]) {
        Sound s = (Sound) clientData;
        int index;

        if (objc < 2) {
            throw new IllegalArgumentException("option ?args?");
            return -1;
        }

        if (Tcl.GetIndexFromObj(interp, objv[1], sndCmdNames, "option", 0, index) != 0) {
            return -1;
        }

        return sndCmdProcs.get(index).exec(s, interp, objc, objv);
    }

    Sound Snack_NewSound(int rate, int encoding, int nchannels) {
        Sound s = new Sound();

        if (s == null) {
            return null;
        }

        // Default sound specifications

        s.samprate = rate;
        s.encoding[0] = encoding;
        if (s.encoding[0] == Sound.LIN16) {
            s.sampsize[0] = 2;
        } else if (s.encoding[0] == Sound.LIN24 || s.encoding[0] == Sound.LIN32 || s.encoding[0] == Sound.SNACK_FLOAT) {
            s.sampsize[0] = 4;
        } else if (s.encoding[0] == Sound.LIN24PACKED) {
            s.sampsize[0] = 3;
        } else {
            s.sampsize[0] = 1;
        }
        if (s.encoding[0] == Sound.LIN8OFFSET) {
            s.maxsamp = 128.0f;
            s.minsamp = 128.0f;
        } else {
            s.maxsamp = 0.0f;
            s.minsamp = 0.0f;
        }
        s.nchannels[0] = nchannels;
        s.length = 0;
        s.maxlength = 0;
        s.abmax = 0.0f;
        s.readStatus = Sound.IDLE;
        s.writeStatus = Sound.IDLE;
        s.firstCB = null;
        s.fileType = Sound.RAW_STRING;
        s.tmpbuf = null;
        s.swap = 0;
        s.headSize = 0;
        s.skipBytes = 0;
        s.storeType = Sound.SOUND_IN_MEMORY;
        s.fcname = null;
        s.interp = null;
        s.cmdPtr = null;
        s.blocks = new Number[Sound.MAXNBLKS][];
        if (s.blocks == null) {
            return null;
        }
        s.blocks[0] = null;
        s.maxblks = Sound.MAXNBLKS;
        s.nblks = 0;
        s.exact = 0;
        s.precision = Sound.SNACK_SINGLE_PREC;
        s.blockingPlay = 0;
        s.debug = 0;
        s.destroy = 0;
        s.guessEncoding = 0;
        s.guessRate = 0;
//        s.rwchan = null;
        s.firstNRead = 0;
        s.buffersize = 0;
        s.forceFormat = 0;
        s.itemRefCnt = 0;
        s.validStart = 0;
        s.linkInfo.linkCh = null;
        s.linkInfo.eof = 0;
        s.inByteOrder = Sound.SnackEndianness.SNACK_NATIVE;
        s.devStr = null;
        s.soundTable = null;
        s.filterName = null;
        s.extHead = null;
        s.extHeadType = 0;
        s.extHead2 = null;
        s.extHead2Type = 0;
        s.loadOffset = 0;
        s.changeCmdPtr = null;
        s.userFlag = 0;
        s.userData = null;

        return s;
    }

    void CleanSound(Sound s, Tcl.Interp interp, String name) {
        Snack_DeleteSound(s);
        s.soundTable.remove(name);
    }

    enum options3 {
        OPTLOAD,
        OPTFILE,
        RATE,
        FREQUENCY,
        CHANNELS,
        ENCODING,
        FORMAT,
        CHANNEL,
        BYTEORDER,
        BUFFERSIZE,
        SKIPHEAD,
        GUESSPROPS,
        FILEFORMAT,
        PRECISION,
        CHGCMD,
        OPTDEBUG
    }

    int ParseSoundCmd(Object cdata, Tcl.Interp interp, int objc, Tcl.Obj objv[], String namep, Sound sp) {
        Sound s;
        int arg, arg1, filearg = 0, flag;
        /* static */int id = 0;
        int samprate = defaultSampleRate;
        int[] nchannels = { 1 };
        int[] encoding = { Sound.LIN16 }, sampsize = { 2 };
        int storeType = -1, guessEncoding = -1, guessRate = -1;
        int forceFormat = -1, skipBytes = -1, buffersize = -1;
        int guessProps = 0, swapIfBE = -1, debug = -1, precision = -1;
        String fileType = null;
        /* static */String ids;
        String name;
        Map hTab = (Map) cdata;
        Object hPtr;
        int length = 0;
        String string = null;
        Tcl.Obj cmdPtr = null;
        final String[] optionStrings = {
            "-load", "-file", "-rate", "-frequency",
            "-channels", "-encoding", "-format", "-channel",
            "-byteorder", "-buffersize", "-skiphead", "-guessproperties",
            "-fileformat", "-precision", "-changecommand", "-debug", null
        };

        if (objc > 1) {
            string = Tcl.GetStringFromObj(objv[1], length);
        }
        if ((objc == 1) || (string.charAt(0) == '-')) {
            do {
                ids = String.format("sound%d", ++id);
            } while (hTab.get(ids) != null);
            name = ids;
            arg1 = 1;
        } else {
            name = string;
            arg1 = 2;
        }
        namep = name;

        hPtr = hTab.get(name);
        if (hPtr != null) {
            Sound t = (Sound) hPtr;
            Snack_StopSound(t, interp);
            Tcl.DeleteCommand(interp, name);
        }

        for (arg = arg1; arg < objc; arg += 2) {
            int index;

            if (Tcl.GetIndexFromObj(interp, objv[arg], optionStrings, "option", 0, index) != 0) {
                return -1;
            }

            if (arg + 1 == objc) {
                Tcl.AppendResult(interp, "No argument given for ", optionStrings[index], " option", (String) null);
                return -1;
            }

            switch (options3.values()[index]) {
            case OPTLOAD: {
                if (arg + 1 == objc) {
                    Tcl.AppendResult(interp, "No filename given", null);
                    return -1;
                }
                filearg = arg + 1;
                storeType = Sound.SOUND_IN_MEMORY;
                break;
            }
            case OPTFILE: {
                if (arg + 1 == objc) {
                    Tcl.AppendResult(interp, "No filename given", null);
                    return -1;
                }
                filearg = arg + 1;
                storeType = Sound.SOUND_IN_FILE;
                break;
            }
            case RATE:
            case FREQUENCY: {
                if (Tcl.GetIntFromObj(interp, objv[arg + 1], samprate) != 0) {
                    return -1;
                }
                guessRate = 0;
                break;
            }
            case CHANNELS: {
                if (GetChannels(interp, objv[arg + 1], nchannels) != 0) {
                    return -1;
                }
                break;
            }
            case ENCODING:
            case FORMAT: {
                if (GetEncoding(interp, objv[arg + 1], encoding, sampsize) != 0) {
                    return -1;
                }
                guessEncoding = 0;
                break;
            }
            case CHANNEL: {
                if (arg + 1 == objc) {
                    Tcl.AppendResult(interp, "No channel name given", null);
                    return -1;
                }
                filearg = arg + 1;
                storeType = Sound.SOUND_IN_CHANNEL;
                break;
            }
            case OPTDEBUG: {
                if (arg + 1 == objc) {
                    Tcl.AppendResult(interp, "No debug flag given", null);
                    return -1;
                }
                if (Tcl.GetIntFromObj(interp, objv[arg + 1], debug) != 0) {
                    return -1;
                }
                break;
            }
            case FILEFORMAT: {
                if (Tcl.GetStringFromObj(objv[arg + 1], null).length() > 0) {
                    if (GetFileFormat(interp, objv[arg + 1], fileType) != 0) {
                        return -1;
                    }
                    forceFormat = 1;
                }
                break;
            }
            case BYTEORDER: {
                String str = Tcl.GetStringFromObj(objv[arg + 1], length);
                if (str.equalsIgnoreCase("littleEndian")) {
                    swapIfBE = 1;
                } else if (str.equalsIgnoreCase("bigEndian")) {
                    swapIfBE = 0;
                } else {
                    Tcl.AppendResult(interp, "-byteorder option should be bigEndian or littleEndian", null);
                    return -1;
                }
                guessEncoding = 0;
                break;
            }
            case BUFFERSIZE: {
                if (Tcl.GetIntFromObj(interp, objv[arg + 1], buffersize) != 0)
                    return -1;
                break;
            }

            case SKIPHEAD: {
                if (Tcl.GetIntFromObj(interp, objv[arg + 1], skipBytes) != 0)
                    return -1;
                break;
            }
            case GUESSPROPS: {
                if (Tcl.GetBooleanFromObj(interp, objv[arg + 1], guessProps) != 0)
                    return -1;
                break;
            }
            case PRECISION: {
                String str = Tcl.GetStringFromObj(objv[arg + 1], length);
                if (str.equalsIgnoreCase("double")) {
                    precision = Sound.SNACK_DOUBLE_PREC;
                } else if (str.equalsIgnoreCase("single")) {
                    precision = Sound.SNACK_SINGLE_PREC;
                } else {
                    Tcl.AppendResult(interp, "-precision option should be single", " or double", null);
                    return -1;
                }
                break;
            }
            case CHGCMD: {
                String str = Tcl.GetStringFromObj(objv[arg + 1], null);

                if (str.length() > 0) {
                    cmdPtr = objv[arg + 1];
                    Tcl.IncrRefCount(cmdPtr);
                }
                break;
            }
            }
        }

        sp = s = Snack_NewSound(samprate, encoding, nchannels);

        hPtr = hTab.put(name, flag);
        hTab.put(hPtr, s);
        s.soundTable = hTab;

        if (guessProps != 0) {
            if (guessEncoding == -1) {
                s.guessEncoding = 1;
            }
            if (guessRate == -1) {
                s.guessRate = 1;
            }
        }
        if (storeType != -1) {
            s.storeType = storeType;
        }
        if (buffersize != -1) {
            s.buffersize = buffersize;
        }
        if (skipBytes != -1) {
            s.skipBytes = skipBytes;
        }
        if (debug != -1) {
            s.debug = debug;
        }
        if (fileType != null) {
            s.fileType = fileType;
        }
        if (forceFormat != -1) {
            s.forceFormat = forceFormat;
        }
        if (precision != -1) {
            s.precision = precision;
        }
        if (swapIfBE == 0) {
            SwapIfLE(s);
        }
        if (swapIfBE == 1) {
            SwapIfBE(s);
        }
        if (cmdPtr != null) {
            s.changeCmdPtr = cmdPtr;
        }

        /* s.fcname = strdup(name); */
        s.interp = interp;

        if (filearg > 0) {
            if (SetFcname(s, interp, objv[filearg]) != 0) {
                CleanSound(s, interp, name);
                return -1;
            }
        }

        if (filearg > 0 && s.fcname.length() > 0) {
            if (s.storeType == Sound.SOUND_IN_MEMORY) {
                String type = LoadSound(s, interp, null, 0, -1);

                if (type == null) {
                    CleanSound(s, interp, name);
                    return -1;
                }
                Snack_UpdateExtremes(s, 0, s.length, Sound.SNACK_NEW_SOUND);
            } else if (s.storeType == Sound.SOUND_IN_FILE) {
                if (GetHeader(s, interp, null) != 0) {
                    s.fileType = NameGuessFileType(s.fcname);
                }
                if (s.encoding[0] == Sound.LIN8OFFSET) {
                    s.maxsamp = 128.0f;
                    s.minsamp = 128.0f;
                } else {
                    s.maxsamp = 0.0f;
                    s.minsamp = 0.0f;
                }
            } else if (s.storeType == Sound.SOUND_IN_CHANNEL) {
                int mode = 0;

                s.rwchan = Tcl.GetChannel(interp, s.fcname, mode);
                if (!(mode & Tcl.READABLE)) {
                    s.rwchan = null;
                }
                if (s.rwchan != null) {
                    Tcl.SetChannelOption(interp, s.rwchan, "-translation", "binary");
                    // #ifdef Tcl.81_API
                    // Tcl.SetChannelOption(interp, s.rwchan, "-encoding",
                    // "binary");
                    // #endif
                }
            }
        }

        return 0;
    }

    private void SoundDeleteCmd(Object clientData) {
        Sound s = (Sound) clientData;
        int i;

        if (s.debug > 1) {
            System.err.printf("  Sound obj cmd deleted\n");
        }
        if (s.destroy == 0) {
            Snack_StopSound(s, s.interp);
        }
        if (s.destroy == 0 || wop == Sound.IDLE) {
            Snack_DeleteSound(s);
        }
    }

    int Snack_SoundCmd(Object cdata, Tcl.Interp interp, int objc, Tcl.Obj objv[]) {
        String name;
        Sound s = null;

        if (ParseSoundCmd(cdata, interp, objc, objv, name, s) != 0) {
            return -1;
        }

        Tcl.CreateObjCommand(interp, name, SoundCmd, (Object) s);

        Tcl.SetObjResult(interp, Tcl.NewStringObj(name, -1));

        return 0;
    }

    Sound Snack_GetSound(Tcl.Interp interp, String name) {
        Tcl.CmdInfo infoPtr;
        Object hPtr = filterHashTable.get(name);

        if (hPtr != null || Tcl.GetCommandInfo(interp, name, infoPtr) == 0) {
            Tcl.AppendResult(interp, name, " : no such sound", (String) null);
            return null;
        }

        return (Sound) infoPtr.objClientData;
    }

    void Snack_SoundDeleteCmd(Object clientData) {
    }

    int Snack_AddSubCmd(int snackCmd, String cmdName, Command cmdProc) {
        int i;

        switch (snackCmd) {
        case Sound.SNACK_SOUND_CMD:
            sndCmdProcs.put(cmdName, cmdProc);
            break;
        case Sound.SNACK_AUDIO_CMD:
            audioCmdProcs.put(cmdName, cmdProc);
            break;
        case Sound.SNACK_MIXER_CMD:
            mixerCmdProcs.put(cmdName, cmdProc);
            break;
        }

        return 0;
    }

    int SetFcname(Sound s, Tcl.Interp interp, Tcl.Obj obj) {
        int length;
        String str = Tcl.GetStringFromObj(obj, length);

        s.fcname = str;

        return 0;
    }

    static int Snack_ProgressCallback(Tcl.Obj cmdPtr, Tcl.Interp interp, String type, double fraction) {
        if (cmdPtr != null) {
            Tcl.Obj cmd = null;
            int res;

            cmd = Tcl.NewListObj(0, null);
            Tcl.ListObjAppendElement(interp, cmd, cmdPtr);
            Tcl.ListObjAppendElement(interp, cmd, Tcl.NewStringObj(type, -1));
            Tcl.ListObjAppendElement(interp, cmd, Tcl.NewDoubleObj(fraction));
//          Tcl.Preserve((ClientData) interp);
            res = Tcl.GlobalEvalObj(interp, cmd);
//          Tcl.Release((ClientData) interp);
            return res;
        }
        return 0;
    }

    int Snack_PlatformIsLittleEndian() {
        return littleEndian;
    }

    int wop;
    static Map filterHashTable;
    static int littleEndian = 0;
    static int defaultSampleRate = 16000;
}

interface Command {
    int exec(Sound s, Tcl.Interp interp, int objc, Tcl.Obj objv[]);
}

/* */
