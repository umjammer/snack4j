/*
 * Copyright (c) 2008 by Naohide Sano, All rights reserved.
 *
 * Programmed by Naohide Sano
 */

package se.kth.snack;

import java.util.List;
import java.util.Map;

import se.kth.snack.Formant.POLE;
import se.kth.snack.Tcl.Obj;


/**
 * Sound.
 * 
 * @author <a href="mailto:vavivavi@yahoo.co.jp">Naohide Sano</a> (nsano)
 * @version 0.00 080707 nsano initial version <br>
 */
public class Sound {
    static final int  SNACK_MONO =  1;
    static final int SNACK_STEREO=2;
    static final int  SNACK_QUAD  = 4;

    static final int  SNACK_SOUND_CMD =1;
    static final int SNACK_AUDIO_CMD =2;
    static final int  SNACK_MIXER_CMD= 3;

    enum SnackEndianness {
        SNACK_NATIVE,
        SNACK_BIGENDIAN,
        SNACK_LITTLEENDIAN
    }

    static final int IDLE = 0;
    static final int READ = 1;
    static final int WRITE = 2;
    static final int PAUSED = 3;
    static final int MAXNBLKS = 200;
    static final int FEXP = 17;
    static final int CEXP = (FEXP + 2);
    static final int DEXP = (FEXP - 1);
    static final int FBLKSIZE = (1 << FEXP);
    static final int DBLKSIZE = (1 << DEXP);
    static final int CBLKSIZE = (1 << CEXP);
    static final int SOUND_IN_MEMORY = 0;
    static final int SOUND_IN_CHANNEL = 1;
    static final int SOUND_IN_FILE = 2;
    /* NFIRSTSAMPLES */
    static final int CHANNEL_HEADER_BUFFER = 20000;
    static final int SNACK_NEW_SOUND = 1;
    static final int SNACK_MORE_SOUND = 2;
    static final int SNACK_DESTROY_SOUND = 3;
    static final int SNACK_SINGLE_PREC = 1;
    static final int SNACK_DOUBLE_PREC = 2;
    static final int LIN16 = 1;
    static final int ALAW = 2;
    static final int MULAW = 3;
    static final int LIN8OFFSET = 4;
    static final int LIN8 = 5;
    static final int LIN24 = 6;
    static final int LIN32 = 7;
    static final int SNACK_FLOAT = 8;
    static final int SNACK_DOUBLE = 9;
    static final int LIN24PACKED = 10;
    static final String QUE_STRING = "QUE";
    static final String RAW_STRING = "RAW";
    static final String WAV_STRING = "WAV";
    static final String AIFF_STRING = "AIFF";
    static final String SMP_STRING = "SMP";
    static final String AU_STRING = "AU";
    static final String SD_STRING = "SD";
    static final String MP3_STRING = "MP3";
    static final String CSL_STRING = "CSL";

    static float FSAMPLE(Sound s, int i) { return (Float) (s.blocks[i >> FEXP][i & (FBLKSIZE - 1)]); }
    static void FSAMPLE(Sound s, int i, float v) { s.blocks[i >> FEXP][i & (FBLKSIZE - 1)] = v; }
    static double DSAMPLE(Sound s, int i) { return (Double) (s.blocks[i >> DEXP][i & (DBLKSIZE - 1)]); }
    static void DSAMPLE(Sound s, int i, double v) { s.blocks[i >> DEXP][i & (DBLKSIZE - 1)] = v; }

    int samprate;
    int[] encoding = new int[1];
    int[] sampsize = new int[1];
    int[] nchannels = new int[1];
    int length;
    int maxlength;
    float maxsamp;
    float minsamp;
    float abmax;
    Object[][] blocks;
    int maxblks;
    int nblks;
    int exact;
    int precision;
    int writeStatus;
    int readStatus;
    short[] tmpbuf;
    int swap;
    int storeType;
    int headSize;
    int skipBytes;
    int buffersize;
    Tcl.Interp     interp;
    Tcl.Obj        cmdPtr;
    String fcname;
    List<Callback> firstCB;
    String fileType;
    int blockingPlay;
    int debug;
    int destroy;
    int guessEncoding;
    SnackEndianness inByteOrder;
    int firstNRead;
    int guessRate;
    int forceFormat;
    int itemRefCnt;
    int validStart;
    SnackLinkedFileInfo linkInfo;
    String devStr;
    Map soundTable;
    String filterName;
    POLE[] extHead;
    byte[] extHead2;
    int extHeadType;
    int extHead2Type;
    int loadOffset;
    Obj changeCmdPtr;
    /* User flags, for new file formats, etc */
    int userFlag;
    /* User data pointer */
    byte[] userData;
}

class SnackLinkedFileInfo {
    float[] buffer;
    int filePos;
    int validSamples;
    int eof;
    Sound sound;
    Tcl.Channel linkCh;
}

class Callback {

    int id;
    Object clientData;
    /** */
    public void exec(Object clientData, int flag) {
        // TODO Auto-generated method stub
        
    }
}

/* */
