/*
 * Copyright (c) 2008 by Naohide Sano, All rights reserved.
 *
 * Programmed by Naohide Sano
 */

package se.kth.snack;


/**
 * Tcl. 
 *
 * @author <a href="mailto:vavivavi@yahoo.co.jp">Naohide Sano</a> (nsano)
 * @version 0.00 2008/07/13 nsano initial version <br>
 */
public class Tcl {

    class Obj {
        
    }

    class Interp {
        
    }

    class Channel {
        
    }

    /** */
    public static void AppendResult(Interp interp, Object... args) {
    }

    /** */
    public static void SetObjResult(Interp interp, Obj list) {
        // TODO Auto-generated method stub
        
    }

    /** */
    public static Obj NewDoubleObj(double d) {
        // TODO Auto-generated method stub
        return null;
    }

    /** */
    public static void ListObjAppendElement(Interp interp, Obj list, Object value) {
        // TODO Auto-generated method stub
        
    }

    /** */
    public static Obj NewListObj(int i, Object object) {
        // TODO Auto-generated method stub
        return null;
    }

    /** */
    public static int GetIndexFromObj(Interp interp, Obj obj, String[] strings, String string, int i, int index) {
        // TODO Auto-generated method stub
        return 0;
    }

    /** */
    public static int GetIntFromObj(Interp interp, Obj obj, int value) {
        // TODO Auto-generated method stub
        return 0;
    }

    /** */
    public static String GetStringFromObj(Obj obj, Object object) {
        // TODO Auto-generated method stub
        return null;
    }

    /** */
    public static void IncrRefCount(Obj obj) {
        // TODO Auto-generated method stub
        
    }

    /** */
    public static void DecrRefCount(Obj obj) {
        // TODO Auto-generated method stub
        
    }

    /** */
    public static int GetDoubleFromObj(Interp interp, Obj obj, double value) {
        // TODO Auto-generated method stub
        return 0;
    }

    /** */
    public static Obj NewStringObj(String type, int i) {
        // TODO Auto-generated method stub
        return null;
    }

    /** */
    public static int GlobalEvalObj(Interp interp, Obj cmd) {
        // TODO Auto-generated method stub
        return 0;
    }

    /** */
    public static Obj NewIntObj(int skipBytes) {
        // TODO Auto-generated method stub
        return null;
    }

    /** */
    public static void DeleteCommand(Interp interp, String string) {
        // TODO Auto-generated method stub
        
    }

    /** */
    public static void AddErrorInfo(Interp interp, String string) {
        // TODO Auto-generated method stub
        
    }

    /** */
    public static void BackgroundError(Interp interp) {
        // TODO Auto-generated method stub
        
    }

    /**  */
    public static void WrongNumArgs(Interp interp, int i, Obj[] objv, String string) {
        // TODO Auto-generated method stub
        
    }

    /** */
    public static int GetBooleanFromObj(Interp interp, Obj obj, int guessProps) {
        // TODO Auto-generated method stub
        return 0;
    }
}

/* */
