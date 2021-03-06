/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 1.3.40
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */

package libtisch;

public class MatcherTUIOInput {
  private long swigCPtr;
  protected boolean swigCMemOwn;

  protected MatcherTUIOInput(long cPtr, boolean cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = cPtr;
  }

  protected static long getCPtr(MatcherTUIOInput obj) {
    return (obj == null) ? 0 : obj.swigCPtr;
  }

  protected void finalize() {
    delete();
  }

  public synchronized void delete() {
    if (swigCPtr != 0) {
      if (swigCMemOwn) {
        swigCMemOwn = false;
        libtischJNI.delete_MatcherTUIOInput(swigCPtr);
      }
      swigCPtr = 0;
    }
  }

  public MatcherTUIOInput(Matcher m) {
    this(libtischJNI.new_MatcherTUIOInput(Matcher.getCPtr(m), m), true);
  }

  public void process_frame() {
    libtischJNI.MatcherTUIOInput_process_frame(swigCPtr, this);
  }

  public void process_blob(BasicBlob b) {
    libtischJNI.MatcherTUIOInput_process_blob(swigCPtr, this, BasicBlob.getCPtr(b), b);
  }

}
