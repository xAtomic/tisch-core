/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 1.3.40
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */

package libtisch;

public class Motion extends FeatureVector {
  private long swigCPtr;

  protected Motion(long cPtr, boolean cMemoryOwn) {
    super(libtischJNI.SWIGMotionUpcast(cPtr), cMemoryOwn);
    swigCPtr = cPtr;
  }

  protected static long getCPtr(Motion obj) {
    return (obj == null) ? 0 : obj.swigCPtr;
  }

  protected void finalize() {
    delete();
  }

  public synchronized void delete() {
    if (swigCPtr != 0) {
      if (swigCMemOwn) {
        swigCMemOwn = false;
        libtischJNI.delete_Motion(swigCPtr);
      }
      swigCPtr = 0;
    }
    super.delete();
  }

  public Motion(long tf) {
    this(libtischJNI.new_Motion__SWIG_0(tf), true);
  }

  public Motion() {
    this(libtischJNI.new_Motion__SWIG_1(), true);
  }

  public FeatureBase clone() {
    long cPtr = libtischJNI.Motion_clone(swigCPtr, this);
    return (cPtr == 0) ? null : new Motion(cPtr, false);
  }

  public void load(InputState state) {
    libtischJNI.Motion_load(swigCPtr, this, InputState.getCPtr(state), state);
  }

  public String name() {
    return libtischJNI.Motion_name(swigCPtr, this);
  }

  public static Motion dynamic_cast(FeatureBase base) {
    long cPtr = libtischJNI.Motion_dynamic_cast(FeatureBase.getCPtr(base), base);
    return (cPtr == 0) ? null : new Motion(cPtr, false);
  }

}
