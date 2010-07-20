/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 1.3.40
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */

package libtisch;

public class Label extends Widget {
  private long swigCPtr;

  protected Label(long cPtr, boolean cMemoryOwn) {
    super(libtischJNI.SWIGLabelUpcast(cPtr), cMemoryOwn);
    swigCPtr = cPtr;
  }

  protected static long getCPtr(Label obj) {
    return (obj == null) ? 0 : obj.swigCPtr;
  }

  protected void finalize() {
    delete();
  }

  public synchronized void delete() {
    if (swigCPtr != 0) {
      if (swigCMemOwn) {
        swigCMemOwn = false;
        libtischJNI.delete_Label(swigCPtr);
      }
      swigCPtr = 0;
    }
    super.delete();
  }

  protected void swigDirectorDisconnect() {
    swigCMemOwn = false;
    delete();
  }

  public void swigReleaseOwnership() {
    swigCMemOwn = false;
    libtischJNI.Label_change_ownership(this, swigCPtr, false);
  }

  public void swigTakeOwnership() {
    swigCMemOwn = true;
    libtischJNI.Label_change_ownership(this, swigCPtr, true);
  }

  public Label(String text, int _w, int _h, int _x, int _y, double angle, int center, int snip, RGBATexture _tex) {
    this(libtischJNI.new_Label__SWIG_0(text, _w, _h, _x, _y, angle, center, snip, libtischJNI.getCPtrAddRef(_tex), _tex), true);
    libtischJNI.Label_director_connect(this, swigCPtr, swigCMemOwn, true);
  }

  public Label(String text, int _w, int _h, int _x, int _y, double angle, int center, int snip) {
    this(libtischJNI.new_Label__SWIG_1(text, _w, _h, _x, _y, angle, center, snip), true);
    libtischJNI.Label_director_connect(this, swigCPtr, swigCMemOwn, true);
  }

  public Label(String text, int _w, int _h, int _x, int _y, double angle, int center) {
    this(libtischJNI.new_Label__SWIG_2(text, _w, _h, _x, _y, angle, center), true);
    libtischJNI.Label_director_connect(this, swigCPtr, swigCMemOwn, true);
  }

  public Label(String text, int _w, int _h, int _x, int _y, double angle) {
    this(libtischJNI.new_Label__SWIG_3(text, _w, _h, _x, _y, angle), true);
    libtischJNI.Label_director_connect(this, swigCPtr, swigCMemOwn, true);
  }

  public Label(String text, int _w, int _h, int _x, int _y) {
    this(libtischJNI.new_Label__SWIG_4(text, _w, _h, _x, _y), true);
    libtischJNI.Label_director_connect(this, swigCPtr, swigCMemOwn, true);
  }

  public Label(String text, int _w, int _h, int _x) {
    this(libtischJNI.new_Label__SWIG_5(text, _w, _h, _x), true);
    libtischJNI.Label_director_connect(this, swigCPtr, swigCMemOwn, true);
  }

  public Label(String text, int _w, int _h) {
    this(libtischJNI.new_Label__SWIG_6(text, _w, _h), true);
    libtischJNI.Label_director_connect(this, swigCPtr, swigCMemOwn, true);
  }

  public void draw() {
    if (getClass() == Label.class) libtischJNI.Label_draw(swigCPtr, this); else libtischJNI.Label_drawSwigExplicitLabel(swigCPtr, this);
  }

  public void action(Gesture gesture) {
    if (getClass() == Label.class) libtischJNI.Label_action(swigCPtr, this, Gesture.getCPtr(gesture), gesture); else libtischJNI.Label_actionSwigExplicitLabel(swigCPtr, this, Gesture.getCPtr(gesture), gesture);
  }

  public void set(String _text) {
    libtischJNI.Label_set(swigCPtr, this, _text);
  }

  protected void setText(String value) {
    libtischJNI.Label_text_set(swigCPtr, this, value);
  }

  protected String getText() {
    return libtischJNI.Label_text_get(swigCPtr, this);
  }

  protected void setCenter(int value) {
    libtischJNI.Label_center_set(swigCPtr, this, value);
  }

  protected int getCenter() {
    return libtischJNI.Label_center_get(swigCPtr, this);
  }

  protected void setSnip(int value) {
    libtischJNI.Label_snip_set(swigCPtr, this, value);
  }

  protected int getSnip() {
    return libtischJNI.Label_snip_get(swigCPtr, this);
  }

}