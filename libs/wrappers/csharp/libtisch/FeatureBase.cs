/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 1.3.40
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */


using System;
using System.Runtime.InteropServices;

public class FeatureBase : IDisposable {
  private HandleRef swigCPtr;
  protected bool swigCMemOwn;

  internal FeatureBase(IntPtr cPtr, bool cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = new HandleRef(this, cPtr);
  }

  internal static HandleRef getCPtr(FeatureBase obj) {
    return (obj == null) ? new HandleRef(null, IntPtr.Zero) : obj.swigCPtr;
  }

  ~FeatureBase() {
    Dispose();
  }

  public virtual void Dispose() {
    lock(this) {
      if (swigCPtr.Handle != IntPtr.Zero) {
        if (swigCMemOwn) {
          swigCMemOwn = false;
          libtischPINVOKE.delete_FeatureBase(swigCPtr);
        }
        swigCPtr = new HandleRef(null, IntPtr.Zero);
      }
      GC.SuppressFinalize(this);
    }
  }

  public FeatureBase(uint _tf) : this(libtischPINVOKE.new_FeatureBase__SWIG_0(_tf), true) {
    if (libtischPINVOKE.SWIGPendingException.Pending) throw libtischPINVOKE.SWIGPendingException.Retrieve();
    SwigDirectorConnect();
  }

  public FeatureBase() : this(libtischPINVOKE.new_FeatureBase__SWIG_1(), true) {
    if (libtischPINVOKE.SWIGPendingException.Pending) throw libtischPINVOKE.SWIGPendingException.Retrieve();
    SwigDirectorConnect();
  }

  public virtual string name() {
    string ret = libtischPINVOKE.FeatureBase_name(swigCPtr);
    if (libtischPINVOKE.SWIGPendingException.Pending) throw libtischPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public virtual FeatureBase clone() {
    IntPtr cPtr = libtischPINVOKE.FeatureBase_clone(swigCPtr);
    FeatureBase ret = (cPtr == IntPtr.Zero) ? null : new FeatureBase(cPtr, false);
    if (libtischPINVOKE.SWIGPendingException.Pending) throw libtischPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public virtual void load(InputState state) {
    libtischPINVOKE.FeatureBase_load(swigCPtr, InputState.getCPtr(state));
    if (libtischPINVOKE.SWIGPendingException.Pending) throw libtischPINVOKE.SWIGPendingException.Retrieve();
  }

  public virtual int next() {
    int ret = libtischPINVOKE.FeatureBase_next(swigCPtr);
    if (libtischPINVOKE.SWIGPendingException.Pending) throw libtischPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public virtual void serialize(SWIGTYPE_p_std__ostream s) {
    libtischPINVOKE.FeatureBase_serialize(swigCPtr, SWIGTYPE_p_std__ostream.getCPtr(s));
    if (libtischPINVOKE.SWIGPendingException.Pending) throw libtischPINVOKE.SWIGPendingException.Retrieve();
  }

  public virtual void unserialize(SWIGTYPE_p_std__istream s) {
    libtischPINVOKE.FeatureBase_unserialize(swigCPtr, SWIGTYPE_p_std__istream.getCPtr(s));
    if (libtischPINVOKE.SWIGPendingException.Pending) throw libtischPINVOKE.SWIGPendingException.Retrieve();
  }

  public int has_result {
    set {
      libtischPINVOKE.FeatureBase_has_result_set(swigCPtr, value);
      if (libtischPINVOKE.SWIGPendingException.Pending) throw libtischPINVOKE.SWIGPendingException.Retrieve();
    } 
    get {
      int ret = libtischPINVOKE.FeatureBase_has_result_get(swigCPtr);
      if (libtischPINVOKE.SWIGPendingException.Pending) throw libtischPINVOKE.SWIGPendingException.Retrieve();
      return ret;
    } 
  }

  protected uint typeflags {
    set {
      libtischPINVOKE.FeatureBase_typeflags_set(swigCPtr, value);
      if (libtischPINVOKE.SWIGPendingException.Pending) throw libtischPINVOKE.SWIGPendingException.Retrieve();
    } 
    get {
      uint ret = libtischPINVOKE.FeatureBase_typeflags_get(swigCPtr);
      if (libtischPINVOKE.SWIGPendingException.Pending) throw libtischPINVOKE.SWIGPendingException.Retrieve();
      return ret;
    } 
  }

  private void SwigDirectorConnect() {
    if (SwigDerivedClassHasMethod("name", swigMethodTypes0))
      swigDelegate0 = new SwigDelegateFeatureBase_0(SwigDirectorname);
    if (SwigDerivedClassHasMethod("clone", swigMethodTypes1))
      swigDelegate1 = new SwigDelegateFeatureBase_1(SwigDirectorclone);
    if (SwigDerivedClassHasMethod("load", swigMethodTypes2))
      swigDelegate2 = new SwigDelegateFeatureBase_2(SwigDirectorload);
    if (SwigDerivedClassHasMethod("next", swigMethodTypes3))
      swigDelegate3 = new SwigDelegateFeatureBase_3(SwigDirectornext);
    if (SwigDerivedClassHasMethod("serialize", swigMethodTypes4))
      swigDelegate4 = new SwigDelegateFeatureBase_4(SwigDirectorserialize);
    if (SwigDerivedClassHasMethod("unserialize", swigMethodTypes5))
      swigDelegate5 = new SwigDelegateFeatureBase_5(SwigDirectorunserialize);
    libtischPINVOKE.FeatureBase_director_connect(swigCPtr, swigDelegate0, swigDelegate1, swigDelegate2, swigDelegate3, swigDelegate4, swigDelegate5);
  }

  private bool SwigDerivedClassHasMethod(string methodName, Type[] methodTypes) {
    System.Reflection.MethodInfo methodInfo = this.GetType().GetMethod(methodName, System.Reflection.BindingFlags.Public | System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance, null, methodTypes, null);
    bool hasDerivedMethod = methodInfo.DeclaringType.IsSubclassOf(typeof(FeatureBase));
    return hasDerivedMethod;
  }

  private string SwigDirectorname() {
    return name();
  }

  private IntPtr SwigDirectorclone() {
    return FeatureBase.getCPtr(clone()).Handle;
  }

  private void SwigDirectorload(IntPtr state) {
    load(new InputState(state, false));
  }

  private int SwigDirectornext() {
    return next();
  }

  private void SwigDirectorserialize(IntPtr s) {
    serialize(new SWIGTYPE_p_std__ostream(s, false));
  }

  private void SwigDirectorunserialize(IntPtr s) {
    unserialize(new SWIGTYPE_p_std__istream(s, false));
  }

  public delegate string SwigDelegateFeatureBase_0();
  public delegate IntPtr SwigDelegateFeatureBase_1();
  public delegate void SwigDelegateFeatureBase_2(IntPtr state);
  public delegate int SwigDelegateFeatureBase_3();
  public delegate void SwigDelegateFeatureBase_4(IntPtr s);
  public delegate void SwigDelegateFeatureBase_5(IntPtr s);

  private SwigDelegateFeatureBase_0 swigDelegate0;
  private SwigDelegateFeatureBase_1 swigDelegate1;
  private SwigDelegateFeatureBase_2 swigDelegate2;
  private SwigDelegateFeatureBase_3 swigDelegate3;
  private SwigDelegateFeatureBase_4 swigDelegate4;
  private SwigDelegateFeatureBase_5 swigDelegate5;

  private static Type[] swigMethodTypes0 = new Type[] {  };
  private static Type[] swigMethodTypes1 = new Type[] {  };
  private static Type[] swigMethodTypes2 = new Type[] { typeof(InputState) };
  private static Type[] swigMethodTypes3 = new Type[] {  };
  private static Type[] swigMethodTypes4 = new Type[] { typeof(SWIGTYPE_p_std__ostream) };
  private static Type[] swigMethodTypes5 = new Type[] { typeof(SWIGTYPE_p_std__istream) };
}
