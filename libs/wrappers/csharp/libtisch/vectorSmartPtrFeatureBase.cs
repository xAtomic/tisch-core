/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 1.3.40
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */


using System;
using System.Runtime.InteropServices;

public class vectorSmartPtrFeatureBase : IDisposable, System.Collections.IEnumerable
#if !SWIG_DOTNET_1
    , System.Collections.Generic.IEnumerable<smartPtrFeatureBase>
#endif
 {
  private HandleRef swigCPtr;
  protected bool swigCMemOwn;

  internal vectorSmartPtrFeatureBase(IntPtr cPtr, bool cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = new HandleRef(this, cPtr);
  }

  internal static HandleRef getCPtr(vectorSmartPtrFeatureBase obj) {
    return (obj == null) ? new HandleRef(null, IntPtr.Zero) : obj.swigCPtr;
  }

  ~vectorSmartPtrFeatureBase() {
    Dispose();
  }

  public virtual void Dispose() {
    lock(this) {
      if (swigCPtr.Handle != IntPtr.Zero) {
        if (swigCMemOwn) {
          swigCMemOwn = false;
          libtischPINVOKE.delete_vectorSmartPtrFeatureBase(swigCPtr);
        }
        swigCPtr = new HandleRef(null, IntPtr.Zero);
      }
      GC.SuppressFinalize(this);
    }
  }

  public vectorSmartPtrFeatureBase(System.Collections.ICollection c) : this() {
    if (c == null)
      throw new ArgumentNullException("c");
    foreach (smartPtrFeatureBase element in c) {
      this.Add(element);
    }
  }

  public bool IsFixedSize {
    get {
      return false;
    }
  }

  public bool IsReadOnly {
    get {
      return false;
    }
  }

  public smartPtrFeatureBase this[int index]  {
    get {
      return getitem(index);
    }
    set {
      setitem(index, value);
    }
  }

  public int Capacity {
    get {
      return (int)capacity();
    }
    set {
      if (value < size())
        throw new ArgumentOutOfRangeException("Capacity");
      reserve((uint)value);
    }
  }

  public int Count {
    get {
      return (int)size();
    }
  }

  public bool IsSynchronized {
    get {
      return false;
    }
  }

#if SWIG_DOTNET_1
  public void CopyTo(System.Array array)
#else
  public void CopyTo(smartPtrFeatureBase[] array)
#endif
  {
    CopyTo(0, array, 0, this.Count);
  }

#if SWIG_DOTNET_1
  public void CopyTo(System.Array array, int arrayIndex)
#else
  public void CopyTo(smartPtrFeatureBase[] array, int arrayIndex)
#endif
  {
    CopyTo(0, array, arrayIndex, this.Count);
  }

#if SWIG_DOTNET_1
  public void CopyTo(int index, System.Array array, int arrayIndex, int count)
#else
  public void CopyTo(int index, smartPtrFeatureBase[] array, int arrayIndex, int count)
#endif
  {
    if (array == null)
      throw new ArgumentNullException("array");
    if (index < 0)
      throw new ArgumentOutOfRangeException("index", "Value is less than zero");
    if (arrayIndex < 0)
      throw new ArgumentOutOfRangeException("arrayIndex", "Value is less than zero");
    if (count < 0)
      throw new ArgumentOutOfRangeException("count", "Value is less than zero");
    if (array.Rank > 1)
      throw new ArgumentException("Multi dimensional array.", "array");
    if (index+count > this.Count || arrayIndex+count > array.Length)
      throw new ArgumentException("Number of elements to copy is too large.");
    for (int i=0; i<count; i++)
      array.SetValue(getitemcopy(index+i), arrayIndex+i);
  }

#if !SWIG_DOTNET_1
  System.Collections.Generic.IEnumerator<smartPtrFeatureBase> System.Collections.Generic.IEnumerable<smartPtrFeatureBase>.GetEnumerator() {
    return new vectorSmartPtrFeatureBaseEnumerator(this);
  }
#endif

  System.Collections.IEnumerator System.Collections.IEnumerable.GetEnumerator() {
    return new vectorSmartPtrFeatureBaseEnumerator(this);
  }

  public vectorSmartPtrFeatureBaseEnumerator GetEnumerator() {
    return new vectorSmartPtrFeatureBaseEnumerator(this);
  }

  // Type-safe enumerator
  /// Note that the IEnumerator documentation requires an InvalidOperationException to be thrown
  /// whenever the collection is modified. This has been done for changes in the size of the
  /// collection but not when one of the elements of the collection is modified as it is a bit
  /// tricky to detect unmanaged code that modifies the collection under our feet.
  public sealed class vectorSmartPtrFeatureBaseEnumerator : System.Collections.IEnumerator
#if !SWIG_DOTNET_1
    , System.Collections.Generic.IEnumerator<smartPtrFeatureBase>
#endif
  {
    private vectorSmartPtrFeatureBase collectionRef;
    private int currentIndex;
    private object currentObject;
    private int currentSize;

    public vectorSmartPtrFeatureBaseEnumerator(vectorSmartPtrFeatureBase collection) {
      collectionRef = collection;
      currentIndex = -1;
      currentObject = null;
      currentSize = collectionRef.Count;
    }

    // Type-safe iterator Current
    public smartPtrFeatureBase Current {
      get {
        if (currentIndex == -1)
          throw new InvalidOperationException("Enumeration not started.");
        if (currentIndex > currentSize - 1)
          throw new InvalidOperationException("Enumeration finished.");
        if (currentObject == null)
          throw new InvalidOperationException("Collection modified.");
        return (smartPtrFeatureBase)currentObject;
      }
    }

    // Type-unsafe IEnumerator.Current
    object System.Collections.IEnumerator.Current {
      get {
        return Current;
      }
    }

    public bool MoveNext() {
      int size = collectionRef.Count;
      bool moveOkay = (currentIndex+1 < size) && (size == currentSize);
      if (moveOkay) {
        currentIndex++;
        currentObject = collectionRef[currentIndex];
      } else {
        currentObject = null;
      }
      return moveOkay;
    }

    public void Reset() {
      currentIndex = -1;
      currentObject = null;
      if (collectionRef.Count != currentSize) {
        throw new InvalidOperationException("Collection modified.");
      }
    }

#if !SWIG_DOTNET_1
    public void Dispose() {
        currentIndex = -1;
        currentObject = null;
    }
#endif
  }

  public void Clear() {
    libtischPINVOKE.vectorSmartPtrFeatureBase_Clear(swigCPtr);
    if (libtischPINVOKE.SWIGPendingException.Pending) throw libtischPINVOKE.SWIGPendingException.Retrieve();
  }

  public void Add(smartPtrFeatureBase x) {
    libtischPINVOKE.vectorSmartPtrFeatureBase_Add(swigCPtr, smartPtrFeatureBase.getCPtr(x));
    if (libtischPINVOKE.SWIGPendingException.Pending) throw libtischPINVOKE.SWIGPendingException.Retrieve();
  }

  private uint size() {
    uint ret = libtischPINVOKE.vectorSmartPtrFeatureBase_size(swigCPtr);
    if (libtischPINVOKE.SWIGPendingException.Pending) throw libtischPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  private uint capacity() {
    uint ret = libtischPINVOKE.vectorSmartPtrFeatureBase_capacity(swigCPtr);
    if (libtischPINVOKE.SWIGPendingException.Pending) throw libtischPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  private void reserve(uint n) {
    libtischPINVOKE.vectorSmartPtrFeatureBase_reserve(swigCPtr, n);
    if (libtischPINVOKE.SWIGPendingException.Pending) throw libtischPINVOKE.SWIGPendingException.Retrieve();
  }

  public vectorSmartPtrFeatureBase() : this(libtischPINVOKE.new_vectorSmartPtrFeatureBase__SWIG_0(), true) {
    if (libtischPINVOKE.SWIGPendingException.Pending) throw libtischPINVOKE.SWIGPendingException.Retrieve();
  }

  public vectorSmartPtrFeatureBase(vectorSmartPtrFeatureBase other) : this(libtischPINVOKE.new_vectorSmartPtrFeatureBase__SWIG_1(vectorSmartPtrFeatureBase.getCPtr(other)), true) {
    if (libtischPINVOKE.SWIGPendingException.Pending) throw libtischPINVOKE.SWIGPendingException.Retrieve();
  }

  public vectorSmartPtrFeatureBase(int capacity) : this(libtischPINVOKE.new_vectorSmartPtrFeatureBase__SWIG_2(capacity), true) {
    if (libtischPINVOKE.SWIGPendingException.Pending) throw libtischPINVOKE.SWIGPendingException.Retrieve();
  }

  private smartPtrFeatureBase getitemcopy(int index) {
    smartPtrFeatureBase ret = new smartPtrFeatureBase(libtischPINVOKE.vectorSmartPtrFeatureBase_getitemcopy(swigCPtr, index), true);
    if (libtischPINVOKE.SWIGPendingException.Pending) throw libtischPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  private smartPtrFeatureBase getitem(int index) {
    smartPtrFeatureBase ret = new smartPtrFeatureBase(libtischPINVOKE.vectorSmartPtrFeatureBase_getitem(swigCPtr, index), false);
    if (libtischPINVOKE.SWIGPendingException.Pending) throw libtischPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  private void setitem(int index, smartPtrFeatureBase val) {
    libtischPINVOKE.vectorSmartPtrFeatureBase_setitem(swigCPtr, index, smartPtrFeatureBase.getCPtr(val));
    if (libtischPINVOKE.SWIGPendingException.Pending) throw libtischPINVOKE.SWIGPendingException.Retrieve();
  }

  public void AddRange(vectorSmartPtrFeatureBase values) {
    libtischPINVOKE.vectorSmartPtrFeatureBase_AddRange(swigCPtr, vectorSmartPtrFeatureBase.getCPtr(values));
    if (libtischPINVOKE.SWIGPendingException.Pending) throw libtischPINVOKE.SWIGPendingException.Retrieve();
  }

  public vectorSmartPtrFeatureBase GetRange(int index, int count) {
    IntPtr cPtr = libtischPINVOKE.vectorSmartPtrFeatureBase_GetRange(swigCPtr, index, count);
    vectorSmartPtrFeatureBase ret = (cPtr == IntPtr.Zero) ? null : new vectorSmartPtrFeatureBase(cPtr, true);
    if (libtischPINVOKE.SWIGPendingException.Pending) throw libtischPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public void Insert(int index, smartPtrFeatureBase x) {
    libtischPINVOKE.vectorSmartPtrFeatureBase_Insert(swigCPtr, index, smartPtrFeatureBase.getCPtr(x));
    if (libtischPINVOKE.SWIGPendingException.Pending) throw libtischPINVOKE.SWIGPendingException.Retrieve();
  }

  public void InsertRange(int index, vectorSmartPtrFeatureBase values) {
    libtischPINVOKE.vectorSmartPtrFeatureBase_InsertRange(swigCPtr, index, vectorSmartPtrFeatureBase.getCPtr(values));
    if (libtischPINVOKE.SWIGPendingException.Pending) throw libtischPINVOKE.SWIGPendingException.Retrieve();
  }

  public void RemoveAt(int index) {
    libtischPINVOKE.vectorSmartPtrFeatureBase_RemoveAt(swigCPtr, index);
    if (libtischPINVOKE.SWIGPendingException.Pending) throw libtischPINVOKE.SWIGPendingException.Retrieve();
  }

  public void RemoveRange(int index, int count) {
    libtischPINVOKE.vectorSmartPtrFeatureBase_RemoveRange(swigCPtr, index, count);
    if (libtischPINVOKE.SWIGPendingException.Pending) throw libtischPINVOKE.SWIGPendingException.Retrieve();
  }

  public static vectorSmartPtrFeatureBase Repeat(smartPtrFeatureBase value, int count) {
    IntPtr cPtr = libtischPINVOKE.vectorSmartPtrFeatureBase_Repeat(smartPtrFeatureBase.getCPtr(value), count);
    vectorSmartPtrFeatureBase ret = (cPtr == IntPtr.Zero) ? null : new vectorSmartPtrFeatureBase(cPtr, true);
    if (libtischPINVOKE.SWIGPendingException.Pending) throw libtischPINVOKE.SWIGPendingException.Retrieve();
    return ret;
  }

  public void Reverse() {
    libtischPINVOKE.vectorSmartPtrFeatureBase_Reverse__SWIG_0(swigCPtr);
    if (libtischPINVOKE.SWIGPendingException.Pending) throw libtischPINVOKE.SWIGPendingException.Retrieve();
  }

  public void Reverse(int index, int count) {
    libtischPINVOKE.vectorSmartPtrFeatureBase_Reverse__SWIG_1(swigCPtr, index, count);
    if (libtischPINVOKE.SWIGPendingException.Pending) throw libtischPINVOKE.SWIGPendingException.Retrieve();
  }

  public void SetRange(int index, vectorSmartPtrFeatureBase values) {
    libtischPINVOKE.vectorSmartPtrFeatureBase_SetRange(swigCPtr, index, vectorSmartPtrFeatureBase.getCPtr(values));
    if (libtischPINVOKE.SWIGPendingException.Pending) throw libtischPINVOKE.SWIGPendingException.Retrieve();
  }

}
