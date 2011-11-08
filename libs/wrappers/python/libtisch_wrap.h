/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 1.3.40
 * 
 * This file is not intended to be easily readable and contains a number of 
 * coding conventions designed to improve portability and efficiency. Do not make
 * changes to this file unless you know what you are doing--modify the SWIG 
 * interface file instead. 
 * ----------------------------------------------------------------------------- */

#ifndef SWIG_libtisch_WRAP_H_
#define SWIG_libtisch_WRAP_H_

#include <map>
#include <string>


class SwigDirector_FeatureBase : public FeatureBase, public Swig::Director {

public:
    SwigDirector_FeatureBase(PyObject *self, unsigned int _tf = 0);
    virtual ~SwigDirector_FeatureBase();
    virtual char const *name() const;
    virtual FeatureBase *clone() const;
    virtual void load(InputState &state);
    virtual int next();
    virtual void serialize(std::ostream &s);
    virtual void unserialize(std::istream &s);
    using FeatureBase::typeflags;


/* Internal Director utilities */
public:
    bool swig_get_inner(const char* name) const {
      std::map<std::string, bool>::const_iterator iv = inner.find(name);
      return (iv != inner.end() ? iv->second : false);
    }

    void swig_set_inner(const char* name, bool val) const
    { inner[name] = val;}

private:
    mutable std::map<std::string, bool> inner;


#if defined(SWIG_PYTHON_DIRECTOR_VTABLE)
/* VTable implementation */
    PyObject *swig_get_method(size_t method_index, const char *method_name) const {
      PyObject *method = vtable[method_index];
      if (!method) {
        swig::SwigVar_PyObject name = SWIG_Python_str_FromChar(method_name);
        method = PyObject_GetAttr(swig_get_self(), name);
        if (method == NULL) {
          std::string msg = "Method in class FeatureBase doesn't exist, undefined ";
          msg += method_name;
          Swig::DirectorMethodException::raise(msg.c_str());
        }
        vtable[method_index] = method;
      };
      return method;
    }
private:
    mutable swig::SwigVar_PyObject vtable[6];
#endif

};


class SwigDirector_Widget : public Widget, public Swig::Director {

public:
    SwigDirector_Widget(PyObject *self, int _w, int _h, int _x = 0, int _y = 0, double _angle = 0.0, RGBATexture *_tex = 0, unsigned int _regflags = (((unsigned int) 1 << INPUT_TYPE_COUNT) -1));
    virtual ~SwigDirector_Widget();
    virtual void outline();
    virtual void update(Widget *target = 0);
    virtual void doUpdate(Widget *target = 0);
    virtual void raise(Widget *widget = 0);
    virtual void lower(Widget *widget = 0);
    virtual void draw();
    virtual void action(Gesture *gesture);
    virtual void enter(double z = 0.0);
    virtual void paint(bool update_stencil = false);
    using Widget::unregister;
    using Widget::mytex;
    using Widget::mycolor;
    using Widget::parent;
    using Widget::m_model;


/* Internal Director utilities */
public:
    bool swig_get_inner(const char* name) const {
      std::map<std::string, bool>::const_iterator iv = inner.find(name);
      return (iv != inner.end() ? iv->second : false);
    }

    void swig_set_inner(const char* name, bool val) const
    { inner[name] = val;}

private:
    mutable std::map<std::string, bool> inner;


#if defined(SWIG_PYTHON_DIRECTOR_VTABLE)
/* VTable implementation */
    PyObject *swig_get_method(size_t method_index, const char *method_name) const {
      PyObject *method = vtable[method_index];
      if (!method) {
        swig::SwigVar_PyObject name = SWIG_Python_str_FromChar(method_name);
        method = PyObject_GetAttr(swig_get_self(), name);
        if (method == NULL) {
          std::string msg = "Method in class Widget doesn't exist, undefined ";
          msg += method_name;
          Swig::DirectorMethodException::raise(msg.c_str());
        }
        vtable[method_index] = method;
      };
      return method;
    }
private:
    mutable swig::SwigVar_PyObject vtable[15];
#endif

};


class SwigDirector_Label : public Label, public Swig::Director {

public:
    SwigDirector_Label(PyObject *self, char const *text, int _w, int _h, int _x = 0, int _y = 0, double angle = 0.0, int center = 0, int snip = 0, RGBATexture *_tex = 0);
    virtual ~SwigDirector_Label();
    virtual void outline();
    virtual void update(Widget *target = 0);
    virtual void doUpdate(Widget *target = 0);
    virtual void raise(Widget *widget = 0);
    virtual void lower(Widget *widget = 0);
    virtual void draw();
    virtual void action(Gesture *gesture);
    virtual void enter(double z = 0.0);
    virtual void paint(bool update_stencil = false);
    using Label::text;
    using Label::center;
    using Label::snip;


/* Internal Director utilities */
public:
    bool swig_get_inner(const char* name) const {
      std::map<std::string, bool>::const_iterator iv = inner.find(name);
      return (iv != inner.end() ? iv->second : false);
    }

    void swig_set_inner(const char* name, bool val) const
    { inner[name] = val;}

private:
    mutable std::map<std::string, bool> inner;


#if defined(SWIG_PYTHON_DIRECTOR_VTABLE)
/* VTable implementation */
    PyObject *swig_get_method(size_t method_index, const char *method_name) const {
      PyObject *method = vtable[method_index];
      if (!method) {
        swig::SwigVar_PyObject name = SWIG_Python_str_FromChar(method_name);
        method = PyObject_GetAttr(swig_get_self(), name);
        if (method == NULL) {
          std::string msg = "Method in class Label doesn't exist, undefined ";
          msg += method_name;
          Swig::DirectorMethodException::raise(msg.c_str());
        }
        vtable[method_index] = method;
      };
      return method;
    }
private:
    mutable swig::SwigVar_PyObject vtable[15];
#endif

};


class SwigDirector_Button : public Button, public Swig::Director {

public:
    SwigDirector_Button(PyObject *self, int _w, int _h, int _x = 0, int _y = 0, double angle = 0.0, RGBATexture *_tex = 0);
    virtual ~SwigDirector_Button();
    virtual void outline();
    virtual void update(Widget *target = 0);
    virtual void doUpdate(Widget *target = 0);
    virtual void raise(Widget *widget = 0);
    virtual void lower(Widget *widget = 0);
    virtual void draw();
    virtual void action(Gesture *gesture);
    virtual void enter(double z = 0.0);
    virtual void paint(bool update_stencil = false);
    virtual void tap(Vector pos, int id);
    virtual void release();
    using Button::active;


/* Internal Director utilities */
public:
    bool swig_get_inner(const char* name) const {
      std::map<std::string, bool>::const_iterator iv = inner.find(name);
      return (iv != inner.end() ? iv->second : false);
    }

    void swig_set_inner(const char* name, bool val) const
    { inner[name] = val;}

private:
    mutable std::map<std::string, bool> inner;


#if defined(SWIG_PYTHON_DIRECTOR_VTABLE)
/* VTable implementation */
    PyObject *swig_get_method(size_t method_index, const char *method_name) const {
      PyObject *method = vtable[method_index];
      if (!method) {
        swig::SwigVar_PyObject name = SWIG_Python_str_FromChar(method_name);
        method = PyObject_GetAttr(swig_get_self(), name);
        if (method == NULL) {
          std::string msg = "Method in class Button doesn't exist, undefined ";
          msg += method_name;
          Swig::DirectorMethodException::raise(msg.c_str());
        }
        vtable[method_index] = method;
      };
      return method;
    }
private:
    mutable swig::SwigVar_PyObject vtable[17];
#endif

};


class SwigDirector_Tile : public Tile, public Swig::Director {

public:
    SwigDirector_Tile(PyObject *self, int _w, int _h, int _x = 0, int _y = 0, double angle = 0.0, RGBATexture *_tex = 0, int _mode = 0xFF);
    virtual ~SwigDirector_Tile();
    virtual void outline();
    virtual void update(Widget *target = 0);
    virtual void doUpdate(Widget *target = 0);
    virtual void raise(Widget *widget = 0);
    virtual void lower(Widget *widget = 0);
    virtual void draw();
    virtual void action(Gesture *gesture);
    virtual void enter(double z = 0.0);
    virtual void paint(bool update_stencil = false);
    virtual void tap(Vector pos, int id);
    virtual void release();
    virtual void apply(Vector delta);
    using Tile::mode;
    using Tile::slide;
    using Tile::vel;


/* Internal Director utilities */
public:
    bool swig_get_inner(const char* name) const {
      std::map<std::string, bool>::const_iterator iv = inner.find(name);
      return (iv != inner.end() ? iv->second : false);
    }

    void swig_set_inner(const char* name, bool val) const
    { inner[name] = val;}

private:
    mutable std::map<std::string, bool> inner;


#if defined(SWIG_PYTHON_DIRECTOR_VTABLE)
/* VTable implementation */
    PyObject *swig_get_method(size_t method_index, const char *method_name) const {
      PyObject *method = vtable[method_index];
      if (!method) {
        swig::SwigVar_PyObject name = SWIG_Python_str_FromChar(method_name);
        method = PyObject_GetAttr(swig_get_self(), name);
        if (method == NULL) {
          std::string msg = "Method in class Tile doesn't exist, undefined ";
          msg += method_name;
          Swig::DirectorMethodException::raise(msg.c_str());
        }
        vtable[method_index] = method;
      };
      return method;
    }
private:
    mutable swig::SwigVar_PyObject vtable[18];
#endif

};


class SwigDirector_Container : public Container, public Swig::Director {

public:
    SwigDirector_Container(PyObject *self, int w, int h, int x, int y, double angle = 0.0, RGBATexture *tex = 0, int mode = 32);
    virtual ~SwigDirector_Container();
    virtual void outline();
    virtual void update(Widget *target = 0);
    virtual void doUpdate(Widget *target = 0);
    virtual void raise(Widget *widget = 0);
    virtual void lower(Widget *widget = 0);
    virtual void draw();
    virtual void action(Gesture *gesture);
    virtual void enter(double z = 0.0);
    virtual void paint(bool update_stencil = false);
    virtual void tap(Vector vec, int id);
    virtual void release();
    virtual void apply(Vector delta);
    using Container::totalHeight;
    using Container::widgets;
    using Container::locked;


/* Internal Director utilities */
public:
    bool swig_get_inner(const char* name) const {
      std::map<std::string, bool>::const_iterator iv = inner.find(name);
      return (iv != inner.end() ? iv->second : false);
    }

    void swig_set_inner(const char* name, bool val) const
    { inner[name] = val;}

private:
    mutable std::map<std::string, bool> inner;


#if defined(SWIG_PYTHON_DIRECTOR_VTABLE)
/* VTable implementation */
    PyObject *swig_get_method(size_t method_index, const char *method_name) const {
      PyObject *method = vtable[method_index];
      if (!method) {
        swig::SwigVar_PyObject name = SWIG_Python_str_FromChar(method_name);
        method = PyObject_GetAttr(swig_get_self(), name);
        if (method == NULL) {
          std::string msg = "Method in class Container doesn't exist, undefined ";
          msg += method_name;
          Swig::DirectorMethodException::raise(msg.c_str());
        }
        vtable[method_index] = method;
      };
      return method;
    }
private:
    mutable swig::SwigVar_PyObject vtable[18];
#endif

};


class SwigDirector_Slider : public Slider, public Swig::Director {

public:
    SwigDirector_Slider(PyObject *self, int _w, int _h, int _x = 0, int _y = 0, double angle = 0.0, RGBATexture *_tex = 0);
    virtual ~SwigDirector_Slider();
    virtual void outline();
    virtual void update(Widget *target = 0);
    virtual void doUpdate(Widget *target = 0);
    virtual void raise(Widget *widget = 0);
    virtual void lower(Widget *widget = 0);
    virtual void draw();
    virtual void action(Gesture *gesture);
    virtual void enter(double z = 0.0);
    virtual void paint(bool update_stencil = false);
    using Slider::pos;


/* Internal Director utilities */
public:
    bool swig_get_inner(const char* name) const {
      std::map<std::string, bool>::const_iterator iv = inner.find(name);
      return (iv != inner.end() ? iv->second : false);
    }

    void swig_set_inner(const char* name, bool val) const
    { inner[name] = val;}

private:
    mutable std::map<std::string, bool> inner;


#if defined(SWIG_PYTHON_DIRECTOR_VTABLE)
/* VTable implementation */
    PyObject *swig_get_method(size_t method_index, const char *method_name) const {
      PyObject *method = vtable[method_index];
      if (!method) {
        swig::SwigVar_PyObject name = SWIG_Python_str_FromChar(method_name);
        method = PyObject_GetAttr(swig_get_self(), name);
        if (method == NULL) {
          std::string msg = "Method in class Slider doesn't exist, undefined ";
          msg += method_name;
          Swig::DirectorMethodException::raise(msg.c_str());
        }
        vtable[method_index] = method;
      };
      return method;
    }
private:
    mutable swig::SwigVar_PyObject vtable[15];
#endif

};


class SwigDirector_Dial : public Dial, public Swig::Director {

public:
    SwigDirector_Dial(PyObject *self, int _r, int _x = 0, int _y = 0, double _angle = 0.0, RGBATexture *_tex = 0);
    virtual ~SwigDirector_Dial();
    virtual void outline();
    virtual void update(Widget *target = 0);
    virtual void doUpdate(Widget *target = 0);
    virtual void raise(Widget *widget = 0);
    virtual void lower(Widget *widget = 0);
    virtual void draw();
    virtual void action(Gesture *gesture);
    virtual void enter(double z = 0.0);
    virtual void paint(bool update_stencil = false);
    using Dial::k_angle;
    using Dial::k_lower;
    using Dial::k_upper;
    using Dial::oldpos;


/* Internal Director utilities */
public:
    bool swig_get_inner(const char* name) const {
      std::map<std::string, bool>::const_iterator iv = inner.find(name);
      return (iv != inner.end() ? iv->second : false);
    }

    void swig_set_inner(const char* name, bool val) const
    { inner[name] = val;}

private:
    mutable std::map<std::string, bool> inner;


#if defined(SWIG_PYTHON_DIRECTOR_VTABLE)
/* VTable implementation */
    PyObject *swig_get_method(size_t method_index, const char *method_name) const {
      PyObject *method = vtable[method_index];
      if (!method) {
        swig::SwigVar_PyObject name = SWIG_Python_str_FromChar(method_name);
        method = PyObject_GetAttr(swig_get_self(), name);
        if (method == NULL) {
          std::string msg = "Method in class Dial doesn't exist, undefined ";
          msg += method_name;
          Swig::DirectorMethodException::raise(msg.c_str());
        }
        vtable[method_index] = method;
      };
      return method;
    }
private:
    mutable swig::SwigVar_PyObject vtable[15];
#endif

};


class SwigDirector_MasterContainer : public MasterContainer, public Swig::Director {

public:
    SwigDirector_MasterContainer(PyObject *self, int w, int h, int defaults = 1);
    virtual ~SwigDirector_MasterContainer();
    virtual void outline();
    virtual void update(Widget *target = 0);
    virtual void doUpdate(Widget *target = 0);
    virtual void raise(Widget *widget = 0);
    virtual void lower(Widget *widget = 0);
    virtual void draw();
    virtual void action(Gesture *gesture);
    virtual void enter(double z = 0.0);
    virtual void paint(bool update_stencil = false);
    virtual void tap(Vector vec, int id);
    virtual void release();
    virtual void apply(Vector delta);
    using MasterContainer::matcher;
    using MasterContainer::input;
    using MasterContainer::inthread;


/* Internal Director utilities */
public:
    bool swig_get_inner(const char* name) const {
      std::map<std::string, bool>::const_iterator iv = inner.find(name);
      return (iv != inner.end() ? iv->second : false);
    }

    void swig_set_inner(const char* name, bool val) const
    { inner[name] = val;}

private:
    mutable std::map<std::string, bool> inner;


#if defined(SWIG_PYTHON_DIRECTOR_VTABLE)
/* VTable implementation */
    PyObject *swig_get_method(size_t method_index, const char *method_name) const {
      PyObject *method = vtable[method_index];
      if (!method) {
        swig::SwigVar_PyObject name = SWIG_Python_str_FromChar(method_name);
        method = PyObject_GetAttr(swig_get_self(), name);
        if (method == NULL) {
          std::string msg = "Method in class MasterContainer doesn't exist, undefined ";
          msg += method_name;
          Swig::DirectorMethodException::raise(msg.c_str());
        }
        vtable[method_index] = method;
      };
      return method;
    }
private:
    mutable swig::SwigVar_PyObject vtable[18];
#endif

};


class SwigDirector_Window : public Window, public Swig::Director {

public:
    SwigDirector_Window(PyObject *self, int w, int h, std::string title, int use_mouse = 0);
    virtual ~SwigDirector_Window();
    virtual void idle();
    virtual void display();
    virtual void reshape(int w, int h);
    virtual void keyboard(int key, int x, int y);
    virtual void mouse(int num, int button, int state, int x, int y);
    virtual void passive(int num, int x, int y);
    virtual void motion(int num, int x, int y);
    virtual void entry(int num, int state);
    virtual void outline();
    virtual void update(Widget *target = 0);
    virtual void doUpdate(Widget *target = 0);
    virtual void raise(Widget *widget = 0);
    virtual void lower(Widget *widget = 0);
    virtual void draw();
    virtual void action(Gesture *gesture);
    virtual void enter(double z = 0.0);
    virtual void paint(bool update_stencil = false);
    virtual void tap(Vector vec, int id);
    virtual void release();
    virtual void apply(Vector delta);


/* Internal Director utilities */
public:
    bool swig_get_inner(const char* name) const {
      std::map<std::string, bool>::const_iterator iv = inner.find(name);
      return (iv != inner.end() ? iv->second : false);
    }

    void swig_set_inner(const char* name, bool val) const
    { inner[name] = val;}

private:
    mutable std::map<std::string, bool> inner;


#if defined(SWIG_PYTHON_DIRECTOR_VTABLE)
/* VTable implementation */
    PyObject *swig_get_method(size_t method_index, const char *method_name) const {
      PyObject *method = vtable[method_index];
      if (!method) {
        swig::SwigVar_PyObject name = SWIG_Python_str_FromChar(method_name);
        method = PyObject_GetAttr(swig_get_self(), name);
        if (method == NULL) {
          std::string msg = "Method in class Window doesn't exist, undefined ";
          msg += method_name;
          Swig::DirectorMethodException::raise(msg.c_str());
        }
        vtable[method_index] = method;
      };
      return method;
    }
private:
    mutable swig::SwigVar_PyObject vtable[26];
#endif

};


#endif
