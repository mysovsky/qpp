#ifndef _QPP_SHAPE_H
#define _QPP_SHAPE_H

//#include <mathf/constants.hpp>
//#include <data/qppdata.hpp>
#include <geom/geom.hpp>
#include <geom/lace3d.hpp>
#include <algorithm>

#ifdef PY_EXPORT
#include <pybind11/pybind11.h>
namespace py = pybind11;
#endif

#define v2d vector2d<VALTYPE>
#define v3d vector3d<VALTYPE>

namespace qpp{

  // forward declarations

  template <class VALTYPE> class shape_union;
  template <class VALTYPE> class shape_intersect;
  template <class VALTYPE> class shape_subtract;
  template <class VALTYPE> class shape_invert;
  template <class VALTYPE> class shape_xor;

  // ------------------ 3D primitives prototype ------------------------

  template <class VALTYPE>
  class shape{ 

  public:

    STRING name;

    shape( const STRING &__name = "")
    {  name = __name; }

    /*
    shape(const STRING & __name = "", qpp_object * __owner = NULL):
      qpp_object(__name,__owner)
      {}*/

    virtual bool within(const v3d & r) const =0;
    // Answers the question whether the r point is situated within the shape

    virtual void scale(VALTYPE s) =0;
    virtual void move(const v3d & v) =0;
    virtual void rotate(const matrix3d<VALTYPE> & Rot) =0;

    // -------------------------------------
    
    virtual VALTYPE volume() const =0;
    virtual v3d rmin() const =0;
    virtual v3d rmax() const =0;
    // Minimal & maximal cartesian coordinates of the shape

    virtual v3d fmin(const periodic_cell<VALTYPE> &v) const =0;
    virtual v3d fmax(const periodic_cell<VALTYPE> &v) const =0;
    // Minimal and maximal fractional coordinates of the shape for given cell

    virtual void write(std::basic_ostream<CHAR,TRAITS> &os, int offset=0) const =0;

    shape<VALTYPE> & operator|(shape<VALTYPE> & sh)
    {  return *new shape_union<VALTYPE>(*this,sh);}

    shape<VALTYPE> & operator&(shape<VALTYPE> & sh)
    {  return *new shape_intersect<VALTYPE>(*this,sh);}

    shape<VALTYPE> & operator-(shape<VALTYPE> & sh)
    {  return *new shape_subtract<VALTYPE>(*this,sh);}

    shape<VALTYPE> & operator^(shape<VALTYPE> & sh)
    {  return *new shape_xor<VALTYPE>(*this,sh);}

    shape<VALTYPE> & operator~()
    {  return *new shape_invert<VALTYPE>(*this);}

#ifdef PY_EXPORT

    shape<VALTYPE> * py_or(shape<VALTYPE> & sh)
    { return new shape_union<VALTYPE>(*this,sh); }

    shape<VALTYPE> * py_and(shape<VALTYPE> & sh)
    { return new shape_intersect<VALTYPE>(*this,sh); }

    shape<VALTYPE> * py_sub(shape<VALTYPE> & sh)
    { return new shape_subtract<VALTYPE>(*this,sh); }

    shape<VALTYPE> * py_xor(shape<VALTYPE> & sh)
    { return new shape_xor<VALTYPE>(*this,sh); }

    shape<VALTYPE> * py_inv()
    { return new shape_invert<VALTYPE>(*this); }

    static void py_export(const char * pyname)
    {
      bp::class_<shape<VALTYPE>, boost::noncopyable>(pyname, bp::no_init)
	.def("__or__",     & shape<VALTYPE>::py_or,  bp::return_value_policy<bp::manage_new_object>())
	.def("__and__",    & shape<VALTYPE>::py_and, bp::return_value_policy<bp::manage_new_object>())
	.def("__sub__",    & shape<VALTYPE>::py_sub, bp::return_value_policy<bp::manage_new_object>())
	.def("__xor__",    & shape<VALTYPE>::py_xor, bp::return_value_policy<bp::manage_new_object>())
	.def("__invert__", & shape<VALTYPE>::py_inv, bp::return_value_policy<bp::manage_new_object>())
	;
    }

#endif

 };
    

  // ------------------ 3D primitives ------------------------

  template <class VALTYPE>
  class shape_box : public shape<VALTYPE>{

    v3d crn, a[3];

    void fill_corners(v3d *corners) const
    {
      corners[0] = crn;
      corners[1] = crn+a[0];
      corners[2] = crn+a[1];
      corners[3] = crn+a[0]+a[1];
      corners[4] = crn+a[2];
      corners[5] = crn+a[0]+a[2];
      corners[6] = crn+a[1]+a[2];
      corners[7] = crn+a[0]+a[1]+a[2];
    }

  public:

    using shape<VALTYPE>::name;

    shape_box(){}

    shape_box(const v3d & a1, const v3d & a2, const v3d & a3, const v3d & r0, 
		       const STRING & __name = "") :
      shape<VALTYPE>(__name)
    {
      crn = r0;
      a[0] = a1;
      a[1] = a2;
      a[2] = a3;
    }

    shape_box(const v3d & a1, const v3d & a2, const v3d & a3, const STRING & __name = "") :
      shape<VALTYPE>(__name)
    {
       crn = 0e0;
      a[0] = a1;
      a[1] = a2;
      a[2] = a3;
    }
    
    shape_box(VALTYPE a1, VALTYPE a2, VALTYPE a3, const STRING & __name = "") :
      shape<VALTYPE>(__name)
    {
      crn = 0e0;
      a[0] = v3d(a1,  0e0, 0e0);
      a[1] = v3d(0e0, a2,  0e0);
      a[2] = v3d(0e0, 0e0, a3);
    }

    shape_box(const shape_box<VALTYPE> & s) :
      shape<VALTYPE>(s.name)
    {
      crn = s.crn;
      a[0] = s.a[0];
      a[1] = s.a[1];
      a[2] = s.a[2];
    }

    virtual void write(std::basic_ostream<CHAR,TRAITS> &os, int offset=0) const
    {
      for (int i=0; i<offset; i++)
	os << " ";
      os << "box";
      if (name != "")
	os << " " << name;
      os << "( a" << a[0] << ", b" << a[1] << ", c" << a[2] << ", corner" << crn << ")";
    }

    // --------------------------------------------------

    virtual bool within(const v3d & r) const
    {
      if ( scal(r-crn,a[1]%a[2])*scal(r-crn-a[0],a[1]%a[2]) > 0e0)
	return false;
      if ( scal(r-crn,a[2]%a[0])*scal(r-crn-a[1],a[2]%a[0]) > 0e0)
	return false;
      if ( scal(r-crn,a[0]%a[1])*scal(r-crn-a[2],a[0]%a[1]) > 0e0)
	return false;
      return true;
    }
    
    virtual v3d rmin() const
    {
      v3d corners[8];
      fill_corners(corners);

      v3d res=crn;
      for (int i=0; i<8; i++)
	for (int j=0; j<3; j++)
	  if ( res(j) > corners[i](j))
	    res(j) = corners[i](j);

      return res;
    }

    virtual v3d rmax() const
    {
      v3d corners[8];
      fill_corners(corners);

      v3d res=crn;
      for (int i=0; i<8; i++)
	for (int j=0; j<3; j++)
	  if ( res(j) < corners[i](j))
	    res(j) = corners[i](j);

      return res;
    }

    virtual v3d fmin(const periodic_cell<VALTYPE> &v) const
    {
      v3d corners[8];
      fill_corners(corners);

      for (int i=0; i<8; i++)
	corners[i] = v.cart2frac(corners[i]);

      v3d res=corners[0];
      for (int i=1; i<8; i++)
	for (int j=0; j<3; j++)
	  if ( res(j) > corners[i](j))
	    res(j) = corners[i](j);
      
      return res;
    }

    virtual v3d fmax(const periodic_cell<VALTYPE> &v) const
    {
      v3d corners[8];
      fill_corners(corners);

      for (int i=0; i<8; i++)
	corners[i] = v.cart2frac(corners[i]);

      v3d res=corners[0];
      for (int i=1; i<8; i++)
	for (int j=0; j<3; j++)
	  if ( res(j) < corners[i](j))
	    res(j) = corners[i](j);
      
      return res;
    }

    // Minimal & maximal fractional coordinates of the shape for given translation vectors v

    virtual VALTYPE volume() const
    {
      return std::abs(det(a[0], a[1], a[2]));
    }

    virtual void scale(VALTYPE s)
    {
      crn *= s;
      a[0] *= s;
      a[1] *= s;
      a[2] *= s;
    }

    virtual void move(const v3d & v)
    {
      crn += v;
    }

    virtual void rotate(const matrix3d<VALTYPE> & Rot)
    {
      crn  = Rot*crn;
      a[0] = Rot*a[0];
      a[1] = Rot*a[1];
      a[2] = Rot*a[2];
    }

#ifdef PY_EXPORT

    static void py_export(const char * pyname)
    {
      bp::class_<shape_box<VALTYPE>, bp::bases<shape<VALTYPE> > >(pyname)
	.def(init<const v3d&, const v3d&, const v3d&, const v3d&, bp::optional<const STRING &> >())
	.def(init<const v3d&, const v3d&, const v3d&, bp::optional<const STRING &> >())
	.def(init<VALTYPE, VALTYPE, VALTYPE, bp::optional<const STRING &> >())        
	;
    }

#endif

  };

#ifdef PY_EXPORT
    
  template <class VALTYPE>
  shape<VALTYPE> * py_shape_box1(const v3d& a, const v3d& b, const v3d& c, const v3d& r)
  {
    return new shape_box<VALTYPE>(a,b,c,r);
  }

  template <class VALTYPE>
  shape<VALTYPE> * py_shape_box2(const v3d& a, const v3d& b, const v3d& c)
  {
    return new shape_box<VALTYPE>(a,b,c);
  }
  
  template <class VALTYPE>
  shape<VALTYPE> * py_shape_box3(VALTYPE a, VALTYPE b, VALTYPE c)
  {
    return new shape_box<VALTYPE>(a,b,c);
  }

#endif


  // ----------------------------------------------------------------

  template <class VALTYPE>
  class shape_sphere : public shape<VALTYPE>{

    VALTYPE R;
    v3d r0;

  public:

    using shape<VALTYPE>::name;

    shape_sphere(){}

    shape_sphere(VALTYPE _R, const STRING & __name = "") : shape<VALTYPE>(__name)
    {
      R = _R;
      r0 = 0e0;
    }

    shape_sphere(VALTYPE _R, const v3d & _r0, const STRING & __name = "") :
      shape<VALTYPE>(__name)
    {
      R = _R;
      r0 = _r0;
    }

    shape_sphere(const shape_sphere<VALTYPE> & s) : shape<VALTYPE>(s.name)
    {
      R = s.R;
      r0 = s.r0;
    }

    virtual void write(std::basic_ostream<CHAR,TRAITS> &os, int offset=0) const
    {
      for (int i=0; i<offset; i++)
	os << " ";
      os << "sphere";
      if (name != "")
	os << " " << name;
      os << "( R=" << R << ", center" << r0 << ")";
    }

    // --------------------------------------------------

    virtual bool within(const v3d & r) const
    {
      return norm(r-r0) <= R;
    }
    
    virtual v3d rmin() const
    {
      return r0 - v3d(R,R,R);
    }

    virtual v3d rmax() const
    {
      return r0 + v3d(R,R,R);
    }

    virtual v3d fmin(const periodic_cell<VALTYPE> &v) const
    {
      matrix3d<VALTYPE> A(v(0),v(1),v(2));
      matrix3d<VALTYPE> B = invert(A);

      v3d res = 0e0;
      for (int i=0; i<3; i++)
	for (int j=0; j<3; j++)
	  res(i) += B(i,j)*B(i,j);
      for (int i=0; i<3; i++)
	res(i) = std::sqrt(res(i));
      return B*r0 - R*res;
    }

    virtual v3d fmax(const periodic_cell<VALTYPE> &v) const
    {      
      matrix3d<VALTYPE> A(v(0),v(1),v(2));
      matrix3d<VALTYPE> B = invert(A);

      v3d res = 0e0;
      for (int i=0; i<3; i++)
	for (int j=0; j<3; j++)
	  res(i) += B(i,j)*B(i,j);
      for (int i=0; i<3; i++)
	res(i) = std::sqrt(res(i));
      return B*r0 + R*res;
    }

    // Minimal & maximal fractional coordinates of the shape for given translation vectors v

    virtual VALTYPE volume() const
    {
      return 4*qpp::pi*R*R*R/3;
    }

    virtual void scale(VALTYPE s)
    {
      R *= s;
    }

    virtual void move(const v3d & v)
    {
      r0 += v;
    }

    virtual void rotate(const matrix3d<VALTYPE> & Rot)
    {
      r0 = Rot*r0;
    }

#ifdef PY_EXPORT

    static void py_export(const char * pyname)
    {
      bp::class_<shape_sphere<VALTYPE>, bp::bases<shape<VALTYPE> > >(pyname)
	.def(bp::init<VALTYPE, bp::optional<const STRING &> >())
	.def(bp::init<VALTYPE, const v3d&, bp::optional<const STRING &> >())
	;
    }

#endif

  };

  // ----------------------------------------------------------------

  template <class VALTYPE>
  class shape_union : public shape<VALTYPE>{

    shape<VALTYPE> *sh1, *sh2;

  public:
    using shape<VALTYPE>::name;

    shape_union(shape<VALTYPE> & __sh1, shape<VALTYPE> &__sh2, const STRING & __name = ""):
      shape<VALTYPE>(__name)
    { sh1 = &__sh1; sh2 = &__sh2; }

    shape_union(const shape_union<VALTYPE> & s) :  shape<VALTYPE>(s.name)
    {
      sh1 = s.sh1;
      sh2 = s.sh2;
    }

    virtual bool within(const v3d & r) const
    { return sh1->within(r) || sh2->within(r); }
 
    virtual VALTYPE volume() const
    { return 0; }

    virtual void scale(VALTYPE s)
    { sh1 -> scale(s); sh2 -> scale(s); }

    virtual void move(const v3d & v)
    { sh1 -> move(v); sh2 -> move(v);}

    virtual void rotate(const matrix3d<VALTYPE> & Rot)
    { sh1->rotate(Rot); sh2->rotate(Rot); }

    virtual v3d rmin() const
    { 
      v3d r, r1 = sh1->rmin(), r2 = sh2->rmin();
      for (int i=0; i<3; i++)
	r(i) = std::min(r1(i),r2(i));
      return r;
    }

    virtual v3d rmax() const
    { 
      v3d r, r1 = sh1->rmax(), r2 = sh2->rmax();
      for (int i=0; i<3; i++)
	r(i) = std::max(r1(i),r2(i));
      return r;
    }

    virtual v3d fmin(const periodic_cell<VALTYPE> &v) const
    {
      v3d f, f1 = sh1->fmin(v), f2 = sh2->fmin(v);
      for (int i=0; i<3; i++)
	f(i) = std::min(f1(i),f2(i));
      return f;
    }

    virtual v3d fmax(const periodic_cell<VALTYPE> &v) const
    {
      v3d f, f1 = sh1->fmax(v), f2 = sh2->fmax(v);
      for (int i=0; i<3; i++)
	f(i) = std::max(f1(i),f2(i));
      return f;
    }

    virtual void write(std::basic_ostream<CHAR,TRAITS> &os, int offset=0) const
    {
      for (int i=0; i<offset; i++) os << " ";
      os << "union";
      if (name != "")
	os << " " << name;
      os << "( shape1 = ";
      sh1->write(os);
      os << ", shape2 = ";
      sh2->write(os);
      os << ")";
    }

  };

  // ----------------------------------------------------------------

  template <class VALTYPE>
  class shape_intersect : public shape<VALTYPE>{

    shape<VALTYPE> *sh1, *sh2;

  public:
    using shape<VALTYPE>::name;

    shape_intersect(shape<VALTYPE> & __sh1, shape<VALTYPE> &__sh2, const STRING & __name = ""):
      shape<VALTYPE>(__name)
    { sh1 = &__sh1; sh2 = &__sh2; }

    shape_intersect(const shape_intersect<VALTYPE> & s) :
      shape<VALTYPE>(s.name)
    {
      sh1 = s.sh1;
      sh2 = s.sh2;
    }

    virtual bool within(const v3d & r) const
    { return sh1->within(r) && sh2->within(r); }
 
    virtual VALTYPE volume() const
    { return 0; }

    virtual void scale(VALTYPE s)
    { sh1 -> scale(s); sh2 -> scale(s); }

    virtual void move(const v3d & v)
    { sh1 -> move(v); sh2 -> move(v);}

    virtual void rotate(const matrix3d<VALTYPE> & Rot)
    { sh1->rotate(Rot); sh2->rotate(Rot); }

    virtual v3d rmin() const
    { 
      v3d r, r1 = sh1->rmin(), r2 = sh2->rmin();
      for (int i=0; i<3; i++)
	r(i) = std::max(r1(i),r2(i));
      return r;
    }

    virtual v3d rmax() const
    { 
      v3d r, r1 = sh1->rmax(), r2 = sh2->rmax();
      for (int i=0; i<3; i++)
	r(i) = std::min(r1(i),r2(i));
      return r;
    }

    virtual v3d fmin(const periodic_cell<VALTYPE> &v) const
    {
      v3d f, f1 = sh1->fmin(v), f2 = sh2->fmin(v);
      for (int i=0; i<3; i++)
	f(i) = std::max(f1(i),f2(i));
      return f;
    }

    virtual v3d fmax(const periodic_cell<VALTYPE> &v) const
    {
      v3d f, f1 = sh1->fmax(v), f2 = sh2->fmax(v);
      for (int i=0; i<3; i++)
	f(i) = std::min(f1(i),f2(i));
      return f;
    }

    virtual void write(std::basic_ostream<CHAR,TRAITS> &os, int offset=0) const
    {
      for (int i=0; i<offset; i++) os << " ";
      os << "intersect";
      if (name != "")
	os << " " << name;
      os << "( shape1 = ";
      sh1->write(os);
      os << ", shape2 = ";
      sh2->write(os);
      os << ")";
    }

  };

  // ----------------------------------------------------------------

  template <class VALTYPE>
  class shape_subtract : public shape<VALTYPE>{

    shape<VALTYPE> *sh1, *sh2;

  public:
    using shape<VALTYPE>::name;

    shape_subtract(shape<VALTYPE> & __sh1, shape<VALTYPE> &__sh2, const STRING & __name = ""):
      shape<VALTYPE>(__name)
    { sh1 = &__sh1; sh2 = &__sh2; }

    shape_subtract(const shape_subtract<VALTYPE> & s) :
      shape<VALTYPE>(s.name)
    {
      sh1 = s.sh1;
      sh2 = s.sh2;
    }

    virtual bool within(const v3d & r) const
    { return sh1->within(r) && (! sh2->within(r)); }
 
    virtual VALTYPE volume() const
    { return 0; }

    virtual void scale(VALTYPE s)
    { sh1 -> scale(s); sh2 -> scale(s); }

    virtual void move(const v3d & v)
    { sh1 -> move(v); sh2 -> move(v);}

    virtual void rotate(const matrix3d<VALTYPE> & Rot)
    { sh1->rotate(Rot); sh2->rotate(Rot); }

    virtual v3d rmin() const
    { return sh1->rmin(); }

    virtual v3d rmax() const
    { return sh1->rmax(); }

    virtual v3d fmin(const periodic_cell<VALTYPE> &v) const
    { return sh1->fmin(v); }

    virtual v3d fmax(const periodic_cell<VALTYPE> &v) const
    { return sh1->fmax(v); }

    virtual void write(std::basic_ostream<CHAR,TRAITS> &os, int offset=0) const
    {
      for (int i=0; i<offset; i++) os << " ";
      os << "subtract";
      if (name != "")
	os << " " << name;
      os << "( shape1 = ";
      sh1->write(os);
      os << ", shape2 = ";
      sh2->write(os);
      os << ")";
    }

  };

  // ----------------------------------------------------------------

  template <class VALTYPE>
  class shape_invert : public shape<VALTYPE>{

    shape<VALTYPE> *sh;

  public:
    using shape<VALTYPE>::name;

    shape_invert(shape<VALTYPE> & __sh, const STRING & __name = ""):
      shape<VALTYPE>(__name)
    { sh = &__sh; }

    shape_invert(const shape_subtract<VALTYPE> & s) :
      shape<VALTYPE>(s.name)
    {
      sh = s.sh;
    }

    virtual bool within(const v3d & r) const
    { return ! sh->within(r); }
 
    virtual VALTYPE volume() const
    { return infty; }

    virtual void scale(VALTYPE s)
    { sh -> scale(s); }

    virtual void move(const v3d & v)
    { sh -> move(v);}

    virtual void rotate(const matrix3d<VALTYPE> & Rot)
    { sh->rotate(Rot); }

    virtual v3d rmin() const
    { return {-infty, -infty, -infty}; }

    virtual v3d rmax() const
    { return {infty, infty, infty}; }

    virtual v3d fmin(const periodic_cell<VALTYPE> &v) const
    { return {-infty, -infty, -infty}; }

    virtual v3d fmax(const periodic_cell<VALTYPE> &v) const
    { return  {infty, infty, infty}; }

    virtual void write(std::basic_ostream<CHAR,TRAITS> &os, int offset=0) const
    {
      for (int i=0; i<offset; i++) os << " ";
      os << "invert";
      if (name != "")
	os << " " << name;
      os << "( shape = ";
      sh->write(os);
      os << ")";
    }

  };

  // ----------------------------------------------------------------

  template <class VALTYPE>
  class shape_xor : public shape<VALTYPE>{

    shape<VALTYPE> *sh1, *sh2;

  public:
    using shape<VALTYPE>::name;

    shape_xor(shape<VALTYPE> & __sh1, shape<VALTYPE> &__sh2, const STRING & __name = ""):
      shape<VALTYPE>(__name)
    { sh1 = &__sh1; sh2 = &__sh2; }

    shape_xor(const shape_subtract<VALTYPE> & s) :
      shape<VALTYPE>(s.name)
    {
      sh1 = s.sh1;
      sh2 = s.sh2;
    }

    virtual bool within(const v3d & r) const
    { 
      bool in1 = sh1->within(r);
      bool in2 = sh2->within(r);
      return (in1 && ! in2) || (in2 && ! in1); 
    }
 
    virtual VALTYPE volume() const
    { return 0; }

    virtual void scale(VALTYPE s)
    { sh1 -> scale(s); sh2 -> scale(s); }

    virtual void move(const v3d & v)
    { sh1 -> move(v); sh2 -> move(v);}

    virtual void rotate(const matrix3d<VALTYPE> & Rot)
    { sh1->rotate(Rot); sh2->rotate(Rot); }

    virtual v3d rmin() const
    { 
      v3d r, r1 = sh1->rmin(), r2 = sh2->rmin();
      for (int i=0; i<3; i++)
	r(i) = std::min(r1(i),r2(i));
      return r;
    }

    virtual v3d rmax() const
    { 
      v3d r, r1 = sh1->rmax(), r2 = sh2->rmax();
      for (int i=0; i<3; i++)
	r(i) = std::max(r1(i),r2(i));
      return r;
    }

    virtual v3d fmin(const periodic_cell<VALTYPE> &v) const
    {
      v3d f, f1 = sh1->fmin(v), f2 = sh2->fmin(v);
      for (int i=0; i<3; i++)
	f(i) = std::min(f1(i),f2(i));
      return f;
    }

    virtual v3d fmax(const periodic_cell<VALTYPE> &v) const
    {
      v3d f, f1 = sh1->fmax(v), f2 = sh2->fmax(v);
      for (int i=0; i<3; i++)
	f(i) = std::max(f1(i),f2(i));
      return f;
    }

    virtual void write(std::basic_ostream<CHAR,TRAITS> &os, int offset=0) const
    {
      for (int i=0; i<offset; i++) os << " ";
      os << "xor";
      if (name != "")
	os << " " << name;
      os << "( shape1 = ";
      sh1->write(os);
      os << ", shape2 = ";
      sh2->write(os);
      os << ")";
    }

  };

  // ----------------------------------------------------------------

#ifdef PY_EXPORT

  template <class VALTYPE>
  struct py_shape : shape<VALTYPE>, bp::wrapper<shape<VALTYPE> >
  {

    bool within(const v3d & r) const
    { this->get_override("within")(r); }

    void scale(VALTYPE s)
    { this->get_override("scale")(s); }

    void move(const v3d & v)
    { this->get_override("move")(v);}

    void rotate(const matrix3d<VALTYPE> & Rot)
    { this->get_override("rotate")(Rot);}

    VALTYPE volume() const
    { this->get_override("volume")();}

    v3d rmin() const
    { this->get_override("rmin")();}

    v3d rmax() const
    { this->get_override("rmax")();}

    v3d fmin(const periodic_cell<VALTYPE> &v) const
    { this->get_override("fmin")(v);}

    v3d fmax(const periodic_cell<VALTYPE> &v) const
    { this->get_override("fmax")(v);}

    void write(std::basic_ostream<CHAR,TRAITS> &os, int offset=0) const
    { this->get_override("write")(os,offset);}

  };

#endif

};

#undef v3d
#undef v2d

#endif
