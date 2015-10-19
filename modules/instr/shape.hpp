#ifndef _QPP_SHAPE_H
#define _QPP_SHAPE_H

#include <mathf/constants.hpp>
#include <data/qppdata.hpp>
#include <geom/geom.hpp>
#include <lace/lace3d.hpp>
#include <algorithm>

#define v2d lace::vector2d<VALTYPE>
#define v3d lace::vector3d<VALTYPE>

namespace qpp{

  // ------------------ 3D primitives prototype ------------------------

  template <class VALTYPE>
  class qpp_shape : public qpp_object{

  public:

    qpp_shape(const STRING & __name = "", qpp_object * __owner = NULL):
      qpp_object(__name,__owner)
    {}

    virtual bool within(const v3d & r) const =0;
    // Answers the question whether the r point is situated within the shape

    virtual void scale(VALTYPE s) =0;
    virtual void move(const v3d & v) =0;
    virtual void rotate(const lace::matrix3d<VALTYPE> & Rot) =0;

    // -------------------------------------
    
    virtual VALTYPE volume() const =0;
    virtual v3d rmin() const =0;
    virtual v3d rmax() const =0;
    // Minimal & maximal cartesian coordinates of the shape

    virtual lace::simple_vector<VALTYPE,3> fmin(const periodic_cell<3,VALTYPE> &v) const =0;
    virtual lace::simple_vector<VALTYPE,3> fmax(const periodic_cell<3,VALTYPE> &v) const =0;
    // Minimal and maximal fractional coordinates of the shape for given cell

    virtual lace::simple_vector<VALTYPE,2> fmin(const periodic_cell<2,VALTYPE> &v) const
    {
      periodic_cell<3,VALTYPE> xv(v(0),v(1),v(0)%v(1));
      v3d f = fmin(xv);
      lace::simple_vector<VALTYPE,2> res;
      res(0) = f(0);
      res(1) = f(1);
      return res;
    }

    virtual lace::simple_vector<VALTYPE,1> fmin(const periodic_cell<1,VALTYPE> &v) const
    {
      v3d a = v(0), arb(1,0,0), b, c;
      if (norm(a%arb)/norm(a)<v3d::tolerance_for_equiv)
	arb=v3d(0,1,0);
      b = a%arb/norm(a%arb);
      c = a%b/norm(a%b);
      periodic_cell<3,VALTYPE> xv(a,b,c);
      v3d f = fmin(xv);
      lace::simple_vector<VALTYPE,1> res;
      res(0) = f(0);
      return res;
    }

    virtual lace::simple_vector<VALTYPE,2> fmax(const periodic_cell<2,VALTYPE> &v) const
    {
      periodic_cell<3,VALTYPE> xv(v(0),v(1),v(0)%v(1));
      v3d f = fmax(xv);
      lace::simple_vector<VALTYPE,2> res;
      res(0) = f(0);
      res(1) = f(1);
      return res;
    }

    virtual lace::simple_vector<VALTYPE,1> fmax(const periodic_cell<1,VALTYPE> &v) const
    {
      v3d a = v(0), arb(1,0,0), b, c;
      if (norm(a%arb)/norm(a)<v3d::tolerance_for_equiv)
	arb=v3d(0,1,0);
      b = a%arb/norm(a%arb);
      c = a%b/norm(a%b);
      periodic_cell<3,VALTYPE> xv(a,b,c);
      v3d f = fmax(xv);
      lace::simple_vector<VALTYPE,1> res;
      res(0) = f(0);
      return res;
    }

    // Minimal & maximal fractional coordinates of the shape for given translation vectors v

    virtual int n_nested() const
    { return 0;}

    virtual qpp_object* nested(int i) const
    { return NULL;}

    virtual int n_param() const
    { 
      //fixme - implement this for different shapes
      return 0; 
    }
    
    virtual void set_n_param(int n)
    {}


    virtual qppobject_type gettype() const
    { return qtype_shape;}
    
  };

  // ------------------ 3D primitives ------------------------

  template <class VALTYPE>
  class shape_box : public qpp_shape<VALTYPE>{

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

    using qpp_shape<VALTYPE>::name;

    shape_box(const v3d & a1, const v3d & a2, const v3d & a3, const v3d & r0, 
		       const STRING & __name = "", qpp_object * __owner = NULL) :
      qpp_shape<VALTYPE>(__name, __owner)
    {
      crn = r0;
      a[0] = a1;
      a[1] = a2;
      a[2] = a3;
    }

    shape_box(const v3d & a1, const v3d & a2, const v3d & a3, 
		       const STRING & __name = "", qpp_object * __owner = NULL):
      qpp_shape<VALTYPE>(__name, __owner)
    {
       crn = 0e0;
      a[0] = a1;
      a[1] = a2;
      a[2] = a3;
    }
    
    shape_box(VALTYPE a1, VALTYPE a2, VALTYPE a3, 
		       const STRING & __name = "", qpp_object * __owner = NULL):
      qpp_shape<VALTYPE>(__name, __owner)
    {
      crn = 0e0;
      a[0] = v3d(a1,  0e0, 0e0);
      a[1] = v3d(0e0, a2,  0e0);
      a[2] = v3d(0e0, 0e0, a3);
    }

    shape_box(const shape_box<VALTYPE> & s)
    {
      crn = s.crn;
      a[0] = s.a[0];
      a[1] = s.a[1];
      a[2] = s.a[2];
    }

    virtual STRING category() const
    { return "box";}    

    virtual qpp_object * copy() const
    {
      return new shape_box<VALTYPE>(*this);
    }

    virtual void write(std::basic_ostream<CHAR,TRAITS> &os, int offset=0) const
    {
      for (int i=0; i<offset; i++)
	os << " ";
      os << "box";
      if (name() != "")
	os << " " << name();
      os << "( a" << a[0] << ", b" << a[1] << ", c" << a[2] << ", corner" << crn << ");\n";
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

    virtual lace::simple_vector<VALTYPE,3> fmin(const periodic_cell<3,VALTYPE> &v) const
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

    virtual lace::simple_vector<VALTYPE,3> fmax(const periodic_cell<3,VALTYPE> &v) const
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

    virtual void rotate(const lace::matrix3d<VALTYPE> & Rot)
    {
      crn  = Rot*crn;
      a[0] = Rot*a[0];
      a[1] = Rot*a[1];
      a[2] = Rot*a[2];
    }

  };

  // ----------------------------------------------------------------

  template <class VALTYPE>
  class shape_sphere : public qpp_shape<VALTYPE>{

    VALTYPE R;
    v3d r0;

  public:

    using qpp_shape<VALTYPE>::name;

    shape_sphere(VALTYPE _R, const STRING & __name = "", qpp_object * __owner = NULL):
      qpp_shape<VALTYPE>(__name, __owner)
    {
      R = _R;
      r0 = 0e0;
    }

    shape_sphere(VALTYPE _R, const v3d & _r0, 
		     const STRING & __name = "", qpp_object * __owner = NULL):
      qpp_shape<VALTYPE>(__name, __owner)
    {
      R = _R;
      r0 = _r0;
    }

    shape_sphere(const shape_sphere<VALTYPE> & s):
      qpp_shape<VALTYPE>(s)
    {
      R = s.R;
      r0 = s.r0;
    }
    /*
    qpp_shape_sphere(qpp_param_array & parm, 
		     const STRING & __name = "", qpp_object * __owner = NULL): 
      qpp_shape<VALTYPE>(__name, __owner)
    {      
    }
    */
    virtual STRING category() const
    { return "sphere";}

    virtual void write(std::basic_ostream<CHAR,TRAITS> &os, int offset=0) const
    {
      for (int i=0; i<offset; i++)
	os << " ";
      os << "sphere";
      if (name() != "")
	os << " " << name();
      os << "( R=" << R << ", center" << r0 << ");\n";
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

    virtual lace::simple_vector<VALTYPE,3> fmin(const periodic_cell<3,VALTYPE> &v) const
    {
      lace::matrix3d<VALTYPE> A(v(0),v(1),v(2));
      lace::matrix3d<VALTYPE> B = invert(A);

      v3d res = 0e0;
      for (int i=0; i<3; i++)
	for (int j=0; j<3; j++)
	  res(i) += B(i,j)*B(i,j);
      for (int i=0; i<3; i++)
	res(i) = std::sqrt(res(i));
      return B*r0 - R*res;
    }

    virtual lace::simple_vector<VALTYPE,3> fmax(const periodic_cell<3,VALTYPE> &v) const
    {      
      lace::matrix3d<VALTYPE> A(v(0),v(1),v(2));
      lace::matrix3d<VALTYPE> B = invert(A);

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

    virtual void rotate(const lace::matrix3d<VALTYPE> & Rot)
    {
      r0 = Rot*r0;
    }

    virtual qpp_object * copy() const
    {
      return new shape_sphere<VALTYPE>(*this);
    }

  };

  // ----------------------------------------------------------------

  template <class VALTYPE>
  class shape_union : public qpp_shape<VALTYPE>{

    qpp_shape<VALTYPE> *sh1, *sh2;

  public:
    using qpp_shape<VALTYPE>::name;

    shape_union(qpp_shape<VALTYPE> & __sh1, qpp_shape<VALTYPE> &__sh2,
		const STRING & __name = "", qpp_object * __owner = NULL):
      qpp_shape<VALTYPE>(__name,__owner)
    { sh1 = &__sh1; sh2 = &__sh2; }

    shape_union(const shape_union<VALTYPE> & s) :
      qpp_shape<VALTYPE>(s)
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

    virtual void rotate(const lace::matrix3d<VALTYPE> & Rot)
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

    virtual lace::simple_vector<VALTYPE,3> fmin(const periodic_cell<3,VALTYPE> &v) const
    {
      v3d f, f1 = sh1->fmin(v), f2 = sh2->fmin(v);
      for (int i=0; i<3; i++)
	f(i) = std::min(f1(i),f2(i));
      return f;
    }

    virtual lace::simple_vector<VALTYPE,3> fmax(const periodic_cell<3,VALTYPE> &v) const
    {
      v3d f, f1 = sh1->fmax(v), f2 = sh2->fmax(v);
      for (int i=0; i<3; i++)
	f(i) = std::max(f1(i),f2(i));
      return f;
    }

    virtual STRING category() const
    { return "union";}

    virtual void write(std::basic_ostream<CHAR,TRAITS> &os, int offset=0) const
    {
      for (int i=0; i<offset; i++) os << " ";
      os << "union";
      if (name() != "")
	os << " " << name();
      os << "( shape1 = ";
      sh1->write(os);
      os << ", shape2 = ";
      sh2->write(os);
      os << ");\n";
    }

    virtual qpp_object * copy() const
    {
      return new shape_union<VALTYPE>(*this);
    }

  };

  // ----------------------------------------------------------------

  template <class VALTYPE>
  class shape_intersect : public qpp_shape<VALTYPE>{

    qpp_shape<VALTYPE> *sh1, *sh2;

  public:
    using qpp_shape<VALTYPE>::name;

    shape_intersect(qpp_shape<VALTYPE> & __sh1, qpp_shape<VALTYPE> &__sh2,
		    const STRING & __name = "", qpp_object * __owner = NULL):
      qpp_shape<VALTYPE>(__name,__owner)
    { sh1 = &__sh1; sh2 = &__sh2; }

    shape_intersect(const shape_intersect<VALTYPE> & s) :
      qpp_shape<VALTYPE>(s)
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

    virtual void rotate(const lace::matrix3d<VALTYPE> & Rot)
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

    virtual lace::simple_vector<VALTYPE,3> fmin(const periodic_cell<3,VALTYPE> &v) const
    {
      v3d f, f1 = sh1->fmin(v), f2 = sh2->fmin(v);
      for (int i=0; i<3; i++)
	f(i) = std::max(f1(i),f2(i));
      return f;
    }

    virtual lace::simple_vector<VALTYPE,3> fmax(const periodic_cell<3,VALTYPE> &v) const
    {
      v3d f, f1 = sh1->fmax(v), f2 = sh2->fmax(v);
      for (int i=0; i<3; i++)
	f(i) = std::min(f1(i),f2(i));
      return f;
    }

    virtual STRING category() const
    { return "intersect";}

    virtual void write(std::basic_ostream<CHAR,TRAITS> &os, int offset=0) const
    {
      for (int i=0; i<offset; i++) os << " ";
      os << "intersect";
      if (name() != "")
	os << " " << name();
      os << "( shape1 = ";
      sh1->write(os);
      os << ", shape2 = ";
      sh2->write(os);
      os << ");\n";
    }

    virtual qpp_object * copy() const
    {
      return new shape_intersect<VALTYPE>(*this);
    }

  };

  // ----------------------------------------------------------------

  template <class VALTYPE>
  class shape_subtract : public qpp_shape<VALTYPE>{

    qpp_shape<VALTYPE> *sh1, *sh2;

  public:
    using qpp_shape<VALTYPE>::name;

    shape_subtract(qpp_shape<VALTYPE> & __sh1, qpp_shape<VALTYPE> &__sh2,
		   const STRING & __name = "", qpp_object * __owner = NULL):
      qpp_shape<VALTYPE>(__name,__owner)
    { sh1 = &__sh1; sh2 = &__sh2; }

    shape_subtract(const shape_subtract<VALTYPE> & s) :
      qpp_shape<VALTYPE>(s)
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

    virtual void rotate(const lace::matrix3d<VALTYPE> & Rot)
    { sh1->rotate(Rot); sh2->rotate(Rot); }

    virtual v3d rmin() const
    { return sh1->rmin(); }

    virtual v3d rmax() const
    { return sh1->rmax(); }

    virtual lace::simple_vector<VALTYPE,3> fmin(const periodic_cell<3,VALTYPE> &v) const
    { return sh1->fmin(v); }

    virtual lace::simple_vector<VALTYPE,3> fmax(const periodic_cell<3,VALTYPE> &v) const
    { return sh1->fmax(v); }

    virtual STRING category() const
    { return "subtract";}

    virtual void write(std::basic_ostream<CHAR,TRAITS> &os, int offset=0) const
    {
      for (int i=0; i<offset; i++) os << " ";
      os << "subtract";
      if (name() != "")
	os << " " << name();
      os << "( shape1 = ";
      sh1->write(os);
      os << ", shape2 = ";
      sh2->write(os);
      os << ");\n";
    }

    virtual qpp_object * copy() const
    {
      return new shape_subtract<VALTYPE>(*this);
    }

  };

  // ----------------------------------------------------------------

};

#undef v3d
#undef v2d

#endif