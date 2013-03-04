#ifndef _LACE3D_H
#define _LACE3D_H

#include <cmath>
#include <ostream>
#include <sstream>
#include <boost/format.hpp>
#include <lace/complex.hpp>

namespace lace{

  template<class VALTYPE = double>
  class vector3d{
    VALTYPE r[3];

  public:

    static double tolerance_for_equiv;

    vector3d(){}

    vector3d(VALTYPE s)
    {
      r[0] = s; r[1] = s; r[2] = s;
    }
    
    vector3d(VALTYPE x, VALTYPE y, VALTYPE z)
    {
      r[0] = x; r[1] = y; r[2] = z;
    }

    vector3d(const vector3d<VALTYPE>& v)
    {
      r[0] = v.r[0]; r[1] = v.r[1]; r[2] = v.r[2];
    }

    inline VALTYPE& operator()(int i)
    {
      return r[i];
    }

    inline VALTYPE& x(){return r[0];}

    inline VALTYPE& y(){return r[1];}

    inline VALTYPE& z(){return r[2];}

    inline VALTYPE norm2(){return r[0]*r[0]+r[1]*r[1]+r[2]*r[2];}

    inline VALTYPE norm(){return std::sqrt(norm2());}

    inline vector3d<VALTYPE> operator+(vector3d<VALTYPE> v)
    {
      vector3d<VALTYPE> res( x()+v.x(), y()+v.y(), z()+v.z() );
      return res;
    }

    inline vector3d<VALTYPE> operator-(vector3d<VALTYPE> v)
    {
      vector3d<VALTYPE> res( x()-v.x(), y()-v.y(), z()-v.z() );
      return res;
    }

    inline vector3d<VALTYPE> operator*(VALTYPE s)
    {
      vector3d<VALTYPE> res( x()*s, y()*s, z()*s );
      return res;
    }

    inline vector3d<VALTYPE> operator/(VALTYPE s)
    {
      vector3d<VALTYPE> res( x()/s, y()/s, z()/s );
      return res;
    }

    inline vector3d<VALTYPE> operator%(vector3d<VALTYPE> v)
    {
      vector3d<VALTYPE> res( y()*v.z() - z()*v.y(), 
			     z()*v.x() - x()*v.z(), 
			     x()*v.y() - y()*v.x() );
      return res;
    }

    inline vector3d<VALTYPE>& operator=(vector3d<VALTYPE> v)
    {
      x() = v.x();
      y() = v.y();
      z() = v.z();
      return *this;
    }

    inline vector3d<VALTYPE>& operator+=(vector3d<VALTYPE> v)
    {
      x() += v.x();
      y() += v.y();
      z() += v.z();
      return *this;
    }

    inline vector3d<VALTYPE>& operator-=(vector3d<VALTYPE> v)
    {
      x() -= v.x();
      y() -= v.y();
      z() -= v.z();
      return *this;
    }

    inline vector3d<VALTYPE>& operator*=(VALTYPE s)
    {
      x() *= s;
      y() *= s;
      z() *= s;
      return *this;
    }

    inline vector3d<VALTYPE> operator-()
    {
      return vector3d<VALTYPE>(-x(), -y(), -z());
    }

    inline bool operator==(vector3d<VALTYPE> b)
    {
      return (*this - b).norm() <=  tolerance_for_equiv;
    }

  };

  template<class VALTYPE>
  inline vector3d<VALTYPE> operator*(VALTYPE s, vector3d<VALTYPE> v)
  {
    return v.operator*(s);
  }

  template<class VALTYPE>
  inline vector3d<VALTYPE> operator*(int s, vector3d<VALTYPE> v)
  {
    return v.operator*((VALTYPE)s);
  }

  template<class VALTYPE>
  inline VALTYPE scal(vector3d<VALTYPE> v1, vector3d<VALTYPE> v2)
  {
    return v1(0)*v2(0) + v1(1)*v2(1) + v1(2)*v2(2);
  }
 
  template<class VALTYPE>
  inline VALTYPE norm2(vector3d<VALTYPE> v)
  {
    return v(0)*v(0) + v(1)*v(1) + v(2)*v(2);
  }

  template<class VALTYPE>
  inline VALTYPE norm(vector3d<VALTYPE> v)
  {
    return std::sqrt(norm2(v));
  }

 /*---------------------------------------------------------*/

  template<class VALTYPE = double>
  class matrix3d{
    VALTYPE a[9];

  public:

    static double tolerance_for_equiv;

    matrix3d(){}

    matrix3d(double axx, double axy, double axz,
	     double ayx, double ayy, double ayz,
	     double azx, double azy, double azz)
    {
      a[0] = axx; a[1] = axy; a[2] = axz;
      a[3] = ayx; a[4] = ayy; a[5] = ayz;
      a[6] = azx; a[7] = azy; a[8] = azz;
    }

    matrix3d(VALTYPE alpha)
    {
      a[0] = a[4] = a[8] = alpha;
      a[1] = a[2] = a[3] = a[5] = a[6] = a[7] = VALTYPE(0);      
    }

    matrix3d(vector3d<VALTYPE> a1, vector3d<VALTYPE> b, vector3d<VALTYPE> c)
    {
      a[0] = a1.x(); a[1] = b.x(); a[2] = c.x();
      a[3] = a1.y(); a[4] = b.y(); a[5] = c.y();
      a[6] = a1.z(); a[7] = b.z(); a[8] = c.z();
    }

    inline VALTYPE & operator()(int i, int j){ return a[i*3+j];}

    //------------------------------------

    inline VALTYPE& xx(){return a[0];}

    inline VALTYPE& xy(){return a[1];}

    inline VALTYPE& xz(){return a[2];}

    inline VALTYPE& yx(){return a[3];}

    inline VALTYPE& yy(){return a[4];}

    inline VALTYPE& yz(){return a[5];}

    inline VALTYPE& zx(){return a[6];}

    inline VALTYPE& zy(){return a[7];}

    inline VALTYPE& zz(){return a[8];}

    inline vector3d<VALTYPE> operator()(int i)
    {
      return vector3d<VALTYPE>(a[i],a[i+3],a[i+6]);
    } 

    //------------------------------------

    inline vector3d<VALTYPE> operator*(vector3d<VALTYPE> & x)
    {
      
      vector3d<VALTYPE> res( a[0]*x.x() + a[1]*x.y() + a[2]*x.z(), 
			     a[3]*x.x() + a[4]*x.y() + a[5]*x.z(),
			     a[6]*x.x() + a[7]*x.y() + a[8]*x.z());
      return res;
    }

    inline matrix3d<VALTYPE> operator*(matrix3d<VALTYPE> b)
    {
      // fixme -- not efficient!!!
      
      matrix3d<VALTYPE> res;

      for(int i=0; i<3; i++)
	for(int j=0; j<3; j++)
	  {
	    res(i,j) = 0e0;
	    for (int k=0; k<3; k++)
	      res(i,j) += (*this)(i,k)*b(k,j);
	  }
      return res;
    }

    inline matrix3d<VALTYPE> T()
    {
      return matrix3d<VALTYPE>( a[0], a[3], a[6],
				a[1], a[4], a[7],
				a[2], a[5], a[8]);
    }

    inline matrix3d<VALTYPE>& operator=(VALTYPE alpha)
    {
      a[0] = a[4] = a[8] = VALTYPE(alpha);
      a[1] = a[2] = a[3] = a[5] = a[6] = a[7] = VALTYPE(0);
      return *this;
    }

    inline matrix3d<VALTYPE> operator+(matrix3d<VALTYPE> b)
    {
      return  matrix3d<VALTYPE>( a[0]+b(0,0), a[1]+b(0,1), a[2]+b(0,2),
				 a[3]+b(1,0), a[4]+b(1,1), a[5]+b(1,2),
				 a[6]+b(2,0), a[7]+b(2,1), a[8]+b(2,2));
    } 

    inline matrix3d<VALTYPE> operator-(matrix3d<VALTYPE> b)
    {
      return  matrix3d<VALTYPE>( a[0]-b(0,0), a[1]-b(0,1), a[2]-b(0,2),
				 a[3]-b(1,0), a[4]-b(1,1), a[5]-b(1,2),
				 a[6]-b(2,0), a[7]-b(2,1), a[8]-b(2,2));
    } 

    inline matrix3d<VALTYPE> operator*(VALTYPE alpha)
    {
      return  matrix3d<VALTYPE>( a[0]*alpha, a[1]*alpha, a[2]*alpha,
				 a[3]*alpha, a[4]*alpha, a[5]*alpha,
				 a[6]*alpha, a[7]*alpha, a[8]*alpha);
    }

    inline matrix3d<VALTYPE> operator*=(VALTYPE alpha)
    {
      for (int i = 0; i<9; i++)
	a[i] *= alpha;
      return *this;
    }

    inline double norm()
    {
      double s = 0e0;
      for (int i = 0; i<9; i++)
	s += a[i]*a[i];
      return std::sqrt(s);
    }

    inline bool operator==(matrix3d<VALTYPE> b)
    {
      return ((*this) - b ).norm() <= tolerance_for_equiv;
    }

    inline VALTYPE det()
    {
      return 
	xx()*( yy()*zz() - yz()*zy() )-
	xy()*( yx()*zz() - zx()*yz() )+
	xz()*( yx()*zy() - zx()*yy() );
    }

  };
  
  template<class VALTYPE>
  inline VALTYPE det(matrix3d<VALTYPE> a)
  {
    return a.det();
  }

  template<class VALTYPE>
  inline VALTYPE det(vector3d<VALTYPE> a, vector3d<VALTYPE> b, vector3d<VALTYPE> c)
  {
    return 
	a.x()*( b.y()*c.z() - b.z()*c.y() )-
	a.y()*( b.x()*c.z() - c.x()*b.z() )+
	a.z()*( b.x()*c.y() - c.x()*b.y() );
  }

  template<class VALTYPE>
  inline double norm(matrix3d<VALTYPE> a)
  {
    return a.norm();
  }

  template<class VALTYPE>
  inline matrix3d<VALTYPE> operator*(VALTYPE alpha, matrix3d<VALTYPE> a)
  {
    return a*alpha;
  }

  template<class VALTYPE>
  double matrix3d<VALTYPE>::tolerance_for_equiv = 1e-6;

  template<class VALTYPE>
  double vector3d<VALTYPE>::tolerance_for_equiv = 1e-6;

  //-------------------------------------------------

  template<class VALTYPE>
  matrix3d<VALTYPE> RotMtrx(vector3d<VALTYPE> n, VALTYPE phi)
  {
    n = n/norm(n);

    VALTYPE c = std::cos(phi), s = std::sin(phi);
    return matrix3d<VALTYPE>
      ( n(0)*n(0)*(1-c) + c,       n(0)*n(1)*(1-c) + n(2)*s,  n(0)*n(2)*(1-c) - n(1)*s,
	n(1)*n(0)*(1-c) - n(2)*s,  n(1)*n(1)*(1-c) + c,       n(1)*n(2)*(1-c) + n(0)*s,
	n(2)*n(0)*(1-c) + n(1)*s,  n(2)*n(1)*(1-c) - n(0)*s,  n(2)*n(2)*(1-c) + c    );
  }

  template<class VALTYPE>
  matrix3d<VALTYPE> Sigma(vector3d<VALTYPE> n)
  {
    n = n/norm(n);
    return matrix3d<VALTYPE>
      ( 1e0 - 2*n(0)*n(0), -2*n(0)*n(1),      -2*n(0)*n(2),
	-2*n(0)*n(1),      1e0 - 2*n(1)*n(1), -2*n(1)*n(2),
	-2*n(0)*n(2),      -2*n(1)*n(2),      1e0 - 2*n(2)*n(2) );
  } 

  template<class VALTYPE>
  vector3d<VALTYPE> solve3d(matrix3d<VALTYPE> A, vector3d<VALTYPE> b)
  {
    VALTYPE 
      D = det(A),
      X = det(b,    A(1), A(2)),
      Y = det(A(0), b,    A(2)),
      Z = det(A(0), A(1), b);
    return vector3d<VALTYPE>(X, Y, Z)/D;
  }

  template<class VALTYPE>
  vector3d<VALTYPE> solve3d(vector3d<VALTYPE> A0, vector3d<VALTYPE> A1, 
			    vector3d<VALTYPE> A2, vector3d<VALTYPE> b)
  {
    VALTYPE 
      D = det(A0, A1, A2),
      X = det(b,  A1, A2),
      Y = det(A0, b,  A2),
      Z = det(A0, A1, b);
    return vector3d<VALTYPE>(X, Y, Z)/D;
  }

  //-------------------------------------------------

  template<class VALTYPE>
  vector3d<typename lace::numeric_type<VALTYPE>::complex> solve_cubeq(VALTYPE a, VALTYPE b, VALTYPE c, VALTYPE d)
  // Solves ax^3 + bx^2 + cx + d = 0
  {
    b /= a;
    c /= a;
    d /= a;
    typename lace::numeric_type<VALTYPE>::complex I(VALTYPE(0),VALTYPE(1));
    VALTYPE aa = c/3 - b*b/9, bb = d/2 - b*c/6 + b*b*b/27;
    VALTYPE discr1 = aa*aa*aa + bb*bb;
    typename lace::numeric_type<VALTYPE>::complex discr2;
    if (discr1 > VALTYPE(0))
      {
	VALTYPE discr23 = std::sqrt(discr1) - bb;
	int sgn = (discr23>0)? 1 : -1;
	discr23 *= sgn;
	discr2 = sgn*std::exp(std::log(discr23)/3);
      }
    else
      {
	discr2 = -bb + std::sqrt(-discr1)*I;
	VALTYPE dabs = lace::abs(discr2), darg = lace::arg(discr2);
	discr2 = lace::polar(std::exp(std::log(dabs)/3),darg/3);
      }
    typename lace::numeric_type<VALTYPE>::complex r1 = discr2 - aa/discr2, r2 = discr2 + aa/discr2;
    return vector3d<typename lace::numeric_type<VALTYPE>::complex>( r1 - b/3, 
								   -r1/2 + std::sqrt(3.0)*r2*I/2 - b/3, 
								   -r1/2 - std::sqrt(3.0)*r2*I/2 - b/3);
  }

  //-------------------------------------------------

  //template<class VALTYPE>
  //void diagon

  //-------------------------------------------------

  template<typename _CharT, class _Traits, class VALTYPE>
  std::basic_ostream<_CharT, _Traits>&
  operator<<(std::basic_ostream<_CharT, _Traits>& __os, vector3d<VALTYPE> v)
  {
    std::basic_ostringstream<_CharT, _Traits> __s;
    __s.flags(__os.flags());
    __s.imbue(__os.getloc());
    __s.precision(__os.precision());
    __s << "(" << v(0) << "," << v(1) << "," << v(2) << ")";
    return __os << __s.str();
  }

  template<typename _CharT, class _Traits, class VALTYPE>
  std::basic_ostream<_CharT, _Traits>&
  operator<<(std::basic_ostream<_CharT, _Traits>& __os, matrix3d<VALTYPE> a)
  {
    std::basic_ostringstream<_CharT, _Traits> __s;
    __s.flags(__os.flags());
    __s.imbue(__os.getloc());
    __s.precision(__os.precision());
    __s  << a(0,0) << " " << a(0,1) << " " << a(0,2) << "\n"
	 << a(1,0) << " " << a(1,1) << " " << a(1,2) << "\n"
	 << a(2,0) << " " << a(2,1) << " " << a(2,2) << "\n" ;
    return __os << __s.str();
  }

};

#endif