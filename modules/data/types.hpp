#ifndef _QPP_TYPES_H
#define _QPP_TYPES_H

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <complex>

#include <stdexcept>

namespace qpp{

#ifdef WCHAR

#define CHAR wchar
#define STRING std::wstring

#else

#define CHAR char
#define STRING std::string

#endif

#undef assert

  typedef std::char_traits<CHAR> TRAITS;

  typedef std::basic_ostream<CHAR,TRAITS> OSTREAM;

  typedef std::basic_istream<CHAR,TRAITS> ISTREAM;

  typedef std::basic_ifstream<CHAR,TRAITS> IFSTREAM;

  typedef std::basic_ofstream<CHAR,TRAITS> OFSTREAM;

  typedef std::basic_stringstream<CHAR,TRAITS> SSTREAM;

  // -----------------------------------------------------

  template<typename VALTYPE>
  struct numeric_type;

  template<>
  struct numeric_type<float>
  {
    typedef float real;
    typedef std::complex<float> complex;
    typedef float norm;
    static bool is_real();
    static bool is_complex();
    static float I();
  };

  template<>
  struct numeric_type<double>
  {
    typedef double real;
    typedef std::complex<double> complex;
    typedef double norm;
    static bool is_real();
    static bool is_complex();
    static double I();
  };

  template<>
  struct numeric_type<long double>
  {
    typedef long double real;
    typedef std::complex<long double> complex;
    typedef long double norm;
    static bool is_real();
    static bool is_complex();
    static long double I();
  };

  template<>
  struct numeric_type<std::complex<float> >
  {
    typedef float real;
    typedef std::complex<float> complex;
    typedef float norm;
    static bool is_real();
    static bool is_complex();
    static std::complex<float> I();
  };

  template<>
  struct numeric_type<std::complex<double> >
  {
    typedef double real;
    typedef std::complex<double> complex;
    typedef double norm;
    static bool is_real();
    static bool is_complex();
    static std::complex<double> I();
  };

  template<>
  struct numeric_type<std::complex<long double> >
  {
    typedef long double real;
    typedef std::complex<long double> complex;
    typedef long double norm;
    static bool is_real();
    static bool is_complex();
    static std::complex<long double> I();
  };

  template<>
  struct numeric_type<int>
  {
    typedef int real;
    typedef std::complex<double> complex;
    typedef double norm;
    static bool is_real();
    static bool is_complex();
    static int I();
  };

  // --------------------------------------------------------

  enum basic_types{type_int, type_real, type_double, type_float, type_bool, type_string};

  // --------------------------------------------------------

  template <class T1, class T2>
  struct convert
  {
    static T1 get(T2 t)
    { 
      // fixme - error handling 
      throw std::invalid_argument("Type error in type cast");
    }
  };

  template<class T>
  struct convert<T,T>
  {
    static T get(T t)
    { return t; }
  };

  template<> struct convert<std::string, const char*>{ static std::string get(const char* t);};
  template<> struct convert<bool,short>{ static bool get(short t);};
  template<> struct convert<bool&,short&>{ static bool & get(short & t);};
  template<> struct convert<double,float>{ static double get(float t);};
  template<> struct convert<double,int>{ static double get(int t);};
  template<> struct convert<float,double>{ static float get(double t);};
  template<> struct convert<float,int>{ static float get(int t);};

  // --------------------------------------------------------

  template <class T>
  struct attributes;

  template <>
  struct attributes<bool>
  {
    static STRING name;
    static basic_types type;
  };

  template <>
  struct attributes<int>
  {
    static STRING name;
    static basic_types type;
  };

  template <>
  struct attributes<float>
  {
    static STRING name;
    static basic_types type;
  };

  template <>
  struct attributes<double>
  {
    static STRING name;
    static basic_types type;
  };

  template <>
  struct attributes<STRING>
  {
    static STRING name;
    static basic_types type;
  };


};

#endif
