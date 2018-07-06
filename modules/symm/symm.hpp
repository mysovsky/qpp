#ifndef _QPP_SYMM_H
#define _QPP_SYMM_H

#include <symm/index.hpp>
#include <vector>

#ifdef PY_EXPORT
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
namespace py = pybind11;
#endif

namespace qpp{
  
  /*!\brief The generators_pack class implements Positionary Generator Form (PGF) for arbitrary finite
    group.
  */
  template <class TRANSFORM>
  class generators_pack{
  public:
    
    std::vector<TRANSFORM> generators;
    index _begin, _end;
    int DIM;
    
    generators_pack(int dim=0)
    {
      DIM=dim;
      generators.resize(DIM);
      _begin = index::D(DIM);
      _end   = index::D(DIM);
    }

    generators_pack(const std::vector<TRANSFORM> & g,
                    const index & __begin, const index & __end)
    {
      DIM = g.size();
      generators.resize(DIM);
      int d=0;
      for (const TRANSFORM & t : g)
        generators[d++]=t;
      _begin = __begin;
      _end   = __end;
    }

    generators_pack(const std::vector<TRANSFORM> & g)
    {
      DIM = g.size();
      generators.resize(DIM);

      int d=0;
      for (const TRANSFORM & t : g)
        generators[d++]=t;

      _begin = index::D(DIM);
      _end   = index::D(DIM);
    }
    
    generators_pack(const generators_pack<TRANSFORM> & G) :
      DIM(G.DIM), generators(G.generators), _begin(G._begin), _end(G._end)
    {}

    int get_dim(){return DIM;}

    void set_dim(int D){
      DIM = D;
      generators.resize(DIM);
      
      _begin = index::D(DIM);
      _end   = index::D(DIM);
    }

    TRANSFORM operator()(const index & n) const{
      if (DIM==0)
        return TRANSFORM::Identity();

      //TRANSFORM A = pow(generators[0],n(0));
      TRANSFORM A = generators[0].pow(n(0));
      for (int d = 1; d<DIM; d++)
        //A = A*pow(generators[d],n(d));
        A = A * generators[d].pow(n(d));
      return A;
    }

    template <class ARRAY>
    void generate(ARRAY & group){
      for (iterator n(_begin, _end); !n.end(); n++)
        group.push_back((*this)(n));
    }

    void auto_order(int d){
      _begin(d) = 0;
      const TRANSFORM & g = generators[d];
      TRANSFORM a = g;
      int n=1;
      while (a != TRANSFORM::Identity()){
          a = a*g;
          n++;
        }
      _end(d) = n-1;
    }

    void auto_orders(){
      for (int d=0; d<DIM; d++)
        auto_order(d);
    }

    inline index begin() const
    { return _begin;}

    inline index end() const
    { return _end;}
    
  };

  // ------------------------------------------------------------------------
  /*
#ifdef PY_EXPORT

  template<class TRANSFORM> class generated_group;

  template <class TRANSFORM>
  int py_group_len(const generated_group<TRANSFORM> & G)
  {
    return G.size();
  }

#endif
  */
  // ------------------------------------------------------------------------

  template<class TRANSFORM>
  class generated_group{
  public:
    
    std::vector<TRANSFORM> group;

    int index(const TRANSFORM & g){
      int i;
      bool result=false;
      for (i=0; i<group.size(); i++)
        if ( group[i] == g ){
            result = true;
            break;
          }
      return result? i : -1;
    }

    generated_group(TRANSFORM E = TRANSFORM::Identity()){
      group.push_back(E);
    }

    generated_group(const generated_group<TRANSFORM> & G):
      group(G.group)
    {}

    inline TRANSFORM & operator[](int i)
    { return group[i]; }

    inline TRANSFORM operator[](int i) const
    { return group[i]; }

    inline int size() const
    { return group.size(); }

    void add(const TRANSFORM & g){
      if ( index(g) >= 0 )
        return;
      int inew = size();
      group.push_back(g);

      while (inew < size()){
          //std::cout << size() << "\n";

          int inewest = size();
          for (int ig1 = 0; ig1 < inewest; ig1++)
            for (int ig2 = inew; ig2 < inewest; ig2++)
              {
                //std::cout << "ig1= " << ig1 << " ig2= " << ig2 << "\n";

                TRANSFORM h1 = group[ig1]*group[ig2];

                //std::cout << "h1= " << h1 << "\n";

                if (index(h1)==-1)
                  group.push_back(h1);
                TRANSFORM h2 = group[ig2]*group[ig1];

                //std::cout << "h2= " << h2 << "\n";

                if (h2 != h1 && index(h2)==-1)
                  group.push_back(h2);
              }
          //std::cout << inew << " " << inewest << "\n";

          inew = inewest;
        }
    }

    virtual void write(std::basic_ostream<CHAR,TRAITS> &os, int offset=0) const{
      // fixme
    }

#ifdef PY_EXPORT

    inline TRANSFORM py_getitem(int i)
    {
      if (i<0)
        i += size();
      if (i<0 || i>=size())
        IndexError("cell: index out of range");
      return group[i];
    }

    inline void py_setitem(int i, const TRANSFORM & t)
    {
      if (i<0)
        i += size();
      if (i<0 || i>=size())
        IndexError("cell: index out of range");
      group[i] = t;
    }

    static void py_export(py::module m, const char * pyname)
    {
      py::class_<generated_group<TRANSFORM> >(m, pyname)
          .def(py::init<>())
          .def(py::init<const generated_group<TRANSFORM> &>())
          .def("index", & generated_group<TRANSFORM>::index )
          .def("add",   & generated_group<TRANSFORM>::add )
          .def("__getitem__",  & generated_group<TRANSFORM>::py_getitem)
          .def("__setitem__",  & generated_group<TRANSFORM>::py_setitem)
          .def("__len__", & generated_group<TRANSFORM>::size)
          ;
      // bp::def("len", py_group_len<TRANSFORM>);
    }

#endif
    
  };
  
};

#endif
