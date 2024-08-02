#include <pyqpp/pyqpp.hpp>
#include <classic/potentials.hpp>
#include <classic/pair_potentials.hpp>
#include <classic/potentials_3b.hpp>
#include <classic/coulomb.hpp>
#include <symm/gen_cell.hpp>

namespace py = pybind11;

template<class REAL, class CELL>
std::vector<int> py_find_core_shells(qpp::xgeometry<REAL,CELL> & geom, REAL maxdistance){
  std::vector<int> res;
  if (qpp::find_core_shells(res,geom,maxdistance))
    return res;
  else
    return std::vector<int>();
}


template<class REAL, class CELL>
void py_potentials_export (py::module m, const char * pyname) {

  // generic potential
  qpp::classical_potential<REAL,CELL>::py_export(m,(std::string("classical_potential_") + pyname).c_str());
  qpp::pair_potential<REAL,CELL>::py_export(m,(std::string("pair_potential_") + pyname).c_str());
  qpp::potential_3body<REAL,CELL>::py_export(m,(std::string("potential_3body") + pyname).c_str());
  
  qpp::mm_calculator<REAL,CELL>::py_export(m,(std::string("mm_calculator_") + pyname).c_str());

  m.def("find_core_shells", &py_find_core_shells<REAL,CELL>);
  
  py::module pp = m.def_submodule("pp");
  
  qpp::buckingham_potential<REAL,CELL>::py_export(pp,(std::string("buckingham_") + pyname).c_str());
  qpp::buckingham4_potential<REAL,CELL>::py_export(pp,(std::string("buckingham4_") + pyname).c_str());
  qpp::morse_potential<REAL,CELL>::py_export(pp,(std::string("morse_") + pyname).c_str());
  qpp::cutcoulomb_potential<REAL,CELL>::py_export(pp,(std::string("cutcoulomb_") + pyname).c_str());
  qpp::spring_potential<REAL,CELL>::py_export(pp,(std::string("spring_") + pyname).c_str());
  qpp::three_harmonic<REAL,CELL>::py_export(pp,(std::string("three_harm_") + pyname).c_str());
}

void pyqpp_potentials_export (py::module m) {
  
  py_potentials_export<float, qpp::periodic_cell<float> >(m, "f");
  py_potentials_export<float, qpp::gen_cell<float, qpp::matrix3<float> > >(m, "pgf");
  py_potentials_export<float, qpp::gen_cell<float, qpp::rotrans<float> > >(m, "cgf");

  qpp::coulomb_point_charges<float>::py_export(m,"f");
   
#ifdef PYTHON_EXP_EXT
  
  py_potentials_export<double, qpp::periodic_cell<double> >(m, "d");
  py_potentials_export<double, qpp::gen_cell<double, qpp::matrix3<double> > >(m, "pgd");
  py_potentials_export<double, qpp::gen_cell<double, qpp::rotrans<double> > >(m, "cgd");

  qpp::coulomb_point_charges<double>::py_export(m,"d");
  
  
#endif
}
