#ifndef _QPP_COULOMB_H
#define  _QPP_COULOMB_H

#include <geom/xgeom.hpp>
#include <symm/gen_cell.hpp>
#include <consts.hpp>
#include <Eigen/Dense>

namespace qpp{

  template<class REAL, class CELL>
  bool find_core_shells(std::vector<int> & pairs, xgeometry<REAL,CELL> & geom, REAL maxdistance){
    bool res = false;
    pairs.resize(geom.nat());
    std::fill(pairs.begin(), pairs.end(), -1);
    geom.build_types();
    int n = geom.n_types();
    std::vector<int> pair_types(n,-1);
    for (int i=0; i<n; i++)
      for (int j=0; j<i; j++){
	STRING_EX s1 = geom.atom_of_type(i);
	STRING_EX s2 = geom.atom_of_type(j);
	int k = common_begin(s1,s2);
	if (k>0 and oneof<STRING_EX>( s1.substr(k),{"cor","core"}))
	  if (oneof<STRING_EX>( s2.substr(k), {"shl","shel","shell"})){
	    pair_types[i] = j;
	    pair_types[j] = i;
	    res = true;
	  }
	if (k>0 and oneof<STRING_EX>( s2.substr(k),{"cor","core"}))
	  if (oneof<STRING_EX>( s1.substr(k), {"shl","shel","shell"})){
	    pair_types[i] = j;
	    pair_types[j] = i;
	    res = true;
	  }
      }

    /*
    std::cout << "find_cores_shells\n";
    for (int t=0; t<geom.n_types(); t++)
      std::cout << t << " " << geom.atom_of_type(t) << "\n";
    for (int p:pair_types) std::cout << p << " ";
    std::cout << "\n";
    */
    
    bonding_table<REAL> bt;
    bt.default_distance = maxdistance;
    neighbours_table<REAL,CELL> nt(geom,bt);
    nt.build();
    
    for (int i = 0; i<geom.nat(); i++)
      if (pair_types[geom.type(i)]!=-1){
	int t = pair_types[geom.type(i)];
	
	//std::cout << i << " ";
	//for (int k=0; k<nt.n(i); k++) std::cout << nt(i,k) << geom.atom(nt(i,k)(0));
	//std::cout << "\n";
	
	for (int k=0; k<nt.n(i); k++){
	  int j = nt(i,k)(0);
	  if (geom.type(j)==t){
	    pairs[i] = j;
	    pairs[j] = i;
	    res = true;
	  }
	}
      }	
    return res;
  }
  
  template<class REAL, class CELL>
  void mmcharges(std::vector<REAL> & charge, xgeometry<REAL,CELL> & geom, const std::vector<int> & regions){
    int i1 = -1, i2 = -1, ir1 = -1, ir2 = -1;
    for (int i = 0; i<geom.nfields(); i++){
      if (geom.field_name(i)=="qmm1")
	i1 = i;
      if (geom.field_name(i)=="qmm2")
	i2 = i;
      if (geom.field_name(i)=="reg1")
	ir1 = i;
      if (geom.field_name(i)=="reg2")
	ir2 = i;
    }
    if (i1==-1 || i2 == -1)
      throw std::runtime_error("\"qmm1\" and \"qmm2\" fields not found\n");
    charge.resize(geom.nat());
    for (int j = 0; j < geom.nat(); j++){
      REAL q = 0e0;
      if (oneof<int>(geom.template xfield<int>(ir1,j),regions))
	q += geom.template xfield<REAL>(i1,j);
      if (oneof<int>(geom.template xfield<int>(ir2,j),regions))
	q += geom.template xfield<REAL>(i2,j);
      charge[j] = q;
    }
  }
    
  template <class REAL>
  struct coulomb_point_charges{

    // length units in which the coordinates of charges are provided
    // energy units in which to return energy
    STRING_EX length_units, energy_units;

    REAL length_scale, energy_scale;
    
    std::vector<REAL> charges;
    std::vector<vector3<REAL> > coords;
    std::vector<int> core_shell;
    
    static REAL too_close;
    REAL core_shell_distance;
    bool shell_model;

    coulomb_point_charges()
    {
      length_units = "bohr";
      energy_units = "au";
    }

    void setscales(){
      STRING_EX lu = tolower(length_units),
	eu = tolower(energy_units);
      if (lu == "bohr")
	length_scale = 1e0;
      else if (lu == "angstrom")
	length_scale = 1e0/ang_to_bohr;
      else if (lu == "nm" || lu == "nanometer")
	length_scale = 10e0/ang_to_bohr;
      if (eu == "au" || eu == "hartree" )
	energy_scale = 1e0;
      else if (eu == "ev")
	energy_scale = hartree_to_ev;
      else if (eu == "kjm")
	energy_scale = 2625.5;
      else if (eu == "kcalm")
	energy_scale = 627.5;
    }

    template<class CELL>
    coulomb_point_charges( xgeometry<REAL,CELL> & geom, const std::vector<int> & reglist = {},
			   REAL _core_shell_distance = 0e0){
      length_units = "bohr";
      energy_units = "au";

      std::vector<int> atoms;
      if (reglist.size()==0)
	for (int i=0; i < geom.nat(); i++) atoms.push_back(i);
      else {
	std::vector<int> iregs;
	for (int i=0; i<geom.nfields(); i++)
	  if ( geom.field_name(i).substr(0,3) == "reg" && geom.field_type(i) == basic_types::type_int )
	    iregs.push_back(i);
	for (int i=0; i<geom.nat(); i++){
	  bool active = false;
	  for (int k:iregs) {
	    int r = geom.template xfield<int>(k,i);
	    if (std::find(reglist.begin(), reglist.end(), r) != reglist.end()) {
	      active = true;
	      break;
	    }
	  }
	  if (active)
	    atoms.push_back(i);
	}	  
      }

      std::vector<REAL> mmcharge;
      mmcharges(mmcharge,geom,reglist);
      
      for (int i:atoms){
	charges.push_back(mmcharge[i]);
	coords.push_back(geom.pos(i));
      }

      core_shell_distance = _core_shell_distance;
      shell_model = find_core_shells(core_shell, geom, core_shell_distance);
    }

    coulomb_point_charges(const std::vector<std::vector<REAL> > &v,
			  const std::vector<int> & _core_shell = {} ){
      length_units = "bohr";
      energy_units = "au";

      shell_model = false;
      
      for (const auto & l:v){
	charges.push_back(l[0]);
	coords.push_back({l[1],l[2],l[3]});
	core_shell.push_back(-1);
      }

      if (_core_shell.size() == coords.size()){
	shell_model = true;
	core_shell = _core_shell;
      }
      
    }
    
    REAL interaction_energy(const coulomb_point_charges<REAL> & c2)
    {
      //std::cout << "interaction_energy 1\n";
      setscales();
      REAL E = 0e0;

      bool selfenergy = (&c2 == this);

      //std::cout << "interaction_energy 2\n";
      
      for (int i = 0; i<charges.size(); i++)
	for (int j = 0; j<c2.charges.size(); j++)
	  {
	    if (selfenergy && core_shell[i] == j)
	      continue;
	    vector3<REAL> r = coords[i] - c2.coords[j];
	    REAL R = r.norm();
	    if (R > too_close)
	      E += charges[i]*c2.charges[j]/R;
	  }
      return energy_scale*E/length_scale;
    }

    void interaction_gradients(Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> & D1E,
			       const coulomb_point_charges<REAL> & c2)
    {
      setscales();
      REAL scale = energy_scale/(length_scale*length_scale);
      bool selfenergy = (&c2 == this);
      for (int i = 0; i<charges.size(); i++)
	{
	  vector3<REAL> g = 0e0;
	  for (int j = 0; j<c2.charges.size(); j++)
	    {
	      if (selfenergy && core_shell[i] == j)
		continue;
	      vector3<REAL> r = coords[i] - c2.coords[j];
	      REAL R = r.norm();
	      if (R > too_close)
		g -= charges[i]*c2.charges[j]*r/(R*R*R);
	    }
	  D1E.template block<1,3>(i,0) =  scale*g;
	}
    }

    void interaction_hessian(Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> & D2E,
			       const coulomb_point_charges<REAL> & c2)
    {
      setscales();
      REAL scale = energy_scale/(length_scale*length_scale*length_scale);
      bool selfenergy = (&c2 == this);
    }
    
#if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)

    Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic>
    py_interaction_gradients(const coulomb_point_charges<REAL> & c2){
      int n = charges.size();
      Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> D1E(n,3);
      interaction_gradients(D1E,c2);
      return D1E;
    }
    
    Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic>
    py_interaction_hessian(const coulomb_point_charges<REAL> & c2){
      int n = charges.size();
      Eigen::Matrix<REAL, Eigen::Dynamic, Eigen::Dynamic> D2E(n*3,n*3);
      interaction_hessian(D2E,c2);
      return D2E;
    }
    
    static void py_export(py::module m, const char * pyname) {
      
      py::class_<coulomb_point_charges<REAL> >(m, (std::string("coulomb_point_charges_")+pyname).c_str())
	.def(py::init<>())
	.def(py::init< xgeometry<REAL,periodic_cell<REAL> > &,
	     const std::vector<int> &, bool >(), py::arg("geom"),
	     py::arg("reglist") = std::vector<int>(), py::arg("shell_model")=0e0)
	.def(py::init< xgeometry<REAL,gen_cell<REAL,matrix3<REAL> > > &,
	     const std::vector<int> &, bool >(), py::arg("geom"),
	     py::arg("reglist") = std::vector<int>(), py::arg("shell_model")=0e0)
	.def(py::init< xgeometry<REAL,gen_cell<REAL,rotrans<REAL> > > &,
	     const std::vector<int> &, bool >(), py::arg("geom"),
	     py::arg("reglist") = std::vector<int>(), py::arg("shell_model")=0e0)
	.def(py::init<const std::vector<std::vector<REAL> > &, const std::vector<int> &>(),
	     py::arg("charges"), py::arg("core_shell") = std::vector<int>() )
	.def("interaction_energy", & coulomb_point_charges<REAL>::interaction_energy )
	.def("interaction_gradients", & coulomb_point_charges<REAL>::py_interaction_gradients )
	.def("interaction_hessian", & coulomb_point_charges<REAL>::py_interaction_hessian )
	.def("setscales", & coulomb_point_charges<REAL>::setscales )
	.def_readwrite_static("too_close", & coulomb_point_charges<REAL>::too_close)
	.def_readwrite("core_shell_distance", & coulomb_point_charges<REAL>::core_shell_distance)
	.def_readwrite("charges", & coulomb_point_charges<REAL>::charges )
	.def_readwrite("coords", & coulomb_point_charges<REAL>::coords )
	.def_readwrite("core_shell", & coulomb_point_charges<REAL>::core_shell )
	.def_readwrite("shell_model", & coulomb_point_charges<REAL>::shell_model )
	.def_readwrite("length_units", & coulomb_point_charges<REAL>::length_units )
	.def_readwrite("length_scale", & coulomb_point_charges<REAL>::length_scale )
	.def_readwrite("energy_units", & coulomb_point_charges<REAL>::energy_units )
	.def_readwrite("energy_scale", & coulomb_point_charges<REAL>::energy_scale )
	;
    }
    
#endif
    
  };

  template <class REAL>
  REAL coulomb_point_charges<REAL>::too_close = 1e-6;
  
  // template<class REAL, class CELL>
  // void coulomb_field( std::vector<REAL> & pot, std::vector<vector3<REAL> > & field, bool do_field,
  // 		      std::vector<vector3<REAL> > & where, const xgeometry<REAL,CELL> & geom)
  // {
  //   for (int i=0; i < where.size(); i++)
  //     for (int j=0; j < geom.size(); j++)
  // 	{
  // 	  vector3<REAL> Rji = where[i] - geom.pos(j);
  // 	  REAL q = geom.charge(j);
  // 	  REAL RRji = Rji.norm();
  // 	  pot[i] += q/RRji;
  // 	  if (do_field)
  // 	    field[i] += q*Rji/(RRji*RRji*RRji);
  // 	}
  // }

  // template<class REAL, class CELL>
  // REAL coulomb_interaction_energy(const xgeometry<REAL,CELL> & geom1, const xgeometry<REAL,CELL> & geom2)
  // {
  //   std::vector<REAL> R;
  //   for (int i = 0; i < geom1.size(); i++)
  //     R.push_back(geom1.pos(i));
  //   std::vector<REAL> pot(geom1.size());
  //   std::vector<vector3<REAL> > empty_field;
  //   coulomb_field(pot,empty_field,false,R,geom2);
  //   REAL E=0;
  //   for (int i = 0; i < geom1.size(); i++)
  //     E += pot[i]*geom1.charge(i);
  //   return E;
  // }
  
};

#endif
