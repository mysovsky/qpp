#ifndef _LACE_H
#define _LACE_H

#include <lace/calculator.hpp>
#include <lace/wizard_instance.hpp>
#include <lace/globals.hpp>
#include <lace/lapack.hpp>
//#include <lace/hlevel.hpp>

namespace lace{

	using _lace_expressions::matrix;
	using _lace_expressions::vector;
	using _lace_storage::rectang;
	using _lace_storage::symmetric;
	using _lace_storage::hermitian;  
	using _lace_storage::banddiag;   
	using _lace_storage::symmband;   
	using _lace_storage::hermband;   
	using _lace_storage::triang;
	using _lace_storage::matrix_type;
	using _lace_storage::vector_type;
	using _lace_expressions::sub;
	using _lace_storage::dense;
	using _lace_storage::sparse;
	using _lace_storage::matrix_shape;

	using _lace_storage::eigvals_valtype;

  /* ----------------------------------------------------------------
                   High level functions
     ---------------------------------------------------------------- */

  template<typename VALTYPE, matrix_type MTP>
  void solve_lu(matrix<VALTYPE,MTP> &A, vector<VALTYPE> &y, vector<VALTYPE> &x)
  {
    assert( A.size(0) == A.size(1) &&
	    A.size(0) == x.size() &&
	    A.size(0) == y.size() );

    matrix<VALTYPE,MTP> _A;
    if (globals::preserv_matrix)
      {
	_A.reshape(A.shape());
	_A.copy(A);
      }
    else
      _A.reference(A);
    // debug
//     std::cout << "_A:\n";
//     prnmtr(_A);
    _lace_storage::_solve_lu<VALTYPE,MTP>(A.shape(),_A.ptr(0,0),y.ptr(0),x.ptr(0));
  }

  template<typename VALTYPE, matrix_type MTP>
  void solve_lu(matrix<VALTYPE,MTP> &A, matrix<VALTYPE> &y, matrix<VALTYPE> &x)
  {
    assert( A.size(0) == A.size(1) &&
	    A.size(0) == x.size(0) &&
	    A.size(0) == y.size(0) &&
	    x.size(1) == y.size(1));
    
    matrix<VALTYPE,MTP> _A;
    if (globals::preserv_matrix)
      {
	_A.reshape(A.shape());
	_A.copy(A);
      }
    else
      _A.reference(A);

    // debug
//     std::cout << "_A:\n";
//     prnmtr(_A);
    _lace_storage::_solve_lu<VALTYPE,MTP>(A.shape(),x.shape(),_A.ptr(0,0),y.ptr(0,0),x.ptr(0,0));
  }


  template <typename VALTYPE, matrix_type MTRTYPE>
  void invert(matrix<VALTYPE,MTRTYPE> &A, matrix<VALTYPE,rectang> &B)
  // returns
  // B = A^(-1)
  // fixme -- replace with (s,d,c,z)getri lapack calls
  {
    int n = A.size(0);
    assert( A.size(1) == n && B.size(0) == n && B.size(1) == n && 
	    "Invalid matrix dimensions in lace::invert");
    matrix<VALTYPE,rectang> I(n);
    I = VALTYPE(0);
    for (int i=0; i<n; i++)
      I(i,i) = VALTYPE(1);
    
    solve_lu(A,I,B);
  }
  
  template <typename VALTYPE, matrix_type MTRTYPE>
  void diagon(matrix<VALTYPE,MTRTYPE> &A, 
	      vector<typename eigvals_valtype<VALTYPE,MTRTYPE>::type> &eig_vals, 
	      matrix<VALTYPE,rectang> & eig_vec)
  {
    matrix<VALTYPE,MTRTYPE> _A;
    if (globals::preserv_matrix)
      {
	_A.reshape(A.shape());
	_A.copy(A);
      }
    else
      _A.reference(A);

    _lace_storage::_diagon<VALTYPE,MTRTYPE>(_A.shape(),eig_vals.shape(),eig_vec.shape(),
					    _A.ptr(0,0),eig_vals.ptr(0),eig_vec.ptr(0,0));
  }


};

#endif