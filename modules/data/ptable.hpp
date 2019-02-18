#ifndef QPP_ptable_H
#define QPP_ptable_H
#include <string>
#include <iostream>
#include <geom/lace3d.hpp>
#include <data/types.hpp>
#include <array>
#include <data/compiler_fallback.hpp>
#include <data/case_icomp.hpp>

#if defined(PY_EXPORT) || defined(QPPCAD_PY_EXPORT)
#pragma push_macro("slots")
#undef slots
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pyqpp/py_indexed_property.hpp>
namespace py = pybind11;
#pragma pop_macro("slots")
#endif

namespace qpp {

    const uint PTABLE_ELEM_N = 100;
    const uint PTABLE_INT_NONE = -1;
    const float PTABLE_FLOAT_NONE = -1.0;
    const STRING PTABLE_NONE = "X";


    struct ptable_atom_ionic_radii {
        float charge;
        float radii;
    };

    struct ptable_atom_record {

        vector3<float> m_color_jmol;
        vector3<float> m_color_gv;
        vector3<float> m_color_cpk;

        std::vector<std::tuple<int,std::string,int> > m_elec_conf;

        STRING       m_name;
        STRING       m_symbol;

        uint         m_number;
        uint         m_n_val_elec;
        float        m_radius; // A
        float        m_mass;
        float        m_covrad_slater;// A
        float        m_covrad_bragg;// A
        float        m_vdw_radius;// A
        float        m_vdw_radius_alvarez;// A
        float        m_ionic_radius;// A
        bool         m_redefined{false};

    };

    class ptable {

    public:

        static ptable *instance;

        std::map<std::string, uint, case_icomp<std::string> > cache_atom_idx;
        std::array<ptable_atom_record, PTABLE_ELEM_N> arecs;

        static ptable* get_inst () {

            if (!instance) {
                    instance = new ptable();
                    instance->init_default();
                }
            return instance;

        }

        void init_default();
        void assembly_ptable_0();
        void assembly_ptable_1();
        void assembly_ptable_2();
        void assembly_ptable_3();
        void assembly_ptable_4();
        void assembly_ptable_5();
        void assembly_ptable_6();
        void assembly_ptable_7();
        void assembly_ptable_8();
        void assembly_ptable_9();

        void init_cache_atom_idx();

        ptable () {}

        static STRING symbol_by_number (const uint number) {
            ptable *table = ptable::get_inst();
            if ((number >= 1) && (number < PTABLE_ELEM_N))
                return table->arecs[number-1].m_symbol;
            return PTABLE_NONE;
        }

        static STRING name_by_number (const int number) {

            ptable *table = ptable::get_inst();
            if (number >= 1 && number < PTABLE_ELEM_N)
                return table->arecs[number-1].m_name;
            return PTABLE_NONE;

        }

        static std::optional<uint> number_by_name (const STRING& name) {

            ptable *table = ptable::get_inst();
            for (int i = 0; i < PTABLE_ELEM_N; i++)
                if (table->arecs[i].m_name == name)
                    return std::optional<uint>(i+1);
            return std::nullopt;

        }

        static std::optional<uint> number_by_symbol (const STRING& symbol) {

            ptable *table = ptable::get_inst();
            if ( table->cache_atom_idx.find(symbol) ==
                 table->cache_atom_idx.end() ) {
              return std::nullopt;
            } else {
              return std::optional<uint>(table->cache_atom_idx[symbol]);
            }
            return std::nullopt;

        }

        static std::optional<float> mass_by_number (const uint number) {

            ptable *table = ptable::get_inst();
            if (number >= 1 && number < PTABLE_ELEM_N)
                return std::optional<float>(table->arecs[number-1].m_mass);
            else return std::nullopt;

        }

        static std::optional<float> ionic_rad_by_number (const uint number) {

            ptable *table = ptable::get_inst();
            if (number >= 1 && number < PTABLE_ELEM_N)
                return std::optional<float>(table->arecs[number-1].m_ionic_radius);
            else return std::nullopt;

        }

        static std::optional<float> cov_rad_by_number (const uint number) {

            ptable *table = ptable::get_inst();
            if (number >= 1 && number < PTABLE_ELEM_N)
                return std::optional<float>(table->arecs[number-1].m_covrad_slater);
            else return std::nullopt;

        }

        static std::optional<float> vdw_rad_by_number (const uint number) {

            ptable *table = ptable::get_inst();
            if (number >= 1 && number < PTABLE_ELEM_N)
                return std::optional<float>(table->arecs[number-1].m_vdw_radius);
            else return std::nullopt;

        }

        static uint nval_elec_by_number (const uint number) {

            ptable *table = ptable::get_inst();
            if (number >= 1 && number < PTABLE_ELEM_N)
                return table->arecs[number-1].m_n_val_elec;
            else return 1.0f;

        }

    };

}
#endif
