#pragma once

#include "Mesh1D.hpp"
#include "dg_mesh1D.hpp"
#include "dg_field.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>
// Suppose the semi-discrete u_t = L(u). We will compute L(u) here.
namespace dg{
    // This part of the code determines whether the source term needs to be computed.
    template<class Equation, class = void>
    struct HasSource : std::false_type {};
    template<class Equation>
    struct HasSource<Equation, std::void_t<decltype(std::declval<const Equation&()>.source(
        double{},double{},std::declval<const typename Equation::State&>()))>>: std::true_type {};
    
    enum class BoundaryType1D {Periodic, Outflow, Reflective};
    // This part of the code determines whether the reflective state needs to be considered.
    template<class Equation,class = void>
    struct HasReflectedState : std::false_type{};

    template<class Equation>
    struct HasReflectedState<Equation, std::void_t<decltype(std::declval<const Equation&>().reflected_state(
        std::declval<const typename Equation::State&>()
    ))>> : std::true_type{};

    template<class Equation, class Space>
    class DGOperator1D{
        using State = typename Equation::State;
        private:
        const Mesh1D& mesh_;
        Equation equation_;
        const Space& space_;
        BoundaryType1D boundary_;
        std::vector<double> flux_;
        std::vector<double> f_;
        std::vector<double> source_;
        std::vector<double> uL_;
        std::vector<double> uR_;

        //Trace/flux arrays: [interface][variable].
        std::size_t face_index(int face, int var) const{
            return static_cast<std::size_t>(face)*equation_.nvar() + var;
        }
        //Volume-flux array: [cell][quadrature point][variable].
        std::size_t volume_index(int cell, int q, int var) const{
            return (static_cast<std::size_t>(cell)*space_.nq() + q)*equation_.nvar() + var;
        }
        void check_field(const DG_Field1D& u) const{
            if(u.ncells()!=mesh_.ncells() || u.nvar()!=equation_.nvar()
              || u.ndof()!=space_.ndof()){
                throw std::invalid_argument("DG_Field1D dimensions do not match mesh, equation and space");
              }
        }
        State quadrature_state(const DG_Field1D& u, int cell ,int q) const{
            State value{};
            for(int var = 0; var < equation_.nvar(); ++var){
                for(int k = 0; k < space_.ndof(); ++k){
                    value[var]+= u(cell,var,k) * space_.phi(q,k);
                }
            }
            return value;
        }
        State trace_state(const DG_Field1D& u, int cell, bool right) const{
            State value{};
            for(int var = 0; var<equation_.nvar() ; ++var){
                for(int k = 0; k < space_.ndof() ; ++k){
                    value[var] += u(cell,var,k)*(right ? space_.phi_right(k) : space_.phi_left(k));
                }
            }
            return value;
        }
    };
}
