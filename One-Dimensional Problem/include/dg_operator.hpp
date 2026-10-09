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
            // |____|____|____| right: ->; left: <-
            State value{};
            for(int var = 0; var<equation_.nvar() ; ++var){
                for(int k = 0; k < space_.ndof() ; ++k){
                    value[var] += u(cell,var,k)*(right ? space_.phi_right(k) : space_.phi_left(k));
                }
            }
            return value;
        }
        double check_speed(const State& state) const{
            for(double value : state){
                if(!std::isfinite(value)){
                    throw std::runtime_error("Non-finite reconstructed DG state");
                }
                const double speed = equation_.max_speed(state);
                if(!std::isfinite(speed) || speed < 0.0){
                    throw std::runtime_error("Invalid characteristic wave speed");
                }
            }
            return speed;
        }
        State reflected_state(const State& interior) const{
            if constexpr(HasReflectedState<Equation>::value){
                return equation_.reflected_state(interior);
            }else{
                throw std::invalid_argument("Reflective boundary requires Equation::reflected_state(interior)");
            }
        }
        public:
        DGOperator1D(const Mesh1D& mesh, Equation equation, const Space& space,
                    BoundaryType1D boundary = BoundaryType1D::Periodic) : mesh_(mesh),
                    equation_(equation),space_(space), boundary_(boundary),
                    flux_((static_cast<std::size_t>(mesh.ncells()) + 1)*(equation.nvar()),0.0),
                    f_(static_cast<std::size_t>(mesh.ncells()) * space.nq() * equation.nvar(),0.0),
                    source_(f_.size(),0.0),uL_(flux_.size(),0.0),uR_(flux_.size(),0.0){
                        if (space.nvar()!=equation.nvar() || static_cast<int>(State{}.size()) != equation_.nvar()){
                            throw std::invalid_argument("Equation and DG space must have the same number of variables");
                        } 
                        if (boundary_ != BoundaryType1D::Periodic && boundary_ != BoundaryType1D::Outflow
                            && boundary_ != BoundaryType1D::Reflective) {
                            throw std::invalid_argument("Unknown one-dimensional boundary type");
                        }
                        if constexpr (!HasReflectedState<Equation>::value) {
                            if (boundary_ == BoundaryType1D::Reflective) {
                                throw std::invalid_argument("This equation has no reflected_state interface");
                            }
                        }
        }

        BoundaryType1D boundary_type() const {return boundary_;}
        
        //Net inward conservative flux from the most recent rhs evalution.

        State net_boundary_flux() const{
            State result{};
            for(int var = 0 ; var < equation_.nvar(); ++var){
                result[var] = flux_[face_index(0,var)] - flux_[face_index(mesh_.ncells(),var)];
            }
            return result;
        }

        // U_t = L(U)
        void rhs(const DG_Field1D& u, double t, DG_Field1D& results) {
            check_field(u);
            check_field(results);
            if(&u == &results) {
                throw std::invalid_argument("rhs input and output must be distince fields");
            }
            assemble_flux(u,t,flux_);
            assemble_volume_flux(u,f_);
        }

        void assemble_flux(const DG_Field1D& u, double /*t*/, std::vector<double>& FLUX){
            FLUX.resize(flux_.size());
            std::fill(uL_.begin(), uL_.end(),0.0);
            std::fill(uR_.begin(), uR_.end(),0.0);
            const int ncell = mesh_.ncells();
            const int nvar = equation_.nvar();

            // Interior faces: cell face - 1 supplies the left trace and cell face
            // supplies the right face.
            for(int face = 1; face < ncell ; ++face){
                const State left = trace_state(u, face - 1, true);
                const State right = trace_state(u, face, false);
                for (int var = 0; var < nvar; ++var){
                    uL_[face_index(face,var)] = left[var];
                    uR_[face_index(face,var)] = right[var];
                }
            }

            // Using boundary conditions

            const State first_left = trace_state(u,0,false); //The first cell, left boundary
            const State last_right = trace_state(u,ncell - 1, true); //The last cell, right boundary
            State exterior_left{};
            State exterior_right{};
            if(boundary == BoundaryType1D::Periodic){ // Periodic Boundary conditions
                exterior_left = last_right;
                exterior_right = first_left;
            }else if(boundary == BoundaryType1D::Outflow){ // Outflow Boundary conditions
                exterior_left = first_left;
                exterior_right = last_right;
            }else{                                         // Reflective Boundary conditions
                exterior_left = reflected_state(first_left);
                exterior_right = reflected_state(last_right);
            }

            for(int var = 0; var < nvar ; ++var){
                uL_[face_index(0,var)] = exterior_left[var];
                uR_[face_index(0,var)] = first_left[var];
                uL_[face_index(ncell,var)] = last_right[var];
                uR_[face_index(ncell,var)] = exterior_right[var];
            }

            // using uL and uR to compute numerical fluxes

            const int last_flux_face = boundary_ ==BoundaryType1D::Periodic
                    ? ncell - 1 : ncell;
            for(int face = 0; face <= last_flux_face; ++face){
                State left{};
                State right{};
                for(int var = 0; var < nvar; ++var){
                    left[var] = uL_[face_index(face, var)];
                    right[var] = uR_[face_index(face,var)];
                }
                const State numerical_flux = LaxFriedrichsFlux(equation_,left,right);
                for(int var = 0; var < nvar; ++var){
                    FLUX[face_index(face,var)] = numerical_flux[var];
                }
            }

            if(boundary_ == BoundaryType1D::Periodic){
                for(int var = 0; var<nvar; ++var){
                    FLUX[face_index(ncell,var)] = FLUX[face_index(0,var)];
                }
            }
        }

        void assemble_volume_flux(const DG_Field1D& u, std::vector<double>& F){
            F.resize(f_.size());
            for (int cell = 0; cell < mesh_.nells(); ++cell){
                for(int q = 0; q < space_.nq() ; ++q){
                    const State physical_flux = equation_.flux(quadrature_state(u,cell,q));
                    for(int var = 0; var < equation_.nvar(); ++var){
                        if(!std::isfinite(physical_flux[var])){
                            throw std::runtime_error("Non-finite volume flux");
                        }
                        F[volume_index(cell,q,var)] = physical_flux[var];
                    }
                }
            }
        }
    };
}
