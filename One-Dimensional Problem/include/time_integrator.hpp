#pragma once
#include "dg_field.hpp"
#include <algorithm>
#include <cstddef>
#include <cmath>
#include <stdexcept>
#include <type_traits>
namespace dg{

struct NoLimiter{
    void operator()(DG_Field1D&) const{}
};
// Reusable wrok fields; no allocation in time loop.
class SSPRK3 {
    DG_Field1D old_;
    DG_Field1D stage_;
    DG_Field1D rhs_;
    template<class Limiter>
    static void apply_limiter(Limiter& limiter, DG_Field1D& field, double time){
        if constexpr (std::is_invocable_v<Limiter&,DG_Field1D&,double>) limiter(field,time);
        else limiter(field);
    }
    public:

    explicit SSPRK3(const DG_Field1D& prototype):old_(prototype.ncells(),prototype.nvar(),prototype.ndof()),
        stage_(prototype.ncells(),prototype.nvar(),prototype.ndof()),
        rhs_(prototype.ncells(),prototype.nvar(),prototype.ndof()){}
    
    template<class Operator>
    void step(DG_Field1D& u, double t, double dt, Operator& op) {
        NoLimiter limiter;
        step(u,t,dt,op,limiter);
    }

    template<class Operator, class Limiter>

    void step(DG_Field1D& u, double t, double dt, Operator& op, Limiter& limiter){
        if(!std::isfinite(t) || !std::isfinite(dt) || dt <= 0.0){
            throw std::invalid_argument("SSPRK3 requires finite time and positive dt");
        }
        if(u.ncells() != old_.ncells() || u.nvar() != old_.nvar() 
            || u.ndof() != old_.ndof()){
                throw std::invalid_argument("SSPRK3 field dimensions changed");
        }
        std::copy(u.data(), u.data() + u.size(), old_.data());

        op.rhs(old_,t,rhs_);
        const std::ptrdiff_t coefficient_count = static_cast<std::ptrdiff_t>(u.size());
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static) if(coefficient_count>=1048576)
        #endif
        for(std::ptrdiff_t j = 0; j < coefficient_count; ++j){
            stage_.data()[j] = old_.data()[j] + dt * rhs_.data()[j];
        }
        apply_limiter(limiter,stage_,t+dt);
        op.rhs(stage_, t+dt, rhs_);
        #ifdef _OPENMP
        #pragma omp parallel for schedule(static) if(coefficient_count>=1048576)
        #endif
        for (std::ptrdiff_t j = 0; j < coefficient_count; ++j) {
            // Equivalent to 3/4 old + 1/4 (stage + dt*rhs).
            stage_.data()[j] = old_.data()[j]
                + 0.25 * (stage_.data()[j] - old_.data()[j] + dt * rhs_.data()[j]);
        }
        apply_limiter(limiter,stage_,t+0.5*dt);

        op.rhs(stage_, t + 0.5 * dt, rhs_);
        int invalid_coefficient=0;
        #ifdef _OPENMP
        #pragma omp parallel for reduction(|:invalid_coefficient) schedule(static) if(coefficient_count>=1048576)
        #endif
        for (std::ptrdiff_t j = 0; j < coefficient_count; ++j) {
            // Equivalent to 1/3 old + 2/3 (stage + dt*rhs).
            u.data()[j] = old_.data()[j]
                + (2.0 / 3.0) * (stage_.data()[j] - old_.data()[j] + dt * rhs_.data()[j]);
            invalid_coefficient |= !std::isfinite(u.data()[j]);
        }
        if (invalid_coefficient) throw std::runtime_error("SSPRK3 produced a non-finite coefficient");
        apply_limiter(limiter,u,t+dt);
    }
};
}