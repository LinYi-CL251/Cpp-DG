#pragma once

#include "dg_field.hpp"
#include "equation.hpp"
#include "basis1D.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>

namespace dg {

// Zhang-Shu type scaling limiter for one-dimensional ideal-gas Euler states.
// It checks every volume quadrature point and both cell endpoints. Cell means
// are unchanged. Density is limited first, followed by a pressure scaling of
// the complete state toward its cell mean.
template<class Space>
class EulerPositivityLimiter1D {
    const Space& space_;
    Euler1D equation_;
    double epsilon_;
    std::size_t calls_ = 0;
    std::size_t limited_cells_ = 0;

    Euler1D::State state_at(const DG_Field1D& field, int cell, double xi) const {
        Euler1D::State state{};
        for (int var = 0; var < 3; ++var) {
            for (int mode = 0; mode < field.ndof(); ++mode) {
                state[var] += field(cell,var,mode)*basis1D(xi,mode);
            }
        }
        return state;
    }
    Euler1D::State state_at_quadrature(const DG_Field1D& field, int cell, int q) const {
        Euler1D::State state{};
        for (int var = 0; var < 3; ++var) {
            for (int mode = 0; mode < field.ndof(); ++mode) {
                state[var] += field(cell,var,mode)*space_.phi(q,mode);
            }
        }
        return state;
    }
    double raw_pressure(const Euler1D::State& state) const {
        return (equation_.gamma()-1.0)
            *(state[2]-0.5*state[1]*state[1]/state[0]);
    }
    template<class Function>
    void for_each_check_state(const DG_Field1D& field, int cell, Function function) const {
        for (int q = 0; q < space_.nq(); ++q) function(state_at_quadrature(field,cell,q));
        function(state_at(field,cell,-1.0));
        function(state_at(field,cell,1.0));
    }

public:
    EulerPositivityLimiter1D(const Space& space, Euler1D equation,
                             double epsilon = 1.0e-13)
        : space_(space), equation_(equation), epsilon_(epsilon) {
        if (!std::isfinite(epsilon_) || epsilon_ <= 0.0 || space_.nvar() != 3) {
            throw std::invalid_argument("Invalid Euler positivity limiter parameters");
        }
    }

    double epsilon() const { return epsilon_; }
    std::size_t calls() const { return calls_; }
    std::size_t limited_cells() const { return limited_cells_; }
    void reset_statistics() { calls_ = limited_cells_ = 0; }

    std::size_t operator()(DG_Field1D& field) {
        if (field.ncells() <= 0 || field.nvar() != 3
            || field.ndof() != space_.ndof()) {
            throw std::invalid_argument("Positivity limiter field dimensions do not match the DG space");
        }
        ++calls_;
        std::size_t changed_cells = 0;
        for (int cell = 0; cell < field.ncells(); ++cell) {
            const Euler1D::State mean = {{field(cell,0,0),field(cell,1,0),field(cell,2,0)}};
            if (!std::isfinite(mean[0]) || !std::isfinite(mean[1])
                || !std::isfinite(mean[2]) || mean[0] <= 0.0) {
                throw std::runtime_error("Euler cell average has nonpositive density");
            }
            const double mean_pressure = equation_.pressure(mean);

            double minimum_density = std::numeric_limits<double>::infinity();
            for_each_check_state(field,cell,[&](const Euler1D::State& state) {
                for (double value : state) {
                    if (!std::isfinite(value)) throw std::runtime_error("Non-finite Euler reconstruction");
                }
                minimum_density = std::min(minimum_density,state[0]);
            });
            const double density_floor = std::min(epsilon_,mean[0]);
            double theta_density = 1.0;
            if (minimum_density < density_floor) {
                theta_density = std::min(1.0,
                    (mean[0]-density_floor)/(mean[0]-minimum_density));
                // Move slightly inside the admissible set to avoid cancellation
                // when a trace lies exactly on the requested floor.
                theta_density *= 1.0-64.0*std::numeric_limits<double>::epsilon();
                for (int mode = 1; mode < field.ndof(); ++mode) {
                    field(cell,0,mode) *= theta_density;
                }
            }

            double minimum_pressure = std::numeric_limits<double>::infinity();
            for_each_check_state(field,cell,[&](const Euler1D::State& state) {
                if (!std::isfinite(state[0]) || state[0] <= 0.0) {
                    throw std::runtime_error("Density scaling failed in positivity limiter");
                }
                minimum_pressure = std::min(minimum_pressure,raw_pressure(state));
            });
            const double pressure_floor = std::min(epsilon_,mean_pressure);
            double theta_pressure = 1.0;
            if (minimum_pressure < pressure_floor) {
                theta_pressure = std::min(1.0,
                    (mean_pressure-pressure_floor)/(mean_pressure-minimum_pressure));
                theta_pressure *= 1.0-64.0*std::numeric_limits<double>::epsilon();
                for (int var = 0; var < 3; ++var) {
                    for (int mode = 1; mode < field.ndof(); ++mode) {
                        field(cell,var,mode) *= theta_pressure;
                    }
                }
            }
            // With very large kinetic/total energies, the pressure floor can be
            // smaller than a double-precision subtraction error. If necessary,
            // contract a little further; never change the mean to hide a failure.
            for (int retry = 0; retry <= 12; ++retry) {
                bool admissible = true;
                for_each_check_state(field,cell,[&](const Euler1D::State& state) {
                    admissible = admissible && std::isfinite(state[0]) && state[0] > 0.0
                        && std::isfinite(raw_pressure(state)) && raw_pressure(state) > 0.0;
                });
                if (admissible) break;
                if (retry == 12) throw std::runtime_error(
                    "Positivity scaling failed in cell " + std::to_string(cell));
                for (int var = 0; var < 3; ++var) {
                    for (int mode = 1; mode < field.ndof(); ++mode) field(cell,var,mode) *= 0.5;
                }
                theta_pressure *= 0.5;
            }
            if (theta_density < 1.0-32.0*std::numeric_limits<double>::epsilon()
                || theta_pressure < 1.0-32.0*std::numeric_limits<double>::epsilon()) {
                ++changed_cells;
            }
        }
        limited_cells_ += changed_cells;
        return changed_cells;
    }
};

template<class FirstLimiter, class SecondLimiter>
class SequentialLimiter {
    FirstLimiter& first_;
    SecondLimiter& second_;
public:
    SequentialLimiter(FirstLimiter& first, SecondLimiter& second)
        : first_(first), second_(second) {}
    void operator()(DG_Field1D& field) {
        first_(field);
        second_(field);
    }
};

} // namespace dg