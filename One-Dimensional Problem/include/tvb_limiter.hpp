#pragma once

#include "characteristic.hpp"
#include "dg_field.hpp"
#include "Mesh1D.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace dg{
    inline double tvb_minmod(double a, double b, double c, double Mh2){
        if (!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(c)
            || !std::isfinite(Mh2) || Mh2 < 0.0) {
            throw std::invalid_argument("TVB minmod requires finite data and a nonnegative threshold");
        }
        if (std::abs(a) <= Mh2) return a;
        if (a > 0.0 && b > 0.0 && c > 0.0) return std::min({a,b,c});
        if (a < 0.0 && b < 0.0 && c < 0.0) return -std::min({-a,-b,-c});
        return 0.0;
    }

    enum class LimitingVariables{ Characteristic, Conservative};
    //Judeg whether the system has characteristic basis
    template<class Equation, class = void>
    struct HasCharacteristicBasis : std::false_type {};
    template<class Equation>
    struct HasCharacteristicBasis<Equation, std::void_t<decltype(
    std::declval<const Equation&>().characteristic_basis(
        std::declval<const typename Equation::State&>()))>> : std::true_type {};

    // A new system supplies characteristic_basis(mean), returning L and R.
    // Scalar equations need no eigenvector interface. Conservative limiting
    // must be explicitly selected for a system that has no characteristic basis.

    template<class Equation>
    class EquationCharacteristics{
        Equation equation_;
        LimitingVariables variables_;
        public:
        using State = typename Equation::State;
        explicit EquationCharacteristics(Equation equation, 
            LimitingVariables variables = LimitingVariables::Characteristic)
            : equation_(equation), variables_(variables) {
                if constexpr (std::tuple_size<State>::value > 1 && !HasCharacteristicBasis<Equation>::value) {
                    if (variables_ == LimitingVariables::Characteristic) {
                        throw std::invalid_argument("System characteristic limiting requires characteristic_basis(mean)");
                    }
                }
        }

        CharacteristicBasis<State> operator()(const State& mean) const{
            if constexpr (HasCharacteristicBasis<Equation>::value) {
                if(variables_ == LimitingVariables::Characteristic){
                    return equation_.characteristic_basis(mean);
                }
            }
            return CharacteristicBasis<State>::identity();
        }
    };

    struct PeriodicMeanNeighbors {
        template<class State>

        State operator()(const std::vector<State>& means, int cell, int offset) const{
            const int n = static_cast<int>(means.size());
            return means[(cell + offset + n) % n];
        }
    };

    struct OutflowMeanNeighbors{
        template<class State>
        State operator()(const std::vector<State>& means, int cell, int offset) const{
            const int neighbor = cell + offset;
            if(neighbor < 0) return means.front();
            if(neighbor >= static_cast<int>(means.size())) return means.back();
            return means[neighbor];
        }
    };

    struct EulerReflectingMeanNeighbors{
        std::array<double,3> operator()(const std::vector<std::array<double,3>>& means,
                                   int cell, int offset) const {
        const int neighbor = cell + offset;
        if (neighbor >= 0 && neighbor < static_cast<int>(means.size())) return means[neighbor];
        auto ghost = means[cell];
        ghost[1] = -ghost[1];
        return ghost;
    }
    };

    // Legendre modal TVB endpoint limiter. Degree 0 is unchanged. Unmarked cells
    // retain every mode. In marked cells, match the two limited endpoint deviations
    // with modes 1 and 2, then remove modes >=3. Mode 0 (the mean) is never modified.
    template<class Equation, class Transform = EquationCharacteristics<Equation>, 
            class MeanNeighbors = PeriodicMeanNeighbors>
    class TVBLimiter1D{
        using State = typename Equation::State;
        const Mesh1D& mesh_;
        int nvar_;
        double M_;
        double threshold_;
        Transform transform_;
        MeanNeighbors neighbors_;
        std::vector<State> means_;
        std::size_t calls_ = 0;
        std::size_t limited_cells_ = 0;
        public:
        // Constructors
        TVBLimiter1D(const Mesh1D& mesh, Equation equation, double M,
                 Transform transform, MeanNeighbors neighbors = {})
        : mesh_(mesh), nvar_(equation.nvar()), M_(M), threshold_(M*mesh.dx()*mesh.dx()),
          transform_(std::move(transform)), neighbors_(std::move(neighbors)),
          means_(mesh.ncells()) {
            if (!std::isfinite(M_) || M_ < 0.0 || mesh.ncells() <= 0
                || !std::isfinite(mesh.dx()) || mesh.dx() <= 0.0
                || !std::isfinite(threshold_) || nvar_ != static_cast<int>(State{}.size())) {
                throw std::invalid_argument("Invalid TVB limiter mesh, M or equation dimensions");
            }
        }
        // Copy constructor
        TVBLimiter1D(const Mesh1D& mesh, Equation equation, double M,
                    LimitingVariables variables = LimitingVariables::Characteristic)
            : TVBLimiter1D(mesh, equation, M, Transform(equation, variables)) {}

        double M() const { return M_; }
        std::size_t calls() const { return calls_; }
        std::size_t limited_cells() const { return limited_cells_; }
        void reset_statistics() { calls_ = limited_cells_ = 0; }

        std::size_t operator()(DG_Field1D& field) {
            if (field.ncells() != mesh_.ncells() || field.nvar() != nvar_ || field.ndof() < 1) {
                throw std::invalid_argument("TVB limiter field dimensions do not match mesh/equation");
            }
            ++calls_;
            if (field.ndof() == 1) return 0;
            for (int cell = 0; cell < field.ncells(); ++cell) {
                for (int v = 0; v < nvar_; ++v) {
                    means_[cell][v] = field(cell,v,0);
                    if (!std::isfinite(means_[cell][v])) throw std::runtime_error("Non-finite cell average");
                }
            }
            std::size_t changed_cells = 0;
            for (int cell = 0; cell < field.ncells(); ++cell) {
                const State left_mean = neighbors_(means_,cell,-1);
                const State right_mean = neighbors_(means_,cell,1);
                State dr{}, dl{}, forward{}, backward{};
                for (int v = 0; v < nvar_; ++v) {
                    forward[v] = right_mean[v] - means_[cell][v]; // \bar{u}_{j+1} - \bar{u}_j
                    backward[v] = means_[cell][v] - left_mean[v]; // \bar{u}_j - \bar{u}_{j-1}
                    for (int k = 1; k < field.ndof(); ++k) {
                        dr[v] += field(cell,v,k); // P_k(+1)=1. // u_{j+\frac12} - \bar{u}_j
                        dl[v] += (k % 2 ? 1.0 : -1.0)*field(cell,v,k); // mean-u(-1). \bar{u}_j - u_{j - \frac12}
                    }
                }
                const auto basis = transform_(means_[cell]);
                // Characteristic space
                const State cr = matrix_vector<State>(basis.left,dr);
                const State cl = matrix_vector<State>(basis.left,dl);
                const State cf = matrix_vector<State>(basis.left,forward);
                const State cb = matrix_vector<State>(basis.left,backward);
                State nr{}, nl{};
                bool changed = false;
                for (int v = 0; v < nvar_; ++v) {
                    nr[v] = tvb_minmod(cr[v],cf[v],cb[v],threshold_);
                    nl[v] = tvb_minmod(cl[v],cf[v],cb[v],threshold_);
                    // Compare in characteristic space; do not let L/R roundoff
                    // alone trigger a rebuild of an otherwise smooth polynomial.
                    const double tolerance = 64.0*std::numeric_limits<double>::epsilon()
                        * std::max({1.0,std::abs(cr[v]),std::abs(cl[v]),std::abs(cf[v]),std::abs(cb[v])});
                    changed = changed || std::abs(nr[v]-cr[v]) > tolerance
                                    || std::abs(nl[v]-cl[v]) > tolerance;
                }
                if (!changed) continue;
                const State limited_right = matrix_vector<State>(basis.right,nr);
                const State limited_left = matrix_vector<State>(basis.right,nl);
                for (int v = 0; v < nvar_; ++v) {
                    field(cell,v,1) = 0.5*(limited_right[v]+limited_left[v]);
                    if (field.ndof() >= 3) field(cell,v,2) = 0.5*(limited_right[v]-limited_left[v]);
                    for (int k = 3; k < field.ndof(); ++k) field(cell,v,k) = 0.0;
                }
                ++changed_cells;
            }
            limited_cells_ += changed_cells;
            return changed_cells;
        }
    };
}