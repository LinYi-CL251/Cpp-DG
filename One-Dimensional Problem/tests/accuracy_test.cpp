// Accuracy Test for the One-Dimensional Problem

//Setup: 
// (i) LinearAdvection: domain [0,1], periodic boundary conditions, u0 = sin(2*pi*x), T = 0.5
//                      Exact solution: u(x,t) = sin(2*pi*(x - 2*t))
// (ii) Burgers: domain [0,1], periodic boundary conditions, u0 = sin(2*pi*x), T = 0.1
//            Exact solution: Using Newton iteration. (before shock at t_b = 1/(2*pi))
//Time stepping: SSPRK3 with dt = cfl*dx/a_max.
// g++ -O3 -I"../include" accuracy_test.cpp -o at.exe
#include "dg_mesh1D.hpp"
#include "dg_operator.hpp"
#include "time_integrator.hpp"
#include "equation.hpp"

#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>
namespace {
    const double PI = std::acos(-1.0);

    double evaluate(const dg::DG_Field1D& u, int cell, double xi, int var = 0){
        double value = 0.0;
        for (int k = 0; k< u.ndof(); ++k){
            value += u(cell,var,k) * dg::basis1D(xi, k);
        }
        return value;
    }

//8-points Gauss-Legendre quadrature
const std::array<double, 8> kErrPts = {{
    -0.9602898564975363, -0.7966664774136267,
    -0.5255324099163290, -0.1834346424956498,
     0.1834346424956498,  0.5255324099163290,
     0.7966664774136267,  0.9602898564975363
}};
const std::array<double, 8> kErrWts = {{
    0.1012285362903763, 0.2223810344533745,
    0.3137066458778873, 0.3626837833783620,
    0.3626837833783620, 0.3137066458778873,
    0.2223810344533745, 0.1012285362903763
}};

template <class Exact>
double measure_L1(const dg::Mesh1D& mesh, const dg::DG_Field1D& u, Exact exact, int var = 0){
    double l1 = 0.0;
    const double jaco = 0.5 * mesh.dx();
    for (int i = 0; i<mesh.ncells(); ++i){
        const double x_mid = mesh.cell_center(i);
        double tmp = 0.0;
        for (int q = 0 ; q < 8 ; ++q){
            const double xi = x_mid + jaco * kErrPts[q];
            tmp += std::abs(evaluate(u,i,kErrPts[q],var) - exact(xi,var)) * kErrWts[q];
        }
        l1 += jaco*tmp;
    }
    return l1;
}

template <class Exact>
double measure_Linf(const dg::Mesh1D& mesh, const dg::DG_Field1D& u, Exact exact, int var = 0){
    double linf = 0.0;
    const double jaco = 0.5 * mesh.dx();
    for (int i = 0; i < mesh.ncells(); ++i) {
        const double x_mid = mesh.cell_center(i);
        for (int q = 0; q < 8; ++q) {
            const double xi = x_mid + jaco * kErrPts[q];
            linf = std::max(linf, std::abs(evaluate(u, i, kErrPts[q], var) - exact(xi, var)));
        }
    }
    return linf;
}

template <class Equation, class Exact>
void run_convergence(const std::string& name, Equation eq, Exact exact,
     double T, const std::vector<int>& cells_list, const std::vector<int>& degrees,
     double x_min, double x_max, double cfl){
        for(int p : degrees){
            std::cout << "\n=== " << name << ", p = " << p << " ===\n";
            std::cout << std::setw(8) << "N" << std::setw(20) << "L1 error"
                    << std::setw(12) << "L1 order" << std::setw(20) << "Linf error"
                    << std::setw(12) << "Linf order" << '\n';
            double prev_l1 = 0.0, prev_linf = 0.0;
            for(int N : cells_list){
                dg::Mesh1D mesh(x_min,x_max,N);
                
                auto init = [&](double x, int /*var*/){return exact(x,0.0);};
                dg::DG_Mesh1D<decltype(init)> space(mesh, p, 6, 1);
                dg::DG_Field1D u(N, 1, p+1);
                space.L2projection(init, u);
                dg::DGOperator1D<Equation, decltype(space)> op(mesh,eq,space,dg::BoundaryType1D::Periodic);
                dg::SSPRK3 rk(u);
                double t = 0.0;
                while (t < T) {
                    const double smax = op.max_speed(u);
                    double dt = cfl * mesh.dx() / std::max(smax, 1e-12);
                    if (t + dt > T) dt = T - t;
                    rk.step(u, t, dt, op);
                    t += dt;
                }
                auto exactT = [&](double x, int /*var*/) {return exact(x,T);};
                const double l1 = measure_L1(mesh, u, exactT);
                const double linf = measure_Linf(mesh, u, exactT);
                const double o1 = (prev_l1 > 0.0 && l1 > 0.0) ? std::log2(prev_l1 / l1)
                                                           : std::numeric_limits<double>::quiet_NaN();
                const double oinf = (prev_linf > 0.0 && linf > 0.0) ? std::log2(prev_linf / linf)
                                                                    : std::numeric_limits<double>::quiet_NaN();
                std::cout << std::setw(8) << N << std::scientific << std::setprecision(6)
                        << std::setw(20) << l1 << std::fixed << std::setprecision(3)
                        << std::setw(12) << o1 << std::scientific << std::setprecision(6)
                        << std::setw(20) << linf << std::fixed << std::setprecision(3)
                        << std::setw(12) << oinf << '\n';
                prev_l1 = l1;
                prev_linf = linf;
            }
        }
     }
}

int main(){
    const std::vector<int> cells = {10, 20, 40 ,80, 160};
    const std::vector<int> degrees = {0,1,2};

    auto exact_adv = [](double x, double t) {return std::sin(2.0*PI*(x-t));};
    run_convergence("Linead advection(a = 1), T = 0.5", dg::LinearAdvection(1.0), exact_adv,0.5,cells,degrees,0.0,1.0,0.16);

    auto exact_burg = [](double x, double t) {
        double x0 = x;
        for (int it = 0; it < 100; ++it) {
            const double f = x0 + t * std::sin(2.0 * PI * x0) - x;
            const double df = 1.0 + t * 2.0 * PI * std::cos(2.0 * PI * x0);
            const double dx = f / df;
            x0 -= dx;
            if (std::abs(dx) < 1e-14) break;
        }
        return std::sin(2.0 * PI * x0);
    };
    run_convergence("Burgers, T=0.1", dg::Burgers(), exact_burg, 0.1, cells, degrees,0.0,1.0,0.16);
    return 0;
}