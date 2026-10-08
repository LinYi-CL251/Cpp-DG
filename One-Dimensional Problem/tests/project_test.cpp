#include "dg_mesh1D.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <utility>


double evaluate(const dg::DG_Field1D& u, int cell, double xi, int var = 0){
    double value = 0.0;
    for(int k = 0; k<u.ndof() ; ++k){
        value += u(cell,var,k)*dg::basis1D(xi, k);
    }
    return value;
}

template<class Function>

double measure_error(const dg::Mesh1D& mesh, const dg::DG_Field1D& u, Function exact, int var = 0){
    double l1_error = 0.0;
    double tmp = 0.0;
    double xi, x_mid;
    const double jaco = 0.5*mesh.dx();
    std::array<double,8> error_points = {{
    -0.9602898564975363, -0.7966664774136267,
    -0.5255324099163290, -0.1834346424956498,
     0.1834346424956498,  0.5255324099163290,
     0.7966664774136267,  0.9602898564975363
    }};
    std::array<double,8> error_weights = {{
    0.1012285362903763, 0.2223810344533745,
    0.3137066458778873, 0.3626837833783620,
    0.3626837833783620, 0.3137066458778873,
    0.2223810344533745, 0.1012285362903763
    }};
    for(int i = 0; i<mesh.ncells();++i){
        tmp = 0.0;
        x_mid = mesh.cell_center(i);
        for(int q = 0; q<8;++q){
            xi = x_mid + jaco*error_points[q];
            tmp+= std::abs(evaluate(u,i,error_points[q],var) - exact(xi,var))*error_weights[q];
        }
        l1_error += jaco*tmp;
    }
    return l1_error;
}
template <class Function>
std::vector<double> project_and_measure(int degree, int cells, double xmin, double xmax, Function exact, int nvar = 1){
    dg::Mesh1D mesh(xmin, xmax,cells);
    const auto initial = [&exact](double x, int var=0){
        return exact(x,var);
    };
    dg::DG_Mesh1D<decltype(initial)> space(mesh, degree,5,nvar);
    dg::DG_Field1D u(mesh.ncells(),nvar,degree + 1);
    space.L2projection(initial ,u);
    std::vector<double> error(nvar,0.0);
    for(int eq = 0; eq < nvar; eq++ ){
        error[eq] = measure_error(mesh,u,exact,eq);
    }
    return error;
}
void test_projection(){
    const double pi = std::acos(-1.0);

    const auto exact = [pi](double x, int /*var*/){return std::sin(2.0*pi*x);};
    const std::array<int,4> cell_counts = {{8,16,32,64}};
    for(int p = 0; p<=2 ; ++p){
        std::cout << "\np = " << p << '\n';
        double previous_l1 = 0.0;
        for(int cell : cell_counts){
            const auto errors = project_and_measure(p,cell,0.0,1.0,exact);
            const double error = errors[0];
            const double l1_order = previous_l1 > 0.0 && error > 0.0
                ? std::log2(previous_l1 / error) : std::numeric_limits<double>::quiet_NaN();
            std::cout << std::setw(6) << cell
                      << std::scientific << std::setprecision(6)
                      << std::setw(15) << 1.0 / cell
                      << std::setw(20) << error;
            if (std::isfinite(l1_order)) {
                std::cout << std::fixed << std::setprecision(4) << std::setw(10) << l1_order;
            } else {
                std::cout << std::setw(10) << "--";
            }
            std::cout<<'\n';
            previous_l1 = error;
        }
    }
};

void test_multi_projection(){
    const double pi = std::acos(-1.0);

    const auto exact = [pi](double x, int eq){
        if(eq == 0){
            return std::sin(2.0*pi*x);
        }else if(eq == 1){
            return x*x*x;
        }else if(eq == 2){
            return std::cos(2.0*pi*x) + 1.0;
        }else{
            throw std::invalid_argument("Not Match");
        }
    };

    const std::array<int,4> cell_counts = {{8,16,32,64}};
    for(int p = 0; p<=2 ; ++p){
        std::cout << "\np = " << p << '\n';
       
        std::vector<double> previous_l1(3,0.0);
        for(int cell : cell_counts){
            const auto errors = project_and_measure(p,cell,0.0,1.0,exact,3);
            for(int eq = 0; eq < 3; ++eq){
                const double error = errors[eq];
                const double l1_order = previous_l1[eq] > 0.0 && error > 0.0
                    ? std::log2(previous_l1[eq] / error) : std::numeric_limits<double>::quiet_NaN();
                std::cout << "Equation"<<eq<<" "<<std::setw(6) << cell
                        << std::scientific << std::setprecision(6)
                        << std::setw(15) << 1.0 / cell
                        << std::setw(20) << error;
                if (std::isfinite(l1_order)) {
                    std::cout << std::fixed << std::setprecision(4) << std::setw(10) << l1_order;
                } else {
                    std::cout << std::setw(10) << "--";
                }
                std::cout<<'\n';
                previous_l1[eq] = error;
            }
            std::cout<<'\n';
        }
    }
};
int main(){
    test_projection();
    test_multi_projection();
}
// g++ -std=c++17 -Wall -Wextra -pedantic -I"One-Dimensional Problem/include" "One-Dimensional Problem/tests/project_test.cpp" -o pt.exe
// g++ -std=c++17 -Wall -Wextra -pedantic -I"../include" "project_test.cpp" -o pt.exe