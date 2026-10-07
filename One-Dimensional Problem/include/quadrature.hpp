#pragma once

#include<vector>
#include<stdexcept>
// Gauss-Legrendre Points
namespace dg{
    struct Quadrature1D
    {
        std::vector<double> points;
        std::vector<double> weights;
    };
    inline Quadrature1D gauss_quadrature(int n){
        Quadrature1D quad;
        if(n <= 2){
            throw std::invalid_argument("Number of points must be greater than 2");
        } 
        else if(n == 3){
            const double a = 0.7745966692414834;
            quad.points = {-a, 0.0, a};
            quad.weights = {5.0/9.0, 8.0/9.0, 5.0/9.0};
        }else if(n == 4){
            const double a = 0.8611363115940526;
            const double b = 0.3399810435848563;
            quad.points = {-a, -b, b, a};
            quad.weights = {0.3478548451374538, 0.6521451548625461, 
                0.6521451548625461, 0.3478548451374538};
        }else if(n == 5){
            const double a = 0.9061798459386640;
            const double b = 0.5384693101056831;
            quad.points = {-a, -b, 0.0, b, a};
            quad.weights = {0.2369268850561891, 0.4786286704993665, 
                0.5688888888888889, 0.4786286704993665, 0.2369268850561891};
        }else if(n == 6){
            const double a = 0.9324695142031521;
            const double b = 0.6612093864662645;
            const double c = 0.2386191860831969;
            quad.points = {-a, -b, -c, c, b, a};
            quad.weights = {0.1713244923791704, 0.3607615730481386,
                0.4679139345726910, 0.4679139345726910,
                0.3607615730481386, 0.1713244923791704};
        }else{
            throw std::invalid_argument("Number of points must be 3, 4, 5 or 6");
        }
        return quad;
    }
}