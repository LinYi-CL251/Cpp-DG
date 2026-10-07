#pragma once
#include<stdexcept>
namespace dg{
    // Legendre basis function
inline double basis1D(double x, int degree){
    switch(degree){
        case 0:
            return 1.0;
        case 1:
            return x;
        case 2:
            return 0.5*(3.0*x*x - 1.0);
        case 3:
            return 0.5*(5.0*x*x*x - 3.0*x);
        case 4:
            return 0.125*(35.0*x*x*x*x - 30.0*x*x + 3.0);
        default:
            throw std::invalid_argument("degree must be between 0 and 4");
    }
}
inline double basis1D_derivate(double x, int degree){
    switch(degree){
        case 0:
            return 0.0;
        case 1:
            return 1.0;
        case 2:
            return 3.0*x;
        case 3:
            return 7.5*x*x - 1.5;
        case 4:
            return (140.0*x*x*x - 60.0*x) / 8.0;
        default:
            throw std::invalid_argument("degree must be between 0 and 4");
    }
}
} //namespace dg