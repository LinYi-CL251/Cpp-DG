#pragma once

#include<algorithm>
#include<cmath>
#include<stdexcept>

namespace dg{
    template<class Equation>
    double LaxFriedrichsFlux(const Equation& eq, double uL, double uR){
        const double alpha = std::max(eq.max_speed(uL), eq.max_speed(uR));
        return 0.5*(eq.flux(uL) + eq.flux(uR)) - 0.5*alpha*(uR - uL);
    }

    //One spectral-radius bound for the WHOLE state, shared by all components
    template<class Equation>
    typename Equation::State LaxFriedrichsFlux(const Equation& eq, const typename Equation::State& uL,
        const typename Equation::State& uR){
            const double alpha = std::max(eq.max_speed(uL), eq.max_speed(uR));
            if(!std::isfinite(alpha) || alpha < 0.0){
                throw std::runtime_error("Invalid numerial-flux wave speed");
            }
            const auto fL = eq.flux(uL);
            const auto fR = eq.flux(uR);
            typename Equation::State result{};
            for(std::size_t var = 0; var < result.size() ; ++ var){
                result[var] = 0.5*(fL[var] + fR[var] - alpha*(uR[var] - uL[var]));
                if(!std::isfinite(result[var])){
                    throw std::runtime_error("Invalid numerical-flux value");
                }
            }
            return result;
        }
}