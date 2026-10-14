#pragma once

#include<cmath>
#include<array>
#include<functional>
#include<stdexcept>

namespace dg{
    
    class LinearAdvection{
        private:
            double speed_;
        public:
            using State = std::array<double,1>;
            int nvar() const {return 1;}
            LinearAdvection(double speed) : speed_(speed){}
            double flux(double u) const{return speed_*u;}
            double max_speed(double /*u*/) const{
                return std::abs(speed_);
            }
            State flux(const State& u) const{return {{flux(u[0])}};}
            double max_speed(const State& u) const{return max_speed(u[0]);}
    };

    class Burgers{
        public:
            using State = std::array<double,1>;
            int nvar() const {return 1;}
            double flux(double u) const{
                return 0.5*u*u;
            }
            double max_speed(double u) const{
                return std::abs(u);
            }
            State flux(const State& u) const {return {{flux(u[0])}};}

            double max_speed(const State& u) const{return max_speed(u[0]);}
    };
}