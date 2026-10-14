#pragma once
#include "characteristic.hpp"
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

    // Homogenous ideal-gas Euler equations, U = [rho, rho*u, E]
    class Euler1D{
        double gamma_;
        public:
            using State = std::array<double,3>;
            explicit Euler1D(double gamma = 1.4):gamma_(gamma){
                if (!std::isfinite(gamma_) || gamma_ <= 1.0){
                    throw std::invalid_argument("gamma must be finite and greater than 1");
                }
            }
            
            int nvar() const {return 3;}
            double gamma() const {return gamma_;}
            State conservative(double rho, double velocity, double pressure) const{
                if (!std::isfinite(rho) || !std::isfinite(velocity)
                    || !std::isfinite(pressure) || rho <= 0.0 || pressure <= 0.0) {
                    throw std::domain_error("Euler primitive state requires positive density and pressure");
                }
                return {{rho, rho * velocity,
                    pressure / (gamma_ - 1.0) + 0.5 * rho * velocity * velocity}};
            }
            double pressure(const State& u) const{
                if (!std::isfinite(u[0]) || !std::isfinite(u[1])
                    || !std::isfinite(u[2]) || u[0] <= 0.0) {
                    throw std::domain_error("Euler state has non-finite data or nonpositive density");
                }
                const double p = (gamma_ - 1.0) * (u[2] - 0.5 * u[1] * u[1] / u[0]);
                if (!std::isfinite(p) || p <= 0.0) {
                    throw std::domain_error("Euler state has nonpositive pressure");
                }
                return p;
            }

            double velocity(const State& u) const{
                return u[1] / u[0];
            }
            State flux(const State& u) const{
                const double p = pressure(u);
                const double v = u[1] / u[0];
                return {{u[1], u[1]*v + p, (u[2] + p)*v}};
            }
            double max_speed(const State& u) const{
                const double p = pressure(u);
                return std::abs(u[1]/u[0]) + std::sqrt(gamma_ * p /u[0]);
            }

            State reflected_state(const State& interior) const{
                State exterior = interior;
                exterior[1] = -exterior[1];
                return exterior;
            }
            // Freeze the Euler eigenvectors at the cell average, as in the
        // characteristic TVB limiter. Ordering: v-c, v, v+c.
        CharacteristicBasis<State> characteristic_basis(const State& mean) const {
            const double p = pressure(mean);
            const double v = mean[1] / mean[0];
            const double H = (mean[2] + p) / mean[0];
            const double c2 = gamma_ * p / mean[0];
            const double c = std::sqrt(c2);
            const double gm1 = gamma_ - 1.0;
            const double q = 0.5 * gm1 * v * v;
            CharacteristicBasis<State> basis;
            basis.right = {{{{1.0, 1.0, 1.0}},
                {{v-c, v, v+c}}, {{H-v*c, 0.5*v*v, H+v*c}}}};
            basis.left = {{{{0.5*(q+c*v)/c2, -0.5*(gm1*v+c)/c2, 0.5*gm1/c2}},
                {{(c2-q)/c2, gm1*v/c2, -gm1/c2}},
                {{0.5*(q-c*v)/c2, -0.5*(gm1*v-c)/c2, 0.5*gm1/c2}}}};
            return basis;
        }
    };
}