#pragma once
#include "Mesh1D.hpp"
#include "basis1D.hpp"
#include "quadrature.hpp"
#include "dg_field.hpp"
#include<vector>
#include<stdexcept>
// Prepare necessary info: 
// The values of basis function(quadrature points, cell boundary, derivate)
namespace dg{
    template <class Function>
    class DG_Mesh1D{
        private:
        const Mesh1D& mesh_; //Computational mesh
        int degree_; //The degree of polynomial
        int ndof_; // degree + 1
        int nvar_; // The number of equations(scalar = 1, system >=2)

        Quadrature1D quad_;
        // Interior Points
        std::vector<double> phi_; 
        std::vector<double> dphi_;
        // Cell boundary values
        std::vector<double> phi_left_; // Cell Left boundary
        std::vector<double> phi_right_; //Cell Right boundary

        std::vector<double> mass_; //We only consider the orthogonal basis. Hence, it will be a diagonal matrix.
        std::size_t index(int q, int k) const {return q*ndof_ + k;}
        public:
        DG_Mesh1D(const Mesh1D& mesh, int degree, int nq, int nvar=1):
        mesh_(mesh),degree_(degree),ndof_(degree+1),nvar_(nvar),quad_(gauss_quadrature(nq)),
        phi_(static_cast<std::size_t>(ndof_) * static_cast<std::size_t>(nq),0.0),
        dphi_(static_cast<std::size_t>(ndof_) * static_cast<std::size_t>(nq),0.0),
        phi_left_(static_cast<std::size_t>(ndof_),0.0),
        phi_right_(static_cast<std::size_t>(ndof_),0.0),
        mass_(static_cast<std::size_t>(ndof_),0.0)
        {
            if(nvar_ <=0){
                throw std::invalid_argument("Number of variables must be positive");
            }
            if(nq < degree_ + 1){
                throw std::invalid_argument("Number of quadrature points must be greater than or equal to degree + 1");
            }

            // Compute the values of basis function
            for(int q = 0; q<nq ; ++q){
                for (int k = 0; k<ndof_ ; ++k){
                    phi_[q*ndof_ + k] = basis1D(quad_.points[q],k);
                    dphi_[q*ndof_+ k] = basis1D_derivate(quad_.points[q],k);
                }
            }
            for(int k = 0; k<ndof_ ; ++k){
                phi_left_[k] = basis1D(-1.0,k);
                phi_right_[k] = basis1D(1.0,k);
            }

            //Compute the mass matrix
            for(int k = 0; k<ndof_; ++k){
                mass_[k] = 0.0;
                for(int q=0 ; q < nq; ++q){
                    mass_[k] += quad_.weights[q]*phi_[q*ndof_ + k]*phi_[q*ndof_ + k]; 
                }
                mass_[k] = 1.0/mass_[k];
            }
        }
        double phi(int q, int k) const{return phi_[index(q,k)];}
        double dphi(int q, int k) const{return dphi_[index(q,k)];}
        double phi_left(int k) const{return phi_left_[k];}
        double phi_right(int k) const{return phi_right_[k];}
        int ndof() const{return ndof_;}
        int nvar() const{return nvar_;}
        int degree() const{return degree_;}
        int nq() const{return quad_.points.size();}
        double point(int q) const{return quad_.points[q];}
        double weight(int q) const{return quad_.weights[q];}
        double inv_ref_mass(int k) const{return mass_[k];}
        void L2projection(Function f, DG_Field1D& u) const{
            if(u.ncells() != mesh_.ncells() || u.ndof()!=ndof_ || u.nvar()!=nvar_){
                throw std::invalid_argument("DGField dimensions do not match DG space");
            }
            int ncells = u.ncells();
            int nq = quad_.points.size();
            double dx = mesh_.dx();
            double x_q, x_mid, integral;
            for(int i = 0; i<ncells ; ++i){
                x_mid = 0.5*(mesh_.cell_left(i) + mesh_.cell_right(i));
                for(int eq = 0; eq < nvar_; ++eq){
                    for(int k = 0; k < ndof_; ++k){
                        integral = 0.0;
                        for(int q = 0; q<nq; ++q){
                            x_q = x_mid + 0.5*dx*quad_.points[q];
                            integral += quad_.weights[q]*phi(q,k)*f(x_q,eq);
                        }
                        u(i,eq,k) = mass_[k]*integral;
                    }
                }
            }
        }
    };
}