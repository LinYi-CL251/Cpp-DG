#pragma once

#include<vector>
#include<limits>

//The data structure of DG coefficient
// It will do not any operation to the solutions.

namespace dg{
    class DG_Field1D{
        private:
        int ncells_;
        int nvar_;
        int ndof_;

        std::vector<double> data_; // We will store the solutions in a one-dimensional vector.

        std::size_t index(int cell, int var, int mode) const{
            return(static_cast<std::size_t>(cell)*nvar_ + var)*ndof_ + mode; // one-by-one relation
        }
        public:

        DG_Field1D(int ncells, int nvar, int ndof):ncells_(ncells), nvar_(nvar), ndof_(ndof),
        data_(static_cast<std::size_t>(ncells)*static_cast<std::size_t>(nvar)*static_cast<std::size_t>(ndof),0.0){}

        double& operator()(int cell, int var, int mode){
            return data_[index(cell, var, mode)];
        }
        const double& operator()(int cell, int var, int mode) const{
            return data_[index(cell, var, mode)];
        }

        std::size_t size() const{return data_.size();}
        double* data(){return data_.data();}
        const double* data() const {return data_.data();};

        int ncells() const{return ncells_;}
        int nvar() const{return nvar_;}
        int ndof() const{return ndof_;}
    };
}