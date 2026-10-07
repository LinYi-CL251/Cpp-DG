#pragma once

// Class Mesh1D: Store the mesh info

namespace dg{

    class Mesh1D {
        private:
        double xmax_, xmin_; // Computational Domain [xmin_, xmax_];
        double dx_; // The Mesh size
        int ncells_;
        public:
            Mesh1D(double xmin, double xmax, int ncells):
            xmin_(xmin),xmax_(xmax),ncells_(ncells){dx_ = (xmax_ - xmin_) / ncells_;}
            double xmax()   const {return xmax_;}
            double xmin()   const {return xmin_;}
            double dx()     const {return dx_;};
            int    ncells() const {return ncells_;}
            double cell_left(int i)   const{
                return xmin_ + static_cast<double>(i) * dx_;  // i start with 0.
            }
            double cell_right(int i)  const{
                return xmin_ + (static_cast<double>(i) + 1.0) * dx_; // i start with 0.
            }
            double cell_center(int i) const{
                return xmin_ + (static_cast<double>(i) + 0.5)*dx_;  // x_mid = 0.5 * (x_left + x_right)
            }

    };
}