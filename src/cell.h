/*
    This file is part of ReadOpenFoam.
    ReadOpenFoam is the programm to read OpenFoam files and create VTK.
    Copyright (C) 2014 TU Bergakademie Freiberg, Institute for Mechanics and Fluid Dynamics

    Author: 2014 Anton Gladky <gladky.anton@gmail.com>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    RheometerAnalyze is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with ReadOpenFoam.  If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once

#include <Eigen/Dense>
#include <memory>
#include <vector>

class cell {
  private:
    Eigen::Vector3d _cC;                            // Center coordinates in cartesian coordinates (r,z,fi)
    std::vector<Eigen::Vector3d> _u;                // Velocity vector
  public:
    cell (Eigen::Vector3d c) {_cC=c;};
    void addU (const Eigen::Vector3d & u) {_u.push_back(u);}
    Eigen::Vector3d c() const {return _cC;}
    const Eigen::Vector3d U(unsigned int i) const {return _u[i];}
    const Eigen::Vector3d Ucyl(unsigned int i) const {return cyl(_u[i]);}  //Returns (rho, z, phi)
    const Eigen::Vector3d Ccyl() const {return cyl(_cC);}                  //Returns (rho, z, phi)
    const Eigen::Vector3d cyl(Eigen::Vector3d c) const {                   //Returns (rho, z, phi)
      double const& x = c(0);
      double const& y = c(1);
      double const& z = c(2);
      
      double rho = sqrt(x*x + y*y);
      
      return Eigen::Vector3d(rho, z, (atan2(y,x)));
    };
};

typedef Eigen::Matrix<Eigen::Vector3d, Eigen::Dynamic, Eigen::Dynamic> SliceMatrix;
