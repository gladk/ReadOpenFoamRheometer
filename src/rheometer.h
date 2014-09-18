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
#include <list>
#include "main.h"
#include "cell.h"

class rheometer {
  private:
    std::string _outputFolder;
    std::string _inputFolder;
    std::string _ccFolder;
    std::vector< boost::filesystem::path > _inputFolders;
    boost::filesystem::path _ccPath;
    std::vector<std::shared_ptr<cell> > _cells;
    std::vector<double> _time;
    std::map<double, unsigned int> _RhoMap, _ZMap;
    SliceMatrix _sliceMatrix;
    std::vector<std::shared_ptr<cell> > _cellVectorTmp;
    
  public:
    rheometer (const std::string & outputFolder, const std::string & inputFolder, const std::string & ccFolder);
    void loadData ();
    void prepareMatrix ();
    void loadCC (std::vector<double> & cc, const std::string & file);
    void calculateAvgU();
};
