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

#include "rheometer.h"

bool sortFileTimeCreate(boost::filesystem::path i, boost::filesystem::path j) {
  return (boost::filesystem::last_write_time(i) < boost::filesystem::last_write_time(j));
}

rheometer::rheometer (const std::string & outputFolder, const std::string & inputFolder, const std::string & ccFolder) {
  
  namespace fs = boost::filesystem;
  _outputFolder = outputFolder;
  _inputFolder = inputFolder;
  _ccFolder = ccFolder;
  
  
  std::vector< fs::path > inputFolders;
  fs::path inpDir (inputFolder);
  fs::directory_iterator end_iter;
  
  if ( fs::exists(inpDir) && fs::is_directory(inpDir))
  {
    for( fs::directory_iterator dir_iter(inpDir) ; dir_iter != end_iter ; ++dir_iter)
    {
      if (fs::is_directory(dir_iter->status()) )
      {
        fs::path curDir (*dir_iter);
        fs::path curFileVel = curDir.string() + "/U";
        
        if(fs::is_regular_file(curFileVel)) {
          inputFolders.push_back(curDir);
        }
      }
    }
  }
  
  std::sort(inputFolders.begin(), inputFolders.end(), sortFileTimeCreate);
  
  std::cout<<inputFolders.size()<<" file(s) found to analyze"<<std::endl;
  
  fs::path ccPath (ccFolder);
  
  if ( fs::exists(ccPath) && fs::is_directory(ccPath)) {
    fs::path ccx = ccPath.string() + "/ccx";
    fs::path ccy = ccPath.string() + "/ccy";
    fs::path ccz = ccPath.string() + "/ccz";
    if(fs::is_regular_file(ccx) and fs::is_regular_file(ccy) and fs::is_regular_file(ccz)) {
        std::cout<<"ccx, ccy and ccz files are found"<<std::endl;
    } else {
      std::cerr<<"Please specify the folder where ccx, ccy and ccz files are available!"<<std::endl;
      exit (EXIT_FAILURE);
    }
  } else {
    std::cerr<<"Please specify the folder where ccx, ccy and ccz files are available!"<<std::endl;
    exit (EXIT_FAILURE);
  }

  //=====================================================
  if (not fs::is_directory(outputFolder)) {
    std::cout<<"The directory " << outputFolder<< " does not exists. Creating."<<std::endl;
    if (fs::create_directory(outputFolder)) {
      std::cout<<"The directory " << outputFolder<< " created."<<std::endl;
    }
  }
  
}

