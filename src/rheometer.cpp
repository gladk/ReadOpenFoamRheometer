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
  boost::filesystem::path ii = i.string() + "/uniform/time";
  boost::filesystem::path jj = j.string() + "/uniform/time";
  return (boost::filesystem::last_write_time(ii) < boost::filesystem::last_write_time(jj));
}

rheometer::rheometer (const std::string & outputFolder, const std::string & inputFolder, const std::string & ccFolder) {
  namespace fs = boost::filesystem;
  
  _outputFolder = outputFolder;
  _inputFolder = inputFolder;
  _ccFolder = ccFolder;
  
  
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
        fs::path curFileTime= curDir.string() + "/uniform/time";
        
        if(fs::is_regular_file(curFileVel) && fs::is_regular_file(curFileTime)) {
          _inputFolders.push_back(curDir);
        }
      }
    }
  }
  
  std::sort(_inputFolders.begin(), _inputFolders.end(), sortFileTimeCreate);
  
  std::cout<<_inputFolders.size()<<" file(s) found to analyze"<<std::endl;
  
  _ccPath = ccFolder;
  
  if ( fs::exists(_ccPath) && fs::is_directory(_ccPath)) {
    fs::path ccx = _ccPath.string() + "/ccx";
    fs::path ccy = _ccPath.string() + "/ccy";
    fs::path ccz = _ccPath.string() + "/ccz";
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
  
  loadData();
}

void rheometer::loadData () {
  namespace fs = boost::filesystem;
  
  //Load coordinates
  std::vector<double> ccx; loadCC(ccx, "ccx");
  std::vector<double> ccy; loadCC(ccy, "ccy");
  std::vector<double> ccz; loadCC(ccz, "ccz");
  
  
  if (ccx.size()==ccy.size() and ccx.size()==ccz.size()) {
    std::cout<<ccx.size()<<" coordinates loaded"<<std::endl;
    for (unsigned int i=0; i<ccx.size(); i++)  {
      _cells.push_back(Eigen::Vector3d(ccx[i],ccy[i],ccz[i]));
    }
  } else {
    std::cerr << "error: ccx!=ccy!=ccz" << std::endl;
    exit (EXIT_FAILURE);
  }
  
  BOOST_FOREACH(fs::path f, _inputFolders) {
    // Time  =======================================================
    fs::path curFileTime = f.string() + "/uniform/time";
    fs::path curFileU = f.string() + "/U";
    if (not (fs::exists(curFileU)) or not (fs::exists(curFileTime))) {continue;}
    
    std::ifstream fileT(curFileTime.string(), std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_istream inT;
    inT.push(fileT);
    unsigned int i = 0;
    for(std::string str; std::getline(inT, str); )
    {
      if (i == 17) {
        std::vector<std::string> tokens;
        std::istringstream iss(str);
        copy(std::istream_iterator<std::string>(iss),
             std::istream_iterator<std::string>(),
             std::back_inserter<std::vector<std::string> >(tokens));
        _time.push_back(stod (tokens[1]));
        break;
      }
      i++;
    }
    std::cout<<"Time "<<_time[_time.size()-1]<<std::endl;
    // U  =======================================================
    
    std::ifstream fileU(curFileU.string(), std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_istream inU;
    inU.push(fileU);
    i = 0;
    const unsigned int ccSize = _cells.size();
    for(std::string str; std::getline(inU, str); )
    {
      if (i > 21 and i < (22 + ccSize)) {
        std::vector<std::string> tokens;
        std::istringstream iss(str);
        copy(std::istream_iterator<std::string>(iss),
             std::istream_iterator<std::string>(),
             std::back_inserter<std::vector<std::string> >(tokens));
        Eigen::Vector3d U(stod(tokens[0].erase(0,1)), stod(tokens[1]), stod(tokens[2].erase(tokens[2].size()-1,1)));
        _cells[i-22].addU(U);
      }
      i++;
    }
  }
}

void rheometer::loadCC (std::vector<double> & cc, const std::string & file) {
  namespace fs = boost::filesystem;
  fs::path ccT = _ccPath.string() + "/"+ file;
  std::ifstream fileT(ccT.string(), std::ios_base::in | std::ios_base::binary);
  boost::iostreams::filtering_istream inT;
  inT.push(fileT);
  unsigned int i = 0;
  unsigned int pointsNumber = 0;
  for(std::string strT; std::getline(inT, strT); ) {
    if (i == 20) {
      pointsNumber = stoi (strT);
    } else if (i > 21 and i < (pointsNumber + 22)) {
      cc.push_back(stod (strT));
    }
    i++;
  }
}
