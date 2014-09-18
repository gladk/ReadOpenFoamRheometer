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

#include "interpolation.h"
void function_cx_1_func(const alglib::real_1d_array &c, const alglib::real_1d_array &x, double &func, void *ptr) 
  {
      /*
       * 
       * The formula (18) in GraMat. Rheology of weakly wetted granular materials - a comparison of experimental and numerical data.
       * Ruediger Schwarze · Anton Gladkyy · Fabian Uhlig · Stefan Luding, 2013
       * 
       */
       func = 0.5 + 0.5*erf((x[0] - c[0])/c[1]);
  }

bool sortFileTimeCreate(boost::filesystem::path i, boost::filesystem::path j) {
  boost::filesystem::path ii = i.string() + "/uniform/time";
  boost::filesystem::path jj = j.string() + "/uniform/time";
  return (boost::filesystem::last_write_time(ii) < boost::filesystem::last_write_time(jj));
}

bool sortCells(const std::shared_ptr<cell> & i, const std::shared_ptr<cell> & j) {
  if (i->Ccyl()(1) < j->Ccyl()(1)) {              // Z
    return true;
  } else if (i->Ccyl()(1) == j->Ccyl()(1)) {      // Rho
    if (i->Ccyl()(0) < j->Ccyl()(0)) {
      return true;
    } else if (i->Ccyl()(0) == j->Ccyl()(0)) {
      if (i->Ccyl()(2) < j->Ccyl()(2)) {          // Phi
        return true;
      } else {
        return false;
      }
    } else {
      return false;
    }
  } else {
    return false;
  }
  return true;
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
  prepareMatrix();
  
  
  /*
  for (size_t r=0; r<RhoValues.size(); r++) {
    for (size_t c=0; c<ZValues.size(); c++) {
      if (not(_sliceMatrix(r,c))) {
        std::cout<<"!!!!!!!r="<<r<<"; c="<<c<<std::endl; 
      } else {
        std::cout<<"!r="<<r<<"; "<<(_sliceMatrix(r,c))->Ccyl()(0)<<"c="<<c<<"; "<<(_sliceMatrix(r,c))->Ccyl()(1)
          <<"   U:["<<(_sliceMatrix(r,c))->U(5)(0)<<"; "<<(_sliceMatrix(r,c))->U(5)(1)<<"; "<<(_sliceMatrix(r,c))->U(5)(2)<<"]"<<std::endl; 
      }
    }
  }
  */
  
  const double maxR = (_sliceMatrix(_sliceMatrix.rows()-1,0))->Ccyl()(0);
  std::cout<<"maxR: "<<maxR<<std::endl;
  double maxOmega = 0;
  
  std::cout<<"_sliceMatrix.cols(): "<<_sliceMatrix.cols()<<std::endl;
  std::cout<<"_sliceMatrix.rows(): "<<_sliceMatrix.rows()<<std::endl;
  
  
  ofstream exportFile ("./gnuplot_daten");
  ofstream exportFileBand ("./gnuplot_daten_band");
  
  exportFile     << "#001_tId\t002_t\t003_rId\t004_r\t005_zId\t006_z\t007_omega\t008_omegaNorm\n";
  exportFileBand << "#001_tId\t002_t\t003_zId\t004_z\005_RZ\t\006_W\n";
  
  for (size_t t=0; t<_time.size(); t++) {
    for (int c=0; c<_sliceMatrix.cols(); c++) {
      if (_sliceMatrix(_sliceMatrix.rows()-1,c)) {
        maxOmega = (_sliceMatrix(_sliceMatrix.rows()-1,c))->U(t)(2) / (2*M_PI*_sliceMatrix(_sliceMatrix.rows()-1,c)->Ccyl()(0));
      }
      alglib::real_2d_array x;
      alglib::real_1d_array y;
      alglib::real_1d_array cA = "[0.08, 0.0075]";
      double epsf = 0;
      double epsx = 0.0000001;
      alglib::ae_int_t maxits = 0;
      alglib::ae_int_t info;
      alglib::lsfitstate state;
      alglib::lsfitreport rep;
      double diffstep = 0.000001;
      
      x.setlength(_sliceMatrix.rows(), 1);
      y.setlength(_sliceMatrix.rows());
      
      /*
      for(unsigned int r=0; r<_cfg->SecRadial(); r++) {
        x(r,0) = (this->getBand(r,0))->midLinedR();
        y(r) = this->getBand(r,h)->omegaNorm();
      }
      */
      for (int r=0; r<_sliceMatrix.rows(); r++) {
        if (not(_sliceMatrix(r,c))) {
          //std::cout<<"!!!!!!!r="<<r<<"; c="<<c<<std::endl; 
        } else {
          const double Omega = (_sliceMatrix(r,c))->U(t)(2) / (2*M_PI*_sliceMatrix(r,c)->Ccyl()(0));
          exportFile << t << "\t";                                // 001_tId
          exportFile << _time[t] << "\t";                         // 002_t
          exportFile << r << "\t";                                // 003_rId
          exportFile << _sliceMatrix(r,c)->Ccyl()(0) << "\t";     // 004_r
          exportFile << c << "\t";                                // 005_zId
          exportFile << _sliceMatrix(r,c)->Ccyl()(1) << "\t";     // 006_z
          exportFile << Omega << "\t";                            // 007_omega
          exportFile << Omega/maxOmega << "\t";                   // 008_omegaNorm
          exportFile << "\n";                                     //
          //std::cerr<<"t-"<<t<<"; r-"<<_sliceMatrix.rows()-1<<"; c-"<<c<<";   OmegaNorm: "<<Omega/maxOmega<<std::endl;
          x(r,0) = _sliceMatrix(r,c)->Ccyl()(0);
          y(r) = Omega/maxOmega;
        }
      }
      
      lsfitcreatef(x, y, cA, diffstep, state);
      lsfitsetcond(state, epsf, epsx, maxits);
      alglib::lsfitfit(state, function_cx_1_func);
      lsfitresults(state, info, cA, rep);
      
      const double Rz = cA(0);
      const double W = cA(1);
      exportFileBand << t << "\t";                                // 001_tId
      exportFileBand << _time[t] << "\t";                         // 002_t
      exportFileBand << c << "\t";                                // 003_zId
      exportFileBand << _sliceMatrix(1,c)->Ccyl()(1) << "\t";     // 004_z
      exportFileBand << Rz << "\t";                               // 005_RZ
      exportFileBand << W << "\t";                                // 006_W
      exportFileBand << "\n";                                     //
    }
  }
  exportFile.close();
  exportFileBand.close();
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
      std::shared_ptr<cell> cellTmp = std::make_shared<cell>(Eigen::Vector3d(ccx[i],ccy[i],ccz[i]));
      _cells.push_back(cellTmp);
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
        _cells[i-22]->addU(U);
      }
      i++;
    }
  }
  
  std::sort(_cells.begin(), _cells.end(), sortCells);
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

void rheometer::prepareMatrix() {
  std::list<double> RhoValues;
  std::list<double> ZValues;
  unsigned int i = 0;
  BOOST_FOREACH(std::shared_ptr <cell> c, _cells) {
    //std::cout<<std::setprecision(10)<<c->Ccyl()(1)<<" "<<c->Ccyl()(0)<<" "<<c->Ccyl()(2)<<std::endl;
    RhoValues.push_back(c->Ccyl()(0));
    ZValues.push_back(c->Ccyl()(1));
    i++;
  }
  RhoValues.sort(); RhoValues.unique();
  ZValues.sort();   ZValues.unique();
  
  i=0;
  BOOST_FOREACH(double rho, RhoValues) {
    _RhoMap[rho]=i;
    i++;
  }
  i=0;
  BOOST_FOREACH(double z, ZValues) {
    _ZMap[z]=i;
    i++;
  }
  
  _sliceMatrix.resize(RhoValues.size(), ZValues.size());
  _sliceMatrix.fill(nullptr);
    
  for (unsigned d=0; d<_cells.size(); d++) {
    if (d==0 or 
        (Eigen::Vector2d(_cells[d]->Ccyl()(0),_cells[d]->Ccyl()(1)) == Eigen::Vector2d(_cells[d-1]->Ccyl()(0),_cells[d-1]->Ccyl()(1)))) {
      _cellVectorTmp.push_back(_cells[d]);
    } else {
      calculateAvgU();
    }
  }
}

void  rheometer::calculateAvgU() {
  std::shared_ptr<cell> cellTmp = std::make_shared<cell>(Eigen::Vector3d((_cellVectorTmp[0])->c()(0), (_cellVectorTmp[0])->c()(1), (_cellVectorTmp[0])->c()(2)));
  for (unsigned int t = 0; t < _time.size(); t++) {
    double avgU = 0;
    for (unsigned int z = 0; z < _cellVectorTmp.size(); z++) {
      avgU+=(_cellVectorTmp[z])->Ucyl(t)(2);
    }
    // We do consider only Phi-direction!
    avgU/=_cellVectorTmp.size();
    cellTmp->addU(Eigen::Vector3d(0.0, 0.0, avgU));
  }
  _sliceMatrix(_RhoMap[(_cellVectorTmp[0])->Ccyl()(0)], _ZMap[(_cellVectorTmp[0])->Ccyl()(1)]) = cellTmp;
  _cellVectorTmp.clear();
}
