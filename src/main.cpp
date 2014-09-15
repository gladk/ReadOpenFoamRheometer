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

namespace po = boost::program_options;

int main(int ac, char* av[])
{
  std::string outputFolder;
  std::string inputFolder;
  std::string ccFolder;
  
  
  
  std::cout<<"\n\
ReadOpenFoam\n\
Copyright (C) 2014 TU Bergakademie Freiberg\nInstitute for Mechanics and Fluid Dynamics\n\
This program comes with ABSOLUTELY NO WARRANTY.\n\
";
  
  try {
    po::options_description desc("Allowed options");
    desc.add_options()
      ("help", "produce help message")
      ("input,i",  po::value<std::string>()->default_value("."), "input folder")
      ("output,o", po::value<std::string>()->default_value("output"), "output folder")
      ("cc,c", po::value<std::string>()->default_value("."), "folder with ccx, ccy and ccz files (cell centers)")
    ;
    
    po::positional_options_description p;
    p.add("input", -1);
    po::variables_map vm;        
    po::store(po::command_line_parser(ac, av).
    options(desc).positional(p).run(), vm);
    po::notify(vm);  
    
    if (vm.count("help")) {
      std::cout << desc ;
      return 0;
    }
    
    if (vm.count("output")) {
      std::cout << "output folder: " << vm["output"].as<std::string>()<<std::endl ;
    }
    outputFolder = vm["output"].as<std::string>();
    
    if (vm.count("input")) {
      std::cout << "input folder: " << vm["input"].as<std::string>()<<std::endl ;
    }
    inputFolder = vm["input"].as<std::string>();
    if (vm.count("cc")) {
      std::cout << "cc folder: " << vm["cc"].as<std::string>()<<std::endl ;
    }
    ccFolder = vm["cc"].as<std::string>();
  }
  
  catch(std::exception& e) {
      std::cerr << "error: " << e.what()<<std::endl;
      exit (EXIT_FAILURE);
  }
  catch(...) {
      std::cerr << "Exception of unknown type!"<<std::endl;
  }
  
  std::shared_ptr<rheometer> rheometerTmp = std::make_shared<rheometer>(outputFolder, inputFolder, ccFolder);
  
  //=====================================================
  
  return 0;
}
