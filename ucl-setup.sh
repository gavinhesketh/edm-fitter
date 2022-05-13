##setup various things:
export PATH=/cvmfs/sft.cern.ch/lcg/releases/LCG_101/CMake/3.20.0/x86_64-centos7-gcc11-opt/bin/:$PATH
export LD_LIBRARY_PATH=/cvmfs/sft.cern.ch/lcg/releases/LCG_101/vdt/0.4.3/x86_64-centos7-gcc11-opt/lib:$LD_LIBRARY_PATH 
export LD_LIBRARY_PATH=/cvmfs/sft.cern.ch/lcg/releases/LCG_101/tbb/2020_U2/x86_64-centos7-gcc11-opt/lib:$LD_LIBRARY_PATH

source /cvmfs/sft.cern.ch/lcg/releases/gcc/11.2.0/x86_64-centos7/setup.sh
source /cvmfs/sft.cern.ch/lcg/releases/LCG_101/ROOT/6.24.06/x86_64-centos7-gcc11-opt/bin/thisroot.sh 

export CMAKE_PREFIX_PATH=/cvmfs/sft.cern.ch/lcg/releases/LCG_101/jsonmcpp/3.9.1/x86_64-centos7-gcc11-opt:$CMAKE_PREFIX_PATH
# cmake -Dnlohmann_json_ROOT=/cvmfs/sft.cern.ch/lcg/releases/LCG_101/jsonmcpp/3.9.1/x86_64-centos7-gcc11-opt ..


echo "*********************** mkedm to compile *******************"
alias mkedm="cd /unix/muons/g-2/hesketh/edm-fitter/build ; make -j8 install ; cd -"
