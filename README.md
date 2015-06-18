# PowerGrid

##  Software for CPU and GPU accelerated iterative magnetic resonance imaging reconstruction. 

##  Subfolders:
*   /PowerGrid - This contains the code directly used by the project
*   /Resources - This contains code that is not directly used by the project for reference. These files are not compiled into the project from this folder.
*   /Support - Convenience headers and support files from other sources. These files are compiled into the project.

## Depedenencies 
*   Xerces-C++ - XML Parser used with the CodeSynthesis Generated files for parsing our config.xml files.

*   Armadillo (http://arma.sourceforge.net) - Templated Linear Algebra library. Gives us MATLAB like syntax in C++

*   MATIO  (http://matio.sourceforge.net) - Library for reading and writing MATLAB .mat files in C/C++

*   FFTW - () -Fastest Fourier Transform in the West - Used for CPU implementations of the FFTs used in Gridding.

### These dependencies are requires for running on Blue Waters
ACML - AMD Core Math Library - A set of math routines including LAPACK and BLAS for running on Opteron processors like on Blue Waters and other Cray systems.



## CPU Build Instructions with PGI

Step 1: Setup and install all required dependencies. 

Step 2: Setup your .bashrc to include the PGI shared libraries

Step 3: Clone the PowerGrid git repository
    
    git clone git@bioe-mrfil-07.bioen.illinois.edu:mrfil/PowerGrid.git

Step 4: Make a build directory in the root of your PowerGrid repository

    cd PowerGrid
    mkdir build
    
Step 5: Use the cmake command to setup your build directory from within your build directory
    
    cd build
    cmake ../ -DCMAKE_CXX_COMPILER=pgc++ -DCXX_FLAGS="-fast"
    
Step 6: Make the software in the buld directory. Use -j# to build with multiple threads.

    make -j2
    
## CPU Build Instructions with GCC

Step 1: Load all required dependencies by loading PowerGridSupport Module

    module load PowerGridSupport/1.0

Step 3: Clone the PowerGrid git repository
    
    git clone git@bioe-mrfil-07.bioen.illinois.edu:mrfil/PowerGrid.git

Step 4: Make a build directory in the root of your PowerGrid repository

    cd PowerGrid
    mkdir build
    
Step 5: Use the cmake command to setup your build directory from within your build directory
    
    cd build
    cmake ../ -DCMAKE_CXX_COMPILER=g++ -DCXX_FLAGS="-O2"
    
Step 6: Make the software in the buld directory. Use -j# to build with multiple threads.

    make -j2
    
## OpenACC Build Instructions with PGI
   
Step 1: Setup and install all required dependencies. 
   
Step 2: Setup your .bashrc to include the PGI shared libraries
  
Step 3: Clone the PowerGrid git repository
       
    git clone git@bioe-mrfil-07.bioen.illinois.edu:mrfil/PowerGrid.git
   
Step 4: Make a build directory in the root of your PowerGrid repository
   
    cd PowerGrid
    mkdir build
       
Step 5: Use the cmake command to setup your build directory from within your build directory
       
    cd build
    cmake ../ -DCMAKE_CXX_COMPILER=pgc++ -DCXX_FLAGS="-fast -ta=tesla,nordc -Minfo=acc"
       
Step 6: Make the software in the buld directory. Use -j# to build with multiple threads.
   
    make -j2