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

## How to setup build environment for PowerGrid with Vagrant
Vagrant (http://www.vagrantup.com) is a tool for creating, managing, and distributing reproducable development environments. In short Vagrant uses an existing virtualization provider (like VirtualBox or AWS) to download, provision, and run a virtual machine (VM).

We use Vagrant with an Ubuntu 15.04 box and a shell script provisioner. The Vagrantfile supports VirtualBox and Parallels provisioners as well as AWS with a few caveats.
 
### To use Vagrant to setup your build environment on a local provider (VMWare, VirtualBox, and Parallels)

Step 1: Install Vagrant from www.vagrantup.com. Binaries exist for OS X, Windows, and Linux. These instructions were made using Vagrant 1.7.2

Step 2: Clone the PowerGrid git repository 

    git clone git@bioe-mrfil-07.bioen.illinois.edu:mrfil/PowerGrid.git
    cd PowerGrid

Step 3: Bring up your VM with the following command

    vagrant up
    
Step 4: Log into your VM. The PowerGrid repo will be mounted in the virtual machine at /vagrant. This folder is the PowerGrid repo directory mounted in the virtual machine.  

    vagrant ssh 
    cd /vagrant
    mkdir build
    cd build
    cmake ../
    make
    
 
### To use Vagrant to setup your build environment on an AWS (Eucalyptus like CUBIC) provider.

These instructions assume your account is setup properly for the cloud and that you have already created AWS keys and a keypair.

Step 1: Install Vagrant from www.vagrantup.com. Binaries exist for OS X, Windows, and Linux. These instructions were made using Vagrant 1.7.2

Step 2: Install the vagrant-aws plugin

    vagrant plugin install vagrant-aws

Step 3: Ensure that your AWS credentials and keys are located in your .bashrc file. Also, ensure that ssh-agent is running for secure key forwarding.

~/.bashrc

    export AWS_ACCESS_KEY=PASTE ACCESS KEY HERE
    export AWS_SECRET_KEY=PASTE SECRET KEY HERE
    export AWS_SSH_KEY=NAME OF SSH KEY IN EUCA
    ssh-agent
    
Step 4: Make sure that you have ssh keys for passwordless ssh login setup on a lab workstation. It is strongly reccomended that you have no passphrase on your ssh key to allow for automated ssh from vagrant.

Step 5: Clone the PowerGrid git repository 

    git clone git@bioe-mrfil-07.bioen.illinois.edu:mrfil/PowerGrid.git
    cd PowerGrid

Step 6: Bring up your VM with the following command

    vagrant up
    
Step 7: Log into your VM. The PowerGrid repo will be mounted in the virtual machine at /vagrant. This folder is the PowerGrid repo directory mounted in the virtual machine. 

    vagrant ssh 
    cd /vagrant
    mkdir build
    cd build
    cmake ../
    make
    
## Troubleshooting:
    ### Vagrant complains that sshfs doesn't work? I can't see my local files in the virtual machine at /vagrant on AWS/Euca?
    
    This happens when the guest (VM) can't login to your host (read physical machine). Check that your SSH Keys are setup correctly. Also, using SSHFS to create a shared file system betwen host and guest relies on SSH connections. If your computer is not reachable at the fully qualified domain name, then the sshfs will fail. Use one of the other providers instead.

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

Step 2: Clone the PowerGrid git repository
    
    git clone git@bioe-mrfil-07.bioen.illinois.edu:mrfil/PowerGrid.git

Step 3: Make a build directory in the root of your PowerGrid repository

    cd PowerGrid
    mkdir build
    
Step 4: Use the cmake command to setup your build directory from within your build directory
    
    cd build
    cmake ../ -DCMAKE_CXX_COMPILER=g++ -DCXX_FLAGS="-O2"
    
Step 5: Make the software in the buld directory. Use -j# to build with multiple threads.

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
