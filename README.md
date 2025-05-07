# SM-FMM
Solid Mechanics Fast Multipole Method based on elastic Green's function and accelerated using FFTs (O4893).

## Compile and run

1. Install miniconda or anaconda. Set up proxy in .condarc (in users folder under your profile).  
    For Windows, it may be necessary to add proxy servers directly in Navigator under 
    File->Preferences->Configure Conda. Example .condarc for LANL win is below:
    ```
    channels:
    - defaults
    
    auto_update_conda: false
    
    proxy servers:
      http: http://proxyout.lanl.gov:8080
      https: https//proxyout.lanl.gov:8080
   ```
2. From anaconda prompt create enviroment and install packages.  
    For osx:
    ```
    conda create -n SMFMM
    conda activate SMFMM
    conda install gfortran -c conda-forge  
    conda install lapack -c conda-forge
    conda install fftw -c conda-forge
    conda install make -c conda-forge
    ```

    Last step is not realy needed but helps set the enviroment variables properly.

    For linux:
   
    It may be necessary to install appropriate version of lapack (and or other packages).
    It was found that on some Linux systems version 3.6.0 from conda forge (oldest) works well,
    while newest does not. We can specify the version by using
    ```
    conda install lapack=3.6.0=* -c conda-forge
    ```
    Packages can be searched on conda using
    ```
    conda search lapack -c conda-forge
    ```
    For linux, following packages appear to be working well (in this order)
    ```
    conda install liblapack=3.9.0=24_linux64_mkl -c conda-forge
    conda install gfortran=14.2.0=h96c4ede_1 -c conda-forge
    conda install fftw=3.3.10=nompi_hf1063bd_110 -c conda-forge
    conda install make=4.4.1=hb9d3cd8_2 -c conda-forge
    ```
    For Win:

    ```
    conda create -n SMFMM
    conda activate SMFMM
    conda install make -c conda-forge
    conda install lapack -c conda-forge
    conda install fftw -c conda-forge
    ```
    Appropriate version of gfortran can be installed as follows.
     - Download archive (tested working version)
         GCC 8.5.0 + MinGW-w64 9.0.0 (MSVCRT) - release 1
           - Win64: 7-Zip archive* | Zip archive
       from https://winlibs.com/.
     - Unzip to desired folder.
     - Add path folder\mingw64\bin to User Enviroment Variable named Path.
   
3. Adjust the FFTW3 variable in makefile to be directory of your enviroment in conda  
    (e.g. /Users/miroslavzecevic/miniconda3/envs/LSEVPFFT/).LAPACK variable in makefile may 
    be left blank or may need to be specified as well (depending on the system).

4. Compile LS-EVPFFT (serial, serial debug, openmp parallelized) using following comands.
    ```
    make debug
    make opt
    make openmp
    ```
5. Run  
    serial  
    osx and linux:
    ```
    ./../src/FMM.out
    ```
    Win:
    ```
    ./../src/LS-EVPFFT.exe
    ```
    
    openmp version  
    osx and linux:
    ```
    ./../src/FMM.out --nthreads n
    ```
    Win:
    ```
    ./../src/LS-EVPFFT.exe --nthreads n
    ```
    where n is number of threads.



## BSD 3-Clause License

Copyright (c) 2025, Los Alamos National Laboratory

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
