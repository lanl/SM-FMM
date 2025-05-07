! This just wraps the FFTW3 C-interface into a Fortran module for ease of use
module FFTW3
  use, intrinsic :: iso_c_binding
  include 'fftw3.f03'
end module
