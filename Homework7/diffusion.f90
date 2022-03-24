program diffusion
implicit none

include "constants.h"

! number of timesteps
integer, parameter :: NSTEP = 30000
 
! time step in seconds
double precision, parameter :: DT = 100000000. ! s
 
! fixed boundary conditions
logical, parameter :: FIXED_BC = .true.
 
! model parameters  (SI)
double precision, parameter :: LENGTH = 3.0d+03 ! m
double precision, parameter :: DENSITY = 2.5d+03 ! kg/m^3
double precision, parameter :: THERMALCONDUCTIVITY = 10.0d-01 ! cal/m/s/K
double precision, parameter :: HEATCAPACITY = 0.3d+03 ! cal/kg/K

integer ispec,i,j,iglob,itime

! Gauss-Lobatto-Legendre points of integration
double precision, dimension(NGLL) :: xigll
 
! weights
double precision, dimension(NGLL) :: wgll
 
! array with derivatives of Lagrange polynomials
double precision, dimension(NGLL,NGLL) :: hprime

! anchors
double precision, dimension(NSPEC) :: x1,x2

! global grid points
double precision, dimension(NGLOB) :: x

! material properties
double precision, dimension(NGLL,NSPEC) :: rho,heat_capacity,thermal_conductivity

! Jacobian `matrix' and Jacobian
double precision, dimension(NGLL,NSPEC) :: dxidx,jacobian

! local mass matrix
double precision mass_local

! global mass matrix
double precision, dimension(NGLOB) :: mass_global

! temperature and temperature time derivative
double precision, dimension(NGLOB) :: temperature,dtemperature_dt

! local to global numbering
integer, dimension(NGLL,NSPEC) :: ibool

! time marching
double precision deltat,deltatover2
double precision dh,diffusivity,time_step

! end fluxes
double precision flux_1,flux_NGLOB

! end temperatures
double precision temperature_1,temperature_NGLOB

! derivatives
double precision dtdx,flux,templ,temp(NGLL)

! movie
character(len=50) moviefile

!++++++++++++++++++++++++++++++++++++++++++++++++++

call define_derivative_matrix(xigll,wgll,hprime)

! evenly spaced anchors between 0 and 1
do ispec = 1,NSPEC
  x1(ispec) = LENGTH*dble(ispec-1)/dble(NSPEC)
  x2(ispec) = LENGTH*dble(ispec)/dble(NSPEC)
enddo

! set up the mesh properties
do ispec = 1,NSPEC
  do i = 1,NGLL
    rho(i,ispec) = DENSITY
    thermal_conductivity(i,ispec) = THERMALCONDUCTIVITY
    heat_capacity(i,ispec) = HEATCAPACITY
    dxidx(i,ispec) = 2. / (x2(ispec)-x1(ispec))
    jacobian(i,ispec) = (x2(ispec)-x1(ispec)) / 2.
  enddo
enddo

! set up local to global numbering
iglob = 1
do ispec = 1,NSPEC
  do i = 1,NGLL
  if(i>1)iglob = iglob+1
    ibool(i,ispec) = iglob
  enddo
enddo

! get the global grid points
do ispec = 1,NSPEC
  do i = 1,NGLL
    iglob = ibool(i,ispec)
    x(iglob) = 0.5*(1.-xigll(i))*x1(ispec)+0.5*(1.+xigll(i))*x2(ispec)
  enddo
enddo

! calculate the global mass matrix 'mass_global'
!---------------------------------------------------
! put your codes here
!---------------------------------------------------

! estimate the time step 'time_step'
!---------------------------------------------------
! put your codes here
!---------------------------------------------------
print *,'time step estimate: ',time_step,' seconds'

! set up the boundary conditions
!---------------------------------------------------
! put your codes here
!---------------------------------------------------

! initialize
temperature(:) = 0.
dtemperature_dt(:) = 0.

do itime = 1,NSTEP

! update temperature and do not forget to set appropriate
! boundary conditions
!---------------------------------------------------
! put your codes here
!---------------------------------------------------


  
! write out snapshots
  
  if(mod(itime-1,1000) == 0) then
    ! WE edit - make a file to store the timestaps outputted
    open(unit=90,file="snapshots/time", status='unknown', position="append")
    write(90,*) itime
    close(90)


    write(moviefile,'(a,i5.5)')'snapshots/snapshot',itime
    open(unit=10,file=trim(moviefile),status='unknown')
    do iglob = 1,NGLOB
      write(10,*) sngl(x(iglob)),sngl(temperature(iglob))
    enddo
    close(10)
  endif

enddo ! end time loop

end program diffusion
!======================================================================
