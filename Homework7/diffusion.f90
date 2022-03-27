program diffusion
implicit none

include "constants.h"

! number of timesteps
integer, parameter :: NSTEP = 4000000
 
! time step in seconds
double precision, parameter :: DT = 10000000. ! s
 
! logicals: 
logical, parameter :: hetero_kappa = .true. 
logical, parameter :: hetero_spacing = .false. 

character(60) :: outdir   = "hetero_k_homo_sp/"
character(60) :: fnameout = "t"
character(60) makedirectory, removedir


! fixed boundary conditions
logical, parameter :: FIXED_BC = .true.
 
! model parameters  (SI)
double precision, parameter :: LENGTH = 3.0d+03 ! m
double precision, parameter :: DENSITY = 2.5d+03 ! kg/m^3
double precision, parameter :: THERMALCONDUCTIVITY = 10.0d-01 , TC2 = 2.0d-01 ! units cal/m/s/K - TC2 used for Right half of domain if hetero_kappa = true
double precision, parameter :: HEATCAPACITY = 0.3d+03 ! cal/kg/K

integer ispec,i,j,iglob,itime,a,b,g, ctr

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
double precision stiffness, val 

! global mass matrix
double precision, dimension(NGLOB) :: mass_global = 0.

! WE - global rhs vector
double precision, dimension(NGLOB) :: RHS = 0.

! temperature and temperature time derivative
double precision, dimension(NGLOB) :: temperature =0. ,dtemperature_dt = 0.

! local to global numbering
integer, dimension(NGLL,NSPEC) :: ibool

! time marching
double precision deltat,deltatover2
double precision dh,diffusivity

! end fluxes
double precision flux_1,flux_NGLOB

! end temperatures
double precision temperature_1,temperature_NGLOB

! derivatives
double precision dtdx,flux,templ,temp(NGLL)

! movie
character(len=50) moviefile

! plot? 
character(60) :: plot  


! CMD for plotting? 
call get_command_argument(1, plot)


! Make output directory
print*, "Removing existing directory:   ", trim(outdir)

removedir = 'rm -r ' // trim(outdir)
call system(removedir)

makedirectory = 'mkdir ' // trim(outdir)
print*, "Creating new directory:        ", trim(outdir)
call system(makedirectory)

!++++++++++++++++++++++++++++++++++++++++++++++++++
print *, "Calculating matrices..."

call define_derivative_matrix(xigll,wgll,hprime)


! set up local to global numbering
iglob = 1
do ispec = 1,NSPEC
  do i = 1,NGLL
  if(i>1)iglob = iglob+1
    ibool(i,ispec) = iglob

  enddo
enddo


! Determine element spatial distribtion depending on if homogeneous or heterogeneous 
if (hetero_spacing) then 
    ! non evenly-spaced: 
    do ispec = 1, NSPEC 
      if (ispec <4) then 
        x1(ispec) = dble(ispec-1) *500.
      x2(ispec) =dble(ispec)*500. 
      else 
        val = 1500.0/9.0
        x1(ispec) = ((dble(ispec-1)-3.0)*val ) + 1500.
        x2(ispec) = ((dble(ispec)-3.0)*val ) + 1500.
      endif 
    enddo 
else 
    ! evenly spaced anchors between 0 and 1
    do ispec = 1,NSPEC
      x1(ispec) = LENGTH*dble(ispec-1)/dble(NSPEC)
      x2(ispec) = LENGTH*dble(ispec)/dble(NSPEC)
    enddo 
endif 


! set up the mesh properties
do ispec = 1,NSPEC
  do a = 1,NGLL
    rho(a,ispec) = DENSITY
    heat_capacity(a,ispec) = HEATCAPACITY
    dxidx(a,ispec) = 2. / (x2(ispec)-x1(ispec))
    jacobian(a,ispec) = (x2(ispec)-x1(ispec)) / 2.

    
    thermal_conductivity(a,ispec) = THERMALCONDUCTIVITY
  enddo
enddo



! get the global grid points & calculate mass matrix 
do ispec = 1,NSPEC
  do a = 1,NGLL
    ! Determine global x values
    iglob = ibool(a,ispec)
    x(iglob) = 0.5*(1.-xigll(a))*x1(ispec)+0.5*(1.+xigll(a))*x2(ispec)

    ! Apply heterogeneity of kappa 
    if (hetero_kappa) then 
      if (x(iglob) > LENGTH/2.0) then
        thermal_conductivity(a,ispec) = TC2
      endif 
    endif

    ! Mass matrix (vector) calculation - do in this loop for efficiency :
    do b = 1, NGLL
      mass_global(iglob) = mass_global(iglob) + wgll(b)*rho(b,ispec)*heat_capacity(b,ispec)*jacobian(b,ispec)
    enddo

  enddo
enddo


! Set initial conditions: 
temperature(:) = 0.
temperature(1) = 10.  
dtemperature_dt(:) = 0.

ctr = 0 

print *, "Running timesteps..."

! Loop in time (time-marching)
do itime = 1, NSTEP

  ! 1.1) PREDICTOR STEP: 
  temperature = temperature + 0.5*DT*dtemperature_dt

  ! 1.2) ENFORCE BOUNDARY CONDITIONS 
  temperature(1) = 10.0
  temperature(NGLOB) = 0.0

  ! 2.1) SOLVE FOR THE RHS 
    RHS(:) = 0. ! reinitialise RHS for each time loop 
    do ispec = 1, NSPEC
        do a = 1, NGLL
          do b = 1, NGLL
            ! reset/initialise the local stiffness (K_ab)
            stiffness = 0.       
      
            do g = 1, NGLL
              ! Loop through summation for quadrature - g is gamma in equations 
              stiffness = stiffness + wgll(g)*thermal_conductivity(g,ispec)*dxidx(g,ispec)*hprime(a, g)*hprime(b, g)
            enddo 
          
            j = ibool(a,ispec)
            RHS(j) = RHS(j)   -(stiffness * temperature(ibool(b,ispec)))
        
          enddo
        enddo
    enddo 

  ! 2.2) Estimate new dtemp_dt 
  do i = 1,NGLOB
    dtemperature_dt(i) = RHS(i)/mass_global(i)
  enddo 

  ! 3) CORRECTOR STEP: 
  temperature = temperature + 0.5*DT*dtemperature_dt


  ! 4) write out snapshots
    
    if(mod(itime-1,20000) == 0) then
      ! WE edit - make a file to store the timestaps outputted
      open(unit=90,file=trim(outdir)//"time", status='unknown', position="append")
      write(90,*) itime
      close(90)
      
      write(moviefile,'(a,i8.8)') trim(outdir)//trim(fnameout), itime
      open(unit=10,file=trim(moviefile),status='unknown')
      do iglob = 1,NGLOB
        write(10,*) sngl(x(iglob)),sngl(temperature(iglob))
      enddo
      close(10)
      ctr = ctr + 1 
    endif

enddo ! end time loop

open(unit=91,file= trim(outdir)//"meta", status='unknown', position="append")
    write(91,*) fnameout
    write(91,*) ctr 
    write(91,*) DT 
    close(91)


! Run python plotter:
if (plot=="plot") then 
  print *, "Finished calculations. Starting plots:"
  call system("python3 plot.py "//trim(outdir))
endif 

end program diffusion
!======================================================================
