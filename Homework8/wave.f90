program wave
implicit none

include "constants.h"

! number of timesteps
integer, parameter :: NSTEP = 15000
integer, parameter :: slice_interval = 10
 
! time step in seconds
double precision, parameter :: DT = 0.01 ! s
double precision, parameter :: DT2 = DT*DT

! Output directory and file names: 
character(60) :: outdir   = "test_wave/"
character(60) :: fnameout = "w"
character(60) makedirectory, removedir

! boundary conditions: if not true then uses neumann
logical, parameter :: dirichlet_BC = .true.
 
! model parameters  (SI)
double precision, parameter :: RHO     = 1 ! kg/m^3
double precision, parameter :: MU      = 1 ! 
double precision, parameter :: LENGTH  = 100 ! m


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
double precision, dimension(NGLL,NSPEC) :: rho_matrix, mu_matrix

! Jacobian inverse `matrix' and Jacobian
double precision, dimension(NGLL,NSPEC) :: dxidx,jacobian

! local stiffness matrix
double precision stiffness

! global mass matrix
double precision, dimension(NGLOB) :: m_glob = 0.

! WE - global rhs vector
double precision, dimension(NGLOB) :: RHS = 0.

! displacement, velocity, acceleration
double precision, dimension(NGLOB) :: dis = 0. , vel = 0., acc = 0.

! local to global numbering
integer, dimension(NGLL,NSPEC) :: ibool


! movie
character(len=50) moviefile

! plot - takes command line input - if contains string "plot", i.e. the cmd line is 
! ./xwave plot  then it will run a python script to animate the results
character(60) :: plot  


! CMD line input for plotting
call get_command_argument(1, plot)

!++++++++++++++++++++++++++++++++++++++++++++++++++
! Delete/re-create directories for slice outputs:  
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
do ispec = 1,NSPEC
  x1(ispec) = LENGTH*dble(ispec-1)/dble(NSPEC)
  x2(ispec) = LENGTH*dble(ispec)/dble(NSPEC)
enddo 


! set up the mesh properties
do ispec = 1,NSPEC
  do a = 1,NGLL
    rho_matrix(a,ispec) = RHO
    mu_matrix(a,ispec)  = MU
    dxidx(a,ispec) = 2. / (x2(ispec)-x1(ispec))
    jacobian(a,ispec) = (x2(ispec)-x1(ispec)) / 2.
  enddo
enddo


m_glob = 0.
! get the global grid points & calculate mass matrix 
do ispec = 1,NSPEC
  do a = 1,NGLL
    iglob = ibool(a,ispec)  ! Determine global x values
    x(iglob) = 0.5*(1.-xigll(a))*x1(ispec)+0.5*(1.+xigll(a))*x2(ispec)    

    ! Mass matrix (vector) calculation 
    m_glob(iglob) = m_glob(iglob) + wgll(a)*rho_matrix(a,ispec)*jacobian(a,ispec)
  enddo
enddo


! Set initial force: 
do i = 1, NGLOB 
  ! I think the HW wants a force of exp(...) so that the acceleration is given as initiall 0
  ! but only seems to work if the displacement is set and not the acc using a = f/m
  !acc(i) = EXP(-((x(i)-50.0)**2.0)*0.1)/m_glob(i)
  dis(i) = EXP(-((x(i)-50.0)**2.0)*0.1)
  !dis(i) = 0.0
  !vel(i) = 0.0

enddo


ctr = 0 

print *, "Running timesteps..."

! Loop in time (time-marching)
do itime = 1, NSTEP


  ! 1) Predictor step: 
  dis = dis + DT*vel + 0.5*DT2*acc
  vel = vel + 0.5*DT*acc 
  acc = 0.

  ! Apply BCs:
  if (dirichlet_BC) then
    dis(1) = 0.
    dis(NGLOB) = 0.
  endif

  ! 2.1) SOLVE FOR THE RHS 
  RHS(:) = 0. ! reinitialise RHS for each time loop 
  do ispec = 1, NSPEC
      do a = 1, NGLL
        do b = 1, NGLL
          ! reset/initialise the local stiffness (K_ab)
          stiffness = 0.       
    
          do g = 1, NGLL
            ! Loop through summation for quadrature - g is gamma in equations 
            stiffness = stiffness + wgll(g)*mu_matrix(g,ispec)*dxidx(g,ispec)*hprime(a, g)*hprime(b, g)
          enddo 
        
          j = ibool(a,ispec)
          RHS(j) = RHS(j) - (stiffness * dis(ibool(b,ispec)))
      
        enddo
      enddo
  enddo 

  ! 2) Solve: 
  acc(:) =  acc(:) + RHS(:)/m_glob(:)

  ! 3) Corrector: 
  vel(:) = vel(:) + 0.5*DT*acc(:)






  ! 4) write out snapshots
    if(mod(itime-1, slice_interval) == 0) then
      ! WE edit - make a file to store the timestaps outputted
      open(unit=90,file=trim(outdir)//"time", status='unknown', position="append")
      write(90,*) itime
      close(90)
      
      write(moviefile,'(a,i8.8)') trim(outdir)//trim(fnameout), itime
      open(unit=10,file=trim(moviefile),status='unknown')
      do iglob = 1,NGLOB
        write(10,*) sngl(x(iglob)),sngl(dis(iglob))
      enddo
      close(10)
      ctr = ctr + 1 
    endif

enddo ! end time loop

open(unit=91,file= trim(outdir)//"meta", status='unknown', position="append")
    write(91,*) fnameout
    write(91,*) ctr 
    write(91,*) DT 
    write(91,*) slice_interval 
    close(91)


! Run python plotter:
if (plot=="plot") then 
  print *, "Finished calculations. Starting plots:"
  call system("python3 plot.py "//trim(outdir))
endif 

end program wave
! =================================================================================================
