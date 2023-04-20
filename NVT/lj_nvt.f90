module parameters
  !ensemble parameters
  integer :: dim_               !number of dimensions (D)
  integer :: N                  !number of particles (N)
  real*8 :: dens                !density
  real*8 :: temp                !temperature
  !simulation parameters
  real*8 :: lbox                !length of box
  real*8 :: volume              !volume of box
  real*8 :: a1,a2               !LJ potential parameters
  integer :: ran                !random seed
  integer :: nseries            !number of series of sampling (MAX 20)
  integer :: nekv               !equilibration length
  !MC parameters
  integer :: ncycl              !number of cycles = N trial moves
  integer :: nsampl             !number of trial moves per sample (same as N)
  real*8 :: delta               !MAX random displacement by sigma
  !MD parameters
  integer :: nsteps             !number of time steps
  real*8 :: tstep,dt2,d2t2      !time step, half step, half squared step
  real*8 :: dt4,dt8             !quarter step, eighth step
  real*8 :: tt                  !temperature time constant
  real*8 :: NfkT                !number of degrees of freedom * temperature
  real*8, dimension(10) :: Q    !heat bath coupling masses (10)
  integer :: nframe             !number of time steps per frame for movie
  !number Pi
  real*8 :: Pi
  PARAMETER (Pi=4.0d0*datan(1.0d0))
  !numberical constant
  real*8 :: c
endmodule parameters

module particles
  real*8, dimension(:,:), allocatable :: pos               !matrix of position vectors
  real*8, dimension(:,:), allocatable :: vel               !matrix of velocity vectors
  real*8, dimension(:,:), allocatable :: acc               !matrix of acceleration vectors
  real*8 :: ekin                                           !kinetic energy
  real*8 :: tem                                            !equipartion theorem temperature
  real*8 :: etot                                           !total internal energy
  real*8, dimension(10) :: vksi                            !heat bath velocities
  real*8, dimension(10) :: G                               !heat bath forces
endmodule particles

module microstate
  real*8, dimension(:,:), allocatable :: pot               !matrix of potential energies
  real*8, dimension(:,:), allocatable :: tla               !matrix of virial pressures
  real*8, dimension(:,:), allocatable :: dis               !matrix of distances
  real*8, dimension(:,:,:), allocatable :: rrr             !matrix of distance vectors
  real*8 :: u,p                                            !potential energy, virial pressure
  real*8 :: rcut                                           !distance of potential energy cutoff
endmodule microstate

module macrostate
  real*8 :: uu,uu2                  !excess internal energy, excess internal energy squared
  real*8 :: pp                      !virial pressure
  integer :: m                      !number of samples
  integer :: accx                   !number of accepted trial moves for displacement
  real*8 :: kx                      !ratio of accepted trial moves for displacement
endmodule macrostate

module correlation
  real*8, dimension(10000) :: gr                !correlation function g(r)
  real*8 :: grinterv                            !g(r) radial interval
  integer :: ngr                                !cycles/steps per g(r)
  integer :: mgr                                !number of samples
endmodule correlation


program core
  use parameters
  use particles
  use microstate
  use macrostate
  use correlation
  implicit none
  real*8 :: time1,time2                         !CPU_time
  integer, dimension(8) :: values1,values2      !value(1)=year, value(2)=month, value(3)=day, value(4)=time difference with UTC in minutes, value(5)=hour, value(6)=minute, value(7)=second, value(8)=milisecond
  character :: tekst                            !reading from input
  character(8) :: fmt                           !string for output format
  integer :: i,j
  real*8 :: cv                                  !excess heat capacity
  real*8 :: therm(25),dtherm(25),thermo(25,20)  !thermodynamic quantities
  real*8 :: dgr(10000),ggr(10000,20)            !correlation function

  call cpu_time(time1)
  call date_and_time(VALUES=values1)
  open(1002,file='log')
  write(1002,*)'-------------------------------------------------------------------------'
  write(1002,*)'start date and time YYYY, MM, DD, UTC, HH, MIN, SS'
  write(1002,'(8i5)')values1

  open(1001,file='param')
  read(1001,*)tekst
  read(1001,*)tekst,dim_
  read(1001,*)tekst,N
  read(1001,*)tekst,dens
  read(1001,*)tekst,temp
  read(1001,*)tekst
  read(1001,*)tekst,a1,a2
  read(1001,*)tekst,ran
  read(1001,*)tekst,nseries
  read(1001,*)tekst,nekv
  read(1001,*)tekst
  read(1001,*)tekst,ncycl
  read(1001,*)tekst,nsampl
  read(1001,*)tekst,delta
  read(1001,*)tekst
  read(1001,*)tekst,nsteps
  read(1001,*)tekst,tstep
  read(1001,*)tekst,tt
  read(1001,*)tekst,nframe
  read(1001,*)tekst
  read(1001,*)tekst,grinterv
  read(1001,*)tekst,ngr
  close(1001)

  allocate( pos(N,dim_) )
  allocate( vel(N,dim_) )
  allocate( acc(N,dim_) )
  allocate( pot(N,N) )
  allocate( tla(N,N) )
  allocate( dis(N,N) )
  allocate( rrr(N,N,dim_) )
  select case(dim_)
  case(1)
    c=2.0d0
  case(2)
    c=2.0d0*Pi
  case(3)
    c=4.0d0*Pi
  endselect
  volume=dble(N)/dens
  lbox=volume**(1.0d0/dim_)
  rcut=lbox/2.0d0
  dt2=tstep/2.0d0
  d2t2=tstep**2/2.0d0
  dt4=tstep/4.0d0
  dt8=tstep/8.0d0
  NfkT=dble(N*dim_)*temp
  Q(1)=NfkT*tt**2
  do i=2,10
    Q(i)=temp*tt**2
  enddo
  vksi=0.0d0
  G(1)=0.0d0
  do i=1,9
    G(i+1)=(Q(i)*vksi(i)**2-temp)/Q(i+1)
  enddo

  write(1002,*)'-------------------------------------------------------------------------'
  write(1002,*)'N=',N
  write(1002,*)'V=',volume
  write(1002,*)'T=',temp
  write(1002,*)'-------------------------------------------------------------------------'
  write(1002,*)'dim_=',dim_
  write(1002,*)'density=',dens
  write(1002,*)'length of box=',lbox
  write(1002,*)'random seed=',ran
  write(1002,*)'MAX random displacement by sigma=',delta
  write(1002,*)'time step=',tstep
  write(1002,*)'-------------------------------------------------------------------------'
  write(1002,*)'Monte Carlo (MC) and Molecular Dynamics (MD) simulation of Lennard-Jones (LJ) particles'

  !----------------------------------------------------------
  !MC
  !----------------------------------------------------------
  call random_pos()
  call snapshot('random')
  write(1002,*)'-------------------------------------------------------------------------'
  write(1002,*)'particles randomly put in box'

  call minimize()
  call snapshot('minpos')
  write(1002,*)'-------------------------------------------------------------------------'
  write(1002,*)'bad contacts removed'
  kx=dble(accx)/dble(10000*N)
  write(1002,*)'acceptance ratio kx=',kx

  call mcequil()
  call snapshot('mcequi')
  uu=uu/dble(m)
  pp=pp/dble(m)
  kx=dble(accx)/dble(nekv*N)
  write(1002,*)'-------------------------------------------------------------------------'
  write(1002,*)'MC equilibration successful'
  write(1002,*)'uu/N=',uu/dble(N)
  write(1002,*)'pp=',pp
  write(1002,*)'acceptance ratio kx=',kx

  open(223,file='thermodynamics')
  write(223,*)'MC'
  write(1002,*)'-------------------------------------------------------------------------'
  write(1002,*)'MC simulation initialized'
  do j=1,nseries
    call mcseries(j)
    write(fmt,'(a2,i4.4)')'mc',j
    call snapshot(fmt)
    uu=uu/dble(m)
    uu2=uu2/dble(m)
    cv=(uu2-uu**2)/temp**2/dble(N)
    pp=pp/dble(m)
    kx=dble(accx)/dble(ncycl*N)
    gr=gr/dble(mgr)
    thermo(1,j)=uu/dble(N)
    thermo(2,j)=pp
    thermo(3,j)=cv
    ggr(:,j)=gr(:)
    write(223,'(5e16.7)')temp,dens,uu/dble(N),pp,cv
    write(1002,*)'-------------------------------------------------------------------------'
    write(1002,*) j,'-th MC run successful'
    write(1002,*)'uu/N=',uu/dble(N)
    write(1002,*)'pp=',pp
    write(1002,*)'cv=',cv
    write(1002,*)'acceptance ratio kx=',kx
  enddo

  open(202,file='mc_thermodynamics')
  write(202,'(2e16.7)')temp,dens
  therm=0.0d0
  dtherm=0.0d0
  do i=1,3
    do j=1,nseries
      therm(i)=therm(i)+thermo(i,j)
    enddo
    therm(i)=therm(i)/dble(nseries)
    do j=1,nseries
      dtherm(i)=dtherm(i)+(thermo(i,j)-therm(i))**2
    enddo
    dtherm(i)=dsqrt(dtherm(i)/dble(nseries))
    write(202,'(2e16.7)')therm(i),dtherm(i)
  enddo
  close(202)

  open(202,file='mc_correlation')
  gr=0.0d0
  dgr=0.0d0
  do i=1,10000
    do j=1,nseries
      gr(i)=gr(i)+ggr(i,j)
    enddo
    gr(i)=gr(i)/dble(nseries)
    do j=1,nseries
      dgr(i)=dgr(i)+(ggr(i,j)-gr(i))**2
    enddo
    dgr(i)=dsqrt(dgr(i)/dble(nseries))
    write(202,'(3f16.7)')(dble(i)-0.5d0)*grinterv,gr(i),dgr(i)
  enddo
  close(202)

  !------------------------------------------------------------------
  !MD
  !------------------------------------------------------------------
  call init_pos()
  call snapshot('mdinit')
  write(1002,*)'-------------------------------------------------------------------------'
  write(1002,*)'positions initialized'

  call init_vel()
  write(1002,*)'-------------------------------------------------------------------------'
  write(1002,*)'velocities initialized'

  call mdequil()
  call snapshot('mdequi')
  uu=uu/dble(m)
  pp=pp/dble(m)
  write(1002,*)'-------------------------------------------------------------------------'
  write(1002,*)'MD equilibration successful'
  write(1002,*)'uu/N=',uu/dble(N)
  write(1002,*)'pp=',pp

  write(223,*)'MD'
  write(1002,*)'-------------------------------------------------------------------------'
  write(1002,*)'MD simulation initialized'
  do j=1,nseries
    call mdseries(j)
    write(fmt,'(a2,i4.4)')'md',j
    call snapshot(fmt)
    write(fmt,'(i4.4)')j
    call velocities(fmt)
    uu=uu/dble(m)
    uu2=uu2/dble(m)
    cv=(uu2-uu**2)/temp**2/dble(N)
    pp=pp/dble(m)
    gr=gr/dble(mgr)
    thermo(1,j)=uu/dble(N)
    thermo(2,j)=pp
    thermo(3,j)=cv
    ggr(:,j)=gr(:)
    write(223,'(5e16.7)')temp,dens,uu/dble(N),pp,cv
    write(1002,*)'-------------------------------------------------------------------------'
    write(1002,*) j,'-th MD run successful'
    write(1002,*)'uu/N=',uu/dble(N)
    write(1002,*)'pp=',pp
    write(1002,*)'cv=',cv
  enddo
  close(223)

  open(202,file='md_thermodynamics')
  write(202,'(2e16.7)')temp,dens
  therm=0.0d0
  dtherm=0.0d0
  do i=1,3
    do j=1,nseries
      therm(i)=therm(i)+thermo(i,j)
    enddo
    therm(i)=therm(i)/dble(nseries)
    do j=1,nseries
      dtherm(i)=dtherm(i)+(thermo(i,j)-therm(i))**2
    enddo
    dtherm(i)=dsqrt(dtherm(i)/dble(nseries))
    write(202,'(2e16.7)')therm(i),dtherm(i)
  enddo
  close(202)

  open(202,file='md_correlation')
  gr=0.0d0
  dgr=0.0d0
  do i=1,10000
    do j=1,nseries
      gr(i)=gr(i)+ggr(i,j)
    enddo
    gr(i)=gr(i)/dble(nseries)
    do j=1,nseries
      dgr(i)=dgr(i)+(ggr(i,j)-gr(i))**2
    enddo
    dgr(i)=dsqrt(dgr(i)/dble(nseries))
    write(202,'(3f16.7)')(dble(i)-0.5d0)*grinterv,gr(i),dgr(i)
  enddo
  close(202)

  call cpu_time(time2)
  call date_and_time(VALUES=values2)
  write(1002,*)'-------------------------------------------------------------------------'
  write(1002,*)'CPU simulation time=',time2-time1
  write(1002,*)'start and finish date and time YYYY, MM, DD, UTC, HH, MIN, SS'
  write(1002,'(8i5)')values1
  write(1002,'(8i5)')values2
  time2=(values2(8)-values1(8))/1000.0d0+values2(7)-values1(7)
  time2=time2+(values2(6)-values1(6))*60.0d0
  time2=time2+(values2(5)-values1(5))*60.0d0*60.0d0
  time2=time2+(values2(3)-values1(3))*60.0d0*60.0d0*24.0d0
  write(1002,*)'real simulation time=',time2
  close(1002)
endprogram core

!----------------------------------------------------------------------------------------------------
!Generates random initial positions for particles inside the box without overlapping
!----------------------------------------------------------------------------------------------------
subroutine random_pos()
  use parameters
  use particles
  implicit none
  real*8 :: ran3,image        !functions
  real*8 :: a                 !random number
  real*8 :: rvec(dim_),r      !coordinates
  integer :: i,j,k,ii,jj,kk

  !first particle
  do ii=1,dim_
    a=ran3(ran)
    pos(1,ii)=(a-0.5d0)*lbox
  enddo

  !prevent overlapping
  do i=2,N
    j=0
    do while (j.lt.0.5)
      do jj=1,dim_
        a=ran3(ran)
        pos(i,jj)=(a-0.5d0)*lbox
      enddo
      j=1
      do k=1,i-1
        do kk=1,dim_
          rvec(kk)=image(pos(i,kk),pos(k,kk),lbox)
        enddo
        r=dsqrt(sum(rvec**2))
        if (r.lt.0.5d0) j=0
      enddo
    enddo
  enddo
endsubroutine random_pos

!----------------------------------------------------------------------------------------------------
!Ran3 random number generator generates uniformly distributed random real numbers [0,1)
!----------------------------------------------------------------------------------------------------
function ran3(idum)
  implicit none
  integer :: idum  !negative seed restarts the initialization
  integer :: MBIG,MSEED,MZ
  real*8 :: ran3,FAC
  PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1.0d0/dble(MBIG))
  integer :: i,iff,ii,inext,inextp,k
  integer :: mj,mk,ma(55)  !the value 55 is special, see Knuth
  SAVE iff,inext,inextp,ma  !static variables
  DATA iff /0/  !first initialization
  
  if (idum.lt.0.or.iff.eq.0) then  !initialization
    iff=1
    mj=MSEED-iabs(idum)
    mj=mod(mj,MBIG)
    ma(55)=mj
    mk=1
    do i=1,54
      ii=mod(21*i,55)
      ma(ii)=mk
      mk=mj-mk
      if (mk.lt.MZ) mk=mk+MBIG  !integer overflow
      mj=ma(ii)
    enddo
    do k=1,4  !warming up the generator
      do i=1,55
        ma(i)=ma(i)-ma(1+mod(i+30,55))
        if (ma(i).lt.MZ) ma(i)=ma(i)+MBIG  !integer overflow
      enddo
    enddo
    inext=0  !prepare the indices
    inextp=31  !the value 31 is special, see Knuth
    idum=1  !the seed is set to one, if set to a negative number the next call will restart the initialization
  endif
  inext=inext+1
  if (inext.eq.56) inext=1  !stay in range
  inextp=inextp+1
  if (inextp.eq.56) inextp=1  !stay in range
  mj=ma(inext)-ma(inextp)
  if (mj.lt.MZ) mj=mj+MBIG  !integer overflow
  ma(inext)=mj
  ran3=mj*FAC  !return a real number [0,1)
  return
endfunction ran3

!----------------------------------------------------------------------------------------------------
!Minimum image convention for periodic boundary conditions
!----------------------------------------------------------------------------------------------------
function image(xa,xb,ll)
  implicit none
  real*8 :: xa,xb,ll  !coordinates
  real*8 :: xc        !distance
  real*8 :: image     !function

  xc=xb-xa
  xc=xc-ll*dnint(xc/ll)  !rounds double to nearest whole double
  image=xc
  return
endfunction image

!----------------------------------------------------------------------------------------------------
!Removes bad contacts
!----------------------------------------------------------------------------------------------------
subroutine minimize()
  use parameters
  use microstate
  use macrostate
  implicit none
  integer :: i

  accx=0

  open(202,file='mc_minimization')
  call interactions()
  do i=1,10000*N
    call mcmove1()
    if (mod(i,nsampl).eq.0) write(202,'(i16,2f16.7)')i,u,p
  enddo
  close(202)
endsubroutine minimize

!----------------------------------------------------------------------------------------------------
!Calculates interactions between all particles
!----------------------------------------------------------------------------------------------------
subroutine interactions()
  use parameters
  use particles
  use microstate
  implicit none
  real*8 :: image,ljpot,tl  !functions
  real*8 :: rvec(dim_),r    !distances
  integer :: i,j,ii

  !microstate module
  pot=0.0d0
  tla=0.0d0
  dis=0.0d0
  rrr=0.0d0
  u=0.0d0
  p=0.0d0

  !interactions
  do i=1,N-1
    do j=i+1,N
      do ii=1,dim_
        rvec(ii)=image(pos(i,ii),pos(j,ii),lbox)
      enddo
      r=dsqrt(sum(rvec**2))
      if (r.lt.rcut) then
        pot(i,j)=ljpot(r,a1,a2,tl)
        tla(i,j)=tl
        dis(i,j)=r
        rrr(i,j,:)=rvec
        pot(j,i)=pot(i,j)
        tla(j,i)=tl
        dis(j,i)=r
        rrr(j,i,:)=-rvec
        u=u+pot(i,j)
        p=p+tla(i,j)
      else
        dis(i,j)=r
        rrr(i,j,:)=rvec
        dis(j,i)=r
        rrr(j,i,:)=-rvec
      endif
    enddo
  enddo
endsubroutine interactions

!----------------------------------------------------------------------------------------------------
!Calculates the LJ potential between two particles
!----------------------------------------------------------------------------------------------------
function ljpot(r,a1,a2,p)
  implicit none
  !using reduced units
  real*8 :: r
  real*8 :: a1,a2
  !a1=epsilon = absolute value of the minimum value of the LJ potential (depth of the potential well)
  !a2=sigma = distance at which the potential becomes positive
  real*8 :: x
  real*8 :: ljpot,p  !LJ potential, virial pressure

  x=a2/r
  x=x**6
  ljpot=4.0d0*a1*x*(x-1.0d0)       ! U
  p=48.0d0*a1*x*(x-0.5d0)          ! -dU/dr * r = F * r
  return
endfunction ljpot

!----------------------------------------------------------------------------------------------------
!Trial move for 1 particle
!----------------------------------------------------------------------------------------------------
subroutine mcmove1()
  use parameters
  use particles
  use microstate
  use macrostate
  implicit none
  real*8 :: ran3,image,ljpot,tl      !functions
  real*8 :: a                        !random number
  real*8 :: pot1(N),tla1(N),dis1(N)  !potential energy
  real*8 :: pos1(dim_),rvec(dim_),r  !coordinates
  real*8 :: du,dp                    !change in potential energy
  integer :: j,k,jj

  !choose random particle
  a=ran3(ran)
  k=int(a*N)+1  !floors double to integer
  
  !new position
  do jj=1,dim_
    a=ran3(ran)
    pos1(jj)=pos(k,jj)+(a-0.5d0)*delta*a2
    pos1(jj)=image(0.0d0,pos1(jj),lbox)
  enddo

  !new interactions
  do j=1,N
    if (j.eq.k) then
      pot1(j)=0.0d0
      tla1(j)=0.0d0
      dis1(j)=0.0d0
    else
      do jj=1,dim_
        rvec(jj)=image(pos1(jj),pos(j,jj),lbox)
      enddo
      r=dsqrt(sum(rvec**2))
      if (r.lt.rcut) then
        pot1(j)=ljpot(r,a1,a2,tl)
        tla1(j)=tl
        dis1(j)=r
      else
        pot1(j)=0.0d0
        tla1(j)=0.0d0
        dis1(j)=r
      endif
    endif
  enddo

  !change in potential energy
  du=0.0d0
  dp=0.0d0
  do j=1,N
    du=du+pot1(j)-pot(k,j)
    dp=dp+tla1(j)-tla(k,j)
  enddo

  !accept or reject
  a=ran3(ran)
  if(a.lt.dexp(-du/temp)) then
    pos(k,:)=pos1(:)
    do j=1,N
      pot(k,j)=pot1(j)
      tla(k,j)=tla1(j)
      dis(k,j)=dis1(j)
      pot(j,k)=pot1(j)
      tla(j,k)=tla1(j)
      dis(j,k)=dis1(j)
    enddo
    u=u+du
    p=p+dp
    accx=accx+1
  endif
endsubroutine mcmove1

!----------------------------------------------------------------------------------------------------
!MC equilibration
!----------------------------------------------------------------------------------------------------
subroutine mcequil()
  use parameters
  use microstate
  use macrostate
  implicit none
  integer :: i

  uu=0.0d0
  pp=0.0d0
  m=0
  accx=0

  open(202,file='mc_equilibration')
  call interactions()
  do i=1,nekv*N
    call mcmove1()
    if (mod(i,nsampl).eq.0) then
      write(202,'(i16,2f16.7)')i,u,p
      uu=uu+u
      pp=pp+p
      m=m+1
    endif
  enddo
  close(202)

  pp=1.0d0+pp/2.0d0/dble(N)/temp
  pp=pp*dens*temp
endsubroutine mcequil

!----------------------------------------------------------------------------------------------------
!MC sampling
!----------------------------------------------------------------------------------------------------
subroutine mcseries(idum)
  use parameters
  use microstate
  use macrostate
  use correlation
  implicit none
  integer :: idum
  character(8) :: fmt
  integer :: i

  uu=0.0d0
  uu2=0.0d0
  pp=0.0d0
  m=0
  accx=0
  gr=0.0d0
  mgr=0

  write(fmt,'(i3.3)')idum
  open(202,file='mc_sam'//trim(fmt))
  call interactions()
  do i=1,ncycl*N
    call mcmove1()
    if (mod(i,nsampl).eq.0) then
      write(202,'(i16,2f16.7)')i,u,p
      uu=uu+u
      uu2=uu2+u**2
      pp=pp+p
      m=m+1
      if (mod(i,ngr*N).eq.0) then
        call corr()
        mgr=mgr+1
      endif
    endif
  enddo
  close(202)

  pp=1.0d0+pp/2.0d0/dble(N)/temp
  pp=pp*dens*temp
  do i=1,10000
    gr(i)=gr(i)/(c*((dble(i)-0.5d0)*grinterv)**(dim_-1)*grinterv)
  enddo
  gr=gr/dble(N-1)/dens
endsubroutine mcseries

!----------------------------------------------------------------------------------------------------
!Pair correlation function g(r) -- Radial distribution function
!----------------------------------------------------------------------------------------------------
subroutine corr()
  use parameters
  use microstate
  use correlation
  implicit none
  integer :: i,j
  integer :: k
  real*8 :: r     !distance

  do i=1,N-1
    do j=i+1,N
      r=dis(i,j)
      k=int(r/grinterv)  !floors double to integer
      if (k.lt.10000) then
        gr(k+1)=gr(k+1)+2.0d0
      endif
    enddo
  enddo
endsubroutine corr

!---------------------------------------------------------------------------------------------
!Initializes positions for a homogeneous fluid
!---------------------------------------------------------------------------------------------
subroutine init_pos()
  use parameters
  use particles
  implicit none
  integer :: i,ii
  real*8 :: image         !function
  real*8 :: sumpos(dim_)  !sum of all positions

  !setting the sum of positions to zero
  sumpos=1.0d0
  do while (sum(sumpos**2).gt.1.0d-5)
    do ii=1,dim_
      sumpos(ii)=sum(pos(:,ii))/dble(N)
    enddo
    do i=1,N
      do ii=1,dim_
        pos(i,ii)=image(sumpos(ii),pos(i,ii),lbox)
      enddo
    enddo
  enddo
endsubroutine init_pos

!-------------------------------------------------------------------------------------------------------------------------------
!Initializes velocities for an isotropic fluid
!Box-Muller transform for Maxwell-Boltzmann distribution
!-------------------------------------------------------------------------------------------------------------------------------
subroutine init_vel()
  use parameters
  use particles
  implicit none
  integer :: i,ii
  real*8 :: ran3                           !function
  real*8 :: sigma,u(dim_+1),s              !sigma, random numbers
  real*8 :: sumvel(dim_)                   !sum of all velocities
  real*8 :: rescale                        !rescaling factor

  !Box-Muller transform
  sigma=dsqrt(temp)
  do i=1,N
    do ii=1,dim_,2
      s=2.0d0
      do while (s.eq.0.0d0.or.s.ge.1.0d0)
        u(ii)=2.0d0*ran3(ran)-1.0d0
        u(ii+1)=2.0d0*ran3(ran)-1.0d0
        s=u(ii)**2+u(ii+1)**2
      enddo
      u(ii)=sigma*u(ii)*dsqrt(-2.0d0*dlog(s)/s)
      u(ii+1)=sigma*u(ii+1)*dsqrt(-2.0d0*dlog(s)/s)
    enddo
    vel(i,:)=u(:dim_)
  enddo

  !setting the sum of velocities to zero
  sumvel=1.0d0
  do while (sum(sumvel**2).gt.1.0d-5)
    do ii=1,dim_
      sumvel(ii)=sum(vel(:,ii))/dble(N)
      vel(:,ii)=vel(:,ii)-sumvel(ii)
    enddo
  enddo

  !calculate the temperature
  call temperature()

  !rescale the particle velocities to the desired temperature
  rescale=dsqrt(temp/tem)
  vel=vel*rescale
endsubroutine init_vel

!----------------------------------------------------------------------------------------------------
!Calculates current temperature using the equipartition theorem
!----------------------------------------------------------------------------------------------------
subroutine temperature()
  use parameters
  use particles
  implicit none
  integer :: i

  ekin=0.0d0  
  do i=1,N
    ekin=ekin+(sum(vel(i,:)**2))/2.0d0
  enddo
  !(x,y,z) -> D degrees of freedom
  !1/2 T for each degree of freedom
  !N*D degrees of freedom (N particles with (x,y,z))
  tem=2.0d0*ekin/dble(N*dim_)
endsubroutine temperature

!----------------------------------------------------------------------------------------------------
!MD equilibration
!----------------------------------------------------------------------------------------------------
subroutine mdequil()
  use parameters
  use particles
  use microstate
  use macrostate
  implicit none
  integer :: i
  real*8 :: t

  t=0.0d0
  uu=0.0d0
  pp=0.0d0
  m=0

  open(202,file='md_equilibration')
  call forces()
  call temperature()
  do i=1,nekv
    call berendsen()
    t=t+tstep
    etot=ekin+u
    write(202,'(6f16.7)')t,u,ekin,etot,p,tem
  enddo
  call forces()
  call temperature()
  do i=1,nekv
    call nose_hoover_chain()
    t=t+tstep
    etot=ekin+u
    write(202,'(6f16.7)')t,u,ekin,etot,p,tem
    uu=uu+u
    pp=pp+p
    m=m+1
  enddo
  close(202)

  pp=1.0d0+pp/2.0d0/dble(N)/temp
  pp=pp*dens*temp
endsubroutine mdequil

!----------------------------------------------------------------------
!Calculates forces from interactions
!----------------------------------------------------------------------
subroutine forces()
  use parameters
  use particles
  use microstate
  implicit none
  integer :: i,j,ii
  real*8 :: f(N,N,dim_)

  call interactions()
  f=0.0d0
  do i=1,N-1
    do j=i+1,N
      do ii=1,dim_
        f(i,j,ii)=tla(i,j)*rrr(i,j,ii)/dis(i,j)**2  
        f(j,i,ii)=-f(i,j,ii)
      enddo
    enddo
  enddo
  do j=1,N
    do ii=1,dim_
      acc(j,ii)=sum(f(:,j,ii))
    enddo
  enddo
endsubroutine forces

!-------------------------------------------------------------------------
!Velocity Verlet algorithm for calculating velocities
!Berendsen rescaling thermostat with temperature time constant tt
!Newton's equations of motion
!-------------------------------------------------------------------------
subroutine berendsen()
  use parameters
  use particles
  implicit none
  integer :: i,ii
  real*8 :: image    !function
  real*8 :: rescale  !rescaling factor

  do i=1,N
    do ii=1,dim_
      !half-update the particle velocities by tstep/2
      vel(i,ii)=vel(i,ii)+acc(i,ii)*dt2
      !update the particle positions by tstep
      pos(i,ii)=pos(i,ii)+vel(i,ii)*tstep
      pos(i,ii)=image(0.0d0,pos(i,ii),lbox)
    enddo
  enddo

  !update the forces
  call forces()

  do i=1,N
    do ii=1,dim_
      !update the particle velocities by tstep/2
      vel(i,ii)=vel(i,ii)+acc(i,ii)*dt2
    enddo
  enddo
  
  !update the temperature
  call temperature()

  !rescale the particle velocities to the desired temperature
  rescale=dsqrt(1.0d0+tstep/tt*(temp/tem-1.0d0))
  vel=vel*rescale
endsubroutine berendsen

!-------------------------------------------------------------------------------------------------------
!Velocity Verlet algorithm for calculating velocities
!Nose-Hoover chain thermostat with 10 heat bath coupling massess Q
!Extended Lagrangian mechanics -- Martyna-Tuckerman-Tobias-Klein equations of motion
!-------------------------------------------------------------------------------------------------------
subroutine nose_hoover_chain()
  use parameters
  use particles
  implicit none
  integer :: i,ii
  real*8 :: image  !function

  !half-update the thermostat
  call nhc_thermostat()
  !half-update the particle velocities
  do i=1,N
    do ii=1,dim_
      vel(i,ii)=vel(i,ii)+acc(i,ii)*dt2
    enddo
  enddo
  !update the particle positions
  do i=1,N
    do ii=1,dim_
      pos(i,ii)=pos(i,ii)+vel(i,ii)*tstep
      pos(i,ii)=image(0.0d0,pos(i,ii),lbox)
    enddo
  enddo
  !update the forces
  call forces()
  !update the particle velocities
  do i=1,N
    do ii=1,dim_
      vel(i,ii)=vel(i,ii)+acc(i,ii)*dt2
    enddo
  enddo
  !update the thermostat
  call temperature()
  call nhc_thermostat()
endsubroutine nose_hoover_chain

!-------------------------------------------------------------------------------------------------------
!Nose-Hoover chain thermostat with 10 heat bath coupling massess Q
!-------------------------------------------------------------------------------------------------------
subroutine nhc_thermostat()
  use parameters
  use particles
  implicit none
  integer :: i,ii
  real*8 :: X
  real*8 :: AA

  G(1)=(2.0d0*ekin-NfkT)/Q(1)
  vksi(10)=vksi(10)+G(10)*dt4
  do i=1,9
    X=dexp(-vksi(11-i)*dt8)
    vksi(10-i)=vksi(10-i)*X
    vksi(10-i)=vksi(10-i)+G(10-i)*dt4
    vksi(10-i)=vksi(10-i)*X
  enddo

  AA=dexp(-vksi(1)*dt2)
  !scale the particle velocities
  do i=1,N
    do ii=1,dim_
      vel(i,ii)=vel(i,ii)*AA
    enddo
  enddo
  !update the temperature
  ekin=ekin*AA**2
  tem=tem*AA**2
  
  G(1)=(2.0d0*ekin-NfkT)/Q(1)
  do i=1,9
    X=dexp(-vksi(i+1)*dt8)
    vksi(i)=vksi(i)*X
    vksi(i)=vksi(i)+G(i)*dt4
    vksi(i)=vksi(i)*X
    G(i+1)=(Q(i)*vksi(i)**2-temp)/Q(i+1)
  enddo
  vksi(10)=vksi(10)+G(10)*dt4
endsubroutine nhc_thermostat

!----------------------------------------------------------------------------------------------------
!MD sampling
!----------------------------------------------------------------------------------------------------
subroutine mdseries(idum)
  use parameters
  use particles
  use microstate
  use macrostate
  use correlation
  implicit none
  integer :: idum
  character(8) :: fmt
  integer :: i
  real*8 :: t

  t=0.0d0
  uu=0.0d0
  uu2=0.0d0
  pp=0.0d0
  m=0
  gr=0.0d0
  mgr=0

  write(fmt,'(i3.3)')idum
  open(202,file='md_sam'//trim(fmt))
  call forces()
  call temperature()
  do i=1,nsteps
    call nose_hoover_chain()
    t=t+tstep
    etot=ekin+u
    write(202,'(6f16.7)')t,u,ekin,etot,p,tem
    uu=uu+u
    uu2=uu2+u**2
    pp=pp+p
    m=m+1
    if (mod(i,ngr).eq.0) then
      call corr()
      mgr=mgr+1
    endif
    if (idum.eq.1.and.mod(i,nframe).eq.0.and.(i/nframe).lt.101) then
      write(fmt,'(i2.2,i6.6)')idum,i/nframe
      call frame(fmt)
    endif
  enddo
  close(202)

  pp=1.0d0+pp/2.0d0/dble(N)/temp
  pp=pp*dens*temp
  do i=1,10000
    gr(i)=gr(i)/(c*((dble(i)-0.5d0)*grinterv)**(dim_-1)*grinterv)
  enddo
  gr=gr/dble(N-1)/dens
endsubroutine mdseries

!-------------------------------------------------------
!Snapshot
!-------------------------------------------------------
subroutine snapshot(fmt)
  use parameters
  use particles
  implicit none
  character(6) :: fmt
  integer :: i,ii

  open(111,file='snapshot_'//fmt)
  do i=1,N
    do ii=1,dim_
      write(111,'(3f16.7)',advance='no')pos(i,ii)/lbox
    enddo
    write(111,'(3f16.7)')a2/2.0d0/lbox
  enddo
  close(111)
endsubroutine snapshot

!-------------------------------------------------------
!Velocities
!-------------------------------------------------------
subroutine velocities(fmt)
  use parameters
  use particles
  implicit none
  character(4) :: fmt
  integer :: i,ii

  open(111,file='md_velocities'//fmt)
  do i=1,N
    do ii=1,dim_
      write(111,'(3f16.7)',advance='no')vel(i,ii)
    enddo
    write(111,*)''
  enddo
  close(111)
endsubroutine velocities

!--------------------------------------------------------
!Movie
!--------------------------------------------------------
subroutine frame(fmt)
  use parameters
  use particles
  implicit none
  character(8) :: fmt
  integer :: i,ii

  open(111,file='frame_'//fmt)
  do i=1,N
    do ii=1,dim_
      write(111,'(3f16.7)',advance='no')pos(i,ii)/lbox
    enddo
    write(111,'(3f16.7)')a2/2.0d0/lbox
  enddo
  close(111)
endsubroutine frame
