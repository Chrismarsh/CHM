!=================================================!
! Distributed Blowing Snow Model                  !
! Richard Essery, University of Wales Aberystwyth !
! 17 July 2006 modified 26 July 2006 TB.                                    !
!=================================================!
      implicit none

      integer
     & n                 ! Grid dimension
     &,nsteps            ! Number of timesteps in met file
     &,sx,sy             ! Grid location of met station
      real 
     & dt                ! Timestep (s)
     &,dx                ! Grid spacing (m)
      character*70         !(len=70)
     & met_file          ! Meteorology file name
     &,topo_file         ! DEM file name
     &,veg_file          ! Vegetation file name
      include 'site_params.inc'
      real     ! Met data 
     & temp              ! Air temperature (C)
     &,rh                ! Relative humidity (fraction)
     &,dir               ! Wind direction (degrees)
     &,Us                ! Windspeed (m/s)
     &,snow              ! Snowfall rate (kg/m2/s)
      real 
     & p(n,n)            ! Probability of blowing snow occurence 
     &,rho(n,n)          ! Snow density (kg/m3)
     &,s(n,n)            ! SWE (kg/m2)
     &,sublim(n,n)       ! Sublimation (kg/m2/s)
     &,trans(n,n)        ! Transport (kg/m/s)
     &,U(n,n,8)          ! Normalized westerly wind component
     &,V(n,n,8)          ! Normalized southerly wind component
     &,W(n,n,8)          ! Normalized wind speed
     &,Ur(n,n)           ! Rotated wind component
     &,Vr(n,n)           ! Rotated wind component
     &,veg_ht(n,n)       ! Vegetation height (m)
     &,wet_snow(n,n)     ! Wet snow mass (kg/m2) 
     &,z(n,n)            ! Elevation (m)
     &,zr(n,n)           ! Rotated topography (m)
      logical 
     & dry_snow          ! Dry / wet snow indicator
      real 
     & r                 ! Rotation angle (radians)
     &,denst
     &,rho_new           ! Fresh snow density (kg/m3)
     &,snow_age          ! Snow age (hours)
     &,sublimation       ! Sublimation windspeed function
     &,subscale          ! Sublimation temperature function
     &,transport         ! Transport function
     &,Wave              ! Average windspeed
      real     ! Transport across cell boundaries (kg/m/s) 
     & tw                ! - west boundary
     &,te                ! - east boundary
     &,ts                ! - south boundary
     &,tn                ! - north boundary
     &,Tin               ! Transport into cell 
     &,Tout              ! Transport out of cell
      integer
     & d                 ! Discretized wind direction
!                        ! 1 = NE, 2 = E, ..., 8 = N
     &,dd(8)             ! Corresponding grid rotations
     &,i                 ! W - E loop counter
     &,j                 ! S - N loop counter
     &,k                 ! Wind direction loop counter
     &,step              ! Timestep counter
      data dd / 6, 5, 4, 3, 2, 1, 8, 7 /
      integer year
	  real day
      real hour

        character*70 base
         character*70 strdir
         character*70 ext
         character*70 fname
         integer rr

! Read DEM and vegetation map
      open_GEM(21,file=topo_file)
      do j=1,6
!        !skip
        read(21,*)
      enddo
      do j=1,n
        read(21,*) (z(i,j),i=1,n)
      enddo
      close(21)
      open_GEM(21,file=veg_file)
      do j=1,6
             !skip
           read(21,*)
      enddo
      do j=1,n
        read(21,*) (veg_ht(i,j),i=1,n)
      enddo
      close(21)

! Blend topography into plane
      call BLEND(n,z,.true.)

! Initialize snow prognostics
      dry_snow = .true.
      snow_age = 0.
      do j=1,n 
      do i=1,n
        rho(i,j) = 100.
        S(i,j) = 0.
        wet_snow(i,j) = 0.
      enddo
      enddo

! Loop over wind directions
      do k=1,8
        r = 45*(k-1) *3.14159/180.
        d = dd(k)
        rr = d
        base= 'Windspeed_Normalized_'
        ext = '.txt'
        write(strdir,'(I0)') rr
        fname = trim(base) // trim(strdir) // trim(ext)

        open_GEM(32,file=fname)

! Rotate topography to wind direction
        call ROTATE(n,z,r,zr)

! Run windflow model
        call WIND(n,dx,zr,U(1,1,d),V(1,1,d))
        call BLEND(n,U(1,1,d),.false.)
        call BLEND(n,V(1,1,d),.false.)

! Rotate wind components back and normalize
        r = -r
        call ROTATE(n,U(1,1,d),r,Ur)
        call ROTATE(n,V(1,1,d),r,Vr)
        Wave = 0.
        do j=1,n
            do i=1,n
              U(i,j,d) = Ur(i,j)*cos(r) + Vr(i,j)*sin(r)
              V(i,j,d) = Vr(i,j)*cos(r) - Ur(i,j)*sin(r)
              W(i,j,d) = sqrt( U(i,j,d)**2 + V(i,j,d)**2 )
              Wave = Wave + W(i,j,d)
            enddo
        enddo
        Wave = Wave/(n*n) 
        do j=1,n
            do i=1,n
              U(i,j,d) = U(i,j,d) / W(i,j,d)
              V(i,j,d) = V(i,j,d) / W(i,j,d)
              W(i,j,d) = W(i,j,d) / Wave
            enddo
        enddo


         do j=1,n
            write(32, '(1f15.6)') (W(i,j,d),i=1,n)
         enddo
         close(32)


        enddo ! End of windflow calculations

! Open met file and start loop over timesteps
      open_GEM(21,file=met_file)
      do step=1,nsteps
        if (mod(step,48).eq.0) 
     &     print 100, 100*step/nsteps
 100	format(I3,'% complete') 
! Read met data and discretize wind direction
        read(21,*) year,day,temp,rh,Us,dir,snow
        rho_new = 67.9 + 51.3*exp(temp/2.6)
        denst = 14.643 - 4000/(273.15 + min(temp,0.)) 
        d = nint(dir/45.)
        if (d.eq.0) d = 8
        Us = Us / W(sx,sy,d)

! Increment snow age, SWE and wet snow mass
        do j=1,n
        do i=1,n 
          wet_snow(i,j) = min( wet_snow(i,j), S(i,j) )
        enddo
        enddo
        if (snow .gt. 0.) then
          snow_age = 1.
          dry_snow = .true.
          do j=1,n
          do i=1,n
            if (s(i,j).gt.0) then
            rho(i,j) = (rho(i,j)*S(i,j) + rho_new*snow*dt) /
     &                 (S(i,j) + snow*dt)
            else
            rho(i,j) = rho_new
            endif
            S(i,j) = S(i,j) + snow*dt
          enddo
          enddo
        else
          snow_age = snow_age + dt / 3600.
        endif 
        if (temp.ge.0.) then
          dry_snow = .false.
          do j=1,n
          do i=1,n
            wet_snow(i,j) = S(i,j) 
          enddo
          enddo
        endif   

! Calulate probabilities of blowing snow occurence in cells
        call PROB(n,snow_age,temp,Us,dry_snow,rho,S,
     &            W(1,1,d),wet_snow,veg_ht,p)

! Calculate transport and sublimation fluxes for fully developed
! blowing snow
        call FLUX(n,p,rh,temp,Us,W(1,1,d),sublim,trans)

! Adjust fluxes for development
        call DEVELOP(n,d,dx,u(1,1,d),v(1,1,d),sublim,trans)

! Increment snow water equivalent and density
        do j=2,n-1
        do i=2,n-1
          if (S(i,j) .le. 0) then
            trans(i,j) = 0.
            sublim(i,j) = 0.
          endif
          tw = trans(i-1,j) * u(i-1,j,d)
          te = trans(i+1,j) * u(i+1,j,d) 
          ts = trans(i,j-1) * v(i,j-1,d) 
          tn = trans(i,j+1) * v(i,j+1,d)
          Tin = max(Tw, 0.) - min(Te, 0.) + max(Ts, 0.) - min(Tn, 0.)
          Tout = trans(i,j) * (abs(u(i,j,d)) + abs(v(i,j,d))) 

          S(i,j) = S(i,j) + dt*( (Tin - Tout)/dx - sublim(i,j) )
          rho(i,j) = rho(i,j) * ( 1 + 0.5e-8*9.81*S(i,j)*dt *
     &                                exp(denst - 0.02*rho(i,j)) )
        enddo
        enddo

! Fill in boundaries
        do i=1,n
          rho(i,1) = rho(i,2)
          rho(i,n) = rho(i,n-1)
          rho(1,i) = rho(2,i)
          rho(n,i) = rho(n-1,i)
          s(i,1) = s(i,2)
          s(i,n) = s(i,n-1)
          s(1,i) = s(2,i)
          s(n,i) = s(n-1,i)
        enddo

      enddo ! End of blowing snow calculations
      close(21)

! Write final SWE and depth grids
      open_GEM(31,file='SWE.txt')
      open_GEM(33,file='depth.txt')
      do j=1,n
        write(31, 101) (S(i,j),i=1,n)
 101  format(5f15.6)
        write(33,102) (S(i,j)/rho(i,j),i=1,n)
 102  format(5f15.6)
      enddo
      close(31)
      close(33)
      
      end

!=================================
      subroutine BLEND(n,z,smooth)
! Blend topography into plane
!=================================
      implicit none
      integer 
     & n                 ! IN Grid dimension
      real
     & z(n,n)            ! INOUT Elevations
      logical smooth     ! IN Set .TRUE. to apply smoothing

      real
     & b                 ! Blending weight
     &,f                 ! Fraction of grid not blended
     &,r                 ! Distance from grid centre
     &,w                 ! Blending width
     &,zm                ! Average elevation
      real, dimension(:,:), allocatable ::  zs !(:,:)           ! Smoothed elevations
      parameter(f=0.6, w=0.316)
      integer
     & i,j,is,js         ! Loop counters
	 
      allocate (zs(n,n))
! Calculate average elevation
      do j=1,n
      do i=1,n
        zm = zm + z(i,j)        
      enddo
      enddo
      zm = zm / (n*n)

! Blend topography into plane
      do j=1,n
      do i=1,n
        r = sqrt((i-n/2.)**2+(j-n/2.)**2) / (0.5*f*n)
        b = 1
        if (r .gt. 1) b = exp(-((r-1)/w)**2)
        z(i,j) = b*z(i,j) + (1 - b)*zm
      enddo
      enddo

! Smooth topography
      if (smooth) then
      do j=3,n-2
      do i=3,n-2
        zs(i,j) = 0.
        do js=j-2,j+2
        do is=i-2,i+2
          zs(i,j) = zs(i,j) + z(is,js)
        enddo
        enddo
      enddo
      enddo
      do j=3,n-2
      do i=3,n-2
        z(i,j) = 0.04*zs(i,j)
      enddo
      enddo
      endif

      open_GEM(31,file='DEMsmooth.txt')
      do j=1,n
        write(31,109) (z(i,j),i=1,n)
109   format(1f15.6)
      enddo
      close(31)

      deallocate (zs)
      return
      end

!================================
      subroutine ROTATE(n,g,r,gr)
! Rotate grids
!================================
      implicit none

      integer 
     & n                 ! IN Grid dimension
      real
     & g(n,n)            ! IN Grid
     &,r                 ! IN Rotation angle (radians)
     &,gr(n,n)           ! OUT Rotated grid

      real
     & x,y,dx,dy         ! Grid coordinates
     &,xr,yr             ! Rotated grid coordinates
      integer
     & i,j               ! Grid indices
     &,ir,jr             ! Rotated grid indices

      do ir=1,n
      do jr=1,n
        xr = ir - (n+1)/2
        yr = jr - (n+1)/2
        x = xr*cos(r) - yr*sin(r)
        y = xr*sin(r) + yr*cos(r)
        i = x + (n+1)/2
        dx = x + (n+1)/2 - i
        j = y + (n+1)/2
        dy = y + (n+1)/2 - j
        if (i .gt. 0 .and. i .lt. n-1 .and.
     &      j .gt. 0 .and. j .lt. n-1) then 
           gr(ir,jr) = g(i,j)*(dx - 1)*(dy - 1) - g(i+1,j)*dx*(dy - 1)
     &              - g(i,j+1)*(dx - 1)*dy + g(i+1,j+1)*dx*dy
        else
          gr(ir,jr) = g(1,1)
        endif
      enddo                                                                 
      enddo

      return
      end

!================================
      subroutine WIND(n,dx,z,u,v)
! Mason-Sykes windflow model
!================================
      implicit none

      integer 
     & n                 ! IN Grid dimension
      real
     & dx                ! IN Grid spacing
     &,z(n,n)            ! IN Elevation map
     &,u(n,n)            ! OUT Westerly wind component
     &,v(n,n)            ! OUT Southerly wind component

!      real
!     & ur(0:n-1,0:n-1)   ! Real part of U and Fourier transform
!     &,ui(0:n-1,0:n-1)   ! Imag. part of U and Fourier transform
!     &,vr(0:n-1,0:n-1)   ! Real part of V and Fourier transform
!     &,vi(0:n-1,0:n-1)   ! Imag. part of V and Fourier transform
!     &,zr(0:n-1,0:n-1)   ! Real part of Z and Fourier transform
!     &,zi(0:n-1,0:n-1)   ! Imag. part of Z and Fourier transform
!     &,rchiu(n/2),ichiu(n/2)
!     &,rchiv(n/2),ichiv(n/2)
	 
	  real, dimension(:,:), allocatable ::  ur, ui, vr, vi, zr, zi

      real, dimension(:), allocatable :: rchiu, ichiu
      real, dimension(:), allocatable ::  rchiv, ichiv
	 
      real 
     & epsil
     &,l                 ! Inner length scale (m)
     &,L0                ! Horizontal length scale (m)
     &,k,m               ! Wavenumbers
     &,wfac              ! Scaled wavenumber
     &,pi2ndx            ! 2*pi / (n*dx)
     &,z0                ! Surface roughness length (m)
     &,z1                ! Windspeed measurement height (m)
       
      parameter ( L0=1e3, z0=1e-2, z1=2. )
	  
      complex 
     & y                 ! Ratio of Bessel functions
     &,K0                ! Modified Bessel function
      integer 
     & i,j                ! Loop counters 

	  allocate(ur(0:n-1,0:n-1),ui(0:n-1,0:n-1))
	  allocate(vr(0:n-1,0:n-1),vi(0:n-1,0:n-1))
	  allocate(zr(0:n-1,0:n-1),zi(0:n-1,0:n-1))
      allocate(rchiu(n/2), ichiu(n/2))
      allocate(rchiv(n/2), ichiv(n/2))
	 
      l = (z0/8.)*(L0/z0)**0.9
      do i=1,100
        l = 2*0.16*L0/log(l/z0)
      enddo
      epsil = 0.66*(log(L0/z0))**2/(log(l/z0)*log(z1/z0))

      pi2ndx = 2*3.14159/(n*dx) 
      do i=1,n/2
        k = pi2ndx*i
        y = 1 - K0(2*sqrt(L0*k*z1/l)) / K0(2*sqrt(L0*k*z0/l))
        rchiu(i) = real(y)
        ichiu(i) = - real((0.,1.)*y)
        y = 1 - K0(2*sqrt(2*L0*k*z1/l)) / K0(2*sqrt(2*L0*k*z0/l))
        rchiv(i) = real(y)
        ichiv(i) = - real((0.,1.)*y)
      enddo

      do j=0,n-1
      do i=0,n-1
        zi(i,j) = 0.
        zr(i,j) = z(i+1,j+1)
      enddo
      enddo

! Fourier transform of elevation
      call FFT2(zr,zi,n,1)

      do i=1,n/2 
        do j=0,n/2     
          k = pi2ndx*i    
          m = pi2ndx*j
          wfac = k**2/sqrt(k**2+m**2)
          ur(i,j) = (rchiu(i)*zr(i,j) + ichiu(i)*zi(i,j))*wfac
          ui(i,j) = (-ichiu(i)*zr(i,j) + rchiu(i)*zi(i,j))*wfac
          wfac = k*m/sqrt(k**2+m**2)
          vr(i,j) = (rchiv(i)*zr(i,j) + ichiv(i)*zi(i,j))*wfac
          vi(i,j) = (-ichiv(i)*zr(i,j) + rchiv(i)*zi(i,j))*wfac
        enddo
        do j=n/2+1,n-1      
          k = pi2ndx*i     
          m = pi2ndx*(j-n)
          wfac = k**2/sqrt(k**2+m**2)
          ur(i,j) = (rchiu(i)*zr(i,j) + ichiu(i)*zi(i,j))*wfac
          ui(i,j) = (-ichiu(i)*zr(i,j) + rchiu(i)*zi(i,j))*wfac
          wfac = k*m/sqrt(k**2+m**2)
          vr(i,j) = (rchiv(i)*zr(i,j) + ichiv(i)*zi(i,j))*wfac
          vi(i,j) = (-ichiv(i)*zr(i,j) + rchiv(i)*zi(i,j))*wfac
        enddo
      enddo
      do i=n/2+1,n-1
        ur(i,0) = ur(n-i,0)
        ui(i,0) = - ui(n-i,0)
      enddo 
      do i=n/2+1,n-1
        do j=1,n-1
          ur(i,j) = ur(n-i,n-j)
          ui(i,j) = - ui(n-i,n-j)
          vr(i,j) = vr(n-i,n-j)
          vi(i,j) = - vi(n-i,n-j)
        enddo
      enddo 

! Inverse Fourier transform of wind components
      call FFT2(ur,ui,n,-1)
      call FFT2(vr,vi,n,-1)

      do j=1,n     
      do i=1,n 
        U(i,j) = 1. + epsil*ur(i-1,j-1)  
        V(i,j) = epsil*vr(i-1,j-1) 
        if (U(i,j) .lt. 0.01) U(i,j) = 0.01
      enddo
      enddo

	  deallocate(ur,ui)
	  deallocate(vr,vi)
	  deallocate(zr,zi)
      deallocate(rchiu, ichiu)
      deallocate(rchiv, ichiv)
	  
      return
      end

!=====================================================
      subroutine FLUX(n,p,rh,temp,Us,W,sublim,trans)
! Calculate transport and sublimation rates for 
! fully-developed blowing snow
!=====================================================
      implicit none
      integer
     & n                 ! IN Grid dimension
     &,d                 ! IN Discretized wind direction
      real 
     & p(n,n)            ! IN Probability of blowing snow occurence 
     &,rh                ! IN Relative humidity (%)
     &,temp              ! IN Air temperature (C)
     &,Us                ! IN Windspeed at met station (m/s)
     &,W(n,n)            ! IN Normalized windspeed
     &,sublim(n,n)       ! OUT Sublimation (kg/m2/s)
     &,trans(n,n)        ! OUT Transport (kg/m/s)

      real
     & cond              ! Thermal conductivity of air (W/m/K)
     &,diff              ! Diffusivity of water vapour in air (m2/s)
     &,Ls                ! Latent heat of sublimation (J/kg)
     &,m                 ! Molecular weight of water (kg/kmole)
     &,R                 ! Universal gas constant (J/kmole/K)
     &,rsat              ! Saturation density of water vapour (kg/m3)
     &,subscale          ! Sublimation scaling function
     &,trnscale
     &,Tk                ! Temperature (K)
     &,Ul                ! Local windspeed (m/s)
      parameter ( Ls=2.838e6, m=18.01, R=8313. )
      integer
     & i,j               ! Loop counters

      Tk = temp + 273.
      diff = 2.06e-5*(Tk/273.)**1.75
      rsat = m*611.15*exp(22.45*temp/Tk)/(R*Tk)
      cond = (0.00063*Tk + 0.0673) / 10.
      subscale = ((Ls*m/(R*Tk)) - 1.)/(cond*(temp+273.))
     &           + 1 /(Ls*diff*rsat)
      subscale = 0.1376*(1 - rh/100.) / subscale
      trnscale = (0.00096*temp**2 + 0.5298*temp + 666.82) * 1e-3
      Us = Us/25.
      do j=1,n
      do i=1,n
        Ul = Us*W(i,j)
        sublim(i,j) = p(i,j)*subscale*Ul**5
        trans(i,j) = p(i,j)*trnscale*Ul**4
      enddo
      enddo

! Write normalized windspeed, sublimation and transport grids
!      open_GEM(32,file='Windspeed_Normalized.txt')
!      open_GEM(34,file='Sublimation.txt')
!      open_GEM(36,file='Transport.txt')
!      do j=1,n
!        write(32, 103) (W(i,j),i=1,n)
!103   format(1f15.6)
!        write(34, 104) (sublim(i,j),i=1,n)
!104   format(5f15.6)
!        write(36, 105) (trans(i,j),i=1,n)
!105   format(5f15.6)
!      enddo
!      close(32)
!      close(34)
!      close(36)

      return
      end

!=========================================================
      subroutine PROB(n,snow_age,temp,Us,dry_snow,rho,S,W,
     &                wet_snow,veg_ht,p)
! Calculate probability of blowing snow occurence
!========================================================= 
      implicit none
      integer
     & n                 ! IN Grid dimension
     &,d                 ! IN Discretized wind direction
      real 
     & snow_age          ! IN Snow age (hours)
     &,temp              ! IN Air temperature (C)
     &,Us                ! IN Windspeed at met station (m/s)
      logical 
     & dry_snow          ! IN Dry snow indicator
      real
     & rho(n,n)          ! IN Snow density (kg/m3)
     &,S(n,n)            ! IN SWE (kg/m2)
     &,W(n,n)            ! IN Normalized windspeed
     &,wet_snow(n,n)     ! IN Wet snow mass (kg/m2)
     &,veg_ht(n,n)       ! IN Vegetation height (m)
     &,p(n,n)            ! OUT Probability of blowing snow occurence 
      real
     & mean              ! Mean of cummulative normal distribution 
     &,var               ! Standard deviation
     &,sd                ! Snow depth (m)
     &,Ul                ! Local windspeed (m/s)
      integer
     & i,j               ! Loop counters

      mean = 0.365*temp + 0.00706*temp**2 + 0.91*log(snow_age) + 11.0
      var = 0.145*temp + 0.00196*temp**2 + 4.23
      if (.not. dry_snow) then
        mean = 21.
        var = 7.
      endif

      do j=1,n
      do i=1,n
        sd = S(i,j) / rho(i,j)
        Ul = Us*W(i,j)
        if (sd. lt. veg_ht(i,j))
     &    Ul = Ul / sqrt(1. + 170.*2*0.1*(veg_ht(i,j) - sd))
        p(i,j) = 1 / ( 1. + exp(1.7*(mean - Ul)/var) )
        if (S(i,j) .le. wet_snow(i,j)) then
          p(i,j) = 1 / ( 1. + exp(1.7*(21. - Ul)/7.) )
          if (Ul.le.7.) p(i,j) = 0.
        endif
        if (sd .le. 0.01) p(i,j) = 0.
        if ( dry_snow ) then
          if (Ul .le. 3.) p(i,j) = 0.
        else
          if (Ul .le. 7.) p(i,j) = 0.
        endif
      enddo
      enddo

      return
      end

!================================================
      subroutine DEVELOP(n,d,dx,U,V,sublim,trans)
! Adjust sublimation and transport fluxes for
! flow development
!================================================
      implicit none
      integer
     & n                 ! IN Grid dimension
     &,d                 ! IN Discretized wind direction
      real
     & dx                ! IN Grid spacing (m)
     &,U(n,n)            ! IN Normalized westerly wind component
     &,V(n,n)            ! IN Normalized southerly wind component
     &,sublim(n,n)       ! INOUT Sublimation (kg/m2/s)
     &,trans(n,n)        ! INOUT Transport (kg/m/s)

      real 
     & a,a1,a2           ! Flux weighting factors
     &,f1                ! Increasing flow fetch requirement (m)
     &,f2                ! Decaying flow fetch requirement (m)
     &,tw,te,ts,tn       ! Transport across cell boundaries
  !   &,Tin(:,:)          ! Transport into cell

      real, dimension(:,:), allocatable :: Tin
      parameter( f1=333., f2=10. )
      integer
     & i,j               ! Loop counters

      allocate(Tin(n,n))
	  
      do j=2,n-1
      do i=2,n-1
        tw = trans(i-1,j) * u(i-1,j)
        te = trans(i+1,j) * u(i+1,j) 
        ts = trans(i,j-1) * v(i,j-1) 
        tn = trans(i,j+1) * v(i,j+1)
        Tin(i,j) = sqrt( (max(Tw, 0.) - min(Te, 0.))**2 + 
     &                   (max(Ts, 0.) - min(Tn, 0.))**2 )
      enddo
      enddo
      do i=1,n
        Tin(i,1) = Trans(i,1)
        Tin(i,n) = Trans(i,n)
        Tin(1,i) = Trans(1,i)
        Tin(n,i) = Trans(n,i)
      enddo

! Cardinal wind directions
      a1 = 1 / (1 + f1/dx)
      a2 = 1 / (1 + f2/dx)
      if (d.eq.2) then
        do j=1,n
        do i=n-1,1,-1
          a = a1
          if (trans(i,j) .lt. Tin(i,j)) a = a2
          trans(i,j) = a*trans(i,j) + (1 - a)*trans(i+1,j)
          sublim(i,j) = a*sublim(i,j) + (1 - a)*sublim(i+1,j)
        enddo
        enddo
        deallocate(Tin)
        return
      elseif (d.eq.4) then
        do j=2,n
        do i=1,n
          a = a1
          if (trans(i,j) .lt. Tin(i,j)) a = a2
          trans(i,j) = a*trans(i,j) + (1 - a)*trans(i,j-1)
          sublim(i,j) = a*sublim(i,j) + (1 - a)*sublim(i,j-1)
        enddo
        enddo
        deallocate(Tin)
        return
      elseif (d.eq.6) then
        do j=1,n
        do i=2,n
          a = a1
          if (trans(i,j) .lt. Tin(i,j)) a = a2
          trans(i,j) = a*trans(i,j) + (1 - a)*trans(i-1,j)
          sublim(i,j) = a*sublim(i,j) + (1 - a)*sublim(i-1,j) 
        enddo
        enddo
      elseif (d.eq.8) then
        do i=1,n
        do j=n-1,1,-1
          a = a1
          if (trans(i,j) .lt. Tin(i,j)) a = a2
          trans(i,j) = a*trans(i,j) + (1 - a)*trans(i,j+1)
          sublim(i,j) = a*sublim(i,j) + (1 - a)*sublim(i,j+1) 
        enddo
        enddo
		
        deallocate(Tin)
        return
      endif

! Diagonal wind directions
      a1 = 1 / (1 + f1/(sqrt(2.)*dx))
      a2 = 1 / (1 + f2/(sqrt(2.)*dx))
      if (d.eq.1) then
        do j=n-1,1,-1
        do i=n-1,1,-1
          a = a1
          if (trans(i,j) .lt. Tin(i,j)) a = a2
          trans(i,j) = a*trans(i,j) + 0.5*(1-a)*trans(i+1,j)
     &                              + 0.5*(1-a)*trans(i,j+1) 
          sublim(i,j) = a*sublim(i,j) + 0.5*(1-a)*sublim(i+1,j)
     &                                + 0.5*(1-a)*sublim(i,j+1)
        enddo
        enddo
        deallocate(Tin)
        return
      elseif (d.eq.3) then
        do j=2,n
        do i=n-1,1,-1
          a = a1
          if (trans(i,j) .lt. Tin(i,j)) a = a2
          trans(i,j) = a*trans(i,j) + 0.5*(1-a)*trans(i+1,j)
     &                              + 0.5*(1-a)*trans(i,j-1)
          sublim(i,j) = a*sublim(i,j) + 0.5*(1-a)*sublim(i+1,j)
     &                                + 0.5*(1-a)*sublim(i,j-1)
        enddo
        enddo
        deallocate(Tin)
        return
      elseif (d.eq.5) then
        do j=2,n
        do i=2,n
          a = a1
          if (trans(i,j) .lt. Tin(i,j)) a = a2
          trans(i,j) = a*trans(i,j) + 0.5*(1-a)*trans(i-1,j)
     &                              + 0.5*(1-a)*trans(i,j-1)
          sublim(i,j) = a*sublim(i,j) + 0.5*(1-a)*sublim(i-1,j)
     &                                + 0.5*(1-a)*sublim(i,j-1) 
        enddo
        enddo
        deallocate(Tin)
        return
      elseif (d.eq.7) then
        do j=n-1,1,-1
        do i=2,n
          a = a1
          if (trans(i,j) .lt. Tin(i,j)) a = a2
          trans(i,j) = a*trans(i,j) + 0.5*(1-a)*trans(i-1,j)
     &                              + 0.5*(1-a)*trans(i,j+1) 
          sublim(i,j) = a*sublim(i,j) + 0.5*(1-a)*sublim(i-1,j)
     &                                + 0.5*(1-a)*sublim(i,j+1)
        enddo
        enddo
      endif

      deallocate(Tin)
      return
      end

!============================================
      complex function K0(x)
! Modified Bessel function K0(xsqrt(i))
! Abramowitz and Stegun 9.9.2 (p379) and 9.11
!============================================
      implicit none 
      real x
      double precision berx,beix,kerx,keix,y,pi
      parameter (pi=3.14159265) 

      x = min(x, 8.) 
      y = x / 8.
      berx = 1d0 - 64d0*y**4 + 113.77777774D0*y**8 
     &       - 32.36345652D0*y**12 + 2.64191397D0*y**16
     &       - 0.08349609D0*y**20 + 0.00122552D0*y**24
     &       - 0.00000901D0*y**28
      beix = 16.D0*y**2 - 113.77777774D0*y**6
     &       + 72.81777742D0*y**10 - 10.56765779D0*y**14
     &       + 0.52185615D0*y**18 - 0.01103667D0*y**22
     &       + 0.00011346D0*y**26
      kerx = - log(x/2)*berx + (pi/4)*beix - 0.57721566
     &       - 59.05819744D0*y**4 + 171.36272133D0*y**8
     &       - 60.60977451D0*y**12 + 5.65539121D0*y**16
     &       - 0.19636347D0*y**20 + 0.00309699D0*y**24
     &       - 0.00002458D0*y**28
      keix = - log(x/2)*beix - (pi/4)*berx + 6.76454936D0*y**2
     &       - 142.91827687D0*y**6 + 124.23569650D0*y**10
     &       - 21.30060904D0*y**14 + 1.17509064D0*y**18
     &       - 0.02695875D0*y**22 + 0.00029532D0*y**26

      K0 = sngl(kerx) + (0.,1.)*sngl(keix)
   
      return
      end

!===================================
      subroutine FFT2(fr,fi,n,isign)
! 2-D Fast Fourier Transform
!===================================
      integer n,isign
      real fr(0:n-1,0:n-1),fi(0:n-1,0:n-1)
      integer i,j,k,nn(2)
      real norm !data(:),

      real, dimension(:), allocatable :: data
	  
      allocate(data(2*n*n))
	  
      nn(1) = n
      nn(2) = n
      k = 1
      do j=0,n-1
      do i=0,n-1
        data(k) = fr(i,j)
        data(k+1) = fi(i,j)
        k = k + 2
      enddo
      enddo
      call fourn(data,nn,2,isign)
      k = 1
      norm = 1.
      if (isign.eq.-1) norm = n*n
      do j=0,n-1
      do i=0,n-1
        fr(i,j) = data(k) / norm
        fi(i,j) = data(k+1) / norm
        k = k + 2
      enddo
      enddo

      deallocate(data)
      return
      end

!=========================================
      subroutine FOURN(data,nn,ndim,isign)
! ndim-dimensional FFT
! Numerical Recipes in Fortran p518
!=========================================
      INTEGER isign,ndim,nn(ndim)
      REAL data(*)
      INTEGER i1,i2,i2rev,i3,i3rev,ibit,idim,ifp1,ifp2,ip1,ip2,ip3,k1,
     *k2,n,nprev,nrem,ntot
      REAL tempi,tempr
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      ntot=1
      do 11 idim=1,ndim
        ntot=ntot*nn(idim)
11    continue
      nprev=1
      do 18 idim=1,ndim
        n=nn(idim)
        nrem=ntot/(n*nprev)
        ip1=2*nprev
        ip2=ip1*n
        ip3=ip2*nrem
        i2rev=1
        do 14 i2=1,ip2,ip1
          if(i2.lt.i2rev)then
            do 13 i1=i2,i2+ip1-2,2
              do 12 i3=i1,ip3,ip2
                i3rev=i2rev+i3-i2
                tempr=data(i3)
                tempi=data(i3+1)
                data(i3)=data(i3rev)
                data(i3+1)=data(i3rev+1)
                data(i3rev)=tempr
                data(i3rev+1)=tempi
12            continue
13          continue
          endif
          ibit=ip2/2
1         if ((ibit.ge.ip1).and.(i2rev.gt.ibit)) then
            i2rev=i2rev-ibit
            ibit=ibit/2
          goto 1
          endif
          i2rev=i2rev+ibit
14      continue
        ifp1=ip1
2       if(ifp1.lt.ip2)then
          ifp2=2*ifp1
          theta=isign*6.28318530717959d0/(ifp2/ip1)
          wpr=-2.d0*sin(0.5d0*theta)**2
          wpi=sin(theta)
          wr=1.d0
          wi=0.d0
          do 17 i3=1,ifp1,ip1
            do 16 i1=i3,i3+ip1-2,2
              do 15 i2=i1,ip3,ifp2
                k1=i2
                k2=k1+ifp1
                tempr=sngl(wr)*data(k2)-sngl(wi)*data(k2+1)
                tempi=sngl(wr)*data(k2+1)+sngl(wi)*data(k2)
                data(k2)=data(k1)-tempr
                data(k2+1)=data(k1+1)-tempi
                data(k1)=data(k1)+tempr
                data(k1+1)=data(k1+1)+tempi
15            continue
16          continue
            wtemp=wr
            wr=wr*wpr-wi*wpi+wr
            wi=wi*wpr+wtemp*wpi+wi
17        continue
          ifp1=ifp2
        goto 2
        endif
        nprev=n*nprev
18    continue

      return
      end
