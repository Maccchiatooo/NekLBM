c-----------------------------------------------------------------------
      subroutine nek_init(comm)
c
      include 'SIZE'
      include 'TOTAL'
      include 'DOMAIN'
c
      include 'OPCTR'
      include 'CTIMER'
      include 'LBMD3Q13'

      integer comm
      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal
  
      common /rdump/ ntdump

      real kwave2
      logical ifemati

      real rtest
      integer itest
      integer*8 itest8
      character ctest
      logical ltest 

      common /c_is1/ glo_num(lx1 * ly1 * lz1, lelt)
      common /ivrtx/ vertex((2 ** ldim) * lelt)
      integer*8 glo_num, ngv
      integer*8 vertex

      ! set word size for REAL
      wdsize = sizeof(rtest)
      ! set word size for INTEGER
      isize = sizeof(itest)
      ! set word size for INTEGER*8
      isize8 = sizeof(itest8) 
      ! set word size for LOGICAL
      lsize = sizeof(ltest) 
      ! set word size for CHARACTER
      csize = sizeof(ctest)

      call setupcomm(comm,newcomm,newcommg,'','')
      intracomm   = newcomm   ! within a session
      nekcomm     = newcomm
      iglobalcomm = newcommg  ! across all sessions
      call iniproc()

      if (nid.eq.nio) call printHeader

      etimes = dnekclock()
      istep  = 0

      call opcount(1)

      call initdim         ! Initialize / set default values.
      call initdat
      call initLBM_13
      call files

      call readat          ! Read .rea +map file

      if (nio.eq.0) then
         write(6,12) 'nelgt/nelgv/lelt:',nelgt,nelgv,lelt
         write(6,12) 'lx1/lx2/lx3/lxd: ',lx1,lx2,lx3,lxd
 12      format(1X,A,4I12)
         write(6,*)
      endif

      call setvar          ! Initialize most variables

      instep=1             ! Check for zero steps
      if (nsteps.eq.0 .and. fintim.eq.0.) instep=0

      igeom = 2
      call setup_topo      ! Setup domain topology  

      call genwz           ! Compute GLL points, weights, etc.
      if(nio.eq.0) write(6,*) 'call usrdat'
      call usrdat
      if(nio.eq.0) write(6,'(A,/)') ' done :: usrdat' 
      call gengeom(igeom)  ! Generate geometry, after usrdat 

      if (ifmvbd) call setup_mesh_dssum ! Set mesh dssum (needs geom)
      if(nio.eq.0) write(6,*) 'call usrdat2'
      call usrdat2
      if(nio.eq.0) write(6,'(A,/)') ' done :: usrdat2' 
      call count_bdry   ! count the number of faces with assigned BCs
      call fix_geom

      call vrdsmsh          ! verify mesh topology
      call mesh_metrics     ! print some metrics

      call setlog(.true.)   ! Initalize logical flags

      if (ifneknekc) call neknek_setup

      call bcmask  ! Set BC masks for Dirichlet boundaries.

      if (fintim.ne.0.0 .or. nsteps.ne.0) 
     $   call geneig(igeom) ! eigvals for tolerances


      if(ifcvode) call cv_setsize
      if(nio.eq.0) write(6,*) 'call usrdat3'
      call usrdat3
      if(nio.eq.0) write(6,'(A,/)') ' done :: usrdat3'
      call setics
      call setprop

      if (ifcvode .and. nsteps.gt.0) call cv_init

      call comment
      call sstest (isss) 

      call dofcnt

      jp = 0  ! Set perturbation field count to 0 for baseline flow
      p0thn = p0th

      call in_situ_init()

      call time00       !     Initalize timers to ZERO
      call opcount(2)

      ntdump=0
      if (timeio.ne.0.0) ntdump = int( time/timeio )

      tinit = dnekclock_sync() - etimes
      if (nio.eq.0) then
        write (6,*) ' '
        if (time.ne.0.0) write (6,'(a,e14.7)') ' Initial time:',time
        write (6,'(a,g13.5,a)') 
     &     ' Initialization successfully completed ', tinit, ' sec'
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine nek_solve

      include 'SIZE'
      include 'TSTEP'
      include 'INPUT'
      include 'CTIMER'
      include 'LBMD3Q13'

      integer :: unit_number, iostat
      character(len=100) :: filename
      call nekgsync()

      if (instep.eq.0) then
        if(nid.eq.0) write(6,'(/,A,/,A,/)') 
     &     ' nsteps=0 -> skip time loop',
     &     ' running solver in post processing mode'
      else
        if(nio.eq.0) write(6,'(/,A,/)') 'Starting time loop ...'
      endif

      isyc  = 0
      if(ifsync) isyc=1
      itime = 0
#ifdef TIMER
      itime = 1
#endif

      ! start measurements
      dtmp = dnekgflops()

      istep  = 0
      msteps = 1

      irstat = int(param(120))
      call lbm_initialization_3
      call outpost(velxn,velyn,velzn,velzn,velzn,'   ')

      unit_number=123
      filename = "energy.txt"
      if(nid.eq.0)then
      open(unit=unit_number, file=filename, status="unknown", 
     $   action="write", iostat=iostat)
      endif


      do kstep=1,nsteps
            
         call nek_advance
         call userchk

         if(mod(kstep,100000).eq.0)then
            call outpost(velxn,velyn,velzn,velzn,velzn,'   ')
            
         endif
        call wall_avg
        call wall_shear
      enddo
      


      if(nio.eq.0) close(unit=123)
 1001 lastep=1

      call comment

c     check for post-processing mode
      if (instep.eq.0) then
         nsteps=0
         istep=0
         if(nio.eq.0) write(6,*) 'call userchk'
         call userchk
         if(nio.eq.0) write(6,*) 'done :: userchk'
         call prepost (.true.,'his')
      else
         if (nio.eq.0) write(6,'(/,A,/)') 
     $      'end of time-step loop' 
      endif


      RETURN
      END

c-----------------------------------------------------------------------
      subroutine nek_advance

      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'

      if (iftran) call settime
      call comment
      call lbm_solver

      return      
      end

c-----------------------------------------------------------------------
      subroutine nek_end

      include 'SIZE'
      include 'TOTAL'
      include 'DPROCMAP'

      if(instep.ne.0) call runstat


#ifdef DPROCMAP
#ifdef MPI
      call MPI_Win_free(dProcmapH, ierr)
#endif
#endif 
      call in_situ_end()
      call exitt0()
      return
      end

      subroutine lbm_solver
  
        call lbm_collision_3
        call lbm_propogation
        call lbm_update
      
      endsubroutine

        subroutine lbm_initialization_3
        implicit none
        include 'SIZE'
        include 'TOTAL'
        include 'LBMD3Q13'

        real tmpvel1,tmpvel2,tmpvel3,tmpw2,xa,ya,za
        real uc,Re,cs2,delta,kappa,machN,taui,R,rr
        integer         i,q,kx1,kx2,ky1,ky2,kz1,kz2
        integer ix,iy,iz,ie,iface,nto,ia
        real zlength,xr,yr,th,zo,rlam
        real amp_z,freq_z,freq_t,amp_tht,amp_clip
        real blt,phase_z,arg_tht,amp_sin
        integer npts
        real big,ieg,rand
        real kz,kx,xl1,zl1,alpha,beta,rand1,rand2
        real rand3,eps
        real yp,up,lp,u_tau

        tim_t=0
        cs2=1.d0
        npts=lx1*ly1*lz1*lelv
        machN=0.1d0/sqrt(cs2)
        uc=machN*sqrt(cs2)
        R=0.5d0
        Re=5300
        zlength=4.d0
        dp=0.00018454966
        taui=uc/Re
        u_tau=sqrt(dp*R/2.d0  )


        do ix=1,lx1
        do iy=1,ly1
        do iz=1,lz1
        do ie=1,nelv
        i=ix+(iy-1)*lx1+(iz-1)*lx1**2+(ie-1)*lx1**3

          xa=xm1(i,1,1,1)
          ya=ym1(i,1,1,1)   
          za=zm1(i,1,1,1)  
      
        xr = xa/R ! convert -R<x<R to -1<xr<1
        yr = ya/R ! convert -R<y<R to -1<yr<1
        rr = xr*xr + yr*yr
        if (rr.gt.0) rr=sqrt(rr)
        
        th = atan2(ya,xa) 
        zo = 2*pi*za/zlength ! convert 0<z<zlength to 0<zo<2*pi


          rlam = sqrt(xa**2 + ya**2)
          lp=(R-rlam)*u_tau/taui
          if(lp.ge.10)then
          velzn(i,1) =(1.d0/0.4*log(lp)+5.d0)*u_tau!uc*2*(1-rlam**2/R**2)!
          else
          velzn(i,1) =lp*u_tau
          endif

          velxn(i,1)=0.d0
          velyn(i,1)=0.d0
          denstN(i)=1.0

        enddo
        enddo
        enddo
        enddo


        do i=1,npts

       tmpvel1=0.5*(velxn(i,1)**2+velyn(i,1)**2+velzn(i,1)**2)/cs2
        do q=1,q_vel
            tmpvel2 = (e_vx(q)*velxn(i,1)
     $                  + e_vy(q)*velyn(i,1)
     $                  + e_vz(q)*velzn(i,1))/cs2

            tmpvel3 = (e_vx(q)*velxn(i,1)
     $                  + e_vy(q)*velyn(i,1)
     $              + e_vz(q)*velzn(i,1))**2/cs2/cs2/2.d0

            tmpw2   = 1.0 - tmpvel1+tmpvel2+tmpvel3

            f_1(i,1,q) = w_q(q)*denstN(i)*tmpw2
        enddo
        enddo

        
        return
      endsubroutine

!-----------------------------------------
      subroutine lbm_collision_3
        implicit none
        include 'SIZE'
        include 'TOTAL'
        include 'LBMD3Q13'
                        
        real tmpvel1,tmpvel2,tmpvel3,tmpw2,f_1_eq
        real uc,Re,cs2,delta,kappa,machN,taui
        integer i,q
        integer npts
        real R,kz,k_th,l0,l1,b0,zo,xa,ya,za
        real grav,zlength,tauii,u_tau,rlam
        real lp,xr,yr,rr,th

        npts=lx1*ly1*lz1*lelv
        Re=5300.d0
        cs2=1.d0
        machN=0.1d0/sqrt(cs2)
        uc=machN*sqrt(cs2)
        taui=0.5+uc/dt/cs2/Re


        dp=0.00018454966!32.d0*uc**2/Re

        R=0.5d0
        kz=3.d0
        k_th=2.d0
        l0=0.2*R
        l1=0.4*R
        b0=50.d0
        kappa=0.5d0
        grav=dp
        zlength=4.0
        tauii=uc/Re
        u_tau=sqrt(grav*R/2.d0)
        
      do i=1,npts
        dpx(i)=0
        dpy(i)=0
        dpz(i)=dp
        dpr(i)=0
        dpth(i)=0
      enddo

      if(tim_t.le.220/dt)then
      do i=1,npts

      xa=xm1(i,1,1,1)
      ya=ym1(i,1,1,1)   
      za=zm1(i,1,1,1)  
      
      xr = xa/R ! convert -R<x<R to -1<xr<1
      yr = ya/R ! convert -R<y<R to -1<yr<1
      rr = xr*xr + yr*yr
      if (rr.gt.0) rr=sqrt(rr)
        
      th = atan2(ya,xa) 
      zo = 2*pi*za/zlength ! convert 0<z<zlength to 0<zo<2*pi

      rlam = sqrt(xa**2 + ya**2)
      lp=(R-rlam)*u_tau/taui

      if(R-rlam.lt.l0+l1.and.R-rlam.gt.l0)then
      dpr(i)=-grav*kappa*b0*R/rlam*kz*l1/4.0*
     $        sin(2.d0*pi*tim_t/(220/dt*8))*
     $        (1.d0-cos(2.d0*pi*(R-rlam-l0)/l1))*
     $        cos(kz*2.d0*pi*za/4.0)*cos(k_th*th)

      dpth(i)=grav*(1.d0-kappa)*b0*kz/k_th*2*pi*R/4.0*
     $        sin(2.d0*pi*tim_t/(220/dt*8))*
     $        sin(2.d0*pi*(R-rlam-l0)/l1)*
     $        cos(kz*2.d0*pi*za/4.0) *sin(k_th*th)

      dpz(i)=dp-grav*b0*R/rlam*sin(2.d0*pi*tim_t/(220/dt*8))*
     $        sin(2.d0*pi*(R-rlam-l0)/l1)*
     $        sin(kz*2.d0*pi*za/4.0)*cos(k_th*th)

      dpx(i)=dpr(i)*cos(th)-dpth(i)*sin(th)
      dpy(i)=dpr(i)*sin(th)+dpth(i)*cos(th)

      endif
      enddo
      endif

        do i=1,npts
       tmpvel1=0.5*(velxn(i,1)**2+velyn(i,1)**2+velzn(i,1)**2)/cs2
        do q=1,q_vel
            tmpvel2 = (e_vx(q)*velxn(i,1)
     $                  + e_vy(q)*velyn(i,1)
     $                  + e_vz(q)*velzn(i,1))/cs2

            tmpvel3 = (e_vx(q)*velxn(i,1)
     $                  + e_vy(q)*velyn(i,1)
     $              + e_vz(q)*velzn(i,1))**2/cs2/cs2/2.d0

            tmpw2   = 1.0 - tmpvel1+tmpvel2+tmpvel3

            force_p=w_q(q)*(((e_vx(q)*velxn(i,1)+e_vy(q)*velyn(i,1)+
     $       +e_vz(q)*velzn(i,1))*e_vx(q)
     $       -cs2*velxn(i,1))/cs2/cs2*dpx(i)+
     $      ((e_vx(q)*velxn(i,1)+e_vy(q)*velyn(i,1)+
     $       +e_vz(q)*velzn(i,1))*e_vy(q)
     $       -cs2*velyn(i,1))/cs2/cs2*dpy(i)+
     $      ((e_vx(q)*velxn(i,1)+e_vy(q)*velyn(i,1)+
     $       +e_vz(q)*velzn(i,1))*e_vz(q)
     $       -cs2*velzn(i,1))/cs2/cs2*dpz(i))



            f_1_eq    = w_q(q)*denstN(i)*tmpw2-0.5*dt*force_p

            f_1(i,1,q)  = f_1(i,1,q)-(f_1(i,1,q)-f_1_eq)/taui+
     $       dt*force_p

        enddo
        enddo

        return
      endsubroutine

!-----------------------------------------
      subroutine lbm_propogation
        implicit none
        include 'SIZE'
        include 'TOTAL'
        include 'LBMD3Q13'
        
        integer q

        do q=1,(q_vel-1)/2
          call convect_q_bb(q)
        enddo

        return
      endsubroutine
! !-----------------------------------------

! ! !-----------------------------------------
          subroutine convect_q_bb(nq)
          implicit none
          include 'SIZE'
          include 'TOTAL'
          include 'LBMD3Q13'
          integer nq,i,npts
          
          npts=lx1*ly1*lz1*lelv
  

          do i=1,npts
            f_q(i)=f_1(i,1,nq)
            f_bq(i)=f_1(i,1,bb(nq))
          enddo
  
          call lbm_rk3_q_bb(nq)
  
          do i=1,npts
            f_1(i,1,nq)=f_q(i)
            f_1(i,1,bb(nq))=f_bq(i)
          enddo
          return
          endsubroutine
          

      subroutine lbm_boundary(cv,cvb,f_tem,f_tem_b,nq)
            implicit none
            include 'SIZE'
            include 'TOTAL'
            include 'LBMD3Q13'

            integer :: i,npts,nq
            real :: xa,ya,f_tem(lx1*ly1*lz1*lelv)
            real :: f_tem_b(lx1*ly1*lz1*lelv)
            real :: flx_b(lx1*ly1*lz1*lelv)
            real :: cv(lx1*ly1*lz1*lelv)
            real :: cvb(lx1*ly1*lz1*lelv)
            integer :: bbq 
            real :: flx_bb(lx1*ly1*lz1*lelv)


            npts=lx1*ly1*lz1*lelv
            bbq = bb(nq)     

            call compute_flx(f_tem,f_tem_b,flx_b,nq)
        call compute_flx(f_tem_b,f_tem,flx_bb,bbq)

        
            do i=1,npts
              cv(i)=cv(i)+flx_b(i)
         cvb(i)=cvb(i)+flx_bb(i)
            enddo
        return
      endsubroutine

                             
      subroutine compute_flx(f_tem,f_tem_b,flx_b,nq)
        ! implicit none
        include 'SIZE'
        include 'TOTAL'
        include 'LBMD3Q13'

        integer :: i,npts,nq,e,nto,ix,iy,iz,f,ia
        real :: xa,ya,za,f_tem(1),f_tem_b(1)
        real :: flx_b(1)
        real :: flx_t(1)
        real :: ex,ey,ez
        real :: norm_x,norm_y,norm_z,a_w,nme
        real :: diff_f

        common /finewts/ wghtc(lx1*lz1)



        npts=lx1*ly1*lz1*lelv

        ex=e_vx(nq)
        ey=e_vy(nq)
        ez=e_vz(nq)

        do i=1,npts
            flx_b(i)=0.d0
        enddo

        do e=1,nelv
        do f=1,2*ndim 

        if (cbc(f,e,1).eq.'W  ') then
        call facind (kx1,kx2,ky1,ky2,kz1,kz2,lx1,ly1,lz1,f)
        ia=0
              do iz=kz1,kz2
              do iy=ky1,ky2
              do ix=kx1,kx2
              ia=ia+1

              nto=ix+(iy-1)*lx1+(iz-1)*lx1**2
     $            +(e-1)*lx1**3
              
      diff_f=(f_tem(nto)-f_tem_b(nto))
              norm_x=unx(ia,1,f,e)
              norm_y=uny(ia,1,f,e)
              norm_z=unz(ia,1,f,e)
              nme=norm_x*ex+norm_y*ey+norm_z*ez
              a_w=area(ia,1,f,e)

              if(nme<0)then
              flx_b(nto)=a_w*nme*diff_f
              endif

              enddo
              enddo
              enddo
          endif

        enddo
        enddo

        return
        endsubroutine


! !-----------------------------------------
      subroutine lbm_update
        implicit none
        include 'SIZE'
        include 'TOTAL'
        include 'LBMD3Q13'
        call q_filter_lbm_3(param(103))   
        call lbm_density_3
        call lbm_velocity_3
        tim_t=tim_t+1

      endsubroutine
! !-----------------------------------------

      subroutine lbm_density_3
        implicit none
        include 'SIZE'
        include 'TOTAL'
        include 'LBMD3Q13'

        integer    i,q
        integer npts
        npts=lx1*ly1*lz1*lelv

        do i =1,npts
            denstN(i)=0.d0
        enddo

        do q = 1,q_vel
        do i = 1,npts
            denstN(i)=denstN(i)+f_1(i,1,q)
        enddo
        enddo

        return  
      endsubroutine
! !-----------------------------------------
      subroutine lbm_velocity_3
        implicit none
        include 'SIZE'
        include 'TOTAL'
        include 'LBMD3Q13'

        integer i,q
        integer npts
        integer f,e,ix,iy,iz
        integer kx1,kx2,ky1,ky2,kz1,kz2,nto

        npts=lx1*ly1*lz1*lelv
        
        do i = 1,npts
            velxn(i,1)=0.d0
            velyn(i,1)=0.d0
            velzn(i,1)=0.d0
        enddo

          do q = 1,q_vel
          do i = 1,npts
           velxn(i,1)=velxn(i,1)+e_vx(q)*f_1(i,1,q)/denstN(i)
           velyn(i,1)=velyn(i,1)+e_vy(q)*f_1(i,1,q)/denstN(i)
           velzn(i,1)=velzn(i,1)+e_vz(q)*f_1(i,1,q)/denstN(i)
        enddo
        enddo

        return
      endsubroutine


      subroutine initLBM_13
C--------------------------------------------------------------------
C
C     Initialize LBM
C
C--------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'
      include 'LBMD3Q13'

!      6  2  5 
!        \|/
!      3--+--1
!        /|\
!      7  4  8

c     Set default logicals
            integer::npts,NTOT1
            real abscir,abscis
            abscir=sqrt((5.d0+sqrt(5.d0))/2.d0)
            abscis=sqrt((5.d0-sqrt(5.d0))/2.d0)

            NTOT1=lx1*ly1*lz1*lelv
            npts=lx1*ly1*lz1*lelv*q_vel

            CALL RZERO(velxn,NTOT1)
            CALL RZERO(velyn,NTOT1)
            CALL RZERO(velzn,NTOT1)
            CALL RZERO(denstN,NTOT1)
            CALL RZERO(f_1,npts)

            e_vx(13) =  0.d0
            e_vy(13) =  0.d0
            e_vz(13) =  0.d0

            e_vx(1) =  abscir
            e_vy(1) =  abscis
            e_vz(1) =  0.d0

            e_vx(7) =  -abscir
            e_vy(7) = -abscis
            e_vz(7) =  0.d0            

            e_vx(2) =  abscir
            e_vy(2) =  -abscis
            e_vz(2) =  0.d0

            e_vx(8) =  -abscir
            e_vy(8) =  abscis
            e_vz(8) =  0.d0

            e_vx(3) =  abscis
            e_vy(3) =  0.d0
            e_vz(3) =  abscir

            e_vx(9) =  -abscis
            e_vy(9) =  0.d0
            e_vz(9) =  -abscir

            e_vx(4) =  -abscis
            e_vy(4) =  0.d0
            e_vz(4) =  abscir

            e_vx(10) =  abscis
            e_vy(10) =  0.d0
            e_vz(10) =  -abscir


            e_vx(5) =  0.d0
            e_vy(5) =  abscir
            e_vz(5) =  abscis


            e_vx(11) = 0.d0
            e_vy(11) =  -abscir
            e_vz(11) =  -abscis


            e_vx(6) =  0.d0
            e_vy(6) =  -abscir
            e_vz(6) =  abscis

            e_vx(12) =  0.d0
            e_vy(12) =  abscir
            e_vz(12) =  -abscis         


            bb(1) = 7
            bb(2) = 8
            bb(3) = 9
            bb(4) = 10
            bb(5) = 11
            bb(6) = 12
            bb(7) = 1
            bb(8) = 2
            bb(9) = 3
            bb(10) = 4
            bb(11) = 5
            bb(12) = 6
            bb(13) =13
           


            w_q(13) = 2.0 / 5.0
            w_q(1) = 1.0 / 20.0
            w_q(2) = 1.0 / 20.0
            w_q(3) = 1.0 / 20.0
            w_q(4) = 1.0 / 20.0
            w_q(5) = 1.0 / 20.0
            w_q(6) = 1.0 / 20.0
            w_q(7) = 1.0 / 20.0
            w_q(8) = 1.0 / 20.0
            w_q(9) = 1.0 / 20.0
            w_q(10) = 1.0 / 20.0
            w_q(11) = 1.0 / 20.0
            w_q(12) = 1.0 / 20.0
            


      RETURN
      END


      subroutine lbm_rk3_q_bb(nq)
        implicit none
        include 'SIZE'
        include 'TOTAL'
        include 'LBMD3Q13'

        integer q,i,npts,nq
        real rk_a,rk_b1,rk_b2,rk_b3,rk_b4
        real ex,ey,ez,ex_b,ey_b,ez_b
        real uc,machN,Re,cs2
        integer f,e,ix,iy,iz,ie
        real R,xa,ya,za,tau_wsum,xr,yr
        real u_avesum,a_w
        integer kx1,kx2,ky1,ky2,kz1,kz2
        integer ia,nto,iia
        real rr,th,zo,lp,rlam
        real kz,k_th,l0,l1,b0,kappa
        real grav


        npts=lx1*ly1*lz1*lelv
        Re=5300.d0
        cs2=1.d0
        machN=0.1d0/sqrt(cs2)
        uc=machN*sqrt(cs2)
        dp=0.00018454966!32.d0*uc**2/Re
        R=0.5d0
        kz=3.d0
        k_th=2.d0
        l0=0.2*R
        l1=0.4*R
        b0=50.d0
        kappa=0.5d0
        grav=dp
!        zlength=4.0
!        taui=uc/Re
!        u_tau=sqrt(grav*R/2.d0)

        ex=e_vx(nq)
        ey=e_vy(nq)    
        ez=e_vz(nq)    

        ex_b=e_vx(bb(nq))
        ey_b=e_vy(bb(nq))    
        ez_b=e_vz(bb(nq))    
        
        call gradm1(dfx,dfy,dfz,f_q)
        call gradm1(dfx_b,dfy_b,dfz_b,f_bq)

      do i=1,npts
      xa=xm1(i,1,1,1)
      ya=ym1(i,1,1,1) 

      th = atan2(ya,xa)   
        dp1(i)=w_q(nq)*(ez*(dpz(i))/cs2+
     $   ey/cs2*dpy(i)+
     $   ex/cs2*dpx(i))
        dp2(i)=w_q(bb(nq))*(ez_b*(dpz(i))/cs2+
     $   ey_b/cs2*dpy(i)+
     $   ex_b/cs2*dpx(i))
      enddo
      do i=1,npts
       conv_b(i)=(-ex*dfx(i)-ey*dfy(i)-ez*dfz(i))
        conv_1b(i)=(-ex_b*dfx_b(i)-ey_b*dfy_b(i)-ez_b*dfz_b(i))
        conv_b(i)=conv_b(i)+dp1(i)
        conv_1b(i)=conv_1b(i)+dp2(i)
        enddo

          call col2 (conv_b,bm1,npts) 
          call col2 (conv_1b,bm1,npts) 
          call lbm_boundary(conv_b,conv_1b,f_q,f_bq,nq)
          call fgslib_gs_op(gsh_fld(1),conv_b,1,1,0)
          call fgslib_gs_op(gsh_fld(1),conv_1b,1,1,0)
          call col2(conv_b,binvm1, npts)
          call col2(conv_1b,binvm1, npts)



        do i =1,npts
          f_temp1(i)=f_q(i)+2.d0*dt*conv_b(i)/3.d0
          f_temp1_b(i)=f_bq(i)+2.d0*dt*conv_1b(i)/3.d0
        enddo  
        call gradm1(dfx,dfy,dfz,f_temp1)
        call gradm1(dfx_b,dfy_b,dfz_b,f_temp1_b)


      do i=1,npts
      conv(i)=(-ex*dfx(i)-ey*dfy(i)-ez*dfz(i))
      conv_2(i)=(-ex_b*dfx_b(i)-ey_b*dfy_b(i)-ez_b*dfz_b(i))

        conv(i)=conv(i)+dp1(i)
        conv_2(i)=conv_2(i)+dp2(i)
        enddo
        call col2 (conv,bm1,npts)
        call col2 (conv_2,bm1,npts)
        call lbm_boundary(conv,conv_2,f_temp1,f_temp1_b,nq)
        call fgslib_gs_op(gsh_fld(1),conv,1,1,0)
        call fgslib_gs_op(gsh_fld(1),conv_2,1,1,0)
        call col2(conv,binvm1, npts)
        call col2(conv_2,binvm1, npts)

      do i =1,npts
        f_temp2(i)=2.d0/3.d0*f_q(i)+f_temp1(i)/3.d0
        f_temp2(i)=f_temp2(i)+4.d0*conv(i)*dt/9.d0

        f_temp2_b(i)=2.d0/3.d0*f_bq(i)+f_temp1_b(i)/3.d0
        f_temp2_b(i)=f_temp2_b(i)+4.d0*conv_2(i)*dt/9.d0
      enddo  

      call gradm1(dfx,dfy,dfz,f_temp2)
      call gradm1(dfx_b,dfy_b,dfz_b,f_temp2_b)

      do i=1,npts
        conv(i)=(-ex*dfx(i)-ey*dfy(i)-ez*dfz(i))
        conv_2(i)=(-ex_b*dfx_b(i)-ey_b*dfy_b(i)-ez_b*dfz_b(i))
        conv(i)=conv(i)+dp1(i)
        conv_2(i)=conv_2(i)+dp2(i)

      enddo
      call col2 (conv,bm1,npts) 
      call col2 (conv_2,bm1,npts) 
      call lbm_boundary(conv,conv_2,f_temp2,f_temp2_b,nq)
      call fgslib_gs_op(gsh_fld(1),conv,1,1,0)
      call fgslib_gs_op(gsh_fld(1),conv_2,1,1,0)
      call col2(conv,binvm1, npts)
      call col2(conv_2,binvm1, npts)

      do i =1,npts
        f_q(i)=37.d0*f_q(i)/64.d0+5.d0/32.d0*dt*conv_b(i)
        f_q(i)=f_q(i)+27.d0*f_temp2(i)/64.d0+9.d0/16.d0*dt*conv(i)

        f_bq(i)=37.d0*f_bq(i)/64.d0+5.d0/32.d0*dt*conv_1b(i)
      f_bq(i)=f_bq(i)+27.d0*f_temp2_b(i)/64.d0+9.d0/16.d0*dt*conv_2(i)
      enddo 

        return
        endsubroutine


      subroutine q_filter_lbm_3(wght)
c
c     filter vx,vy,vz, and p by simple interpolation
c
      include 'SIZE'
      include 'TOTAL'
c
c
c     These are the dimensions that we interpolate onto for v and p:
      parameter(lxv=lx1-1)
      parameter(lxp=lx2-1)
c
      real intdv(lx1,lx1)
      real intuv(lx1,lx1)
      real intdp(lx1,lx1)
      real intup(lx1,lx1)
      real intv(lx1,lx1)
      real intp(lx1,lx1)
c
      save intdv
      save intuv
      save intdp
      save intup
      save intv
      save intp

      common /ctmp0/ intw,intt
     $             , wk1,wk2
     $             , zgmv,wgtv,zgmp,wgtp,tmax(100),omax(103)

      real intw(lx1,lx1)
      real intt(lx1,lx1)
      real wk1  (lx1,lx1,lx1,lelt)
      real wk2  (lx1,lx1,lx1)
      real zgmv(lx1),wgtv(lx1),zgmp(lx1),wgtp(lx1)
c
c     outpost arrays
      parameter (lt=lx1*ly1*lz1*lelv)
      common /scruz/ w1(lt),w2(lt),w3(lt),wt(lt)

      character*18 sfmt

      integer icalld
      save    icalld
      data    icalld /0/

      logical if_fltv

      ncut = param(101)+1

      if(wght.le.0) return
      if(ifaxis) call exitti('Filtering not supported w/ IFAXIS!$',1)
      if(nid.eq.0 .and. loglevel.gt.2) write(6,*) 'apply q_filter ',
     $                                            ncut, wght

      imax = nid
      imax = iglmax(imax,1)
      jmax = iglmax(imax,1)

      if (icalld.eq.0) then
         icalld = 1
         call build_new_filter(intv,zgm1,lx1,ncut,wght,nio)
      elseif (icalld.lt.0) then   ! old (std.) filter
         icalld = 1
         call zwgll(zgmv,wgtv,lxv)
         call igllm(intuv,intw,zgmv,zgm1,lxv,lx1,lxv,lx1)
         call igllm(intdv,intw,zgm1,zgmv,lx1,lxv,lx1,lxv)
c
         call zwgl (zgmp,wgtp,lxp)
         call iglm (intup,intw,zgmp,zgm2,lxp,lx2,lxp,lx2)
         call iglm (intdp,intw,zgm2,zgmp,lx2,lxp,lx2,lxp)
c
c        Multiply up and down interpolation into single operator
c
         call mxm(intup,lx2,intdp,lxp,intp,lx2)
         call mxm(intuv,lx1,intdv,lxv,intv,lx1)
c
c        Weight the filter to make it a smooth (as opposed to truncated)
c        decay in wave space
         w0 = 1.-wght
         call ident(intup,lx2)
         call add2sxy(intp,wght,intup,w0,lx2*lx2)

         call ident   (intuv,lx1)
         call add2sxy (intv ,wght,intuv,w0,lx1*lx1)
      endif
      call filterq_lbm_3(intv,lx1,lz1,wk1,wk2,intt)
      return
      end
! c-----------------------------------------------------------------------
      subroutine filterq_lbm_3(f,nx,nz,w1,w2,ft)
c filterq

      include 'SIZE'
      include 'TSTEP'
      include 'LBMD3Q13'

      real w1(1),w2(1)
c
      real f(nx,nx),ft(nx,nx)

      integer iii,q
c
      integer e
c
      call transpose(ft,nx,f,nx)
c
      nxyz=nx*nx*nz
      dmax = 0.
      emax = 0.

      nel = nelfld(1)

         do q=1,q_vel
         do e=1,nel
          call copy(w2,f_1(1,e,q),nxyz)
          call mxm(f ,nx,w2,nx,w1,nx*nx)
          i=1
          j=1
          do k=1,nx
             call mxm(w1(i),nx,ft,nx,w2(j),nx)
             i = i+nx*nx
             j = j+nx*nx
          enddo
          call mxm (w2,nx*nx,ft,nx,w1,nx)
          call sub3(w2,f_1(1,e,q),w1,nxyz)
          call copy(f_1(1,e,q),w1,nxyz)
        enddo
        enddo

      return
      end   
