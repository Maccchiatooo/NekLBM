c-----------------------------------------------------------------------
      subroutine nek_init(comm)
c
      include 'SIZE'
      include 'TOTAL'
      include 'DOMAIN'
c
      include 'OPCTR'
      include 'CTIMER'
      include 'LBMD3Q19'

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
      call initLBM_19
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

      call gengeom(igeom)  ! Generate geometry, after usrdat 

      if (ifmvbd) call setup_mesh_dssum ! Set mesh dssum (needs geom)

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
      include 'LBMD3Q19'
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



      do kstep=1,63660
            
         call nek_advance
         call userchk
         if(kstep.eq.38546)then
            call outpost(velxn,velyn,velzn,velzn,velzn,'un')
         endif
         if(mod(kstep,6366).eq.0)then
            call outpost(velxn,velyn,velzn,velzn,velzn,'   ')
          endif
        
      enddo

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
        include 'LBMD3Q19'

        real tmpvel2,tmpvel3,tmpw2,xa,ya,za
        real uc,Re,cs2,delta,kappa,machN,taui,R,rr
        integer         i,q,kx1,kx2,ky1,ky2,kz1,kz2
        integer ix,iy,iz,e,f,nto,ia
        real zlength,xr,yr,th,zo,rlam
        real amp_z,freq_z,freq_t,amp_tht,amp_clip
        real blt,phase_z,arg_tht,amp_sin
        integer npts
        real big,ieg,rand
        real kz,kx,xl1,zl1,alpha,beta,rand1,rand2
        real tmpvel1(lx1*ly1*lz1*lelv)
        real eps
        cs2=1.d0/3.d0
        npts=lx1*ly1*lz1*lelv
        machN=0.04d0
        uc=0.1!machN*sqrt(cs2)
        Re=1600.d0
        taui=0.5+uc/dt/cs2/Re/2/pi


        do i=1,npts
          xa=xm1(i,1,1,1)
          ya=ym1(i,1,1,1)   
          za=zm1(i,1,1,1)  
            
          velxn(i,1)= uc*sin(2.d0*pi*xa)*cos(2.d0*pi*ya)
     $              * cos(2.d0*pi*za)
          velyn(i,1)=-uc*cos(2.d0*pi*xa)*sin(2.d0*pi*ya)
     $               *cos(2.d0*pi*za)
          velzn(i,1)=0.d0
          denstN(i)=1.0*cs2+1.d0/16.d0*0.01*(cos(4.d0*pi*xa)
     $     +cos(4.d0*pi*ya))*(cos(4.d0*pi*za)+2.d0)
        enddo
        do i=1,npts
      tmpvel1(i)=0.5*(velxn(i,1)**2+velyn(i,1)**2+velzn(i,1)**2)/cs2
            enddo


        do q=1,q_vel
            do i=1,npts
            tmpvel2 = (e_vx(q)*velxn(i,1)
     $                  + e_vy(q)*velyn(i,1)
     $                  + e_vz(q)*velzn(i,1))/cs2

            tmpvel3 = (e_vx(q)*velxn(i,1)
     $                  + e_vy(q)*velyn(i,1)
     $              + e_vz(q)*velzn(i,1))**2/cs2/cs2/2.d0

            tmpw2   =  tmpvel1(i)-tmpvel2-tmpvel3

            f_1(i,1,q) = w_q(q)*denstN(i)/cs2-w_q(q)*tmpw2
        enddo
        enddo

        
        return
      endsubroutine

!-----------------------------------------
      subroutine lbm_collision_3
        implicit none
        include 'SIZE'
        include 'TOTAL'
        include 'LBMD3Q19'

                        
        real tmpvel2,tmpvel3,tmpw2,f_1_eq
        real uc,Re,cs2,delta,kappa,machN,taui
        integer i,q
        integer npts
        real tmpvel1(lx1*ly1*lz1*lelv)


        npts=lx1*ly1*lz1*lelv
        Re=1600.d0
        cs2=1.d0/3.d0
        machN=0.04d0
        uc=0.1!machN*sqrt(cs2)
        taui=0.5+uc/dt/cs2/Re/2/pi


        do i=1,npts
       tmpvel1(i)=0.5*(velxn(i,1)**2+velyn(i,1)**2+velzn(i,1)**2)/cs2
       enddo
        do q=1,q_vel
        do i=1,npts
            tmpvel2 = (e_vx(q)*velxn(i,1)
     $                  + e_vy(q)*velyn(i,1)
     $                  + e_vz(q)*velzn(i,1))/cs2

            tmpvel3 = (e_vx(q)*velxn(i,1)
     $                  + e_vy(q)*velyn(i,1)
     $              + e_vz(q)*velzn(i,1))**2/cs2/cs2/2.d0

            tmpw2   = tmpvel1(i)-tmpvel2-tmpvel3

            f_1_eq    = w_q(q)*denstN(i)/cs2-w_q(q)*tmpw2
            f_1(i,1,q)  = f_1(i,1,q)-(f_1(i,1,q)-f_1_eq)/taui
        enddo
        enddo

        return
      endsubroutine

!-----------------------------------------
      subroutine lbm_propogation
        implicit none
        include 'SIZE'
        include 'TOTAL'
        include 'LBMD3Q19'
        
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
          include 'LBMD3Q19'
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
          



! !-----------------------------------------
      subroutine lbm_update
        implicit none
        include 'SIZE'
        include 'TOTAL'
        include 'LBMD3Q19'
        call q_filter_lbm_3(param(103))   
        call lbm_density_3
        call lbm_velocity_3
      endsubroutine
! !-----------------------------------------

      subroutine lbm_density_3
        implicit none
        include 'SIZE'
        include 'TOTAL'
        include 'LBMD3Q19'

        integer    i,q
        integer npts

        npts=lx1*ly1*lz1*lelv

        do i =1,npts
            denstN(i)=0.d0
        enddo

      do q = 1,q_vel
        do i = 1,npts

            denstN(i)=denstN(i)+f_1(i,1,q)/3.d0
        enddo
        enddo



        return  
      endsubroutine
! !-----------------------------------------
      subroutine lbm_velocity_3
      !   implicit none
        include 'SIZE'
        include 'TOTAL'
        include 'CTIMER'
        include 'LBMD3Q19'

        integer i,q
        real ya,xa,za
        real ek(lx1*ly1*lz1*lelv)
        real dissip(lx1*ly1*lz1*lelv)
        real eksum,dissipsum
        real uc
        integer npts
        npts=lx1*ly1*lz1*lelv

        do i = 1,npts
            velxn(i,1)=0.d0
            velyn(i,1)=0.d0
            velzn(i,1)=0.d0
        enddo

        do q = 1,q_vel
        do i = 1,npts
           velxn(i,1)=velxn(i,1)+e_vx(q)*f_1(i,1,q)
           velyn(i,1)=velyn(i,1)+e_vy(q)*f_1(i,1,q)
           velzn(i,1)=velzn(i,1)+e_vz(q)*f_1(i,1,q)
        enddo
        enddo

        return
      endsubroutine

      subroutine initLBM_19
C--------------------------------------------------------------------
C
C     Initialize LBM
C
C--------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'
      include 'LBMD3Q19'

!      6  2  5 
!        \|/
!      3--+--1
!        /|\
!      7  4  8

c     Set default logicals
            integer::npts,NTOT1

            NTOT1=lx1*ly1*lz1*lelv
            npts=lx1*ly1*lz1*lelv*q_vel

            CALL RZERO(velxn,NTOT1)
            CALL RZERO(velyn,NTOT1)
            CALL RZERO(velzn,NTOT1)
            CALL RZERO(denstN,NTOT1)
            CALL RZERO(f_1,npts)

            e_vx(19) =  0.d0
            e_vy(19) =  0.d0
            e_vz(19) =  0.d0

            e_vx(1) =  1.d0
            e_vy(1) =  0.d0
            e_vz(1) =  0.d0

            e_vx(10) =  -1.d0
            e_vy(10) =  0.d0
            e_vz(10) =  0.d0            

            e_vx(2) =  0.d0
            e_vy(2) =  1.d0
            e_vz(2) =  0.d0

            e_vx(11) =  0.d0
            e_vy(11) =  -1.d0
            e_vz(11) =  0.d0

            e_vx(3) =  0.d0
            e_vy(3) =  0.d0
            e_vz(3) =  1.d0

            e_vx(12) =  0.d0
            e_vy(12) =  0.d0
            e_vz(12) =  -1.d0


            e_vx(4) =  -1.d0
            e_vy(4) =  -1.d0
            e_vz(4) =  0.d0

            e_vx(13) =  1.d0
            e_vy(13) =  1.d0
            e_vz(13) =  0.d0


            e_vx(5) =  1.d0
            e_vy(5) =  -1.d0
            e_vz(5) =  0.d0


            e_vx(14) =  -1.d0
            e_vy(14) =  1.d0
            e_vz(14) =  0.d0


            e_vx(6) =  -1.d0
            e_vy(6) =  0.d0
            e_vz(6) =  1.d0

            e_vx(15) =  1.d0
            e_vy(15) =  0.d0
            e_vz(15) =  -1.d0            

            e_vx(7) =  0.d0
            e_vy(7) =  -1.d0
            e_vz(7) =  1.d0

            e_vx(16) =  0.d0
            e_vy(16) =  1.d0
            e_vz(16) =  -1.d0

            e_vx(8) =  1.d0
            e_vy(8) =  0.d0
            e_vz(8) =  1.d0

            e_vx(17) =  -1.d0
            e_vy(17) =  0.d0
            e_vz(17) =  -1.d0

            e_vx(9) =  0.d0
            e_vy(9) =  1.d0
            e_vz(9) =  1.d0

            e_vx(18) =  0.d0
            e_vy(18) =  -1.d0
            e_vz(18) =  -1.d0

            bb(1) = 10
            bb(2) = 11
            bb(3) = 12
            bb(4) = 13
            bb(5) = 14
            bb(6) = 15
            bb(7) = 16
            bb(8) = 17
            bb(9) = 18
            bb(10) = 1
            bb(11) = 2
            bb(12) = 3
            bb(13) = 4
            bb(14) = 5
            bb(15) = 6
            bb(16) = 7
            bb(17) = 8
            bb(18) = 9
            bb(19) = 19



            w_q(19) = 1.0 / 3.0
            w_q(1) = 1.0 / 18.0
            w_q(10) = 1.0 / 18.0
            w_q(2) = 1.0 / 18.0
            w_q(11) = 1.0 / 18.0
            w_q(3) = 1.0 / 18.0
            w_q(12) = 1.0 / 18.0
            w_q(4) = 1.0 / 36.0
            w_q(13) = 1.0 / 36.0
            w_q(5) = 1.0 / 36.0
            w_q(14) = 1.0 / 36.0
            w_q(6) = 1.0 / 36.0
            w_q(15) =1.0 / 36.0
            w_q(7) = 1.0 / 36.0
            w_q(16) = 1.0 / 36.0
            w_q(8) = 1.0 / 36.0
            w_q(17) =1.0 / 36.0
            w_q(9) =1.0 / 36.0
            w_q(18) =1.0 / 36.0
            


      RETURN
      END

      subroutine lbm_rk3_q_bb(nq)
        implicit none
        include 'SIZE'
        include 'TOTAL'
        include 'LBMD3Q19'

        integer q,i,npts,nq
        real rk_a,rk_b1,rk_b2,rk_b3,rk_b4
        real ex,ey,ez,ex_b,ey_b,ez_b
        real uc,machN,Re,cs2,dp

        
        npts=lx1*ly1*lz1*lelv

        ex=e_vx(nq)
        ey=e_vy(nq)    
        ez=e_vz(nq)    

        ex_b=e_vx(bb(nq))
        ey_b=e_vy(bb(nq))    
        ez_b=e_vz(bb(nq))    
        call gradm1(dfx,dfy,dfz,f_q)
        call gradm1(dfx_b,dfy_b,dfz_b,f_bq)

      do i=1,npts
       conv_b(i)=(-ex*dfx(i)-ey*dfy(i)-ez*dfz(i))
        conv_1b(i)=(-ex_b*dfx_b(i)-ey_b*dfy_b(i)-ez_b*dfz_b(i))
        enddo

          call col2 (conv_b,bm1,npts) 
          call col2 (conv_1b,bm1,npts) 
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
        enddo
        call col2 (conv,bm1,npts)
        call col2 (conv_2,bm1,npts)
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
      enddo
      call col2 (conv,bm1,npts) 
      call col2 (conv_2,bm1,npts) 
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
      include 'LBMD3Q19'

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
