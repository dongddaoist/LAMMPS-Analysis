!
!    ***  H-bond analysis by D.D.
!
      PROGRAM hbond_all
      implicit none
!
!
      integer nptypes,maxshots,MaxInpFiles,nread,nskip,nion,nss
      integer maxat,maxnch
      parameter (nion=125)
      parameter (nss=12000)
      parameter (maxshots=100000) !maximum number of snaps, see dim of xl(3,77,*)
      parameter (MaxInpFiles=40)       ! # of fort.77 files
      parameter (nread=1)        ! how often to read snapshots from fort.77
      parameter (nskip=0)
      parameter (maxat=74375)
      parameter (maxnch=18000)
!
      real*4 xx(3,maxat),rboxgen(3),el2gen(3)
      integer ifn(nion),iln(nion),ifw(nss),nshift
      integer nfiles,nat_ion,nat_w,ll,idi,imol
      integer ichain,iat,itype,itmp,iflag,i,jj,icross1,icross2,ipt
      integer ifiles,natoms,isnap,j,ipos,ich,nshots,  &
       npos,ik,i1,j1,iii,istart,Oi,Hi(2),h2o_type
      integer ii,h1_1,h1_2,h2_1,h2_2,kk,o1,o2,num_r,loci
      integer id_d(10),id_a(10),ndch,nach,id_di(10),id_ai(10),id_w(3)
      integer id_dc(10),id_dci(10),donated(10),jjj,indx_ss(nss)
      integer indx_p(2000),natp,indx_s(1000),ns,hb(500,2,20000)
      integer knt_hb(20000),nat,knt,typei,ns2,indx_s2(1000),natps
      integer indx_ps(2000),iend,imax,jmax
      character(len=256) fname(MaxInpFiles),fname1,fname2
      character(len=10) str1
      character(len=5)typea
      real*4 dis_oo,a_oho,dis_c,a_c,h_bond_mat(10,10),d1,dd2
      real*4 kount_mat(maxnch,2),N_D,N_A,doo(3),shift_xyz(3),ai
      real*4 xun1(3,maxat),xun(3,maxat),pi,kount,dx,dy,dz
      real*4 voh(3),box(3),r_min,r_max,r_bin,mass(100),mass_tot
      real*4 kount_d(500),kount_ch(500),r_array(500),com(3)
      real*4 v_h21_o1(3),v_h21_o2(3),v_h22_o1(3),v_h22_o2(3)
      real*4 kount_a(500),ang,rox,rxo,rxh1,rxh2,roh1,roh2
      real*4 vox(3),vxo(3),vxh1(3),vxh2(3),voh1(3),voh2(3)    
      real*4 disij,hba,hbd,dr(3),voo(3),vh1o(3),vh2o(3),ro,tt
      real*4 rh1o,rh2o,xlo,xhi,ylo,yhi,zlo,zhi,chgi,dchg,deltt
      real*4 r_max2,dr2(3),vj(3),acfi,nass,ttime,sumi
      integer ki,kj,jstart,jend,nhb0,nhb,k1,k2,knti_tot,used
      integer,allocatable::acf(:,:)
!
!
!    ***  initializations
!
!    Files for analysis
!
      knt_hb(:)=0_4
      r_max=5.0
      nat_w=3
      deltt=2.0
      open(12,file='lammps.dat')
      read(12,*) str1
      read(12,*) nat,str1
      write(6,*) 'nat',nat
      open(16,file='acf.dat')
      open(13,file='test2')
      natp=0
      do while (1>0)
          read(12,*) str1
          write(6,*) str1
          if (str1== "Atoms") then
              !read(12,*) str1
              goto 2000
          endif
      enddo
2000  continue
      knt=0
      natps=0
      do ii=1,nat
          read(12,*)   idi,imol,typei,chgi
          if (imol<=2) then
              natp=natp+1
              indx_p(natp)=ii
          endif
          if (typei==13 .and. abs(chgi+1.0424)< 0.00001) then
              knt=knt+1
              indx_ss(knt)=ii
              write(6,*) ii
          endif
      enddo
      write(6,*) 'knt',knt
      allocate(acf(knt,maxshots))
      acf(:,:)=0_4
      dis_c=3.45
      a_c=30.0
      pi=3.141592654
      open(22,file="files",status="old")
      read(22,*) nfiles
      do i=1,nfiles
        read(22,'(a160)') fname(i)
        write(6,'(a160)') fname(i)
        open (30+i,file=fname(i),status='old')
      enddo
      write(6,*) "Opened ",nfiles
      close(22)
!       
!
      ipos=0
      iii=0
      kount_ch(:)=0_8
      do ifiles = 31,30+nfiles
         do while ( 1>0) !.not. IS_IOSTAT_END(ifiles))
         !do while (.not. feof(ifiles))
            read(ifiles,'(A10)',END=200) str1
            do ii=1,4
                read(ifiles,*) str1
            enddo
            read(ifiles,*) xlo,xhi
            read(ifiles,*) ylo,yhi
            read(ifiles,*) zlo,zhi
            read(ifiles,*)   str1
            box(1)=xhi-xlo
            box(2)=yhi-ylo
            box(3)=zhi-zlo
            write(6,*) str1
            do ii=1,nat
               read(ifiles,*) idi,imol, &
       xx(1,ii),xx(2,ii),xx(3,ii)    
               xx(1,ii)=xx(1,ii)-xlo
               xx(2,ii)=xx(2,ii)-ylo
               xx(3,ii)=xx(3,ii)-zlo 
            enddo
            rboxgen(1)=real(box(1))
            rboxgen(2)=real(box(2))
            rboxgen(3)=real(box(3))
            !write(6,*) rboxgen
            do j=1,3
               el2gen(j)=0.5*rboxgen(j)
            enddo
            ipos=ipos+1
            write(6,*) "ipos= ",ipos
            indx_s(:)=0_4
            do ii=1,knt
               do jj=1,natp
                  do kk=1,3
                     dr(kk)=xx(kk,indx_ss(ii))-xx(kk,indx_p(jj))
                     if (dr(kk)>el2gen(kk)) then
                         dr(kk)=dr(kk)-rboxgen(kk)
                     elseif (dr(kk)<-el2gen(kk)) then
                         dr(kk)=dr(kk)+rboxgen(kk)    
                     endif
                  enddo
                  disij=sqrt(dr(1)**2+dr(2)**2+dr(3)**2)
                  used=0
                  if (disij<r_max) then
                      acf(ii,ipos)=1 
                      goto 1000 
                  endif
               enddo
1000           continue               
            enddo
!  
!
         enddo  !  do while
200      continue
      enddo  ! do  ii=1,nfiles
      write(6,*) 'kng',knt,ipos
      imax=ipos-1
      do ii=1,imax
          ns=0 
          ns2=0
          jmax=ipos-ii
          !write(6,*) 'here2',imax,jmax
          do jj=1,jmax
              istart=jj
              iend=jj+ii
              do kk=1,knt
                  
                  if (acf(kk,istart)==1) then
                      ns=ns+1   
                      do ll=istart,iend 
                        if (acf(kk,ll)==0)  goto 300 
                      enddo
300                   continue                      
                      if (ll==iend) ns2=ns2+1
                  endif
              enddo
          enddo
          acfi=dble(ns2)/dble(ns)
          ttime=dble(ii)*2.0
          write(16,*) ttime,acfi 
      enddo
!           
      end  PROGRAM hbond_all
