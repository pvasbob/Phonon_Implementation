!------------------------------------------------------------------------------
! pnfam_solver.f90
!
! M.T. Mustonen and T. Shafer, UNC Chapel Hill, 2013-15
!------------------------------------------------------------------------------
module pnfam_solver
   use logger
   use blockmatrix_type
   use external_field_type
   use density_set_type
   implicit none
   integer, parameter, private :: dp = kind(1d0)
   real(dp), parameter :: pi = 3.141592653589793238462643_dp
   
   type(blockmatrix), private :: greenXre, greenXim, greenYre, greenYim, f20, f02,              &
      h20re, h20im, h02re, h02im, rex, imx, rey, imy, rerho_np, imrho_np, rerho_pn, imrho_pn,   &
      rekp, imkp, rekm, imkm, redp, imdp, redm, imdm, reh_np, imh_np, reh_pn, imh_pn
   type(blockmatrix), dimension(:), allocatable, private :: g20, g02
   type(density_set), private :: rhox
   
   ! Finite temperature and odd-mass nuclei
   type(blockmatrix), private :: greenPre, greenPim, greenQre, greenQim, f11, f11t, &
         h11re, h11im, h11tre, h11tim, rep, imp, req, imq, mat_fbfa_p, mat_fbfa_q,  &
         mat_1fafb_x, mat_1fafb_y
   type(blockmatrix), dimension(:), allocatable, private :: g11, g11t
   
   integer, private :: nxy
   
   !----------------------------------QQst      modify18
   character(100), public,  parameter :: qq_io_path = "/21dayscratch/scr/q/u/qunqun/pn_io_1/"
   !----------------------------------QQend     modify18
contains
   
   !---------------------------------------------------------------------------
   ! Subroutine for solving the FAM equations for one value of energy. Before
   ! calling this the first time, the module data structures must be initialized
   ! with init_pnfam_solver.
   !---------------------------------------------------------------------------
   subroutine pnfam_solve(omega, width, Ep, En, vp, up, vn, un, max_iter, &
      convergence_epsilon, quench, bminus, str, iter_used, cstr)
      use constants, only : iu
      use broyden_mixer
      use hamiltonian
      use hfbtho_basis, only : ft_active, ft_temp, hfb_blo_active
      !----------------------
      use pnfam_WX
      !----------------------
      implicit none
      real(dp), intent(in) :: omega, width, Ep(:), En(:), convergence_epsilon
      integer, intent(in) :: max_iter
      type(blockmatrix), intent(in) :: vp, vn, up, un
      logical, intent(in) :: bminus
      integer, intent(out) :: iter_used
      real(dp),    dimension(:), intent(out) :: str
      complex(dp), dimension(:), intent(out), optional :: cstr ! complex strength
      
      ! Finite temperature
      real(dp)    :: ft_real_factor
      complex(dp) :: ft_cmplx_factor
      
      integer :: iter, ixterm
      logical :: use_diagonal_blocks

      ! Adjustments to solution
      real(dp), intent(in) :: quench
      !------------------------------------QQst         modify19
      ! variable for assigning X,Y to matX and matY
      !---------------------------------------
     !integer, parameter :: dp = kind(1.0d0)
      integer, parameter :: dc = kind(cmplx(1.0,1.0,kind=dp))
      integer, allocatable :: Ph_readin_select(:,:)
      integer, allocatable :: U_select(:,:),V_select(:,:)
      integer, allocatable :: X_select(:,:),Y_select(:,:)
      integer, allocatable :: Ph_db(:)
      integer :: Ph_nb, Ph_NumberofX
      integer :: Ph_sum_db
      integer, allocatable :: Ph_sum_db_colm1(:)          !!!Q: dimension is iHFB_NB, record sum of nd of first column minus 1 blocks.
      Complex(dc), allocatable  :: Ph_matX(:,:),Ph_matY(:,:)
      Complex(dc), allocatable  :: Ph_X(:), Ph_Y(:)       !!!Q: maybe not necessary.
      !---------------------------------------
      ! variable for assigning matN11 and matP11
      !---------------------------------------
      !
      integer, allocatable :: WX_readin_select(:,:)
      integer, allocatable :: WX_id(:)
      integer :: WX_iHFB_NB, WX_NumberofPH11,WX_NumberofPH11_maxval
      integer :: WX_sum_nd
      integer, allocatable :: WX_sum_nd_colm1(:)          !!!Q: dimension is iHFB_NB, record sum of nd of first column minus 1 blocks.
      real(dp) :: WX_realomega 
      real(dp) :: W_omega 
      !---------------------------------------
      Complex(Kind(1.000D0)), allocatable  :: WX_NH11_E101(:),WX_PH11_E101(:)
      Complex(Kind(1.000D0)), allocatable  :: WX_BNH11_E101(:),WX_BPH11_E101(:)  !!!Q: bmod01
      !
      Complex(Kind(1.000D0)), allocatable  :: WX_matN11_E101(:,:),WX_matP11_E101(:,:)
      Complex(Kind(1.000D0)), allocatable  :: WX_matBN11_E101(:,:),WX_matBP11_E101(:,:)  !!!Q: bmod01
      !
      Complex, allocatable  :: WX_matpnbasisN11_E101(:,:),WX_matpnbasisP11_E101(:,:) !!!Q:bmod01    
     !real(dp), allocatable  :: WX_matpnbasisN11_E101(:,:),WX_matpnbasisP11_E101(:,:)
      integer, allocatable :: WX_pnbasis_select(:,:)   !modifyX01
      !---------------------------------------
 
      ! variable used by both matX, matY and matN11, matP11
      !---------------------------------------
      integer :: sum_block
      integer :: qqj, qqi,qqk, ct1, ct2, ct3      !!!Q: ct means count, for any
      integer :: inbcc                        !!!Q: inside block column count, used to decide which column to be assigned.
      integer :: inbrc                        !!!Q: inside block row count, used to decide whcih row element to be assigned.
      integer :: inbrr                        !!!Q: used to record how many elements already assigned to a row inside a block.
      logical :: Ph_test=.True.               !!!Q: logical to open or close phonon calculation.
      complex,parameter :: cunit=cmplx(0.0_dp,1.0_dp,kind=dp)
      complex,parameter :: runit=cmplx(1.0_dp,0.0_dp,kind=dp)
      integer :: nph, nph_rec   !!! number of phonons 

     !Complex(Kind(1.000D0)), allocatable  :: WX_NH11_store(:,:),WX_PH11_store(:,:)
     !Complex(Kind(1.000D0)), allocatable  :: WX_BNH11_store(:,:),WX_BPH11_store(:,:)  !!!Q: bmod01
      real(dp), allocatable :: WX_realomega_store(:) 
      integer, allocatable :: WX_iHFB_NB_store(:), WX_NumberofPH11_store(:)
      integer, allocatable :: WX_id_store(:,:)
      integer, allocatable :: WX_readin_select_store(:,:,:)
      integer, allocatable :: iso_store(:), L_store(:), k_store(:)



      Complex(Kind(1.000D0)), allocatable  :: WX_NH11_store1D(:)
      Complex(Kind(1.000D0)), allocatable  :: WX_BNH11_store1D(:)
      Complex(Kind(1.000D0)), allocatable  :: WX_PH11_store1D(:)
      Complex(Kind(1.000D0)), allocatable  :: WX_BPH11_store1D(:)


      !---------------------------------------------
      ! matrix PX, PPX,  
      !---------------------------------------------
      !interm variable matrix, used for both WX and WY
   !  Complex(Kind(1.000D0)), allocatable  :: Ph_matPdF(:,:),Ph_matPdFg(:,:), Ph_matFNc(:,:),Ph_matFNcg(:,:)
      ! matrix distinguished.
      !---------------------------------------------
   !  Complex(Kind(1.000D0)), allocatable  :: Ph_matPPdXg(:,:), Ph_matXNcgN(:,:), Ph_matPXNcg(:,:), Ph_matPdXgN(:,:)
   !  Complex(Kind(1.000D0)), allocatable  :: Ph_matPPdYg(:,:), Ph_matYNcgN(:,:), Ph_matPYNcg(:,:), Ph_matPdYgN(:,:)
   !  Complex(Kind(1.000D0)), allocatable  :: Ph_matWX(:,:), Ph_matWY(:,:)             !modify32
      !---------------------------------------------
      real, allocatable :: Ph_ReWX(:),Ph_ImWX(:)  
      real, allocatable :: Ph_ReWY(:),Ph_ImWY(:)  
      !---------------------------------------------
      Complex, allocatable  :: WgreenX(:), WgreenY(:)                                             
      !---------------------------------------------
 
      real :: Ph_mult=1.0
      real :: Delta = +0.1
      !-------------------------QQst    modify30
      real :: si_mini=100.0
      integer :: iter_start, iter_end
      integer, parameter :: iter_endmax = 100
      !-------------------------QQend   modify30
      complex, dimension(1) :: cstr_mini
      !--------------------------
      integer :: which_iso
      integer :: which_L
      integer :: which_k
      integer :: which_ph
      character(len=10)  :: iso_id
      character(len=10)  :: L_id
      character(len=10)  :: k_id
      character(len=10)  :: ph_id
      logical :: file_exists


      open(48,file='Ph_converge',access='append')    ! modify18
     !open(49,file='Ph_from_simplex_to_pnfam_dH11')
      open(50,file='50_checke')    ! modify18
     !
      !------------------------------------QQend
      !
      ! The "diagonal" blocks of the generalized density (P and -Q) are only
      ! active when starting from a finite-temperature or odd-mass HFB solution
      if (ft_active .or. hfb_blo_active) then
         use_diagonal_blocks = .true.
      else
         use_diagonal_blocks = .false.
      end if
      
      ! Cannot handle both T>0 and odd A
      if (hfb_blo_active .and. ft_active) then
         write(*,'(a)') 'ERROR: this code cannot handle T>0 and odd-Z/odd-N simultaneously.'
         stop
      end if
      
      
      ! Initial guess for X and Y
      ! We keep the real and complex part of X and Y separate to allow
      ! optimizing various matrix operations, where often only one of the
      ! matrices is complex and lifting the real ones to be complex is a waste
      ! of computing time
      ReX%elem = 0; ImX%elem = 0 ; ReY%elem = 0 ; ImY%elem = 0
      !-------------------------------------QQst         modify19
     !write(33,*) 'size(ReX%elem)',' ', size(ReX%elem)
     !write(33,*) 'size(En)',' ',size(En)
     !write(33,*) 'size(Ep)',' ',size(Ep)
     !write(33,*) 'omega',' ',omega,' ','width',' ',width,' ','dp',' ',dp
     !write(33,*) 'cunit',' ',17.26*cunit,' ','runit',' ',17.26*runit
      !-------------------------------------QQend   
      if (use_diagonal_blocks) then
         ReP%elem = 0; ImP%elem = 0 ; ReQ%elem = 0 ; ImQ%elem = 0
      end if
      
      ! Compute the denominators in the FAM equations with minus sign attached,
      ! e.g. -1/(E_p + E_n - \omega)
      if (bminus) then
         call greenx(Ep, En, -cmplx(omega,width,dp), greenXre, greenXim)
         call greenx(Ep, En,  cmplx(omega,width,dp), greenYre, greenYim)
         
         if (use_diagonal_blocks) then
            call greenx(Ep, -En, -cmplx(omega,width,dp), greenPre, greenPim)
            call greenx(Ep, -En,  cmplx(omega,width,dp), greenQre, greenQim)            
         end if
      else
         call greenx(En, Ep, -cmplx(omega,width,dp), greenXre, greenXim)
         call greenx(En, Ep,  cmplx(omega,width,dp), greenYre, greenYim)
         
         if (use_diagonal_blocks) then
            call greenx(En, -Ep, -cmplx(omega,width,dp), greenPre, greenPim)
            call greenx(En, -Ep,  cmplx(omega,width,dp), greenQre, greenQim)            
         end if
      end if
   
      ! The iteration loop
      si = 1
      iter_used = -1  ! in case of no convergence, this value will remain
      !---------------------------------------------------------------QQst  modify20
      ! allocate for Phonon calculation.
      !---------------------
      if(Ph_test) then			!!!Q: prepare before iter loop
          !--------------------------------
	  !--------------------------------
	  ! readin binfo from exterfield, assign array X and Y to matX and matY
	  !--------------------------------
	  open(300,file='Ph_binary_output_extfield_GT0',status='unknown',form='unformatted')
	  Read(300) Ph_nb
	  Read(300) Ph_NumberofX                    !!!Q: read
	  !
	  if(allocated(Ph_db)) deallocate(Ph_sum_db_colm1,Ph_db)
	  Allocate(Ph_sum_db_colm1(Ph_nb),Ph_db(Ph_nb))
	  ! 
	  Read(300) Ph_db                 !!!Q: pnfam_solver.f90 does't have db(:)
	  Read(300) Ph_sum_db
	  Read(300) Ph_sum_db_colm1
	  !
	  !
	  if(allocated(Ph_readin_select)) deallocate(Ph_readin_select)
	  Allocate(Ph_readin_select(Ph_nb,Ph_nb))
	  !
	  Read(300) Ph_readin_select                   !!!Q: read
	  close(300)

	  Allocate(U_select(Ph_nb,Ph_nb), V_select(Ph_nb,Ph_nb), &
                   X_select(Ph_nb,Ph_nb), Y_select(Ph_nb,Ph_nb))

         U_select = 0
         V_select = 0
         X_select = 0
         Y_select = 0
         do qqi = 1, Ph_nb  
             do qqj = 1, Ph_nb
                 if(qqi.eq.(qqj-Ph_nb/2)) then
                     V_select(qqi,qqj) = 1
                 end if
                 if((qqi-Ph_nb/2).eq.qqj) then
                     V_select(qqi,qqj) = 1
                 end if
             end do 
         end do
         !
         X_select= matmul(Ph_readin_select, V_select)
         Y_select= matmul(V_select, Ph_readin_select)
         open(20,file='X_select')
         do qqi = 1, Ph_nb
            do qqj = 1, Ph_nb
              write(20,*)qqi, qqj, U_select(qqi,qqj), V_select(qqi,qqj), &
                        X_select(qqi,qqj)
            end do
         end do
         close(20)


          !-----------------------
          ! Wgreens         
          !----------------------
	  if(allocated(WgreenX)) deallocate(WgreenX, WgreenY)
	  Allocate(WgreenX(Ph_sum_db**2), WgreenY(Ph_sum_db**2))
          WgreenX = 0
          WgreenY = 0
	  !
         !Call WgreenX_func(Ph_sum_db, Ep, En, omega, WX_realomega, width, WgreenX) 
         !Call WgreenY_func(Ph_sum_db, Ep, En, omega, WX_realomega, width, WgreenY) 
          !
          Allocate(Ph_ReWX(Ph_NumberofX), Ph_ImWX(Ph_NumberofX))
          Ph_ReWX=0
          Ph_ImWX=0
 
          Allocate(Ph_ReWY(Ph_NumberofX), Ph_ImWY(Ph_NumberofX))
          Ph_ReWY=0
          Ph_ImWY=0
 
	  !----------------------------
	  ! preserve below to check
	  !----------------------------
	 !sum_block=0
	 !do qqj=1,Ph_nb
	 !   do qqi=1,Ph_nb
	!if(Ph_readin_select(qqj,qqi)==1) then 
	!   sum_block=sum_block+Ph_db(qqj)*Ph_db(qqi)             !!!Q: number of elements sum block 
	!  !write(34,*) qqj, qqi, Ph_db(qqj), Ph_db(qqi), sum_block
	!end if
	!    end do
	! end do
          !
          !-----------------------------
          ! allocate Ph_matX and Ph_matY 
          !-----------------------------
          Allocate(Ph_matX(Ph_sum_db,Ph_sum_db),Ph_matY(Ph_sum_db,Ph_sum_db))
          Ph_matX=0
          Ph_matY=0



         nph = 0
         !=================================================
         !=================================================
         ! find out how many phonons, and allocate storage memory.
         !=================================================
        do which_iso=  0,1
           write(iso_id,'(i0)') which_iso
             do which_L = 0,7
               write(L_id,'(i0)') which_L
               do which_k = 0, which_L
                 do which_ph = 1,15
                   write(k_id,'(i0)') which_k
                   write(ph_id,'(i0)') which_ph
                   !if(which_k.eq.0) then
                   !   Ph_mult = 1                   !!!Q: apply for deformed nuclear.          
                   !else
                   !   Ph_mult = 2  
                   !end if
                 
                

	  inquire(file=trim(adjustl(qq_io_path))//'QQ_binary_output_binfo_'//trim(adjustl(iso_id))//trim(adjustl(L_id))//trim(adjustl(k_id))//trim(adjustl(ph_id)),exist=file_exists)    !modify18
      if(file_exists) then
         nph = nph + 1
      endif
    
                 end do
               end do
             end do
        end do

        nph_rec = nph

        write(50, *) 'nph:', nph


	  if(allocated(WX_realomega_store)) deallocate(WX_realomega_store, WX_iHFB_NB_store, WX_NumberofPH11_store)
	  allocate(WX_realomega_store(nph), WX_iHFB_NB_store(nph), WX_NumberofPH11_store(nph))
      if(allocated(iso_store)) deallocate(iso_store, L_store, k_store)
      allocate(iso_store(nph), L_store(nph), k_store(nph))

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        nph = 0
        do which_iso=  0,1
        write(iso_id,'(i0)') which_iso
        do which_L = 0,7
        write(L_id,'(i0)') which_L
        do which_k = 0, which_L
           do which_ph = 1,15
              write(k_id,'(i0)') which_k
              write(ph_id,'(i0)') which_ph
                if(which_k.eq.0) then
                   Ph_mult = 1                   !!!Q: apply for deformed nuclear.          
                else
                   Ph_mult = 2  
                end if
                

	  inquire(file=trim(adjustl(qq_io_path))//'QQ_binary_output_binfo_'//trim(adjustl(iso_id))//trim(adjustl(L_id))//trim(adjustl(k_id))//trim(adjustl(ph_id)),exist=file_exists)    !modify18
      if(file_exists) then
      nph = nph  + 1
	  open(110,file=trim(adjustl(qq_io_path))//'QQ_binary_output_binfo_'//trim(adjustl(iso_id))//trim(adjustl(L_id))//trim(adjustl(k_id))//trim(adjustl(ph_id)),status='unknown',form='unformatted')       !modify18
          Read(110) WX_realomega
	  Read(110) WX_iHFB_NB                         !!!Q: read


	  if(allocated(WX_id)) deallocate(WX_id, WX_sum_nd_colm1)
	  allocate(WX_id(WX_iHFB_NB), WX_sum_nd_colm1(WX_iHFB_NB))

          if(allocated(WX_readin_select))  deallocate(WX_readin_select)
	  Allocate(WX_readin_select(WX_iHFB_NB,WX_iHFB_NB))


	  Read(110) WX_NumberofPH11                    !!!Q: read
	  Read(110) WX_id                              !!!Q: read
	  Read(110) WX_readin_select                   !!!Q: read
	  !
	  close(110)
 
      WX_realomega_store(nph) = WX_realomega
      WX_iHFB_NB_store(nph) = WX_iHFB_NB
      WX_NumberofPH11_store(nph) = WX_NumberofPH11
 
      iso_store(nph) =  which_iso
      L_store(nph) =  which_L 
      k_store(nph) =  which_k 

      end if    ! file_exists
         end do ! which_ph
         end do ! which_k
         end do ! which_L
         end do ! which_iso

    WX_NumberofPH11_maxval = maxval(WX_NumberofPH11_store)
    write(50,*) 'WX_realomega_store' ,WX_realomega_store
    write(50,*) 'WX_iHFB_NB_store' ,WX_iHFB_NB_store 
    write(50,*) 'WX_NumberofPH11_store',WX_NumberofPH11_store
    write(50,*) 'WX_NumberofPH11_maxval',WX_NumberofPH11_maxval


	  if(allocated(WX_id_store)) deallocate(WX_id_store, WX_readin_select_store)
	  allocate(WX_id_store(nph, WX_iHFB_NB_store(1)), WX_readin_select_store(nph, WX_iHFB_NB_store(1),WX_iHFB_NB_store(1)))


	 !if(allocated(WX_NH11_store)) deallocate(WX_NH11_store, WX_BNH11_store, WX_PH11_store, WX_BPH11_store)
	 !allocate(WX_NH11_store(nph, WX_NumberofPH11_maxval), WX_BNH11_store(nph, WX_NumberofPH11_maxval), & 
     !          WX_PH11_store(nph, WX_NumberofPH11_maxval), WX_BPH11_store(nph, WX_NumberofPH11_maxval))

      if(allocated(WX_NH11_store1D)) deallocate(WX_NH11_store1D,WX_BNH11_store1D,WX_PH11_store1D,WX_BPH11_store1D)
      allocate(WX_NH11_store1D(nph*WX_NumberofPH11_maxval))
     !allocate(WX_BNH11_store1D(nph*WX_NumberofPH11_maxval))
      allocate(WX_PH11_store1D(nph*WX_NumberofPH11_maxval))
     !allocate(WX_BPH11_store1D(nph*WX_NumberofPH11_maxval))

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! loop phonon binfo once more to populate WX_id_store,
    ! WX_readin_select_store,  WX_?H11_store
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        nph = 0
        do which_iso=  0,1
        write(iso_id,'(i0)') which_iso
        do which_L = 0,7
        write(L_id,'(i0)') which_L
        do which_k = 0, which_L
           do which_ph = 1,15
              write(k_id,'(i0)') which_k
              write(ph_id,'(i0)') which_ph
                if(which_k.eq.0) then
                   Ph_mult = 1                   !!!Q: apply for deformed nuclear.          
                else
                   Ph_mult = 2  
                end if
                

	  inquire(file=trim(adjustl(qq_io_path))//'QQ_binary_output_binfo_'//trim(adjustl(iso_id))//trim(adjustl(L_id))//trim(adjustl(k_id))//trim(adjustl(ph_id)),exist=file_exists)    !modify18
      if(file_exists) then

      nph = nph + 1

	  open(110,file=trim(adjustl(qq_io_path))//'QQ_binary_output_binfo_'//trim(adjustl(iso_id))//trim(adjustl(L_id))//trim(adjustl(k_id))//trim(adjustl(ph_id)),status='unknown',form='unformatted')       !modify18
          Read(110) WX_realomega
	      Read(110) WX_iHFB_NB                         !!!Q: read


	      if(allocated(WX_id)) deallocate(WX_id, WX_sum_nd_colm1)
	      allocate(WX_id(WX_iHFB_NB), WX_sum_nd_colm1(WX_iHFB_NB))

          if(allocated(WX_readin_select))  deallocate(WX_readin_select)
	  Allocate(WX_readin_select(WX_iHFB_NB,WX_iHFB_NB))


	  Read(110) WX_NumberofPH11                    !!!Q: read
	  Read(110) WX_id                              !!!Q: read
	  Read(110) WX_readin_select                   !!!Q: read
	  !
	  close(110)






          if(allocated(WX_NH11_E101)) deallocate(WX_NH11_E101, & 
                                                 WX_PH11_E101)
	  Allocate(WX_NH11_E101(WX_NumberofPH11),WX_PH11_E101(WX_NumberofPH11))
	  WX_NH11_E101=0
	  WX_PH11_E101=0

          !--------------------------------------------------
    !     if(allocated(WX_BNH11_E101)) deallocate(WX_BNH11_E101, & 
    !                                            WX_BPH11_E101)
	! Allocate(WX_BNH11_E101(WX_NumberofPH11),WX_BPH11_E101(WX_NumberofPH11))
	! WX_BNH11_E101=0
	! WX_BPH11_E101=0
          !--------------------------------------------------

          !----------------------------------------E101
          ! 
	  open(201,file=trim(adjustl(qq_io_path))//'QQ_binary_NH11_'//trim(adjustl(iso_id))//trim(adjustl(L_id))//trim(adjustl(k_id))//trim(adjustl(ph_id)),& 
                                                                status='unknown',form='unformatted')        !modify18
	  Read(201) WX_NH11_E101
	! Read(201) WX_BNH11_E101
	  close(201)
	  !
	  !
	  open(201,file=trim(adjustl(qq_io_path))//'QQ_binary_PH11_'//trim(adjustl(iso_id))//trim(adjustl(L_id))//trim(adjustl(k_id))//trim(adjustl(ph_id)),&
                                                                status='unknown',form='unformatted')        !modify18
	  Read(201) WX_PH11_E101
	! Read(201) WX_BPH11_E101
	  close(201)
          ! 

     WX_id_store(nph,:) = WX_id
     WX_readin_select_store(nph,:,:) = WX_readin_select
!
    do qqi = 1, size(WX_NH11_E101)
     ! WX_NH11_store(nph, qqi) = WX_NH11_E101 (qqi)
     ! WX_BNH11_store(nph,qqi) = WX_BNH11_E101(qqi)  
     ! WX_PH11_store(nph, qqi) = WX_PH11_E101(qqi)
     ! WX_BPH11_store(nph,qqi) = WX_BPH11_E101(qqi)


       WX_NH11_store1D((nph-1)*WX_NumberofPH11_maxval + qqi) = WX_NH11_E101(qqi)
     ! WX_BNH11_store1D((nph-1)*WX_NumberofPH11_maxval + qqi) = WX_BNH11_E101(qqi)
       WX_PH11_store1D((nph-1)*WX_NumberofPH11_maxval + qqi) = WX_PH11_E101(qqi)
     ! WX_BPH11_store1D((nph-1)*WX_NumberofPH11_maxval + qqi) = WX_BPH11_E101(qqi)
    enddo
!
    


      end if    ! file_exists
         end do ! which_ph
         end do ! which_k
         end do ! which_L
         end do ! which_iso


   !do qqi = 1, size(WX_NH11_E101)
   !  if(WX_NH11_store(nph, qqi).Ne. WX_NH11_E101(qqi) .and. &  
   !    WX_NH11_store(nph, qqi).Ne. WX_NH11_E101(qqi) .and.&
   !    WX_NH11_store(nph, qqi).Ne. WX_NH11_E101(qqi) .and.& 
   !    WX_NH11_store(nph, qqi).Ne. WX_NH11_E101(qqi))  then 
   !      write(50,*) 'Wrong element assignment in WX_NH11_store(nph))'
   !  endif
   !enddo
   !write(50,*) 'WX_NH11_store(nph) has same element as WX_NH11_E101'
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         ! 
      end if			!!!Q: preparation before iter loop end.
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! end preparation
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
      ! 
      !-------------------------------------------------------------------------QQend 
      ! iter loop begin
      !--------------------------------------
      ! 
      do iter=0,max_iter
         if (abs(quench) < 1e-10) then
            h20re%elem = f20%elem ; h20im%elem = 0  ! only the external field
            h02re%elem = f02%elem ; h02im%elem = 0
            
            if (use_diagonal_blocks) then
               h11re%elem  = f11%elem  ; h11im%elem  = 0
               h11tre%elem = f11t%elem ; h11tim%elem = 0
            end if
         else
         ! Compute the perturbed densities
         if (bminus) then
            ! rho_np = -Un X(-)' Vp' + Vn Y(-)' Up'
            call triprod('n', un, 't', rex, 't', vp, -1d0, 0d0, rerho_np)
            call triprod('n', un, 't', imx, 't', vp, -1d0, 0d0, imrho_np)
            call triprod('n', vn, 't', rey, 't', up, 1d0, 1d0, rerho_np)
            call triprod('n', vn, 't', imy, 't', up, 1d0, 1d0, imrho_np)
            
            ! rho_pn = Up X(-) Vn' - Vp Y(-) Un'
            call triprod('n', up, 'n', rex, 't', vn, 1d0, 0d0, rerho_pn)
            call triprod('n', up, 'n', imx, 't', vn, 1d0, 0d0, imrho_pn)
            call triprod('n', vp, 'n', rey, 't', un, -1d0, 1d0, rerho_pn)
            call triprod('n', vp, 'n', imy, 't', un, -1d0, 1d0, imrho_pn)
            
            ! kappa(+)_np = -un X(-)' up' + vn Y(-)' vp'
            call triprod('n', un, 't', rex, 't', up, -1d0, 0d0, rekp)
            call triprod('n', un, 't', imx, 't', up, -1d0, 0d0, imkp)
            call triprod('n', vn, 't', rey, 't', vp, 1d0, 1d0, rekp)
            call triprod('n', vn, 't', imy, 't', vp, 1d0, 1d0, imkp)
            
            ! kappa(-)_np = vn X(-)*' vp' - un Y(-)'* up'
            call triprod('n', vn, 't', rex, 't', vp, 1d0, 0d0, rekm)
            call triprod('n', vn, 't', imx, 't', vp, -1d0, 0d0, imkm)
            call triprod('n', un, 't', rey, 't', up, -1d0, 1d0, rekm)
            call triprod('n', un, 't', imy, 't', up, 1d0, 1d0, imkm)
            
            if (use_diagonal_blocks) then
               ! rho_pn += up P(-) un' - vp Q(-) vn'
               call triprod('n', up, 'n', rep, 't', un, 1d0, 1d0, rerho_pn)
               call triprod('n', up, 'n', imp, 't', un, 1d0, 1d0, imrho_pn)
               call triprod('n', vp, 'n', req, 't', vn,-1d0, 1d0, rerho_pn)
               call triprod('n', vp, 'n', imq, 't', vn,-1d0, 1d0, imrho_pn)
               
               ! rho_np += Un Q(-)' Up' - Vn P(-)' Vp'
               call triprod('n', un, 't', req, 't', up, 1d0, 1d0, rerho_np)
               call triprod('n', un, 't', imq, 't', up, 1d0, 1d0, imrho_np)
               call triprod('n', vn, 't', rep, 't', vp,-1d0, 1d0, rerho_np)
               call triprod('n', vn, 't', imp, 't', vp,-1d0, 1d0, imrho_np)
               
               ! kappa(+)_np += un Q(-)' vp' - vn P(-)' up'
               call triprod('n', un, 't', req, 't', vp, 1d0, 1d0, rekp)
               call triprod('n', un, 't', imq, 't', vp, 1d0, 1d0, imkp)
               call triprod('n', vn, 't', rep, 't', up,-1d0, 1d0, rekp)
               call triprod('n', vn, 't', imp, 't', up,-1d0, 1d0, imkp)
               
               ! kappa(-)_np += -vn Q(-)'* up' + un P(-)'* vp'
               call triprod('n', vn, 't', req, 't', up,-1d0, 1d0, rekm)
               call triprod('n', vn, 't', imq, 't', up, 1d0, 1d0, imkm)
               call triprod('n', un, 't', rep, 't', vp, 1d0, 1d0, rekm)
               call triprod('n', un, 't', imp, 't', vp,-1d0, 1d0, imkm)
            end if
            
         else
            
            ! rho_np = Un X(+) Vp' - Vn Y(+) Up'
            call triprod('n', un, 'n', rex, 't', vp, 1d0, 0d0, rerho_np)
            call triprod('n', un, 'n', imx, 't', vp, 1d0, 0d0, imrho_np)
            call triprod('n', vn, 'n', rey, 't', up, -1d0, 1d0, rerho_np)
            call triprod('n', vn, 'n', imy, 't', up, -1d0, 1d0, imrho_np)
            
            ! rho_pn = -Up X(+)' Vn' + Vp Y(+)' Un'
            call triprod('n', up, 't', rex, 't', vn, -1d0, 0d0, rerho_pn)
            call triprod('n', up, 't', imx, 't', vn, -1d0, 0d0, imrho_pn)
            call triprod('n', vp, 't', rey, 't', un, 1d0, 1d0, rerho_pn)
            call triprod('n', vp, 't', imy, 't', un, 1d0, 1d0, imrho_pn)
            
            ! kappa(+)_pn = -up X(+)' un' + vp Y(+)' vn'
            call triprod('n', up, 't', rex, 't', un, -1d0, 0d0, rekp)
            call triprod('n', up, 't', imx, 't', un, -1d0, 0d0, imkp)
            call triprod('n', vp, 't', rey, 't', vn, 1d0, 1d0, rekp)
            call triprod('n', vp, 't', imy, 't', vn, 1d0, 1d0, imkp)
            
            ! kappa(-)_pn = vp X(+)*' vn' - up Y(+)'* un'
            call triprod('n', vp, 't', rex, 't', vn, 1d0, 0d0, rekm)
            call triprod('n', vp, 't', imx, 't', vn, -1d0, 0d0, imkm)
            call triprod('n', up, 't', rey, 't', un, -1d0, 1d0, rekm)
            call triprod('n', up, 't', imy, 't', un, 1d0, 1d0, imkm)
            
            if (use_diagonal_blocks) then
               ! rho_np += un P(+) up' - vn Q(+) vp'
               call triprod('n', un, 'n', rep, 't', up, 1d0, 1d0, rerho_np)
               call triprod('n', un, 'n', imp, 't', up, 1d0, 1d0, imrho_np)
               call triprod('n', vn, 'n', req, 't', vp,-1d0, 1d0, rerho_np)
               call triprod('n', vn, 'n', imq, 't', vp,-1d0, 1d0, imrho_np)
               
               ! rho_pn += Up Q(+)' Un' - Vp P(+)' Vn'
               call triprod('n', up, 't', req, 't', un, 1d0, 1d0, rerho_pn)
               call triprod('n', up, 't', imq, 't', un, 1d0, 1d0, imrho_pn)
               call triprod('n', vp, 't', rep, 't', vn,-1d0, 1d0, rerho_pn)
               call triprod('n', vp, 't', imp, 't', vn,-1d0, 1d0, imrho_pn)
               
               ! kappa(+)_pn += up Q(+)' vn' - vp P(+)' un'
               call triprod('n', up, 't', req, 't', vn, 1d0, 1d0, rekp)
               call triprod('n', up, 't', imq, 't', vn, 1d0, 1d0, imkp)
               call triprod('n', vp, 't', rep, 't', un,-1d0, 1d0, rekp)
               call triprod('n', vp, 't', imp, 't', un,-1d0, 1d0, imkp)
               
               ! kappa(-)_pn += -vp Q(+)'* un' + up P(+)'* vn'
               call triprod('n', vp, 't', req, 't', un,-1d0, 1d0, rekm)
               call triprod('n', vp, 't', imq, 't', un, 1d0, 1d0, imkm)
               call triprod('n', up, 't', rep, 't', vn, 1d0, 1d0, rekm)
               call triprod('n', up, 't', imp, 't', vn,-1d0, 1d0, imkm)
            end if            
         end if

         ! Compute the perturbed fields (functional derivatives of the EDF)
         ! in the coordinate space and then the mean fields h_pn and h_np
         call density(rerho_pn, imrho_pn, rekp, imkp, rhox)
         call meanfield(rhox, reh_pn, imh_pn)
         call pairingfield(rhox, redp, imdp)    ! Delta+ from rho-tilde+(r)
         
         call density(rerho_np, imrho_np, rekm, imkm, rhox)
         call meanfield(rhox, reh_np, imh_np)
         call pairingfield(rhox, redm, imdm)    ! Delta- from rho-tilde-(r)
         
         ! Transform the Hamiltonian into quasiparticle basis
         if (bminus) then
            ! H20 = up' h(pn) vn - vp' h(np)' un
            !     + up' Delta+(pn) un - vp' Delta-* vn
            call triprod('t', up, 'n', reh_pn, 'n', vn, 1d0, 0d0, h20re)
            call triprod('t', up, 'n', imh_pn, 'n', vn, 1d0, 0d0, h20im)
            call triprod('t', vp, 't', reh_np, 'n', un, -1d0, 1d0, h20re)
            call triprod('t', vp, 't', imh_np, 'n', un, -1d0, 1d0, h20im)
            call triprod('t', up, 'n', redp, 'n', un, 1d0, 1d0, h20re)
            call triprod('t', up, 'n', imdp, 'n', un, 1d0, 1d0, h20im)
            call triprod('t', vp, 'n', redm, 'n', vn, -1d0, 1d0, h20re)
            call triprod('t', vp, 'n', imdm, 'n', vn, 1d0, 1d0, h20im)
            
            ! H02 = up' h(np)' vn - vp' h(pn) un
            !     + up' Delta-* un - vp' Delta+ vn
            call triprod('t', up, 't', reh_np, 'n', vn, 1d0, 0d0, h02re)
            call triprod('t', up, 't', imh_np, 'n', vn, 1d0, 0d0, h02im)
            call triprod('t', vp, 'n', reh_pn, 'n', un, -1d0, 1d0, h02re)
            call triprod('t', vp, 'n', imh_pn, 'n', un, -1d0, 1d0, h02im)
            call triprod('t', up, 'n', redm, 'n', un, 1d0, 1d0, h02re)
            call triprod('t', up, 'n', imdm, 'n', un, -1d0, 1d0, h02im)
            call triprod('t', vp, 'n', redp, 'n', vn, -1d0, 1d0, h02re)
            call triprod('t', vp, 'n', imdp, 'n', vn, -1d0, 1d0, h02im)
            
            if (use_diagonal_blocks) then
               ! H11  =  up' h(pn) un - vp' h(np)' vn + up' d(+)  vn - vp' d(-)*  un
               call triprod('t', up, 'n', reh_pn, 'n', un, 1d0, 0d0, h11re)
               call triprod('t', up, 'n', imh_pn, 'n', un, 1d0, 0d0, h11im)
               call triprod('t', vp, 't', reh_np, 'n', vn,-1d0, 1d0, h11re)
               call triprod('t', vp, 't', imh_np, 'n', vn,-1d0, 1d0, h11im)
               call triprod('t', up, 'n', redp,   'n', vn, 1d0, 1d0, h11re)
               call triprod('t', up, 'n', imdp,   'n', vn, 1d0, 1d0, h11im)
               call triprod('t', vp, 'n', redm,   'n', un,-1d0, 1d0, h11re)
               call triprod('t', vp, 'n', imdm,   'n', un, 1d0, 1d0, h11im)
               
               ! H11~ = -vp' h(pn) vn + up' h(np)' un - vp' d(+)   un + up' d(-)* vn
               call triprod('t', vp, 'n', reh_pn, 'n', vn,-1d0, 0d0, h11tre)
               call triprod('t', vp, 'n', imh_pn, 'n', vn,-1d0, 0d0, h11tim)
               call triprod('t', up, 't', reh_np, 'n', un, 1d0, 1d0, h11tre)
               call triprod('t', up, 't', imh_np, 'n', un, 1d0, 1d0, h11tim)
               call triprod('t', vp, 'n', redp,   'n', un,-1d0, 1d0, h11tre)
               call triprod('t', vp, 'n', imdp,   'n', un,-1d0, 1d0, h11tim)
               call triprod('t', up, 'n', redm,   'n', vn, 1d0, 1d0, h11tre)
               call triprod('t', up, 'n', imdm,   'n', vn,-1d0, 1d0, h11tim)
            end if
            
         else
            
            ! H20 = un' h(np) vn - vn' h(pn)' up
            !     + un' Delta+(np) un - vn' Delta-* vp
            call triprod('t', un, 'n', reh_np, 'n', vp, 1d0, 0d0, h20re)
            call triprod('t', un, 'n', imh_np, 'n', vp, 1d0, 0d0, h20im)
            call triprod('t', vn, 't', reh_pn, 'n', up, -1d0, 1d0, h20re)
            call triprod('t', vn, 't', imh_pn, 'n', up, -1d0, 1d0, h20im)
            call triprod('t', un, 'n', redp, 'n', up, 1d0, 1d0, h20re)
            call triprod('t', un, 'n', imdp, 'n', up, 1d0, 1d0, h20im)
            call triprod('t', vn, 'n', redm, 'n', vp, -1d0, 1d0, h20re)
            call triprod('t', vn, 'n', imdm, 'n', vp, 1d0, 1d0, h20im)
            
            ! H02 = un' h(pn)' vp - vn' h(np) up
            !     + un' Delta-* up - vn' Delta+ vp
            call triprod('t', un, 't', reh_pn, 'n', vp, 1d0, 0d0, h02re)
            call triprod('t', un, 't', imh_pn, 'n', vp, 1d0, 0d0, h02im)
            call triprod('t', vn, 'n', reh_np, 'n', up, -1d0, 1d0, h02re)
            call triprod('t', vn, 'n', imh_np, 'n', up, -1d0, 1d0, h02im)
            call triprod('t', un, 'n', redm, 'n', up, 1d0, 1d0, h02re)
            call triprod('t', un, 'n', imdm, 'n', up, -1d0, 1d0, h02im)
            call triprod('t', vn, 'n', redp, 'n', vp, -1d0, 1d0, h02re)
            call triprod('t', vn, 'n', imdp, 'n', vp, -1d0, 1d0, h02im)
            
            if (use_diagonal_blocks) then
               ! H11  =  un' h(np) up - vn' h(pn)' vp + un' d(+) vp - vn' d(-)* up
               call triprod('t', un, 'n', reh_np, 'n', up, 1d0, 0d0, h11re)
               call triprod('t', un, 'n', imh_np, 'n', up, 1d0, 0d0, h11im)
               call triprod('t', vn, 't', reh_pn, 'n', vp,-1d0, 1d0, h11re)
               call triprod('t', vn, 't', imh_pn, 'n', vp,-1d0, 1d0, h11im)
               call triprod('t', un, 'n', redp,   'n', vp, 1d0, 1d0, h11re)
               call triprod('t', un, 'n', imdp,   'n', vp, 1d0, 1d0, h11im)
               call triprod('t', vn, 'n', redm,   'n', up,-1d0, 1d0, h11re)
               call triprod('t', vn, 'n', imdm,   'n', up, 1d0, 1d0, h11im)
               
               ! H11~ = -vn' h(np) vp + un' h(pn)' up - vn' d(+) up + un' d(-)* vp
               call triprod('t', vn, 'n', reh_np, 'n', vp,-1d0, 0d0, h11tre)
               call triprod('t', vn, 'n', imh_np, 'n', vp,-1d0, 0d0, h11tim)
               call triprod('t', un, 't', reh_pn, 'n', up, 1d0, 1d0, h11tre)
               call triprod('t', un, 't', imh_pn, 'n', up, 1d0, 1d0, h11tim)
               call triprod('t', vn, 'n', redp,   'n', up,-1d0, 1d0, h11tre)
               call triprod('t', vn, 'n', imdp,   'n', up,-1d0, 1d0, h11tim)
               call triprod('t', un, 'n', redm,   'n', vp, 1d0, 1d0, h11tre)
               call triprod('t', un, 'n', imdm,   'n', vp,-1d0, 1d0, h11tim)
            end if
         end if

         ! Quench the residual interaction
         h20re%elem  = h20re%elem(:)*quench ; h20im%elem  = h20im%elem(:)*quench
         h02re%elem  = h02re%elem(:)*quench ; h02im%elem  = h02im%elem(:)*quench

         !------------------------------------------QQst        modify21
         ! assign array X and Y to Ph_matX and Ph_matY to construct product 
         !----------------------------
        if(Ph_test) then

          !--------------------------------
          ! readin binfo, N11 and P11 from lpfam
          !--------------------------------
	  !
        do nph = 1, nph_rec
           if(k_store(nph) .eq. 0) then 
              Ph_mult = 1
           else
              Ph_mult = 1
           end if
     !   do which_iso=  0,1
     !   write(iso_id,'(i0)') which_iso
     !   do which_L = 0,7
     !   write(L_id,'(i0)') which_L
     !   do which_k = 0, which_L
     !      do which_ph = 1,15
     !         write(k_id,'(i0)') which_k
     !         write(ph_id,'(i0)') which_ph
     !           if(which_k.eq.0) then
     !              Ph_mult = 1                   !!!Q: apply for deformed nuclear.          
     !           else
     !              Ph_mult = 2  
     !           end if
     !           
!
	  !inquire(file=trim(adjustl(qq_io_path))//'QQ_binary_output_binfo_'//trim(adjustl(iso_id))//trim(adjustl(L_id))//trim(adjustl(k_id))//trim(adjustl(ph_id)),exist=file_exists)    !modify18
     ! if(file_exists) then
     ! 
     ! nph = nph + 1
    
	 !open(110,file=trim(adjustl(qq_io_path))//'QQ_binary_output_binfo_'//trim(adjustl(iso_id))//trim(adjustl(L_id))//trim(adjustl(k_id))//trim(adjustl(ph_id)),status='unknown',form='unformatted')       !modify18
     !    Read(110) WX_realomega
	 !Read(110) WX_iHFB_NB                         !!!Q: read


	  if(allocated(WX_sum_nd_colm1)) deallocate(WX_sum_nd_colm1)
	  allocate(WX_sum_nd_colm1(WX_iHFB_NB))

     !    if(allocated(WX_readin_select))  deallocate(WX_readin_select)
	 !Allocate(WX_readin_select(WX_iHFB_NB,WX_iHFB_NB))


	 !Read(110) WX_NumberofPH11                    !!!Q: read
	 !Read(110) WX_id                              !!!Q: read
	 !Read(110) WX_readin_select                   !!!Q: read
	  !
	 !close(110)

     W_omega = 1.0
     if(WX_realomega_store(nph) .lt.0.5) then
          write(50,*) iter, 'TLK: ', iso_store(nph), L_store(nph), k_store(nph), 'cycle', WX_realomega_store(nph)
          cycle
     end if
 


         !open(600,file='600_Ph_db_WX_id')
         !write(600, *) Ph_db
         !write(600, *) WX_id_store(nph,:)
         !close(600)
!-------------------------------------------------------------------
          if(allocated(WX_pnbasis_select)) deallocate(WX_pnbasis_select)
          Allocate(WX_pnbasis_select(2*WX_iHFB_NB, 2*WX_iHFB_NB))
          !
	  WX_sum_nd=0
	  WX_sum_nd_colm1=0
	  do qqj=1,WX_iHFB_NB
	      WX_sum_nd=WX_sum_nd + WX_id_store(nph, qqj)
	      WX_sum_nd_colm1(qqj)=WX_sum_nd - WX_id_store(nph, qqj)
	  end do
	  ! 
	  !
          if(allocated(WX_matN11_E101)) deallocate(WX_matN11_E101, &
                                        WX_matP11_E101)
	  Allocate(WX_matN11_E101(WX_sum_nd,WX_sum_nd),WX_matP11_E101(WX_sum_nd,WX_sum_nd))
	  WX_matN11_E101=0
	  WX_matP11_E101=0
          !
          ! 
          !-----------------------------------------------------------
          if(allocated(WX_matBN11_E101)) deallocate(WX_matBN11_E101, &
                                        WX_matBP11_E101)
	  Allocate(WX_matBN11_E101(WX_sum_nd,WX_sum_nd),WX_matBP11_E101(WX_sum_nd,WX_sum_nd))
	  WX_matBN11_E101=0
	  WX_matBP11_E101=0
          !-----------------------------------------------------------
          !
          !
          if(allocated(WX_matpnbasisN11_E101)) deallocate(WX_matpnbasisN11_E101,&
                                                WX_matpnbasisP11_E101) 
	  Allocate(WX_matpnbasisN11_E101(2*WX_sum_nd,2*WX_sum_nd),WX_matpnbasisP11_E101(2*WX_sum_nd,2*WX_sum_nd))
	  WX_matpnbasisN11_E101=0
	  WX_matpnbasisP11_E101=0


!---------------------------------------------------------------------------------



          write(50,*) iter, trim(adjustl(iso_id))//trim(adjustl(L_id))//trim(adjustl(k_id))//trim(adjustl(ph_id)), & 
                 Ph_sum_db, 2*WX_sum_nd, Ph_mult, WX_realomega_store(nph) 

     !    if(allocated(WX_NH11_E101)) deallocate(WX_NH11_E101, & 
     !                                           WX_PH11_E101)
	 !Allocate(WX_NH11_E101(WX_NumberofPH11),WX_PH11_E101(WX_NumberofPH11))
	 !WX_NH11_E101=0
	 !WX_PH11_E101=0

     !    !--------------------------------------------------
     !    if(allocated(WX_BNH11_E101)) deallocate(WX_BNH11_E101, & 
     !                                           WX_BPH11_E101)
	 !Allocate(WX_BNH11_E101(WX_NumberofPH11),WX_BPH11_E101(WX_NumberofPH11))
	 !WX_BNH11_E101=0
	 !WX_BPH11_E101=0
     !    !--------------------------------------------------


          Call WgreenX_func(Ph_sum_db, Ep, En, omega, WX_realomega_store(nph), width, WgreenX) 
          Call WgreenY_func(Ph_sum_db, Ep, En, omega, WX_realomega_store(nph), width, WgreenY) 
 
          do qqi = 1, 2*WX_iHFB_NB
             do qqj = 1, 2*WX_iHFB_NB
                if ((qqi.le.WX_iHFB_NB).and.(qqj.le.WX_iHFB_NB)) then
                        WX_pnbasis_select(qqi, qqj) = WX_readin_select_store(nph, qqi, qqj)
                end if
                !
                if ((qqi.le.WX_iHFB_NB).and.(qqj.gt.WX_iHFB_NB)) then
                        WX_pnbasis_select(qqi, qqj) = WX_readin_select_store(nph, qqi, qqj-WX_iHFB_NB)
                end if
                !
                if ((qqi.gt.WX_iHFB_NB).and.(qqj.le.WX_iHFB_NB)) then
                        WX_pnbasis_select(qqi, qqj) = WX_readin_select_store(nph, qqi-WX_iHFB_NB, qqj)
                end if
                !
                if ((qqi.gt.WX_iHFB_NB).and.(qqj.gt.WX_iHFB_NB)) then
                        WX_pnbasis_select(qqi, qqj) = WX_readin_select_store(nph, qqi-WX_iHFB_NB, qqj-WX_iHFB_NB)
                end if
             end do
          end do
 


          !----------------------------------------E101
          ! 
	  !open(201,file=trim(adjustl(qq_io_path))//'QQ_binary_NH11_'//trim(adjustl(iso_id))//trim(adjustl(L_id))//trim(adjustl(k_id))//trim(adjustl(ph_id)),& 
     !                                                           status='unknown',form='unformatted')        !modify18
	  !Read(201) WX_NH11_E101
	  !Read(201) WX_BNH11_E101
	  !close(201)
	  !!
	  !!
	  !open(201,file=trim(adjustl(qq_io_path))//'QQ_binary_PH11_'//trim(adjustl(iso_id))//trim(adjustl(L_id))//trim(adjustl(k_id))//trim(adjustl(ph_id)),&
     !                                                           status='unknown',form='unformatted')        !modify18
	  !Read(201) WX_PH11_E101
	  !Read(201) WX_BPH11_E101
	  !close(201)
          ! 
	  !-------------------------
	  ! asssign NH11 to matN11
	  ! asssign NP11 to matP11
	  !-------------------------
	  !
	  ct1=1
	  ct3=0
	  inbcc=0
	  inbrc=0
	  inbrr=0
	  do qqj=1,WX_iHFB_NB
	     do qqi=1,WX_iHFB_NB
		if(WX_readin_select_store(nph, qqi,qqj)==1)  then
		   inbcc=WX_sum_nd_colm1(qqj)+1
		   inbrc=WX_sum_nd_colm1(qqi)
		   ct3=ct3+WX_id_store(nph, qqj)*WX_id_store(nph, qqi)
		   do ct2=ct1, ct3 
		      inbrc=inbrc+1
		      inbrr=inbrr+1
                      !----------------------------E101
		      WX_matN11_E101(inbrc, inbcc)=WX_NH11_store1D((nph-1)*WX_NumberofPH11_maxval+ct2)              
		     !WX_matBN11_E101(inbrc, inbcc)=WX_BNH11_store1D((nph-1)*WX_NumberofPH11_maxval+ct2)              
                      !
		      WX_matP11_E101(inbrc, inbcc)=WX_PH11_store1D((nph-1)*WX_NumberofPH11_maxval+ct2)              
		     !WX_matBP11_E101(inbrc, inbcc)=WX_BPH11_store1D((nph-1)*WX_NumberofPH11_maxval+ct2)              
                      !
		      if(inbrr.eq. WX_id_store(nph, qqi)) then
			 inbcc=inbcc+1        !!!Q: col increase by 1
			 inbrc=WX_sum_nd_colm1(qqi)   !!!Q: reset inbrc   
			 inbrr=0              !!!Q: reset inbrr
		      end if
		   end do
		   ct1=ct1+WX_id_store(nph, qqj)*WX_id_store(nph, qqi)
		end if
	     end do
	  end do
          !
          !----------------------------------------------------------
	

          !-----------------------------------QQst      modify22
          do qqj = 1, Ph_sum_db   
             do qqi = 1, Ph_sum_db
                ! (row positive Omega, col positive Omega)
                if((qqi.Le.WX_sum_nd) .and.(qqj.Le.WX_sum_nd)) then 
                        !-----------------------------E101
                        WX_matpnbasisN11_E101(qqi,qqj) = &
                             !0.5*(WX_matN11_E101(qqi,qqj)+WX_matBN11_E101(qqi,qqj))
                              real(WX_matN11_E101(qqi,qqj))*2
                        WX_matpnbasisP11_E101(qqi,qqj) = &
                             !0.5*(WX_matP11_E101(qqi,qqj)+WX_matBP11_E101(qqi,qqj))
                              real(WX_matP11_E101(qqi,qqj))*2
                                 
                end if
                ! (row positive Omega, col negative Omega) 
                if((qqi.Le.WX_sum_nd) .and.(qqj.gt.WX_sum_nd)) then 
                        !-----------------------------E101
                        WX_matpnbasisN11_E101(qqi,qqj) = &
                            ! 0.5*cunit*(WX_matN11_E101(qqi,qqj-WX_sum_nd)-WX_matBN11_E101(qqi,qqj-WX_sum_nd))
                             +aimag(WX_matN11_E101(qqi,qqj-WX_sum_nd))*2
                        WX_matpnbasisP11_E101(qqi,qqj) = &
                            ! 0.5*cunit*(WX_matP11_E101(qqi,qqj-WX_sum_nd)-WX_matBP11_E101(qqi,qqj-WX_sum_nd))
                             +aimag(WX_matP11_E101(qqi,qqj-WX_sum_nd))*2
                end if
                ! (row negative Omega, col positive Omega) 
                if((qqi.gt.WX_sum_nd) .and.(qqj.Le.WX_sum_nd)) then 
                        !-----------------------------E101
                        WX_matpnbasisN11_E101(qqi,qqj) = &
                            !-0.5*cunit*(WX_matN11_E101(qqi-WX_sum_nd,qqj)-WX_matBN11_E101(qqi-WX_sum_nd,qqj))
                             -aimag(WX_matN11_E101(qqi-WX_sum_nd,qqj))*2
                        WX_matpnbasisP11_E101(qqi,qqj) = &
                            !-0.5*cunit*(WX_matP11_E101(qqi-WX_sum_nd,qqj)-WX_matBP11_E101(qqi-WX_sum_nd,qqj))
                             -aimag(WX_matP11_E101(qqi-WX_sum_nd,qqj))*2
                end if
                ! (row negative Omega, col negative Omega) 
                if((qqi.gt.WX_sum_nd) .and.(qqj.gt.WX_sum_nd)) then 
                        !-----------------------------E101
                        WX_matpnbasisN11_E101(qqi,qqj) = &
                            ! 0.5*(WX_matN11_E101(qqi-WX_sum_nd,qqj-WX_sum_nd)+WX_matBN11_E101(qqi-WX_sum_nd,qqj-WX_sum_nd))
                              real(WX_matN11_E101(qqi-WX_sum_nd,qqj-WX_sum_nd))*2
                        WX_matpnbasisP11_E101(qqi,qqj) = &
                            ! 0.5*(WX_matP11_E101(qqi-WX_sum_nd,qqj-WX_sum_nd)+WX_matBP11_E101(qqi-WX_sum_nd,qqj-WX_sum_nd))
                              real(WX_matP11_E101(qqi-WX_sum_nd,qqj-WX_sum_nd))*2
                end if
             end do
          end do
          !----------------------------------QQend             


         !=================================================
         !=================================================
         !=================================================

         !Call WX_cal          
          call WX_cal(nxy,Ph_nb, Ph_db, Ph_readin_select, &
                     Ph_sum_db,Ph_sum_db_colm1, Ph_NumberofX, WX_sum_nd, & 
                     ReX%elem, ImX%elem, ReY%elem, ImY%elem, &
                     WX_matpnbasisN11_E101, WX_matpnbasisP11_E101, &
                     Ph_ReWX, Ph_ImWX, Ph_ReWY, Ph_ImWY, &
                     WgreenX, WgreenY,&
                     WX_iHFB_NB,WX_pnbasis_select,X_select,Y_select)       

             h20re%elem = h20re%elem - Ph_mult*Ph_ReWX*W_omega 
             h20Im%elem = h20Im%elem - Ph_mult*Ph_ImWX*W_omega 
             !
             h02re%elem = h02re%elem - Ph_mult*Ph_ReWY*W_omega 
             h02Im%elem = h02Im%elem + Ph_mult*Ph_ImWY*W_omega    !!!Q: need '-' as (WY*)* conjg.

           ! h20re%elem = h20re%elem + 5.000*Ph_ReWX 
           ! h20Im%elem = h20Im%elem + 5.000*Ph_ImWX 
           ! !
           ! h02re%elem = h02re%elem + 5.000*Ph_ReWY 
           ! h02Im%elem = h02Im%elem - 5.000*Ph_ImWY    !!!Q: need '-' as (WY*)* conjg.
 
 

      !end if    ! file_exists
      !   end do ! which_ph
      !   end do ! which_k
      !   end do ! which_L
      !   end do ! which_iso
      end do    ! nph
        end if  ! Ph_test
      
        !------------------------------------------QQend 
         ! Add the external field (assuming real f20, f02 here)
         h20re%elem = f20%elem + h20re%elem
         h02re%elem = f02%elem + h02re%elem
         !-----------------------------------------QQst       modify21 
         ! implement Ph_WX and Ph_WY to FAM equation for X and Y.
         !---------------------
       ! if(Ph_test) then
       !    !write(*,*) iter, si, Ph_mult ,convergence_epsilon, Ph_test            !xxx
       !     !rite(48,*) iter, si, Ph_mult              !xxx
       !     h20re%elem = h20re%elem + Ph_mult*Ph_ReWX 
       !     h20Im%elem = h20Im%elem + Ph_mult*Ph_ImWX 
       !     !
       !     h02re%elem = h02re%elem + Ph_mult*Ph_ReWY 
       !    !h02Im%elem = h02Im%elem - Ph_mult*Ph_ImWY    !!!Q: '-' instead of '+' as we need take conjg for WY
       !     h02Im%elem = h02Im%elem + Ph_mult*Ph_ImWY    !!!Q: no longer need '-' as no longer need gake conjg.
       ! end if
 
         
         !-----------------------------------QQend  
         if (use_diagonal_blocks) then
            h11re%elem  = h11re%elem(:)*quench ; h11im%elem  = h11im%elem(:)*quench
            h11tre%elem = h11tre%elem(:)*quench; h11tim%elem = h11tim%elem(:)*quench
            h11re%elem  = f11%elem  + h11re%elem
            h11tre%elem = f11t%elem + h11tre%elem
         end if

         end if  ! End computation of the residual interaction
 

         ! Solve the new X and Y
         call complex_multiply(greenXre, greenXim, h20re, h20im, rex, imx)
         call complex_multiply(greenYre, greenYim, h02re, h02im, rey, imy)
 
         ! As well as P and Q
         if (use_diagonal_blocks) then
            call complex_multiply(greenPre, greenPim, h11re,  h11im,  rep, imp)
            call complex_multiply(greenQre, greenQim, h11tre, h11tim, req, imq)
            
            ! Statistical factors (T > 0) or q.p. occupations (A odd)
            rex%elem(:) = mat_1fafb_x%elem(:)*rex%elem(:)
            imx%elem(:) = mat_1fafb_x%elem(:)*imx%elem(:)
            rey%elem(:) = mat_1fafb_y%elem(:)*rey%elem(:)
            imy%elem(:) = mat_1fafb_y%elem(:)*imy%elem(:)

            rep%elem(:) = mat_fbfa_p%elem(:)*rep%elem(:)
            imp%elem(:) = mat_fbfa_p%elem(:)*imp%elem(:)
            req%elem(:) = mat_fbfa_q%elem(:)*req%elem(:)
            imq%elem(:) = mat_fbfa_q%elem(:)*imq%elem(:)
         end if
         
         ! Use Broyden mixing to stabilize convergence
         ! T > 0 / odd-A have twice the number of vectors to handle P and Q
         if (use_diagonal_blocks) then
            call qrpa_broyden(iter, nxy, rex%elem, rey%elem, imx%elem, &
               imy%elem, rep%elem, imp%elem, req%elem, imq%elem)
         else
            call qrpa_broyden(iter, nxy, rex%elem, rey%elem, imx%elem, imy%elem)            
         end if
         
         
         ! Strength function and cross-terms
         ! Real strength
         str(1) = -1.0_dp/pi*(dot_product(f20%elem, imx%elem) + dot_product(f02%elem, imy%elem))
         
         ! Complex strength
         if (present(cstr)) then
            cstr(1) = iu*str(1) - &
               1.0_dp/pi*(dot_product(f20%elem, rex%elem) + dot_product(f02%elem, rey%elem))
         end if
         
         ! Cross-terms
         if (allocated(g20)) then
            do ixterm=1, size(g20)
               str(1+ixterm) = -1.0_dp/pi * &
                  (dot_product(g20(ixterm)%elem, imx%elem) + &
                   dot_product(g02(ixterm)%elem, imy%elem))
      
               if (present(cstr)) then
                  cstr(1+ixterm) = iu*str(1+ixterm) - 1.0_dp/pi * &
                     (dot_product(g20(ixterm)%elem, rex%elem)   + &
                      dot_product(g02(ixterm)%elem, rey%elem))
               end if
            end do
         end if
         
         ! For odd A or T > 0 we also have traces involving diagonal blocks
         if (use_diagonal_blocks) then
            ! Real
            str(1) = str(1) &
               - dot_product(f11%elem,  imp%elem)/pi &
               - dot_product(f11t%elem, imq%elem)/pi
            
            ! Complex
            if (present(cstr)) then
               cstr(1) = iu*str(1) &
                  - dot_product(f20%elem,  rex%elem)/pi &
                  - dot_product(f02%elem,  rey%elem)/pi &
                  - dot_product(f11%elem,  rep%elem)/pi &
                  - dot_product(f11t%elem, req%elem)/pi
            end if
            
            ! Cross terms
            if (allocated(g20)) then
               do ixterm=1, size(g20)
                  str(1+ixterm) = str(1+ixterm) &
                     - dot_product(g11(ixterm)%elem,  imp%elem)/pi &
                     - dot_product(g11t(ixterm)%elem, imq%elem)/pi

                  if (present(cstr)) then
                     cstr(1+ixterm) = iu*str(1+ixterm) &
                        - dot_product(g20(ixterm)%elem,  rex%elem)/pi &
                        - dot_product(g02(ixterm)%elem,  rey%elem)/pi &
                        - dot_product(g11(ixterm)%elem,  rep%elem)/pi &
                        - dot_product(g11t(ixterm)%elem, req%elem)/pi
                  end if
               end do
            end if   
         end if
         
         ! There is an additional prefactor when T > 0
         if (ft_active) then
            ! Note that this has a pole at E=0, and as a result we need to
            ! abandon the computation of the point E = 0 + i\gamma when the
            ! finite-temperature is active. If we are absolutely at (0,0),
            ! the computation should be dropped due to the pole.
            if (abs(cmplx(omega,width,dp)) < 1.0d-6) then
               write(*,'(a)') 'ERROR: cannot compute the strength function at&
                  & finite temperature for E = 0 + i0 MeV.'
               iter_used = -2
               str(:) = 0
               if (present(cstr)) cstr(:) = 0
               stop
            else if (abs(omega) < 1.0d-6) then
               if (iter == 0) write(*,'(a)') 'Warning: cannot compute the real-valued&
                   & strength function at finite temperature for Re(E) = 0 MeV.'
               str(:) = 0
               if (.not.present(cstr)) then
                  iter_used = -2
                  write(*,'(a)') 'Abandoning this point since a complex strength&
                      & function was not requested.'
                  exit
               end if
            end if
         
            ! Actually compute the prefactors
            ! 1) Compute the prefactor which needs Re(E) only
            ! 2) Compute the full complex prefactor
            ! If Re(E)/T is too large, we have problems. In such a case, fall
            ! back to factor = 1/(1-e(-E/T)) ~ 1/(1-1/e(<huge>)) ~ 1
            if (omega/ft_temp > log(huge(1.0_dp))) then
               ft_real_factor  = 1.0_dp
               ft_cmplx_factor = 1.0_dp
               if (iter == 0) then
                  write(*,'(a)') 'Warning: setting ft_factor = 1 to avoid overflows.'
               end if
            else
               ft_real_factor  = (1-exp(-omega/ft_temp))**(-1)
               ft_cmplx_factor = (1-exp(-cmplx(omega,width,dp)/ft_temp))**(-1)
            end if
         
            ! Multiply in the prefactors
            str(:) = str(:)*ft_real_factor
            if (present(cstr)) cstr(:) = cstr(:)*ft_cmplx_factor
         end if
         
         !-------------------------------QQst   modify30
         if(aimag(cstr(1)).gt.0.0) then
                 cstr_mini(:) = cstr(:)
         end if

         if (iter==max_iter) then
                 cstr(:) = cstr_mini(:)
                 exit
         end if
 
         !-------------------------------QQend modify30

         ! Check for convergence
         if (si < convergence_epsilon) then
            iter_used = iter
            exit
         end if
         write(*,*) iter,',', si,',', aimag(cstr),',',omega,',', Ph_test,',', Ph_mult,',',si_mini,',', aimag(cstr_mini)
         write(48,*) iter,',', si,',', aimag(cstr),',',omega,',', Ph_test,',', Ph_mult,',',si_mini,',', aimag(cstr_mini)
      end do
      !rite(48,*) omega,',', aimag(cstr)
  !-----------------------------------QQst modify19    
     !close(33)
     !close(34)
     !close(35)
     !close(36)
     !close(37)
     !close(38)
     !close(39)
     !close(40)
     !close(41)
     !close(42)
     !close(43)
     !close(44)
     !close(45)
     !close(46)
     !close(47)
      close(48)
     !close(49)
      close(50)
    !

   end subroutine pnfam_solve
   
   
   !---------------------------------------------------------------------------
   ! Allocates and initializes the data structures used by pnfam_solve.
   !---------------------------------------------------------------------------
   subroutine init_pnfam_solver(f, vp, up, vn, un, bminus, g)
      
      use constants,    only : IT_NEUTRON, IT_PROTON
      use hfbtho_basis, only : ft_active, hfb_blo_active
      
      implicit none
      
      type(blockmatrix), intent(in) :: f, vp, up, vn, un
      type(external_field), dimension(:), allocatable, intent(in) :: g
      logical, intent(in) :: bminus
      
      integer :: ixterm
      logical :: use_diagonal_blocks
      
      
      ! While the different FAM matrices may have a different non-trivial block
      ! structure, they still share the same number of elements in the non-trivial
      ! blocks:
      nxy = size(f%elem)
      
      
      ! The "diagonal" blocks of the generalized density (P and -Q) are only
      ! active when starting from a finite-temperature or odd-mass HFB solution
      if (ft_active .or. hfb_blo_active) then
         use_diagonal_blocks = .true.
      else
         use_diagonal_blocks = .false.
      end if
      
      ! Cannot handle both T>0 and odd A
      if (hfb_blo_active .and. ft_active) then
         write(*,'(a)') 'ERROR: this code cannot handle T>0 and odd-Z/odd-N simultaneously.'
         stop
      end if
      
      
      ! Allocate needed arrays
      call allocate_blockmatrix(f20, nxy)
      call allocate_blockmatrix(f02, nxy)
      
      ! Special allocation for the cross-terms
      if (allocated(g)) then
         allocate(g20(1:size(g)), g02(1:size(g)))
         
         do ixterm=1, size(g)       
            call allocate_blockmatrix(g20(ixterm), nxy)
            call allocate_blockmatrix(g02(ixterm), nxy)
         end do
         
         if (use_diagonal_blocks) then
            allocate(g11(1:size(g)), g11t(1:size(g)))
            
            do ixterm=1, size(g)       
               call allocate_blockmatrix(g11(ixterm),  nxy)
               call allocate_blockmatrix(g11t(ixterm), nxy)
            end do
         end if
      end if
   
      ! Usual FAM anti-diagonal blocks
      call allocate_blockmatrix(h20re, nxy)
      call allocate_blockmatrix(h20im, nxy)
      call allocate_blockmatrix(h02re, nxy)
      call allocate_blockmatrix(h02im, nxy)
   
      call allocate_blockmatrix(rex, nxy)
      call allocate_blockmatrix(imx, nxy)
      call allocate_blockmatrix(rey, nxy)
      call allocate_blockmatrix(imy, nxy)
   
      call allocate_blockmatrix(greenXre, nxy)
      call allocate_blockmatrix(greenXim, nxy)
      call allocate_blockmatrix(greenYre, nxy)
      call allocate_blockmatrix(greenYim, nxy)

      call allocate_blockmatrix(reh_np, nxy)
      call allocate_blockmatrix(imh_np, nxy)
      call allocate_blockmatrix(reh_pn, nxy)
      call allocate_blockmatrix(imh_pn, nxy)
      call allocate_blockmatrix(rerho_pn, nxy)
      call allocate_blockmatrix(imrho_pn, nxy)
      call allocate_blockmatrix(rerho_np, nxy)
      call allocate_blockmatrix(imrho_np, nxy)
   
      call allocate_blockmatrix(redp, nxy)
      call allocate_blockmatrix(imdp, nxy)
      call allocate_blockmatrix(redm, nxy)
      call allocate_blockmatrix(imdm, nxy)
      call allocate_blockmatrix(rekp, nxy)
      call allocate_blockmatrix(imkp, nxy)
      call allocate_blockmatrix(rekm, nxy)
      call allocate_blockmatrix(imkm, nxy)

      ! Diagonal blocks for T > 0 and odd-A nuclei
      if (use_diagonal_blocks) then
         call allocate_blockmatrix(f11,         nxy)    ! Diagonal blocks of F
         call allocate_blockmatrix(f11t,        nxy)
         call allocate_blockmatrix(rep,         nxy)    ! Diagonal blocks of \delta R
         call allocate_blockmatrix(imp,         nxy)
         call allocate_blockmatrix(req,         nxy)
         call allocate_blockmatrix(imq,         nxy)
         call allocate_blockmatrix(h11re,       nxy)    ! Diagonal blocks of \delta H
         call allocate_blockmatrix(h11im,       nxy)
         call allocate_blockmatrix(h11tre,      nxy)
         call allocate_blockmatrix(h11tim,      nxy)
         call allocate_blockmatrix(greenPre,    nxy)    ! Diagonal blocks of Green's function
         call allocate_blockmatrix(greenPim,    nxy)
         call allocate_blockmatrix(greenQre,    nxy)
         call allocate_blockmatrix(greenQim,    nxy)
         call allocate_blockmatrix(mat_fbfa_p,  nxy)    ! Statistical factors
         call allocate_blockmatrix(mat_fbfa_q,  nxy)
         call allocate_blockmatrix(mat_1fafb_x, nxy)
         call allocate_blockmatrix(mat_1fafb_y, nxy)
      end if
   
      ! Allocate the coordinate-space densities
      call allocate_density(rhox)
   
      ! Transform the external operator to the quasiparticle basis
      if (bminus) then
         call triprod('t', up, 'n', f, 'n', vn,  1d0, 0d0, f20)      !  up' f- vn
         call triprod('t', vp, 'n', f, 'n', un, -1d0, 0d0, f02)      ! -vp' f- un
         
         if (use_diagonal_blocks) then
            call triprod('t', up, 'n', f, 'n', un, 1d0, 0d0, f11)    !  up' f- un
            call triprod('t', vp, 'n', f, 'n', vn,-1d0, 0d0, f11t)   ! -vp' f- vn
         end if
      else
         call triprod('t', un, 'n', f, 'n', vp,  1d0, 0d0, f20)      !  un' f+ vp
         call triprod('t', vn, 'n', f, 'n', up, -1d0, 0d0, f02)      ! -vn' f+ un
         
         if (use_diagonal_blocks) then
            call triprod('t', un, 'n', f, 'n', up, 1d0, 0d0, f11)    !  un' f+ up
            call triprod('t', vn, 'n', f, 'n', vp,-1d0, 0d0, f11t)   ! -vn' f+ vp
         end if
      end if
      
      ! Special transformation for the cross-terms
      if (allocated(g)) then
         do ixterm=1, size(g)
            if (bminus) then
               call triprod('t', up, 'n', g(ixterm)%mat, 'n', vn,  1d0, 0d0, g20(ixterm))
               call triprod('t', vp, 'n', g(ixterm)%mat, 'n', un, -1d0, 0d0, g02(ixterm))
               
               if (use_diagonal_blocks) then
                  call triprod('t', up, 'n', g(ixterm)%mat, 'n', un,  1d0, 0d0, g11(ixterm))
                  call triprod('t', vp, 'n', g(ixterm)%mat, 'n', vn, -1d0, 0d0, g11t(ixterm))
               end if
            else
               call triprod('t', un, 'n', g(ixterm)%mat, 'n', vp,  1d0, 0d0, g20(ixterm))
               call triprod('t', vn, 'n', g(ixterm)%mat, 'n', up, -1d0, 0d0, g02(ixterm))
               
               if (use_diagonal_blocks) then
                  call triprod('t', un, 'n', g(ixterm)%mat, 'n', up,  1d0, 0d0, g11(ixterm))
                  call triprod('t', vn, 'n', g(ixterm)%mat, 'n', vp, -1d0, 0d0, g11t(ixterm))
               end if
            end if
         end do
      end if
      
      ! Initialize the block structure of various matrices which
      ! are not constructed with the routine triprod
      
      ! Anti-diagonal blocks
      call copy_block_structure(f20, rex)  ! X has the same structure as F20
      call copy_block_structure(f20, imx)
      call copy_block_structure(f02, rey)  ! Y has the same structure as F02
      call copy_block_structure(f02, imy)
   
      call copy_block_structure(f, reh_pn)  ! h_pn has the same structure as f
      call copy_block_structure(f, imh_pn)
      call copy_block_structure(rex, redp)  ! Delta(+) has the same structure as X
      call copy_block_structure(imx, imdp)
      
      call copy_block_structure(f, reh_np)  ! h_np has the same structure as f
      call copy_block_structure(f, imh_np)
      call copy_block_structure(rey, redm)  ! Delta(-) has the same structure as Y
      call copy_block_structure(imy, imdm)
      
      call copy_block_structure(rex, greenXre)  ! -1/(Ep+En-omega) has the same structure as X
      call copy_block_structure(imx, greenXim)
      call copy_block_structure(rey, greenYre)  ! -1/(Ep+En+omega*) has the same structure as Y
      call copy_block_structure(imy, greenYim)
      
      call copy_block_structure(f20, h20re)  ! delta-H20 has the same structure as F20
      call copy_block_structure(f20, h20im)
      call copy_block_structure(f02, h02re)  ! delta-H02 has the same structure as F02
      call copy_block_structure(f02, h02im)
      
      ! Diagonal blocks
      if (use_diagonal_blocks) then
         call copy_block_structure(f11,  ReP)
         call copy_block_structure(f11,  ImP)
         call copy_block_structure(f11,  greenPre)
         call copy_block_structure(f11,  greenPim)
         call copy_block_structure(f11,  h11re)
         call copy_block_structure(f11,  h11im)
         call copy_block_structure(f11t, ReQ)
         call copy_block_structure(f11t, ImQ)
         call copy_block_structure(f11t, greenQre)
         call copy_block_structure(f11t, greenQim)
         call copy_block_structure(f11t, h11tre)
         call copy_block_structure(f11t, h11tim)
         call copy_block_structure(f11,  mat_fbfa_p)
         call copy_block_structure(f11t, mat_fbfa_q)
         call copy_block_structure(f20,  mat_1fafb_x)
         call copy_block_structure(f02,  mat_1fafb_y)
         
         ! The (1-fa-fb) and (fb-fa) matrices are different
         ! depending on whether A is odd or T > 0.
         mat_1fafb_x%elem(:) = 0
         mat_1fafb_y%elem(:) = 0
         mat_fbfa_p%elem(:)  = 0
         mat_fbfa_q%elem(:)  = 0
         
         if (ft_active) then
            if (bminus .eqv. .true.) then
               call ft_1fafb_matrix(ta=IT_PROTON, tb=IT_NEUTRON, mat=mat_1fafb_x)
               call ft_1fafb_matrix(ta=IT_PROTON, tb=IT_NEUTRON, mat=mat_1fafb_y)
               call ft_fbfa_matrix (ta=IT_PROTON, tb=IT_NEUTRON, mat=mat_fbfa_p)
               call ft_fbfa_matrix (ta=IT_PROTON, tb=IT_NEUTRON, mat=mat_fbfa_q)
            else
               call ft_1fafb_matrix(ta=IT_NEUTRON, tb=IT_PROTON, mat=mat_1fafb_x)
               call ft_1fafb_matrix(ta=IT_NEUTRON, tb=IT_PROTON, mat=mat_1fafb_y)
               call ft_fbfa_matrix (ta=IT_NEUTRON, tb=IT_PROTON, mat=mat_fbfa_p)
               call ft_fbfa_matrix (ta=IT_NEUTRON, tb=IT_PROTON, mat=mat_fbfa_q)
            end if
         else if (hfb_blo_active) then
            call odd_efa_fafb_matrix(class='XY', bminus=bminus, mat=mat_1fafb_x)   
            call odd_efa_fafb_matrix(class='XY', bminus=bminus, mat=mat_1fafb_y)   
            call odd_efa_fafb_matrix(class='PQ', bminus=bminus, mat=mat_fbfa_p)   
            call odd_efa_fafb_matrix(class='PQ', bminus=bminus, mat=mat_fbfa_q)
         end if
      end if
      
   end subroutine init_pnfam_solver
   
   
   subroutine deallocate_pnfam_solver(g)
      use hfbtho_basis, only : ft_active, hfb_blo_active
      implicit none
      type(external_field), dimension(:), allocatable :: g
      
      integer :: ixterm
      logical :: use_diagonal_blocks
      
      if (ft_active .or. hfb_blo_active) then
         use_diagonal_blocks = .true.
      else
         use_diagonal_blocks = .false.
      end if
      
      call deallocate_blockmatrix(f20)
      call deallocate_blockmatrix(f02)
      
      if (use_diagonal_blocks) then
         call deallocate_blockmatrix(f11)
         call deallocate_blockmatrix(f11t)         
      end if
      
      ! Special allocation for the cross-terms
      if (allocated(g)) then
         
         do ixterm=1, size(g)       
            call deallocate_blockmatrix(g20(ixterm))
            call deallocate_blockmatrix(g02(ixterm))
         end do
         
         if (use_diagonal_blocks) then
            do ixterm=1, size(g)       
               call deallocate_blockmatrix(g11(ixterm))
               call deallocate_blockmatrix(g11t(ixterm))
            end do
            
            deallocate(g11, g11t)
         end if
         
         deallocate(g20, g02)
      end if

      call deallocate_blockmatrix(h20re)
      call deallocate_blockmatrix(h20im)
      call deallocate_blockmatrix(h02re)
      call deallocate_blockmatrix(h02im)
   
      call deallocate_blockmatrix(rex)
      call deallocate_blockmatrix(imx)
      call deallocate_blockmatrix(rey)
      call deallocate_blockmatrix(imy)
   
      call deallocate_blockmatrix(greenXre)
      call deallocate_blockmatrix(greenXim)
      call deallocate_blockmatrix(greenYre)
      call deallocate_blockmatrix(greenYim)
   
      call deallocate_blockmatrix(reh_np)
      call deallocate_blockmatrix(imh_np)
      call deallocate_blockmatrix(reh_pn)
      call deallocate_blockmatrix(imh_pn)
      call deallocate_blockmatrix(rerho_pn)
      call deallocate_blockmatrix(imrho_pn)
      call deallocate_blockmatrix(rerho_np)
      call deallocate_blockmatrix(imrho_np)
   
      call deallocate_blockmatrix(redp)
      call deallocate_blockmatrix(imdp)
      call deallocate_blockmatrix(redm)
      call deallocate_blockmatrix(imdm)
      call deallocate_blockmatrix(rekp)
      call deallocate_blockmatrix(imkp)
      call deallocate_blockmatrix(rekm)
      call deallocate_blockmatrix(imkm)
      
      if (use_diagonal_blocks) then
         call deallocate_blockmatrix(f11)    ! Diagonal blocks of F
         call deallocate_blockmatrix(f11t)
         call deallocate_blockmatrix(rep)    ! Diagonal blocks of \delta R
         call deallocate_blockmatrix(imp)
         call deallocate_blockmatrix(req)
         call deallocate_blockmatrix(imq)
         call deallocate_blockmatrix(h11re)    ! Diagonal blocks of \delta H
         call deallocate_blockmatrix(h11im)
         call deallocate_blockmatrix(h11tre)
         call deallocate_blockmatrix(h11tim)
         call deallocate_blockmatrix(greenPre)    ! Diagonal blocks of Green's function
         call deallocate_blockmatrix(greenPim)
         call deallocate_blockmatrix(greenQre)
         call deallocate_blockmatrix(greenQim)
         call deallocate_blockmatrix(mat_fbfa_p)   ! Statistical factors
         call deallocate_blockmatrix(mat_fbfa_q)
         call deallocate_blockmatrix(mat_1fafb_x)
         call deallocate_blockmatrix(mat_1fafb_y)
      end if
      
   end subroutine deallocate_pnfam_solver
   
   
   !---------------------------------------------------------------------------
   ! A wrapper for the more general modified Broyden mixing routine
   !---------------------------------------------------------------------------
   subroutine qrpa_broyden(niter, nuv, x1, x2, x3, x4, x5, x6, x7, x8)
      use broyden_mixer, only : qrpa_alphamix, qrpa_broyden_method, si, &
         broyden_history_size, qrpa_bbroyden
      implicit none
      integer, intent(in) :: nuv
      real(dp), intent(inout)           :: x1(nuv), x2(nuv), x3(nuv), x4(nuv)
      real(dp), intent(inout), optional :: x5(nuv), x6(nuv), x7(nuv), x8(nuv)
      integer, intent(in)  :: niter
      integer :: i, nvec
      real(dp), allocatable, save :: qrpa_broin(:), qrpa_broout(:)
      
      ! Calculate the size of the broyden vectors.
      ! 5, 6, 7, 8 are present at finite temperature
      i = 4*nuv
      if (present(x5)) i = i+nuv
      if (present(x6)) i = i+nuv
      if (present(x7)) i = i+nuv
      if (present(x8)) i = i+nuv
      nvec = i/nuv

      ! On the first iteration
      if (niter==0) then
         if (allocated(qrpa_broin)) deallocate(qrpa_broin,qrpa_broout)
         allocate(qrpa_broin(i),qrpa_broout(i))
         qrpa_broin = 0
      end if

      ! The Broyden mixing routine operates on a real array; Therefore we need
      ! to unpack the complex values to a real array, call the Broyden routine,
      ! and then repack the real and imaginary parts to the complex arrays
      qrpa_broout(1:nuv) = x1
      qrpa_broout(nuv+1:2*nuv) = x2
      qrpa_broout(2*nuv+1:3*nuv) = x3
      qrpa_broout(3*nuv+1:4*nuv) = x4
      
      ! nvec should be either 4 (T=0) or 8 (T>0)
      if (nvec /= 4 .and. nvec /= 8) then
         write(*,'(a)') 'Error: nvec /= 4 or 8 in qrpa_broyden'
         stop
      end if
      
      ! Set up the Broyden vectors
      if (nvec == 8) then
         qrpa_broout(4*nuv+1:5*nuv) = x5
         qrpa_broout(5*nuv+1:6*nuv) = x6
         qrpa_broout(6*nuv+1:7*nuv) = x7
         qrpa_broout(7*nuv+1:8*nuv) = x8
      end if
      
      ! Broyden
      call qrpa_broyden_method(i,qrpa_broout,qrpa_broin,qrpa_alphamix,si,niter,broyden_history_size,qrpa_bbroyden)
      if (niter==0) si = 1

      x1 = qrpa_broin(1:nuv)
      x2 = qrpa_broin(nuv+1:2*nuv)
      x3 = qrpa_broin(2*nuv+1:3*nuv)
      x4 = qrpa_broin(3*nuv+1:4*nuv)

      ! Set up the Broyden vectors
      if (nvec == 8) then
         x5 = qrpa_broin(4*nuv+1:5*nuv)
         x6 = qrpa_broin(5*nuv+1:6*nuv)
         x7 = qrpa_broin(6*nuv+1:7*nuv)
         x8 = qrpa_broin(7*nuv+1:8*nuv)
      end if

   end subroutine qrpa_broyden


   !---------------------------------------------------------------------------
   ! The "Green functions" in the FAM equations.  The minus sign has been
   ! absorbed here.  This version has been extended for general K and parity.
   ! Note omega is complex, RE(omega)=omega, IM(omega)=half_width.
   !---------------------------------------------------------------------------
   subroutine greenx(E1, E2, omega, greenRe, greenIm)
      use hfbtho_basis
      implicit none
      real(dp), intent(in) :: E1(:), E2(:)
      complex(dp), intent(in) :: omega
      type(blockmatrix), intent(inout) :: greenRe, greenIm

      integer :: ir, ic, ipt, i1, i2, i1a, i1b, i2a, i2b
      !
      !------------------------------QQst
      integer :: qqi
    ! open(45,file='QQ_45_REqpN')
    ! open(46,file='QQ_46_REqpP')
      !--------
    ! do qqi = 1, sizeof(E1)
    !    write(45,*) qqi, E1(qqi)
    ! end do
      !-----------------
    ! do qqi = 1, sizeof(E2)
    !    write(46,*) qqi, E2(qqi)
    ! end do
      !
      !------------------------------QQend
 
     ipt = 1
      do ir=1,nb   ! block row index
         ic = greenRe%ir2c(ir)   ! block column index
         if (ic == 0) cycle
         i1a = isstart(ir) ; i1b = i1a + db(ir) - 1
         i2a = isstart(ic) ; i2b = i2a + db(ic) - 1
         do i2 = i2a,i2b   ! true column index
         do i1 = i1a,i1b   ! true row index
            greenRe%elem(ipt) = real(-1/(E1(i1) + E2(i2) + omega), dp)
            greenIm%elem(ipt) = aimag(-1/(E1(i1) + E2(i2) + omega))
            ipt = ipt + 1
         end do
         end do
      end do

      !
     !close(45)
     !close(46)
   end subroutine greenx
   
      
   ! The moment of inertia of the HFB solution in the cranking approximation
   ! Call this separately for protons and neutrons, add the results up
   function moment_of_inertia(E, U, V) result(moi)
      use external_field_type
      use extfield
      use hfbtho_basis
      implicit none
      real(dp), intent(in) :: E(:)
      type(blockmatrix), intent(in) :: U, V  ! qp energies and amplitudes
      real(dp) :: moi
      
      type(blockmatrix) :: M1, M2
      type(external_field) :: Jplus, Jminus
      integer :: ipt, ibx1, ibx2, i1, i2, nd1, nd2, ix1, ix2
      
      ! Finite temperature / odd A
      if (ft_active) then
         call writelog('Warning: moment_of_inertia has not been updated for finite-temperature.')
      end if
      if (hfb_blo_active) then
         call writelog('Warning: moment_of_inertia has not been updated for odd-mass nuclei.')
      end if
      call writelog('Warning: moment_of_inertia has not been checked for using the PWI.')
      
      ! Because the blockmatrix type requires there to be only one nonzero block
      ! in each row and column, we need to separate the matrix into two:
      ! M1 = U^dagger J_+ V^* - V^dagger J_- U^*
      ! M2 = U^dagger J_- V^* - V^dagger J_+ U^*
      
      ! Construct raising and lowering operators of the total angular momentum
      Jplus%k = 1 ; Jplus%parity_even = .true. ; Jplus%rank = 1
      call init_fam_mapping(Jplus)
      ipt = 1
      do ibx1 = 1,nb
         ibx2 = Jplus%mat%ir2c(ibx1)
         if (ibx2 == 0) cycle
         nd1 = db(ibx1) ; nd2 = db(ibx2)
         if (nd1 == 0 .or. nd2 == 0) cycle
         
         do i2=1, nd2
            ! Calculate the true index of the state |2>
            ix2 = i2 + isstart(ibx2) - 1

            do i1=1, nd1
               ! Calculate the true index of the state <1|
               ix1 = i1 + isstart(ibx1) - 1
               
               !if (nr(ix1)==nr(ix2).and.nz(ix1)==nz(ix2)) then
                  if (nl(ix1)==nl(ix2) .and. ns(ix1)==ns(ix2)+2) then
                     Jplus%mat%elem(ipt) = dot_product(wf(:,ix1), wf(:,ix2))  ! ok
                  end if
                  if (nl(ix1)==nl(ix2)+1 .and. ns(ix1)==ns(ix2)) then
                     Jplus%mat%elem(ipt) = dot_product(-wf(:,ix1)/y(:), wfdz(:,ix2)) &
                        + dot_product(wf(:,ix1)*z(:), wfdr(:,ix2)) &
                        - nl(ix2)*dot_product(wf(:,ix1)*z(:)*y(:), wf(:,ix2))
                  end if
               !end if
               
               ipt = ipt + 1
            end do
         end do
      end do
      
      Jminus%k = -1 ; Jminus%parity_even = .true. ; Jminus%rank = 1
      call init_fam_mapping(Jminus)
      ipt = 1
      do ibx1 = 1,nb
         ibx2 = Jminus%mat%ir2c(ibx1)
         if (ibx2 == 0) cycle
         nd1 = db(ibx1) ; nd2 = db(ibx2)
         if (nd1 == 0 .or. nd2 == 0) cycle
         
         do i2=1, nd2
            ! Calculate the true index of the state |2>
            ix2 = i2 + isstart(ibx2) - 1

            do i1=1, nd1
               ! Calculate the true index of the state <1|
               ix1 = i1 + isstart(ibx1) - 1
               
               !if (nr(ix1)==nr(ix2).and.nz(ix1)==nz(ix2)) then
                  if (nl(ix1)==nl(ix2) .and. ns(ix1)==ns(ix2)-2) then
                     Jminus%mat%elem(ipt) = dot_product(wf(:,ix1)/y(:), wfdz(:,ix2)) &
                        - dot_product(z(:)*wf(:,ix1), wfdr(:,ix2)) &
                        - nl(ix2)*dot_product(wf(:,ix1)*z(:)*y(:), wf(:,ix2))
                  end if
                  if (nl(ix1)==nl(ix2)-1 .and. ns(ix1)==ns(ix2)) then
                     Jminus%mat%elem(ipt) = dot_product(wf(:,ix1), wf(:,ix2))
                  end if
               !end if

               ipt = ipt + 1
            end do
         end do
      end do
      
      ! Prep the matrices M1 and M2
      call allocate_blockmatrix(M1,size(Jplus%mat%elem))
      call allocate_blockmatrix(M2,size(Jplus%mat%elem))
      
      ! Transform them into the qp basis
      call triprod('T', U, 'N', Jplus%mat, 'N', V, 1d0, 0d0, M1)
      call triprod('T', V, 'N', Jminus%mat, 'N', U, -1d0, 1d0, M1)
      call triprod('T', U, 'N', Jminus%mat, 'N', V, 1d0, 0d0, M2)
      call triprod('T', V, 'N', Jplus%mat, 'N', U, -1d0, 1d0, M2)
      
      ! Sum over all qp pairs, block by block
      moi = 0
      ipt = 1
      do ibx1 = 1,nb
         ibx2 = M1%ir2c(ibx1)
         if (ibx2 == 0) cycle
         nd1 = db(ibx1) ; nd2 = db(ibx2)
         if (nd1 == 0 .or. nd2 == 0) cycle
         
         do i2=1, nd2
            ix2 = i2 + isstart(ibx2) - 1
            do i1=1, nd1
               ix1 = i1 + isstart(ibx1) - 1
               moi = moi + M1%elem(ipt)**2/(E(ix1) + E(ix2))
               ipt = ipt + 1
            end do
         end do
      end do
      ipt = 1
      do ibx1 = 1,nb
         ibx2 = M2%ir2c(ibx1)
         if (ibx2 == 0) cycle
         nd1 = db(ibx1) ; nd2 = db(ibx2)
         if (nd1 == 0 .or. nd2 == 0) cycle
         
         do i2=1, nd2
            ix2 = i2 + isstart(ibx2) - 1
            do i1=1, nd1
               ix1 = i1 + isstart(ibx1) - 1
               moi = moi + M2%elem(ipt)**2/(E(ix1) + E(ix2))
               ipt = ipt + 1
            end do
         end do
      end do
      
      moi = moi/4
      
   end function moment_of_inertia

   
   !----------------------------------------------------------------------------
   ! Compute the finite-temperature statistical factor matrices (1 - f_a - f_b)
   !----------------------------------------------------------------------------
   subroutine ft_1fafb_matrix(ta, tb, mat)
      use constants, only : IT_NEUTRON, IT_PROTON
      use hfbtho_basis, only : pwi_start_n, pwi_start_p, pwi_dim_n, pwi_dim_p, &
         pwi_qp_n, pwi_qp_p, ft_active, ft_fp, ft_fn, isstart
      
      implicit none
      
      integer, intent(in) :: ta, tb
      type(blockmatrix), intent(inout) :: mat
      
      integer :: im, ibr, ibc, ia, ib, ik, ka, kb
      logical :: active_a, active_b
      
      integer,  pointer :: kstarta(:), kstartb(:), kdima(:), kdimb(:), kqpa(:), kqpb(:)
      real(dp), pointer :: fa(:), fb(:)
      
      
      mat%elem(:) = 0
      
      ! Check that finite temperature should be on; if not, exit with the defaults
      if (.not.ft_active) then
         write(*,'(a)') 'WARNING: finite temperature is not active, but ft_1fafb_matrix was called.'
         return
      end if
      
      ! Check isospin integers
      if  ((ta /= IT_NEUTRON .and. ta /= IT_PROTON)   &
      .or. (tb /= IT_NEUTRON .and. tb /= IT_PROTON)) then
         write(*,'(A)') 'ERROR: bad isospin integers in ft_1fafb_matrix.'
         return
      end if
      
      ! Set up pointers
      nullify(kstarta, kdima, kqpa, fa)
      if (ta == IT_NEUTRON) then
         kstarta => pwi_start_n;  kdima => pwi_dim_n;  kqpa => pwi_qp_n;  fa => ft_fn
      else if (ta == IT_PROTON) then
         kstarta => pwi_start_p;  kdima => pwi_dim_p;  kqpa => pwi_qp_p;  fa => ft_fp
      end if
      
      nullify(kstartb, kdimb, kqpb, fb)
      if (tb == IT_NEUTRON) then
         kstartb => pwi_start_n;  kdimb => pwi_dim_n;  kqpb => pwi_qp_n;  fb => ft_fn
      else if (tb == IT_PROTON) then
         kstartb => pwi_start_p;  kdimb => pwi_dim_p;  kqpb => pwi_qp_p;  fb => ft_fp
      end if

      im = 0; ka = 0; kb = 0
      
      do ibr=1, nb
         
         ibc = mat%ir2c(ibr)
         if (ibc == 0) cycle
         
         do ib=isstart(ibc),isstart(ibc)+db(ibc)-1
            ! Is the q.p. with index ib active?
            active_b = .false.
            do ik=kstartb(ibc),kstartb(ibc)+kdimb(ibc)-1
               if (kqpb(ik) == ib) then
                  active_b = .true.
                  kb = ik
                  exit
               end if
            end do
            
            do ia=isstart(ibr),isstart(ibr)+db(ibr)-1
               ! Is the q.p. with index ia active?
               active_a = .false.
               do ik=kstarta(ibr),kstarta(ibr)+kdima(ibr)-1
                  if (kqpa(ik) == ia) then
                     active_a = .true.
                     ka = ik
                     exit
                  end if
               end do

               im = im + 1
               if (active_a .and. active_b) mat%elem(im) = 1 - fa(ka) - fb(kb)
            end do
         end do
      end do

   end subroutine ft_1fafb_matrix
   
   !----------------------------------------------------------------------------
   ! Compute the finite-temperature statistical factor matrices (f_b - f_a)
   !----------------------------------------------------------------------------
   subroutine ft_fbfa_matrix(ta, tb, mat)
      use constants, only : IT_NEUTRON, IT_PROTON
      use hfbtho_basis, only : pwi_start_n, pwi_start_p, pwi_dim_n, pwi_dim_p, &
         pwi_qp_n, pwi_qp_p, ft_active, ft_fp, ft_fn, isstart
      
      implicit none
      
      integer, intent(in) :: ta, tb
      type(blockmatrix), intent(inout) :: mat
      
      integer :: im, ibr, ibc, ia, ib, ik, ka, kb
      logical :: active_a, active_b
      
      integer,  pointer :: kstarta(:), kstartb(:), kdima(:), kdimb(:), kqpa(:), kqpb(:)
      real(dp), pointer :: fa(:), fb(:)
      
      
      mat%elem(:) = 0
      
      ! Check that finite temperature should be on; if not, exit with the defaults
      if (.not.ft_active) then
         write(*,'(a)') 'WARNING: finite temperature is not active, but ft_1fafb_matrix was called.'
         return
      end if
      
      ! Check isospin integers
      if  ((ta /= IT_NEUTRON .and. ta /= IT_PROTON)   &
      .or. (tb /= IT_NEUTRON .and. tb /= IT_PROTON)) then
         write(*,'(A)') 'ERROR: bad isospin integers in ft_1fafb_matrix.'
         return
      end if
      
      ! Set up pointers
      nullify(kstarta, kdima, kqpa, fa)
      if (ta == IT_NEUTRON) then
         kstarta => pwi_start_n;  kdima => pwi_dim_n;  kqpa => pwi_qp_n;  fa => ft_fn
      else if (ta == IT_PROTON) then
         kstarta => pwi_start_p;  kdima => pwi_dim_p;  kqpa => pwi_qp_p;  fa => ft_fp
      end if
      
      nullify(kstartb, kdimb, kqpb, fb)
      if (tb == IT_NEUTRON) then
         kstartb => pwi_start_n;  kdimb => pwi_dim_n;  kqpb => pwi_qp_n;  fb => ft_fn
      else if (tb == IT_PROTON) then
         kstartb => pwi_start_p;  kdimb => pwi_dim_p;  kqpb => pwi_qp_p;  fb => ft_fp
      end if

      im = 0; ka = 0; kb = 0
      
      do ibr=1, nb
         
         ibc = mat%ir2c(ibr)
         if (ibc == 0) cycle
         
         do ib=isstart(ibc),isstart(ibc)+db(ibc)-1
            ! Is the q.p. with index ib active?
            active_b = .false.
            do ik=kstartb(ibc),kstartb(ibc)+kdimb(ibc)-1
               if (kqpb(ik) == ib) then
                  active_b = .true.
                  kb = ik
                  exit
               end if
            end do
            
            do ia=isstart(ibr),isstart(ibr)+db(ibr)-1
               ! Is the q.p. with index ia active?
               active_a = .false.
               do ik=kstarta(ibr),kstarta(ibr)+kdima(ibr)-1
                  if (kqpa(ik) == ia) then
                     active_a = .true.
                     ka = ik
                     exit
                  end if
               end do

               im = im + 1
               if (active_a .and. active_b) mat%elem(im) = fb(kb) - fa(ka)
            end do
         end do
      end do
   
   end subroutine ft_fbfa_matrix
   
   !----------------------------------------------------------------------------
   ! In the equal-filling approximation there are finite-temperature-style
   ! q.p. occupation factors to multiply X, Y, P, and Q
   !----------------------------------------------------------------------------
   subroutine odd_efa_fafb_matrix(class, bminus, mat)
      
      use constants,        only : translate_uppercase, IT_NEUTRON, IT_PROTON
      use blockmatrix_type, only : blockmatrix, nb, db
      use hfbtho_basis,     only : isstart, hfb_blo_active, hfb_blo_qp, qp_t_reverse
      
      implicit none
      
      character(len=2),  intent(in)    :: class
      logical,           intent(in)    :: bminus
      type(blockmatrix), intent(inout) :: mat
         
      character(len=2) :: class_
      
      logical  :: l1, l2
      integer  :: ibr, ibc, ic, ir, iqpc, iqpr, im, blo_n, blo_p, blo_nr, blo_pr
      real(dp) :: fr, fc
      
      
      ! Check class and rowcol
      class_  = class
      call translate_uppercase(class_)
      
      if (class_ /= 'PQ' .and. class_ /= 'XY') then
         write(*,'("ERROR: unrecognized class ",a," in odd_efa_fafb_matrix")') class_
         stop
      end if

      ! Odd q.p. indices
      blo_nr = 0;  blo_pr = 0
      
      blo_n = hfb_blo_qp(IT_NEUTRON)
      blo_p = hfb_blo_qp(IT_PROTON)
      
      if (blo_n /= 0) blo_nr = qp_t_reverse(blo_n)
      if (blo_p /= 0) blo_pr = qp_t_reverse(blo_p)
      
      ! Multiply in the factors
      mat%elem = 0
      
      do ibr=1, nb
         ibc = mat%ir2c(ibr)
         if (ibc == 0) cycle
         
         do ic=1, db(ibc)
            iqpc = isstart(ibc) + ic - 1
            
            ! Beta-minus & neutron / Beta-plus & proton
            l1 = ((hfb_blo_active .and. bminus .and. blo_n /= 0) .and. &
                  (iqpc == blo_n .or. iqpc == blo_nr))
            l2 = ((hfb_blo_active .and. (.not. bminus) .and. blo_p /= 0) .and. &
                  (iqpc == blo_p .or. iqpc == blo_pr))
            
            if (l1 .or. l2) then
               fc = 0.5_dp
            else
               fc = 0
            end if

            do ir=1, db(ibr)
               iqpr = isstart(ibr) + ir - 1
               
               ! Beta-minus & proton / Beta-plus & neutron
               l1 = ((hfb_blo_active .and. bminus .and. blo_p /= 0) .and. &
                     (iqpr == blo_p .or. iqpr == blo_pr))
               l2 = ((hfb_blo_active .and. (.not. bminus) .and. blo_n /= 0) .and. &
                     (iqpr == blo_n .or. iqpr == blo_nr))
            
               if (l1 .or. l2) then
                  fr = 0.5_dp
               else
                  fr = 0
               end if
               
               im = mat%ir2m(ibr) + (ic-1)*db(ibr) + ir - 1
               
               if (class_ == 'PQ') then
                  mat%elem(im) = fc - fr
               else if (class_ == 'XY') then
                  mat%elem(im) = 1.0_dp - fr - fc
               end if
            end do
         end do
      end do

   end subroutine odd_efa_fafb_matrix

end module pnfam_solver
