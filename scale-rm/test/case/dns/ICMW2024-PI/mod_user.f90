!-------------------------------------------------------------------------------
!> module User
!!
!! @par Description
!!          
!!          ICMW 2021, PI chamber case, 
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2021-6-5 (S. Shima)   [new]
!! @li      2021-6-8 (S. Shima)   [add] set the boundary condition
!!
!<
!-------------------------------------------------------------------------------
module mod_user
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_tracer
  use scale_time
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: USER_config
  public :: USER_setup
  public :: USER_resume0
  public :: USER_resume
  public :: USER_step

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical,  private, save :: USER_do  = .true.
  real(RP), private, save :: Cm_max = 784.0_RP
  real(RP), private, save :: TEMP_BOTTOM = 277.15_RP ! bottom wall temperature [K]
  real(RP), private, save :: TEMP_TOP    = 257.15_RP ! top wall temperature [K]
  real(RP), private, save :: TEMP_SIDE   = 261.15_RP ! side wall temperature [K]
  real(RP), private, save :: RHICE_SIDE  = 0.3_RP    ! side wall relative humidity over ice []
 
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Config
  subroutine USER_config

    return
  end subroutine USER_config

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine USER_setup
    use scale_process, only: &
       PRC_MPIstop
    implicit none
    integer :: ierr

    namelist / PARAM_USER / &
       USER_do,        &
       Cm_max, &
       TEMP_BOTTOM, &
       TEMP_TOP, &
       TEMP_SIDE, &
       RHICE_SIDE

    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[USER]/Categ[MAIN]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_USER,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_USER. Check!'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_USER)

    if ( USER_do ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Enable User step'

    endif

    return
  end subroutine USER_setup

  !-----------------------------------------------------------------------------
  !> Resuming operation, before calculating tendency
  subroutine USER_resume0
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine USER_resume0

  !-----------------------------------------------------------------------------
  !> Resuming operation
  subroutine USER_resume
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine USER_resume

  !-----------------------------------------------------------------------------
  !> Step
  !! This subroutine is called every TIME_DT
  subroutine USER_step
    use scale_history, only: &
       HIST_in
    use scale_process, only: &
       PRC_MPIstop, mype => PRC_myrank
    use scale_rm_process, only: &
       PRC_HAS_W,   &
       PRC_HAS_E,   &
       PRC_HAS_S,   &
       PRC_HAS_N
    use scale_const, only: &
       P00   => CONST_PRE00
    use mod_atmos_vars, only: &
       DENS, &
       MOMX, &
       MOMY, &
       MOMZ, &
       RHOT, &
       QTRC, &
       DENS_tp, &
       MOMZ_tp, &
       MOMX_tp, &
       MOMY_tp, &
       RHOT_tp, &
       RHOQ_tp
    use scale_grid, only: &
       CX => GRID_CX, &
       CY => GRID_CY, &
       CZ => GRID_CZ, &
       FX => GRID_FX, &
       FY => GRID_FY, &
       FZ => GRID_FZ, &
       RCDZ => GRID_RCDZ, &
       RCDX => GRID_RCDX, &
       RCDY => GRID_RCDY
    use scale_atmos_hydrometeor, only: &
       I_QV
    use scale_atmos_thermodyn, only: &
       THERMODYN_qd => ATMOS_THERMODYN_qd, &
       THERMODYN_cp => ATMOS_THERMODYN_cp, &
       THERMODYN_r  => ATMOS_THERMODYN_r
    use scale_atmos_saturation, only : &
       moist_pres2qsat_liq => ATMOS_SATURATION_pres2qsat_liq, &
       moist_pres2qsat_ice => ATMOS_SATURATION_pres2qsat_ice

    implicit none

    real(RP) :: Uabs, Cm, Ch, Ce
    
    integer :: k, i, j, iq
    character(len=H_LONG) :: basename_sd_out
    integer :: fid_sdm_0
    integer :: fid_sdm_1
    integer :: fid_sdm_2
    integer :: fid_sdm_3
    integer :: fid_sdm_4
    integer :: ierr
    logical :: exist    

    ! surface flux for bottom/top boundary
    real(RP) :: SFLX_BT_VELX(IA,JA)
    real(RP) :: SFLX_BT_VELY(IA,JA)
    real(RP) :: SFLX_BT_POTT(IA,JA)
    real(RP) :: SFLX_BT_QV  (IA,JA)
    real(RP) :: SFLX_BT_SH  (IA,JA)

    ! surface flux for west/east boundary
    real(RP) :: SFLX_WE_VELY(KA,JA)
    real(RP) :: SFLX_WE_VELZ(KA,JA)
    real(RP) :: SFLX_WE_POTT(KA,JA)
    real(RP) :: SFLX_WE_QV  (KA,JA)
    real(RP) :: SFLX_WE_SH  (KA,JA)

    ! surface flux for south/north boundary
    real(RP) :: SFLX_SN_VELZ(KA,IA)
    real(RP) :: SFLX_SN_VELX(KA,IA)
    real(RP) :: SFLX_SN_POTT(KA,IA)
    real(RP) :: SFLX_SN_QV  (KA,IA)
    real(RP) :: SFLX_SN_SH  (KA,IA)

    ! surface flux for west/east/south/north
    real(RP) :: total_SFLX_W_SH
    real(RP) :: total_SFLX_W_QV
    real(RP) :: total_SFLX_E_SH
    real(RP) :: total_SFLX_E_QV
    real(RP) :: total_SFLX_S_SH
    real(RP) :: total_SFLX_S_QV
    real(RP) :: total_SFLX_N_SH
    real(RP) :: total_SFLX_N_QV
    real(RP) :: total_SFLX_T_SH
    real(RP) :: total_SFLX_T_QV    
    !integer  :: histitemid

    real(RP) :: q(QA)
    real(RP) :: qdry
    real(RP) :: Rtot
    real(RP) :: CPtot
    real(RP) :: CPovCV
    real(RP) :: RovCP

    real(RP),parameter :: k_karman=0.4
    real(RP) :: z0m ! momentum roughness length [m]
    real(RP) :: z0h ! sensivle heat roughness length [m]
    real(RP) :: z0e ! vapor roughness length [m]

    real(RP) :: WALL_TEMP, WALL_PRES, WALL_PT, WALL_QV
    real(RP) :: D_REFP

    z0m = 0.75e-3_RP   ! [m] specified in the ICMW2024 PI case description  
    z0h = 0.619*z0m    ! [m] specified in the ICMW2024 PI case description
    z0e = 0.756*z0m    ! [m] specified in the ICMW2024 PI case description  
  
    !---------------------------------------------------------------------------

    if ( USER_do ) then
      ! open the file0
      fid_sdm_0 = IO_get_available_fid()
      write(basename_sd_out,'(A22,I6.6)') 'Surface_flux_west.pe',mype

      if( IO_L ) write(IO_FID_LOG,*) '*** Atmos physics  step: Surface flux (USER_step)'
      
      inquire(file=trim(basename_sd_out), exist=exist)
      if (exist) then
         open(fid_sdm_0, file=trim(basename_sd_out), status="old", position="append", &
                        action="write", form="formatted", iostat=ierr)
      else
         open(fid_sdm_0, file=trim(basename_sd_out), status="new", position="append", &
                        action="write", form="formatted", iostat=ierr)
         write(fid_sdm_0, '(A)') ' time  total_SFLX_W_SH  total_SFLX_W_QV '
      end if

      ! open the file1
      fid_sdm_1 = IO_get_available_fid()
      write(basename_sd_out,'(A22,I6.6)') 'Surface_flux_east.pe',mype

      if( IO_L ) write(IO_FID_LOG,*) '*** Atmos physics  step: Surface flux (USER_step)'
      
      inquire(file=trim(basename_sd_out), exist=exist)
      if (exist) then
         open(fid_sdm_1, file=trim(basename_sd_out), status="old", position="append", &
                        action="write", form="formatted", iostat=ierr)
      else
         open(fid_sdm_1, file=trim(basename_sd_out), status="new", position="append", &
                        action="write", form="formatted", iostat=ierr)
         write(fid_sdm_1, '(A)') ' time   total_SFLX_E_SH  total_SFLX_E_QV '
      end if

      ! open the file2
      fid_sdm_2 = IO_get_available_fid()
      write(basename_sd_out,'(A22,I6.6)') 'Surface_flux_south.pe',mype

      if( IO_L ) write(IO_FID_LOG,*) '*** Atmos physics  step: Surface flux (USER_step)'
      
      inquire(file=trim(basename_sd_out), exist=exist)
      if (exist) then
         open(fid_sdm_2, file=trim(basename_sd_out), status="old", position="append", &
                        action="write", form="formatted", iostat=ierr)
      else
         open(fid_sdm_2, file=trim(basename_sd_out), status="new", position="append", &
                        action="write", form="formatted", iostat=ierr)
         write(fid_sdm_2, '(A)') ' time   total_SFLX_S_SH  total_SFLX_S_QV '
      end if

      ! open the file3
      fid_sdm_3 = IO_get_available_fid()
      write(basename_sd_out,'(A22,I6.6)') 'Surface_flux_north.pe',mype

      if( IO_L ) write(IO_FID_LOG,*) '*** Atmos physics  step: Surface flux (USER_step)'
      
      inquire(file=trim(basename_sd_out), exist=exist)
      if (exist) then
         open(fid_sdm_3, file=trim(basename_sd_out), status="old", position="append", &
                        action="write", form="formatted", iostat=ierr)
      else
         open(fid_sdm_3, file=trim(basename_sd_out), status="new", position="append", &
                        action="write", form="formatted", iostat=ierr)
         write(fid_sdm_3, '(A)') ' time   total_SFLX_N_SH  total_SFLX_N_QV '
      end if
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!! reset the verlocity boundary condition at the side walls
      if( .NOT. PRC_HAS_W )then !west boundary
         do i = 1,IS-1
            MOMX(:,i,:) = 0.0_RP
            MOMY(:,i,:) = MOMY(:,IS,:)
            MOMZ(:,i,:) = MOMZ(:,IS,:)
            DENS(:,i,:) = DENS(:,IS,:)
            RHOT(:,i,:) = RHOT(:,IS,:)
            QTRC(:,i,:,:) = QTRC(:,IS,:,:)
         end do
      endif

      if( .NOT. PRC_HAS_E )then !east boundary
         do i = IE,IA
            MOMX(:,i,:) = 0.0_RP
            MOMY(:,i,:) = MOMY(:,IE,:)
            MOMZ(:,i,:) = MOMZ(:,IE,:)
            DENS(:,i,:) = DENS(:,IE,:)
            RHOT(:,i,:) = RHOT(:,IE,:)
            QTRC(:,i,:,:) = QTRC(:,IE,:,:)
         end do
      endif

      if( .NOT. PRC_HAS_S )then !south boundary
         do j = 1,JS-1
            MOMX(:,:,j) = MOMX(:,:,JS)
            MOMY(:,:,j) = 0.0_RP
            MOMZ(:,:,j) = MOMZ(:,:,JS)
            DENS(:,:,j) = DENS(:,:,JS)
            RHOT(:,:,j) = RHOT(:,:,JS)
            QTRC(:,:,j,:) = QTRC(:,:,JS,:)
         end do
      endif

      if( .NOT. PRC_HAS_N )then !north boundary
         do j = JE,JA
            MOMX(:,:,j) = MOMX(:,:,JE)
            MOMY(:,:,j) = 0.0_RP
            MOMZ(:,:,j) = MOMZ(:,:,JE)
            MOMZ(:,:,j) = MOMZ(:,:,JE)
            DENS(:,:,j) = DENS(:,:,JE)
            RHOT(:,:,j) = RHOT(:,:,JE)
            QTRC(:,:,j,:) = QTRC(:,:,JE,:)
         end do
      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!! Bottom !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      k = KS
      WALL_TEMP=TEMP_BOTTOM
      D_REFP=CZ(KS)-FZ(KS-1)

      !! momx
      !$omp parallel do private(i,j,k,Uabs,Cm) collapse(2)
      do j = JS, JE
        do i = IS, IE
          ! at (u, y, layer)
          Uabs = sqrt( &
                  + ( 2.0_RP *   MOMX(k,i,j) )**2 &
                  + ( 0.5_RP * ( MOMY(k,i,j-1) + MOMY(k,i,j) + MOMY(k,i+1,j-1) + MOMY(k,i+1,j) ) )**2 &
                  ) / ( DENS(k,i,j) + DENS(k,i+1,j) )
          Cm = k_karman**2 / (log(D_REFP/z0m))**2
          Cm = min( Cm, Cm_max )

          SFLX_BT_VELX(i,j) = - 2 * Cm * Uabs * MOMX(k,i,j) / ( DENS(k,i,j) + DENS(k,i+1,j) )
          
          MOMX_tp(k,i,j) = MOMX_tp(k,i,j) + 0.5_RP*( DENS(k,i,j) + DENS(k,i+1,j) )*SFLX_BT_VELX(i,j)*RCDZ(k)
        enddo
      enddo
      !$omp end parallel do

      !! momy
      !$omp parallel do private(i,j,k,Uabs,Cm) collapse(2)
      do j = JS, JE
        do i = IS, IE
          ! at (x, v, layer)
          Uabs = sqrt( &
            + ( 0.5_RP * ( MOMX(k,i-1,j) + MOMX(k,i,j) + MOMX(k,i-1,j+1) + MOMX(k,i,j+1) ) )**2 &
            + ( 2.0_RP *   MOMY(k,i,j) )**2 &
            ) / ( DENS(k,i,j) + DENS(k,i,j+1) )
          Cm = k_karman**2 / (log(D_REFP/z0m))**2
          Cm = min( Cm, Cm_max )

          SFLX_BT_VELY(i,j) = - 2 * Cm * Uabs *MOMY(k,i,j) / ( DENS(k,i,j) + DENS(k,i,j+1) )

          MOMY_tp(k,i,j) = MOMY_tp(k,i,j) + 0.5_RP*( DENS(k,i,j) + DENS(k,i,j+1) )*SFLX_BT_VELY(i,j)*RCDZ(k)
             
        enddo
      enddo
      !$omp end parallel do

      !! rhot, dens, and rhoqv
      do j = JS, JE
        do i = IS, IE
          ! at cell center
          Uabs = sqrt( &
                  + ( MOMX(k,i-1,j) + MOMX(k,i,j) )**2 &
                  + ( MOMY(k,i,j-1) + MOMY(k,i,j) )**2 &
                  ) / DENS(k,i,j) * 0.5_RP

          Ch = k_karman**2 / (log(D_REFP/z0m)*log(D_REFP/z0h))
          Ce = k_karman**2 / (log(D_REFP/z0m)*log(D_REFP/z0e))

          !--- saturation at surface
          do iq = 1, QA
             q(iq) = QTRC(k,i,j,iq)
          enddo

          call THERMODYN_qd( qdry,  q, TRACER_MASS )
          call THERMODYN_cp( CPtot, q, TRACER_CP, qdry )
          call THERMODYN_r ( Rtot,  q, TRACER_R,  qdry )
          
          CPovCV = CPtot / ( CPtot - Rtot )
          WALL_PRES   = P00 * ( RHOT(k,i,j) * Rtot / P00 )**CPovCV
          
          RovCP = Rtot / CPtot
          WALL_PT = WALL_TEMP * ( P00 / WALL_PRES )**RovCP 

          call moist_pres2qsat_liq( WALL_QV, WALL_TEMP, WALL_PRES )

          SFLX_BT_POTT(i,j) = Ch * Uabs * ( WALL_PT - RHOT(k,i,j)/DENS(k,i,j) )
          SFLX_BT_QV  (i,j) = Ce * Uabs * ( WALL_QV - QTRC(k,i,j,I_QV) )

          RHOT_tp(k,i,j) = RHOT_tp(k,i,j) + ( DENS(k,i,j)*SFLX_BT_POTT(i,j) + SFLX_BT_QV(i,j)*RHOT(k,i,j) )*RCDZ(k)
          DENS_tp(k,i,j) = DENS_tp(k,i,j) + DENS(k,i,j)*SFLX_BT_QV(i,j) *RCDZ(k) 
          RHOQ_tp(k,i,j,I_QV) = RHOQ_tp(k,i,j,I_QV) + DENS(k,i,j)*SFLX_BT_QV(i,j)*RCDZ(k)

          SFLX_BT_SH(i,j) = CPtot*DENS(k,i,j)*SFLX_BT_POTT(i,j)

        enddo
      enddo
      
      call HIST_in( SFLX_BT_SH(:,:), 'SFLX_B_SH', 'surface sensible heat flux at the bottom', 'W/m2')
      call HIST_in( SFLX_BT_QV(:,:), 'SFLX_B_QV', 'surface moisuture flux at the bottom', 'kg/kg m/s')
   
      !!!! Top !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      k = KE
      WALL_TEMP=TEMP_TOP
      D_REFP=FZ(KE)-CZ(KE)

      !! momx
      !$omp parallel do private(i,j,k,Uabs,Cm) collapse(2)
      do j = JS, JE
        do i = IS, IE
          ! at (u, y, layer)
          Uabs = sqrt( &
                  + ( 2.0_RP *   MOMX(k,i,j) )**2 &
                  + ( 0.5_RP * ( MOMY(k,i,j-1) + MOMY(k,i,j) + MOMY(k,i+1,j-1) + MOMY(k,i+1,j) ) )**2 &
                  ) / ( DENS(k,i,j) + DENS(k,i+1,j) )
          Cm = k_karman**2 / (log(D_REFP/z0m))**2
          Cm = min( Cm, Cm_max )

          SFLX_BT_VELX(i,j) = - 2 * Cm * Uabs * MOMX(k,i,j) / ( DENS(k,i,j) + DENS(k,i+1,j) )
          
          MOMX_tp(k,i,j) = MOMX_tp(k,i,j) + 0.5_RP*( DENS(k,i,j) + DENS(k,i+1,j) )*SFLX_BT_VELX(i,j)*RCDZ(k)
        enddo
      enddo
      !$omp end parallel do

      !! momy
      !$omp parallel do private(i,j,k,Uabs,Cm) collapse(2)
      do j = JS, JE
        do i = IS, IE
          ! at (x, v, layer)
          Uabs = sqrt( &
            + ( 0.5_RP * ( MOMX(k,i-1,j) + MOMX(k,i,j) + MOMX(k,i-1,j+1) + MOMX(k,i,j+1) ) )**2 &
            + ( 2.0_RP *   MOMY(k,i,j) )**2 &
            ) / ( DENS(k,i,j) + DENS(k,i,j+1) )
          Cm = k_karman**2 / (log(D_REFP/z0m))**2
          Cm = min( Cm, Cm_max )

          SFLX_BT_VELY(i,j) = - 2 * Cm * Uabs *MOMY(k,i,j) / ( DENS(k,i,j) + DENS(k,i,j+1) )

          MOMY_tp(k,i,j) = MOMY_tp(k,i,j) + 0.5_RP*( DENS(k,i,j) + DENS(k,i,j+1) )*SFLX_BT_VELY(i,j)*RCDZ(k)
             
        enddo
      enddo
      !$omp end parallel do

      !! rhot, dens, and rhoqv
      do j = JS, JE
        do i = IS, IE
          ! at cell center
          Uabs = sqrt( &
                  + ( MOMX(k,i-1,j) + MOMX(k,i,j) )**2 &
                  + ( MOMY(k,i,j-1) + MOMY(k,i,j) )**2 &
                  ) / DENS(k,i,j) * 0.5_RP

          Ch = k_karman**2 / (log(D_REFP/z0m)*log(D_REFP/z0h))
          Ce = k_karman**2 / (log(D_REFP/z0m)*log(D_REFP/z0e))

          !--- saturation at surface
          do iq = 1, QA
             q(iq) = QTRC(k,i,j,iq)
          enddo

          call THERMODYN_qd( qdry,  q, TRACER_MASS )
          call THERMODYN_cp( CPtot, q, TRACER_CP, qdry )
          call THERMODYN_r ( Rtot,  q, TRACER_R,  qdry )
          
          CPovCV = CPtot / ( CPtot - Rtot )
          WALL_PRES   = P00 * ( RHOT(k,i,j) * Rtot / P00 )**CPovCV
          
          RovCP = Rtot / CPtot
          WALL_PT = WALL_TEMP * ( P00 / WALL_PRES )**RovCP 

          call moist_pres2qsat_ice( WALL_QV, WALL_TEMP, WALL_PRES )

          SFLX_BT_POTT(i,j) = Ch * Uabs * ( WALL_PT - RHOT(k,i,j)/DENS(k,i,j) )
          SFLX_BT_QV  (i,j) = Ce * Uabs * ( WALL_QV - QTRC(k,i,j,I_QV) )

          RHOT_tp(k,i,j) = RHOT_tp(k,i,j) + ( DENS(k,i,j)*SFLX_BT_POTT(i,j) + SFLX_BT_QV(i,j)*RHOT(k,i,j) )*RCDZ(k)
          DENS_tp(k,i,j) = DENS_tp(k,i,j) + DENS(k,i,j)*SFLX_BT_QV(i,j) *RCDZ(k) 
          RHOQ_tp(k,i,j,I_QV) = RHOQ_tp(k,i,j,I_QV) + DENS(k,i,j)*SFLX_BT_QV(i,j)*RCDZ(k)

          SFLX_BT_SH(i,j) = CPtot*DENS(k,i,j)*SFLX_BT_POTT(i,j)

        enddo
      enddo

      call HIST_in( SFLX_BT_SH(:,:), 'SFLX_T_SH', 'surface sensible heat flux at the top', 'W/m2')
      call HIST_in( SFLX_BT_QV(:,:), 'SFLX_T_QV', 'surface moisuture flux at the top', 'kg/kg m/s')

      !!!! West !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SFLX_WE_SH(:,:) = 0.0_RP
      SFLX_WE_QV(:,:) = 0.0_RP

      if( .NOT. PRC_HAS_W )then !west boundary

         i = IS
         WALL_TEMP=TEMP_SIDE
         D_REFP=CX(IS)-FX(IS-1)

         !! momy
         !$omp parallel do private(i,j,k,Uabs,Cm) collapse(2)
         do j = JS, JE
            do k = KS, KE
               ! at (v, z, layer)
               Uabs = sqrt( &
                    + ( 2.0_RP *   MOMY(k,i,j) )**2 &
                    + ( 0.5_RP * ( MOMZ(k-1,i,j) + MOMZ(k,i,j) + MOMZ(k-1,i,j+1) + MOMZ(k,i,j+1) ) )**2 &
                    ) / ( DENS(k,i,j) + DENS(k,i,j+1) )
               Cm = k_karman**2 / (log(D_REFP/z0m))**2
               Cm = min( Cm, Cm_max )
               
               SFLX_WE_VELY(k,j) = - 2 * Cm * Uabs * MOMY(k,i,j) / ( DENS(k,i,j) + DENS(k,i,j+1) )
          
               MOMY_tp(k,i,j) = MOMY_tp(k,i,j) + 0.5_RP*( DENS(k,i,j) + DENS(k,i,j+1) )*SFLX_WE_VELY(k,j)*RCDX(i)
            enddo
         enddo
         !$omp end parallel do
         
         !! momz
         !$omp parallel do private(i,j,k,Uabs,Cm) collapse(2)
         do j = JS, JE
            do k = KS, KE
               ! at (y, w, layer)
               Uabs = sqrt( &
                    + ( 0.5_RP * ( MOMY(k,i,j-1) + MOMY(k,i,j) + MOMY(k+1,i,j-1) + MOMY(k+1,i,j) ) )**2 &
                    + ( 2.0_RP *   MOMZ(k,i,j) )**2 &
                    ) / ( DENS(k,i,j) + DENS(k+1,i,j) )
               Cm = k_karman**2 / (log(D_REFP/z0m))**2
               Cm = min( Cm, Cm_max )
               
               SFLX_WE_VELZ(k,j) = - 2 * Cm * Uabs *MOMZ(k,i,j) / ( DENS(k,i,j) + DENS(k+1,i,j) )
               
               MOMZ_tp(k,i,j) = MOMZ_tp(k,i,j) + 0.5_RP*( DENS(k,i,j) + DENS(k+1,i,j) )*SFLX_WE_VELZ(k,j)*RCDX(i)           
            enddo
         enddo
         !$omp end parallel do

         !! rhot, dens, and rhoqv
         do j = JS, JE
            do k = KS, KE
               ! at cell center
               Uabs = sqrt( &
                    + ( MOMY(k,i,j-1) + MOMY(k,i,j) )**2 &
                    + ( MOMZ(k-1,i,j) + MOMZ(k,i,j) )**2 &
                    ) / DENS(k,i,j) * 0.5_RP

               Ch = k_karman**2 / (log(D_REFP/z0m)*log(D_REFP/z0h))
               Ce = k_karman**2 / (log(D_REFP/z0m)*log(D_REFP/z0e))

               !--- saturation at surface
               do iq = 1, QA
                  q(iq) = QTRC(k,i,j,iq)
               enddo

               call THERMODYN_qd( qdry,  q, TRACER_MASS )
               call THERMODYN_cp( CPtot, q, TRACER_CP, qdry )
               call THERMODYN_r ( Rtot,  q, TRACER_R,  qdry )
          
               CPovCV = CPtot / ( CPtot - Rtot )
               WALL_PRES   = P00 * ( RHOT(k,i,j) * Rtot / P00 )**CPovCV
               
               RovCP = Rtot / CPtot
               WALL_PT = WALL_TEMP * ( P00 / WALL_PRES )**RovCP 

               call moist_pres2qsat_ice( WALL_QV, WALL_TEMP, WALL_PRES )
               WALL_QV = WALL_QV * RHICE_SIDE

               SFLX_WE_POTT(k,j) = Ch * Uabs * ( WALL_PT - RHOT(k,i,j)/DENS(k,i,j) )
               SFLX_WE_QV  (k,j) = Ce * Uabs * ( WALL_QV - QTRC(k,i,j,I_QV) )

               RHOT_tp(k,i,j) = RHOT_tp(k,i,j) + ( DENS(k,i,j)*SFLX_WE_POTT(k,j) + SFLX_WE_QV(k,j)*RHOT(k,i,j) )*RCDX(i)
               DENS_tp(k,i,j) = DENS_tp(k,i,j) + DENS(k,i,j)*SFLX_WE_QV(k,j) *RCDX(i)
               RHOQ_tp(k,i,j,I_QV) = RHOQ_tp(k,i,j,I_QV) + DENS(k,i,j)*SFLX_WE_QV(k,j)*RCDX(i)

               SFLX_WE_SH(k,j) = CPtot*DENS(k,i,j)*SFLX_WE_POTT(k,j)

            enddo
         enddo
         total_SFLX_W_SH = sum(SFLX_WE_SH(KS+4:KE-4,:))
         total_SFLX_W_QV = sum(SFLX_WE_QV(KS+4:KE-4,:))
         if (modulo(TIME_NOWSEC, 30.0) == 0.0) then
            write(fid_sdm_0,'(F10.1, E30.14, E30.14)') TIME_NOWSEC, total_SFLX_W_SH, total_SFLX_W_QV
         end if
      end if

      !call HIST_in( SFLX_WE_SH(:,:), 'SFLX_W_SH', 'surface sensible heat flux at the west side', 'W/m2')
      !call HIST_in( SFLX_WE_QV(:,:), 'SFLX_W_QV', 'surface moisuture flux at the west side', 'kg/kg m/s')

      !!!! East !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SFLX_WE_SH(:,:) = 0.0_RP
      SFLX_WE_QV(:,:) = 0.0_RP

      if( .NOT. PRC_HAS_E )then !east boundary

         i = IE
         WALL_TEMP=TEMP_SIDE
         D_REFP=FX(IE)-CX(IE)

         !! momy
         !$omp parallel do private(i,j,k,Uabs,Cm) collapse(2)
         do j = JS, JE
            do k = KS, KE
               ! at (v, z, layer)
               Uabs = sqrt( &
                    + ( 2.0_RP *   MOMY(k,i,j) )**2 &
                    + ( 0.5_RP * ( MOMZ(k-1,i,j) + MOMZ(k,i,j) + MOMZ(k-1,i,j+1) + MOMZ(k,i,j+1) ) )**2 &
                    ) / ( DENS(k,i,j) + DENS(k,i,j+1) )
               Cm = k_karman**2 / (log(D_REFP/z0m))**2
               Cm = min( Cm, Cm_max )
               
               SFLX_WE_VELY(k,j) = - 2 * Cm * Uabs * MOMY(k,i,j) / ( DENS(k,i,j) + DENS(k,i,j+1) )
          
               MOMY_tp(k,i,j) = MOMY_tp(k,i,j) + 0.5_RP*( DENS(k,i,j) + DENS(k,i,j+1) )*SFLX_WE_VELY(k,j)*RCDX(i)
            enddo
         enddo
         !$omp end parallel do
         
         !! momz
         !$omp parallel do private(i,j,k,Uabs,Cm) collapse(2)
         do j = JS, JE
            do k = KS, KE
               ! at (y, w, layer)
               Uabs = sqrt( &
                    + ( 0.5_RP * ( MOMY(k,i,j-1) + MOMY(k,i,j) + MOMY(k+1,i,j-1) + MOMY(k+1,i,j) ) )**2 &
                    + ( 2.0_RP *   MOMZ(k,i,j) )**2 &
                    ) / ( DENS(k,i,j) + DENS(k+1,i,j) )
               Cm = k_karman**2 / (log(D_REFP/z0m))**2
               Cm = min( Cm, Cm_max )
               
               SFLX_WE_VELZ(k,j) = - 2 * Cm * Uabs *MOMZ(k,i,j) / ( DENS(k,i,j) + DENS(k+1,i,j) )
               
               MOMZ_tp(k,i,j) = MOMZ_tp(k,i,j) + 0.5_RP*( DENS(k,i,j) + DENS(k+1,i,j) )*SFLX_WE_VELZ(k,j)*RCDX(i)           
            enddo
         enddo
         !$omp end parallel do

         !! rhot, dens, and rhoqv
         do j = JS, JE
            do k = KS, KE
               ! at cell center
               Uabs = sqrt( &
                    + ( MOMY(k,i,j-1) + MOMY(k,i,j) )**2 &
                    + ( MOMZ(k-1,i,j) + MOMZ(k,i,j) )**2 &
                    ) / DENS(k,i,j) * 0.5_RP

               Ch = k_karman**2 / (log(D_REFP/z0m)*log(D_REFP/z0h))
               Ce = k_karman**2 / (log(D_REFP/z0m)*log(D_REFP/z0e))

               !--- saturation at surface
               do iq = 1, QA
                  q(iq) = QTRC(k,i,j,iq)
               enddo

               call THERMODYN_qd( qdry,  q, TRACER_MASS )
               call THERMODYN_cp( CPtot, q, TRACER_CP, qdry )
               call THERMODYN_r ( Rtot,  q, TRACER_R,  qdry )
          
               CPovCV = CPtot / ( CPtot - Rtot )
               WALL_PRES   = P00 * ( RHOT(k,i,j) * Rtot / P00 )**CPovCV
               
               RovCP = Rtot / CPtot
               WALL_PT = WALL_TEMP * ( P00 / WALL_PRES )**RovCP 

               call moist_pres2qsat_ice( WALL_QV, WALL_TEMP, WALL_PRES )
               WALL_QV = WALL_QV * RHICE_SIDE

               SFLX_WE_POTT(k,j) = Ch * Uabs * ( WALL_PT - RHOT(k,i,j)/DENS(k,i,j) )
               SFLX_WE_QV  (k,j) = Ce * Uabs * ( WALL_QV - QTRC(k,i,j,I_QV) )

               RHOT_tp(k,i,j) = RHOT_tp(k,i,j) + ( DENS(k,i,j)*SFLX_WE_POTT(k,j) + SFLX_WE_QV(k,j)*RHOT(k,i,j) )*RCDX(i)
               DENS_tp(k,i,j) = DENS_tp(k,i,j) + DENS(k,i,j)*SFLX_WE_QV(k,j) *RCDX(i)
               RHOQ_tp(k,i,j,I_QV) = RHOQ_tp(k,i,j,I_QV) + DENS(k,i,j)*SFLX_WE_QV(k,j)*RCDX(i)

               SFLX_WE_SH(k,j) = CPtot*DENS(k,i,j)*SFLX_WE_POTT(k,j)

            enddo
         enddo
         total_SFLX_E_SH = sum(SFLX_WE_SH(KS+4:KE-4,:))
         total_SFLX_E_QV = sum(SFLX_WE_QV(KS+4:KE-4,:))
         if (modulo(TIME_NOWSEC, 30.0) == 0.0) then
            write(fid_sdm_1,'(F10.1, E30.14, E30.14)') TIME_NOWSEC, total_SFLX_E_SH, total_SFLX_E_QV
         end if
      end if
      !call HIST_in( SFLX_WE_SH(:,:), 'SFLX_E_SH', 'surface sensible heat flux at the east side', 'W/m2')
      !call HIST_in( SFLX_WE_QV(:,:), 'SFLX_E_QV', 'surface moisuture flux at the east side', 'kg/kg m/s')

      !!!! South !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SFLX_SN_SH(:,:) = 0.0_RP
      SFLX_SN_QV(:,:) = 0.0_RP

      if( .NOT. PRC_HAS_S )then !south boundary

         j = JS
         WALL_TEMP=TEMP_SIDE
         D_REFP=CY(JS)-FY(JS-1)

         !! momz
         !$omp parallel do private(i,j,k,Uabs,Cm) collapse(2)
         do i = IS, IE
            do k = KS, KE
               ! at (w, x, layer)
               Uabs = sqrt( &
                    + ( 2.0_RP *   MOMZ(k,i,j) )**2 &
                    + ( 0.5_RP * ( MOMX(k,i-1,j) + MOMZ(k,i,j) + MOMZ(k+1,i-1,j) + MOMZ(k+1,i,j) ) )**2 &
                    ) / ( DENS(k+1,i,j) + DENS(k,i,j) )
               Cm = k_karman**2 / (log(D_REFP/z0m))**2
               Cm = min( Cm, Cm_max )
               
               SFLX_SN_VELZ(k,i) = - 2 * Cm * Uabs * MOMZ(k,i,j) / ( DENS(k,i,j) + DENS(k+1,i,j) )
          
               MOMZ_tp(k,i,j) = MOMZ_tp(k,i,j) + 0.5_RP*( DENS(k,i,j) + DENS(k+1,i,j) )*SFLX_SN_VELZ(k,i)*RCDY(j)
            enddo
         enddo
         !$omp end parallel do
         
         !! momx
         !$omp parallel do private(i,j,k,Uabs,Cm) collapse(2)
         do i = IS, IE
            do k = KS, KE
               ! at (z, u, layer)
               Uabs = sqrt( &
                    + ( 0.5_RP * ( MOMZ(k-1,i,j) + MOMZ(k,i,j) + MOMZ(k-1,i+1,j) + MOMZ(k,i+1,j) ) )**2 &
                    + ( 2.0_RP *   MOMX(k,i,j) )**2 &
                    ) / ( DENS(k,i,j) + DENS(k,i+1,j) )
               Cm = k_karman**2 / (log(D_REFP/z0m))**2
               Cm = min( Cm, Cm_max )
               
               SFLX_SN_VELX(k,i) = - 2 * Cm * Uabs *MOMX(k,i,j) / ( DENS(k,i,j) + DENS(k,i+1,j) )
               
               MOMX_tp(k,i,j) = MOMX_tp(k,i,j) + 0.5_RP*( DENS(k,i,j) + DENS(k,i+1,j) )*SFLX_SN_VELX(k,i)*RCDY(j)
            enddo
         enddo
         !$omp end parallel do

         !! rhot, dens, and rhoqv
         do i = IS, IE
            do k = KS, KE
               ! at cell center
               Uabs = sqrt( &
                    + ( MOMZ(k-1,i,j) + MOMZ(k,i,j) )**2 &
                    + ( MOMX(k,i-1,j) + MOMX(k,i,j) )**2 &
                    ) / DENS(k,i,j) * 0.5_RP

               Ch = k_karman**2 / (log(D_REFP/z0m)*log(D_REFP/z0h))
               Ce = k_karman**2 / (log(D_REFP/z0m)*log(D_REFP/z0e))

               !--- saturation at surface
               do iq = 1, QA
                  q(iq) = QTRC(k,i,j,iq)
               enddo

               call THERMODYN_qd( qdry,  q, TRACER_MASS )
               call THERMODYN_cp( CPtot, q, TRACER_CP, qdry )
               call THERMODYN_r ( Rtot,  q, TRACER_R,  qdry )
          
               CPovCV = CPtot / ( CPtot - Rtot )
               WALL_PRES   = P00 * ( RHOT(k,i,j) * Rtot / P00 )**CPovCV
               
               RovCP = Rtot / CPtot
               WALL_PT = WALL_TEMP * ( P00 / WALL_PRES )**RovCP 

               call moist_pres2qsat_ice( WALL_QV, WALL_TEMP, WALL_PRES )
               WALL_QV = WALL_QV * RHICE_SIDE

               SFLX_SN_POTT(k,i) = Ch * Uabs * ( WALL_PT - RHOT(k,i,j)/DENS(k,i,j) )
               SFLX_SN_QV  (k,i) = Ce * Uabs * ( WALL_QV - QTRC(k,i,j,I_QV) )

               RHOT_tp(k,i,j) = RHOT_tp(k,i,j) + ( DENS(k,i,j)*SFLX_SN_POTT(k,i) + SFLX_SN_QV(k,i)*RHOT(k,i,j) )*RCDY(j)
               DENS_tp(k,i,j) = DENS_tp(k,i,j) + DENS(k,i,j)*SFLX_SN_QV(k,i) *RCDY(j)
               RHOQ_tp(k,i,j,I_QV) = RHOQ_tp(k,i,j,I_QV) + DENS(k,i,j)*SFLX_SN_QV(k,i)*RCDY(j)

               SFLX_SN_SH(k,i) = CPtot*DENS(k,i,j)*SFLX_SN_POTT(k,i)

            enddo
         enddo
         total_SFLX_S_SH = sum(SFLX_SN_SH(KS+4:KE-4,:))
         total_SFLX_S_QV = sum(SFLX_SN_QV(KS+4:KE-4,:))
         if (modulo(TIME_NOWSEC, 30.0) == 0.0) then
            write(fid_sdm_2,'(F10.1, E30.14, E30.14)') TIME_NOWSEC, total_SFLX_S_SH, total_SFLX_S_QV
         end if
      end if
      !call HIST_in( SFLX_SN_SH(:,:), 'SFLX_S_SH', 'surface sensible heat flux at the south side', 'W/m2')
      !call HIST_in( SFLX_SN_QV(:,:), 'SFLX_S_QV', 'surface moisuture flux at the south side', 'kg/kg m/s')

      !!!! North !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SFLX_SN_SH(:,:) = 0.0_RP
      SFLX_SN_QV(:,:) = 0.0_RP

      if( .NOT. PRC_HAS_N )then !north boundary

         j = JE
         WALL_TEMP=TEMP_SIDE
         D_REFP=FY(JE)-CY(JE)

         !! momz
         !$omp parallel do private(i,j,k,Uabs,Cm) collapse(2)
         do i = IS, IE
            do k = KS, KE
               ! at (w, x, layer)
               Uabs = sqrt( &
                    + ( 2.0_RP *   MOMZ(k,i,j) )**2 &
                    + ( 0.5_RP * ( MOMX(k,i-1,j) + MOMZ(k,i,j) + MOMZ(k+1,i-1,j) + MOMZ(k+1,i,j) ) )**2 &
                    ) / ( DENS(k+1,i,j) + DENS(k,i,j) )
               Cm = k_karman**2 / (log(D_REFP/z0m))**2
               Cm = min( Cm, Cm_max )
               
               SFLX_SN_VELZ(k,i) = - 2 * Cm * Uabs * MOMZ(k,i,j) / ( DENS(k,i,j) + DENS(k+1,i,j) )
          
               MOMZ_tp(k,i,j) = MOMZ_tp(k,i,j) + 0.5_RP*( DENS(k,i,j) + DENS(k+1,i,j) )*SFLX_SN_VELZ(k,i)*RCDY(j)
            enddo
         enddo
         !$omp end parallel do
         
         !! momx
         !$omp parallel do private(i,j,k,Uabs,Cm) collapse(2)
         do i = IS, IE
            do k = KS, KE
               ! at (z, u, layer)
               Uabs = sqrt( &
                    + ( 0.5_RP * ( MOMZ(k-1,i,j) + MOMZ(k,i,j) + MOMZ(k-1,i+1,j) + MOMZ(k,i+1,j) ) )**2 &
                    + ( 2.0_RP *   MOMX(k,i,j) )**2 &
                    ) / ( DENS(k,i,j) + DENS(k,i+1,j) )
               Cm = k_karman**2 / (log(D_REFP/z0m))**2
               Cm = min( Cm, Cm_max )
               
               SFLX_SN_VELX(k,i) = - 2 * Cm * Uabs *MOMX(k,i,j) / ( DENS(k,i,j) + DENS(k,i+1,j) )
               
               MOMX_tp(k,i,j) = MOMX_tp(k,i,j) + 0.5_RP*( DENS(k,i,j) + DENS(k,i+1,j) )*SFLX_SN_VELX(k,i)*RCDY(j)
            enddo
         enddo
         !$omp end parallel do

         !! rhot, dens, and rhoqv
         do i = IS, IE
            do k = KS, KE
               ! at cell center
               Uabs = sqrt( &
                    + ( MOMZ(k-1,i,j) + MOMZ(k,i,j) )**2 &
                    + ( MOMX(k,i-1,j) + MOMX(k,i,j) )**2 &
                    ) / DENS(k,i,j) * 0.5_RP

               Ch = k_karman**2 / (log(D_REFP/z0m)*log(D_REFP/z0h))
               Ce = k_karman**2 / (log(D_REFP/z0m)*log(D_REFP/z0e))

               !--- saturation at surface
               do iq = 1, QA
                  q(iq) = QTRC(k,i,j,iq)
               enddo

               call THERMODYN_qd( qdry,  q, TRACER_MASS )
               call THERMODYN_cp( CPtot, q, TRACER_CP, qdry )
               call THERMODYN_r ( Rtot,  q, TRACER_R,  qdry )
          
               CPovCV = CPtot / ( CPtot - Rtot )
               WALL_PRES   = P00 * ( RHOT(k,i,j) * Rtot / P00 )**CPovCV
               
               RovCP = Rtot / CPtot
               WALL_PT = WALL_TEMP * ( P00 / WALL_PRES )**RovCP 

               call moist_pres2qsat_ice( WALL_QV, WALL_TEMP, WALL_PRES )
               WALL_QV = WALL_QV * RHICE_SIDE

               SFLX_SN_POTT(k,i) = Ch * Uabs * ( WALL_PT - RHOT(k,i,j)/DENS(k,i,j) )
               SFLX_SN_QV  (k,i) = Ce * Uabs * ( WALL_QV - QTRC(k,i,j,I_QV) )

               RHOT_tp(k,i,j) = RHOT_tp(k,i,j) + ( DENS(k,i,j)*SFLX_SN_POTT(k,i) + SFLX_SN_QV(k,i)*RHOT(k,i,j) )*RCDY(j)
               DENS_tp(k,i,j) = DENS_tp(k,i,j) + DENS(k,i,j)*SFLX_SN_QV(k,i) *RCDY(j)
               RHOQ_tp(k,i,j,I_QV) = RHOQ_tp(k,i,j,I_QV) + DENS(k,i,j)*SFLX_SN_QV(k,i)*RCDY(j)

               SFLX_SN_SH(k,i) = CPtot*DENS(k,i,j)*SFLX_SN_POTT(k,i)

            enddo
         enddo
         total_SFLX_N_SH = sum(SFLX_SN_SH(KS+4:KE-4,:))
         total_SFLX_N_QV = sum(SFLX_SN_QV(KS+4:KE-4,:))
         if (modulo(TIME_NOWSEC, 30.0) == 0.0) then
            write(fid_sdm_3,'(F10.1, E30.14, E30.14)') TIME_NOWSEC, total_SFLX_N_SH, total_SFLX_N_QV
         end if
      end if
      !call HIST_in( SFLX_SN_SH(:,:), 'SFLX_N_SH', 'surface sensible heat flux at the north side', 'W/m2')
      !call HIST_in( SFLX_SN_QV(:,:), 'SFLX_N_QV', 'surface moisuture flux at the north side', 'kg/kg m/s')
      close(fid_sdm_0)
      close(fid_sdm_1)
      close(fid_sdm_2)
      close(fid_sdm_3)   
   endif
   return
  end subroutine USER_step

end module mod_user
