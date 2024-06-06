MODULE GLOBAL
    USE PARAM
    IMPLICIT NONE
    SAVE

# if defined (PARALLEL)
! MPI variables
    INTEGER :: myid,ier
    INTEGER :: comm2d
    INTEGER :: n_west, n_east, n_suth, n_nrth
    INTEGER :: npx,npy
    INTEGER :: ndims=2
    INTEGER :: NumberProcessor
    INTEGER,DIMENSION(:,:),ALLOCATABLE :: ProcessorID 
    INTEGER, DIMENSION(2) :: dims, coords
    LOGICAL, DIMENSION(2) :: periods
    LOGICAL :: reorder = .true.
    INTEGER :: nprocs                 
    INTEGER :: iista, iiend, jjsta, jjend 
# endif
    INTEGER :: px,py

! TIMER
    REAL(SP) :: tbegin,tend
    INTEGER :: ISTAGE


! CHARACTER AND LOGICAL VARIABLES
    CHARACTER(LEN=80) INPUT_FILE_NAME
    CHARACTER(LEN=80) TITLE
    CHARACTER(LEN=80) ReadType
    CHARACTER(LEN=80) DEPTH_TYPE
    CHARACTER(LEN=80) RESULT_FOLDER
    CHARACTER(LEN=80) WaveMaker
    CHARACTER(LEN=80) Time_Scheme
    CHARACTER(LEN=80) CONSTR
    CHARACTER(LEN=80) HIGH_ORDER
    CHARACTER(LEN=80) FIELD_IO_TYPE 
    LOGICAL :: NO_MASK_FILE = .TRUE.
    LOGICAL :: NO_UV_FILE = .TRUE.
    CHARACTER(LEN=16)::FORMAT_LEN=' '


    REAL(SP),PARAMETER,DIMENSION(3)::alpha=(/0.0_SP,3.0_SP/4.0_SP,1.0_SP/3.0_SP/)
    REAL(SP),PARAMETER,DIMENSION(3)::beta=(/1.0_SP,1.0_SP/4.0_SP,2.0_SP/3.0_SP/)

    ! MUSCL Parameter
    REAL(SP)::Kappa  

! PHYSICAL DOMAIN VARIABLES
    INTEGER :: Mglob,Nglob,Mloc,Nloc,Mloc1,Nloc1
    INTEGER, PARAMETER :: Nghost = 3
    INTEGER :: Ibeg,Iend,Jbeg,Jend,Iend1,Jend1
    REAL(SP):: DX,DY
    ! Wavemaker Coordinates
    REAL(SP) :: DXg,DYg  

! TIME DOMAIN VARIABLES
    REAL(SP)::  DT
    REAL(SP)::  TIME
    REAL(SP)::  TOTAL_TIME
    REAL(SP)::  PLOT_INTV
    REAL(SP)::  PLOT_COUNT
    REAL(SP)::  SCREEN_INTV
    REAL(SP)::  SCREEN_COUNT
    REAL(SP)::  DT_fixed
    REAL(SP)::  PLOT_START_TIME
    

! OUTPUT `var_XXXXX` X count
    INTEGER :: icount=-1  ! for output file number
    
! NUMERICS 
    REAL(SP) :: MinDepth=0.001
    REAL(SP) :: MinDepthFrc=0.5
    REAL(SP) :: ArrTimeMin=0.001
    REAL(SP) :: CFL=0.15_SP
    REAL(SP) :: FroudeCap=10.0_SP
    REAL(SP) :: SWE_ETA_DEP=0.7_SP

    LOGICAL :: BATHY_CORRECTION = .FALSE.
    LOGICAL :: BREAKWATER = .FALSE.

! CARTESIDAN COORDINATE
   REAL(SP),DIMENSION(:),ALLOCATABLE ::  Xco,Yco

! DISPERSION
    ! gamma1 is for extra M term
    REAL(SP) :: gamma1
    REAL(SP) :: gamma2
    REAL(SP) :: gamma3 ! linear shallow water
    REAL(SP) :: Beta_ref=-0.531_SP
    ! KENNEDY EQUATION
    REAL(SP) :: Beta_1,Beta_2

    ! a1=beta_ref^2/2-1/6, a2=beta_ref+1/2, b1=beta_ref^2, b2=beta_ref
    REAL(SP) :: a1,a2,b1,b2
    LOGICAL :: DISPERSION=.FALSE.
    LOGICAL :: DISP_TIME_LEFT=.FALSE. 
    LOGICAL :: StretchGrid = .FALSE.
    ! put time derivative dispersion term on left
    LOGICAL :: FIXED_DT = .FALSE.

! LOCAL VARIABLES (?) TODO
    REAL(SP),DIMENSION(:,:),ALLOCATABLE :: &
       U4xL,U4xR,V4yL,V4yR, &
       V4xL,V4xR,U4yL,U4yR, &   
       U4,V4,U1p,V1p, &
       U1pp,V1pp, &
       U2,V2,U3,V3, &
       DelxU,DelxHU,DelxV,DelxEtar,&
       DelxHV, DelyHU, &
       DelyU,DelyHV,DelyV,DelyEtar,&
       UxL,UxR,VxL,VxR,&
       HUxL,HUxR,HUyL,HUyR,HxL,HxR, &
       EtaRxL,EtaRxR,&
       UyL,UyR,VyL,VyR,&
       HVxL,HVxR,HVyL,HVyR,HyL,HyR, &
       EtaRyL,EtaRyR, &
       PL,PR,QL,QR, &
       FxL,FxR,FyL,FyR, &
       GxL,GxR,GyL,GyR, &
       SxL,SxR,SyL,SyR, &

! DERIVATIVES
    ! cross-derivatives in space
       Vxy,DVxy,Uxy,DUxy, &
    ! second-derivatives in space
       Uxx,DUxx,Vyy,DVyy, &
    ! first-derivatives in space
       Ux,Vx,Uy,Vy,DUx,DUy,DVx,DVy, &
       ETAx,ETAy, ETAT, ETATx,ETATy, &
    ! time-derivatives
       U0,V0,Ut,Vt,Utx,Vty,Utxx,Utxy,Vtxy,Vtyy,&
       DUtxx,DUtxy,DVtxy,DVtyy,DUtx,DVty,&
    ! original variables
       Fx,Fy,U,V,HU,HV,&
       Gx,Gy,P,Q,SourceX,SourceY,Int2Flo, &
       tmp4preview,HeightMax,HeightMin,VelocityMax,&
       MomentumFluxMax,VorticityMax,ARRTIME

! WETTING AND DRYING (TODO)
     INTEGER,DIMENSION(:,:),ALLOCATABLE :: MASK,MASK_STRUC,MASK9
     REAL(SP) :: Dmass,WetArea,DwetEta

! WAVEMAKER
     LOGICAL :: EqualEnergy = .FALSE.
     LOGICAL :: WaveMakerCurrentBalance = .FALSE.
     REAL(SP) :: WaveMakerCd
     REAL(SP)::AMP_SOLI,DEP_SOLI,LAG_SOLI, CPH_SOLI,XWAVEMAKER, &
               Xc,Yc, WID,Xc_WK,Tperiod,AMP_WK,DEP_WK,Theta_WK, &
               rlamda,Time_ramp,D_gen,Beta_gen,Width_WK,Delta_WK,&
               Ywidth_WK,Yc_WK
     LOGICAL :: SolitaryPositiveDirection = .TRUE.
     REAL(SP),DIMENSION(:,:),ALLOCATABLE :: D_gen_ir,rlamda_ir,phase_ir
     REAL(SP),DIMENSION(:),ALLOCATABLE :: Beta_gen_ir,omgn_ir
     REAL(SP),DIMENSION(:),ALLOCATABLE :: omgn2D
     REAL(SP),DIMENSION(:,:),ALLOCATABLE :: Wavemaker_Mass
     REAL(SP) :: FreqMin,FreqMax,FreqPeak,GammaTMA,Hmo,ThetaPeak,&
                 Sigma_Theta
     REAL(SP),DIMENSION(:,:,:),ALLOCATABLE ::Cm,Sm
     INTEGER :: Nfreq,Ntheta 
     REAL(SP)::x1_Nwave = 5.0, &
               x2_Nwave = 5.0, &
               a0_Nwave = 1.0, &
               gamma_Nwave = -3.0, &
               dep_Nwave = 1.0

! FRICTION TERMS
     LOGICAL :: IN_Cd=.FALSE.
     REAL(SP):: Cd_fixed
     CHARACTER(LEN=80) CD_FILE
     REAL(SP),DIMENSION(:,:),ALLOCATABLE :: Cd

     LOGICAL :: ROLLER=.FALSE.
     REAL(SP) :: ROLLER_SWITCH 

! SPONGELAYER
     REAL(SP),DIMENSION(:,:),ALLOCATABLE :: SPONGE,SpongeMaker
     REAL(SP)::Sponge_west_width,Sponge_east_width, &
               Sponge_south_width,Sponge_north_width, &
               R_sponge,A_sponge


! smagorinsky and wave height
   REAL(SP),DIMENSION(:,:),ALLOCATABLE :: Umean,Vmean,&
               ETAmean,Usum,Vsum,ETAsum, nu_smg, &
               UUsum,UUmean,UVsum,UVmean,VVsum,VVmean, &
               WWsum,WWmean,FRCXsum,FRCXmean,FRCYsum,FRCYmean,Wsurf, &
               DxSxx,DySxy,DySyy,DxSxy,PgrdX,PgrdY,DxUUH,DyUVH,DyVVH,DxUVH, &
               P_center,Q_center,U_davg,V_davg,U_davg_sum,V_davg_sum, &
               U_davg_mean,V_davg_mean,P_sum,Q_sum, &
               P_mean,Q_mean, &
               BreakDissX,BreakDissY,BreakDissX_sum,BreakDissY_sum
   REAL(SP)::T_INTV_mean = 20.0,T_sum=0.0,C_smg=0.25
   REAL(SP),DIMENSION(:,:),ALLOCATABLE :: &
             WaveHeightRMS,WaveHeightAve,Emax,Emin,& 
             HrmsSum,HavgSum
   INTEGER, DIMENSION(:,:),ALLOCATABLE :: Num_Zero_Up


! H = ETA + DEPTH
    REAL(SP),DIMENSION(:,:),ALLOCATABLE :: Depth,H,&
         DepthNode,Depthx,Depthy
    REAL(SP)::Depth_Flat, SLP,Xslp

! updating variables
    REAL(SP),DIMENSION(:,:),ALLOCATABLE::Ubar0,Vbar0,Eta0,&
                            Ubar,Vbar,Eta

! WATER LEVEL
    REAL(SP) :: WaterLevel = 0.0

! PARAMETERS TO OUTPUT
    INTEGER :: OUTPUT_RES = 1
    LOGICAL :: OUT_U=.FALSE., OUT_V=.FALSE.,OUT_ETA=.FALSE., &
               OUT_MASK=.FALSE.



! dispersion control
    LOGICAL :: OBSTACLE=.FALSE., HOT_START=.FALSE.
! sponge
    LOGICAL :: SPONGE_ON=.FALSE.
! breaking
    LOGICAL :: SHOW_BREAKING=.TRUE.

! slope control, use it in caution, should set false unless for a spherical
! ocean basin domain and slope > 1:5
    LOGICAL :: SLOPE_CTR = .FALSE.
    REAL(SP) :: MAX_SLOPE


!       LOGICAL :: WAVE_DATA = .FALSE.
    CHARACTER(LEN=80) WAVE_DATA_TYPE

! eddy viscosity breaking
     REAL(SP),DIMENSION(:,:),ALLOCATABLE :: AGE_BREAKING
     REAL(SP),DIMENSION(:,:),ALLOCATABLE :: ROLLER_FLUX,UNDERTOW_U,UNDERTOW_V
     REAL(SP),DIMENSION(:,:),ALLOCATABLE :: BreakSourceX,BreakSourceY
     REAL(SP),DIMENSION(:,:),ALLOCATABLE :: FrcInsX,FrcInsY
     REAL(SP) :: Cbrk1=0.65,Cbrk2=0.35,T_brk
     ! use T_brk to judge breakers
     LOGICAL :: INI_UVZ=.FALSE.
     LOGICAL :: BED_DEFORMATION = .FALSE.

    REAL(SP),DIMENSION(:,:),ALLOCATABLE :: nu_break,nu_sponge
    REAL(SP) :: nu_bkg
    LOGICAL :: DIFFUSION_SPONGE = .FALSE.
    LOGICAL :: DIRECT_SPONGE = .FALSE.
    LOGICAL :: FRICTION_SPONGE = .FALSE.
    LOGICAL :: VISCOSITY_BREAKING = .FALSE.
    REAL(SP),DIMENSION(:,:),ALLOCATABLE :: CD_4_SPONGE
    REAL(SP) :: Csp = 0.15
    REAL(SP) :: CDsponge = 0.0
    REAL(SP) :: WAVEMAKER_Cbrk
    LOGICAL :: ETA_LIMITER = .FALSE.
    REAL(SP) :: CrestLimit, TroughLimit

! ETA BLOW UP 
    REAL(SP) :: EtaBlowVal

! STEADY_TIME
    REAL(SP) :: STEADY_TIME

END MODULE GLOBAL
