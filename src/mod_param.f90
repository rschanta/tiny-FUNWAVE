MODULE PARAM  

! MPI SETUP
# if defined PARALLEL
        USE MPI
# endif
        IMPLICIT NONE    

! SINGLE PRECISION
    INTEGER, PARAMETER::SP=SELECTED_REAL_KIND(6,30)
! MPI Precision
# if defined PARALLEL
       INTEGER, PARAMETER::MPI_SP=MPI_REAL
# endif

! CONSTANTS
        REAL(SP), PARAMETER::pi=3.141592653
        REAL(SP), PARAMETER::R_earth = 6371000.0_SP
        REAL(SP), PARAMETER::SMALL=0.000001_SP
        REAL(SP), PARAMETER::LARGE=999999.0_SP
        REAL(SP), PARAMETER:: grav=9.81_SP
        REAL(SP), PARAMETER:: zero = 0.0_SP
        REAL(SP), PARAMETER:: RHO_AW = 0.0012041_SP  ! relative to water
        REAL(SP), PARAMETER:: RHO_AIR = 1.15_SP  ! absolute value
        REAL(SP), PARAMETER:: RHO_WATER = 1000.0_SP
        REAL(SP), PARAMETER:: DEG2RAD = 0.0175_SP

! GLOBAL VARS
        INTEGER :: I,J,K, RES
        INTEGER :: itmp1,itmp2,itmp3,itmp4,itmp5
        REAL(SP):: tmp1,tmp2,tmp3,tmp4,tmp5

END MODULE PARAM
    