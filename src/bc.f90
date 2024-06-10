SUBROUTINE BOUNDARY_CONDITION
    USE GLOBAL
    IMPLICIT NONE
    REAL(SP)::Xi,Deps


! DEAL WITH THE EDGES OF THE DOMAIN FIRSt

! WESTERN EDGE
#   if defined (PARALLEL)
        if ( n_west .eq. MPI_PROC_NULL ) then
#   endif
  IF(WaveMaker(1:3)=='ABS'.OR.WaveMaker(1:11)=='LEFT_BC_IRR')THEN
    ! do nothing
  ELSE
    DO J=Jbeg,Jend
     P(Ibeg,J)=ZERO
     Xi=EtaRxR(Ibeg,J)
     Deps=Depthx(Ibeg,J)
     Fx(Ibeg,J)=0.5_SP*GRAV*(Xi*Xi*Gamma3+2.0_SP*Xi*Deps)
     Gx(Ibeg,J)=ZERO
     ENDDO
  ENDIF ! left bc wavemaker


# if defined (PARALLEL)
     endif
# endif

! EASTERN EDGE
# if defined (PARALLEL)
       if ( n_east .eq. MPI_PROC_NULL ) then
# endif


    DO J=Jbeg,Jend
     P(Iend1,J)=ZERO
     Xi=EtaRxL(Iend1,J)
     Deps=Depthx(Iend1,J)
     Fx(Iend1,J)=0.5_SP*GRAV*(Xi*Xi*Gamma3+2.0_SP*Xi*Deps)
     Gx(Iend1,J)=ZERO
    ENDDO


# if defined (PARALLEL)
     endif
# endif

# if defined (CARTESIAN)
! y direction
  IF(PERIODIC)THEN
!   do nothing
  ELSE
# endif

! SOUTHERN EDGE
# if defined (PARALLEL)
     if ( n_suth .eq. MPI_PROC_NULL ) then
# endif

    DO I=Ibeg,Iend
     Q(I,Jbeg)=ZERO
     Fy(I,Jbeg)=ZERO
     Xi=EtaRyR(I,Jbeg)
     Deps=Depthy(I,Jbeg)
     Gy(I,Jbeg)=0.5_SP*GRAV*(Xi*Xi*Gamma3+2.0_SP*Xi*Deps)
     ENDDO

# if defined (PARALLEL)
     endif
# endif

! NORTHERN EDGE
# if defined (PARALLEL)
     if ( n_nrth .eq. MPI_PROC_NULL ) then
# endif

    DO I=Ibeg,Iend
     Q(I,Jend1)=ZERO
     Fy(I,Jend1)=ZERO
     Xi=EtaRyL(I,Jend1)
     Deps=Depthy(I,Jend1)
     Gy(I,Jend1)=0.5_SP*GRAV*(Xi*Xi*Gamma3+2.0_SP*Xi*Deps)
    ENDDO

# if defined (PARALLEL)
    endif
# endif

# if defined (CARTESIAN)
   ENDIF
# endif

! mask points
! Jeff pointed out the loop should be Jbeg-1, Jend+1
! The problem is that the fluxes on the inter-processor boundaries may be
!modified if the point next to the boundary (e.g., in the ghost cells,
!managed by a different processor) is land, but as is the routine doesnt
!check for this. 

    DO j=Jbeg-1,Jend+1
    DO i=Ibeg-1,Iend+1
     IF(MASK(I,J)<1)THEN
       P(I,J)=ZERO
! Jeff reported a bug here for parallel version
# if defined(PARALLEL)
       IF((I/=Ibeg).or.(n_west.ne.MPI_PROC_NULL))THEN
# else
       IF(I/=Ibeg)THEN
# endif
!         Fx(I,J)=0.5_SP*GRAV*HxL(I,J)*HxL(I,J)*MASK(I-1,J)
!new splitting method
     Xi=EtaRxL(I,J)
     Deps=Depthx(I,J)
        Fx(I,J)=0.5_SP*GRAV*(Xi*Xi*Gamma3+2.0_SP*Xi*Deps)*MASK(I-1,J)
       ELSE
        Fx(I,J)=ZERO
       ENDIF
       Gx(I,J)=ZERO

       P(I+1,J)=ZERO
! Jeff also here
# if defined(PARALLEL)
       IF((I/=Iend).or.(n_east.ne.MPI_PROC_NULL))THEN
# else
       IF(I/=Iend)THEN
# endif
!         Fx(I+1,J)=0.5_SP*GRAV*HxR(I+1,J)*HxR(I+1,J)*MASK(I+1,J)
! new splitting method
     Xi=EtaRxR(I+1,J)
     Deps=Depthx(I+1,J)
        Fx(I+1,J)=0.5_SP*GRAV*(Xi*Xi*Gamma3+2.0_SP*Xi*Deps)*MASK(I+1,J)
       ELSE
        Fx(I+1,J)=ZERO
       ENDIF
       Gx(I+1,J)=ZERO

       Q(I,J)=ZERO
       Fy(I,J)=ZERO
! Jeff also here
# if defined(PARALLEL)
       IF((J/=Jbeg).or.(n_suth.ne.MPI_PROC_NULL))THEN
# else
       IF(J/=Jbeg)THEN
# endif
!         Gy(I,J)=0.5_SP*GRAV*HyL(I,J)*HyL(I,J)*MASK(I,J-1)
! new splitting method
     Xi=EtaRyL(I,J)
     Deps=Depthy(I,J)
        Gy(I,J)=0.5_SP*GRAV*(Xi*Xi*Gamma3+2.0_SP*Xi*Deps)*MASK(I,J-1)
       ELSE
        Gy(I,J)=ZERO
       ENDIF

       Q(I,J+1)=ZERO
       Fy(I,J+1)=ZERO
! Jeff also here
# if defined(PARALLEL)
       IF((J/=Jend).or.(n_nrth.ne.MPI_PROC_NULL))THEN
# else
       IF(J/=Jend)THEN
# endif
!         Gy(I,J+1)=0.5_SP*GRAV*HyR(I,J+1)*HyR(I,J+1)*MASK(I,J+1)
! new splitting method
     Xi=EtaRyR(I,J+1)
     Deps=Depthy(I,J+1)
        Gy(I,J+1)=0.5_SP*GRAV*(Xi*Xi*Gamma3+2.0_SP*Xi*Deps)*MASK(I,J+1)
       ELSE
        Gy(I,J+1)=ZERO
       ENDIF
     ENDIF
    ENDDO
    ENDDO

END SUBROUTINE BOUNDARY_CONDITION


!-------------------------------------------------------------------
!
!   This subroutine is used to collect data into ghost cells                                                         
!
!   HISTORY:
!   07/09/2010 Fengyan Shi, use dummy variables 2) add vtype=3
!
!-------------------------------------------------------------------
SUBROUTINE EXCHANGE_DISPERSION
   USE GLOBAL
   IMPLICIT NONE
   INTEGER :: VTYPE

# if defined (CARTESIAN)
   VTYPE=2
   CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Uxx,VTYPE,PERIODIC)
   CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DUxx,VTYPE,PERIODIC)
   VTYPE=3
   CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Vyy,VTYPE,PERIODIC)
   CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DVyy,VTYPE,PERIODIC)

   VTYPE=1
   CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Uxy,VTYPE,PERIODIC)
   CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DUxy,VTYPE,PERIODIC)
   CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Vxy,VTYPE,PERIODIC)
   CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DVxy,VTYPE,PERIODIC)

   IF(Gamma2>ZERO)THEN

     IF(DISP_TIME_LEFT)THEN
       VTYPE=1 ! symetric in both direction
       CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,ETAT,VTYPE,PERIODIC)
       VTYPE=2  ! like u
       CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,ETATx,VTYPE,PERIODIC)
       VTYPE=3  ! like v
       CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,ETATy,VTYPE,PERIODIC) 
     ELSE
       VTYPE=2
       CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Ut,VTYPE,PERIODIC)
       VTYPE=3
       CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Vt,VTYPE,PERIODIC)

       VTYPE=1
       CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Utx,VTYPE,PERIODIC)
       VTYPE=1
       CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Vty,VTYPE,PERIODIC)

       VTYPE=2
       CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Utxx,VTYPE,PERIODIC)
       VTYPE=3
       CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Vtyy,VTYPE,PERIODIC)

       VTYPE=1
       CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Utxy,VTYPE,PERIODIC) 
       CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Vtxy,VTYPE,PERIODIC) 

       VTYPE=2
       CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DUtxx,VTYPE,PERIODIC)
       VTYPE=3
       CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DVtyy,VTYPE,PERIODIC)

       VTYPE=1
       CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DUtxy,VTYPE,PERIODIC) 
       CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DVtxy,VTYPE,PERIODIC) 
   
     ENDIF

     VTYPE=1  ! symetric in both direction
     CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Ux,VTYPE,PERIODIC)
     CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DUx,VTYPE,PERIODIC)
     CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Vy,VTYPE,PERIODIC)
     CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DVy,VTYPE,PERIODIC)
     VTYPE=3  !like v
     CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Uy,VTYPE,PERIODIC)
     CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DUy,VTYPE,PERIODIC)
     Vtype=2  !like u
     CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Vx,VTYPE,PERIODIC)
     CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DVx,VTYPE,PERIODIC)
     VTYPE=2  ! like u
     CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,ETAx,VTYPE,PERIODIC)
     VTYPE=3  ! like v
     CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,ETAy,VTYPE,PERIODIC)   

   ENDIF

# else
   VTYPE=2
   CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Uxx,VTYPE)
   CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DUxx,VTYPE)
   VTYPE=3
   CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Vyy,VTYPE)
   CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DVyy,VTYPE)

   VTYPE=1
   CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Uxy,VTYPE)
   CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DUxy,VTYPE)
   CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Vxy,VTYPE)
   CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DVxy,VTYPE)

     VTYPE=1  ! symetric in both direction
     CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Ux,VTYPE)
     CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DUx,VTYPE)
     CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Vy,VTYPE)
     CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DVy,VTYPE)
     VTYPE=3  !like v
     CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Uy,VTYPE)
     CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DUy,VTYPE)
     Vtype=2  !like u
     CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Vx,VTYPE)
     CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,DVx,VTYPE)
     VTYPE=2  ! like u
     CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,ETAx,VTYPE)
     VTYPE=3  ! like v
     CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,ETAy,VTYPE)   
# endif
   
END SUBROUTINE EXCHANGE_DISPERSION

!---------------------------------------------------------------------------------------
!
!   EXCHANGE subroutine is used to collect data into ghost cells                                                         
!
!   HISTORY:
!   07/09/2010 Fengyan Shi 
!     1) use dummy variables 2) add vtype=3
!   08/19/2015 Choi, corrected segmentation fault, maybe memory leaking                                       
!
!---------------------------------------------------------------------------------------
SUBROUTINE EXCHANGE
   USE GLOBAL
   IMPLICIT NONE
   INTEGER :: VTYPE
# if defined (CARTESIAN)
   REAL(SP),DIMENSION(Mloc,Nloc) :: rMASK
# endif

# if defined (CARTESIAN)
    ! MIRRORING ALL SIDES
    VTYPE=1
        CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Eta,VTYPE,PERIODIC)
        ! for radiation stress
        CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Wsurf,VTYPE,PERIODIC)
        CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,P_center,VTYPE,PERIODIC)
        CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,Q_center,VTYPE,PERIODIC)
        CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,U_davg,VTYPE,PERIODIC)
        CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,V_davg,VTYPE,PERIODIC)

        IF(VISCOSITY_BREAKING) THEN
        CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,AGE_BREAKING,VTYPE,PERIODIC)
        ENDIF
        
        !ykchoi (08.19.2015) :: new variable (WAVEMAKER_VIS) 	
        IF(VISCOSITY_BREAKING .OR. WAVEMAKER_VIS) THEN
        CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,nu_break,VTYPE,PERIODIC)
        ENDIF

        IF(VISCOSITY_BREAKING .OR. DIFFUSION_SPONGE .OR. WAVEMAKER_VIS)THEN
            CALL PHI_COLL_VARIABLE_LENGTH(Mloc1,Nloc,Ibeg,Iend1,Jbeg,Jend,Nghost,P,VTYPE)
            CALL PHI_COLL_VARIABLE_LENGTH(Mloc,Nloc1,Ibeg,Iend,Jbeg,Jend1,Nghost,Q,VTYPE)
        ENDIF

        rMASK = MASK ! for periodic boundary condition
        CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,rMASK,VTYPE,PERIODIC)  
        MASK = rMASK  
    
    ! NEGATIVE MIRRORING ON WEST/EAST, MIRRORING ON NORTH/SOUTH
    VTYPE=2
        CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,U,VTYPE,PERIODIC)
        CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,HU,VTYPE,PERIODIC)
    
    ! MIRRORING ON EAST/WEST, NEGATIVE MIRRORING ON NORTH/SOUTH
   VTYPE=3
   CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,V,VTYPE,PERIODIC)
   CALL PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,HV,VTYPE,PERIODIC)

# endif

! etaR x mask is a wrong idea
!    Eta=Eta*MASK

   U=U*MASK
   V=V*MASK
   HU=HU*MASK
   HV=HV*MASK

   
END SUBROUTINE EXCHANGE

!-----------------------------------------------------------------------------------
!
!   PHI_COLL_VARIABLE_LENGTH subroutine is used to collect data into ghost cells                                                         
!
!   HISTORY:
!   07/09/2010 Fengyan Shi                                      
!
!-----------------------------------------------------------------------------------
SUBROUTINE PHI_COLL_VARIABLE_LENGTH(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,PHI,VTYPE)

   USE PARAM
   USE GLOBAL,ONLY : WaveMaker

# if defined (PARALLEL)
   USE GLOBAL, ONLY : n_east,n_west,n_suth,n_nrth
# endif

   IMPLICIT NONE
   INTEGER,INTENT(IN) :: VTYPE
   INTEGER,INTENT(IN) :: Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost
   REAL(SP),INTENT(INOUT) :: PHI(Mloc,Nloc)

     ! x-direction
# if defined (PARALLEL)
   if ( n_west .eq. MPI_PROC_NULL ) then
# endif

   IF(WaveMaker(1:3)=='ABS'.OR.WaveMaker(1:11)=='LEFT_BC_IRR')THEN
     ! do nothing
   ELSE
     DO J=Jbeg,Jend  
     DO K=1,Nghost
       PHI(K,J)=PHI(Ibeg+Nghost-K,J)
     ENDDO
     ENDDO
   ENDIF

# if defined (PARALLEL)
   endif
# endif

# if defined (PARALLEL)
   if ( n_east .eq. MPI_PROC_NULL ) then
# endif

     DO J=Jbeg,Jend  
     DO K=1,Nghost
       PHI(Iend+K,J)=PHI(Iend-K+1,J)
     ENDDO
     ENDDO

# if defined (PARALLEL)
   endif
# endif


     ! y-direction and corners
# if defined (PARALLEL)
   if ( n_suth .eq. MPI_PROC_NULL ) then
# endif

     DO I=1,Mloc
     DO K=1,Nghost
       PHI(I,K)=PHI(I,Jbeg+Nghost-K)
     ENDDO
     ENDDO

# if defined (PARALLEL)
   endif
# endif

# if defined (PARALLEL)
   if ( n_nrth .eq. MPI_PROC_NULL ) then
# endif

     DO I=1,Mloc
     DO K=1,Nghost
       PHI(I,Jend+K)=PHI(I,Jend-K+1)
     ENDDO
     ENDDO

# if defined (PARALLEL)
   endif
# endif


END SUBROUTINE PHI_COLL_VARIABLE_LENGTH

!-------------------------------------------------------------------------------------
!
!   PHI_COLL COLLECT DATA INTO GHOST CELLS
!
!   HISTORY:
!     05/01/2010  Fengyan Shi
!     09/07/2010 Fengyan Shi, fix:
!       1) u v symmetric problem, 2) remove use global 3) fix bug
!     05/27/2010 Gangfeng Ma, corrected some bugs
!
!-------------------------------------------------------------------------------------



SUBROUTINE PHI_COLL(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost,PHI,VTYPE,PERIODIC)
    USE PARAM
    USE GLOBAL,ONLY : WaveMaker

! Set up parallel stuff
    USE GLOBAL, ONLY : n_east,n_west,n_suth,n_nrth, &
                    comm2d, ier,myid,PX,PY,&
                      NumberProcessor,ProcessorID

! Variable Declarations
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: VTYPE
    INTEGER,INTENT(IN) :: Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,Nghost
    REAL(SP),INTENT(INOUT) :: PHI(Mloc,Nloc)


! Works for Cartesian/Parallel only
    LOGICAL :: PERIODIC
    INTEGER,DIMENSION(1) :: req
    INTEGER :: nreq,len,ll,l,II,JJ
    integer,dimension(MPI_STATUS_SIZE,1) :: status
    REAL(SP),DIMENSION(Mloc,Nghost) :: xx,send2d
    REAL(sp) :: myvar,mybeta


!@! PERIODIC BOUNDARY CONDITIONS : WRAPS NORTH AND SOUTH
IF(PERIODIC)THEN
    IF (PY>1)THEN
        len=Mloc*Nghost

        ! LOOP ACROSS HORIZONTAL DOMAIN
        DO II = 1,PX
            ! IF WE ARE AT SOUTHERN EDGE
            if(myid==ProcessorID(II,1)) then
                DO I = 1,Mloc
                    DO J = 1,Nghost
                        send2d(I,J)=PHI(I,J+Nghost)
                    ENDDO
                ENDDO
                ! SEND TO NORTHERN EDGE
                call MPI_SEND(send2d,len,MPI_SP,ProcessorID(II,PY),0,MPI_COMM_WORLD,ier)
            endif 

            ! IF WE ARE AT NORTHERN EDGE
            if(myid==ProcessorID(II,PY))then
                ! RECEIVE DATA FROM SOUTHERN EDGE
                call MPI_IRECV(xx,len,MPI_SP,ProcessorID(II,1),0,MPI_COMM_WORLD,req(1),ier)
                ! WAIT UNTIL THE LAST PIECE OF DATA WAS RECEIVED
                call MPI_WAITALL(1,req,status,ier)
                    DO I=1,Mloc
                        DO J=1,Nghost
                            PHI(I,Jend+J)=xx(I,J)
                        ENDDO
                    ENDDO
            endif

            ! IF WE ARE AT NORTHERN EDGE
            if(myid==ProcessorID(II,PY)) then
                DO I = 1,Mloc
                    DO J = 1,Nghost
                        send2d(I,J)=PHI(I,Jend-Nghost+J)
                    ENDDO
                ENDDO
                ! send from master
                call MPI_SEND(send2d,len,MPI_SP,ProcessorID(II,1),0,MPI_COMM_WORLD,ier)
            endif 

            ! IF WE ARE AT SOUTHERN EDGE
            if(myid==ProcessorID(II,1))then
                ! RECEIVE DATA FROM NORTHERN EDGE
                call MPI_IRECV(xx,len,MPI_SP,ProcessorID(II,PY),0,MPI_COMM_WORLD,req(1),ier)
                ! WAIT UNTIL THE LAST PIECE OF DATA WAS RECEIVED
                call MPI_WAITALL(1,req,status,ier)
                DO I=1,Mloc
                    DO J=1,Nghost
                        PHI(I,J)=xx(I,J)
                    ENDDO
                ENDDO
            endif
        ENDDO  

    ! SIMPLER CASE OF Y = 1
    ELSE  
        DO I = Ibeg,Iend
            DO J = 1,Nghost
                PHI(I,J) = PHI(I,Jend-Nghost+J)
                PHI(I,Jend+J) = PHI(I,J+Nghost)
            ENDDO
        ENDDO
   ENDIF
ENDIF ! end periodic 

!@! VTYPE == 1 (eta)
IF(VTYPE==1) THEN  
    ! If at western edge
    if ( n_west .eq. MPI_PROC_NULL ) then
        ! Wavemaker
        IF(WaveMaker(1:3)=='ABS'.OR.WaveMaker(1:11)=='LEFT_BC_IRR')THEN
            ! do nothing
        ELSE
            DO J=Jbeg,Jend  
                DO K=1,Nghost
                    PHI(K,J)=PHI(Ibeg+Nghost-K,J)
                ENDDO
            ENDDO
        ENDIF  

    ! If at eastern edge
    if ( n_east .eq. MPI_PROC_NULL ) then
        DO J=Jbeg,Jend  
            DO K=1,Nghost
                PHI(Iend+K,J)=PHI(Iend-K+1,J)
            ENDDO
        ENDDO

    IF(.NOT.PERIODIC)THEN
        ! If at southern edge
        if ( n_suth .eq. MPI_PROC_NULL ) then
            DO I=1,Mloc
                DO K=1,Nghost
                    PHI(I,K)=PHI(I,Jbeg+Nghost-K)
                ENDDO
            ENDDO
        ENDIF

        ! If at northern edge
        if ( n_nrth .eq. MPI_PROC_NULL ) then
            DO I=1,Mloc
                DO K=1,Nghost
                    PHI(I,Jend+K)=PHI(I,Jend-K+1)
                ENDDO
            ENDDO
        ENDIF
ENDIF ! end vtype=1

!@! VTYPE == 2 (u)
IF(VTYPE==2) THEN  ! for u (x-mirror condition)
    ! If at western edge
    if ( n_west .eq. MPI_PROC_NULL ) then
        IF(WaveMaker(1:3)=='ABS'.OR.WaveMaker(1:11)=='LEFT_BC_IRR')THEN
            ! do nothing
        ELSE
            DO J=Jbeg,Jend
                DO K=1,Nghost
                    PHI(K,J)=-PHI(Ibeg+Nghost-K,J)
                ENDDO
            ENDDO
        ENDIF 
   endif

    ! If at eastern edge
    if ( n_east .eq. MPI_PROC_NULL ) then
        DO J=Jbeg,Jend
            DO K=1,Nghost
                PHI(Iend+K,J)=-PHI(Iend-K+1,J)
            ENDDO
        ENDDO
    endif



    ! y-direction and corners
    IF(.NOT.PERIODIC)THEN
        ! If at southern edge
        if ( n_suth .eq. MPI_PROC_NULL ) then
            DO I=1,Mloc
                DO K=1,Nghost
                    PHI(I,K)=PHI(I,Jbeg+Nghost-K)
                ENDDO
            ENDDO
        endif

        ! If at northern edge
        if ( n_nrth .eq. MPI_PROC_NULL ) then
            DO I=1,Mloc
                DO K=1,Nghost
                    PHI(I,Jend+K)=PHI(I,Jend-K+1)
                ENDDO
            ENDDO
        endif
    ENDIF
ENDIF ! end vtype=2

!@! VTYPE == 3 (v) (y-mirror)
IF(VTYPE==3) THEN ! for v (y-mirror condition)
    ! If at western edge
    if ( n_west .eq. MPI_PROC_NULL ) then
        IF(WaveMaker(1:3)=='ABS'.OR.WaveMaker(1:11)=='LEFT_BC_IRR')THEN
            ! do nothing
        ELSE
            DO J=Jbeg,Jend
                DO K=1,Nghost
                    PHI(K,J)=PHI(Ibeg+Nghost-K,J)
                ENDDO
            ENDDO
        ENDIF ! end left_bc wavemaker
   endif


   ! If at eastern edge
    if ( n_east .eq. MPI_PROC_NULL ) then
        DO J=Jbeg,Jend
            DO K=1,Nghost
                PHI(Iend+K,J)=PHI(Iend-K+1,J)
            ENDDO
        ENDDO
    endif

    IF(.NOT.PERIODIC)THEN
        ! If at southern edge
        if ( n_suth .eq. MPI_PROC_NULL ) then
            DO I=1,Mloc
                DO K=1,Nghost
                    PHI(I,K)=-PHI(I,Jbeg+Nghost-K)
                ENDDO
            ENDDO
        endif

        ! If at northern edge
        if ( n_nrth .eq. MPI_PROC_NULL ) then
            DO I=1,Mloc
                DO K=1,Nghost
                    PHI(I,Jend+K)=-PHI(I,Jend-K+1)
                ENDDO
            ENDDO
        endif
    ENDIF
ENDIF 

!@! VTYPE == 4 (cross derivatives?) uv cross-mirror ! VTYPE==4 for u and v cross-mirror
IF(VTYPE==4) THEN 
    ! If at western edge
    if ( n_west .eq. MPI_PROC_NULL ) then
        DO J=Jbeg,Jend
            DO K=1,Nghost
                PHI(K,J)=ZERO
            ENDDO
        ENDDO
    endif


    ! If at eastern edge
    if ( n_east .eq. MPI_PROC_NULL ) then
        DO J=Jbeg,Jend
            DO K=1,Nghost
                PHI(Iend+K,J)=ZERO
            ENDDO
        ENDDO
    endif

    ! If at southern edge (not exact)
    if ( n_suth .eq. MPI_PROC_NULL ) then
        DO I=1,Mloc
            DO K=1,Nghost
                PHI(I,K)=ZERO
            ENDDO
        ENDDO
    endif

    ! If at northern edge
    if ( n_nrth .eq. MPI_PROC_NULL ) then
        DO I=1,Mloc
            DO K=1,Nghost
                PHI(I,Jend+K)=ZERO
            ENDDO
        ENDDO
    endif
ENDIF ! end vtype=4

!@! VTYPE == 5 (symmetric?)
IF(VTYPE==5)THEN
    ! If at western edge
    if ( n_west .eq. MPI_PROC_NULL ) then
        DO J=Jbeg,Jend
            DO K=1,Nghost
                PHI(K,J)=PHI(Ibeg+Nghost-K,J)
            ENDDO
        ENDDO
    endif

    ! If at eastern edge
    if ( n_east .eq. MPI_PROC_NULL ) then
        DO J=Jbeg,Jend
            DO K=1,Nghost
                PHI(Iend+K,J)=PHI(Iend-K+1,J)
            ENDDO
        ENDDO
    endif

    ! If at southern edge
    if ( n_suth .eq. MPI_PROC_NULL ) then
        DO I=1,Mloc
            DO K=1,Nghost
                PHI(I,K)=PHI(I,Jbeg+Nghost-K)
            ENDDO
        ENDDO
    endif

    ! If at northern edge
    if ( n_nrth .eq. MPI_PROC_NULL ) then
        DO I=1,Mloc
            DO K=1,Nghost
                PHI(I,Jend+K)=PHI(I,Jend-K+1)
            ENDDO
        ENDDO
    endif
ENDIF 

!@! VTYPE == 6 (anti-symmetric?)
IF(VTYPE==6)THEN

    ! If at western edge
    if ( n_west .eq. MPI_PROC_NULL ) then
        DO J=Jbeg,Jend
            DO K=1,Nghost
                PHI(K,J)=-PHI(Ibeg+Nghost-K,J)
            ENDDO
        ENDDO 
    endif

    ! If at eastern edge
    if ( n_east .eq. MPI_PROC_NULL ) then
        DO J=Jbeg,Jend
            DO K=1,Nghost
                PHI(Iend+K,J)=-PHI(Iend-K+1,J)
            ENDDO
        ENDDO 
    endif

    ! If at southern edge
    if ( n_suth .eq. MPI_PROC_NULL ) then
        DO I=1,Mloc
            DO K=1,Nghost
                PHI(I,K)=-PHI(I,Jbeg+Nghost-K)
            ENDDO
        ENDDO   
    endif

    ! If at northern edge
    if ( n_nrth .eq. MPI_PROC_NULL ) then
        DO I=1,Mloc
            DO K=1,Nghost
                PHI(I,Jend+K)=-PHI(I,Jend-K+1)
            ENDDO
        ENDDO     
    endif
ENDIF ! end vtype=6

   call phi_exch (PHI)


END SUBROUTINE PHI_COLL



