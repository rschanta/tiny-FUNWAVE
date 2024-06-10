SUBROUTINE FLUXES
    USE GLOBAL
    IMPLICIT NONE

    ! ONLY APPLICABLE TO 4th ORDER
    IF(HIGH_ORDER(1:3)=='FOU') THEN   
        CALL CONSTRUCTION_HO   ! note choi used the name construction_ho_vanleer
        CALL WAVE_SPEED(Mloc,Nloc,Mloc1,Nloc1,UxL,UxR,VyL,VyR,HxL,HxR,HyL,HyR, &
            SxL,SxR,SyL,SyR)                       
    ENDIF

    ! ONLY USE HLL CASE
    IF(CONSTR(1:3)=='HLL')THEN
        CALL FLUX_AT_INTERFACE_HLL
    ENDIF

    ! UPDATE BOUNDARY CONDITION
    CALL BOUNDARY_CONDITION

END SUBROUTINE FLUXES

!-------------------------------------------------------------------------------------
!
!    FLUX_AT_INTERFACE_HLL is subroutine to do HLL scheme
!
!    HISTORY: 
!    05/06/2010 Fengyan Shi
!
!-------------------------------------------------------------------------------------
SUBROUTINE FLUX_AT_INTERFACE_HLL
    USE GLOBAL
    REAL(SP)::SR,SL,FL,FR,UL,UR

    CALL HLL(Mloc1,Nloc,SxL,SxR,PL,PR,EtaRxL,EtaRxR,P)
    CALL HLL(Mloc,Nloc1,SyL,SyR,QL,QR,EtaRyL,EtaRyR,Q)
    CALL HLL(Mloc1,Nloc,SxL,SxR,FxL,FxR,HUxL,HUxR,Fx)
    CALL HLL(Mloc,Nloc1,SyL,SyR,FyL,FyR,HUyL,HUyR,Fy)
    CALL HLL(Mloc1,Nloc,SxL,SxR,GxL,GxR,HVxL,HVxR,Gx)
    CALL HLL(Mloc,Nloc1,SyL,SyR,GyL,GyR,HVyL,HVyR,Gy)

END SUBROUTINE FLUX_AT_INTERFACE_HLL

!-------------------------------------------------------------------------------------
!
!    HLL is subroutine for HLL scheme
!
!    HISTORY: 
!      05/06/2010 Fengyan Shi
!
!-------------------------------------------------------------------------------------
SUBROUTINE HLL(M,N,SL,SR,FL,FR,UL,UR,FOUT)
    USE PARAM
    ! Dimension of Input and Output Arrays
    INTEGER,INTENT(IN)::M,N
    ! SL: Left wave speeds
    ! SR: Right wave speeds
    ! FL: Flux values on left
    ! FR: Flux values on right (also called !)
    ! UL: state values on left
    ! UR state values on right 
    ! FOUT: Output flux
    REAL(SP),INTENT(IN),DIMENSION(M,N)::SL,SR,FL,FR,UL,UR
    REAL(SP),INTENT(OUT),DIMENSION(M,N)::FOUT

    ! Loop over entire domain
    DO J=1,N
        DO I=1,M     
            ! Consider left wave speed > 0 (info propagating to right: what comes left must go right)
            IF(SL(I,J)>=ZERO) THEN
                FOUT(I,J)=FL(I,J)
            ! Consider right wave speed < 0 (info propagating to left: what comes right must go left)
            ELSEIF(SR(I,J)<=ZERO) THEN
                FOUT(I,J)=FR(I,J)
            ELSE
                ! Left wave < 0 and Right wave > 0
                ! Weighted combination of left and right by difference in state
                FOUT(I,J)=SR(I,J)*FL(I,J)-SL(I,J)*FR(I,J)+SL(I,J)*SR(I,J)*(UR(I,J)-UL(I,J))

                ! Prevent division by 0
                IF((ABS(SR(I,J)-SL(I,J)))<SMALL)THEN
                    FOUT(I,J)=FOUT(I,J)/SMALL
                ELSE
                    FOUT(I,J)=FOUT(I,J)/(SR(I,J)-SL(I,J))
                ENDIF
            ENDIF
        ENDDO
    ENDDO

END SUBROUTINE HLL



! ---------------------------------------------------
!
!    WAVE_SPEED is subroutine to calculate wave speed
!    no shear wave is calculated yet
!
!    HISTORY: 
!        01/21/2012 Fengyan Shi
!        10/29/2012 Fengyan Shi
!             Steve Brandt mentioned HxL and HxR not calculated 
!             inside ghost cells, it is not an issue though corrected
!             There's a bug found in Dmitry's case. In x direction, 
!             should use M1 and N1 for y direction.
!
! --------------------------------------------------
SUBROUTINE WAVE_SPEED(M,N,M1,N1,UL,UR,VL,VR,HxL,HxR,HyL,HyR,&
     SxL,SxR,SyL,SyR)
    USE PARAM
    USE GLOBAL, ONLY : Nghost
    IMPLICIT NONE
    INTEGER,INTENT(IN)::M,N,M1,N1
    REAL(SP),INTENT(IN),DIMENSION(M1,N)::UL,UR,HxL,HxR
    REAL(SP),INTENT(IN),DIMENSION(M,N1)::VL,VR,HyL,HyR
    REAL(SP),INTENT(OUT),DIMENSION(M1,N)::SxL,SxR
    REAL(SP),INTENT(OUT),DIMENSION(M,N1)::SyL,SyR         
    REAL(SP)::SQR_PHI_L,SQR_PHI_R,SQR_PHI_S,U_S


! X INTERFACE
    DO J=1+Nghost,N-Nghost
        DO I=1+Nghost,M1-Nghost
            ! Equation (40) components
                SQR_PHI_L=SQRT(GRAV*ABS(HxL(I,J)))
                SQR_PHI_R=SQRT(GRAV*ABS(HxR(I,J)))
                SQR_PHI_S=0.5*(SQR_PHI_L+SQR_PHI_R)+0.25*(UL(I,J)-UR(I,J))  
            ! Equation (39)
                U_S=0.5*(UL(I,J)+UR(I,J))+SQR_PHI_L-SQR_PHI_R
            ! Equations (37) and (38)
                SxL(I,J)=MIN(UL(I,J)-SQR_PHI_L,U_S-SQR_PHI_S)
                SxR(I,J)=MAX(UR(I,J)+SQR_PHI_R,U_S+SQR_PHI_S)
        ENDDO
    ENDDO

! UPDATE GHOST CELLS WITH CORRECT SPEEDS
    DO J=1+Nghost,N-Nghost
        DO I=1,Nghost
            SxL(I,J)=SxL(Nghost+1,J)
            SxR(I,J)=SxR(Nghost+1,J)
        ENDDO
        DO I=M1-Nghost+1,M1
            SxL(I,J)=SxL(M1-Nghost,J)
            SxR(I,J)=SxR(M1-Nghost,J)       
        ENDDO
    ENDDO

    DO I=1,M1
        DO J=1,Nghost
            SxL(I,J)=SxL(I,Nghost+1)
            SxR(I,J)=SxR(I,Nghost+1)
        ENDDO
        DO J=N-Nghost+1,N
            SxL(I,J)=SxL(I,N-Nghost)
            SxR(I,J)=SxR(I,N-Nghost)
        ENDDO
    ENDDO

! Y INTERFACE
    DO J=1+Nghost,N1-Nghost
        DO I=1+Nghost,M-Nghost
            ! Equation (40) components
                SQR_PHI_L=SQRT(GRAV*ABS(HyL(I,J)))
                SQR_PHI_R=SQRT(GRAV*ABS(HyR(I,J)))
                SQR_PHI_S=0.5*(SQR_PHI_L+SQR_PHI_R)+0.25*(VL(I,J)-VR(I,J)) 
            ! Equation (39)
                U_S=0.5*(VL(I,J)+VR(I,J))+SQR_PHI_L-SQR_PHI_R
            ! Equations (37) and (38)
                SyL(I,J)=MIN(VL(I,J)-SQR_PHI_L,U_S-SQR_PHI_S)
                SyR(I,J)=MAX(VR(I,J)+SQR_PHI_R,U_S+SQR_PHI_S)
        ENDDO
    ENDDO


! UPDATE GHOST CELLS WITH CORRECT SPEEDS
    DO I=1+Nghost,M-Nghost
        DO J=1,Nghost
            SyL(I,J)=SyL(I,Nghost+1)
            SyR(I,J)=SyR(I,Nghost+1)
        ENDDO
        DO J=N1-Nghost+1,N1
            SyL(I,J)=SyL(I,N1-Nghost)
            SyR(I,J)=SyR(I,N1-Nghost)       
        ENDDO
    ENDDO

    DO J=1,N1
        DO I=1,Nghost
            SyL(I,J)=SyL(Nghost+1,J)
            SyR(I,J)=SyR(Nghost+1,J)
        ENDDO
        DO I=M-Nghost+1,M
            SyL(I,J)=SyL(M-Nghost,J)
            SyR(I,J)=SyR(M-Nghost,J)
        ENDDO
    ENDDO

END SUBROUTINE WAVE_SPEED

! --------------- ykchoi (08/28/2016)
! This is subroutine for the fourth order MUSCL-TVD scheme.
! Van-Leer limiter is used for the third-order part, 
! and the minmod limiter is used for the fourth order part in the fourth order scheme.
!
! Reference
! Hybrid finite-volume finite-difference scheme for the solution of Boussinesq equations
! Erduran et al. (2005)
!   08/31/21016  fyshi changed the subroutine name HO_vanleer to HO
!                to make a consistency
! ---------------------------------------------------------------------------------------
SUBROUTINE CONSTRUCTION_HO
    USE GLOBAL
    IMPLICIT NONE

    REAL(SP),ALLOCATABLE :: MASK9u(:,:), MASK9v(:,:)
    ALLOCATE( MASK9u(1:Mloc1,1:Nloc), MASK9v(1:Mloc,1:Nloc1) )
   
    MASK9u(1:Mloc,1:Nloc)=MASK9(1:Mloc,1:Nloc)
    MASK9v(1:Mloc,1:Nloc)=MASK9(1:Mloc,1:Nloc)
    
    MASK9u(Mloc1,1:Nloc)=MASK9(Mloc,1:Nloc)
    MASK9v(1:Mloc,Nloc1)=MASK9(1:Mloc,Nloc)

! construct in x-direction
    CALL CONSTRUCT_HO_X(Mloc,Nloc,Mloc1,Ibeg,Iend,Jbeg,Jend,MASK,U,UxL,UxR)
    CALL CONSTRUCT_HO_X(Mloc,Nloc,Mloc1,Ibeg,Iend,Jbeg,Jend,MASK,V,VxL,VxR)
    CALL CONSTRUCT_HO_X(Mloc,Nloc,Mloc1,Ibeg,Iend,Jbeg,Jend,MASK,HU,HUxL,HUxR)
    CALL CONSTRUCT_HO_X(Mloc,Nloc,Mloc1,Ibeg,Iend,Jbeg,Jend,MASK,HV,HVxL,HVxR)
    CALL CONSTRUCT_HO_X(Mloc,Nloc,Mloc1,Ibeg,Iend,Jbeg,Jend,MASK,Eta,EtaRxL,EtaRxR)


    HxL=EtaRxL+Depthx
    HxR=EtaRxR+Depthx

    ! DEAL WITH DISPERSION
    IF(DISPERSION)THEN
        CALL CONSTRUCT_HO_X(Mloc,Nloc,Mloc1,Ibeg,Iend,Jbeg,Jend,MASK,U4,U4xL,U4xR)
        CALL CONSTRUCT_HO_X(Mloc,Nloc,Mloc1,Ibeg,Iend,Jbeg,Jend,MASK,V4,V4xL,V4xR)  

    ENDIF

    ! MOMENTUM
        
        PL(1:Mloc1,1:Nloc)=HUxL(1:Mloc1,1:Nloc) &
        + Gamma1*MASK9u(1:Mloc1,1:Nloc)*HxL(1:Mloc1,1:Nloc)*U4xL(1:Mloc1,1:Nloc)	

        PR(1:Mloc1,1:Nloc)=HUxR(1:Mloc1,1:Nloc) &
        + Gamma1*MASK9u(1:Mloc1,1:Nloc)*HxR(1:Mloc1,1:Nloc)*U4xR(1:Mloc1,1:Nloc)

        FxL(1:Mloc1,1:Nloc)=Gamma3*PL(1:Mloc1,1:Nloc)*(UxL(1:Mloc1,1:Nloc)  &
        +Gamma1*MASK9u(1:Mloc1,1:Nloc)*U4xL(1:Mloc1,1:Nloc)) &
        +0.5*GRAV*((EtaRxL(1:Mloc1,1:Nloc))*(EtaRxL(1:Mloc1,1:Nloc))*Gamma3  &
        +2.0_SP*(EtaRxL(1:Mloc1,1:Nloc))*(Depthx(1:Mloc1,1:Nloc)))

        FxR(1:Mloc1,1:Nloc)=Gamma3*PR(1:Mloc1,1:Nloc)*(UxR(1:Mloc1,1:Nloc)  &
        +Gamma1*MASK9u(1:Mloc1,1:Nloc)*U4xR(1:Mloc1,1:Nloc)) &
        +0.5*GRAV*((EtaRxR(1:Mloc1,1:Nloc))*(EtaRxR(1:Mloc1,1:Nloc))*Gamma3  &
        +2.0_SP*(EtaRxR(1:Mloc1,1:Nloc))*(Depthx(1:Mloc1,1:Nloc)))


    ! CROSS TERMS OF MOMENTUM
        GxL(1:Mloc1,1:Nloc)=HxL(1:Mloc1,1:Nloc)*    &
            ( UxL(1:Mloc1,1:Nloc) + Gamma1*MASK9u(1:Mloc1,1:Nloc)*U4xL(1:Mloc1,1:Nloc) )*   &
            ( VxL(1:Mloc1,1:Nloc) + Gamma1*MASK9u(1:Mloc1,1:Nloc)*V4xL(1:Mloc1,1:Nloc) )*   &
            Gamma3

        GxR(1:Mloc1,1:Nloc)=HxR(1:Mloc1,1:Nloc)*    &
            ( UxR(1:Mloc1,1:Nloc) + Gamma1*MASK9u(1:Mloc1,1:Nloc)*U4xR(1:Mloc1,1:Nloc) )*   &
            ( VxR(1:Mloc1,1:Nloc) + Gamma1*MASK9u(1:Mloc1,1:Nloc)*V4xR(1:Mloc1,1:Nloc) )*   &
            Gamma3


     
! construct in y-direction
    CALL CONSTRUCT_HO_Y(Mloc,Nloc,Nloc1,Ibeg,Iend,Jbeg,Jend,MASK,U,UyL,UyR)
    CALL CONSTRUCT_HO_Y(Mloc,Nloc,Nloc1,Ibeg,Iend,Jbeg,Jend,MASK,V,VyL,VyR)
    CALL CONSTRUCT_HO_Y(Mloc,Nloc,Nloc1,Ibeg,Iend,Jbeg,Jend,MASK,HV,HVyL,HVyR)
    CALL CONSTRUCT_HO_Y(Mloc,Nloc,Nloc1,Ibeg,Iend,Jbeg,Jend,MASK,HU,HUyL,HUyR)
# if defined (UseEtaScreen)
    CALL CONSTRUCT_HO_Y(Mloc,Nloc,Nloc1,Ibeg,Iend,Jbeg,Jend,MASK,EtaScreen,EtaRyL,EtaRyR)
# else
    CALL CONSTRUCT_HO_Y(Mloc,Nloc,Nloc1,Ibeg,Iend,Jbeg,Jend,MASK,Eta,EtaRyL,EtaRyR)
# endif

    HyL=EtaRyL+Depthy
    HyR=EtaRyR+Depthy

    IF(DISPERSION)THEN

      CALL CONSTRUCT_HO_Y(Mloc,Nloc,Nloc1,Ibeg,Iend,Jbeg,Jend,MASK,V4,V4yL,V4yR)
    CALL CONSTRUCT_HO_Y(Mloc,Nloc,Nloc1,Ibeg,Iend,Jbeg,Jend,MASK,U4,U4yL,U4yR)

    ENDIF

    QL(1:Mloc,1:Nloc1)=HVyL(1:Mloc,1:Nloc1)   &
      + Gamma1*MASK9v(1:Mloc,1:Nloc1)*HyL(1:Mloc,1:Nloc1)*V4yL(1:Mloc,1:Nloc1)

    QR(1:Mloc,1:Nloc1)=HVyR(1:Mloc,1:Nloc1)    &
      + Gamma1*MASK9v(1:Mloc,1:Nloc1)*HyR(1:Mloc,1:Nloc1)*V4yR(1:Mloc,1:Nloc1)

    GyL(1:Mloc,1:Nloc1)=Gamma3*QL(1:Mloc,1:Nloc1)*(VyL(1:Mloc,1:Nloc1)   &
       +Gamma1*MASK9v(1:Mloc,1:Nloc1)*V4yL(1:Mloc,1:Nloc1)) &
       +0.5*GRAV*((EtaRyL(1:Mloc,1:Nloc1))*(EtaRyL(1:Mloc,1:Nloc1))*Gamma3   &
       +2.0_SP*(EtaRyL(1:Mloc,1:Nloc1))*(Depthy(1:Mloc,1:Nloc1)))

    GyR(1:Mloc,1:Nloc1)=Gamma3*QR(1:Mloc,1:Nloc1)*(VyR(1:Mloc,1:Nloc1)    &
       +Gamma1*MASK9v(1:Mloc,1:Nloc1)*V4yR(1:Mloc,1:Nloc1)) &
       +0.5*GRAV*((EtaRyR(1:Mloc,1:Nloc1))*(EtaRyR(1:Mloc,1:Nloc1))*Gamma3    &
       +2.0_SP*(EtaRyR(1:Mloc,1:Nloc1))*(Depthy(1:Mloc,1:Nloc1)))





   FyL(1:Mloc,1:Nloc1)=HyL(1:Mloc,1:Nloc1)*   &
    ( UyL(1:Mloc,1:Nloc1) + Gamma1*MASK9v(1:Mloc,1:Nloc1)*U4yL(1:Mloc,1:Nloc1) )*   &
    ( VyL(1:Mloc,1:Nloc1) + Gamma1*MASK9v(1:Mloc,1:Nloc1)*V4yL(1:Mloc,1:Nloc1) )*   &
    Gamma3
      
   FyR(1:Mloc,1:Nloc1)=HyR(1:Mloc,1:Nloc1)*   &
    ( UyR(1:Mloc,1:Nloc1) + Gamma1*MASK9v(1:Mloc,1:Nloc1)*U4yR(1:Mloc,1:Nloc1) )*   &
    ( VyR(1:Mloc,1:Nloc1) + Gamma1*MASK9v(1:Mloc,1:Nloc1)*V4yR(1:Mloc,1:Nloc1) )*   &
    Gamma3


   deallocate( MASK9u, MASK9v )
    
END SUBROUTINE CONSTRUCTION_HO

! --------------- ykchoi (08/28/2016)
! Reference
! Hybrid finite-volume finite-difference scheme for the solution of Boussinesq equations
! Erduran et al. (2005)
! ---------------------------------------------------------------------------------------
SUBROUTINE CONSTRUCT_HO_X(M,N,M1,Ibeg,Iend,Jbeg,Jend,MASK,Vin,OutL,OutR)
    USE PARAM
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: M,N,M1,Ibeg,Iend,Jbeg,Jend

    REAL(SP),INTENT(IN),DIMENSION(M,N) :: Vin
    INTEGER,INTENT(IN),DIMENSION(M,N) :: MASK
    ! Output on left and right
        REAL(SP),INTENT(OUT),DIMENSION(M1,N) :: OutL,OutR
    
    
        REAL(SP),DIMENSION(M,N) :: Din
        ! Gradients of input field and minmods
        REAL(SP) :: TXP1,TXP2,TXP3,DVP1,DVP2,DVP3

        REAL(SP) :: VAN1, VAN2, RAT

    Din=0.0_SP
    DO J=Jbeg,Jend
       
     ! Fourth order MUSCL TVD with van Leer and minmod limiter
     ! The minmod limiter is used for the fourth order part in the fourth order scheme
    DO I=Ibeg-1,Iend+2
        ! Delta terms that appear in Equation 31
        TXP1=Vin(I-1,J)-Vin(I-2,J)
        TXP2=Vin(I,J)-Vin(I-1,J)
        TXP3=Vin(I+1,J)-Vin(I,J)

        ! minmod(-0.5,0.5,1.5)
            if (TXP1.ge.0.0_SP) then
                DVP1=MAX(0.0_SP,MIN(TXP1,2.0*TXP2,2.0*TXP3))          
            else
                DVP1=MIN(0.0_SP,MAX(TXP1,2.0*TXP2,2.0*TXP3))          
            endif
        ! minmod(0.5,1.5,-0.5)
            if (TXP2.ge.0.0_SP) then
                DVP2=MAX(0.0_SP,MIN(TXP2,2.0*TXP3,2.0*TXP1))          
            else
                DVP2=MIN(0.0_SP,MAX(TXP2,2.0*TXP3,2.0*TXP1))          
            endif
        ! minmod(1.5,-0.5,0.5)
            if (TXP3.ge.0.0_SP) then
                DVP3=MAX(0.0_SP,MIN(TXP3,2.0*TXP1,2.0*TXP2))          
            else
                DVP3=MIN(0.0_SP,MAX(TXP3,2.0*TXP1,2.0*TXP2))          
            endif

            ! DEAL WITH DRY CELLS (gradient in direction of wet region)
            IF(MASK(I-2,J)==0.OR.MASK(I+1,J)==0)THEN
                TXP2=Vin(I,J)-Vin(I-1,J)
                TXP1=TXP2
                TXP3=TXP2

                if (TXP1.ge.0.0_SP) then
                    DVP1=MAX(0.0_SP,MIN(TXP1,2.0*TXP2,2.0*TXP3))
                    else
                    DVP1=MIN(0.0_SP,MAX(TXP1,2.0*TXP2,2.0*TXP3))
                endif
                if (TXP2.ge.0.0_SP) then
                    DVP2=MAX(0.0_SP,MIN(TXP2,2.0*TXP3,2.0*TXP1))
                    else
                    DVP2=MIN(0.0_SP,MAX(TXP2,2.0*TXP3,2.0*TXP1))
                endif
                if (TXP3.ge.0.0_SP) then
                    DVP3=MAX(0.0_SP,MIN(TXP3,2.0*TXP1,2.0*TXP2))          
                    else
                    DVP3=MIN(0.0_SP,MAX(TXP3,2.0*TXP1,2.0*TXP2))          
                endif
            ! here actually DVP1=DVP2=DVP3=TXP2
            ENDIF

            ! ZERO GRADIENT CONDITION
            ! dry I-2 D-1 D I+1, zero gradient
            IF(MASK(I-1,J)==0.OR.MASK(I,J)==0)THEN
                DVP1=ZERO
                DVP2=ZERO
                DVP3=ZERO
            ENDIF

            ! Not sure
            Din(I,J)=TXP2-1.0_SP/6.0_SP*(DVP3-2.0_SP*DVP2+DVP1)
    ENDDO

    ! Fourth order MUSCL TVD with van Leer and minmod limiter
    ! The van-Leer limiter is used for the third order part in the fourth order scheme
    DO I=Ibeg,Iend+1
        ! OUT ON THE LEFT
            ! Setting up for Equations 29 and 30
                TMP1 = Din(I-1,J);   TMP2 = Din(I,J);
            ! Deal specially with the case of small gradients
                IF( ABS(TMP1).le.SMALL ) TMP1 = SMALL*SIGN( 1.0_SP, TMP1 )
                IF( ABS(TMP2).le.SMALL ) TMP2 = SMALL*SIGN( 1.0_SP, TMP2 )

            ! Equation 33: upwind ratio
                RAT = TMP2/TMP1
            ! Van Leer from Eqn (33)
                VAN1 = 0.0_SP;
                IF(abs(1.0 + RAT) > SMALL ) VAN1 = ( RAT + ABS(RAT) )/( 1.0_SP + RAT)

            ! Equation 33: upwind ratio SWITCH DIRECTION
                RAT = TMP1/TMP2
            ! Van Leer from Eqn (33)
                VAN2 = 0.0_SP;
                IF( abs(1.0 + RAT) > SMALL ) VAN2 = ( RAT + ABS(RAT) )/( 1.0_SP + RAT)

            ! EQUATION 29
            OutL(I,J) = Vin(I-1,J) + ( 1.0_SP/6.0_SP )*( VAN1*TMP1 + 2.0_SP*VAN2*TMP2 )

        ! OUT ON THE RIGHT
            ! Setting up for Equations 29 and 30
                TMP1 = Din(I,J);   TMP2 = Din(I+1,J);
            ! Deal specially with the case of small gradients
                IF( ABS(TMP1).le.SMALL ) TMP1 = SMALL*SIGN( 1.0_SP, TMP1 )
                IF( ABS(TMP2).le.SMALL ) TMP2 = SMALL*SIGN( 1.0_SP, TMP2 )
            ! Equation 33: upwind ratio
                RAT = TMP2/TMP1
            ! Van Leer from Eqn (33)
                VAN1 = 0.0_SP;
                IF( abs(1.0 + RAT) > SMALL ) VAN1 = ( RAT + ABS(RAT) )/( 1.0_SP + RAT)
            ! Equation 33: upwind ratio SWITCH DIRECTION
                RAT = TMP1/TMP2
            ! Van Leer from Eqn (33)
                VAN2 = 0.0_SP;
                IF( abs(1.0 + RAT) > SMALL ) VAN2 = ( RAT + ABS(RAT) )/( 1.0_SP + RAT)
            ! EQUATION 30
            OutR(I,J) = Vin(I,J) - ( 1.0_SP/6.0_SP )*( 2.0_SP*VAN1*TMP1 + VAN2*TMP2 )
        ENDDO

    ENDDO 

END SUBROUTINE CONSTRUCT_HO_X

! --------------- ykchoi (08/28/2016)
! Reference
! Hybrid finite-volume finite-difference scheme for the solution of Boussinesq equations
! Erduran et al. (2005)
! ---------------------------------------------------------------------------------------
SUBROUTINE CONSTRUCT_HO_Y(M,N,N1,Ibeg,Iend,Jbeg,Jend,MASK,Vin,OutL,OutR)
    USE PARAM
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: M,N,N1,Ibeg,Iend,Jbeg,Jend

    REAL(SP),INTENT(IN),DIMENSION(M,N) :: Vin
    INTEGER,INTENT(IN),DIMENSION(M,N) :: MASK
    REAL(SP),INTENT(OUT),DIMENSION(M,N1) :: OutL,OutR

    REAL(SP),DIMENSION(M,N) :: Din
    REAL(SP) :: TYP1,TYP2,TYP3,DVP1,DVP2,DVP3

    REAL(SP) :: VAN1, VAN2, RAT

    Din=0.0_SP
    DO I=Ibeg,Iend

     ! Fourth order MUSCL TVD with van Leer and minmod limiter
     ! The minmod limiter is used for the fourth order part in the fourth order scheme
          DO J=Jbeg-1,Jend+2
          TYP1=Vin(I,J-1)-Vin(I,J-2)
          TYP2=Vin(I,J)-Vin(I,J-1)
          TYP3=Vin(I,J+1)-Vin(I,J)

          if (TYP1.ge.0.0_SP) then
             DVP1=MAX(0.0_SP,MIN(TYP1,2.0*TYP2,2.0*TYP3))          
          else
             DVP1=MIN(0.0_SP,MAX(TYP1,2.0*TYP2,2.0*TYP3))          
          endif
          if (TYP2.ge.0.0_SP) then
             DVP2=MAX(0.0_SP,MIN(TYP2,2.0*TYP3,2.0*TYP1))          
          else
             DVP2=MIN(0.0_SP,MAX(TYP2,2.0*TYP3,2.0*TYP1))          
          endif
          if (TYP3.ge.0.0_SP) then
             DVP3=MAX(0.0_SP,MIN(TYP3,2.0*TYP1,2.0*TYP2))          
          else
             DVP3=MIN(0.0_SP,MAX(TYP3,2.0*TYP1,2.0*TYP2))          
          endif

          ! dry D-2 J-1 J D+1, lower-order
          IF(MASK(I,J-2)==0.OR.MASK(I,J+1)==0)THEN
            TYP2=Vin(I,J)-Vin(I,J-1)
            TYP1=TYP2
            TYP3=TYP2
            if (TYP1.ge.0.0_SP) then
               DVP1=MAX(0.0_SP,MIN(TYP1,2.0*TYP2,2.0*TYP3))
            else
               DVP1=MIN(0.0_SP,MAX(TYP1,2.0*TYP2,2.0*TYP3))
            endif
            if (TYP2.ge.0.0_SP) then
               DVP2=MAX(0.0_SP,MIN(TYP2,2.0*TYP3,2.0*TYP1))
            else
               DVP2=MIN(0.0_SP,MAX(TYP2,2.0*TYP3,2.0*TYP1))
            endif
            if (TYP3.ge.0.0_SP) then
               DVP3=MAX(0.0_SP,MIN(TYP3,2.0*TYP1,2.0*TYP2))
            else
               DVP3=MIN(0.0_SP,MAX(TYP3,2.0*TYP1,2.0*TYP2))
            endif
          ! here actually DVP1=DVP2=DVP3=TYP2
          ENDIF    
          ! dry J-2 D-1 D J+1, zero gradient
          IF(MASK(I,J-1)==0.OR.MASK(I,J)==0)THEN
            DVP1=ZERO
            DVP2=ZERO
            DVP3=ZERO
          ENDIF

          Din(I,J)=TYP2-1.0_SP/6.0_SP*(DVP3-2.0_SP*DVP2+DVP1)
          ENDDO
    
     ! Fourth order MUSCL TVD with van Leer and minmod limiter
     ! The van-Leer limiter is used for the third order part in the fourth order scheme
     DO J=Jbeg,Jend+1
          TMP1 = Din(I,J-1);   TMP2 = Din(I,J);

          IF( ABS(TMP1).le.SMALL ) TMP1 = SMALL*SIGN( 1.0_SP, TMP1 )
          IF( ABS(TMP2).le.SMALL ) TMP2 = SMALL*SIGN( 1.0_SP, TMP2 )
          RAT = TMP2/TMP1	   
          VAN1 = 0.0_SP;
          IF( abs(1.0 + RAT) > SMALL ) VAN1 = ( RAT + ABS(RAT) )/( 1.0_SP + RAT)

          RAT = TMP1/TMP2
          VAN2 = 0.0_SP;
          IF( abs(1.0 + RAT) > SMALL ) VAN2 = ( RAT + ABS(RAT) )/( 1.0_SP + RAT)

          OutL(I,J) = Vin(I,J-1) + ( 1.0_SP/6.0_SP )*( VAN1*TMP1 + 2.0_SP*VAN2*TMP2 )

      !!!!!!!!!!!!!!!

          TMP1 = Din(I,J);   TMP2 = Din(I,J+1);

          IF( ABS(TMP1).le.SMALL ) TMP1 = SMALL*SIGN( 1.0_SP, TMP1 )
          IF( ABS(TMP2).le.SMALL ) TMP2 = SMALL*SIGN( 1.0_SP, TMP2 )
          RAT = TMP2/TMP1
          VAN1 = 0.0_SP;
          IF( abs(1.0 + RAT) > SMALL ) VAN1 = ( RAT + ABS(RAT) )/( 1.0_SP + RAT)

          RAT = TMP1/TMP2
          VAN2 = 0.0_SP;
          IF( abs(1.0 + RAT) > SMALL ) VAN2 = ( RAT + ABS(RAT) )/( 1.0_SP + RAT)

          OutR(I,J) = Vin(I,J) - ( 1.0_SP/6.0_SP )*( 2.0_SP*VAN1*TMP1 + VAN2*TMP2 )
       ENDDO

    ENDDO     

END SUBROUTINE CONSTRUCT_HO_Y


