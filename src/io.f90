!-------------------------------------------------------------------------------------
!     OUTPUT
!           - Calculate statistics to display (Statistics)
!           - Outputs timesteps to individual files (Preview)
!-------------------------------------------------------------------------------------
SUBROUTINE OUTPUT
    USE GLOBAL
    IMPLICIT NONE

      ! Update screen
      SCREEN_COUNT=SCREEN_COUNT+DT

      ! Calculate statistics to display
      IF(SCREEN_COUNT>=SCREEN_INTV)THEN
            SCREEN_COUNT=SCREEN_COUNT-SCREEN_INTV
            CALL STATISTICS
      ENDIF

      ! Output individual eta_XXXXX files
      IF(TIME>=PLOT_START_TIME)THEN
            PLOT_COUNT=PLOT_COUNT+DT
            IF(PLOT_COUNT>=PLOT_INTV)THEN
                  PLOT_COUNT=PLOT_COUNT-PLOT_INTV
                  CALL PREVIEW
            ENDIF
      ENDIF 
END SUBROUTINE OUTPUT

!-------------------------------------------------------------------------------------
!     READ_INPUT 
!           - Read in all data from input.txt file
!-------------------------------------------------------------------------------------

SUBROUTINE READ_INPUT
      ! Import necessary modules
      USE GLOBAL
      USE INPUT_READ

      !! Variable declarations
      IMPLICIT NONE
      CHARACTER(LEN=80) FILE_NAME
      CHARACTER(LEN=80) MKFOLDER
      INTEGER::LINE
      INTEGER :: ierr
      INTEGER :: I_comp
      LOGICAL :: INPUT_PHASE = .FALSE.
      ! ìnput.txt` FILE NAME
            CHARACTER(LEN=80)::INPUT_NAME=''
      ! RESULT_FOLDER
            CHARACTER(LEN=80)::FDIR=' '

      ! Set up MPI if needed
#       if defined (PARALLEL)
            CALL MPI_COMM_SIZE (MPI_COMM_WORLD, nprocs, ier)  
            CALL MPI_COMM_RANK (MPI_COMM_WORLD, myid, ier)
#       endif

      ! Set up RESULT_FOLDER, time_dt.out, and LOG.txt file
      FDIR=TRIM(RESULT_FOLDER)
      OPEN(10000,FILE='time_dt.out',STATUS='UNKNOWN')
      OPEN(3,FILE='LOG.txt')   

      ! Get name of `input.txt` file
      CALL GETARG(1,INPUT_NAME) 
            ! defaults to `input.txt` if not given
            if (INPUT_NAME .eq. '') Then
                  FILE_NAME='input.txt'
            ! uses argument name otherwise
            Else
                  FILE_NAME=INPUT_NAME
            endif
      INPUT_FILE_NAME=FILE_NAME

      !@! START LOG FILE
      CALL READ_STRING(TITLE,FILE_NAME,'TITLE',ierr)
      ! Default to 'TEST RUN' if no name given
      IF(ierr==1)THEN
        TITLE='---TEST RUN---'
      ENDIF
      ! Start up log file (parallel)
#       if defined (PARALLEL)
            if (myid.eq.0) WRITE(3,*)'-------------- LOG FILE -----------------'
            if (myid.eq.0) WRITE(3,*)TITLE
            if (myid.eq.0) WRITE(3,*)' --------------input start --------------'
      ! Start up log file (not parallel)
#       else
            WRITE(3,*)'-------------- LOG FILE -----------------'
            WRITE(3,*)TITLE
            WRITE(3,*)' --------------input start --------------'  
            WRITE(3,*)'                                         '   
#       endif

      !@! Parallel Info
#             if defined (PARALLEL)
                  if (myid.eq.0) WRITE(3,*)'                                         '
                  if (myid.eq.0) WRITE(3,*)'-------------- PARALLEL -----------------'    
#             endif
            ! Print PX and PY Parallel Information
#             if defined (PARALLEL)
                  ! PX Information
                  CALL READ_INTEGER(PX,FILE_NAME,'PX',ierr)
                  IF(ierr == 1) THEN
                        PX = 1
                  ENDIF
                  ! PY Information
                  CALL READ_INTEGER(PY,FILE_NAME,'PY',ierr)  
                  IF(ierr == 1) THEN
                        PY = 1
                  ENDIF
                  ! Print PX and PY info to log
                  if (myid.eq.0) WRITE(3,'(A7,I3,A7,I3)') 'PX   =',PX,'PY   =', PY
#             endif

      !@! Grid Info
#             if defined (PARALLEL)
                  if (myid.eq.0) WRITE(3,*)'                                         '
                  if (myid.eq.0) WRITE(3,*)'-------------- GRID INFO -----------------'
#             else
                  WRITE(3,*)'                                         '
                  WRITE(3,*)'-------------- GRID INFO -----------------'   
#             endif

            ! Get and check Mglob dimension
            CALL READ_INTEGER(Mglob,FILE_NAME,'Mglob',ierr)
            IF(ierr==1)THEN
#                   if defined (PARALLEL)
                        if (myid.eq.0) THEN
                              WRITE(*,'(A40,A40)')'Mglob:', 'NOT DEFINED, STOP'
                              WRITE(3,'(A40,A40)')'Mglob:', 'NOT DEFINED, STOP'
                        endif
                        call MPI_FINALIZE ( ier )
#                   else
                        WRITE(*,'(A40,A40)')'Mglob:', 'NOT DEFINED, STOP'
                        WRITE(3,'(A40,A40)')'Mglob:', 'NOT DEFINED, STOP'
#                   endif
                        STOP
            ENDIF

            ! Get and check Nglob dimension
            CALL READ_INTEGER(Nglob,FILE_NAME,'Nglob',ierr)
            IF(ierr==1)THEN
#                   if defined (PARALLEL)
                        if (myid.eq.0) THEN
                              WRITE(*,'(A40,A40)')'Nglob:', 'NOT DEFINED, STOP'
                              WRITE(3,'(A40,A40)')'Nglob:', 'NOT DEFINED, STOP'
                        endif
                        call MPI_FINALIZE ( ier )
#                   else
                        WRITE(*,'(A40,A40)')'Nglob:', 'NOT DEFINED, STOP'
                        WRITE(3,'(A40,A40)')'Nglob:', 'NOT DEFINED, STOP'
#                   endif
                        STOP
            ENDIF

            ! Print grid dimensions
#             if defined (PARALLEL)
                  if (myid.eq.0) WRITE(3,'(A7,I8,A7,I8)') 'Mglob=',Mglob,'Nglob=', Nglob
#             else
                  WRITE(3,'(A7,I8,A7,I8)') 'Mglob=',Mglob,'Nglob=', Nglob
#             endif

            ! Cartesian Settings
            CALL READ_FLOAT(DX,FILE_NAME,'DX',ierr)
            ! Error message: Did you intend to use Spherical?
            IF(ierr==1)THEN
                PRINT *,"Did you intend to use Spherical Coordinates?"
#                         if defined (PARALLEL)
                        if (myid.eq.0) THEN
                            WRITE(*,'(A40,A40)')'DX:', 'NOT DEFINED, STOP'
                            WRITE(3,'(A40,A40)')'DX:', 'NOT DEFINED, STOP'
                        endif
                        call MPI_FINALIZE ( ier )
#                         else
                        WRITE(*,'(A40,A40)')'DX:', 'NOT DEFINED, STOP'
                        WRITE(3,'(A40,A40)')'DX:', 'NOT DEFINED, STOP'
#                         endif
                STOP
            ENDIF

            CALL READ_FLOAT(DY,FILE_NAME,'DY',ierr)
            ! Error message: Undefined DX and DY
            IF(ierr==1)THEN
#                         if defined (PARALLEL)
                        if (myid.eq.0) THEN
                            WRITE(*,'(A40,A40)')'DY:', 'NOT DEFINED, STOP'
                            WRITE(3,'(A40,A40)')'DY:', 'NOT DEFINED, STOP'
                        endif
                        call MPI_FINALIZE ( ier )
#                         else
                        WRITE(*,'(A40,A40)')'DY:', 'NOT DEFINED, STOP'
                        WRITE(3,'(A40,A40)')'DY:', 'NOT DEFINED, STOP'
#                         endif
                STOP
            ENDIF
            
            ! Print DX and DY out
#                   if defined (PARALLEL)
                if (myid.eq.0) WRITE(3,'(A4,F12.2,A4,F12.2)')'DX=',DX,'DY=',DY
#                   else
                WRITE(3,'(A4,F12.2,A4,F12.2)')'DX=',DX,'DY=',DY
#                   endif


      !@! DEPTH INFO
            CALL READ_STRING(DEPTH_TYPE,FILE_NAME,'DEPTH_TYPE',ierr)

            ! Error handling for DEPTH_TYPE > default to flat                     
            IF(ierr==1)THEN
                  DEPTH_TYPE = 'FLAT'
            ENDIF

            ! Write out DEPTH_TYPE
#             if defined (PARALLEL)
                  if (myid.eq.0) WRITE(3,'(A12,A50)')'DEPTH_TYPE:', DEPTH_TYPE
#             else
                  WRITE(3,'(A12,A50)')'DEPTH_TYPE:', DEPTH_TYPE
#             endif

            ! DEPTH_TYPE = SLOPE
            IF(DEPTH_TYPE(1:3)=='SLO')THEN
                  CALL READ_FLOAT(DEPTH_FLAT,FILE_NAME,'DEPTH_FLAT',ierr) 
                  ! Error handling for DEPTH_TYPE=SLOPE DEPTH_FLAT variable
                  IF(ierr==1)THEN
                        DEPTH_FLAT = 10.0_SP
                  ENDIF

                  CALL READ_FLOAT(SLP,FILE_NAME,'SLP',ierr) 
                  ! Error handling for DEPTH_TYPE=SLOPE SLP variable
                  IF(ierr==1)THEN
                        SLP = 0.1_SP
                  ENDIF 

                  CALL READ_FLOAT(Xslp,FILE_NAME,'Xslp',ierr) 
                  ! Error handling for DEPTH_TYPE=SLOPE Xslp variable
                  IF(ierr==1)THEN
                        Xslp = 0.0_SP
                  ENDIF 

                  ! Write out DEPTH_FLAT = SLP parameters
#                   if defined (PARALLEL)
                        if (myid.eq.0) WRITE(3,'(A10,F12.2)')'DEPTH_FLAT=', DEPTH_FLAT 
                        if (myid.eq.0) WRITE(3,'(A5,F12.2)')'SLP=', SLP
                        if (myid.eq.0) WRITE(3,'(A6,F12.2)')'Xslp=', Xslp  
#                   else
                        WRITE(3,'(A10,F12.2)')'DEPTH_FLAT=', DEPTH_FLAT 
                        WRITE(3,'(A5,F12.2)')'SLP=', SLP
                        WRITE(3,'(A6,F12.2)')'Xslp=', Xslp  
#                   endif
            ENDIF  

      !@! DEPTH CORRECTION
            BATHY_CORRECTION = .FALSE. 

      !@! TIME INFORMATION
#             if defined (PARALLEL)
                  if (myid.eq.0) WRITE(3,*)'                                         '
                  if (myid.eq.0) WRITE(3,*)'-------------- TIME INFO -----------------'
#             else
                  WRITE(3,*)'                                         '
                  WRITE(3,*)'-------------- TIME INFO -----------------'   
#             endif

            CALL READ_FLOAT(TOTAL_TIME,FILE_NAME,'TOTAL_TIME',ierr)
            ! Error handling for TOTAL_TIME variable
            IF(ierr==1)THEN
#                   if defined (PARALLEL)
                        if (myid.eq.0) THEN
                              WRITE(*,'(A40,A40)')'TOTAL_TIME:', 'NOT FOUND, STOP'
                              WRITE(3,'(A40,A40)')'TOTAL_TIME:', 'NOT FOUND, STOP'
                        endif
                        call MPI_FINALIZE ( ier )
#                   else
                        WRITE(*,'(A40,A40)')'TOTAL_TIME:', 'NOT FOUND, STOP'
                        WRITE(3,'(A40,A40)')'TOTAL_TIME:', 'NOT FOUND, STOP'
#                   endif
                  STOP
            ENDIF

            CALL READ_FLOAT(PLOT_START_TIME,FILE_NAME,'PLOT_START_TIME',ierr)
            ! Error handling for PLOT_START_TIME variable
            IF(ierr==1)THEN
                  PLOT_START_TIME = 0.0
            ENDIF

            CALL READ_FLOAT(PLOT_INTV,FILE_NAME,'PLOT_INTV',ierr)
            ! Error handling for PLOT_INTV variable
            IF(ierr==1)THEN
                  PLOT_INTV = 1.0
            ENDIF

            ! SCREEN_INTV set to 1
            SCREEN_INTV = 1.0

            ! Write out time information
#             if defined (PARALLEL)
                  if (myid.eq.0) WRITE(3,'(A12,F12.2)')'TOTAL_TIME=', TOTAL_TIME
                  if (myid.eq.0) WRITE(3,'(A12,F12.2)')'PLOT_INTV= ', PLOT_INTV
#             else
                  WRITE(3,'(A12,F12.2)')'TOTAL_TIME=', TOTAL_TIME
                  WRITE(3,'(A12,F12.2)')'PLOT_INTV= ', PLOT_INTV
#             endif


      !@! WAVEMAKER INFORMATION
#             if defined (PARALLEL)
                  if (myid.eq.0) WRITE(3,*)'                                         '
                  if (myid.eq.0) WRITE(3,*)'-------------- WAVEMAKER -----------------'
#             else
                  WRITE(3,*)'                                         '
                  WRITE(3,*)'-------------- WAVEMAKER -----------------'   
#             endif

            ! Get Wavemaker type
            CALL READ_STRING(WaveMaker,FILE_NAME,'WAVEMAKER',ierr)
            ! Error handling for no wavemaker
            IF(ierr==1)THEN
                  WaveMaker = 'nothing'
#                   if defined (PARALLEL)
                        if (myid.eq.0) THEN
                        WRITE(*,'(A40)')'No WaveMaker'
                        WRITE(3,'(A40)')'No WaveMaker'
                        endif
#                   else
                        WRITE(*,'(A40)')'No WaveMaker'
                        WRITE(3,'(A40)')'No WaveMaker'
#                   endif
            ENDIF

            ! Print out wavemaker type
#             if defined (PARALLEL)
                  if (myid.eq.0) WRITE(3,'(A11,A50)')'WAVEMAKER:', WAVEMAKER
#             else
                  WRITE(3,'(A11,A50)')'WAVEMAKER:', WAVEMAKER
#             endif
            
            

            ! WK_REG Wavemaker type
            IF(WaveMaker(1:6)=='WK_REG')THEN
                  CALL READ_FLOAT(Xc_WK,FILE_NAME,'Xc_WK',ierr)
                  ! Error handling for no Xc_WK variable
                  IF(ierr==1)THEN
#                         if defined (PARALLEL)
                              if (myid.eq.0) THEN
                                    WRITE(*,'(A60)')'Xc_WK NOT FOUND, specify Xc_WK in input.txt'
                                    WRITE(3,'(A60)')'Xc_WK NOT FOUND, specify Xc_WK in input.txt'
                              endif
                              call MPI_FINALIZE ( ier )
#                         else
                              WRITE(*,'(A60)')'Xc_WK NOT FOUND, specify Xc_WK in input.txt'
                              WRITE(3,'(A60)')'Xc_WK NOT FOUND, specify Xc_WK in input.txt'
#                         endif
                        STOP
                  ENDIF

                  ! Yc_WK for 1D
                  Yc_WK = ZERO


                  CALL READ_FLOAT(Tperiod,FILE_NAME,'Tperiod',ierr)
                  ! Error handling for no Tperiod variable
                  IF(ierr==1)THEN
#                         if defined (PARALLEL)
                              if (myid.eq.0) THEN
                                    WRITE(*,'(A60)')'Tperiod NOT FOUND, specify Tperiod in input.txt'
                                    WRITE(3,'(A60)')'Tperiod NOT FOUND, specify Tperiod in input.txt'
                              endif
                              call MPI_FINALIZE ( ier )
#                         else
                              WRITE(*,'(A60)')'Tperiod NOT FOUND, specify Tperiod in input.txt'
                              WRITE(3,'(A60)')'Tperiod NOT FOUND, specify Tperiod in input.txt'
#                         endif
                        STOP
                  ENDIF

                  CALL READ_FLOAT(AMP_WK,FILE_NAME,'AMP_WK',ierr)
                  ! Error handling for no AMP_WK variable
                  IF(ierr==1)THEN
#                         if defined (PARALLEL)
                              if (myid.eq.0) THEN
                                    WRITE(*,'(A60)')'AMP_WK NOT FOUND, specify AMP_WK in input.txt'
                                    WRITE(3,'(A60)')'AMP_WK NOT FOUND, specify AMP_WK in input.txt'
                              endif
                              call MPI_FINALIZE ( ier )
#                         else
                              WRITE(*,'(A60)')'AMP_WK NOT FOUND, specify AMP_WK in input.txt'
                              WRITE(3,'(A60)')'AMP_WK NOT FOUND, specify AMP_WK in input.txt'
#                         endif
                        STOP
                  ENDIF

                  CALL READ_FLOAT(DEP_WK,FILE_NAME,'DEP_WK',ierr)
                  ! Error handling for no DEP_WK variable
                  IF(ierr==1)THEN
#                         if defined (PARALLEL)
                              if (myid.eq.0) THEN
                                    WRITE(*,'(A60)')'DEP_WK NOT FOUND, specify DEP_WK in input.txt'
                                    WRITE(3,'(A60)')'DEP_WK NOT FOUND, specify DEP_WK in input.txt'
                              endif
                              call MPI_FINALIZE ( ier )
#                         else
                              WRITE(*,'(A60)')'DEP_WK NOT FOUND, specify DEP_WK in input.txt'
                              WRITE(3,'(A60)')'DEP_WK NOT FOUND, specify DEP_WK in input.txt'
#                         endif
                        STOP
                  ENDIF

                  ! Specify Theta_WK as 0
                  Theta_WK = 0.0_SP

                  ! Time_ramp as 0
                  Time_ramp = 0.0_SP

                

                  CALL READ_FLOAT(Delta_WK,FILE_NAME,'Delta_WK',ierr)
                  ! Error handling for no Delta_WK variable
                  IF(ierr==1)THEN
                        Delta_WK = 0.5_SP
                  ENDIF

                  ! Write out parameters
#                   if defined (PARALLEL)
                        if (myid.eq.0) WRITE(3,'(A10,F12.2)')'Xc_WK   =', Xc_WK
                        if (myid.eq.0) WRITE(3,'(A10,F12.2)')'Tperiod =', Tperiod
                        if (myid.eq.0) WRITE(3,'(A10,F12.2)')'AMP_WK  =', AMP_WK
                        if (myid.eq.0) WRITE(3,'(A10,F12.2)')'DEP_WK  =', DEP_WK
                        if (myid.eq.0) WRITE(3,'(A10,F12.2)')'Theta_WK=', Theta_WK
                        if (myid.eq.0) WRITE(3,'(A11,F12.2)')'Delta_WK=', Delta_WK
#                   else
                        WRITE(3,'(A10,F12.2)')'Xc_WK   =', Xc_WK
                        WRITE(3,'(A10,F12.2)')'Tperiod =', Tperiod
                        WRITE(3,'(A10,F12.2)')'AMP_WK  =', AMP_WK
                        WRITE(3,'(A10,F12.2)')'DEP_WK  =', DEP_WK
                        WRITE(3,'(A10,F12.2)')'Theta_WK=', Theta_WK
                        WRITE(3,'(A11,F12.2)')'Delta_WK=', Delta_WK
#                   endif
                        ENDIF  


      !@! SPONGE INFORMATION
#             if defined (PARALLEL)
                  if (myid.eq.0) WRITE(3,*)'                                         '
                  if (myid.eq.0) WRITE(3,*)'-------------- SPONGE -----------------'
#             else
                  WRITE(3,*)'                                         '
                  WRITE(3,*)'-------------- SPONGE -----------------'   
#             endif

            ! Get the type of sponge
                  CALL READ_LOGICAL(DIFFUSION_SPONGE,FILE_NAME,'DIFFUSION_SPONGE',ierr)
                  IF(ierr==1)THEN
                        DIFFUSION_SPONGE = .FALSE.
                  ENDIF

                  CALL READ_LOGICAL(DIRECT_SPONGE,FILE_NAME,'DIRECT_SPONGE',ierr)
                  IF(ierr==1)THEN
                        DIRECT_SPONGE = .FALSE.
                  ENDIF

                  CALL READ_LOGICAL(FRICTION_SPONGE,FILE_NAME,'FRICTION_SPONGE',ierr)
                  IF(ierr==1)THEN
                        FRICTION_SPONGE = .FALSE.
                  ENDIF

            ! Direct Sponge
                  IF(DIRECT_SPONGE)THEN

#                         if defined (PARALLEL)
                              if (myid.eq.0) THEN
                                    WRITE(*,'(A40)')'DIRECT_SPONGE IS USED'
                                    WRITE(3,'(A40)')'DIRECT_SPONGE IS USED'
                              endif
#                         else
                              WRITE(*,'(A40)')'DIRECT_SPONGE IS USED'
                              WRITE(3,'(A40)')'DIRECT_SPONGE IS USED'
#                         endif
                  ENDIF

            ! Diffusion Sponge PARAMETERS
            IF(DIFFUSION_SPONGE)THEN
                   Csp = 0.1_SP
            ENDIF

            ! Friction Sponge
            IF(FRICTION_SPONGE)THEN
#                   if defined (PARALLEL)
                        if (myid.eq.0) THEN
                              WRITE(*,'(A40)')'FRICTION_SPONGE IS USED'
                              WRITE(3,'(A40)')'FRICTION_SPONGE IS USED'
                        endif
#                   else
                        WRITE(*,'(A40)')'FRICTION_SPONGE IS USED'
                        WRITE(3,'(A40)')'FRICTION_SPONGE IS USED'
#                   endif
                        ! Set CDSponge
                        CALL READ_FLOAT(CDsponge,FILE_NAME,'CDsponge',ierr)
                        IF(ierr==1)THEN
                              CDsponge = 5.0_SP
                        ENDIF

#                               if defined (PARALLEL)
                                    if (myid.eq.0) WRITE(3,'(A26,F12.2)')'FRICTION_SPONGE CDsponge=', CDsponge
#                               else
                                    WRITE(3,'(A26,F12.2)')'FRICTION_SPONGE CDsponge=', CDsponge
#                               endif
            ENDIF  

            ! Common sponge parameters
            IF(DIFFUSION_SPONGE.OR.DIRECT_SPONGE.OR.FRICTION_SPONGE)THEN
                  ! WEST WIDTH
                  CALL READ_FLOAT(Sponge_west_width,FILE_NAME,'Sponge_west_width',ierr)
                        IF(ierr==1)THEN
                              Sponge_west_width = 0.0_SP
                        ENDIF
                  ! EAST WIDTH
                  CALL READ_FLOAT(Sponge_east_width,FILE_NAME,'Sponge_east_width',ierr)
                        IF(ierr==1)THEN
                              Sponge_east_width = 0.0_SP
                        ENDIF
                  ! SOUTH WIDTH
                  CALL READ_FLOAT(Sponge_south_width,FILE_NAME,'Sponge_south_width',ierr)
                        IF(ierr==1)THEN
                              Sponge_south_width = 0.0_SP
                        ENDIF
                  ! NORTH WIDTH
                  CALL READ_FLOAT(Sponge_north_width,FILE_NAME,'Sponge_north_width',ierr)
                  IF(ierr==1)THEN
                        Sponge_north_width = 0.0_SP
                  ENDIF



                  ! Write out parameters
#                   if defined (PARALLEL)
                        if (myid.eq.0) WRITE(3,'(A20,F12.2)')'Sponge_west_width =', Sponge_west_width
                        if (myid.eq.0) WRITE(3,'(A20,F12.2)')'Sponge_east_width =', Sponge_east_width
                        if (myid.eq.0) WRITE(3,'(A20,F12.2)')'Sponge_south_width=', Sponge_south_width
                        if (myid.eq.0) WRITE(3,'(A20,F12.2)')'Sponge_north_width=', Sponge_north_width
                        if (myid.eq.0) WRITE(3,'(A20,F12.2)')'R_sponge          =', R_sponge
                        if (myid.eq.0) WRITE(3,'(A20,F12.2)')'A_sponge          =', A_sponge
#                   else
                        WRITE(3,'(A20,F12.2)')'Sponge_west_width =', Sponge_west_width
                        WRITE(3,'(A20,F12.2)')'Sponge_east_width =', Sponge_east_width
                        WRITE(3,'(A20,F12.2)')'Sponge_south_width=', Sponge_south_width
                        WRITE(3,'(A20,F12.2)')'Sponge_north_width=', Sponge_north_width
                        WRITE(3,'(A20,F12.2)')'R_sponge          =', R_sponge
                        WRITE(3,'(A20,F12.2)')'A_sponge          =', A_sponge
#                   endif
            ENDIF 

! WAVEMAKER CURRENT BALANCE
        WaveMakerCurrentBalance=.FALSE.
# if defined (PARALLEL)
      if (myid.eq.0) WRITE(3,'(A40)')'No WavemakerCurrentBalance'
# else
       WRITE(3,'(A40)')'No WavemakerCurrentBalance'
# endif

     

      !@! PHYSICS INFORMATION
#             if defined (PARALLEL)
                  if (myid.eq.0) WRITE(3,*)'                                         '
                  if (myid.eq.0) WRITE(3,*)'-------------- PHYSICS -----------------'
#             else
                  WRITE(3,*)'                                         '
                  WRITE(3,*)'-------------- PHYSICS -----------------'   
#             endif

      !@! DISPERSION INFORMATION
            ! Dispersion type
            DISPERSION = .TRUE.
            Gamma1 = 1.0_SP
            Gamma2 = 1.0_SP
            Beta_ref = - 0.531_SP
            Gamma3 = 1.0_SP



      !@! VISCOSITY BREAKING INFORMATION
            VISCOSITY_BREAKING = .TRUE.
            IF(ROLLER) VISCOSITY_BREAKING = .TRUE.
            SWE_ETA_DEP = 0.80_SP
            Cd_fixed = 0.0_SP


      !@! NUMERICS INFORMATION
            Time_Scheme = 'Runge_Kutta'
            CONSTR='HLLC'
            HIGH_ORDER='FOURTH'  
            CFL = 0.5_SP
  
      !@! Froude Number Cap
            FroudeCap = 3.0_SP

      !@! MinDepth
            CALL READ_FLOAT(MinDepth,FILE_NAME,'MinDepth',ierr)
            IF(ierr==1)THEN
                  MinDepth = 0.1_SP
            ENDIF

            MinDepthFrc = 0.1_SP

            !  merge two parameters into the minimum one, change to min according to Harris 11/13/2023
                  MinDepthFrc=MIN(MinDepthFrc,MinDepth)
                  MinDepth=MinDepthFrc

#             if defined (PARALLEL)
                  if (myid.eq.0) WRITE(3,'(A40)')'USE MIN(MinDepthFrc, MinDepth)'
#             else
                  WRITE(3,'(A40)')'USE MIN(MinDepthFrc, MinDepth)'
#             endif

            ! Print out parameters
#             if defined (PARALLEL)
                  if (myid.eq.0) WRITE(3,'(A10,F12.6)')'MinDepth=', MinDepth
#             else
                  WRITE(3,'(A10,F12.6)')'MinDepth=', MinDepth
#             endif

#             if defined (PARALLEL)
                  if (myid.eq.0) WRITE(3,'(A13,F12.6)')'MinDepthFrc=', MinDepthFrc
#             else
                  WRITE(3,'(A13,F12.2)')'MinDepthFrc=', MinDepthFrc
#             endif

      !@! WAVE BREAKING
#             if defined (PARALLEL)
                  if (myid.eq.0) WRITE(3,*)'                                         '
                  if (myid.eq.0) WRITE(3,*)'--------- WAVE BREAKING -----------------'
#             else
                  WRITE(3,*)'                                         '
                  WRITE(3,*)'-------- WAVE BREAKING -----------------'   
#             endif

	
            IF(VISCOSITY_BREAKING) SHOW_BREAKING = .TRUE.
            IF(SHOW_BREAKING)THEN
                  !$! Cbrk1 (Float): Breaking parameter 1
                  CALL READ_FLOAT(Cbrk1,FILE_NAME,'Cbrk1',ierr)
                  IF(ierr==1)THEN
                        Cbrk1 = 0.65_SP
            ENDIF

            IF(VISCOSITY_BREAKING)THEN
#                   if defined (PARALLEL)
                        if (myid.eq.0) WRITE(3,'(A8,F12.6)')'Cbrk1 =', Cbrk1
#                   else
                        WRITE(3,'(A8,F12.6)')'Cbrk1 =', Cbrk1
#                   endif
            ENDIF

            !$! Cbrk2 (Float): Breaking parameter 2
            CALL READ_FLOAT(Cbrk2,FILE_NAME,'Cbrk2',ierr)
            IF(ierr==1)THEN
                  Cbrk2 = 0.35_SP
            ENDIF

            ! Write out Cbrk2 parameter
            IF(VISCOSITY_BREAKING)THEN
#                   if defined (PARALLEL)
                        if (myid.eq.0) WRITE(3,'(A8,F12.6)')'Cbrk2 =', Cbrk2
#                   else
                        WRITE(3,'(A8,F12.6)')'Cbrk2 =', Cbrk2
#                   endif
            ENDIF



            ! Error handling for viscosity breaking and wavemaker viscosity
	      IF( VISCOSITY_BREAKING .AND. WAVEMAKER_VIS ) THEN
#                   if defined (PARALLEL)
                        IF (myid.eq.0) then
                              WRITE(*,*) "==============================================="
                              WRITE(*,*)  "STOP :: VISCOSITY_BREAKING=T, WAVEMAKER_VIS=T"
                              WRITE(*,*) "==============================================="
                        ENDIF
                        call MPI_FINALIZE ( ier )
#                   else
                        WRITE(*,*) "==============================================="
                        WRITE(*,*)  "STOP :: VISCOSITY_BREAKING=T, WAVEMAKER_VIS=T"
                        WRITE(*,*) "==============================================="
                        STOP
#                   endif
            ENDIF


      !@! WAVE AVERAGED PROPERTIES
#             if defined (PARALLEL)
                  if (myid.eq.0) WRITE(3,*)'                                         '
                  if (myid.eq.0) WRITE(3,*)'-------- WAVE-AVERAGED PROPERTY -----------------'
#             else
                  WRITE(3,*)'                                         '
                  WRITE(3,*)'------- WAVE-AVERAGED PROPERTY -----------------'   
#             endif

            !#! T_INTV_mean (Float): Time to take mean
            CALL READ_FLOAT(T_INTV_mean,FILE_NAME,'T_INTV_mean',ierr)
            IF(ierr==1)THEN
                  T_INTV_mean = LARGE
#                   if defined (PARALLEL)
                        if (myid.eq.0) THEN
                              WRITE(*,'(A40)')'T_INTV_mean Default:  LARGE'
                              WRITE(3,'(A40)')'T_INTV_mean Default:  LARGE'
                        endif
#                   else
                        WRITE(*,'(A40)')'T_INTV_mean Default:  LARGE'
                        WRITE(3,'(A40)')'T_INTV_mean Default:  LARGE'
#                   endif
            ENDIF

            !#! STEADY_TIME (Float): Time until steady
	      CALL READ_FLOAT(STEADY_TIME,FILE_NAME,'STEADY_TIME',ierr)
            IF(ierr==1)THEN
                  STEADY_TIME = LARGE
#                   if defined (PARALLEL)
                        if (myid.eq.0) THEN
                              WRITE(*,'(A40)')'STEADY_TIME Default:  LARGE'
                              WRITE(3,'(A40)')'STEADY_TIME Default:  LARGE'
                        endif
#                   else
                        WRITE(*,'(A40)')'STEADY_TIME Default:  LARGE'
                        WRITE(3,'(A40)')'STEADY_TIME Default:  LARGE'
#                   endif
            ENDIF


            ! Print out parameters
#             if defined (PARALLEL)
                  if (myid.eq.0) WRITE(3,'(A14,F12.6)')'T_INTV_mean =', T_INTV_mean
                  if (myid.eq.0) WRITE(3,'(A14,F12.6)')'STEADY_TIME =', STEADY_TIME
#             else
                  WRITE(3,'(A14,F12.6)')'T_INTV_mean =', T_INTV_mean
                  WRITE(3,'(A14,F12.6)')'STEADY_TIME =', STEADY_TIME
#             endif
	
      
      !@! OUTPUT INFO
#             if defined (PARALLEL)
                  if (myid.eq.0) WRITE(3,*)'                                         '
                  if (myid.eq.0) WRITE(3,*)'-------------- OUTPUT INFO -----------------'
#             else
                  WRITE(3,*)'                                         '
                  WRITE(3,*)'-------------- OUTPUT INFO -----------------'   
#             endif

            !$! RESULT_FOLDER (String): folder where results will go
            CALL READ_STRING(RESULT_FOLDER,FILE_NAME,'RESULT_FOLDER',ierr)
            ! Default to 'output` name
            IF(ierr==1)THEN
                  RESULT_FOLDER = './output/'
            ENDIF
            ! Print out RESULT_FOLDER
#             if defined (PARALLEL)
                  if (myid.eq.0) WRITE(3,'(A15,A50)')'RESULT_FOLDER:', RESULT_FOLDER
#             else
                  WRITE(3,'(A15,A50)')'RESULT_FOLDER:', RESULT_FOLDER
#             endif

      !@! BINARY/ASCII OUTPUTS
#             if defined (PARALLEL)
                  !$! FIELD_IO_TYPE (String): type for output (binary of ascii)
                  CALL READ_STRING(FIELD_IO_TYPE,FILE_NAME,'FIELD_IO_TYPE',ierr)
                  ! Default to ASCII
                  IF(ierr.EQ.1) FIELD_IO_TYPE = 'ASCII'
                  
                  IF (myid.EQ.0) WRITE(3,*) 'FIELD_IO_TYPE = ' , FIELD_IO_TYPE
#             endif
 
      !@! CREATION OF RESULT FOLDER
            MKFOLDER = "mkdir -p "//TRIM(RESULT_FOLDER)
#             if defined (PARALLEL)
                  IF (myid.eq.0) THEN
#                         if defined (INTEL)
                              RES = SYSTEM(TRIM(MKFOLDER))
#                         else
                              CALL SYSTEM(TRIM(MKFOLDER))
#                         endif
                  ENDIF
#             else
#                   if defined(INTEL)
                        RES = SYSTEM(TRIM(MKFOLDER))
#                   else
                        CALL SYSTEM(TRIM(MKFOLDER))
#                   endif
#             endif


      !@! OUTPUT PARAMETERS
            !$! OUTPUT_RES
            CALL READ_INTEGER(OUTPUT_RES,FILE_NAME,'OUTPUT_RES',ierr)
            IF(ierr==1)THEN
                  OUTPUT_RES = 1
#                   if defined (PARALLEL)
                        if (myid.eq.0) THEN
                              WRITE(*,'(A60)')'OUTPUT_RES NOT FOUND, OUTPUT_RES=1: full resolution'
                              WRITE(3,'(A60)')'OUTPUT_RES NOT FOUND, OUTPUT_RES=1: full resolution'
                        endif
#                   else
                        WRITE(*,'(A60)')'OUTPUT_RES NOT FOUND, OUTPUT_RES=1: full resolution'
                        WRITE(3,'(A60)')'OUTPUT_RES NOT FOUND, OUTPUT_RES=1: full resolution'
#                   endif
            ENDIF

            ! Print out OUTPUT_RES
#             if defined (PARALLEL)
                  if (myid.eq.0) WRITE(3,'(A15,I10)')'OUTPUT_RES',OUTPUT_RES
#             else
                  WRITE(3,'(A15,I10)')'OUTPUT_RES',OUTPUT_RES
#             endif 

      !@! FINAL ROLL CALL OF VARIABLES
            CALL READ_LOGICAL(OUT_DEPTH,FILE_NAME,'DEPTH_OUT',ierr)
            CALL READ_LOGICAL(OUT_U,FILE_NAME,'U',ierr)
            CALL READ_LOGICAL(OUT_V,FILE_NAME,'V',ierr)
            CALL READ_LOGICAL(OUT_ETA,FILE_NAME,'ETA',ierr)

      !@! EtaBlowUp Value
      IF(ierr==1)THEN
            EtaBlowVal = 10.0_SP
#             if defined (PARALLEL)
                  if (myid.eq.0) THEN
                        WRITE(*,'(A40)')'EtaBlowVal Default:  100xmax_depth'
                        WRITE(3,'(A40)')'EtaBlowVal Default:  100xmax_depth'
                  endif
#             else
                  WRITE(*,'(A40)')'EtaBlowVal Default:  100xmax_depth'
                  WRITE(3,'(A40)')'EtaBlowVal Default:  100xmax_depth'
#             endif
       ENDIF

      !@! FINAL WRITING OF ALL VARIABLES
#             if defined (PARALLEL)
                  if (myid.eq.0)   then
#             endif
                  WRITE(3,'(A15,L2)')'OUT_DEPTH',OUT_DEPTH
                  WRITE(3,'(A15,L2)')'OUT_U',OUT_U
                  WRITE(3,'(A15,L2)')'OUT_V',OUT_V
                  WRITE(3,'(A15,L2)')'OUT_ETA',OUT_ETA

#             if defined (PARALLEL)
                  endif
#             endif

      !@! INPUT ENED
#             if defined (PARALLEL)
                  if (myid.eq.0) WRITE(3,*)'                                         '
                  if (myid.eq.0) WRITE(3,*)' --------------input end --------------' 
                  if (myid.eq.0) WRITE(3,*)'                                         '
#             else
                  WRITE(3,*)'                                         '
                  WRITE(3,*)' --------------input end --------------' 
                  WRITE(3,*)'                                         '
#             endif
END SUBROUTINE READ_INPUT



!-------------------------------------------------------------------------------------
!
!    PREVIEW is subroutine for print-out of field data
!
!  HISTORY:
!    05/01/2010  Fengyan Shi
!    06/01/2015  Young-Kwang Choi, change file number to 5 digits, 
!                        such as eta_00001
!
!-------------------------------------------------------------------------------------
SUBROUTINE PREVIEW
      ! Import necessary modules
      USE GLOBAL
      IMPLICIT NONE

      !! Variable Declarations
      CHARACTER(LEN=80)::FILE_NAME=' '
      CHARACTER(LEN=80)::FILE_NAME_MEAN=' '
      CHARACTER(LEN=80)::TMP_NAME=' '
      CHARACTER(LEN=80)::FDIR=' '

      ! RESULT_FOLDER from earlier
      FDIR=TRIM(RESULT_FOLDER)
      ! Iteration counter
      ICOUNT=ICOUNT+1

      !! Display progress of iteration
#             if defined (PARALLEL)
                  if (myid.eq.0)then
                  WRITE(3,102)'PRINTING FILE NO.', icount, ' TIME/TOTAL: ', TIME,'/',Total_Time
                  WRITE(*,102)'PRINTING FILE NO.', icount, ' TIME/TOTAL: ', TIME,'/',Total_Time        
                  endif
#             else
                  WRITE(*,102)'PRINTING FILE NO.', icount, ' TIME/TOTAL: ', TIME,'/',Total_Time
#             endif

            102     FORMAT(A20,I6,A14,F12.3,A2,F12.3)


            itmp1=mod(icount/10000,10)
            itmp2=mod(icount/1000,10)
            itmp3=mod(icount/100,10)
            itmp4=mod(icount/10,10)
            itmp5=mod(icount,10)

            write(file_name(1:1),'(I1)')itmp1
            write(file_name(2:2),'(I1)')itmp2
            write(file_name(3:3),'(I1)')itmp3
            write(file_name(4:4),'(I1)')itmp4
            write(file_name(5:5),'(I1)')itmp5   !ykchoi
      
      
      !! Things to do on first iteration
            IF(ICOUNT==1)THEN
            IF(OUT_DEPTH.OR.BREAKWATER)THEN
                  TMP_NAME = TRIM(FDIR)//'dep.out'
                  call PutFile(TMP_NAME,DEPTH)
                  TMP_NAME = TRIM(FDIR)//'cd_breakwater.out'
                  call PutFile(TMP_NAME,CD_breakwater)
            ENDIF
            ENDIF

      !! Write out time to time_dt file
            write(10000,*)time, dt
      
      !! Write out variables of iterest at each timestep
            IF(OUT_ETA)THEN
                  TMP_NAME = TRIM(FDIR)//'eta_'//TRIM(FILE_NAME)
                  call PutFile(TMP_NAME,Eta)
            ENDIF

#             if defined (UseEtaScreen)
            IF(OUT_ETAscreen)THEN
                  TMP_NAME = TRIM(FDIR)//'etasrn_'//TRIM(FILE_NAME)
                  call PutFile(TMP_NAME,EtaScreen)
            ENDIF
#             endif

            IF(OUT_Hmax)THEN
                  TMP_NAME = TRIM(FDIR)//'hmax_'//TRIM(FILE_NAME)
                  call PutFile(TMP_NAME,HeightMax)
            ENDIF

            IF(OUT_Hmin)THEN
                  TMP_NAME = TRIM(FDIR)//'hmin_'//TRIM(FILE_NAME)
                  call PutFile(TMP_NAME,HeightMin)
            ENDIF

            IF(OUT_Umax)THEN
                  TMP_NAME = TRIM(FDIR)//'umax_'//TRIM(FILE_NAME)
                  call PutFile(TMP_NAME,VelocityMax)
            ENDIF
            
            IF(OUT_MFmax)THEN                                                                                            
                  TMP_NAME = TRIM(FDIR)//'MFmax_'//TRIM(FILE_NAME)                                                          
                  call PutFile(TMP_NAME,MomentumFluxMax)                                                                              
            ENDIF      
            
            IF(OUT_VORmax)THEN                                                                                            
                  TMP_NAME = TRIM(FDIR)//'VORmax_'//TRIM(FILE_NAME)                                                          
                  call PutFile(TMP_NAME,VorticityMax)                                                                              
            ENDIF            
            
            IF(OUT_U)THEN
                  TMP_NAME = TRIM(FDIR)//'u_'//TRIM(FILE_NAME)
                  call PutFile(TMP_NAME,U)
            ENDIF

            IF(OUT_V)THEN
                  TMP_NAME = TRIM(FDIR)//'v_'//TRIM(FILE_NAME)
                  call PutFile(TMP_NAME,V)
            ENDIF

            IF(OUT_MASK)THEN
                  TMP_NAME = TRIM(FDIR)//'mask_'//TRIM(FILE_NAME)
                  Int2Flo=MASK
                  call PutFile(TMP_NAME,Int2Flo)
            ENDIF

            IF(OUT_MASK9)THEN
                  TMP_NAME = TRIM(FDIR)//'mask9_'//TRIM(FILE_NAME)
                  Int2Flo=MASK9
                  call PutFile(TMP_NAME,Int2Flo)
            ENDIF

      210   FORMAT(5000I3)

            IF(OUT_P)THEN
                  TMP_NAME = TRIM(FDIR)//'p_'//TRIM(FILE_NAME)
                  call PutFile(TMP_NAME,P(1:Mloc,1:Nloc))
            ENDIF

            IF(OUT_Q)THEN
                  TMP_NAME = TRIM(FDIR)//'q_'//TRIM(FILE_NAME)
                  call PutFile(TMP_NAME,Q(1:Mloc,1:Nloc))
            ENDIF


            IF(OUT_AGE)THEN
                  IF(SHOW_BREAKING)THEN
                  TMP_NAME = TRIM(FDIR)//'age_'//TRIM(FILE_NAME)
                  call PutFile(TMP_NAME,AGE_BREAKING)
                  ENDIF
            ENDIF

            IF(OUT_ROLLER)THEN
                  TMP_NAME = TRIM(FDIR)//'roller_'//TRIM(FILE_NAME)
                  call PutFile(TMP_NAME,ROLLER_FLUX)
            ENDIF

            IF(OUT_UNDERTOW)THEN
                  TMP_NAME = TRIM(FDIR)//'U_undertow_'//TRIM(FILE_NAME)
                  call PutFile(TMP_NAME,UNDERTOW_U)
                  TMP_NAME = TRIM(FDIR)//'V_undertow_'//TRIM(FILE_NAME)
                  call PutFile(TMP_NAME,UNDERTOW_V)
            ENDIF

                  IF(VISCOSITY_BREAKING)THEN
                  IF(OUT_NU)THEN
                  TMP_NAME = TRIM(FDIR)//'nubrk_'//TRIM(FILE_NAME)
                  call PutFile(TMP_NAME,nu_break)
                  ENDIF

      ENDIF

       IF(OUT_FrcX)THEN
            TMP_NAME = TRIM(FDIR)//'FrcInsX_'//TRIM(FILE_NAME)
            call PutFile(TMP_NAME,FrcInsX)
       ENDIF
       IF(OUT_FrcY)THEN
            TMP_NAME = TRIM(FDIR)//'FrcInsY_'//TRIM(FILE_NAME)
            call PutFile(TMP_NAME,FrcInsY)
       ENDIF
       IF(OUT_BrkdisX)THEN
            TMP_NAME = TRIM(FDIR)//'BrkSrcX_'//TRIM(FILE_NAME)
            call PutFile(TMP_NAME,BreakSourceX)
       ENDIF
       IF(OUT_BrkdisY)THEN
            TMP_NAME = TRIM(FDIR)//'BrkSrcY_'//TRIM(FILE_NAME)
            call PutFile(TMP_NAME,BreakSourceY)
       ENDIF


      IF(OUT_Time)THEN
            TMP_NAME = TRIM(FDIR)//'time_'//TRIM(FILE_NAME)
            call PutFile(TMP_NAME,ARRTIME)
      ENDIF


101   continue

END SUBROUTINE PREVIEW

!-------------------------------------------------------------------------------------
!
!    PREVIEW_MEAN is subroutine for print-out of mean field data
!
!  HISTORY:
!    03/22/2016  Fengyan Shi
!-------------------------------------------------------------------------------------
SUBROUTINE PREVIEW_MEAN
      USE GLOBAL
      IMPLICIT NONE
      REAL(SP),DIMENSION(Mloc,Nloc) :: tmpout 

      CHARACTER(LEN=80)::FILE_NAME=' '
      CHARACTER(LEN=80)::FDIR=' '
      CHARACTER(LEN=80)::TMP_NAME=' '

      FDIR=TRIM(RESULT_FOLDER)

      ICOUNT_MEAN=ICOUNT_MEAN+1

#       if defined (PARALLEL)
            if (myid.eq.0)then
            WRITE(3,102)'PRINTING MEAN FILE', icount_mean
            WRITE(*,102)'PRINTING MEAN FILE', icount_mean
            endif
#       else
            WRITE(*,102)'PRINTING MEAN FILE', icount_mean
#       endif

      102     FORMAT(A20,I6)

            itmp1=mod(icount_mean/10000,10)
            itmp2=mod(icount_mean/1000,10)
            itmp3=mod(icount_mean/100,10)
            itmp4=mod(icount_mean/10,10)
            itmp5=mod(icount_mean,10)

            write(file_name(1:1),'(I1)')itmp1
            write(file_name(2:2),'(I1)')itmp2
            write(file_name(3:3),'(I1)')itmp3
            write(file_name(4:4),'(I1)')itmp4
            write(file_name(5:5),'(I1)')itmp5  

      IF(OUT_Umean)THEN
            TMP_NAME = TRIM(FDIR)//'umean_'//TRIM(FILE_NAME)
            call PutFile(TMP_NAME,Umean)
            tmpout = P_mean / Max(Depth+ETAmean,MinDepthFrc)
            TMP_NAME = TRIM(FDIR)//'ulagm_'//TRIM(FILE_NAME)
            call PutFile(TMP_NAME,tmpout)
      ENDIF
      
      IF(OUT_Vmean)THEN
            TMP_NAME = TRIM(FDIR)//'vmean_'//TRIM(FILE_NAME)
            call PutFile(TMP_NAME,Vmean)
            tmpout = Q_mean / Max(Depth+ETAmean,MinDepthFrc)
            TMP_NAME = TRIM(FDIR)//'vlagm_'//TRIM(FILE_NAME)
            call PutFile(TMP_NAME,tmpout)
      ENDIF

      IF(OUT_ETAmean)THEN
          TMP_NAME = TRIM(FDIR)//'etamean_'//TRIM(FILE_NAME)
          call PutFile(TMP_NAME,ETAmean)
      ENDIF

      IF(OUT_WaveHeight)THEN
            TMP_NAME = TRIM(FDIR)//'Hrms_'//TRIM(FILE_NAME)
            call PutFile(TMP_NAME,WaveHeightRMS)
            TMP_NAME = TRIM(FDIR)//'Havg_'//TRIM(FILE_NAME)
            call PutFile(TMP_NAME,WaveHeightAve)
            TMP_NAME = TRIM(FDIR)//'Hsig_'//TRIM(FILE_NAME)
            call PutFile(TMP_NAME,SigWaveHeight)
      ENDIF


      IF(OUT_Radiation)THEN
            TMP_NAME = TRIM(FDIR)//'Sxx_'//TRIM(FILE_NAME)
            tmpout = UUmean-WWmean+0.5*9.8*ETA2mean
            call PutFile(TMP_NAME,tmpout)
            TMP_NAME = TRIM(FDIR)//'Sxy_'//TRIM(FILE_NAME)
            call PutFile(TMP_NAME,UVmean)
            TMP_NAME = TRIM(FDIR)//'Syy_'//TRIM(FILE_NAME)
            tmpout = VVmean-WWmean+0.5*9.8*ETA2mean
            call PutFile(TMP_NAME,tmpout)

            TMP_NAME = TRIM(FDIR)//'DxSxx_'//TRIM(FILE_NAME)
            call PutFile(TMP_NAME,DxSxx)
            TMP_NAME = TRIM(FDIR)//'DySxy_'//TRIM(FILE_NAME)
            call PutFile(TMP_NAME,DySxy)
            TMP_NAME = TRIM(FDIR)//'DySyy_'//TRIM(FILE_NAME)
            call PutFile(TMP_NAME,DySyy)
            TMP_NAME = TRIM(FDIR)//'DxSxy_'//TRIM(FILE_NAME)
            call PutFile(TMP_NAME,DxSxy)
            TMP_NAME = TRIM(FDIR)//'PgrdX_'//TRIM(FILE_NAME)
            call PutFile(TMP_NAME,PgrdX)
            TMP_NAME = TRIM(FDIR)//'PgrdY_'//TRIM(FILE_NAME)
            call PutFile(TMP_NAME,PgrdY)
            TMP_NAME = TRIM(FDIR)//'DxUUH_'//TRIM(FILE_NAME)
            call PutFile(TMP_NAME,DxUUH)
            TMP_NAME = TRIM(FDIR)//'DyUVH_'//TRIM(FILE_NAME)
            call PutFile(TMP_NAME,DyUVH)
            TMP_NAME = TRIM(FDIR)//'DyVVH_'//TRIM(FILE_NAME)
            call PutFile(TMP_NAME,DyVVH)
            TMP_NAME = TRIM(FDIR)//'DxUVH_'//TRIM(FILE_NAME)
            call PutFile(TMP_NAME,DxUVH)
            TMP_NAME = TRIM(FDIR)//'FRCX_'//TRIM(FILE_NAME)
            call PutFile(TMP_NAME,FRCXmean)
            TMP_NAME = TRIM(FDIR)//'FRCY_'//TRIM(FILE_NAME)
            call PutFile(TMP_NAME,FRCYmean)
            TMP_NAME = TRIM(FDIR)//'BrkDissX_'//TRIM(FILE_NAME)
            call PutFile(TMP_NAME,BreakDissX)
            TMP_NAME = TRIM(FDIR)//'BrkDissY_'//TRIM(FILE_NAME)
            call PutFile(TMP_NAME,BreakDissY)

        ENDIF
                    

END SUBROUTINE PREVIEW_MEAN


# if defined (PARALLEL)
!-------------------------------------------------------------------------------------
!
!    GetFile is subroutine for reading field data
!
!    HISTORY:
!    05/01/2010  Fengyan Shi
!    05/08/2017  Young-Kwang Choi
!-------------------------------------------------------------------------------------

SUBROUTINE GetFile(FILE,PHI)
      USE GLOBAL
      IMPLICIT NONE

      REAL(SP),DIMENSION(MGlob+2*Nghost,NGlob+2*Nghost) :: PHIGLOB
      CHARACTER(LEN=80) FILE
      REAL(SP),DIMENSION(Mloc,Nloc),INTENT(OUT) :: PHI

![-------ykchoi (08/May/2017)
      INTEGER :: irank, lenx, leny, lenxy, ireq
      INTEGER :: Nista, Niend, Njsta, Njend
      INTEGER :: istanum, iendnum, jstanum, jendnum
      INTEGER, ALLOCATABLE :: Nistas(:), Niends(:), Njstas(:), Njends(:)
      INTEGER :: istatus(mpi_status_size)
      REAL(SP), ALLOCATABLE :: xx(:,:)
! -------ykchoi (08/May/2017) ]

! TEMP

      if (myid.eq.0) then
            OPEN(1,FILE=TRIM(FILE))
            DO J=Nghost+1,NGlob+NGhost
            READ(1,*)(PHIGLOB(I,J),I=Nghost+1,MGlob+Nghost)
            ENDDO
            CLOSE(1)
      ! ghost cells
            DO I=Nghost+1,MGlob+Nghost
            DO J=1,Nghost
                  PHIGLOB(I,J)=PHIGLOB(I,Nghost+1)
            ENDDO
            DO J=NGlob+Nghost+1,NGlob+2*Nghost
                  PHIGLOB(I,J)=PHIGLOB(I,NGlob+Nghost)
            ENDDO
            ENDDO
            DO J=1,NGlob+2*Nghost
            DO I=1,Nghost
                  PHIGLOB(I,J)=PHIGLOB(Nghost+1,J)
            ENDDO
            DO I=MGlob+Nghost+1,MGlob+2*Nghost
                  PHIGLOB(I,J)=PHIGLOB(MGlob+Nghost,J)
            ENDDO
            ENDDO
      endif

      ![-------ykchoi (08/May/2017)
      Nista = iista + Nghost;
      Niend = iiend + Nghost;
      Njsta = jjsta + Nghost;
      Njend = jjend + Nghost;

      allocate( Nistas(nprocs), Niends(nprocs), Njstas(nprocs), Njends(nprocs) )

      call MPI_Gather( Nista, 1, MPI_INTEGER, Nistas, 1, MPI_INTEGER, &
                        0, MPI_COMM_WORLD, ier )
      call MPI_Gather( Niend, 1, MPI_INTEGER, Niends, 1, MPI_INTEGER, &
                        0, MPI_COMM_WORLD, ier )
      call MPI_Gather( Njsta, 1, MPI_INTEGER, Njstas, 1, MPI_INTEGER, &
                        0, MPI_COMM_WORLD, ier )
      call MPI_Gather( Njend, 1, MPI_INTEGER, Njends, 1, MPI_INTEGER, &
                        0, MPI_COMM_WORLD, ier )

      if( myid == 0 )then
            PHI = PHIGLOB( 1:Mloc, 1:Nloc )
      endif

      do irank=1, px*py-1
            if( myid == 0 ) then
            istanum = Nistas(irank+1) - Nghost
            iendnum = Niends(irank+1) + Nghost
            jstanum = Njstas(irank+1) - Nghost
            jendnum = Njends(irank+1) + Nghost

            lenx = iendnum - istanum + 1
            leny = jendnum - jstanum + 1
            lenxy = lenx*leny
            allocate( xx(lenx, leny) )

            xx = PHIGLOB( istanum:iendnum, jstanum:jendnum )
            call mpi_isend( xx, lenxy, mpi_sp, irank, 1, mpi_comm_world, ireq, ier )
            call mpi_wait( ireq, istatus, ier )
            deallocate( xx )

            elseif( myid == irank ) then
            
            lenx = Niend-Nista+1+2*Nghost
            leny = Njend-Njsta+1+2*Nghost
            lenxy = lenx*leny

            call mpi_irecv( PHI, lenxy, mpi_sp, 0, 1, mpi_comm_world, ireq, ier )
            call mpi_wait( ireq, istatus, ier )

            endif
      enddo

      deallocate( Nistas, Niends, Njstas, Njends )

! -------ykchoi (08/May/2017) ]

END SUBROUTINE Getfile

# endif

# if defined (PARALLEL)
!-------------------------------------------------------------------------------------
!
!    PutFile is subroutine for print-out of field data
!
!    HISTORY:
!      05/01/2010  Fengyan Shi
!      05/06/2017  Young-Kwang Choi 
!-------------------------------------------------------------------------------------

SUBROUTINE PutFile(FILE_NAME,PHI)
      USE GLOBAL
      USE PARALLEL_FIELD_IO
      IMPLICIT NONE

      CHARACTER(LEN=80) FILE_NAME
      REAL(SP),DIMENSION(Mloc,Nloc),INTENT(IN) :: PHI

      CHARACTER(LEN=80)::TMP_NAME=' '

      SELECT CASE (TRIM(FIELD_IO_TYPE))
      CASE ('ASCII' , 'ascii')
            CALL PutFileASCII(FILE_NAME,PHI)
      CASE ('BINARY' , 'binary' )
            Call PutFileBinary(FILE_NAME,PHI)
      CASE DEFAULT
            !Defaults to ASCII case for non-valid input
            CALL PutFileASCII(FILE_NAME,PHI)
      END SELECT

END SUBROUTINE Putfile

# else
!-------------------------------------------------------------------------------------
!
!    PutFile is subroutine for print-out of field data
!
!    HISTORY:
!      05/01/2010  Fengyan Shi
!
!-------------------------------------------------------------------------------------
SUBROUTINE PutFile(FILE,PHI)
      USE PARAM
      USE GLOBAL
      IMPLICIT NONE
      REAL(SP),DIMENSION(Mloc,Nloc),INTENT(IN) :: PHI
      CHARACTER(LEN=80) FILE
      LOGICAL :: FirstCallPutFile = .TRUE.
      SAVE  FirstCallPutFile

! first time call 
      IF(FirstCallPutFile)THEN
            FirstCallPutFile = .FALSE.
      ! format length
            write(FORMAT_LEN(1:1),'(A1)') '('
            write(FORMAT_LEN(2:8),'(I7)') Mglob
            write(FORMAT_LEN(9:15),'(A7)') 'E16.6E4'
            write(FORMAT_LEN(16:16),'(A1)') ')'
      ENDIF

            OPEN(1,FILE=TRIM(FILE))
#       if defined(DEBUG)
            DO J=1,Nloc
            WRITE(1,100)(real(PHI(I,J)),I=1,Mloc)
            ENDDO
#       else
            DO J=Nghost+1,Nloc-Nghost,OUTPUT_RES
            WRITE(1,FORMAT_LEN)(real(PHI(I,J)),I=Nghost+1,Mloc-Nghost,OUTPUT_RES)
            ENDDO
#       endif
      100  FORMAT(5000E16.6)
      !100   FORMAT(FORMAT_LEN)
            CLOSE(1)
END SUBROUTINE PutFile

# endif







