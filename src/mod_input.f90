MODULE INPUT_READ

    USE PARAM,ONLY : SP
    #if defined (PARALLEL)
      USE GLOBAL,ONLY : myid,ier
      USE MPI
    # endif
      IMPLICIT NONE
      LOGICAL :: FILE_EXIST
      CHARACTER(LEN=80) :: VAL_TYPE 
      CHARACTER(LEN=80) :: VAL_READ 
    
    INTERFACE READ_FOUR_TYPE_VALUES
      Module Procedure READ_FLOAT
      Module Procedure READ_STRING
      Module Procedure READ_LOGICAL
      Module Procedure READ_INTEGER
    END INTERFACE
    
    CONTAINS
      
    ! INTEGER
    
    SUBROUTINE READ_INTEGER(INTEGER_VAL,InputFileName,MatchCharacter,ERR_IF_ONE)  
      IMPLICIT NONE
      INTEGER,          INTENT(OUT) :: INTEGER_VAL 
      CHARACTER(LEN=*), INTENT(IN ) :: InputFileName
      CHARACTER(LEN=*), INTENT(IN ) :: MatchCharacter
      INTEGER           :: ERR_IF_ONE
    
      INQUIRE(FILE=InputFileName,EXIST=FILE_EXIST)
      IF(.NOT.FILE_EXIST)THEN
    # if defined (PARALLEL)
       IF(MYID==0) WRITE(*,*) TRIM(InputFileName), ' CAN NOT BE FOUND. STOP'
       CALL MPI_FINALIZE (ier)
       STOP
    # else
        WRITE(*,*) TRIM(InputFileName), 'CAN NOT BE FOUND. STOP'
        STOP
    # endif
      ENDIF
    
      CALL READ_VALUE(InputFileName,MatchCharacter,VAL_READ,ERR_IF_ONE)
      IF(ERR_IF_ONE == 1)THEN
    # if defined (PARALLEL)
       IF(MYID==0) WRITE(*,*) TRIM(MatchCharacter), ' DOES NOT EXIST. USE DEFAULT VALUE'
    # else
        WRITE(*,*) TRIM(MatchCharacter), ' DOES NOT EXIST. USE DEFAULT VALUE'
    # endif
         RETURN
      ENDIF
    
      !check argument type
      CALL CHECK_CONSISTENCY_TYPE(VAL_READ,VAL_TYPE)
      IF(TRIM(VAL_TYPE) .NE. 'INTEGER')THEN
    # if defined (PARALLEL)
       IF(MYID==0) WRITE(*,*) TRIM(MatchCharacter), ' IS NOT INTEGER, STOP!'
       CALL MPI_FINALIZE (ier)
       STOP
    # else
        WRITE(*,*) TRIM(MatchCharacter), ' IS NOT INTEGER, STOP!'
        STOP
    # endif
      ENDIF
    
      READ(VAL_READ,*)INTEGER_VAL   
    
      
      END SUBROUTINE READ_INTEGER
    
    
    ! Float
    
      SUBROUTINE READ_FLOAT(FLOAT_VAL,InputFileName,MatchCharacter,ERR_IF_ONE)  
      IMPLICIT NONE
      REAL(SP),         INTENT(OUT) :: FLOAT_VAL
      CHARACTER(LEN=*), INTENT(IN ) :: InputFileName
      CHARACTER(LEN=*), INTENT(IN ) :: MatchCharacter
    
      INTEGER           :: ERR_IF_ONE
    
      INQUIRE(FILE=InputFileName,EXIST=FILE_EXIST)
      IF(.NOT.FILE_EXIST)THEN
    # if defined (PARALLEL)
       IF(MYID==0) WRITE(*,*) TRIM(InputFileName), ' CAN NOT BE FOUND. STOP'
       CALL MPI_FINALIZE (ier)
       STOP
    # else
        WRITE(*,*) TRIM(InputFileName), 'CAN NOT BE FOUND. STOP'
        STOP
    # endif
      ENDIF
    
      CALL READ_VALUE(InputFileName,MatchCharacter,VAL_READ,ERR_IF_ONE)
      IF(ERR_IF_ONE == 1)THEN
    # if defined (PARALLEL)
       IF(MYID==0) WRITE(*,*) TRIM(MatchCharacter), ' DOES NOT EXIST. USE DEFAULT VALUE'
    # else
        WRITE(*,*) TRIM(MatchCharacter), ' DOES NOT EXIST. USE DEFAULT VALUE'
    # endif
         RETURN
      ENDIF
    
      CALL CHECK_CONSISTENCY_TYPE(VAL_READ,VAL_TYPE)
      IF(TRIM(VAL_TYPE) .NE. 'float')THEN
    # if defined (PARALLEL)
       IF(MYID==0) WRITE(*,*) TRIM(MatchCharacter), ' IS NOT FLOAT, STOP!'
       CALL MPI_FINALIZE (ier)
       STOP
    # else
        WRITE(*,*) TRIM(MatchCharacter), ' IS NOT FLOAT, STOP!'
        STOP
    # endif
      ENDIF
    
      READ(VAL_READ,*)FLOAT_VAL  
        
      
      End SUBROUTINE READ_FLOAT
    
    
    ! Logical
      SUBROUTINE READ_LOGICAL(LOGICAL_VAL,InputFileName,MatchCharacter,ERR_IF_ONE)  
      IMPLICIT NONE
      LOGICAL,          INTENT(OUT) :: LOGICAL_VAL
      CHARACTER(LEN=*), INTENT(IN ) :: InputFileName
      CHARACTER(LEN=*), INTENT(IN ) :: MatchCharacter
    
      INTEGER           :: ERR_IF_ONE
    
      INQUIRE(FILE=InputFileName,EXIST=FILE_EXIST)
      IF(.NOT.FILE_EXIST)THEN
    # if defined (PARALLEL)
       IF(MYID==0) WRITE(*,*) TRIM(InputFileName), ' CAN NOT BE FOUND. STOP'
       CALL MPI_FINALIZE (ier)
       STOP
    # else
        WRITE(*,*) TRIM(InputFileName), 'CAN NOT BE FOUND. STOP'
        STOP
    # endif
      ENDIF
    
      call READ_VALUE(InputFileName,MatchCharacter,VAL_READ,ERR_IF_ONE)
      IF(ERR_IF_ONE == 1)THEN
    # if defined (PARALLEL)
       IF(MYID==0) WRITE(*,*) TRIM(MatchCharacter), ' DOES NOT EXIST. USE DEFAULT VALUE'
    # else
        WRITE(*,*) TRIM(MatchCharacter), ' DOES NOT EXIST. USE DEFAULT VALUE'
    # endif
         RETURN
      ENDIF
      !check argument type
      CALL CHECK_CONSISTENCY_TYPE(VAL_READ,VAL_TYPE)
      IF(TRIM(VAL_TYPE) .NE. 'logical')THEN
    # if defined (PARALLEL)
       IF(MYID==0) WRITE(*,*) TRIM(MatchCharacter), ' IS NOT LOGICAL, SET FALSE EXCEPT FIRST LETTER IS T'
    # else
        WRITE(*,*) TRIM(MatchCharacter), ' IS NOT LOGICAL, SET FALSE EXCEPT FIRST LETTER IS T'
    # endif
      ENDIF
    
      LOGICAL_VAL = .FALSE.
      IF(VAL_READ(1:1) == 'T') LOGICAL_VAL = .TRUE.
    
      End SUBROUTINE READ_LOGICAL
    
    
      SUBROUTINE READ_STRING(STRING_VAL,InputFileName,MatchCharacter,ERR_IF_ONE)  
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(OUT) :: STRING_VAL
      CHARACTER(LEN=*), INTENT(IN ) :: InputFileName
      CHARACTER(LEN=*), INTENT(IN ) :: MatchCharacter
      INTEGER, INTENT(OUT) :: ERR_IF_ONE
    
      INQUIRE(FILE=InputFileName,EXIST=FILE_EXIST)
      IF(.NOT.FILE_EXIST)THEN
    # if defined (PARALLEL)
       IF(MYID==0) WRITE(*,*) TRIM(InputFileName), ' CAN NOT BE FOUND. STOP'
       CALL MPI_FINALIZE (ier)
       STOP
    # else
        WRITE(*,*) TRIM(InputFileName), 'CAN NOT BE FOUND. STOP'
        STOP
    # endif
      ENDIF
    
    
      call READ_VALUE(InputFileName,MatchCharacter,VAL_READ,ERR_IF_ONE)
      IF(ERR_IF_ONE == 1)THEN
    # if defined (PARALLEL)
       IF(MYID==0) WRITE(*,*) TRIM(MatchCharacter), ' DOES NOT EXIST. USE DEFAULT VALUE'
    # else
        WRITE(*,*) TRIM(MatchCharacter), ' DOES NOT EXIST. USE DEFAULT VALUE'
    # endif
         RETURN
      ENDIF
    
      !check argument type
      call CHECK_CONSISTENCY_TYPE(VAL_READ,VAL_TYPE)
      IF(TRIM(VAL_TYPE) .NE. 'string')THEN
    # if defined (PARALLEL)
       IF(MYID==0) WRITE(*,*) TRIM(MatchCharacter), ' IS NOT STRING, STOP'
       CALL MPI_FINALIZE (ier)
       STOP
    # else
        WRITE(*,*) TRIM(MatchCharacter), ' IS NOT STRING, STOP'
        STOP
    # endif
      ENDIF
    
      STRING_VAL = ADJUSTL(TRIM(VAL_READ))
      
    
      End SUBROUTINE READ_STRING
    
    
    !  read value from file 
    
      SUBROUTINE READ_VALUE(FileName,MatchCharacter,VAL,ERR_IF_ONE)
      IMPLICIT NONE
      CHARACTER(LEN=* ), INTENT(IN ) :: FileName 
      CHARACTER(LEN=* ), INTENT(IN ) :: MatchCharacter 
      CHARACTER(LEN=80), INTENT(OUT) :: VAL
      INTEGER,           INTENT(OUT) :: ERR_IF_ONE
      
      CHARACTER(LEN=80)  :: TMPV,TMPV_NEXT,TEXT,TEXT_READ
      INTEGER            :: LOCATION,FRONT,TEXT_LENGTH,I
      CHARACTER(LEN=1), PARAMETER :: COM_CHAR = '!'
    
      ERR_IF_ONE   = 0
      OPEN(10,FILE=TRIM(FileName)) ; REWIND(10) 
    
       DO WHILE(.TRUE.)
    
         READ(10,'(A)',END=100) TEXT_READ 
    
    
        LOCATION = index(TEXT_READ,com_char)
        if(LOCATION == 0)LOCATION = len_trim(TEXT_READ)+1
        TMPV = TEXT_READ(1:LOCATION-1)
        TEXT = TMPV
        TEXT_LENGTH = len_trim(TEXT)
        FRONT = index(TEXT,"=")
        TMPV_NEXT = TEXT(1:FRONT-1)
        IF(ADJUSTL(TRIM(TMPV_NEXT)) == TRIM(MatchCharacter))THEN
          VAL  = ADJUSTL(TEXT(FRONT+1:TEXT_LENGTH))
          RETURN
        ENDIF
    
      ENDDO
    
      100 CONTINUE
      ERR_IF_ONE = 1
      CLOSE(10)
    
      End SUBROUTINE READ_VALUE  
    
      SUBROUTINE CHECK_CONSISTENCY_TYPE(VREAD,VTYPE)
      IMPLICIT NONE
      CHARACTER(LEN=80), INTENT(OUT) :: VTYPE 
      CHARACTER(LEN=80), INTENT(IN ) :: VREAD  
    
      CHARACTER(LEN=16) :: NUMBERS
      INTEGER           :: LOGICAL_VAL_READ
      INTEGER           :: DotLocation
      INTEGER           :: I
      numbers = "0123456789+-Ee. " 
    
      LOGICAL_VAL_READ = LEN_TRIM(ADJUSTL(VREAD))
    
    ! check it it is logical
      IF(LOGICAL_VAL_READ == 1)THEN
      IF(VREAD(1:1) == "T" .or. VREAD(1:1) == "F")then 
        VTYPE = 'logical' 
        RETURN
      ENDIF
      ENDIF
        
      ! check if it is a string
      DO I=1,LOGICAL_VAL_READ
        IF(INDEX(NUMBERS,VREAD(i:i)) == 0)THEN
          VTYPE = "string" 
          RETURN
        ENDIF
      ENDDO
    
      ! check if it is a float
      DotLocation = index(VREAD,".")
      IF(DotLocation .NE. 0) then
        VTYPE = "float"
      ELSE
        VTYPE = "INTEGER"
      ENDIF
    
      End SUBROUTINE CHECK_CONSISTENCY_TYPE
    
    END MODULE INPUT_READ
    
    
    