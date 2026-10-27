MODULE History_Write_Fixture
  USE History_Fixture_Types
  IMPLICIT NONE
CONTAINS
  ! Include the complete production routine, not a copy of its control flow.
  INCLUDE 'history_write_under_test.inc'

  SUBROUTINE History_Update(input,state_diag,rc)
    TYPE(OptInput), INTENT(IN) :: input
    TYPE(DgnState), INTENT(IN) :: state_diag
    INTEGER, INTENT(OUT) :: rc
    update_calls=update_calls+1
    rc=0
  END SUBROUTINE
END MODULE

PROGRAM History_Startup_Test
  USE History_Write_Fixture
  IMPLICIT NONE
  TYPE(OptInput) :: input
  TYPE(ChmState) :: chemistry
  TYPE(DgnState) :: diagnostics
  TYPE(GrdState) :: grid
  TYPE(DgnMap), TARGET :: mapping
  TYPE(HistContainer), TARGET :: containers(4)
  TYPE(MetaHistContainer), TARGET :: nodes(4)
  CHARACTER(LEN=32), PARAMETER :: names(4)=[CHARACTER(LEN=32) :: &
       'BoundaryConditions','SolarCache','Processes','Restart']
  INTEGER :: a,b,c,d,i,rc,permutations,assertions

  ALLOCATE(chemistry%Species(2))
  DO i=1,2
    ALLOCATE(chemistry%Species(i)%Conc(1,1,1))
  ENDDO
  chemistry%Species(1)%Conc=1.25_fp
  chemistry%Species(2)%Conc=3.5_fp
  mapping%nSlots=2
  mapping%slot2id=[2,1]
  diagnostics%Map_SpeciesBC=>mapping
  ALLOCATE(diagnostics%SpeciesBC(1,1,1,2))
  assertions=0
  permutations=0

  ! Exercise every traversal order, including BC first, middle and last.
  DO a=1,4
    DO b=1,4
      IF(b==a) CYCLE
      DO c=1,4
        IF(c==a.OR.c==b) CYCLE
        DO d=1,4
          IF(d==a.OR.d==b.OR.d==c) CYCLE
          CALL setup([a,b,c,d])
          CALL History_Write(input,chemistry,diagnostics,grid,rc)
          CALL expect(rc==0,'startup succeeds')
          CALL expect(containers(1)%writes==1,'exactly one startup BC write')
          CALL expect(containers(1)%defines==1,'BC file defined at startup')
          CALL expect(ALL(containers(1)%captured==[3.5_fp,1.25_fp]),'mapped initial state captured')
          CALL expect(containers(1)%FileWriteAlarm==3600._fp,'BC next alarm advanced')
          CALL expect(ALL(containers(2:4)%writes==0),'no other startup writes')
          CALL expect(ALL(containers(2:4)%defines==0),'no other startup file definitions')
          CALL expect(ALL(containers(2:4)%FileCloseAlarm==0._fp),'other startup alarms untouched')
          CALL expect(update_calls==1,'one startup update call')
          CALL expect(chemistry%Species(1)%Conc(1,1,1)==1.25_fp,'first species unchanged')
          CALL expect(chemistry%Species(2)%Conc(1,1,1)==3.5_fp,'second species unchanged')

          ! Normal writes still visit every collection and advance alarms.
          containers%ElapsedSec=3600._fp
          containers%writes=0
          containers%defines=0
          containers%closes=0
          update_calls=0
          CALL History_Write(input,chemistry,diagnostics,grid,rc)
          CALL expect(rc==0,'post-start succeeds')
          CALL expect(ALL(containers%writes==1),'all due collections written')
          CALL expect(ALL(containers%defines==1),'all due collections defined')
          CALL expect(ALL(containers%FileWriteAlarm==7200._fp),'write alarms advanced')
          CALL expect(update_calls==0,'no extra post-start update call')
          CALL expect(containers(4)%closes==2,'restart closes after write')
          permutations=permutations+1
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  CALL expect(permutations==24,'all permutations exercised')

  CALL setup([2,3,4])
  CALL History_Write(input,chemistry,diagnostics,grid,rc)
  CALL expect(rc==0.AND.ALL(containers%writes==0),'absent BC startup is a no-op')
  CALL expect(update_calls==0,'absent BC does not update')

  CALL setup([INTEGER ::])
  CALL History_Write(input,chemistry,diagnostics,grid,rc)
  CALL expect(rc==0.AND.update_calls==0,'empty list returns cleanly')

  CALL setup([1])
  CALL History_Write(input,chemistry,diagnostics,grid,rc)
  CALL expect(rc==0.AND.containers(1)%writes==1,'BC-only list captures startup')

  CALL setup([2,1,3,4])
  containers%ElapsedSec=600._fp
  containers%FileCloseAlarm=3600._fp
  CALL History_Write(input,chemistry,diagnostics,grid,rc)
  CALL expect(rc==0.AND.ALL(containers%writes==0),'not-due post-start collections skipped')
  CALL expect(ALL(containers%defines==0),'not-due files not reopened')

  CALL setup([2,1,3,4])
  chemistry%Species(2)%Units=2
  CALL History_Write(input,chemistry,diagnostics,grid,rc)
  CALL expect(rc/=0,'invalid startup species units still rejected')
  CALL expect(ALL(containers%writes==0),'error prevents startup writes')

  WRITE(*,'(a,i0,a,i0,a)') 'PASS: ',permutations,' collection orders, ',assertions,' assertions'
CONTAINS
  SUBROUTINE setup(order)
    INTEGER, INTENT(IN) :: order(:)
    INTEGER :: n,j
    containers=HistContainer()
    DO j=1,4
      containers(j)%Name=names(j)
      nodes(j)%Container=>containers(j)
      NULLIFY(nodes(j)%Next)
    ENDDO
    NULLIFY(CollectionList)
    n=SIZE(order)
    IF(n>0) THEN
      CollectionList=>nodes(order(1))
      DO j=1,n-1
        nodes(order(j))%Next=>nodes(order(j+1))
      ENDDO
    ENDIF
    diagnostics%SpeciesBC=-1._fp
    update_calls=0
  END SUBROUTINE
  SUBROUTINE expect(ok,message)
    LOGICAL, INTENT(IN) :: ok
    CHARACTER(*), INTENT(IN) :: message
    assertions=assertions+1
    IF(.NOT.ok) THEN
      WRITE(*,'(a)') 'FAIL: '//message
      STOP 1
    ENDIF
  END SUBROUTINE
END PROGRAM
