!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: aqueous_mod
!
! !DESCRIPTION: Module Aqueous\_Mod contines arrays and routines for the
!  aqueous chemistry scheme.
!\\
!\\
! !INTERFACE: 
!
MODULE Aqueous_Mod
!
! !USES:
!
  USE Precision_Mod            ! For GEOS-Chem Precision (fp)

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Partition_Soluble
!    
! !REVISION HISTORY:
!  01 Jan 2016 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: partition_soluble
!
! !DESCRIPTION: Subroutine Partition\_Soluble calculates partitioning of species
!  between gas and aqueous phase. Concentrations in liquid cloud water, ice
!  cloud water and (when relevant) aerosol are calculated together.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Partition_Soluble( am_I_Root, Input_Opt, State_Met, State_Chm, RC  )
!
! !USES:
!
    USE CHEMGRID_MOD,         ONLY : ITS_IN_THE_CHEMGRID, ITS_IN_THE_TROP
#if !defined( NO_ISORROPIA )
    Use IsoropiaII_Mod,       Only : Get_ISRInfo
#endif

    USE CMN_SIZE_MOD,         ONLY : IIPAR, JJPAR, LLPAR
    USE ErrCode_Mod
    USE ERROR_MOD
    USE PhysConstants,        ONLY : AirMW
    USE State_Chm_Mod,        ONLY : ChmState
    USE State_Chm_Mod,        ONLY : Ind_
    USE State_Met_Mod,        ONLY : MetState
    USE Henry_Mod,            ONLY : Calc_KH
    USE Henry_Mod,            ONLY : Calc_Heff
    USE Input_Opt_Mod,        ONLY : OptInput
    USE Species_Mod,          ONLY : Species
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root  ! Is this the root CPU?
    TYPE(OptInput), INTENT(IN)    :: Input_Opt  ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetState), INTENT(INOUT) :: State_Met  ! Meteorology State object
    TYPE(ChmState), INTENT(INOUT) :: State_Chm  ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC         ! Success or failure
! 
! !REVISION HISTORY:
!  01 Dec 2016 - S. D. Eastham - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                :: F, M, SpcID
    LOGICAL                :: prtDebug

    ! Species indexing (constants)
    Integer                :: id_HBr,  id_HCl,  id_HI
    Integer                :: id_BrAq, id_ClAq, id_IAq

    ! Species indexing (variable)
    Integer                :: id_Other,  id_Fine
    Integer                :: X

    ! Species properties
    Real(fp)               :: pKa,     CR
    Real(fp)               :: KH,      K0

    ! Aerosol volume calculation
    Real(fp)               :: QL,      QI,      T2L,      DR_Ratio
    Real(fp)               :: VFine,   VLCld,   VICld
    Real(fp)               :: VSulf,   VSSA

    ! Aerosol H+ calculation
    Real(fp)               :: pHFine,  pHLCld,  pHICld
    Real(fp)               :: HEffFine,HEffLCld,HEffICld
    Real(fp)               :: rFine,   rLCld,   rICld

    ! Aerosol water content (mol H2O/m3 air)
    Real(fp)               :: H2OFine

    ! Concentrations in molec/cm3
    Real(fp)               :: CGas,    CLiq,    CFine,    COther
    Real(fp)               :: CICld,   CLCld,   CTot

    ! Fractions in each phase/state
    Real(fp)               :: xGas,    xLiq,    xFine,    xOther
    Real(fp)               :: xICld,   xLCld

    ! Temporary variables
    Real(fp)               :: xDenom
    Integer                :: RCLocal

    ! Grid box properties
    Real(fp)               :: TK

    ! Index variables
    Integer                :: I, J, L, N

    ! Flag to calculate partitioning but not overwrite data
    Logical                :: lockCFine

    ! Flag to include in-aerosol quantity for partitioning
    Logical                :: useCFine

    ! Flag to calculate partitioning into aerosol phase by assuming equilibrium
    Logical                :: calcHFine

    ! Objects
    TYPE(Species), POINTER :: SpcInfo

    !=================================================================
    ! Partition_Soluble begins here!
    !=================================================================

    ! Assume success
    RC = GC_SUCCESS

    ! Initialize pointer
    SpcInfo => NULL()

    ! Partition soluble halogen species between liquid and gas phase throughout
    ! the troposphere (SDE 2016-12-30)
    id_HBr  = ind_('HBr' )
    id_HCl  = ind_('HCl' )
    id_HI   = ind_('HI'  )
    id_BrAq = ind_('BrAq')
    id_ClAq = ind_('ClAq')
    id_IAq  = ind_('IAq' )

    ! Setup for OpenMP
    ! Note: There is a not expected to be significant variation in loop run
    ! time, so an experiment should be run to see if STATIC scheduling makes
    ! much sense. Note that OpenMP will only parallelize over the first loop, in
    ! this case I
    !$OMP PARALLEL DO                                           &
    !$OMP DEFAULT  ( SHARED                                   ) &
    !$OMP PRIVATE  ( I,        J,        L,        TK         ) &
    !$OMP PRIVATE  ( VLCld,    VICld,    VSulf,    VSSA       ) &
    !$OMP PRIVATE  ( QL,       QI                             ) &
    !$OMP PRIVATE  ( xGas,     xLCld,    xFine,    xICld      ) &
    !$OMP PRIVATE  (           rLCld,    rFine,    rICld      ) &
    !$OMP PRIVATE  ( xOther,   CFine,    COther,   CTot       ) &
    !$OMP PRIVATE  ( CR,       K0,       pKa,      KH         ) &
    !$OMP PRIVATE  ( HEffLCld, HEffICld, HEffFine             ) &
    !$OMP PRIVATE  ( pHLCld,   pHICld,   pHFine,   H2OFine    ) &
    !$OMP PRIVATE  ( X,        T2L,      DR_Ratio, xDenom     ) &
    !$OMP PRIVATE  ( id_Other, id_Fine,  SpcInfo,  RCLocal    ) &
    !$OMP PRIVATE  ( useCFine, lockCFine,calcHFine            ) &
    !$OMP SCHEDULE ( DYNAMIC                                  )
    Do I=1,IIPar
    Do J=1,JJPar
    Do L=1,LLPar

       ! If not in the troposphere, return all material dissolved into "clouds"
       ! into the gas phase but leave that dissolved in fine aerosol where it is
       If ( .not. ITS_IN_THE_TROP(I,J,L,State_Met) ) Then
          Do X=1,3
             ! 1: Gas, 2: liquid cloud, 3: ice cloud
             State_Chm%AqPtn(I,J,L,X,1) = 1.0e+0_fp
             State_Chm%AqPtn(I,J,L,X,2) = 0.0e+0_fp 
             State_Chm%AqPtn(I,J,L,X,3) = 0.0e+0_fp 
          End Do
          Cycle
       End If

       ! Safety
       VLCld = 0.0e+0_fp
       VICld = 0.0e+0_fp
       VSulf = 0.0e+0_fp
       VSSA  = 0.0e+0_fp

       ! Convert temperature to REAL*8?
       !TK_8 = State_Met%T(I,J,L)
       TK = State_Met%T(I,J,L)

       ! Defined QL and QI as m3 H2O/m3 air in-cloud
#if defined( GEOS_5 ) || defined( MERRA) || defined( GEOS_FP) || defined( MERRA2 )
       ! In GEOS-5, QL is [kg H2O/kg air]; convert to m3/m3, assuming a density
       ! of 1000 kg/m3 for water
       QL = State_Met%QL(I,J,L) * State_Met%AirDen(I,J,L) * 1.0e-3_fp 
       QI = State_Met%QI(I,J,L) * State_Met%AirDen(I,J,L) * 1.0e-3_fp 
#else
       ! Otherwise, QL is [m3 H2O/m3 air]; no conversion necessary
       QL = State_Met%QL(I,J,L)
       QI = State_Met%QI(I,J,L)
#endif

       ! Convert from in-cloud to grid-box quantities
       QL = State_Met%CLDF(I,J,L) * QL
       QI = State_Met%CLDF(I,J,L) * QI

       ! Get the total volume of cloud water, including a thin surface layer
       ! on cloud ice (2% of the radius)
       DR_RATIO = 2e-2_fp
       T2L = 1.0e0_fp / ( 1.0e0_fp - (1.0e0_fp - DR_RATIO)**3.0e0_fp )
       VLCld = QL 
       VICld = QI/T2L

       ! For sulfate, use the calculated data which takes into account RH
       VSulf = State_Chm%AeroArea(I,J,L,8)*State_Chm%AeroRadi(I,J,L,8)/3.0e+0_fp

       ! Add in the volume of fine sea salt. These are assumed to be mixed for
       ! the purposes of aqueous-phase chemistry, but coarse sea salt is assumed
       ! to be uninvolved
       VSSA  = State_Chm%AeroArea(I,J,L,11)*State_Chm%AeroRadi(I,J,L,11)/3.0e+0_fp

       ! Total volume of fine aerosol
       VFine = VSSA + VSulf

       ! Use calculated liquid cloud pH. Assume 4.5 for ice cloud water.
       pHLCld = State_Chm%phCloud(I,J,L)
       pHICld = 4.5e+0_fp

       If (VFine.gt.1.0e-20_fp) Then
          ! For sulfate, use the pH calculated by ISORROPIA II, if present.
          ! Otherwise, use a pH of 0.0 when dominated by sulfate, and 5.0 when
          ! dominated by fine sea salt (by volume)
          ! NOTE: Also get the aerosol liquid water content (mol/L)
#if defined( NO_ISORROPIA )
          If (VSSA.gt.VSulf) Then
             pHFine = 5.0e+0_fp
          Else
             pHFine = 0.0e+0_fp
          End If
          ! Temporary estimate; assume 25 mol/L aerosol and convert to mol/m3 air
          H2OFine = 25.0e+0_fp * (1.0e+3_fp * VFine)
#else
          ! Retrieve pH and [H2O] from ISORROPIA
          pHFine  = Get_ISRInfo(I,J,L,1)
          ! ISORROPIA uses -999 to indicate no aerosol
          If (pHFine < -100.0e+0_fp) Then
             pHFine = 0.0e+0_fp
             H2OFine = 0.0e+0_fp
          Else
             ! ISORROPIA returns water in mol/m3 air
             H2OFine = Get_ISRInfo(I,J,L,3)
          End If
#endif
       Else
          ! No fine aerosol
          H2OFine = 0.0e+0_fp
          ! Assume a pH of 0.0
          pHFine = 0.0e+0_fp
       End If

       ! Store the fine aerosol pH
       State_Chm%pHFine(I,J,L) = pHFine

       Do X=1,3
          ! Get the total quantity of the soluble species
          ! Options:
          ! lockCFine: If TRUE, do not allow any changes to fine aerosol conc.
          ! calcHFine: If TRUE, calculate fine aerosol HEff based on K0, pH etc.
          ! useCFine : If TRUE, include the fraction already considered to be
          !            dissolved in fine aerosol when calculating partitioning.
          If (X==1) Then
             ! =====BROMINE AND CHLORINE*======================================
             ! Calculate HEff for fine aerosol assuming Henry's law equilibrium
             ! between gas, cloud and aerosol.
             id_Other  = id_HBr
             id_Fine   = id_BrAq
             lockCFine = .False.
             calcHFine = .True.
             useCFine  = .True.
          ElseIf (X==2) Then
             ! =====CHLORINE===================================================
             ! ISORROPIA already calculated partitioning of Cl between HCl(g)
             ! and aerosol. As such, we won't allow Henry's law to move any  
             ! additional Cl into or out of the fine aerosol.
             ! NOTE: If ISORROPIA is disabled, we partition as usual.
             id_Other = id_HCl
             id_Fine  = id_ClAq
#if !defined( NO_ISSOROPIA )
             ! Only calculate partitioning of non-aerosol Cl-
             ! Enforce the existing partitioning between Cl- in fine aerosol and
             ! all other forms by calculating an effective Henry's law constant
             ! for fine aerosol based on the ISORROPIA results. This allows
             ! competition between fine aerosol and cloud water to be
             ! represented.
             lockCFine = .True.
             calcHFine = .False.
             useCFine  = .True.
#else
             ! Standard treatment, same as bromine
             lockCFine = .False.
             calcHFine = .True.
             useCFine  = .True.
#endif
          ElseIf (X==3) Then
             ! =====IODINE=====================================================
             ! Uptake to fine aerosol is considered to be an irreversible
             ! process, represented through pseudoreaxtion. As such, we remove
             ! the already-dissolved aerosol from consideration and do not
             ! allow equilibrium partitioning of HI into fine aerosol.
             ! Instead, act as if the only iodine present is that in the HI
             ! tracer.
             id_Other  = id_HI
             id_Fine   = id_IAq
             lockCFine = .True.
             calcHFine = .False.
             useCFine  = .False.
             HEffFine  = 0.0e+0_fp
          End If

          ! Check - do we actually have the fine tracer?
          If (id_Fine .eq. 0) Then
             ! Only calculate partitioning of HX into gas and cloud water
             useCFine  = .False. 
             lockCFine = .True. 
             calcHFine = .False.
             HEffFine  = 0.0e+0_fp
          End If
                
          ! Total (molec/cm3)
          COther = State_Chm%Species(I,J,L,id_Other)
          If (useCFine) Then
             CFine = State_Chm%Species(I,J,L,id_Fine)
             CTot  = CFine + COther
          Else
             CFine = 0.0e+0_fp
             CTot  = COther
          End If

          ! Shunt - don't bother if the concentration is too low
          ! If concentration is less than 100 molec/cm3, assume 100% gas phase
          If (CTot.gt.100.0e0_fp) Then

             ! Get info about this species from the species database
             SpcInfo => State_Chm%SpcData(id_Other)%Info
             ! KPP species ID
             K0  = SpcInfo%Henry_K0
             CR  = SpcInfo%Henry_CR
             pKa = SpcInfo%Henry_pKa
      
             SpcInfo => NULL()

             ! Assume equilibrium between X in:
             !  - Gas phase
             !  - Aqueous phase (dissolved in liquid cloud water)
             !  - Aqueous phase (dissolved in surface layer on cloud ice)
             !  - Aqueous phase (dissolved in fine aerosol)
             ! Need an effective Henry's law factor for each (unitless)
             Call Calc_KH( K0, CR, TK, KH, RCLocal )

             ! Get HEff for the liquid available in warm clouds and in the
             ! surface layer of cloud ice
             Call Calc_HEff( pKa, pHLCld, KH, HEffLCld, RCLocal)
             Call Calc_HEff( pKa, pHICld, KH, HEffICld, RCLocal)

             ! Calculate the ratio of the dissolved species to the gas
             ! phase concentration (NOT to the total concentration!)
             rLCld = HEffLCld*VLCld
             rICld = HEffICld*VICld

             ! Do we want to include fine aersol X- in the calculation?
             If (lockCFine.and.useCFine) Then
                ! Calculate a Henry's law constant which will result in no
                ! change in the fine aerosol concentration
                If (CFine.gt.1.0e+0_fp) Then
                   If (cOther.gt.100.0e+0_fp) Then
                      rFine = (1.0e+0_fp + rLCld + rICld)*cTot/cOther
                   Else
                      ! Effectively no gas-phase? This will cause trouble
                      rFine = (1.0e+0_fp + rLCld + rICld)
                   End If
                Else
                   rFine = 0.0e+0_fp
                End If
             ElseIf (calcHFine) Then
                ! Calculate aerosol concentration as just another equilibrium
                ! problem, based on the species data
                Call Calc_HEff( pKa, pHFine, KH, HEffFine, RCLocal)
                rFine = HEffFine*VFine
             Else
                ! Exclude existing dissolved X entirely; no exchange
                rFine = 0.0e+0_fp
             End If

             ! Now get the denominator
             xDenom = 1.0e+0_fp + rFine + rLCld + rICld

             ! xGas is the fraction of the available total present in the gas phase
             xGas = Min(Max(1.0e+0_fp / xDenom,0.0e+0_fp),1.0e+0_fp)

             ! How much is present in fine aerosol, and how much is in "other"
             xFine = xGas * rFine
             If (useCFine .and. (xFine.gt.1.0e-20)) Then
                xOther = 1.0e+0_fp - xFine
             Else
                ! Either no X- made it into the fine aerosol, or the fine
                ! aerosol was excluded from the calculation
                xOther = 1.0e+0_fp
                xFine  = 0.0e+0_fp
             End If

             ! Calculate the fraction of the total in cloud water
             xLCld = xGas * rLCld
             xICld = xGas * rICld

             ! Store the result back in State_Chm
             If (.not.lockCFine) Then
                State_Chm%Species(I,J,L,id_Other) = xOther * CTot
                State_Chm%Species(I,J,L,id_Fine ) = xFine * CTot
             End If
          Else
             ! Nothing here - assume 100% gas phase
             xGas   = 1.0e+0_fp
             xOther = 0.0e+0_fp
             xLCld  = 0.0e+0_fp
             xICld  = 0.0e+0_fp
             xFine  = 0.0e+0_fp
          End If

          ! Store the phase partitioning of remaining HX in State_Chm
          If (xOther .lt. 1.0e-20_fp) Then
             xGas  = 1.0e+0_fp
             xLCld = 0.0e+0_fp
             xICld = 0.0e+0_fp
          Else
             ! Convert to being fractions of total non-fine-aerosol X-
             xGas  = xGas/xOther
             xLCld  = xLCld/xOther
             xICld  = xICld/xOther
          End If
          State_Chm%AqPtn(I,J,L,X,1) = xGas
          State_Chm%AqPtn(I,J,L,X,2) = xLCld
          State_Chm%AqPtn(I,J,L,X,3) = xICld
       End Do
 
    End Do
    End Do
    End Do
    !$OMP END PARALLEL DO

  END SUBROUTINE Partition_Soluble
!EOC
END MODULE Aqueous_Mod
