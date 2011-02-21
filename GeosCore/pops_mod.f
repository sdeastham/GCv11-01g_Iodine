!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: pops_mod.f
!
! !DESCRIPTION: This module contains variables and routines for the 
!  GEOS-Chem peristent organic pollutants (POPs) simulation. 
!\\
!\\
! !INTERFACE: 
!
      MODULE POPS_MOD
! 
! !USES:
!
      IMPLICIT NONE
! Make everything Private ...
      PRIVATE
!
! !PUBLIC TYPES:
!
      PUBLIC :: EMISSPOPS

! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: CHEMPOPS
      PUBLIC :: INIT_POPS

!
! !PUBLIC DATA MEMBERS:
!

! !REVISION HISTORY:
!  20 September 2010 N.E. Selin - Initial Version
!
! !REMARKS:
! Under construction
!
!EOP

!******************************************************************************
!Comment header:
!
!  Module Variables:
!  ===========================================================================
!  (1 ) TCOSZ    (REAL*8) : Sum of COS(Solar Zenith Angle ) [unitless]
!  (2 ) TTDAY    (REAL*8) : Total daylight time at location (I,J) [minutes]
!  (3 ) ZERO_DVEL(REAL*8) : Array with zero dry deposition velocity [cm/s]
!  (4 ) COSZM    (REAL*8) : Max daily value of COS(S.Z. angle) [unitless]  
!
!  Module Routines:
!  ===========================================================================
!  (1 ) CHEMPOPS
!  (2 ) INIT_POPS
!  (3 ) CHEM_POPGP
!  (4 ) EMISSPOPS
!  (5 ) EMITPOP
!  (6 ) OHNO3TIME
!  (7 ) CLEANUP_POPS
!
!  Module Functions:
!  ===========================================================================
!
!
!  GEOS-CHEM modules referenced by pops_mod.f
!  ===========================================================================
!
!
!  Nomenclature: 
!  ============================================================================
!
!
!  POPs Tracers
!  ============================================================================
!  (1 ) POPG               : Gaseous POP - total tracer  
!  (2 ) POPPOC             : OC-sorbed POP  - total tracer
!  (3 ) POPPBC             : BC-sorbed POP  - total tracer
!
!
!  References:
!  ============================================================================
!
!
!  Notes:
!  ============================================================================
!  (1) 20 September 2010 N.E. Selin - Initial version
!  (2) 4 January 2011 C.L. Friedman - Expansion on initial version
!
!
!******************************************************************************
!
      ! References to F90 modules

      

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Parameters
      REAL*8,  PARAMETER   :: SMALLNUM = 1D-20
      ! Arrays
      REAL*8,  ALLOCATABLE :: TCOSZ(:,:)
      REAL*8,  ALLOCATABLE :: TTDAY(:,:)
      REAL*8,  ALLOCATABLE :: ZERO_DVEL(:,:)
      REAL*8,  ALLOCATABLE :: COSZM(:,:)
      REAL*8,  ALLOCATABLE :: EPOP_G(:,:,:)
      REAL*8,  ALLOCATABLE :: EPOP_OC(:,:,:)
      REAL*8,  ALLOCATABLE :: EPOP_BC(:,:,:)
      REAL*8,  ALLOCATABLE :: EPOP_P_TOT(:,:,:)
      REAL*4,  ALLOCATABLE :: POP_TOT_EM(:,:,:)
     
      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  CHEMPOPS
!
! !DESCRIPTION: This routine is the driver routine for POPs chemistry 
!  (eck, 9/20/10)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CHEMPOPS
!
! !INPUT PARAMETERS: 
!
!
! !INPUT/OUTPUT PARAMETERS: 
!
!
! !OUTPUT PARAMETERS:
!
 
      ! References to F90 modules
      USE DRYDEP_MOD,    ONLY : DEPSAV
      USE ERROR_MOD,     ONLY : DEBUG_MSG
      USE GLOBAL_OH_MOD, ONLY : GET_GLOBAL_OH
      USE GLOBAL_OC_MOD, ONLY : GET_GLOBAL_OC !clf, 1/20/2011
      USE GLOBAL_BC_MOD, ONLY : GET_GLOBAL_BC !clf, 1/20/2011
      USE PBL_MIX_MOD,   ONLY : GET_PBL_MAX_L
      USE LOGICAL_MOD,   ONLY : LPRT, LGTMM, LNLPBL !CDH added LNLPBL
      USE TIME_MOD,      ONLY : GET_MONTH, ITS_A_NEW_MONTH
      USE TRACER_MOD,    ONLY : N_TRACERS
!      USE TRACERID_MOD,  ONLY : N_HG_CATS
      USE DRYDEP_MOD,    ONLY : DRYPOPG, DRYPOPP_OC, DRYPOPP_BC

#     include "CMN_SIZE"      ! Size parameters

!
! !REVISION HISTORY: 
!  20 September 2010 - N.E. Selin - Initial Version
!
! !REMARKS:
! (1) Based initially on CHEMMERCURY from MERCURY_MOD (eck, 9/20/10)
!
!EOP
!------------------------------------------------------------------------------
!******************************************************************************
!Comment header
!  Subroutine CHEMPOPS is the driver routine for POPs chemistry
!  in the GEOS-CHEM module. (eck, clf, 1/4/2011)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) 
!
!
!  Local variables:
!  ============================================================================
!  (1 )


!  NOTES:
!  (1 )
!******************************************************************************

!BOC

      ! Local variables
      LOGICAL, SAVE          :: FIRST = .TRUE.
      INTEGER                :: I, J, L, MONTH, N, PBL_MAX


      IF (FIRST) THEN 
          CALL INIT_POPS
          FIRST = .FALSE.
      END IF

      !=================================================================
      ! CHEMPOPS begins here!
      !
      ! Read monthly mean OH fields for oxidation and monthly OC and BC
      ! fields for gas-particle partioning
      !=================================================================
      IF ( ITS_A_NEW_MONTH() ) THEN 

         ! Get the current month
         MONTH = GET_MONTH()

         ! Read monthly mean OH and O3 from disk
         CALL GET_GLOBAL_OH( MONTH )
         IF ( LPRT ) CALL DEBUG_MSG( '### CHEMPOPS: a GET_GLOBAL_OH' )

         ! Read monthly OC from disk
         CALL GET_GLOBAL_OC( MONTH )
         IF ( LPRT ) CALL DEBUG_MSG( '### CHEMPOPS: a GET_GLOBAL_OC' )

         ! Read monthly BC from disk
         CALL GET_GLOBAL_BC( MONTH )
         IF ( LPRT ) CALL DEBUG_MSG( '### CHEMPOPS: a GET_GLOBAL_BC' ) 

      ENDIF
     
      ! If it's a new 6-hr mean, then get the current average 3-D temperature    

      !=================================================================
      ! Perform chemistry on POPs tracers
      !=================================================================
      
      ! Compute diurnal scaling for OH
      CALL OHNO3TIME
      IF ( LPRT ) CALL DEBUG_MSG( 'CHEMPOPS: a OHNO3TIME' )

      !-------------------------
      ! GAS AND PARTICLE PHASE chemistry
      !-------------------------
      IF ( LPRT ) CALL DEBUG_MSG( 'CHEMPOPS: b CHEM_GASPART' )
      
      ! Add option for non-local PBL (cdh, 08/27/09)
      IF ( LNLPBL ) THEN

         ! Dry deposition occurs with PBL mixing,
         ! pass zero deposition frequency
         CALL CHEM_POPGP( ZERO_DVEL, ZERO_DVEL, ZERO_DVEL)
         
      ELSE

         IF ( DRYPOPG > 0 .and. DRYPOPP_OC > 0 .and. DRYPOPP_BC > 0 )
     &      THEN
         
            ! Dry deposition active for both POP-Gas and POP-Particle; 
            ! pass drydep frequency to CHEM_POPGP (NOTE: DEPSAV has units 1/s)
            CALL CHEM_POPGP(DEPSAV(:,:,DRYPOPG), DEPSAV(:,:,DRYPOPP_OC),
     &          DEPSAV(:,:,DRYPOPP_BC) )

           ELSEIF (DRYPOPG > 0 .and. DRYPOPP_OC > 0 .and. 
     &          DRYPOPP_BC .le. 0 ) THEN

            ! Only POPG and POPP_OC dry deposition are active
            CALL CHEM_POPGP(DEPSAV(:,:,DRYPOPG), DEPSAV(:,:,DRYPOPP_OC), 
     &          ZERO_DVEL) 

           ELSEIF (DRYPOPG > 0 .and. DRYPOPP_OC .le. 0 .and. 
     &          DRYPOPP_BC > 0 ) THEN

            ! Only POPG and POPP_BC dry deposition are active
            CALL CHEM_POPGP(DEPSAV(:,:,DRYPOPG), ZERO_DVEL, 
     &          DEPSAV(:,:,DRYPOPP_BC)) 
         
           ELSEIF (DRYPOPG > 0 .and. DRYPOPP_OC .le. 0 .and. 
     &          DRYPOPP_BC .le. 0 ) THEN

            ! Only POPG dry deposition is active
            CALL CHEM_POPGP( DEPSAV(:,:,DRYPOPG), ZERO_DVEL, ZERO_DVEL) 
            
           ELSEIF (DRYPOPG <= 0 .and. DRYPOPP_OC > 0 .and. 
     &          DRYPOPP_BC > 0) THEN

            ! Only POPP dry deposition is active
            CALL CHEM_POPGP( ZERO_DVEL , DEPSAV(:,:,DRYPOPP_OC), 
     &           DEPSAV(:,:,DRYPOPP_BC))

           ELSEIF (DRYPOPG <= 0 .and. DRYPOPP_OC > 0 .and. 
     &          DRYPOPP_BC <= 0) THEN

            ! Only POPP_OC dry deposition is active
            CALL CHEM_POPGP( ZERO_DVEL , DEPSAV(:,:,DRYPOPP_OC), 
     &           ZERO_DVEL)

           ELSEIF (DRYPOPG <= 0 .and. DRYPOPP_OC <= 0 .and. 
     &          DRYPOPP_BC > 0) THEN

            ! Only POPP_OC dry deposition is active
            CALL CHEM_POPGP( ZERO_DVEL , ZERO_DVEL, 
     &           DEPSAV(:,:,DRYPOPP_BC))            
         ELSE

            ! No dry deposition, pass zero deposition frequency
            CALL CHEM_POPGP( ZERO_DVEL, ZERO_DVEL, ZERO_DVEL)

         ENDIF

      ENDIF      

      IF ( LPRT ) CALL DEBUG_MSG( 'CHEMPOPS: a CHEM_GASPART' )
   
    
      ! Return to calling program
      END SUBROUTINE CHEMPOPS

!EOC
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  CHEM_POPGP
!
! !DESCRIPTION: This routine does chemistry for POPs gas and particles
!  (eck, 9/20/10)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CHEM_POPGP (V_DEP_G, V_DEP_P_OC, V_DEP_P_BC)

!     References to F90 modules
      USE TRACER_MOD,   ONLY : STT,        XNUMOL
      USE TRACERID_MOD, ONLY : IDTPOPG,    IDTPOPPOC,  IDTPOPPBC
      USE DIAG53_MOD,   ONLY : AD53_PG_PP, AD53_POPG_OH
      USE DIAG53_MOD,   ONLY : ND53,       LD53
      USE TIME_MOD,     ONLY : GET_TS_CHEM
      USE DIAG_MOD,     ONLY : AD44
      USE LOGICAL_MOD,  ONLY : LNLPBL,     LGTMM
      USE PBL_MIX_MOD,  ONLY : GET_FRAC_UNDER_PBLTOP
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE DAO_MOD,      ONLY : T,          AIRVOL
      USE ERROR_MOD,    ONLY : DEBUG_MSG

#     include "CMN_SIZE" ! Size parameters
#     include "CMN_DIAG" ! ND44

!
! !INPUT PARAMETERS: 
!
      REAL*8, INTENT(IN)    :: V_DEP_G(IIPAR,JJPAR)
      REAL*8, INTENT(IN)    :: V_DEP_P_OC(IIPAR,JJPAR)
      REAL*8, INTENT(IN)    :: V_DEP_P_BC(IIPAR,JJPAR)

!
! !INPUT/OUTPUT PARAMETERS: 
!
!
! !OUTPUT PARAMETERS:
!
    
!
!
! !REVISION HISTORY: 
!  20 September 2010 - N.E. Selin - Initial Version
!
! !REMARKS:
! (1) Based initially on CHEM_HG0_HG2 from MERCURY_MOD (eck, 9/20/10)
!
!EOP
!------------------------------------------------------------------------------
!******************************************************************************
!Comment header
!  Subroutine CHEM_POPGP is the chemistry subroutine for the oxidation,
!  gas-particle partitioning, and deposition of POPs.
!  (eck, clf, 1/4/2011)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) V_DEP_G (REAL*8)    : Dry deposition frequency for gaseous POP [/s]
!  (2 ) V_DEP_P_OC (REAL*8) : Dry deposition frequency for OC-POP [/s]
!  (3 ) V_DEP_P_BC (REAL*8) : Dry deposition frequency for BC-POP [/s]
!
!  Local variables:
!  ============================================================================
!  (1 )
!
!
!     
!  NOTES:
!  (1 ) 
!  
!  REFS:
!  (1 ) For OH rate constant: Brubaker & Hites. 1998. OH reaction kinetics of
!  PAHs and PCDD/Fs. J. Phys. Chem. A. 102:915-921. 
!
!******************************************************************************
!BOC
      ! Local variables
      INTEGER               :: I, J, L, NN
      REAL*8                :: DTCHEM
      REAL*8                :: KOA_T,        KBC_T
      REAL*8                :: KOC_BC_T,     KBC_OC_T
      REAL*8                :: TK
      REAL*8                :: AREA_CM2
      REAL*8                :: F_PBL
      REAL*8                :: C_OH,         C_OC,        C_BC
      REAL*8                :: K_OH
      REAL*8                :: K_OX
      REAL*8                :: E_KOX_T
      REAL*8                :: K_DEPG,       K_DEPP_OC,   K_DEPP_BC
      REAL*8                :: OLD_POPG,     OLD_POPP_OC, OLD_POPP_BC
      REAL*8                :: NEW_POPG,     NEW_POPP_OC, NEW_POPP_BC
      REAL*8                :: POPG_BL,      POPP_OC_BL,  POPP_BC_BL
      REAL*8                :: POPG_FT,      POPP_OC_FT,  POPP_BC_FT
      REAL*8                :: TMP_POPG,     TMP_OX!,     TMP_POPP
      REAL*8                :: GROSS_OX,     GROSS_OX_OH, NET_OX
      REAL*8                :: DEP_POPG,     DEP_POPP_OC, DEP_POPP_BC
      REAL*8                :: DEP_POPG_DRY, DEP_POPP_OC_DRY
      REAL*8                :: DEP_POPP_BC_DRY
      REAL*8                :: DEP_DRY_FLXG, DEP_DRY_FLXP_OC
      REAL*8                :: DEP_DRY_FLXP_BC
      REAL*8                :: OLD_POP_T
      REAL*8                :: VR_OC_AIR,    VR_BC_AIR
      REAL*8                :: VR_OC_BC,     VR_BC_OC
      REAL*8                :: F_POP_OC,     F_POP_BC,    F_POP_G
      REAL*8                :: MPOP_OC,      MPOP_BC,     MPOP_G

      ! Delta H for POP [kJ/mol]. Delta H is enthalpy of phase transfer
      ! from gas phase to OC. For now we use Delta H for phase transfer 
      ! from the gas phase to the pure liquid state. For PHENANTHRENE, 
      ! this is taken as the negative of the Delta H for phase transfer
      ! from the pure liquid state to the gas phase (Schwarzenbach,
      !  Gschwend, Imboden, 2003, pg 200, Table 6.3), or -74000 [J/mo]l.
      REAL*8, PARAMETER     :: DEL_H      = -74d3

      ! R = universal gas constant for adjusting KOA for temp: 8.3145 [J/mol/K]
      REAL*8, PARAMETER     :: R          = 8.31d0  

      ! KOA_298 for partitioning of gas phase POP to atmospheric OC
      ! KOA_298 = Cpop in octanol/Cpop in atmosphere at 298 K 
      ! For phenanthrene, KOA_298 = log 7.64, or 4.37*10^7 [unitless]
      ! (Ma et al., J. Chem. Eng. Data, 2010, 55:819-825).
      REAL*8, PARAMETER     :: KOA_298    = 4.37d7

      ! KBC_298 for partitioning of gas phase POP to atmospheric BC
      ! KBC_298 = Cpop in black carbon/Cpop in atmosphere at 298 K
      REAL*8, PARAMETER     :: KBC_298    = 1d10

      ! DENS_OCT = density of octanol, needed for partitioning into OC
      ! 820 [kg/m^3]
      REAL*8, PARAMETER     :: DENS_OCT   = 82d1

      ! DENS_BC = density of BC, needed for partitioning onto BC
      ! 1 [kg/L] or 1000 [kg/m^3] 
      ! From Lohmann and Lammel, Environ. Sci. Technol., 2004, 38:3793-3803.
      REAL*8, PARAMETER     :: DENS_BC    = 1d3

      ! K for reaction POPG + OH  [cm3 /molecule /s]
      ! Currently set for phenanthrene (Source: Brubaker & Hites, 1998)
      ! Could potentially set this to change with temmperature (clf)
      REAL*8, PARAMETER     :: K_POPG_OH  = 2.70d-11 !(Gas phase)

      ! K for reaction POPP + NO3 could be added here someday

      !=================================================================
      ! CHEM_POPGP begins here!
      !=================================================================

      ! Chemistry timestep [s]
      DTCHEM = GET_TS_CHEM() * 60d0


! Here, put some code that reacts popg into popp and 
! oxidizes popg and
! deposits both 

      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Save local temperature in TK for convenience [K]
         TK = T(I,J,L)

         ! Get monthly mean OH, OC, and BC concentrations
         C_OH        = GET_OH( I, J, L )
         C_OC        = GET_OC( I, J, L )
         C_BC        = GET_BC( I, J, L ) 

         ! Define KOA_T, the octanol-air partition coeff at temp T [K]
         KOA_T = KOA_298 * EXP((-DEL_H/R) * ((1d0/TK) - (1d0/298d0)))

         ! Define KBC_T, the BC-air partition coeff at temp T [K]
         KBC_T = KBC_298 * EXP((-DEL_H/R) * ((1d0/TK) - (1d0/298d0)))

         ! Fraction of box (I,J,L) underneath the PBL top [dimensionless]
         F_PBL = GET_FRAC_UNDER_PBLTOP( I, J, L )

         ! Define K for the oxidation reaction with POPG [/s]
         K_OH        = K_POPG_OH * C_OH
         
         ! Could add K for oxidation by NO3 here one day [/s]

         ! Total K for oxidation [/s]
         K_OX        = K_OH !+ ...

         ! Add a K for a POPP oxidation rxn here someday [/s]
         ! (Data from J. Kroll?)

         ! Define Ks for dry deposition of gas phase POP [/s]
         K_DEPG = V_DEP_G(I,J)

         ! Define Ks for dry deposition of particle phase POP [/s]
         K_DEPP_OC = V_DEP_P_OC(I,J)

         ! Define Ks for dry deposition of particle phase POP [/s]
         K_DEPP_BC = V_DEP_P_BC(I,J)

         ! Precompute exponential factors [dimensionless]
         E_KOX_T  = EXP( -K_OX  * DTCHEM )

         !==============================================================
         ! GAS-PARTICLE PARTITIONING
         !==============================================================

         OLD_POPG = MAX( STT(I,J,L,IDTPOPG), SMALLNUM )  ![kg]

         OLD_POPP_OC = MAX( STT(I,J,L,IDTPOPPOC), SMALLNUM )  ![kg]

         OLD_POPP_BC = MAX( STT(I,J,L,IDTPOPPBC), SMALLNUM )  ![kg]

         ! Total POPs in box I,J,L 
         OLD_POP_T = OLD_POPG + OLD_POPP_OC + OLD_POPP_BC

         ! Get monthly mean OC and BC concentrations [kg/box]
         C_OC        = GET_OC( I, J, L )
         C_BC        = GET_BC( I, J, L )

         ! Convert C_OC and C_BC units to volume per box 
         ! [m^3 OC or BC/box]
         C_OC        = C_OC / DENS_OCT
         C_BC        = C_BC / DENS_BC


         ! Define volume ratios:
         ! VR_OC_AIR = volume ratio of OC to air [unitless]     

         VR_OC_AIR   = C_OC / AIRVOL(I,J,L)

         ! VR_OC_BC  = volume ratio of OC to BC [unitless]
         VR_OC_BC    = C_OC / C_BC

         ! VR_BC_AIR = volume ratio of BC to air [unitless]
         VR_BC_AIR   = VR_OC_AIR / VR_OC_BC

         ! VR_BC_OC  = volume ratio of BC to OC [unitless]
         VR_BC_OC    = 1d0 / VR_OC_BC 

         ! Define temperature-dependent partition coefficients:
         ! KOA_T, the octanol-air partition coeff at temp T [unitless]
         KOA_T = KOA_298 * EXP((-DEL_H/R) * ((1d0/TK) - (1d0/298d0)))

         ! Define KBC_T, the BC-air partition coeff at temp T [unitless]
         KBC_T = KBC_298 * EXP((-DEL_H/R) * ((1d0/TK) - (1d0/298d0)))

         ! Define KOC_BC_T, the theoretical OC-BC part coeff at temp T [unitless]
         KOC_BC_T = KOA_T / KBC_T

         ! Define KBC_OC_T, the theoretical BC_OC part coeff at temp T [unitless]
         KBC_OC_T = 1d0 / KOC_BC_T

         ! Redefine fractions of total POPs in box (I,J,L) that are OC-phase, 
         ! BC-phase, and gas phase 
         ! OC-phase:
         F_POP_OC  = 1d0 / (1d0 + (1d0 / (KOA_T * VR_OC_AIR)) + 
     &               (1d0 / (KOC_BC_T * VR_OC_BC)))
            
         ! BC-phase:
         F_POP_BC  = 1d0 / (1d0 + (1d0 / (KBC_T * VR_BC_AIR)) + 
     &               (1d0 / (KBC_OC_T * VR_BC_OC)))

         ! Gas-phase:
         F_POP_G   = 1d0 - F_POP_OC - F_POP_BC 

         ! Calculate new masses of POP in each phase [kg/s]
         ! OC-phase:
         MPOP_OC    = F_POP_OC * OLD_POP_T

         ! BC-phase
         MPOP_BC     = F_POP_BC * OLD_POP_T

         ! Gas-phase
         MPOP_G      = F_POP_G  * OLD_POP_T        

         !==============================================================
         ! CHEMISTRY AND DEPOSITION REACTIONS
         !==============================================================
         
         IF ( F_PBL < 0.05D0 .OR. 
     &           K_DEPG < SMALLNUM ) THEN

               !==============================================================
               ! Entire box is in the free troposphere
               ! or deposition is turned off, so use RXN without deposition
               ! for gas phase POPs
               ! For particle POPs, no rxn and no deposition
               !==============================================================

               CALL RXN_OX_NODEP( MPOP_G, K_OX, 
     &              DTCHEM, E_KOX_T, NEW_POPG,
     &              GROSS_OX )

               NEW_POPP_OC = MPOP_OC
               NEW_POPP_BC = MPOP_BC

               ! No deposition occurs [kg]
               DEP_POPG = 0D0
               DEP_POPP_OC = 0D0
               DEP_POPP_BC = 0D0
               

            ELSE IF ( F_PBL > 0.95D0 ) THEN 

               !==============================================================
               ! Entire box is in the boundary layer
               ! so use RXN with deposition for gas phase POPs
               ! Deposition only (no rxn) for particle phase POPs
               !==============================================================

               CALL RXN_OX_WITHDEP( MPOP_G,   K_OX,
     &              K_DEPG,   DTCHEM,  E_KOX_T, NEW_POPG,
     &              GROSS_OX,  DEP_POPG )

               CALL NO_RXN_WITHDEP( MPOP_OC, K_DEPP_OC, DTCHEM,
     &              NEW_POPP_OC, DEP_POPP_OC )

               CALL NO_RXN_WITHDEP( MPOP_BC, K_DEPP_BC, DTCHEM,
     &              NEW_POPP_BC, DEP_POPP_BC )

            ELSE

               !==============================================================
               ! Box spans the top of the boundary layer
               ! Part of the mass is in the boundary layer and subject to 
               ! deposition while part is in the free troposphere and
               ! experiences no deposition.
               !
               ! We apportion the mass between the BL and FT according to the
               ! volume fraction of the box in the boundary layer.
               ! Arguably we should assume uniform mixing ratio, instead of
               ! uniform density but if the boxes are short, the air density
               ! doesn't change much.
               ! But assuming uniform mixing ratio across the inversion layer
               ! is a poor assumption anyway, so we are just using the
               ! simplest approach.
               !==============================================================

               ! Boundary layer portion of POPG [kg]
               POPG_BL = MPOP_G * F_PBL 

               ! Boundary layer portion of POPP_OC [kg]
               POPP_OC_BL = MPOP_OC * F_PBL

               ! Boundary layer portion of POPP_BC [kg]
               POPP_BC_BL = MPOP_BC * F_PBL

               ! Free troposphere portion of POPG [kg]
               POPG_FT = MPOP_G - POPG_BL

               ! Free troposphere portion of POPP_OC [kg]
               POPP_OC_FT = MPOP_OC - POPP_OC_BL

               ! Free troposphere portion of POPP_BC [kg]
               POPP_BC_FT = MPOP_BC - POPP_BC_BL
               
               ! Do chemistry with deposition on BL fraction for gas phase
               CALL RXN_OX_WITHDEP( POPG_BL,  K_OX,
     &              K_DEPG,   DTCHEM, 
     &              E_KOX_T, NEW_POPG, 
     &              GROSS_OX,  DEP_POPG )

               ! Do chemistry without deposition on the FT fraction for gas phase
               CALL RXN_OX_NODEP(  POPG_FT, K_OX,
     &              DTCHEM, E_KOX_T,
     &              TMP_POPG, TMP_OX )

               ! Do deposition (no chemistry) on BL fraction for particulate phase
               ! No deposition (and no chem) on the FT fraction
               ! for the particulate phase

               CALL NO_RXN_WITHDEP( POPP_OC_BL, K_DEPP_OC, DTCHEM,  
     &              NEW_POPP_OC ,DEP_POPP_OC)

               CALL NO_RXN_WITHDEP( POPP_BC_BL, K_DEPP_BC, DTCHEM,  
     &              NEW_POPP_BC, DEP_POPP_BC)
               
               ! Recombine the boundary layer and free troposphere parts [kg]
               NEW_POPG = NEW_POPG + TMP_POPG
               NEW_POPP_OC = NEW_POPP_OC
               NEW_POPP_BC = NEW_POPP_BC
               
               ! Total gross oxidation of gas phase in the BL and FT [kg]
               GROSS_OX = GROSS_OX + TMP_OX

            ENDIF

            ! Ensure positive concentration [kg]
            NEW_POPG = MAX( NEW_POPG, SMALLNUM )
            NEW_POPP_OC = MAX( NEW_POPP_OC, SMALLNUM )
            NEW_POPP_BC = MAX( NEW_POPP_BC, SMALLNUM )

            ! Archive new POPG and POPP values [kg]
            STT(I,J,L,IDTPOPG) = NEW_POPG
            STT(I,J,L,IDTPOPPOC) = NEW_POPP_OC
            STT(I,J,L,IDTPOPPBC) = NEW_POPP_BC

            ! Net oxidation [kg] (equal to gross ox for now)
            NET_OX = OLD_POPG - NEW_POPG - DEP_POPG

            ! Error check on gross oxidation [kg]
            IF ( GROSS_OX < 0D0 ) 
     &          CALL DEBUG_MSG('CHEM_POPGP: negative gross oxidation')

            ! Apportion gross oxidation between OH and possibly
            ! NO3 someday [kg]
            IF ( (K_OX     < SMALLNUM) .OR. 
     &           (GROSS_OX < SMALLNUM) ) THEN
               GROSS_OX_OH = 0D0
!               GROSS_OX_NO3 = 0D0
            ELSE
               GROSS_OX_OH = GROSS_OX ! * K_OH / K_OX
!               GROSS_OX_NO3 = GROSS_OX * K_NO3 / K_OX
            ENDIF

            ! Apportion deposition [kg]
            ! Right now only using dry deposition (no sea salt) (clf, 1/27/11)
            IF ( (K_DEPG  < SMALLNUM) .OR. 
     &           (DEP_POPG < SMALLNUM) ) THEN
               DEP_POPG_DRY  = 0D0
            ELSE
               DEP_POPG_DRY  = DEP_POPG ! If ever use dep with sea
            ! salt aerosols, will need to multiply DEP_POPG by the ratio 
            ! of K_DRYG (rate of dry dep) to K_DEPG (total dep rate).  
            ENDIF

            IF ( (K_DEPP_OC  < SMALLNUM) .OR. 
     &           (DEP_POPP_OC < SMALLNUM) ) THEN
               DEP_POPP_OC_DRY  = 0D0
            ELSE
               DEP_POPP_OC_DRY  = DEP_POPP_OC 
            ENDIF

            IF ( (K_DEPP_BC  < SMALLNUM) .OR. 
     &           (DEP_POPP_BC < SMALLNUM) ) THEN
               DEP_POPP_BC_DRY  = 0D0
            ELSE
               DEP_POPP_BC_DRY  = DEP_POPP_BC 
            ENDIF

            !=================================================================
            ! ND44 diagnostic: drydep flux of POPG and POPP [molec/cm2/s]
            !=================================================================
            IF ( ( ND44 > 0 .OR. LGTMM ) .AND. (.NOT. LNLPBL) ) THEN
            ! Not using LGTMM right now (logical switch for using GTMM soil model)
            ! Also not using non-local PBL mode yet (clf, 1/27/2011)

               ! Grid box surface area [cm2]
               AREA_CM2 = GET_AREA_CM2( J )

               ! Amt of POPG lost to drydep [molec/cm2/s]
               DEP_DRY_FLXG  = DEP_POPG_DRY * XNUMOL(IDTPOPG) / 
     &              ( AREA_CM2 * DTCHEM )

               ! Archive POPG drydep flux in AD44 array [molec/cm2/s]
               AD44(I,J,IDTPOPG,1) = AD44(I,J,IDTPOPG,1) +
     &              DEP_DRY_FLXG

               ! Amt of POPPOC lost to drydep [molec/cm2/s]
               DEP_DRY_FLXP_OC = DEP_POPP_OC_DRY * 
     &                 XNUMOL(IDTPOPPOC)/( AREA_CM2 * DTCHEM )        

               ! Archive POPPOC drydep flux in AD44 array [molec/cm2/s]
               AD44(I,J,IDTPOPPOC,1) = 
     &              AD44(I,J,IDTPOPPOC,1) + DEP_DRY_FLXP_OC

               ! Amt of POPPBC lost to drydep [molec/cm2/s]
               DEP_DRY_FLXP_BC = DEP_POPP_BC_DRY * 
     &                 XNUMOL(IDTPOPPBC)/( AREA_CM2 * DTCHEM )        

               ! Archive POPPBC drydep flux in AD44 array [molec/cm2/s]
               AD44(I,J,IDTPOPPBC,1) = 
     &              AD44(I,J,IDTPOPPBC,1) + DEP_DRY_FLXP_BC


            ENDIF
           

            !==============================================================
            ! ND53 diagnostic: Oxidized POPG (OH-POPG) production [kg]
            !==============================================================

            IF ( ND53 > 0 .AND. L <= LD53 ) THEN ! LD53 is max level

               ! Store chemistry diagnostics only for total tracer
               !IF ( ID_HG0(NN) == ID_HG_TOT) THEN
               ! Don't think we need this - POP only has one regional category
               ! right now (clf, 1/27/11)

               !AD53_POPG_OX(I,J,L)= AD53_POPG_OX(I,J,L) + NET_OX
               ! Can add this later if we have more than one oxidant
               AD53_POPG_OH(I,J,L) = AD53_POPG_OH(I,J,L)  + 
     &              GROSS_OX_OH
                  


               !ENDIF
            ENDIF

      ENDDO
      ENDDO
      ENDDO


! save into stt
! END OMP stuff here if added

      END SUBROUTINE CHEM_POPGP 

!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  RXN_OX_NODEP
!
! !DESCRIPTION: Subroutine RXN_OX_NODEP calculates new mass of POPG for given
! oxidation rates, without any deposition. This is for the free troposphere, or
! simulations with deposition turned off. (clf, 1/27/11, based on RXN_REDOX_NODEP
! in mercury_mod.f).
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE RXN_OX_NODEP( OLD_POPG, K_OX, DT, E_KOX_T,
     &     NEW_POPG, GROSS_OX )
!

! !INPUT PARAMETERS: 
      REAL*8,  INTENT(IN)  :: OLD_POPG,  DT
      REAL*8,  INTENT(IN)  :: K_OX
      REAL*8,  INTENT(IN)  :: E_KOX_T
      
!
! !INPUT/OUTPUT PARAMETERS:   
!
!
! !OUTPUT PARAMETERS:
      REAL*8,  INTENT(OUT) :: NEW_POPG,  GROSS_OX
!
! !REVISION HISTORY: 
!  27 January 2011 - CL Friedman - Initial Version
!
! !REMARKS:
! (1) Based on RXN_REDOX_NODEP in mercury_mod.f
!
!EOP
!------------------------------------------------------------------------------
!******************************************************************************
!Comment header
!  Subroutine RXN_OX_NODEP calculates new mass of POPG for given
! oxidation rates, without any deposition. This is for the free troposphere, or
! simulations with deposition turned off. (clf, 1/27/11, based on RXN_REDOX_NODEP
! in mercury_mod.f).
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) OLD_POPG (REAL*8) : 
!  (2 ) DT       (REAL*8) : 
!  (3 ) K_OX     (REAL*8) :
!  (4 ) E_KOX_T  (REAL*8) :
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) NEW_POPG (REAL*8) :
!  (2 ) GROSS_OX (REAL*8) :
!
!  Local variables:
!  ============================================================================
!  (1 )
!
!  NOTES:
!  (1 ) 
!  
!  REFS:
!  (1 )  
!
!******************************************************************************
!BOC
      
      ! Local variables
      ! None

      !=================================================================
      ! RXN_OX_NODEP begins here!
      !=================================================================

         !=================================================================
         ! Oxidation
         !=================================================================

         ! New concentration of POPG
         NEW_POPG = OLD_POPG * E_KOX_T

         ! Gross oxidation is the same as net oxidation
         GROSS_OX = OLD_POPG - NEW_POPG

      END SUBROUTINE RXN_OX_NODEP
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  RXN_OX_WITHDEP
!
! !DESCRIPTION: Subroutine RXN_OX_WITHDEP calculates new mass of POPG for given
! rates of oxidation and deposition. This is for the boundary layer.
! (clf, 1/27/11, based on RXN_REDOX_NODEP in mercury_mod.f).
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE RXN_OX_WITHDEP( OLD_POPG, K_OX, K_DEPG, DT, E_KOX_T,
     &     NEW_POPG, GROSS_OX, DEP_POPG )
!
      ! References to F90 modules
      USE ERROR_MOD,    ONLY : ERROR_STOP


! !INPUT PARAMETERS: 
      REAL*8,  INTENT(IN)  :: OLD_POPG,  DT
      REAL*8,  INTENT(IN)  :: K_OX, K_DEPG
      REAL*8,  INTENT(IN)  :: E_KOX_T
      
!
! !INPUT/OUTPUT PARAMETERS:   
!
!
! !OUTPUT PARAMETERS:
      REAL*8,  INTENT(OUT) :: NEW_POPG,  GROSS_OX
      REAL*8,  INTENT(OUT) :: DEP_POPG
!
! !REVISION HISTORY: 
!  27 January 2011 - CL Friedman - Initial Version
!
! !REMARKS:
! (1) Based on RXN_REDOX_WITHDEP in mercury_mod.f
!
!EOP
!------------------------------------------------------------------------------
!******************************************************************************
!Comment header
!  Subroutine RXN_OX_WITHDEP calculates new mass of POPG for given
! rates of oxidation and deposition. This is for the boundary layer.
! (clf, 1/27/11, based on RXN_REDOX_WITHDEP in mercury_mod.f).
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) OLD_POPG (REAL*8) : 
!  (2 ) DT       (REAL*8) : 
!  (3 ) K_OX     (REAL*8) :
!  (4 ) K_DEPG   (REAL*8) :
!  (5 ) E_KOX_T  (REAL*8) :
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) NEW_POPG (REAL*8) :
!  (2 ) GROSS_OX (REAL*8) :
!  (3 ) DEP_POPG (REAL*8) :
!
!  Local variables:
!  ============================================================================
!  (1 )
!
!  NOTES:
!  (1 ) 
!  
!  REFS:
!  (1 )  
!
!******************************************************************************
!BOC
      
      ! Local variables
      REAL*8               :: E_KDEPG_T

      !=================================================================
      ! RXN_OX_WITHDEP begins here!
      !=================================================================

      ! Precompute exponential factor for deposition [dimensionless]
      E_KDEPG_T = EXP( -K_DEPG * DT )

      IF (K_OX < SMALLNUM) THEN      

         !=================================================================
         ! No Chemistry, Deposition only
         !=================================================================

         ! New mass of POPG [kg]
         NEW_POPG = OLD_POPG * E_KDEPG_T

         ! Oxidation of POPG [kg]
         GROSS_OX = 0D0

         ! Deposited POPG [kg]
         DEP_POPG = OLD_POPG - NEW_POPG

      ELSE

         !=================================================================
         ! Oxidation and Deposition 
         !=================================================================

         ![POPG](t) = [POPG](0) exp( -(kOx + kDPOPG) t)

         !Ox(t)     = ( [POPG](0) - [POPG](t) ) * kOx / ( kOx + kDPOPG )

         !Dep_POPG(t)   = ( [POPG](0) - [POPG](t) - Ox(t) ) 

         ! New concentration of POPG [kg]
         NEW_POPG = OLD_POPG * E_KOX_T * E_KDEPG_T 

         ! Gross oxidized gas phase mass [kg]
         GROSS_OX = ( OLD_POPG - NEW_POPG ) * K_OX / ( K_OX + K_DEPG )
         GROSS_OX = MAX( GROSS_OX, 0D0 )

         ! POPG deposition [kg]
         DEP_POPG = ( OLD_POPG - NEW_POPG - GROSS_OX )
         DEP_POPG = MAX( DEP_POPG, 0D0 )

      ENDIF

      END SUBROUTINE RXN_OX_WITHDEP
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  NO_RXN_WITHDEP
!
! !DESCRIPTION: Subroutine NO_RXN_WITHDEP calculates new mass of POPP for given
! rate of deposition. No oxidation of POPP. This is for the boundary layer.
! (clf, 2/9/11)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE NO_RXN_WITHDEP( OLD_POPP, K_DEPP, DT,
     &     NEW_POPP, DEP_POPP )
!
      ! References to F90 modules
      USE ERROR_MOD,    ONLY : ERROR_STOP


! !INPUT PARAMETERS: 
      REAL*8,  INTENT(IN)  :: OLD_POPP
      REAL*8,  INTENT(IN)  :: K_DEPP
      REAL*8,  INTENT(IN)  :: DT
      
!
! !INPUT/OUTPUT PARAMETERS:   
!
!
! !OUTPUT PARAMETERS:
      REAL*8,  INTENT(OUT) :: NEW_POPP
      REAL*8,  INTENT(OUT) :: DEP_POPP
!
! !REVISION HISTORY: 
!  9 February 2011 - CL Friedman - Initial Version
!
! !REMARKS:
!
!EOP
!------------------------------------------------------------------------------
!******************************************************************************
!Comment header
!  Subroutine NO_RXN_WITHDEP calculates new mass of POPP for given
! rate of deposition. This is for the boundary layer.
! (clf, 1/27/11, based on RXN_REDOX_NODEP in mercury_mod.f).
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) OLD_POPP (REAL*8) : 
!  (2 ) K_DEPP   (REAL*8) :
!  (3 ) DT       (REAL*8) :
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) NEW_POPP (REAL*8) :
!  (2 ) DEP_POPP (REAL*8) :
!
!  Local variables:
!  ============================================================================
!  (1 )
!
!  NOTES:
!  (1 ) 
!  
!  REFS:
!  (1 )  
!
!******************************************************************************
!BOC
      
      ! Local variables
      REAL*8               :: E_KDEPP_T

      !=================================================================
      ! NO_RXN_WITHDEP begins here!
      !=================================================================

      ! Precompute exponential factors [dimensionless]
      E_KDEPP_T = EXP( -K_DEPP * DT )     

      !=================================================================
      ! No Chemistry, Deposition only
      !=================================================================

      ! New mass of POPP [kg]
      NEW_POPP = OLD_POPP * E_KDEPP_T

      ! Oxidation of POPP [kg]
      ! GROSS_OX = 0D0

      ! POPP deposition [kg]
      DEP_POPP = OLD_POPP - NEW_POPP
      DEP_POPP = MAX( DEP_POPP, 0D0 )


      END SUBROUTINE NO_RXN_WITHDEP
!EOC

!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  EMISSPOPS
!
! !DESCRIPTION: This routine is the driver routine for POPs emissions
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE EMISSPOPS
!
! !INPUT PARAMETERS: 
!

!
! !INPUT/OUTPUT PARAMETERS: 
!
!
! !OUTPUT PARAMETERS:
!
!
! !REVISION HISTORY: 
!  20 September 2010 - N.E. Selin - Initial Version
!
! !REMARKS:
! (1) Based initially on EMISSMERCURY from MERCURY_MOD (eck, 9/20/10)
!
!EOP
!------------------------------------------------------------------------------
!******************************************************************************
!Comment header
!  Subroutine EMISSPOPS is the driver subroutine for POPs emissions.
!  !
!  Arguments as Input:
!  ============================================================================
!  (1 )
!
!  Local variables:
!  ============================================================================
!  (1 ) I, J, L   (INTEGER) : Long, lat, level
!  (2 ) N         (INTEGER) : Tracer ID
!  (3 ) PBL_MAX   (INTEGER) : Maximum extent of boundary layer [level]
!  (4 ) DTSRCE    (REAL*8)  : Emissions time step  [s]
!  (5 ) T_POP     (REAL*8)  : POP emission rate [kg/s]
!  (6 ) E_POP     (REAL*8)  : POPs emitted into box [kg]
!  (7 ) F_OF_PBL  (REAL*8)  : Fraction of box within boundary layer [unitless]
!     
!  NOTES:
!  (1 ) 
!  
!  REFS:
!  (1 )
!******************************************************************************
!BOC
      ! References to F90 modules
      USE ERROR_MOD,         ONLY : DEBUG_MSG, ERROR_STOP
      USE LOGICAL_MOD,       ONLY : LPRT, LNLPBL !CDH added LNLPBL
      USE TIME_MOD,          ONLY : GET_MONTH, ITS_A_NEW_MONTH
      USE TRACER_MOD,        ONLY : STT
      USE VDIFF_PRE_MOD,     ONLY : EMIS_SAVE !cdh for LNLPBL
      USE GRID_MOD,          ONLY : GET_XMID, GET_YMID
      USE DAO_MOD,           ONLY : T, AIRVOL
      ! Reference to diagnostic arrays
      USE DIAG53_MOD,   ONLY : AD53, ND53!, AD53_PG_PP
      USE PBL_MIX_MOD,  ONLY : GET_FRAC_OF_PBL, GET_PBL_MAX_L
      USE TIME_MOD,     ONLY : GET_TS_EMIS
      USE TRACERID_MOD, ONLY : IDTPOPG, IDTPOPPOC,  IDTPOPPBC
      
#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_DEP"      ! FRCLND

      ! Local variables
      INTEGER               :: I,   J,    L,    N,    PBL_MAX
      REAL*8                :: DTSRCE, F_OF_PBL, TK
      REAL*8                :: E_POP, T_POP
      REAL*8                :: C_OC,        C_BC
      REAL*4                :: F_POP_OC, F_POP_BC 
      REAL*8                :: F_POP_G
      REAL*8                :: KOA_T, KBC_T, KOC_BC_T, KBC_OC_T
      REAL*8                :: VR_OC_AIR, VR_OC_BC
      REAL*8                :: VR_BC_AIR, VR_BC_OC
      REAL*8                :: SUM_OC_EM, SUM_BC_EM, SUM_G_EM
      LOGICAL, SAVE         :: FIRST = .TRUE.

      ! Delta H for POP [kJ/mol]. Delta H is enthalpy of phase transfer
      ! from gas phase to OC. For now we use Delta H for phase transfer 
      ! from the gas phase to the pure liquid state. For PHENANTHRENE, 
      ! this is taken as the negative of the Delta H for phase transfer
      ! from the pure liquid state to the gas phase (Schwarzenbach,
      !  Gschwend, Imboden, 2003, pg 200, Table 6.3), or -74000 [J/mol].
      REAL*8, PARAMETER     :: DEL_H      = -74d3

      ! R = universal gas constant for adjusting KOA for temp: 8.3145 [J/mol/K]
      REAL*8, PARAMETER     :: R          = 8.31d0  

      ! KOA_298 for partitioning of gas phase POP to atmospheric OC
      ! KOA_298 = Cpop in octanol/Cpop in atmosphere at 298 K 
      ! For phenanthrene, KOA_298 = log 7.64, or 4.37*10^7 [unitless]
      ! (Ma et al., J. Chem. Eng. Data, 2010, 55:819-825).
      REAL*8, PARAMETER     :: KOA_298    = 4.37d7

      ! KBC_298 for partitioning of gas phase POP to atmospheric BC
      ! KBC_298 = Cpop in black carbon/Cpop in atmosphere at 298 K
      REAL*8, PARAMETER     :: KBC_298    = 1d10

      ! DENS_OCT = density of octanol, needed for partitioning into OC
      ! 820 [kg/m^3]
      REAL*8, PARAMETER     :: DENS_OCT   = 82d1

      ! DENS_BC = density of BC, needed for partitioning onto BC
      ! 1 [kg/L] or 1000 [kg/m^3]
      ! From Lohmann and Lammel, Environ. Sci. Technol., 2004, 38:3793-3803.
      REAL*8, PARAMETER     :: DENS_BC    = 1d3


      !=================================================================
      ! EMISSPOPS begins here!
      !=================================================================

      ! First-time initialization
      IF ( FIRST ) THEN

         ! Read anthro emissions from disk
         CALL POPS_READYR

         ! Reset first-time flag
         FIRST = .FALSE.
      ENDIF

      ! If we are using the non-local PBL mixing,
      ! we need to initialize the EMIS_SAVE array (cdh, 08/27/09)
      IF (LNLPBL) EMIS_SAVE = 0d0

      ! Emission timestep [s]
      DTSRCE  = GET_TS_EMIS() * 60d0

      ! Maximum extent of the PBL [model levels]
      PBL_MAX = GET_PBL_MAX_L() 

      ! Save local temperature in TK for convenience [K]
      TK = T(I,J,L)

      ! Loop over grid boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N,  T_POP, F_OF_PBL, EPOP_G, EPOP_OC )
!$OMP+PRIVATE( EPOP_BC, EPOP_P_TOT)
      DO J = 1, JJPAR
      DO I = 1, IIPAR          

         !Here, save the total from the anthro emissions array
         !into the T_POP variable [kg/s]
         T_POP = POP_TOT_EM(I,J,1)
 
         !==============================================================
         ! Apportion total POPs emitted to gas phase, OC-bound, and BC-bound
         ! emissions (clf, 2/1/2011)         
         ! Then partition POP throughout PBL; store into STT [kg]
         ! Now make sure STT does not underflow (cdh, bmy, 4/6/06; eck 9/20/10)
         !==============================================================

         ! Loop up to max PBL level
         DO L = 1, PBL_MAX

            ! Get monthly mean OC and BC concentrations [kg/box]
            C_OC        = GET_OC( I, J, L )
            C_BC        = GET_BC( I, J, L )

            ! Convert C_OC and C_BC units to volume per box 
            ! [m^3 OC or BC/box]
            C_OC        = C_OC / DENS_OCT
            C_BC        = C_BC / DENS_BC


            ! Define volume ratios:
            ! VR_OC_AIR = volume ratio of OC to air [unitless]

            ! AIRVOL     = AIRVOL(I,J,L)     

            VR_OC_AIR   = C_OC / AIRVOL(I,J,L)

            ! VR_OC_BC  = volume ratio of OC to BC [unitless]
            VR_OC_BC    = C_OC / C_BC

            ! VR_BC_AIR = volume ratio of BC to air [unitless]
            VR_BC_AIR   = VR_OC_AIR / VR_OC_BC

            ! VR_BC_OC  = volume ratio of BC to OC [unitless]
            VR_BC_OC    = 1d0 / VR_OC_BC 


            ! Fraction of PBL that box (I,J,L) makes up [unitless]
            F_OF_PBL    = GET_FRAC_OF_PBL( I, J, L )


            ! Define temperature-dependent partition coefficients:
            ! KOA_T, the octanol-air partition coeff at temp T [unitless]
            KOA_T = KOA_298 * EXP((-DEL_H/R) * ((1d0/TK) - (1d0/298d0)))

            ! Define KBC_T, the BC-air partition coeff at temp T [unitless]
            KBC_T = KBC_298 * EXP((-DEL_H/R) * ((1d0/TK) - (1d0/298d0)))

            ! Define KOC_BC_T, the theoretical OC-BC part coeff at temp T [unitless]
            KOC_BC_T = KOA_T / KBC_T

            ! Define KBC_OC_T, the theoretical BC_OC part coeff at temp T [unitless]
            KBC_OC_T = 1d0 / KOC_BC_T


            ! Define fractions of total emissions that are OC-phase, 
            ! BC-phase, and gas phase 
            ! OC-phase:
            F_POP_OC  = 1d0 / (1d0 + (1d0 / (KOA_T * VR_OC_AIR)) + 
     &                  (1d0 / (KOC_BC_T * VR_OC_BC)))
            
            ! BC-phase:
            F_POP_BC  = 1d0 / (1d0 + (1d0 / (KBC_T * VR_BC_AIR)) + 
     &                  (1d0 / (KBC_OC_T * VR_BC_OC)))

            ! Gas-phase:
            F_POP_G   = 1d0 - F_POP_OC - F_POP_BC


            ! Calculate rates of POP emissions in each phase [kg/s]
            ! OC-phase:
            EPOP_OC    = F_POP_OC * T_POP

            ! BC-phase
            EPOP_BC    = F_POP_BC * T_POP

            ! Gas-phase
            EPOP_G     = F_POP_G * T_POP

            ! Total particulate
            EPOP_P_TOT = EPOP_OC + EPOP_BC
            

           ! Removed IDTPOPTOT as a tracer (clf, 2/14/2011)
           ! !-----------------
           ! ! TOTAL POP EMISSIONS (gas plus particulate emission)
           ! !-----------------
           ! N           = IDTPOPTOT
           ! E_POP       = F_OF_PBL * T_POP * DTSRCE
           ! CALL EMITPOP( I, J, L, N, E_POP )

            !-----------------
            ! OC-PHASE EMISSIONS
            !-----------------
            N           = IDTPOPPOC
            E_POP       = F_OF_PBL * EPOP_OC(I,J,L) * DTSRCE
            CALL EMITPOP( I, J, L, N, E_POP )

            !-----------------
            ! BC-PHASE EMISSIONS
            !-----------------
            N           = IDTPOPPBC
            E_POP       = F_OF_PBL * EPOP_BC(I,J,L) * DTSRCE
            CALL EMITPOP( I, J, L, N, E_POP )

            !-----------------
            ! GASEOUS EMISSIONS
            !-----------------
            N           = IDTPOPG
            E_POP       = F_OF_PBL * EPOP_G(I,J,L) * DTSRCE
            CALL EMITPOP( I, J, L, N, E_POP )
         ENDDO

         !==============================================================
         ! Sum different POPs emissions phases (OC, BC, and gas phase)
         ! through bottom layer to top of PBL for storage in ND53 diagnostic
         !==============================================================

         SUM_OC_EM = SUM ( EPOP_OC( I, J, 1:PBL_MAX ) ) 
         SUM_BC_EM = SUM ( EPOP_BC( I, J, 1:PBL_MAX ) )
         SUM_G_EM  = SUM ( EPOP_G ( I, J, 1:PBL_MAX ) )

         !==============================================================
         ! ND53 diagnostic: POP emissions [kg]
         ! 1=total; 2=OC; 3=BC; 4=gas phase
         !==============================================================
         IF ( ND53 > 0 ) THEN
            AD53(I,J,1) = AD53(I,J,1) + (T_POP     * DTSRCE)
            AD53(I,J,2) = AD53(I,J,2) + (SUM_OC_EM * DTSRCE)
            AD53(I,J,3) = AD53(I,J,3) + (SUM_BC_EM * DTSRCE)
            AD53(I,J,4) = AD53(I,J,4) + (SUM_G_EM  * DTSRCE)
         ENDIF

         
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE EMISSPOPS

!EOC
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  EMITPOP
!
! !DESCRIPTION: This routine directs emission either to STT directly or to EMIS_SAVE
!  for use by the non-local PBL mixing. This is a programming convenience.
!  (cdh, 08/27/09, modified for pops by eck, 9/20/10)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE EMITPOP( I, J, L, ID, E_POP )
!! 
! !USES:
      ! Reference to diagnostic arrays
      USE TRACER_MOD,   ONLY : STT
      USE LOGICAL_MOD,  ONLY : LNLPBL
      USE VDIFF_PRE_MOD,ONLY : EMIS_SAVE
! !INPUT PARAMETERS: 
      INTEGER, INTENT(IN)   :: I, J, L, ID
      REAL*8,  INTENT(IN)   :: E_POP
!
!
! !INPUT/OUTPUT PARAMETERS: 
!
!
! !OUTPUT PARAMETERS:
!
!
! !REVISION HISTORY: 
!  20 September 2010 - N.E. Selin - Initial Version
!
! !REMARKS:
! (1) Based initially on EMITHG from MERCURY_MOD (eck, 9/20/10)
!
!EOP
!******************************************************************************
!Comment header
!  Subroutine EMITPOP directs emission either to STT directly or to EMIS_SAVE
!  for use by the non-local PBL mixing. This is a programming convenience.
!  (cdh, 08/27/09, modified for pops by eck, 9/20/10)
!  
!  Arguments as Input:
!  ============================================================================
!  (1 ) I, J, L            INTEGERS  Grid box dimensions
!  (2 ) ID                 INTEGER   Tracer ID
!  (3 ) E_POP              REAL*8    POP emissions [kg/s]
!
!  Local variables:
!  ============================================================================
!  (1 ) 
!     
!  NOTES:
!  (1 ) Based on EMITHG in mercury_mod.f
!  
!  REFS:
!  (1 )

!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! EMITPOP begins here!
      !=================================================================

      ! Save emissions [kg/s] for non-local PBL mixing or emit directly.
      ! Make sure that emitted mass is non-negative
      ! This is here only for consistency with old code which warned of
      ! underflow error (cdh, 08/27/09, modified for POPs 9/20/10)
      IF (LNLPBL) THEN
         EMIS_SAVE(I,J,ID) = EMIS_SAVE(I,J,ID) + MAX( E_POP, 0D0 )
      ELSE
         STT(I,J,L,ID) = STT(I,J,L,ID) + MAX( E_POP, 0D0 )
      ENDIF

      END SUBROUTINE EMITPOP
!EOC
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  POPS_READYR
!
! !DESCRIPTION: 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE POPS_READYR
!! 
! !USES:
      ! References to F90 modules
      USE BPCH2_MOD,         ONLY : READ_BPCH2, GET_TAU0
      USE DIRECTORY_MOD,     ONLY : DATA_DIR_1x1
    
! !INPUT PARAMETERS: 
!
!
! !INPUT/OUTPUT PARAMETERS: 
!
!
! !OUTPUT PARAMETERS:
!
!
! !REVISION HISTORY: 
!  3 February 2011 - CL Friedman - Initial Version
!
! !REMARKS:
! (1) Based initially on MERCURY_READYR from MERCURY_MOD (clf, 2/3/2011)
!
!EOP
!******************************************************************************
!Comment header
!  Subroutine POP_READYR read the year-invariant emissions for POPs (PAHs) from 
!  all sources combined
!  
!  Arguments as Input:
!  ============================================================================
!  (1 ) 
!
!  Local variables:
!  ============================================================================
!  (1 ) 
!     
!  NOTES:
!  (1 ) Based on MERCURY_READYR in mercury_mod.f
!  
!  REFS:
!  (1 ) Zhang, Y. and Tao, S. 2009. Global atmospheric emission inventory
!  of polycyclic aromatic hydrocarbons (PAHs) for 2004. Atm Env. 43:812-819.

!------------------------------------------------------------------------------
!BOC
#     include "CMN_SIZE"       ! Size parameters

      ! Local variables
      !REAL*4               :: POP_TOT_EM(72,46,1)    
      REAL*8               :: XTAU
      REAL*8, PARAMETER    :: SEC_PER_YR = 365.35d0 * 86400d0  
      CHARACTER(LEN=225)   :: FILENAME 

      !=================================================================
      ! POP_READYR begins here!
      !=================================================================

      
      ! POLYCYCLIC AROMATIC HYDROCARBONS (PAHS):
      ! PAH emissions are for the year 2004
      ! Each PAH congener is emitted individually and contained in separate
      ! files
 
      ! Filename for congener you wish to model:
      !FILENAME = TRIM( DATA_DIR_1x1 )       // 
!     &           'PAHs_2004/PHE_EM_4x5.bpch' 
      FILENAME = '/net/fs05/d1/geosdata/data/GEOS_4x5/PAHs_2004/' //
     &           'PHE_EM_4x5.bpch'

      
      ! Timestamp for emissions
      ! All PAH emissions are for the year 2004
      XTAU = GET_TAU0( 1, 1, 2004)

      ! Add year to the filename (?)
      ! CALL EXPAND_DATE( FILENAME, NYMD, 000000 )

      ! Echo info (?)
      ! WRITE( 6, 100 ) TRIM( FILENAME )
! 100        FORMAT( '     - MERCURY_READYR: Reading ', a )

       ! Read data in [1000 kg/yr]
       CALL READ_BPCH2( FILENAME, 'PG-SRCE', 1, 
     &           XTAU,      72,     46,    
     &           1,         POP_TOT_EM,   QUIET=.FALSE. )


      WRITE ( 6, '(a)') 'Yikes'
      STOP

     

       ! Convert from [1000 kg/yr] to [kg/s]
       POP_TOT_EM = POP_TOT_EM / 1000d0 / SEC_PER_YR

      !=================================================================
      ! Print totals to the screen in [Gg/yr]
      !=================================================================
      ! Need to figure out how to code this
      
      ! Return to calling program
      END SUBROUTINE POPS_READYR

!EOC
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  GET_OH
!
! !DESCRIPTION: Function GET_OH returns monthly mean OH and imposes a diurnal
! variation. 
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_OH( I, J, L ) RESULT( OH_MOLEC_CM3 )

      ! References to F90 modules
      USE DAO_MOD,       ONLY : SUNCOS 
      USE GLOBAL_OH_MOD, ONLY : OH
      USE TIME_MOD,      ONLY : GET_TS_CHEM

#     include "CMN_SIZE"  ! Size parameters

! !INPUT PARAMETERS: 

      INTEGER, INTENT(IN) :: I, J, L
!
! !OUTPUT PARAMETERS:
!
!
! !REVISION HISTORY: 
!  03 February 2011 - CL Friedman - Initial Version
!
! !REMARKS:
! Copied GET_OH function from mercury_mod.f - CLF
!
!EOP
!------------------------------------------------------------------------------
!BOC

       ! Local variables
       INTEGER       :: JLOOP
       REAL*8        :: OH_MOLEC_CM3
    
       !=================================================================
       ! GET_OH begins here!
       !=================================================================

       ! 1-D grid box index for SUNCOS
       JLOOP = ( (J-1) * IIPAR ) + I

       ! Test for sunlight...
       IF ( SUNCOS(JLOOP) > 0d0 .and. TCOSZ(I,J) > 0d0 ) THEN

         ! Impose a diurnal variation on OH during the day
         OH_MOLEC_CM3 = OH(I,J,L)                      *           
     &                  ( SUNCOS(JLOOP) / TCOSZ(I,J) ) *
     &                  ( 1440d0        / GET_TS_CHEM() )

         ! Make sure OH is not negative
         OH_MOLEC_CM3 = MAX( OH_MOLEC_CM3, 0d0 )
               
       ELSE

         ! At night, OH goes to zero
         OH_MOLEC_CM3 = 0d0

       ENDIF

       ! Return to calling program
       END FUNCTION GET_OH
!EOC

!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  GET_OC
!
! !DESCRIPTION: Function GET_OC returns monthly mean organic carbon 
! concentrations [kg/box]
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_OC( I, J, L) RESULT( C_OC )
!
! !INPUT PARAMETERS:

!     References to F90 modules
      USE GLOBAL_OC_MOD, ONLY : OC 

      INTEGER, INTENT(IN) :: I, J, L 
!
! !OUTPUT PARAMETERS:
!
!
! !REVISION HISTORY: 
!  03 February 2011 - CL Friedman - Initial Version
!
! !REMARKS:
! Test
!
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Local variables
      REAL*8            :: C_OC

      !=================================================================
      ! GET_OC begins here!
      !=================================================================

      ! Get organic carbon concentration [kg/box] for this gridbox and month
      C_OC = OC(I,J,L)

      ! Return to calling program
      END FUNCTION GET_OC
!EOC

!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  GET_BC
!
! !DESCRIPTION: Function GET_BC returns monthly mean black carbon concentrations
! [kg/box]
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_BC( I, J, L) RESULT( C_BC )

!
! !INPUT PARAMETERS: 
!
      ! References to F90 modules
      USE GLOBAL_BC_MOD, ONLY : BC
   
      INTEGER, INTENT(IN) :: I, J, L   
!
! !OUTPUT PARAMETERS:
!
!
! !REVISION HISTORY: 
!  03 February 2011 - CL Friedman - Initial Version
!
! !REMARKS:
! Test
!
!EOP
!------------------------------------------------------------------------------
!BOC

      ! Local variables
      REAL*8      :: C_BC    
    
      !=================================================================
      ! GET_BC begins here!
      !=================================================================

      ! Get black carbon concentration [kg/box] for this gridbox and month
      C_BC = BC(I,J,L)

      ! Return to calling program

      END FUNCTION GET_BC
!EOC

!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  OHNO3TIME
!
! !DESCRIPTION: Subroutine OHNO3TIME computes the sum of cosine of the solar zenith
!  angle over a 24 hour day, as well as the total length of daylight. 
!  This is needed to scale the offline OH and NO3 concentrations.
!  (rjp, bmy, 12/16/02, 12/8/04)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE OHNO3TIME
!! 
! !USES:
      ! References to F90 modules
      USE GRID_MOD, ONLY : GET_XMID,    GET_YMID_R
      USE TIME_MOD, ONLY : GET_NHMSb,   GET_ELAPSED_SEC
      USE TIME_MOD, ONLY : GET_TS_CHEM, GET_DAY_OF_YEAR, GET_GMT


#     include "CMN_SIZE"  ! Size parameters
#     include "CMN_GCTM"  ! Physical constants
! !INPUT PARAMETERS: 
!
!
! !INPUT/OUTPUT PARAMETERS: 
!
!
! !OUTPUT PARAMETERS:
!
!
! !REVISION HISTORY: 
!  20 September 2010 - N.E. Selin - Initial Version for POPS_MOD
!
! !REMARKS:
!  (1 ) Copy code from COSSZA directly for now, so that we don't get NaN
!        values.  Figure this out later (rjp, bmy, 1/10/03)
!  (2 ) Now replace XMID(I) with routine GET_XMID from "grid_mod.f".  
!        Now replace RLAT(J) with routine GET_YMID_R from "grid_mod.f". 
!        Removed NTIME, NHMSb from the arg list.  Now use GET_NHMSb,
!        GET_ELAPSED_SEC, GET_TS_CHEM, GET_DAY_OF_YEAR, GET_GMT from 
!        "time_mod.f". (bmy, 3/27/03)
!  (3 ) Now store the peak SUNCOS value for each surface grid box (I,J) in 
!        the COSZM array. (rjp, bmy, 3/30/04)
!  (4 ) Also added parallel loop over grid boxes (eck, bmy, 12/8/04)
!  (5 ) copied from mercury_mod by eck (9/20/10)
!******************************************************************************
!
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Local variables
      LOGICAL, SAVE       :: FIRST = .TRUE.
      INTEGER             :: I, IJLOOP, J, L, N, NT, NDYSTEP
      REAL*8              :: A0, A1, A2, A3, B1, B2, B3
      REAL*8              :: LHR0, R, AHR, DEC, TIMLOC, YMID_R
      REAL*8              :: SUNTMP(MAXIJ)
      
      !=================================================================
      ! OHNO3TIME begins here!
      !=================================================================

      !  Solar declination angle (low precision formula, good enough for us):
      A0 = 0.006918
      A1 = 0.399912
      A2 = 0.006758
      A3 = 0.002697
      B1 = 0.070257
      B2 = 0.000907
      B3 = 0.000148
      R  = 2.* PI * float( GET_DAY_OF_YEAR() - 1 ) / 365.

      DEC = A0 - A1*cos(  R) + B1*sin(  R)
     &         - A2*cos(2*R) + B2*sin(2*R)
     &         - A3*cos(3*R) + B3*sin(3*R)

      LHR0 = int(float( GET_NHMSb() )/10000.)

      ! Only do the following at the start of a new day
      IF ( FIRST .or. GET_GMT() < 1e-5 ) THEN 
      
         ! Zero arrays
         TTDAY(:,:) = 0d0
         TCOSZ(:,:) = 0d0
         COSZM(:,:) = 0d0

         ! NDYSTEP is # of chemistry time steps in this day
         NDYSTEP = ( 24 - INT( GET_GMT() ) ) * 60 / GET_TS_CHEM()         

         ! NT is the elapsed time [s] since the beginning of the run
         NT = GET_ELAPSED_SEC()

         ! Loop forward through NDYSTEP "fake" timesteps for this day 
         DO N = 1, NDYSTEP
            
            ! Zero SUNTMP array
            SUNTMP(:) = 0d0

            ! Loop over surface grid boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, YMID_R, IJLOOP, TIMLOC, AHR )
            DO J = 1, JJPAR

               ! Grid box latitude center [radians]
               YMID_R = GET_YMID_R( J )

            DO I = 1, IIPAR

               ! Increment IJLOOP
               IJLOOP = ( (J-1) * IIPAR ) + I
               TIMLOC = real(LHR0) + real(NT)/3600.0 + GET_XMID(I)/15.0
         
               DO WHILE (TIMLOC .lt. 0)
                  TIMLOC = TIMLOC + 24.0
               ENDDO

               DO WHILE (TIMLOC .gt. 24.0)
                  TIMLOC = TIMLOC - 24.0
               ENDDO

               AHR = abs(TIMLOC - 12.) * 15.0 * PI_180

            !===========================================================
            ! The cosine of the solar zenith angle (SZA) is given by:
            !     
            !  cos(SZA) = sin(LAT)*sin(DEC) + cos(LAT)*cos(DEC)*cos(AHR) 
            !                   
            ! where LAT = the latitude angle, 
            !       DEC = the solar declination angle,  
            !       AHR = the hour angle, all in radians. 
            !
            ! If SUNCOS < 0, then the sun is below the horizon, and 
            ! therefore does not contribute to any solar heating.  
            !===========================================================

               ! Compute Cos(SZA)
               SUNTMP(IJLOOP) = sin(YMID_R) * sin(DEC) +
     &                          cos(YMID_R) * cos(DEC) * cos(AHR)

               ! TCOSZ is the sum of SUNTMP at location (I,J)
               ! Do not include negative values of SUNTMP
               TCOSZ(I,J) = TCOSZ(I,J) + MAX( SUNTMP(IJLOOP), 0d0 )

               ! COSZM is the peak value of SUMTMP during a day at (I,J)
               ! (rjp, bmy, 3/30/04)
               COSZM(I,J) = MAX( COSZM(I,J), SUNTMP(IJLOOP) )

               ! TTDAY is the total daylight time at location (I,J)
               IF ( SUNTMP(IJLOOP) > 0d0 ) THEN
                  TTDAY(I,J) = TTDAY(I,J) + DBLE( GET_TS_CHEM() )
               ENDIF
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

            ! Increment elapsed time [sec]
            NT = NT + ( GET_TS_CHEM() * 60 )             
         ENDDO

         ! Reset first-time flag
         FIRST = .FALSE.
      ENDIF

      ! Return to calling program
      END SUBROUTINE OHNO3TIME

!EOC
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  INIT_POPS
!
! !DESCRIPTION: Subroutine INIT_POPS allocates and zeroes all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_POPS
!
      ! References to F90 modules
      USE DRYDEP_MOD,   ONLY : DEPNAME,   NUMDEP
      USE ERROR_MOD,    ONLY : ALLOC_ERR, ERROR_STOP
      USE LOGICAL_MOD,  ONLY : LSPLIT,    LDRYD,     LNLPBL
      USE TRACER_MOD,   ONLY : N_TRACERS
      USE PBL_MIX_MOD,  ONLY : GET_PBL_MAX_L
c$$$
#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_DIAG"     ! ND44

! !INPUT PARAMETERS: 
!
!
! !INPUT/OUTPUT PARAMETERS: 
!
!
! !OUTPUT PARAMETERS:
!
!
! !REVISION HISTORY: 
!  20 September 2010 - N.E. Selin - Initial Version
!
! !REMARKS:
! (1) Based initially on INIT_MERCURY from MERCURY_MOD (eck, 9/20/10)
!
!EOP
!------------------------------------------------------------------------------
!BOC

      ! Local variables
      LOGICAL, SAVE         :: IS_INIT = .FALSE. 
      INTEGER               :: AS, N!, PBL_MAX
      !=================================================================
      ! INIT_POPS begins here!
      !=================================================================

      ! Maximum extent of the PBL
      !PBL_MAX = GET_PBL_MAX_L()

      ! Return if we have already allocated arrays
      IF ( IS_INIT ) RETURN

      !=================================================================
      ! Allocate and initialize arrays
      !=================================================================
      ALLOCATE( COSZM( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'COSZM' )
      COSZM = 0d0

      ALLOCATE( TCOSZ( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TCOSZ' )
      TCOSZ = 0d0

      ALLOCATE( TTDAY( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TTDAY' )
      TTDAY = 0d0

      ALLOCATE( EPOP_G( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPOP_G' )
      EPOP_G = 0d0

      ALLOCATE( EPOP_OC( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPOP_OC' )
      EPOP_OC = 0d0

      ALLOCATE( EPOP_BC( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPOP_BC' )
      EPOP_BC = 0d0

      ALLOCATE( EPOP_P_TOT( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPOP_P_TOT' )
      EPOP_P_TOT = 0d0

      ALLOCATE( EPOP_P_TOT( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EPOP_P_TOT' )
      EPOP_P_TOT = 0d0

      ALLOCATE( POP_TOT_EM( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'POP_TOT_EM' )
      POP_TOT_EM = 0d0

      ! Allocate ZERO_DVEL if we use non-local PBL mixing or
      ! if drydep is turned off 
      IF ( LNLPBL .OR. (.not. LDRYD) ) THEN
         ALLOCATE( ZERO_DVEL( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'ZERO_DVEL' )
         ZERO_DVEL = 0d0
      ENDIF

      !=================================================================
      ! Done
      !=================================================================

      ! Reset IS_INIT, since we have already allocated arrays
      IS_INIT = .TRUE.
      
      ! Return to calling program
      END SUBROUTINE INIT_POPS

!EOC
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  CLEANUP_POPS
!
! !DESCRIPTION: Subroutine CLEANUP_POPS deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_POPS
!

! !INPUT PARAMETERS: 
!
!
! !INPUT/OUTPUT PARAMETERS: 
!
!
! !OUTPUT PARAMETERS:
!
!
! !REVISION HISTORY: 
!  20 September 2010 - N.E. Selin - Initial Version
!
! !REMARKS:
! (1) Based initially on INIT_MERCURY from MERCURY_MOD (eck, 9/20/10)
!
!EOP
!------------------------------------------------------------------------------
!BOC

      IF ( ALLOCATED( COSZM    ) ) DEALLOCATE( COSZM   )     
      IF ( ALLOCATED( TCOSZ    ) ) DEALLOCATE( TCOSZ   )
      IF ( ALLOCATED( TTDAY    ) ) DEALLOCATE( TTDAY   )
      IF ( ALLOCATED( ZERO_DVEL) ) DEALLOCATE( ZERO_DVEL )
      IF ( ALLOCATED( EPOP_G   ) ) DEALLOCATE( EPOP_G   )
      IF ( ALLOCATED( EPOP_OC  ) ) DEALLOCATE( EPOP_OC  )
      IF ( ALLOCATED( EPOP_BC  ) ) DEALLOCATE( EPOP_P_TOT )
      IF ( ALLOCATED( POP_TOT_EM ) ) DEALLOCATE( POP_TOT_EM )

      END SUBROUTINE CLEANUP_POPS
!EOC
!------------------------------------------------------------------------------
      END MODULE POPS_MOD

