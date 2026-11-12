module Reaction_Sandbox_Mineral_Armoring_class

#include "petsc/finclude/petscsys.h"
  use petscsys

  use Reaction_Sandbox_Base_class
  use PFLOTRAN_Constants_module

  implicit none

  private

! This sandbox couples the kinetic surface area of one mineral (e.g. a
! sulfide such as Pyrrhotite) to the extent of dissolution of a second,
! "armoring" mineral (e.g. Calcite) occupying the same grid cell:
!
!   f_dissolved = (volfrac0_armor - volfrac_armor) / volfrac0_armor
!   area = area0 * f_dissolved ** SURFACE_AREA_POWER
!
! where f_dissolved is clipped to [0,1] and area0 is the initial surface
! area of MINERAL (rt_auxvar%mnrl_area0). f_dissolved = 0 when the armoring
! mineral is unaltered (MINERAL is fully armored/unreactive) and 1 once the
! armoring mineral has been completely dissolved (MINERAL exposed at its
! full initial surface area).
!
! Optionally, the surface area of MINERAL is also reduced as MINERAL itself
! is consumed by dissolution (shrinking-particle behavior):
!
!   area = area * (volfrac / volfrac0) ** DEPLETION_POWER
!
! where volfrac/volfrac0 are the current/initial volume fractions of
! MINERAL in the grid cell. This factor is applied to dissolution only, so
! that the rate goes to zero smoothly as MINERAL is exhausted. The default
! DEPLETION_POWER of 0 disables the factor (area depends only on armoring).
!
! Regardless of DEPLETION_POWER, the dissolution rate is capped so that a
! single transport timestep cannot dissolve more MINERAL than remains in
! the grid cell.
!
! The dissolution/precipitation rate for MINERAL is otherwise computed with
! a standard single-term TST-like rate law leveraging the stoichiometry and
! logK for MINERAL already present in the thermodynamic database, with an
! optional single-species power-law prefactor and Arrhenius temperature
! correction:
!
!   rate = -area * sign(1-QK) * |1-QK| * RATE_CONSTANT * &
!           {PREFACTOR_SPECIES activity}**PREFACTOR_ALPHA * arrhenius_factor
!
! Optional low-concentration cutoff on the prefactor species (e.g. O2(aq)):
!
!   rate = rate * m / (PREFACTOR_HALF_SATURATION + m)
!
! where m is the molality of PREFACTOR_SPECIES. With a half-saturation
! constant far below normal concentrations (e.g. 1e-9 M for O2) the factor
! is ~1 wherever the species is present, but the rate goes to zero smoothly
! (order alpha+1) at trace levels, so a sub-linear sink cannot drive the
! species toward zero in anoxic cells (log-formulation underflow).
! Default 0 disables the factor.
!
! NO_ARMORING sets f_dissolved = 1 everywhere (ARMORING_MINERAL is then not
! needed). With DEPLETION_POWER 2/3 this reproduces the standard continuum
! rate law (MINERAL_KINETICS with SURFACE_AREA_FUNCTION VOLUME_FRACTION_RATIO,
! power 2/3), so armored and continuum runs can share one rate law and one
! O2 cutoff and differ only in the armoring factor.
!
! For REACTION_SANDBOX,MINERAL_ARMORING to function correctly, the
! corresponding entry for MINERAL within the default MINERAL_KINETICS block
! must have its RATE_CONSTANT (and any PREFACTOR RATE_CONSTANT) set to zero,
! as this sandbox fully replaces that default kinetic rate law. The default
! MINERAL_KINETICS entry for MINERAL is still required (to supply the
! molar volume, logK/stoichiometry from the database, and the initial
! surface area/volume fraction read from the CONSTRAINT).
!
! No analytical Jacobian is provided (the prefactor and Arrhenius terms
! make hand-derivation error-prone); NUMERICAL_JACOBIAN must be specified
! within NUMERICAL_METHODS TRANSPORT, NEWTON_SOLVER when this sandbox is
! used.
!
! Example Input:
!
! CHEMISTRY
!   ...
!   REACTION_SANDBOX
!     MINERAL_ARMORING
!       MINERAL Pyrrhotite
!       ARMORING_MINERAL Calcite
!       RATE_CONSTANT 2.81d-5 mol/m^2-sec
!       PREFACTOR_SPECIES O2(aq)
!       PREFACTOR_ALPHA 0.5
!       ACTIVATION_ENERGY 5.d4 J/mol
!       SURFACE_AREA_POWER 0.66666667
!       DEPLETION_POWER 0.66666667
!     END
!   END
!   ...
! END

  type, public, &
    extends(reaction_sandbox_base_type) :: reaction_sandbox_mineral_armoring_type
    character(len=MAXWORDLENGTH) :: mineral_name
    character(len=MAXWORDLENGTH) :: armor_mineral_name
    character(len=MAXWORDLENGTH) :: prefactor_species_name
    PetscInt :: auxiliary_offset
    PetscInt :: mineral_id
    PetscInt :: armor_mineral_id
    PetscInt :: prefactor_species_id
    PetscReal :: rate_constant
    PetscReal :: activation_energy
    PetscReal :: prefactor_alpha
    PetscReal :: surface_area_power
    PetscReal :: depletion_power
    PetscReal :: affinity_threshold
    PetscReal :: prefactor_half_sat
    PetscBool :: no_armoring
  contains
    procedure, public :: ReadInput => MineralArmoringReadInput
    procedure, public :: Setup => MineralArmoringSetup
    procedure, public :: AuxiliaryPlotVariables => &
                           MineralArmoringAuxiliaryPlotVariables
    procedure, public :: Evaluate => MineralArmoringEvaluate
    procedure, public :: UpdateKineticState => MineralArmoringUpdateKineticState
  end type reaction_sandbox_mineral_armoring_type

  public :: MineralArmoringCreate

contains

! ************************************************************************** !

function MineralArmoringCreate()
  !
  ! Allocates mineral armoring reaction object.
  !
  ! Author: Aaron Sitchler
  ! Date: 08/24/26
  !
  implicit none

  class(reaction_sandbox_mineral_armoring_type), pointer :: MineralArmoringCreate

  allocate(MineralArmoringCreate)
  MineralArmoringCreate%mineral_name = ''
  MineralArmoringCreate%armor_mineral_name = ''
  MineralArmoringCreate%prefactor_species_name = ''
  MineralArmoringCreate%auxiliary_offset = UNINITIALIZED_INTEGER
  MineralArmoringCreate%mineral_id = UNINITIALIZED_INTEGER
  MineralArmoringCreate%armor_mineral_id = UNINITIALIZED_INTEGER
  MineralArmoringCreate%prefactor_species_id = 0
  MineralArmoringCreate%rate_constant = UNINITIALIZED_DOUBLE
  MineralArmoringCreate%activation_energy = 0.d0
  MineralArmoringCreate%prefactor_alpha = UNINITIALIZED_DOUBLE
  MineralArmoringCreate%surface_area_power = UNINITIALIZED_DOUBLE
  MineralArmoringCreate%depletion_power = 0.d0
  MineralArmoringCreate%affinity_threshold = 0.d0
  MineralArmoringCreate%prefactor_half_sat = 0.d0
  MineralArmoringCreate%no_armoring = PETSC_FALSE
  nullify(MineralArmoringCreate%next)

end function MineralArmoringCreate

! ************************************************************************** !

subroutine MineralArmoringReadInput(this,input,option)
  !
  ! Reads mineral armoring reaction parameters
  !
  ! Author: Aaron Sitchler
  ! Date: 08/24/26
  !
  use Option_module
  use Input_Aux_module
  use String_module

  implicit none

  class(reaction_sandbox_mineral_armoring_type) :: this
  type(input_type), pointer :: input
  type(option_type) :: option

  character(len=MAXWORDLENGTH) :: word
  character(len=MAXSTRINGLENGTH) :: error_string

  error_string = 'CHEMISTRY,REACTION_SANDBOX,MINERAL_ARMORING'
  call InputPushBlock(input,option)
  do
    call InputReadPflotranString(input,option)
    if (InputError(input)) exit
    if (InputCheckExit(input,option)) exit

    call InputReadCard(input,option,word)
    call InputErrorMsg(input,option,'keyword',error_string)
    call StringToUpper(word)

    select case(word)
      case('MINERAL')
        call InputReadWord(input,option,this%mineral_name,PETSC_TRUE)
        call InputErrorMsg(input,option,word,error_string)
      case('ARMORING_MINERAL')
        call InputReadWord(input,option,this%armor_mineral_name,PETSC_TRUE)
        call InputErrorMsg(input,option,word,error_string)
      case('RATE_CONSTANT')
        call InputReadDouble(input,option,this%rate_constant)
        call InputErrorMsg(input,option,word,error_string)
        call InputReadAndConvertUnits(input,this%rate_constant, &
                                      'mol/m^2-sec', &
                                      trim(error_string)//',RATE_CONSTANT', &
                                      option)
      case('ACTIVATION_ENERGY')
        call InputReadDouble(input,option,this%activation_energy)
        call InputErrorMsg(input,option,word,error_string)
        call InputReadAndConvertUnits(input,this%activation_energy,'J/mol', &
                              trim(error_string)//',ACTIVATION_ENERGY',option)
      case('PREFACTOR_SPECIES')
        call InputReadWord(input,option,this%prefactor_species_name, &
                           PETSC_TRUE)
        call InputErrorMsg(input,option,word,error_string)
      case('PREFACTOR_ALPHA')
        call InputReadDouble(input,option,this%prefactor_alpha)
        call InputErrorMsg(input,option,word,error_string)
      case('SURFACE_AREA_POWER')
        call InputReadDouble(input,option,this%surface_area_power)
        call InputErrorMsg(input,option,word,error_string)
      case('DEPLETION_POWER')
        call InputReadDouble(input,option,this%depletion_power)
        call InputErrorMsg(input,option,word,error_string)
      case('AFFINITY_THRESHOLD')
        call InputReadDouble(input,option,this%affinity_threshold)
        call InputErrorMsg(input,option,word,error_string)
      case('PREFACTOR_HALF_SATURATION')
        call InputReadDouble(input,option,this%prefactor_half_sat)
        call InputErrorMsg(input,option,word,error_string)
      case('NO_ARMORING')
        this%no_armoring = PETSC_TRUE
      case default
        call InputKeywordUnrecognized(input,word,error_string,option)
    end select
  enddo
  call InputPopBlock(input,option)

  if (len_trim(this%mineral_name) == 0 .or. &
      (len_trim(this%armor_mineral_name) == 0 .and. &
       .not.this%no_armoring)) then
    option%io_buffer = 'MINERAL and ARMORING_MINERAL must be set for &
      &REACTION_SANDBOX,MINERAL_ARMORING.'
    call PrintErrMsg(option)
  endif
  if (.not.this%no_armoring .and. &
      StringCompareIgnoreCase(this%mineral_name, &
                              this%armor_mineral_name)) then
    option%io_buffer = 'MINERAL and ARMORING_MINERAL must differ for &
      &REACTION_SANDBOX,MINERAL_ARMORING ('//trim(this%mineral_name)// &
      &'), as a mineral armoring itself would never react.'
    call PrintErrMsg(option)
  endif
  if (Uninitialized(this%rate_constant) .or. &
      Uninitialized(this%surface_area_power)) then
    option%io_buffer = 'RATE_CONSTANT and SURFACE_AREA_POWER must be set &
      &for REACTION_SANDBOX,MINERAL_ARMORING.'
    call PrintErrMsg(option)
  endif
  if (len_trim(this%prefactor_species_name) > 0 .and. &
      Uninitialized(this%prefactor_alpha)) then
    option%io_buffer = 'PREFACTOR_ALPHA must be set when PREFACTOR_SPECIES &
      &is specified for REACTION_SANDBOX,MINERAL_ARMORING.'
    call PrintErrMsg(option)
  endif
  if (this%prefactor_half_sat < 0.d0 .or. &
      (this%prefactor_half_sat > 0.d0 .and. &
       len_trim(this%prefactor_species_name) == 0)) then
    option%io_buffer = 'PREFACTOR_HALF_SATURATION must be non-negative and &
      &requires PREFACTOR_SPECIES for REACTION_SANDBOX,MINERAL_ARMORING.'
    call PrintErrMsg(option)
  endif
  if (this%depletion_power < 0.d0) then
    option%io_buffer = 'DEPLETION_POWER must be non-negative for &
      &REACTION_SANDBOX,MINERAL_ARMORING.'
    call PrintErrMsg(option)
  endif

end subroutine MineralArmoringReadInput

! ************************************************************************** !

subroutine MineralArmoringSetup(this,reaction,option)
  !
  ! Sets up the mineral armoring reaction with parameters read from input
  !
  ! Author: Aaron Sitchler
  ! Date: 08/24/26
  !
  use Reaction_Aux_module
  use Reaction_Mineral_Aux_module
  use Option_module

  implicit none

  class(reaction_sandbox_mineral_armoring_type) :: this
  class(reaction_rt_type) :: reaction
  type(option_type) :: option

  ! rt_auxvar%auxiliary_data(:) is allocated to reaction%nauxiliary
  ! the offset points this sandbox to the correct entries for storing the
  ! rate and surface area for plotting
  this%auxiliary_offset = reaction%nauxiliary
  reaction%nauxiliary = reaction%nauxiliary + 2

  this%mineral_id = &
    ReactionMnrlGetKinMnrlIDFromName(this%mineral_name,reaction%mineral, &
                                     option)
  if (len_trim(this%armor_mineral_name) > 0) then
    this%armor_mineral_id = &
      ReactionMnrlGetKinMnrlIDFromName(this%armor_mineral_name, &
                                       reaction%mineral,option)
  endif
  if (len_trim(this%prefactor_species_name) > 0) then
    this%prefactor_species_id = &
      ReactionAuxGetPriSpecIDFromName(this%prefactor_species_name,reaction, &
                                      option)
  endif

end subroutine MineralArmoringSetup

! ************************************************************************** !

subroutine MineralArmoringAuxiliaryPlotVariables(this,list,reaction,option)
  !
  ! Adds mineral armoring auxiliary plot variables to output list
  !
  ! Author: Aaron Sitchler
  ! Date: 08/24/26
  !
  use Option_module
  use Reaction_Aux_module
  use Output_Aux_module
  use Variables_module, only : REACTION_AUXILIARY

  class(reaction_sandbox_mineral_armoring_type) :: this
  type(output_variable_list_type), pointer :: list
  type(option_type) :: option
  class(reaction_rt_type) :: reaction

  character(len=MAXWORDLENGTH) :: word
  character(len=MAXWORDLENGTH) :: units

  word = trim(this%mineral_name)//' Armoring Sandbox Rate'
  units = 'mol/m^3-sec'
  call OutputVariableAddToList(list,word,OUTPUT_RATE,units, &
                                REACTION_AUXILIARY, &
                                this%auxiliary_offset+1)

  word = trim(this%mineral_name)//' Armoring Sandbox Area'
  units = 'm^2/m^3'
  call OutputVariableAddToList(list,word,OUTPUT_GENERIC,units, &
                                REACTION_AUXILIARY, &
                                this%auxiliary_offset+2)

end subroutine MineralArmoringAuxiliaryPlotVariables

! ************************************************************************** !

subroutine MineralArmoringEvaluate(this,Residual,Jacobian,compute_derivative, &
                                   rt_auxvar,global_auxvar,material_auxvar, &
                                   reaction,option)
  !
  ! Evaluates the mineral armoring reaction storing residual but no Jacobian
  !
  ! Author: Aaron Sitchler
  ! Date: 08/24/26
  !
  use Option_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Material_Aux_module
  use Reaction_Mineral_Aux_module
  use Utility_module, only : Arrhenius

  implicit none

  class(reaction_sandbox_mineral_armoring_type) :: this
  type(option_type) :: option
  class(reaction_rt_type) :: reaction
  PetscBool :: compute_derivative
  PetscReal :: Residual(reaction%ncomp) ! [mole / sec]
  PetscReal :: Jacobian(reaction%ncomp,reaction%ncomp)
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar

  PetscReal, parameter :: TREF = 25.d0
  type(mineral_type), pointer :: mineral

  PetscReal :: ln_conc(reaction%ncomp)
  PetscReal :: ln_act(reaction%ncomp)
  PetscReal :: volfrac0_armor, volfrac_armor, f_dissolved
  PetscReal :: volfrac0, volfrac
  PetscReal :: area
  PetscReal :: rate
  PetscReal :: lnQK, QK
  PetscReal :: affinity_factor, sign_
  PetscReal :: prefactor_term
  PetscReal :: arrhenius_factor

  PetscInt :: ncomp
  PetscInt :: i, icomp
  PetscInt :: imnrl, imnrl_armor
  PetscInt :: iauxiliary
  PetscBool :: calculate_rate

  mineral => reaction%mineral
  iauxiliary = this%auxiliary_offset + 1

  imnrl = this%mineral_id
  imnrl_armor = this%armor_mineral_id

  ! for REACTION_SANDBOX MINERAL_ARMORING to function correctly, the
  ! RATE_CONSTANT (and any PREFACTOR RATE_CONSTANT) in the default
  ! MINERAL_KINETICS block for MINERAL must be set to zero, as this sandbox
  ! fully replaces the default kinetic rate law
  if (dabs(rt_auxvar%mnrl_rate(imnrl)) > 1.d-40) then
    option%io_buffer = 'For REACTION_SANDBOX MINERAL_ARMORING to function &
      &correctly, the RATE_CONSTANT(S) in the default MINERAL_KINETICS &
      &block for '//trim(this%mineral_name)//' must be set to zero.'
    call PrintErrMsg(option)
  endif

  ln_conc = log(rt_auxvar%pri_molal)
  ln_act = ln_conc+log(rt_auxvar%pri_act_coef)

  ! fraction of the armoring mineral that has dissolved, normalized by the
  ! LOCAL initial volume fraction in this grid cell
  if (this%no_armoring) then
    ! continuum: MINERAL fully exposed from t = 0
    f_dissolved = 1.d0
  else if (rt_auxvar%mnrl_volfrac0(imnrl_armor) > 0.d0) then
    volfrac0_armor = rt_auxvar%mnrl_volfrac0(imnrl_armor)
    volfrac_armor = rt_auxvar%mnrl_volfrac(imnrl_armor)
    f_dissolved = (volfrac0_armor-volfrac_armor)/volfrac0_armor
    f_dissolved = max(0.d0,min(1.d0,f_dissolved))
  else
    f_dissolved = 0.d0
  endif

  ! m^2 mnrl/m^3 bulk
  area = rt_auxvar%mnrl_area0(imnrl)*f_dissolved**this%surface_area_power

  ! compute ion activity product from the database stoichiometry/logK for
  ! MINERAL
  lnQK = -mineral%kinmnrl_logK(imnrl)*LOG_TO_LN
  if (mineral%kinmnrlh2oid(imnrl) > 0) then
    lnQK = lnQK + mineral%kinmnrlh2ostoich(imnrl)*rt_auxvar%ln_act_h2o
  endif
  ncomp = mineral%kinmnrlspecid(0,imnrl)
  do i = 1, ncomp
    icomp = mineral%kinmnrlspecid(i,imnrl)
    lnQK = lnQK + mineral%kinmnrlstoich(i,imnrl)*ln_act(icomp)
  enddo
  QK = exp(lnQK)
  affinity_factor = 1.d0-QK
  sign_ = sign(1.d0,affinity_factor)

  rate = 0.d0
  ! only calculate rate if mineral is present or precipitation (sign<0)
  calculate_rate = rt_auxvar%mnrl_volfrac(imnrl) > 0 .or. sign_ < 0.d0
  ! optional supersaturation threshold, disabling precipitation below it
  if (calculate_rate .and. sign_ < 0.d0 .and. &
      this%affinity_threshold > 0.d0 .and. &
      QK < this%affinity_threshold) then
    calculate_rate = PETSC_FALSE
  endif

  if (calculate_rate) then
    prefactor_term = 1.d0
    if (this%prefactor_species_id > 0) then
      prefactor_term = exp(this%prefactor_alpha* &
                           ln_act(this%prefactor_species_id))
      ! optional smooth cutoff at trace prefactor-species concentration
      if (this%prefactor_half_sat > 0.d0) then
        prefactor_term = prefactor_term* &
          rt_auxvar%pri_molal(this%prefactor_species_id)/ &
          (this%prefactor_half_sat+ &
           rt_auxvar%pri_molal(this%prefactor_species_id))
      endif
    endif
    arrhenius_factor = 1.d0
    if (this%activation_energy > 0.d0) then
      arrhenius_factor = Arrhenius(this%activation_energy, &
                                   global_auxvar%temp,TREF)
    endif
    ! optional shrinking of MINERAL surface area as MINERAL is consumed
    ! (dissolution only)
    volfrac0 = rt_auxvar%mnrl_volfrac0(imnrl)
    volfrac = rt_auxvar%mnrl_volfrac(imnrl)
    if (sign_ > 0.d0 .and. this%depletion_power > 0.d0 .and. &
        volfrac0 > 0.d0) then
      area = area*(max(0.d0,volfrac)/volfrac0)**this%depletion_power
    endif
    ! mol/sec/m^3 bulk
    rate = -area*sign_*dabs(affinity_factor)*this%rate_constant* &
           prefactor_term*arrhenius_factor
    ! limit dissolution to the MINERAL remaining in the grid cell so that
    ! the final timestep cannot release more solute than the mineral holds
    ! (volume fraction is only updated after the timestep converges)
    if (sign_ > 0.d0 .and. option%tran_dt > 0.d0) then
      ! mol/sec/m^3 bulk
      rate = max(rate,-max(0.d0,volfrac)/ &
                 (mineral%kinmnrl_molar_vol(imnrl)*option%tran_dt))
    endif
  endif

  rt_auxvar%auxiliary_data(iauxiliary) = rate
  rt_auxvar%auxiliary_data(iauxiliary+1) = area
  ! mol/sec
  rate = rate*material_auxvar%volume
  do i = 1, ncomp
    icomp = mineral%kinmnrlspecid(i,imnrl)
    Residual(icomp) = Residual(icomp) + mineral%kinmnrlstoich(i,imnrl)*rate
  enddo

  if (compute_derivative) then
    option%io_buffer = 'NUMERICAL_JACOBIAN must be specified within the &
      &NEWTON_SOLVER TRANSPORT block of NUMERICAL_METHODS TRANSPORT, as &
      &an analytical Jacobian is not provided for REACTION_SANDBOX,&
      &MINERAL_ARMORING.'
    call PrintErrMsg(option)
  endif

end subroutine MineralArmoringEvaluate

! ************************************************************************** !

subroutine MineralArmoringUpdateKineticState(this,rt_auxvar,global_auxvar, &
                                             material_auxvar,reaction,option)
  !
  ! Updates mineral volume fraction at end of converged timestep based on
  ! latest rate
  !
  ! Author: Aaron Sitchler
  ! Date: 08/24/26
  !
  use Option_module
  use Reaction_Aux_module
  use Reactive_Transport_Aux_module
  use Global_Aux_module
  use Material_Aux_module

  implicit none

  class(reaction_sandbox_mineral_armoring_type) :: this
  type(reactive_transport_auxvar_type) :: rt_auxvar
  type(global_auxvar_type) :: global_auxvar
  type(material_auxvar_type) :: material_auxvar
  class(reaction_rt_type) :: reaction
  type(option_type) :: option

  PetscInt :: imnrl
  PetscReal :: delta_volfrac

  imnrl = this%mineral_id
  ! rate = mol/m^3/sec
  ! dvolfrac = m^3 mnrl/m^3 bulk = rate (mol mnrl/m^3 bulk/sec) *
  !                                mol_vol (m^3 mnrl/mol mnrl)
  delta_volfrac = rt_auxvar%auxiliary_data(this%auxiliary_offset+1)* &
                  reaction%mineral%kinmnrl_molar_vol(imnrl)* &
                  option%tran_dt
  ! m^3 mnrl/m^3 bulk
  rt_auxvar%mnrl_volfrac(imnrl) = rt_auxvar%mnrl_volfrac(imnrl) + &
                                  delta_volfrac
  ! zero to avoid negative volume fractions
  if (rt_auxvar%mnrl_volfrac(imnrl) < 0.d0) &
    rt_auxvar%mnrl_volfrac(imnrl) = 0.d0

end subroutine MineralArmoringUpdateKineticState

end module Reaction_Sandbox_Mineral_Armoring_class
