PROGRAM main
  use Parameter_Variables
  use System_Variables
  use Dimensional_Variables
  use Static_Variables
  use Structure_module
  use Dynamics_Variables
  use Writting
  implicit  none
  call Hard_sphere_dynamics(.True.,.True.,.True.)

CONTAINS
  SUBROUTINE Hard_sphere_dynamics(VW_op,Gamm_op,Sk_writting_op)
    implicit  none
    ! Array of two dimension type Float with 8 presicion points, this is dynamic memory    
    real * 8, dimension(:,:), Allocatable :: gamma_values
    real * 8:: gam
    ! Vector(Array one Dimension, whatever who the F"#% cares)
    integer, Allocatable, Dimension(:) :: i_k_test ![X]
    integer :: i_k_test_number  ![X]
    logical, Intent (in) :: VW_op,Gamm_op,Sk_writting_op ![X]
    ! ################################################
    !Points assignement
    ! ################################################
    !Wave vector points //this is ommited cause this it is specified in the list of initialization
    kpoints=2**12               ![X]
    print *, 'kpoints value:', kpoints
    short_times=2**4
    medium_times=2**7
    decimations=2**5
    !this is ommited cause this it is specified in the list of initialization
    ! [ ] Testing this variable Value 
    i_k_test_number= 7 ![X]
    ! ################################################
    !Differential variables assignement
    ! ################################################
    !this is ommited cause this it is specified in the list of initialization 
    dk=1.0d-2                   ! [X]
    print *, 'value of dk', dk
    dtau=1.0d-7
    ! ################################################
    !Parameter System variables assignement and allocation
    ! ################################################
    !Space dimension
    SDimen=3  ! [X]
    ! Number of species
    Species=1 ![X]

    
    allocate(i_k_test(i_k_test_number)) ![X]
    i_k_test(1)=200  ![X]
    i_k_test(2)=350  ![X]
    i_k_test(3)=500  ![X]
    i_k_test(4)=710  ![X]
    i_k_test(5)=900  ![X]
    i_k_test(6)=1100 ![X]
    i_k_test(7)=1300 ![X]
    !Memory Allocation for System Variables
    call Sys_Variables_Alloc(species) ! [X]
    !System diameter
    sigma(1)=distu              ![X]
    !Diffusion coefficient of the system particle
    call D0M_ini(species,sigma) ![X]
    ! ################################################
    !Independant System Variables
    ! ################################################
    !Dimensionless Density
    eta(1)=0.3d0    ![X]//********
    !Calculation of m_rho and m_rhoi matrixes
    call rho_ini()  ![X]
    ! ################################################
    !Static Variables Memory Allocation
    !( ################################################
    call Static_Variables_mem_alloc(kpoints,species) ![X]
    call Calc_static_k(dk,kpoints)
    ! ################################################
    !Structure Factor Calculation
    ! ################################################
    call Calc_Sk_hs_py_mono(VW_op) !True for VW Correction, False Otherwise
    Call Sk_writting ("sk.dat") !Structure Factor Writting
  END SUBROUTINE Hard_sphere_dynamics
END PROGRAM main
