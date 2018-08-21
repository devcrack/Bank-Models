! Main Program !
!It will only call a subroutine to do something
Program Main
  Use Parameter_Variables
  Use System_Variables
  Use Dimensional_Variables
  Use Static_Variables
  Use Structure_module
  Use Writting
  Use SCGLE_Variables
  Use SCGLE
  Use NESCGLE_variables
  Use NESCGLE
  Implicit None
  ! ########################################################################
  ! ############ <The choose one> ##########################################
  ! ########################################################################
  !Call Hard_sphere_dynamics(.True.,.True.,.True.)
! ########################################################################

  !Call Hard_sphere_mixture_dynamics(.True.,.True.,.True.)
  !Call Hard_Disk_mixture_dynamics(.True.,.True.)
  !Call  Hard_Disk_mixture_arrest()
  !Call SW_santos_structure(.True.,.False.)
  !Call HD_AY_with_sharma_sharma_dynamics(.True.,.True.)
  !Call Hard_Disk_mono_dynamics(.True.,.True.)
  !Call Hard_Disk_iso_diffusivity_mixture(.True.,.True.)
  !Call Hard_Disk_mixture_multi_dynamics(.True.,.True.)
  !Call Calc_isochoric_ua_Tf_diagram_mono()
  Call Hard_Sphere_attractive_yukawa_mono_dynamics(.True.,.True.)
Contains
! Subroutines !
Subroutine Hard_sphere_dynamics(VW_op,Gamm_op,Sk_writting_op)
  Implicit None
  Real * 8, Dimension(:,:), Allocatable :: gamma_values
  Real * 8:: gam
  Integer, Allocatable, Dimension(:) :: i_k_test
  Integer :: i_k_test_number
  Logical, Intent (in) :: VW_op,Gamm_op,Sk_writting_op
  !Points assignement
  kpoints=2**12 !Wave vector points
  short_times=2**4
  medium_times=2**7
  decimations=2**5
  i_k_test_number=7
  !Differential variables assignement
  dk=1.0d-2
  dtau=1.0d-7
  !Parameter System variables assignement and allocation
  SDimen=3  !Space dimension
  Species=1 !Number of species
  Allocate(i_k_test(i_k_test_number))
  i_k_test(1)=200
  i_k_test(2)=350
  i_k_test(3)=500
  i_k_test(4)=710
  i_k_test(5)=900
  i_k_test(6)=1100
  i_k_test(7)=1300
  Call Sys_Variables_Alloc(species) !Memory Allocation for System Variables
  sigma(1)=distu  !System diameter
  Call D0M_ini(species,sigma) !Diffusion coefficient of the system particle
  !Independant System Variables
  eta(1)=0.3d0  !Dimensionless Density
  Call rho_ini() !Calculation of m_rho and m_rhoi matrixes
  !Static Variables Memory Allocation
  Call Static_Variables_mem_alloc(kpoints,species)
  Call Calc_static_k(dk,kpoints)
  !Structure Factor Calculation
  Call Calc_Sk_hs_py_mono(VW_op) !True for VW Correction, False Otherwise
  Call calc_ck_hk_ski()!Calculation for the direct correlation function and total correlation function
  If (Sk_writting_op.Eqv..True.) Then
    Call Sk_writting ("sk.dat") !Structure Factor Writting
  End If
  Call kc_ini(species,1.305d0)  !Initialization of kc, SCGLE parameter
  Call Calc_lambdak(Species,kpoints)  !Calculation of lambda(k), SCGLE parameter function
  Allocate(gamma_values(Species,Species)) !Allocation of gamma variables
  If (Gamm_op.Eqv..True.) Then !
    Call Gamma_mixture(gamma_values)  !Call of gamma subroutine
    gam=Gamma_mono()
  End If
  Call Long_time_dynamics(i_k_test,i_k_test_number,1)
  Call SCGLE_run_info_file("Run_info.dat",gamma_values,"Hard Sphere Monocomponent",i_k_test(1))
  Deallocate(i_k_test)
End Subroutine

Subroutine Hard_sphere_mixture_dynamics(VW_op,Gamm_op,Sk_writting_op)
  Implicit None
  Real * 8, Dimension(:,:), Allocatable :: gamma_values
  Real * 8:: gam
  Integer :: i_k_test_number
  Integer, Dimension(:), Allocatable :: i_k_test
  Logical, Intent (in) :: VW_op,Gamm_op,Sk_writting_op
  !Points assignement
  kpoints=2**12 !Wave vector points
  short_times=10
  medium_times=100
  decimations=2**5
  i_k_test_number=1
  Allocate(i_k_test(i_k_test_number))
  i_k_test(1)=710
  !Differential variables assignement
  dk=1.0d-2
  dtau=1.0d-7
  !Parameter System variables assignement and allocation
  SDimen=3  !Space dimension
  Species=2 !Number of species
  Call Sys_Variables_Alloc(species) !Memory Allocation for System Variables
  sigma(1)=distu  !System diameter
  sigma(2)=0.5d0*sigma(1)
  Call D0M_ini(species,sigma) !Diffusion coefficient of the system particle
  !Independant System Variables
  eta(1)=0.2d0  !Dimensionless Density
  eta(2)=0.2d0
  Print *, "Holis"
  Call rho_ini() !Calculation of m_rho and m_rhoi matrixes
  !Static Variables Memory Allocation
  Call Static_Variables_mem_alloc(kpoints,species)
  Call Calc_static_k(dk,kpoints)
  !Structure Factor Calculation
  Call Calc_sk_HS_py_mixture()  !Calculation of Hard Sphere Structure Factor of a mixture using Percus-Yevick Approximation
  Call calc_ck_hk_ski() !Calculation for the direct correlation function, total correlation function and Structure Factor Inverse
  If (Sk_writting_op.Eqv..True.) Then
    Call Sk_writting ("sk.dat") !Structure Factor Writting
    stop
  
  End If
  Call kc_ini(species,1.545d0)  !Initialization of kc, SCGLE parameter
  Call Calc_lambdak(Species,kpoints)  !Calculation of lambda(k), SCGLE parameter function
  Allocate(gamma_values(Species,Species)) !Allocation of gamma variables
  If (Gamm_op.Eqv..True.) Then
    Call Gamma_mixture(gamma_values)  !Call of gamma subroutine
  End If
  Call Long_time_dynamics(i_k_test,i_k_test_number,1)
  Call SCGLE_run_info_file("Run_info.dat",gamma_values,"Hard Sphere Mixture",i_k_test(1))
End Subroutine


Subroutine Hard_Disk_mixture_dynamics(Gamm_op,Sk_writting_op)
  Implicit None
  Real * 8, Dimension(:,:), Allocatable :: gamma_values
  Real * 8:: gam
  Integer :: i_k_test_number
  Integer, Dimension(:), Allocatable :: i_k_test
  Logical, Intent (in) :: Gamm_op,Sk_writting_op
  !Points assignement
  kpoints=2**13 !Wave vector points
  short_times=10!2**4
  medium_times=2**7
  decimations=2**5
  i_k_test_number=7
  Allocate(i_k_test(i_k_test_number))
  i_k_test(1)=200
  i_k_test(2)=350
  i_k_test(3)=500
  i_k_test(4)=710
  i_k_test(5)=900
  i_k_test(6)=1100
  i_k_test(7)=1300
  !Differential variables assignement
  dk=1.0d-2
  dtau=1.0d-7
  !Parameter System variables assignement and allocation
  SDimen=2  !Space dimension
  Species=2 !Number of species
  Call Sys_Variables_Alloc(species) !Memory Allocation for System Variables
  sigma(1)=distu  !System diameter
  sigma(2)=0.69d0*sigma(1)
  Call D0M_ini(species,sigma) !Diffusion coefficient of the system particle
  !Independant System Variables
  eta(1)=0.01d0  !Dimensionless Density
  eta(2)=0.665d0
  Print *, "Holis"
  Call rho_ini() !Calculation of m_rho and m_rhoi matrixes
  !Static Variables Memory Allocation
  Call Static_Variables_mem_alloc(kpoints,species)
  Call Calc_static_k(dk,kpoints)
  !Structure Factor Calculation
  Call calc_ck_hd_binary_mixture_roth() !Calculation of the Structure Factor of a Hard Disk mixture System using Roth Expressions
  Call calc_ck_hk_ski()!Calculation for the direct correlation function and total correlation function
  If (Sk_writting_op.Eqv..True.) Then
    Call Sk_writting ("sk.dat") !Structure Factor Writting
  End If
!  Call kc_ini(species,1.15d0)  !Initialization of kc, SCGLE parameter, 1.15 for arrest at 0.70
  Call kc_ini(species,0.7d0)  !Initialization of kc, SCGLE parameter, 0.7 for arrest at 0.72
  Call Calc_lambdak(Species,kpoints)  !Calculation of lambda(k), SCGLE parameter function
  Allocate(gamma_values(Species,Species)) !Allocation of gamma variables
  If (Gamm_op.Eqv..True.) Then !
    Call Gamma_mixture(gamma_values)  !Call of gamma subroutine
    !Stop
  End If
  Call Long_time_dynamics(i_k_test,i_k_test_number,1)
  Call SCGLE_run_info_file("Run_info.dat",gamma_values,"Hard Disk Mixture",i_k_test(1))
  Print*, "mean rho Dl=",(((m_rho(1,1)**2.d0)*Dl(1))+((m_rho(2,2)**2.d0)*Dl(2)))/((m_rho(1,1)**2.d0)+(m_rho(2,2)**2.d0))
  Print*, "mean phi Dl=",(((eta(1))*Dl(1))+((eta(2))*Dl(2)))/(sum(eta))
  DEAllocate(i_k_test)
End Subroutine

Subroutine Hard_Disk_mono_dynamics(Gamm_op,Sk_writting_op)
  Implicit None
  Real * 8, Dimension(:,:), Allocatable :: gamma_values
  Real * 8:: gam
  Integer :: i_k_test_number
  Integer, Dimension(:), Allocatable :: i_k_test
  Logical, Intent (in) :: Gamm_op,Sk_writting_op
  !Points assignement
  kpoints=2**13 !Wave vector points
  short_times=10!2**4
  medium_times=2**7
  decimations=2**5
  i_k_test_number=7
  Allocate(i_k_test(i_k_test_number))
  i_k_test(1)=200
  i_k_test(2)=350
  i_k_test(3)=500
  i_k_test(4)=710
  i_k_test(5)=900
  i_k_test(6)=1100
  i_k_test(7)=1300
  !Differential variables assignement
  dk=1.0d-2
  dtau=1.0d-7
  !Parameter System variables assignement and allocation
  SDimen=2  !Space dimension
  Species=1 !Number of species
  Call Sys_Variables_Alloc(species) !Memory Allocation for System Variables
  sigma(1)=distu  !System diameter
  !sigma(2)=0.69d0*sigma(1)
  Call D0M_ini(species,sigma) !Diffusion coefficient of the system particle
  !Independant System Variables
  eta(1)=0.67d0  !Dimensionless Density
  !eta(2)=0.53823d0
  Print *, "Holis"
  Call rho_ini() !Calculation of m_rho and m_rhoi matrixes
  !Static Variables Memory Allocation
  Call Static_Variables_mem_alloc(kpoints,species)
  Call Calc_static_k(dk,kpoints)
  !Structure Factor Calculation
  Call calc_ck_hd_mono_roth() !Calculation of the Structure Factor of a Hard Disk monocomponent System using Roth Expressions
  Call calc_ck_hk_ski()!Calculation for the direct correlation function and total correlation function
  If (Sk_writting_op.Eqv..True.) Then
    Call Sk_writting ("sk.dat") !Structure Factor Writting
  End If
  !Call kc_ini(species,1.15d0)  !Initialization of kc, SCGLE parameter, 1.15 for arrest at 0.70
  Call kc_ini(species,0.7d0)  !Initialization of kc, SCGLE parameter, 0.7 for arrest at 0.72
  Call Calc_lambdak(Species,kpoints)  !Calculation of lambda(k), SCGLE parameter function
  Allocate(gamma_values(Species,Species)) !Allocation of gamma variables
  If (Gamm_op.Eqv..True.) Then !
    Call Gamma_mixture(gamma_values)  !Call of gamma subroutine
    !Stop
  End If
  Call Long_time_dynamics(i_k_test,i_k_test_number,1)
  Call SCGLE_run_info_file("Run_info.dat",gamma_values,"Hard Disk Mixture",i_k_test(1))
  DEAllocate(i_k_test)
End Subroutine


Subroutine SW_santos_structure(Sk_writting_op,Gamm_op)
  Implicit None
  Integer :: i_k_test_number
  Integer, Dimension(:), Allocatable:: i_k_test
  Logical, Intent (in) :: Gamm_op,Sk_writting_op
  Real * 8 :: Temp,compressibility
  Real * 8, Dimension(:,:), Allocatable :: gamma_values
  !Points assignement
  kpoints=2**12 !Wave vector points
  short_times=10
  medium_times=100
  decimations=2**5
  i_k_test_number=1
  Allocate(i_k_test(i_k_test_number))
  i_k_test(1)=710
  !Differential variables assignement
  dk=1.0d-2
  dtau=1.0d-7
  !Parameter System variables assignement and allocation
  SDimen=3  !Space dimension
  Species=1 !Number of species
  Call Sys_Variables_Alloc(species) !Memory Allocation for System Variables
  sigma(1)=distu  !System diameter
  Call D0M_ini(species,sigma) !Diffusion coefficient of the system particle
  !Independant System Variables
  eta(1)=0.2d0  !Dimensionless Density
  Temp=0.428447705d0 !Dimensionless Temperature
  Print *, "Holis"
  Call rho_ini() !Calculation of m_rho and m_rhoi matrixes
  !Static Variables Memory Allocation
  Call Static_Variables_mem_alloc(kpoints,species)
  Call Calc_static_k(dk,kpoints)
  !Structure Factor Calculation
  Call calc_sk_hs_sw_mono_santos (Temp) !Calculation of the Structure Factor of a Hard Disk Square Well System using Santos Expressions
  Call calc_ck_hk_ski()!Calculation for the direct correlation function and total correlation function
  If (Sk_writting_op.Eqv..True.) Then
    Call Sk_writting ("sk.dat") !Structure Factor Writting
    Call Calc_S0_HS_SW_santos(Temp,compressibility)
    Print*, compressibility, "Holis"
  End If
  Stop
  Call kc_ini(species,1.15d0)  !Initialization of kc, SCGLE parameter
  Call Calc_lambdak(Species,kpoints)  !Calculation of lambda(k), SCGLE parameter function
  Allocate(gamma_values(Species,Species)) !Allocation of gamma variables
  If (Gamm_op.Eqv..True.) Then !
    Call Gamma_mixture(gamma_values)  !Call of gamma subroutin
    !gam=Gamma_mono()
  End If
  Call Long_time_dynamics(i_k_test,i_k_test_number,1)
  Call SCGLE_run_info_file("Run_info.dat",gamma_values,"Hard Disk Mixture",i_k_test(1))
End Subroutine

Subroutine HD_AY_with_sharma_sharma_dynamics(Sk_writting_op,Gamm_op)
  Implicit None
  Integer :: i_k_test_number
  Integer, Dimension(:), Allocatable:: i_k_test
  Logical, Intent (in) :: Gamm_op,Sk_writting_op
  Real * 8 :: Temp,compressibility
  Real * 8, Dimension(:,:), Allocatable :: gamma_values
  !Points assignement
  kpoints=2**12 !Wave vector points
  short_times=10
  medium_times=100
  decimations=2**5
  i_k_test_number=1
  Allocate(i_k_test(i_k_test_number))
  i_k_test(1)=710
  !Differential variables assignement
  dk=1.0d-2
  dtau=1.0d-7
  !Parameter System variables assignement and allocation
  SDimen=2    !Space dimension
  Species=1 !Number of species
  Call Sys_Variables_Alloc(species) !Memory Allocation for System Variables
  sigma(1)=distu  !System diameter
  Call D0M_ini(species,sigma) !Diffusion coefficient of the system particle
  !Independant System Variables
  eta(1)=0.2d0  !Dimensionless Density
  Temp=0.5!Dimensionless Temperature
  Print *, "Holis"
  Call rho_ini() !Calculation of m_rho and m_rhoi matrixes
  !Static Variables Memory Allocation
  Call Static_Variables_mem_alloc(kpoints,species)
  Call Calc_static_k(dk,kpoints)
  Call  z_yuk_mono(2.d0,1)
  !Structure Factor Calculation
  Call sharma_sharma_attractive_yukawa(Temp) !Calculation of the Structure Factor of a Hard Disk Square Well System using Santos Expressions
  Call calc_ck_hk_ski()!Calculation for the direct correlation function and total correlation function
  If (Sk_writting_op.Eqv..True.) Then
    Call Sk_writting ("sk.dat") !Structure Factor Writting
  End If
  Call kc_ini(species,1.305d0)  !Init ialization of kc, SCGLE parameter
  Call Calc_lambdak(Species,kpoints)  !Calculation of lambda(k), SCGLE parameter function
  Allocate(gamma_values(Species,Species)) !Allocation of gamma variables
  If (Gamm_op.Eqv..True.) Then
    Call Gamma_mixture(gamma_values)  !Call of gamma subroutine
    !gam=Gamma_mono()
  End If
  stop
  Call Long_time_dynamics(i_k_test,i_k_test_number,1)
  Call SCGLE_run_info_file("Run_info.dat",gamma_values,"Hard Disk Mixture",i_k_test(1))
End Subroutine
!
Subroutine Hard_Disk_mixture_arrest()
  Implicit None
  Real * 8, Dimension(:,:), Allocatable :: gamma_values
  Integer :: i1
  Logical :: condition
  Real * 8 :: deta,eta_max,eta_min,error,gam1,gam2
  !Points assignement
  kpoints=2**12 !Wave vector points
  !Differential variables assignement
  dk=1.0d-2
  deta=1.0d-2
  !Parameter System variables assignement and allocation
  SDimen=2  !Space dimension
  Species=2 !Number of species
  Call Sys_Variables_Alloc(species) !Memory Allocation for System Variables
  sigma(1)=distu  !System diameter
  sigma(2)=0.69d0*sigma(1)
  Call D0M_ini(species,sigma) !Diffusion coefficient of the system particle
  Allocate(gamma_values(Species,Species)) !Allocation of gamma variables
  Call Static_Variables_mem_alloc(kpoints,species)
  Call Calc_static_k(dk,kpoints)
  Print*, "holis"

  Print*, "holis"
!  Call kc_ini(species,1.15d0)  !Initialization of kc, SCGLE parameter, 1.15 for arrest at 0.70    Call kc_ini(species,0.7d0)  !Initialization of kc, SCGLE parameter, 0.7 for arrest at 0.72
  Call kc_ini(species,0.7d0)  !Initialization of kc, SCGLE parameter, 0.7 for arrest at 0.72
  Call Calc_lambdak(Species,kpoints)  !Calculation of lambda(k), SCGLE parameter function
  Open (unit=50,File="arrest.dat",Status="Replace")
  Do i1=1,70
    !Independant System Variables
    eta(1)=i1*0.01d0  !Dimensionless Density
    eta(2)=0.71d0-(0.01*i1)
    Call rho_ini() !Calculation of m_rho and m_rhoi matrixes
    !Structure Factor Calculation
    Call calc_ck_hd_binary_mixture_roth() !Calculation of the Structure Factor of a Hard Disk mixture System using Roth Expressions
    Call calc_ck_hk_ski()!Calculation for the direct correlation function and total correlation function
    Call Gamma_mixture(gamma_values)  !Call of gamma subroutine
    condition=.False.
    If (gamma_values(1,1)>gamma_max.And.gamma_values(2,2)>gamma_max) Then
      Do While (condition.Eqv..False.)
        eta_min=eta(2)
        eta(2)=eta(2)+deta
        !Print *, "Holis"
        Call rho_ini() !Calculation of m_rho and m_rhoi matrixes
        !Structure Factor Calculation
        Call calc_ck_hd_binary_mixture_roth() !Calculation of the Structure Factor of a Hard Disk mixture System using Roth Expressions
        Call calc_ck_hk_ski()!Calculation for the direct correlation function and total correlation function
        Call Gamma_mixture(gamma_values)  !Call of gamma subroutine
        If (gamma_values(1,1)<gamma_max.And.gamma_values(2,2)<gamma_max) Then
          condition=.True.
          eta_max=eta(2)
        End If
      End Do
    Else
      Do While (condition.Eqv..False.)
        eta_max=eta(2)
        eta(2)=eta(2)-deta
        !Print *, "Holis"
        Call rho_ini() !Calculation of m_rho and m_rhoi matrixes
        !Structure Factor Calculation
        Call calc_ck_hd_binary_mixture_roth() !Calculation of the Structure Factor of a Hard Disk mixture System using Roth Expressions
        Call calc_ck_hk_ski()!Calculation for the direct correlation function and total correlation function
        Call Gamma_mixture(gamma_values)  !Call of gamma subroutine
        If (gamma_values(1,1)>gamma_max.And.gamma_values(2,2)>gamma_max) Then
          condition=.True.
          eta_min=eta(2)
        End If
      End Do
    End If
    error=1.d0
    Do While (error>1.d-6)
      eta(2)=(eta_max+eta_min)/2.d0
      !Print *, "Holis"
      Call rho_ini() !Calculation of m_rho and m_rhoi matrixes
      !Structure Factor Calculation
      Call calc_ck_hd_binary_mixture_roth() !Calculation of the Structure Factor of a Hard Disk mixture System using Roth Expressions
      Call calc_ck_hk_ski()!Calculation for the direct correlation function and total correlation function
      Call Gamma_mixture(gamma_values)  !Call of gamma subroutine
      If (gamma_values(1,1)>gamma_max.And.gamma_values(2,2)>gamma_max) Then
        eta_min=eta(2)
      Else If(gamma_values(1,1)<gamma_max.And.gamma_values(2,2)<gamma_max) Then
        eta_max=eta(2)
        gam1=gamma_values(1,1)
        gam2=gamma_values(2,2)
      End If
      error=Abs(eta_max-eta_min)
      print*, "error=", error
    End Do
    Write(50,*) eta(1), eta_max, gam1,gam2
    Flush(50)
  End Do
  Close(50)
End Subroutine

Subroutine Hard_Disk_iso_diffusivity_mixture(Gamm_op,Sk_writting_op)
  Implicit None
  Real * 8, Dimension(:,:), Allocatable :: gamma_values
  Real * 8:: gam,mean_DL,Wanted_Dl,Dl_error,eta_min,eta_max
  Integer :: i_k_test_number,i1
  Integer, Dimension(:), Allocatable :: i_k_test
  Logical, Intent (in) :: Gamm_op,Sk_writting_op
  Character (LEN=2) :: char_run
  Real * 8, Parameter :: Dl_tolerance=1.d-3
  !Points assignement
  kpoints=2**12 !Wave vector points
  short_times=10!2**4
  medium_times=2**7
  decimations=2**5
  i_k_test_number=7
  Allocate(i_k_test(i_k_test_number))
  i_k_test(1)=200
  i_k_test(2)=350
  i_k_test(3)=500
  i_k_test(4)=710
  i_k_test(5)=900
  i_k_test(6)=1100
  i_k_test(7)=1300
  !Differential variables assignement
  dk=1.0d-2
  dtau=1.0d-7
  !Parameter System variables assignement and allocation
  SDimen=2  !Space dimension
  Species=2 !Number of species
  Call Sys_Variables_Alloc(species) !Memory Allocation for System Variables
  sigma(1)=distu  !System diameter
  sigma(2)=0.69d0*sigma(1)
  Call D0M_ini(species,sigma) !Diffusion coefficient of the system particle
  !Static Variables Memory Allocation
  Call Static_Variables_mem_alloc(kpoints,species)
  Call Calc_static_k(dk,kpoints)
  !  Call kc_ini(species,1.15d0)  !Initialization of kc, SCGLE parameter, 1.15 for arrest at 0.70
  Call kc_ini(species,0.7d0)  !Initialization of kc, SCGLE parameter, 0.7 for arrest at 0.72
  Call Calc_lambdak(Species,kpoints)  !Calculation of lambda(k), SCGLE parameter function
  Allocate(gamma_values(Species,Species)) !Allocation of gamma variables
  Wanted_Dl=1.d-1
  !Independant System Variables
  Open (unit=50,File="isodifussive.dat",Status="Replace")
  Do i1=1,34
    write(char_run,"(I2.1)") i1
    dtau=1.0d-7
    eta(1)=0.67d0-DBLE(i1-1)*2.d-2  !Dimensionless Density
    eta(2)=0.005d0+(DBLE(i1-1)*2.d-2)
    Call rho_ini() !Calculation of m_rho and m_rhoi matrixes
    !Structure Factor Calculation
    Call calc_ck_hd_binary_mixture_roth() !Calculation of the Structure Factor of a Hard Disk mixture System using Roth Expressions
    Call calc_ck_hk_ski()!Calculation for the direct correlation function and total correlation function
    If (Sk_writting_op.Eqv..True.) Then
      Call Sk_writting("sk_"//char_run//".dat")  !Structure Factor Writting
    End If
    If (Gamm_op.Eqv..True.) Then !
      Call Gamma_mixture(gamma_values)  !Call of gamma subroutine
    End If
    Call Long_time_dynamics(i_k_test,i_k_test_number,i1)
    mean_DL=Dl(1)*eta(1)+Dl(2)*eta(2)
    mean_DL=mean_DL/sum(eta)
    Dl_error=Abs(mean_Dl-Wanted_Dl)/Wanted_Dl
    Call Dynamics_Variables_mem_deallocation()
    dtau=1.0d-7
    If (mean_DL>Wanted_Dl) Then
      eta_min=eta(2)
      Do while (mean_DL>Wanted_Dl)
        eta(2)=eta(2)+0.005d0
        Call rho_ini() !Calculation of m_rho and m_rhoi matrixes
        !Structure Factor Calculation
        Call calc_ck_hd_binary_mixture_roth() !Calculation of the Structure Factor of a Hard Disk mixture System using Roth Expressions
        Call calc_ck_hk_ski()!Calculation for the direct correlation function and total correlation function
        If (Sk_writting_op.Eqv..True.) Then
          Call Sk_writting ("sk_"//char_run//".dat") !Structure Factor Writting
        End If
        If (Gamm_op.Eqv..True.) Then !
          Call Gamma_mixture(gamma_values)  !Call of gamma subroutine
        End If
        Call Long_time_dynamics(i_k_test,i_k_test_number,i1)
        mean_DL=Dl(1)*eta(1)+Dl(2)*eta(2)
        mean_DL=mean_DL/sum(eta)
        Dl_error=Abs(mean_Dl-Wanted_Dl)/Wanted_Dl
        If(mean_DL>Wanted_Dl) Then
          eta_min=eta(2)
        End If
        Call Dynamics_Variables_mem_deallocation()
        dtau=1.0d-7
      End Do
      eta_max=eta(2)
    Else
      eta_max=eta(2)
      Do while (mean_DL<Wanted_Dl)
        eta(2)=eta(2)-0.005d0
        Call rho_ini() !Calculation of m_rho and m_rhoi matrixes
        !Structure Factor Calculation
        Call calc_ck_hd_binary_mixture_roth() !Calculation of the Structure Factor of a Hard Disk mixture System using Roth Expressions
        Call calc_ck_hk_ski()!Calculation for the direct correlation function and total correlation function
        If (Sk_writting_op.Eqv..True.) Then
          Call Sk_writting ("sk_"//char_run//".dat") !Structure Factor Writting
        End If
        If (Gamm_op.Eqv..True.) Then !
          Call Gamma_mixture(gamma_values)  !Call of gamma subroutine
        End If
        Call Long_time_dynamics(i_k_test,i_k_test_number,i1)
        mean_DL=Dl(1)*eta(1)+Dl(2)*eta(2)
        mean_DL=mean_DL/sum(eta)
        Dl_error=Abs(mean_Dl-Wanted_Dl)/Wanted_Dl
        If(mean_DL<Wanted_Dl) Then
          eta_max=eta(2)
        End If
        Call Dynamics_Variables_mem_deallocation()
        dtau=1.0d-7
      End Do
      eta_min=eta(2)
    End If
    Do While(Dl_error>Dl_tolerance)
      eta(2)=(eta_min+eta_max)/2.d0
      Call rho_ini() !Calculation of m_rho and m_rhoi matrixes
      !Structure Factor Calculation
      Call calc_ck_hd_binary_mixture_roth() !Calculation of the Structure Factor of a Hard Disk mixture System using Roth Expressions
      Call calc_ck_hk_ski()!Calculation for the direct correlation function and total correlation function
      If (Sk_writting_op.Eqv..True.) Then
        Call Sk_writting ("sk_"//char_run//".dat") !Structure Factor Writting
      End If
      If (Gamm_op.Eqv..True.) Then !
        Call Gamma_mixture(gamma_values)  !Call of gamma subroutine
      End If
      Call Long_time_dynamics(i_k_test,i_k_test_number,i1)
      mean_DL=Dl(1)*eta(1)+Dl(2)*eta(2)
      mean_DL=mean_DL/sum(eta)
      Dl_error=Abs(mean_Dl-Wanted_Dl)/Wanted_Dl
      If(mean_DL<Wanted_Dl) Then
        eta_max=eta(2)
      Else
        eta_min=eta(2)
      End If
      If (Dl_error>Dl_tolerance) Then
        Call Dynamics_Variables_mem_deallocation()
        dtau=1.0d-7
      End If
    End Do
    Call SCGLE_run_info_file("Run_info_"//char_run//".dat",gamma_values,"Hard Disk Mixture",i_k_test(1))
    Write(50,*) eta(1),eta(2),Dl(1),Dl(2),mean_DL
    flush(50)
    Call Dynamics_Variables_mem_deallocation()
    dtau=1.0d-7
  End Do
  Close(50)
  DEAllocate(i_k_test)
End Subroutine

Subroutine Hard_Disk_iso_diffusivity_mixture_2(Gamm_op,Sk_writting_op)
  Implicit None
  Real * 8, Dimension(:,:), Allocatable :: gamma_values
  Real * 8:: gam,mean_DL,Wanted_Dl,Dl_error,eta_min,eta_max
  Integer :: i_k_test_number,i1
  Integer, Dimension(:), Allocatable :: i_k_test
  Logical, Intent (in) :: Gamm_op,Sk_writting_op
  Character (LEN=2) :: char_run
  Real * 8, Parameter :: Dl_tolerance=1.d-3
  !Points assignement
  kpoints=2**12 !Wave vector points
  short_times=10!2**4
  medium_times=2**7
  decimations=2**5
  i_k_test_number=7
  Allocate(i_k_test(i_k_test_number))
  i_k_test(1)=200
  i_k_test(2)=350
  i_k_test(3)=500
  i_k_test(4)=710
  i_k_test(5)=900
  i_k_test(6)=1100
  i_k_test(7)=1300
  !Differential variables assignement
  dk=1.0d-2
  dtau=1.0d-7
  !Parameter System variables assignement and allocation
  SDimen=2  !Space dimension
  Species=2 !Number of species
  Call Sys_Variables_Alloc(species) !Memory Allocation for System Variables
  sigma(1)=distu  !System diameter
  sigma(2)=0.69d0*sigma(1)
  Call D0M_ini(species,sigma) !Diffusion coefficient of the system particle
  !Static Variables Memory Allocation
  Call Static_Variables_mem_alloc(kpoints,species)
  Call Calc_static_k(dk,kpoints)
  !  Call kc_ini(species,1.15d0)  !Initialization of kc, SCGLE parameter, 1.15 for arrest at 0.70
  Call kc_ini(species,0.7d0)  !Initialization of kc, SCGLE parameter, 0.7 for arrest at 0.72
  Call Calc_lambdak(Species,kpoints)  !Calculation of lambda(k), SCGLE parameter function
  Allocate(gamma_values(Species,Species)) !Allocation of gamma variables
  Wanted_Dl=1.d-1
  !Independant System Variables
  Open (unit=50,File="isodifussive.dat",Status="Replace")
  Do i1=1,34
    write(char_run,"(I2.1)") i1
    dtau=1.0d-7
    eta(1)=0.67d0-DBLE(i1-1)*2.d-2  !Dimensionless Density
    eta(2)=0.005d0+(DBLE(i1-1)*2.d-2)
    Call rho_ini() !Calculation of m_rho and m_rhoi matrixes
    !Structure Factor Calculation
    Call calc_ck_hd_binary_mixture_roth() !Calculation of the Structure Factor of a Hard Disk mixture System using Roth Expressions
    Call calc_ck_hk_ski()!Calculation for the direct correlation function and total correlation function
    If (Sk_writting_op.Eqv..True.) Then
      Call Sk_writting("sk_"//char_run//".dat")  !Structure Factor Writting
    End If
    If (Gamm_op.Eqv..True.) Then !
      Call Gamma_mixture(gamma_values)  !Call of gamma subroutine
    End If
    Call Long_time_dynamics(i_k_test,i_k_test_number,i1)
    mean_DL=Dl(1)*(m_rho(1,1)**2.d0)+Dl(2)*(m_rho(2,2)**2.d0)
    mean_DL=mean_DL/((m_rho(1,1)**2.d0)+(m_rho(2,2)**2.d0))
    Dl_error=Abs(mean_Dl-Wanted_Dl)/Wanted_Dl
    Call Dynamics_Variables_mem_deallocation()
    dtau=1.0d-7
    If (mean_DL>Wanted_Dl) Then
      eta_min=eta(2)
      Do while (mean_DL>Wanted_Dl)
        eta(2)=eta(2)+0.005d0
        Call rho_ini() !Calculation of m_rho and m_rhoi matrixes
        !Structure Factor Calculation
        Call calc_ck_hd_binary_mixture_roth() !Calculation of the Structure Factor of a Hard Disk mixture System using Roth Expressions
        Call calc_ck_hk_ski()!Calculation for the direct correlation function and total correlation function
        If (Sk_writting_op.Eqv..True.) Then
          Call Sk_writting ("sk_"//char_run//".dat") !Structure Factor Writting
        End If
        If (Gamm_op.Eqv..True.) Then !
          Call Gamma_mixture(gamma_values)  !Call of gamma subroutine
        End If
        Call Long_time_dynamics(i_k_test,i_k_test_number,i1)
        mean_DL=Dl(1)*(m_rho(1,1)**2.d0)+Dl(2)*(m_rho(2,2)**2.d0)
        mean_DL=mean_DL/((m_rho(1,1)**2.d0)+(m_rho(2,2)**2.d0))
        Dl_error=Abs(mean_Dl-Wanted_Dl)/Wanted_Dl
        If(mean_DL>Wanted_Dl) Then
          eta_min=eta(2)
        End If
        Call Dynamics_Variables_mem_deallocation()
        dtau=1.0d-7
      End Do
      eta_max=eta(2)
    Else
      eta_max=eta(2)
      Do while (mean_DL<Wanted_Dl)
        eta(2)=eta(2)-0.005d0
        Call rho_ini() !Calculation of m_rho and m_rhoi matrixes
        !Structure Factor Calculation
        Call calc_ck_hd_binary_mixture_roth() !Calculation of the Structure Factor of a Hard Disk mixture System using Roth Expressions
        Call calc_ck_hk_ski()!Calculation for the direct correlation function and total correlation function
        If (Sk_writting_op.Eqv..True.) Then
          Call Sk_writting ("sk_"//char_run//".dat") !Structure Factor Writting
        End If
        If (Gamm_op.Eqv..True.) Then !
          Call Gamma_mixture(gamma_values)  !Call of gamma subroutine
        End If
        Call Long_time_dynamics(i_k_test,i_k_test_number,i1)
        mean_DL=Dl(1)*(m_rho(1,1)**2.d0)+Dl(2)*(m_rho(2,2)**2.d0)
        mean_DL=mean_DL/((m_rho(1,1)**2.d0)+(m_rho(2,2)**2.d0))
        Dl_error=Abs(mean_Dl-Wanted_Dl)/Wanted_Dl
        If(mean_DL<Wanted_Dl) Then
          eta_max=eta(2)
        End If
        Call Dynamics_Variables_mem_deallocation()
        dtau=1.0d-7
      End Do
      eta_min=eta(2)
    End If
    Do While(Dl_error>Dl_tolerance)
      eta(2)=(eta_min+eta_max)/2.d0
      Call rho_ini() !Calculation of m_rho and m_rhoi matrixes
      !Structure Factor Calculation
      Call calc_ck_hd_binary_mixture_roth() !Calculation of the Structure Factor of a Hard Disk mixture System using Roth Expressions
      Call calc_ck_hk_ski()!Calculation for the direct correlation function and total correlation function
      If (Sk_writting_op.Eqv..True.) Then
        Call Sk_writting ("sk_"//char_run//".dat") !Structure Factor Writting
      End If
      If (Gamm_op.Eqv..True.) Then !
        Call Gamma_mixture(gamma_values)  !Call of gamma subroutine
      End If
      Call Long_time_dynamics(i_k_test,i_k_test_number,i1)
      mean_DL=Dl(1)*(m_rho(1,1)**2.d0)+Dl(2)*(m_rho(2,2)**2.d0)
      mean_DL=mean_DL/((m_rho(1,1)**2.d0)+(m_rho(2,2)**2.d0))
      Dl_error=Abs(mean_Dl-Wanted_Dl)/Wanted_Dl
      If(mean_DL<Wanted_Dl) Then
        eta_max=eta(2)
      Else
        eta_min=eta(2)
      End If
      If (Dl_error>Dl_tolerance) Then
        Call Dynamics_Variables_mem_deallocation()
        dtau=1.0d-7
      End If
    End Do
    Call SCGLE_run_info_file("Run_info_"//char_run//".dat",gamma_values,"Hard Disk Mixture",i_k_test(1))
    Write(50,*) eta(1),eta(2),Dl(1),Dl(2),mean_DL
    flush(50)
    Call Dynamics_Variables_mem_deallocation()
    dtau=1.0d-7
  End Do
  Close(50)
  DEAllocate(i_k_test)
End Subroutine

Subroutine Hard_Disk_mixture_multi_dynamics(Gamm_op,Sk_writting_op)
  Implicit None
  Real * 8, Dimension(:,:), Allocatable :: gamma_values
  Real * 8:: gam,mean_DL,eta_total
  Integer :: i_k_test_number,i1
  Integer, Dimension(:), Allocatable :: i_k_test
  Logical, Intent (in) :: Gamm_op,Sk_writting_op
  Character (LEN=2) :: char_run
  Real * 8, Parameter :: Dl_tolerance=1.d-3
  !Points assignement
  kpoints=2**13 !Wave vector points
  short_times=10!2**4
  medium_times=2**7
  decimations=2**5
  i_k_test_number=7
  Allocate(i_k_test(i_k_test_number))
  i_k_test(1)=200
  i_k_test(2)=350
  i_k_test(3)=500
  i_k_test(4)=710
  i_k_test(5)=900
  i_k_test(6)=1100
  i_k_test(7)=1300
  !Differential variables assignement
  dk=1.0d-2
  dtau=1.0d-7
  !Parameter System variables assignement and allocation
  SDimen=2  !Space dimension
  Species=2 !Number of species
  Call Sys_Variables_Alloc(species) !Memory Allocation for System Variables
  sigma(1)=distu  !System diameter
  sigma(2)=0.37d0*sigma(1)
  Call D0M_ini(species,sigma) !Diffusion coefficient of the system particle
  !Static Variables Memory Allocation
  Call Static_Variables_mem_alloc(kpoints,species)
  Call Calc_static_k(dk,kpoints)
  !  Call kc_ini(species,1.15d0)  !Initialization of kc, SCGLE parameter, 1.15 for arrest at 0.70
  Call kc_ini(species,0.7d0)  !Initialization of kc, SCGLE parameter, 0.7 for arrest at 0.72
  Call Calc_lambdak(Species,kpoints)  !Calculation of lambda(k), SCGLE parameter function
  Allocate(gamma_values(Species,Species)) !Allocation of gamma variables
  !Independant System Variables
  Open (unit=50,File="isodifussive.dat",Status="Replace")
  Do i1=1,12
    write(char_run,"(I2.1)") i1
    eta_total=0.05d0*DBLE(1+i1)
    dtau=1.0d-7
    eta(1)=(1.d0-0.5d0)*eta_total
    eta(2)=0.5d0*eta_total
    Call rho_ini() !Calculation of m_rho and m_rhoi matrixes
    !Structure Factor Calculation
    Call calc_ck_hd_binary_mixture_roth() !Calculation of the Structure Factor of a Hard Disk mixture System using Roth Expressions
    Call calc_ck_hk_ski()!Calculation for the direct correlation function and total correlation function
    If (Sk_writting_op.Eqv..True.) Then
      Call Sk_writting("sk_"//char_run//".dat")  !Structure Factor Writting
    End If
    If (Gamm_op.Eqv..True.) Then !
      Call Gamma_mixture(gamma_values)  !Call of gamma subroutine
    End If
    Call Long_time_dynamics(i_k_test,i_k_test_number,i1)
    Call SCGLE_run_info_file("Run_info_"//char_run//".dat",gamma_values,"Hard Disk Mixture",i_k_test(1))
    mean_DL=Dl(1)*(m_rho(1,1)**2.d0)+Dl(2)*(m_rho(2,2)**2.d0)
    mean_DL=mean_DL/((m_rho(1,1)**2.d0)+(m_rho(2,2)**2.d0))
    Write(50,*) eta(1),eta(2),Dl(1),Dl(2),mean_DL
    flush(50)
    Call Dynamics_Variables_mem_deallocation()
  End Do
  Close(50)
  DEAllocate(i_k_test)
End Subroutine

Subroutine Hard_Disk_iso_diffusivity_mixture_3(Gamm_op,Sk_writting_op)
  Implicit None
  Real * 8, Dimension(:,:), Allocatable :: gamma_values
  Real * 8:: gam,mean_DL,Wanted_Dl,Dl_error,eta_min,eta_max
  Integer :: i_k_test_number,i1
  Integer, Dimension(:), Allocatable :: i_k_test
  Logical, Intent (in) :: Gamm_op,Sk_writting_op
  Character (LEN=2) :: char_run
  Real * 8, Parameter :: Dl_tolerance=1.d-3
  !Points assignement
  Print*, "Hard_Disk_iso_diffusivity_mixture_3"
  kpoints=2**12 !Wave vector points
  short_times=10!2**4
  medium_times=2**7
  decimations=2**5
  i_k_test_number=7
  Allocate(i_k_test(i_k_test_number))
  i_k_test(1)=200
  i_k_test(2)=350
  i_k_test(3)=500
  i_k_test(4)=710
  i_k_test(5)=900
  i_k_test(6)=1100
  i_k_test(7)=1300
  !Differential variables assignement
  dk=1.0d-2
  dtau=1.0d-7
  !Parameter System variables assignement and allocation
  SDimen=2  !Space dimension
  Species=2 !Number of species
  Call Sys_Variables_Alloc(species) !Memory Allocation for System Variables
  sigma(1)=distu  !System diameter
  sigma(2)=0.69d0*sigma(1)
  Call D0M_ini(species,sigma) !Diffusion coefficient of the system particle
  !Static Variables Memory Allocation
  Call Static_Variables_mem_alloc(kpoints,species)
  Call Calc_static_k(dk,kpoints)
  !  Call kc_ini(species,1.15d0)  !Initialization of kc, SCGLE parameter, 1.15 for arrest at 0.70
  Call kc_ini(species,0.7d0)  !Initialization of kc, SCGLE parameter, 0.7 for arrest at 0.72
  Call Calc_lambdak(Species,kpoints)  !Calculation of lambda(k), SCGLE parameter function
  Allocate(gamma_values(Species,Species)) !Allocation of gamma variables
  Wanted_Dl=1.d-1
  !Independant System Variables
  Open (unit=50,File="isodifussive.dat",Status="Replace")
  Do i1=10,34
    write(char_run,"(I2.1)") i1
    dtau=1.0d-7
    eta(1)=0.67d0-DBLE(i1-1)*2.d-2  !Dimensionless Density
    eta(2)=0.005d0+(DBLE(i1-1)*2.d-2)
    Call rho_ini() !Calculation of m_rho and m_rhoi matrixes
    !Structure Factor Calculation
    Call calc_ck_hd_binary_mixture_roth() !Calculation of the Structure Factor of a Hard Disk mixture System using Roth Expressions
    Call calc_ck_hk_ski()!Calculation for the direct correlation function and total correlation function
    If (Sk_writting_op.Eqv..True.) Then
      Call Sk_writting("sk_"//char_run//".dat")  !Structure Factor Writting
    End If
    If (Gamm_op.Eqv..True.) Then !
      Call Gamma_mixture(gamma_values)  !Call of gamma subroutine
    End If
    Call Long_time_dynamics(i_k_test,i_k_test_number,i1)
    mean_DL=Dl(1)*(m_rho(1,1)**2.d0)+Dl(2)*(m_rho(2,2)**2.d0)
    mean_DL=mean_DL/((m_rho(1,1)**2.d0)+(m_rho(2,2)**2.d0))
    mean_DL=mean_Dl*(eta(1)+(eta(2)*(1.d0/(sigma(2)**3.d0))))/(eta(1)+(eta(2)*(1.d0/(sigma(2)**2.d0))))
    Dl_error=Abs(mean_Dl-Wanted_Dl)/Wanted_Dl
    Call Dynamics_Variables_mem_deallocation()
    dtau=1.0d-7
    If (mean_DL>Wanted_Dl) Then
      eta_min=eta(2)
      Do while (mean_DL>Wanted_Dl)
        eta(2)=eta(2)+0.005d0
        Call rho_ini() !Calculation of m_rho and m_rhoi matrixes
        !Structure Factor Calculation
        Call calc_ck_hd_binary_mixture_roth() !Calculation of the Structure Factor of a Hard Disk mixture System using Roth Expressions
        Call calc_ck_hk_ski()!Calculation for the direct correlation function and total correlation function
        If (Sk_writting_op.Eqv..True.) Then
          Call Sk_writting ("sk_"//char_run//".dat") !Structure Factor Writting
        End If
        If (Gamm_op.Eqv..True.) Then !
          Call Gamma_mixture(gamma_values)  !Call of gamma subroutine
        End If
        Call Long_time_dynamics(i_k_test,i_k_test_number,i1)
        mean_DL=Dl(1)*(m_rho(1,1)**2.d0)+Dl(2)*(m_rho(2,2)**2.d0)
        mean_DL=mean_DL/((m_rho(1,1)**2.d0)+(m_rho(2,2)**2.d0))
        mean_DL=mean_Dl*(eta(1)+(eta(2)*(1.d0/(sigma(2)**3.d0))))/(eta(1)+(eta(2)*(1.d0/(sigma(2)**2.d0))))
        Dl_error=Abs(mean_Dl-Wanted_Dl)/Wanted_Dl
        If(mean_DL>Wanted_Dl) Then
          eta_min=eta(2)
        End If
        If (Dl_error>Dl_tolerance) Then
          Call Dynamics_Variables_mem_deallocation()
        End If
        dtau=1.0d-7
      End Do
      eta_max=eta(2)
    Else
      eta_max=eta(2)
      Do while (mean_DL<Wanted_Dl)
        eta(2)=eta(2)-0.005d0
        Call rho_ini() !Calculation of m_rho and m_rhoi matrixes
        !Structure Factor Calculation
        Call calc_ck_hd_binary_mixture_roth() !Calculation of the Structure Factor of a Hard Disk mixture System using Roth Expressions
        Call calc_ck_hk_ski()!Calculation for the direct correlation function and total correlation function
        If (Sk_writting_op.Eqv..True.) Then
          Call Sk_writting ("sk_"//char_run//".dat") !Structure Factor Writting
        End If
        If (Gamm_op.Eqv..True.) Then !
          Call Gamma_mixture(gamma_values)  !Call of gamma subroutine
        End If
        Call Long_time_dynamics(i_k_test,i_k_test_number,i1)
        mean_DL=Dl(1)*(m_rho(1,1)**2.d0)+Dl(2)*(m_rho(2,2)**2.d0)
        mean_DL=mean_DL/((m_rho(1,1)**2.d0)+(m_rho(2,2)**2.d0))
        mean_DL=mean_Dl*(eta(1)+(eta(2)*(1.d0/(sigma(2)**3.d0))))/(eta(1)+(eta(2)*(1.d0/(sigma(2)**2.d0))))
        Dl_error=Abs(mean_Dl-Wanted_Dl)/Wanted_Dl
        If(mean_DL<Wanted_Dl) Then
          eta_max=eta(2)
        End If
        If (Dl_error>Dl_tolerance) Then
          Call Dynamics_Variables_mem_deallocation()
        End If
        dtau=1.0d-7
      End Do
      eta_min=eta(2)
    End If
    Do While(Dl_error>Dl_tolerance)
      eta(2)=(eta_min+eta_max)/2.d0
      Call rho_ini() !Calculation of m_rho and m_rhoi matrixes
      !Structure Factor Calculation
      Call calc_ck_hd_binary_mixture_roth() !Calculation of the Structure Factor of a Hard Disk mixture System using Roth Expressions
      Call calc_ck_hk_ski()!Calculation for the direct correlation function and total correlation function
      If (Sk_writting_op.Eqv..True.) Then
        Call Sk_writting ("sk_"//char_run//".dat") !Structure Factor Writting
      End If
      If (Gamm_op.Eqv..True.) Then !
        Call Gamma_mixture(gamma_values)  !Call of gamma subroutine
      End If
      Call Long_time_dynamics(i_k_test,i_k_test_number,i1)
      mean_DL=Dl(1)*(m_rho(1,1)**2.d0)+Dl(2)*(m_rho(2,2)**2.d0)
      mean_DL=mean_DL/((m_rho(1,1)**2.d0)+(m_rho(2,2)**2.d0))
      mean_DL=mean_Dl*(eta(1)+(eta(2)*(1.d0/(sigma(2)**3.d0))))/(eta(1)+(eta(2)*(1.d0/(sigma(2)**2.d0))))
      Dl_error=Abs(mean_Dl-Wanted_Dl)/Wanted_Dl
      If(mean_DL<Wanted_Dl) Then
        eta_max=eta(2)
      Else
        eta_min=eta(2)
      End If
      If (Dl_error>Dl_tolerance) Then
        Call Dynamics_Variables_mem_deallocation()
        dtau=1.0d-7
      End If
    End Do
    Print*, "punto", i1
    Print*, "eta 1=",eta(1),"eta 2=",eta(2)
    Call SCGLE_run_info_file("Run_info_"//char_run//".dat",gamma_values,"Hard Disk Mixture",i_k_test(1))
    Write(50,*) eta(1),eta(2),Dl(1),Dl(2),mean_DL
    Print*, "runfile", i1
    flush(50)
    Call Dynamics_Variables_mem_deallocation()
    dtau=1.0d-7
    Print*, "points calculated:", i1
  End Do
  Close(50)
  DEAllocate(i_k_test)
End Subroutine

Subroutine Hard_Disk_attractive_yukawa_mono_dynamics(Gamm_op,Sk_writting_op)
  Implicit None
  Real * 8, Dimension(:,:), Allocatable :: gamma_values
  Real * 8:: gam
  Integer :: i_k_test_number
  Integer, Dimension(:), Allocatable :: i_k_test
  Logical, Intent (in) :: Gamm_op,Sk_writting_op
  !Points assignement
  kpoints=2**13 !Wave vector points
  short_times=10!2**4
  medium_times=2**7
  decimations=2**5
  i_k_test_number=7
  Allocate(i_k_test(i_k_test_number))
  i_k_test(1)=200
  i_k_test(2)=350
  i_k_test(3)=500
  i_k_test(4)=710
  i_k_test(5)=900
  i_k_test(6)=1100
  i_k_test(7)=1300
  !Differential variables assignement
  dk=1.0d-2
  dtau=1.0d-7
  !Parameter System variables assignement and allocation
  SDimen=2  !Space dimension
  Species=1 !Number of species
  Call Sys_Variables_Alloc(species) !Memory Allocation for System Variables
  sigma(1)=distu  !System diameter
  !sigma(2)=0.69d0*sigma(1)
  Call D0M_ini(species,sigma) !Diffusion coefficient of the system particle
  !Independant System Variables
  eta(1)=0.02d0  !Dimensionless Density
  !eta(2)=0.53823d0
  Print *, "Holis"
  Call rho_ini() !Calculation of m_rho and m_rhoi matrixes
  !Static Variables Memory Allocation
  Call Static_Variables_mem_alloc(kpoints,species)
  Call Calc_static_k(dk,kpoints)
  !Structure Factor Calculation
  Call z_yuk_mono(2.d0,1) !Yukawa potential range value assignement
  Call sharma_sharma_attractive_yukawa_2D(0.078d0) !Calculation of the Structure Factor of a Hard Disk+Attractive Yukawa monocomponent System using Roth Expressions
  Call calc_ck_hk_ski()!Calculation for the direct correlation function and total correlation function
  If (Sk_writting_op.Eqv..True.) Then
    Call Sk_writting ("sk.dat") !Structure Factor Writting
  End If
  !Call kc_ini(species,1.15d0)  !Initialization of kc, SCGLE parameter, 1.15 for arrest at 0.70
  Call kc_ini(species,0.7d0)  !Initialization of kc, SCGLE parameter, 0.7 for arrest at 0.72
  Call Calc_lambdak(Species,kpoints)  !Calculation of lambda(k), SCGLE parameter function
  Allocate(gamma_values(Species,Species)) !Allocation of gamma variables
  If (Gamm_op.Eqv..True.) Then !
    Call Gamma_mixture(gamma_values)  !Call of gamma subroutine
    !Stop
  End If
  !Call Long_time_dynamics(i_k_test,i_k_test_number,1)
  !Call SCGLE_run_info_file("Run_info.dat",gamma_values,"Hard Disk Mixture",i_k_test(1))
  DEAllocate(i_k_test)
End Subroutine Hard_Disk_attractive_yukawa_mono_dynamics

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!This is the GOOD!!


Subroutine Hard_Sphere_attractive_yukawa_mono_dynamics(Gamm_op,Sk_writting_op)
  Implicit None
  Real * 8, Dimension(:,:), Allocatable :: gamma_values
  Real * 8:: gam
  Integer :: i_k_test_number
  Integer, Dimension(:), Allocatable :: i_k_test
  Logical, Intent (in) :: Gamm_op,Sk_writting_op
  !Points assignement
  kpoints=199 !Wave vector points
  short_times=10!2**4
  medium_times=2**7
  decimations=2**5
  i_k_test_number=7
  Allocate(i_k_test(i_k_test_number))
  i_k_test(1)=200
  i_k_test(2)=350
  i_k_test(3)=500
  i_k_test(4)=710
  i_k_test(5)=900
  i_k_test(6)=1100
  i_k_test(7)=1300
  !Differential variables assignement
  dk=1.0d-1
  dtau=1.0d-7
  !Parameter System variables assignement and allocation
  SDimen=3  !Space dimension
  Species=1 !Number of species
  Call Sys_Variables_Alloc(species) !Memory Allocation for System Variables
  sigma(1)=distu  !System diameter
  !sigma(2)=0.69d0*sigma(1)
  Call D0M_ini(species,sigma) !Diffusion coefficient of the system particle
  !Independant System Variables
  ! eta(1)=0.5d0  !Dimensionless Density
  eta(1)=0.6d0  !Dimensionless Density
  !eta(2)=0.53823d0
  Print *, "Holis"
  Call rho_ini() !Calculation of m_rho and m_rhoi matrixes
  !Static Variables Memory Allocation
  Call Static_Variables_mem_alloc(kpoints,species) !This is just Mem Alloc
  Call Calc_static_k(dk,kpoints)                   !Make a simple calulate
  !Structure Factor Calculation
  Call z_yuk_mono(2.d0,1) !Yukawa potential range value assignement
  Call sharma_sharma_attractive_yukawa_2D(1.d0) !Calculation of the Structure Factor of a Hard Disk+Attractive Yukawa monocomponent System using Roth Expressions
  Call calc_ck_hk_ski()!Calculation for the direct correlation function and total correlation function
  If (Sk_writting_op.Eqv..True.) Then
    Call Sk_writting ("sk_eta_050.dat") !Structure Factor Writting
  End If
  !Call kc_ini(species,1.15d0)  !Initialization of kc, SCGLE parameter, 1.15 for arrest at 0.70
  Call kc_ini(species,1.305d0)  !Initialization of kc, SCGLE parameter, 0.7 for arrest at 0.72
  Call Calc_lambdak(Species,kpoints)  !Calculation of lambda(k), SCGLE parameter function
  Allocate(gamma_values(Species,Species)) !Allocation of gamma variables
  If (Gamm_op.Eqv..True.) Then !
    Call Gamma_mixture(gamma_values)  !Call of gamma subroutine
    !Stop
  End If
  !Call Long_time_dynamics(i_k_test,i_k_test_number,1)
  !Call SCGLE_run_info_file("Run_info.dat",gamma_values,"Hard Disk Mixture",i_k_test(1))
  DEAllocate(i_k_test)
End Subroutine

End program
