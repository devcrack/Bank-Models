!Module for Parameter Variables
!pi;
Module Parameter_Variables
  Implicit None
  Real * 8, Parameter :: pi=4.d0*DATAN(1.d0)
End Module Parameter_Variables
!!!!!!!!!!
!Module for Dimensional Variables:
!distu;	D0u;
Module Dimensional_Variables
  Implicit None
  Real * 8, Parameter :: distu=1.d0 !Distance unit
  Real * 8, Parameter :: D0u=1.d0 !D0 unit
End Module Dimensional_Variables
!!!!!!!!!!
!Module for System Variables
!SDimen;	Species; eta; sigma; D0M;
Module System_Variables
  Implicit None
  Integer :: SDimen	!Space Dimension
  Integer :: Species !Number of species used
  Real * 8, Dimension (:), Allocatable :: eta	!Dimensionless density of each species
  Real * 8, Dimension (:), Allocatable :: sigma !Dimensionless characteristic distance of each species
  Real * 8, Dimension (:,:), Allocatable :: D0M !Dimensionless Short Time Diffusion Coefficient Diagonal Matrix (species,species)
  Real * 8, Dimension (:,:), Allocatable :: m_rho,m_rhoi !root density matrixes
Contains
  Subroutine Sys_Variables_Alloc(Ialloc)
    Implicit None
    Integer, Intent(in) :: Ialloc
    Allocate(eta(Ialloc))
    Allocate(sigma(Ialloc))
    Allocate(D0M(Ialloc,Ialloc))
    Allocate(m_rho(Ialloc,Ialloc))
    Allocate(m_rhoi(Ialloc,Ialloc))
  End Subroutine Sys_Variables_Alloc

  !Initialization of the D0M variable
  ! call D0M_ini(species,sigma) 
  Subroutine D0M_ini(Iallocated_s,species_distance)
    Implicit None
    Integer :: i1,i2
    Integer, Intent(in) :: Iallocated_s
    Real * 8, Dimension(:), Intent(in) :: species_distance
    Do i1=1, Iallocated_s
       Do i2=1, Iallocated_s
          If (i1==i2) Then
             ! ####################################################
             ! Good Line
             ! print *,' Good Line #1'
             D0M(i1,i2)=1.d0/species_distance(i1)
             ! print *, 'Good value: = ', D0M(i1,i2)
             ! ####################################################
          Else
             D0M(i1,i2)=0.d0
             print *,' Good Line #2'
          End If
       End Do
    End Do
  End Subroutine D0M_ini

  Subroutine rho_ini ()
    Use Parameter_Variables
    Implicit None
    Integer :: i1,i2
    Do i1=1,Species
       Do i2=1,Species
          If (i1==i2)Then
             If (SDimen==3) Then
                !Here execute in hard sphere
                print *, 'Entra aqui__1'
                ! ####################################################
                ! Good Line
                m_rho(i1,i2)=6.d0*eta(i1)/(pi*(sigma(i1)**3))
                ! print *, 'Another Good Value ', m_rho(i1, i2) 
                ! ###################################################
             Else If (SDimen==2) Then
                print *, 'Entra aqui__2'
                m_rho(i1,i2)=4.d0*eta(i1)/(pi*(sigma(i1)**2))
             End if
             print *, 'Obviusly'
             m_rho(i1,i2)=sqrt(m_rho(i1,i2))
             m_rhoi(i1,i2)=1.d0/m_rho(i1,i2)
          Else
             print *, 'Entra aqui__3'
             m_rho(i1,i2)=0.d0
             m_rhoi(i1,i2)=0.d0
          End If
       End Do
    End Do
  End Subroutine rho_ini

End Module System_Variables
!!!!!!!!!!
!Module for Potential Variables
!sw_lambda; z_yuk
Module Potential_Variables
  Use System_Variables
  Implicit None
  Real * 8, Dimension(:), Allocatable :: sw_lambda !Square Well Potential Reach
  Real * 8, Dimension(:), Allocatable :: z_yuk	!Yukawa Potential Reach

Contains
  Subroutine z_yuk_mono(dummy,species_number)
    Implicit None
    Real * 8, Intent(in) :: dummy
    Integer, Intent(in) :: species_number
    Integer :: i1
    Allocate(z_yuk(species_number))
    Do i1=1,species_number
       z_yuk(i1)=dummy
    End Do
  End Subroutine z_yuk_mono
  Subroutine lambda_sw_mono(dummy)
    Implicit None
    Real * 8, Intent(in) :: dummy
    Allocate(sw_lambda(1))
    sw_lambda(1)=dummy
  End Subroutine lambda_sw_mono
End Module Potential_Variables
!!!!!!!!!!
!Module for Static Variables
!kpoints; dk; Sk; Ski; Ck; Hk;
Module Static_Variables
  Implicit None
  Integer :: kpoints !number of wave vectors used
  Real * 8 :: dk !wave vector differential
  Real * 8, Dimension(:), Allocatable :: k 
  Real * 8, Dimension (:,:,:), Allocatable :: Sk,Ski,ck,hk !Structure Factor/Inverse, Direct Correlation Function,Total Correlation Function

Contains
  Subroutine Static_Variables_mem_alloc(Ialloc_k,Ialloc_s)
    Integer, Intent (in) :: Ialloc_k, Ialloc_s
    Allocate(k(Ialloc_k))
    Allocate(Sk(Ialloc_k,Ialloc_s,Ialloc_s))
    Allocate(Ski(Ialloc_k,Ialloc_s,Ialloc_s))
    Allocate(ck(Ialloc_k,Ialloc_s,Ialloc_s))
    Allocate(hk(Ialloc_k,Ialloc_s,Ialloc_s))
  End Subroutine Static_Variables_mem_alloc
End Module Static_Variables
!!!!!!!!!!
!Module for Dynamic Variables
!Fc; Fs; Del_z;	short_times; medium_times; decimations; m_tau
Module Dynamics_Variables
  Implicit None
  Real * 8, Dimension(:,:,:,:), Allocatable :: Fc,Fs !Intermediate scattering function/self...
  Real * 8, Dimension(:,:,:), Allocatable :: Del_z !Time dependant friction function
  Real * 8, Dimension(:,:), Allocatable :: DiffCoeff, MSD
  Real * 8, Dimension(:), Allocatable :: Dl
  Integer :: short_times
  Integer :: medium_times
  Integer :: decimations
  Real * 8, Dimension(:), Allocatable :: m_tau
  Real * 8 :: dtau

Contains
  Subroutine Dynamics_Variables_mem_allocation(time_points,k_points,species_number)
    Implicit None
    Integer, Intent (in) :: time_points,k_points,species_number
    Allocate(Fc(k_points,time_points,species_number,species_number))
    Allocate(Fs(k_points,time_points,species_number,species_number))
    Allocate(Del_z(time_points,species_number,species_number))
    Allocate(m_tau(time_points))
    Allocate(DiffCoeff(time_points,species_number))
    Allocate(MSD(time_points,species_number))
    Allocate(Dl(species_number))
  End Subroutine Dynamics_Variables_mem_allocation

  Subroutine Dynamics_Variables_mem_deallocation()
    Implicit None
    DEAllocate(Fc)
    DEAllocate(Fs)
    DEAllocate(Del_z)
    DEAllocate(m_tau)
    DEAllocate(DiffCoeff)
    DEAllocate(MSD)
    DEAllocate(Dl)
  End Subroutine Dynamics_Variables_mem_deallocation

End Module Dynamics_Variables
!!!!!!!!!!
Module SCGLE_Variables
  Implicit None
  Real * 8, Parameter :: gamma_max=1.0d10 !Max permitted value of gamma
  Real * 8, Dimension(:), Allocatable :: kc
  Real * 8, Dimension(:,:,:), Allocatable :: lambdak
Contains
  !Subroutine for the inicialization of kc, parameter of the SCGLE theory
  Subroutine kc_ini(Ialloc,kcalpha)
    Use System_Variables
    Use Parameter_Variables
    Implicit None
    Integer, Intent(in) :: Ialloc
    Real * 8, Intent(in) :: kcalpha
    Integer :: i1
    Allocate(kc(Ialloc))
    Do i1=1,Ialloc
       kc(i1)=2.0d0*pi*kcalpha/sigma(i1)!the distance is going to be in terms of the big one
    End Do
  End Subroutine kc_ini
  !Subroutine to calculate the parameter function lambda(k) of the SCGLE theory
  Subroutine Calc_lambdak(Ialloc_s,Ialloc_k)
    Use Static_Variables
    Implicit None
    Integer, Intent(in) :: Ialloc_s
    Integer, Intent(in) :: Ialloc_k
    Integer :: i1,i2,i3
    Allocate(lambdak(Ialloc_k,Ialloc_s,Ialloc_s))
    Do i1=1,Ialloc_s
       Do i2=1,Ialloc_s
          Do i3=1,Ialloc_k
             If (i1==i2) Then
                lambdak(i3,i1,i2)=1.d0/(1.d0+((k(i3)/kc(i1))**2.d0))
             Else
                lambdak(i3,i1,i2)=0.d0
             End If
          End Do
       End Do
    End Do
  End Subroutine Calc_lambdak

End Module SCGLE_Variables
!!!!!!!!!!
Module NESCGLE_variables
  Implicit None
  Real * 8, Parameter :: du=1.d-2
  Real * 8, Dimension (:,:,:), Allocatable :: Sk0,Skf !Initial Structure Factor and Final Structure Factor

Contains
  Subroutine NESCGLE_Static_Variables_mem_alloc(Ialloc_k,Ialloc_s)
    Integer, Intent (in) :: Ialloc_k, Ialloc_s
    Allocate(Sk0(Ialloc_k,Ialloc_s,Ialloc_s))
    Allocate(Skf(Ialloc_k,Ialloc_s,Ialloc_s))
  End Subroutine NESCGLE_Static_Variables_mem_alloc

  Subroutine NESCGLE_Static_Variables_mem_dealloc()
    DEAllocate(Sk0)
    DEAllocate(Skf)
  End Subroutine NESCGLE_Static_Variables_mem_dealloc

End Module NESCGLE_variables
