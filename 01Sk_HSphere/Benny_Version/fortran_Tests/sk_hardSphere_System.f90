program test
  Implicit None
    Integer :: i1               !Define an integer variable
    Logical :: vw_option        !Define a boolean variable 
    Real * 8, Dimension(:), Allocatable :: k_vw !Declare  a Dynamic arwray  with tag : k_vw
    Real * 8 :: eta_vw,dum1,dum2,dum3,dum4,dumsin,dumcos !Declare a number of a floating point varibales 
    If (SDimen/= 3 .or. Species /=1 ) Then               !Making Decision, compare values
        print*, ('Error, cannot compute the system direct correlation function') 
        Stop

     End If                     !End of Desicion
     !End of parameter checks
     ! ------------------------------------------------------------------------------------------------------------
     Allocate(k_vw(kpoints))    !Reserve memory for the dynamic Array
                                ! Sintax of this is: Allocatable(identefier(size_array))
     !Start of VW Correction
     If (vw_option .EQV. .TRUE.) Then !Desicion Compare vw_option with true. This instruction can improve
        ! Notes :
        !        EQV logical equivalence
        !        TRUE logical constant
        
        ! --------------------------------------------------------------------
        ! eta(#) is dynamic array that stores the Dimensionless density of each species
        ! To this point eta_have some thing
        ! --------------------------------------------------------------------
        eta_vw=eta(1)*(1.d0-(eta(1)/16.d0)) ! 1.d0 is a double presicion Real
        k_vw(:)=k(:)*((eta_vw/eta(1))**(1.d0/3.d0))
        ! --------------------------------------------------------------------
        ! --------------------------------------------------------------------
     Else
        eta_vw=eta(1)
        k_vw(:)=k(:)
     End If
     


!Subroutine for the calculation of the Structure factor of a Hard Sphere System
  !Using the Percus-Yevick Closure Relationship
  !vw_option Logical value, True for Verlet-Weiss Correction, False Otherwise
  ! Subroutine Calc_Sk_hs_py_mono(vw_option)
  !   Implicit None
  !   Integer :: i1
  !   Logical :: vw_option
  !   Real * 8, Dimension(:), Allocatable :: k_vw
  !   Real * 8 :: eta_vw,dum1,dum2,dum3,dum4,dumsin,dumcos
  !   !Start of parameter checks
  !   If (SDimen/= 3 .or. Species /=1 ) Then
  !       print*, ('Error, cannot compute the system direct correlation function')
  !       Stop
  !   End If
    !!End of parameter checks
    ! Allocate(k_vw(kpoints))     
    !!Start of VW Correction
    If (vw_option .EQV. .TRUE.) Then
      eta_vw=eta(1)*(1.d0-(eta(1)/16.d0))
      k_vw(:)=k(:)*((eta_vw/eta(1))**(1.d0/3.d0))
    Else
      eta_vw=eta(1)
      k_vw(:)=k(:)
    End If
    !End of VW correction
    dum1=-((1.0+2.0*eta_vw)**2)
    dum2=6.0*eta_vw*((1.0+(eta_vw/2.0))**2)
    dum3=-eta_vw*((1.0+(2.0*eta_vw))**2)/2.0
    dum4=((1.0-eta_vw)**4)
    dum1=dum1/dum4
    dum2=dum2/dum4
    dum3=dum3/dum4
    Do i1=1, kpoints
      dumsin=sin(k_vw(i1))
      dumcos=cos(k_vw(i1))
      ck(i1,1,1)= (dum1*(dumsin-k_vw(i1)*dumcos)/(k_vw(i1)**2.d0))+&
      &(dum2*(((2.d0*k_vw(i1))*dumsin)+((-(k_vw(i1)**2.d0)+2.d0)*dumcos)-2.d0)/&
      &(k_vw(i1)**3.d0))+(dum3*(((4.d0*(k_vw(i1)**3.d0)-24.0*k_vw(i1))*dumsin)&
      &+((-(k_vw(i1)**4.d0)+12.d0*(k_vw(i1)**2.d0)-24.d0)*dumcos)+24.d0)/(k_vw(i1)**5))
      Ck(i1,1,1)=  24.d0*eta_vw*ck(i1,1,1)/k_vw(i1)
      Sk(i1,1,1)=1.d0/(1.d0-ck(i1,1,1))
      Ski(i1,1,1)=1.d0-Ck(i1,1,1)
    End Do
  End Subroutine
