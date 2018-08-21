Module Gamma
  Use Variables
  Implicit None

  Contains
  Real * 8 Function Gamma_mono()
    Use Variables
    Use Integration_methods
    Implicit None
    Integer :: i1
    Real * 8 :: gamma_test,dimen_dummy,error,gamma
    Real * 8, Dimension(:), Allocatable :: dum1, dum2, dum3, dum4
    Real * 8, Dimension (:), Allocatable :: integral_vector
    Real * 8, Parameter :: tolerance=1.e-7;
    Allocate(integral_vector(kpoints))
    Allocate(dum1(kpoints))
    Allocate(dum2(kpoints))
    Allocate(dum3(kpoints))
    Allocate(dum4(kpoints))
    If (Dimen==3)Then
      dimen_dummy=36.0*pi*eta(1) !/* 3D */
    Else If (Dimen==2)Then
      dimen_dummy=16.0*eta(1) !/* 2D */
    End If

    Do i1=1,kpoints
      dum1(i1)=(k(i1)*(sk(i1,1,1)-1.d0)*lambdak(i1,1,1))**2
      If (Dimen==3)Then
        dum1(i1)=k(i1)*k(i1)*dum1(i1) !/* 3D */
      Else If (Dimen==2)Then
        dum1(i1)=k(i1)*dum1(i1) !/* 2D */
      End If
      dum2(i1)=lambdak(i1,1,1)*sk(i1,1,1)
      dum3(i1)=k(i1)*k(i1)
    End Do
    gamma_test=1.d-5
    error=1.d0
    gamma=0.d0
    Do While (error>tolerance .AND. gamma<gamma_max)
      Do i1=1,kpoints
        dum4(i1)=dum3(i1)*gamma_test
        integral_vector(i1)=dum1(i1)/((dum2(i1)+dum4(i1))*(lambdak(i1,1,1)+dum4(i1)))
      End Do
      !/*gamma=simpson_3_4_fixed_differential(integral_vector, k_points, dk);*/
      gamma=rectangle_fixed_differential(integral_vector, dk)
      gamma= dimen_dummy/(gamma*dimen_dist*dimen_dist)  ! /* Note, IT IS A HARD SPHERE TERM be care of eta */
      error=abs(gamma_test-gamma)/gamma
      gamma_test=gamma
    !  /*printf("Gamma=%e, error=%e \n",gamma,error);*/
      !Print*, "Gamma=",gamma,"error=",error
    End Do
    Print*,"Gamma=",gamma,"  error=%e \n",error
    Gamma_mono=gamma
    Deallocate (dum1)
    Deallocate (dum2)
    Deallocate (dum3)
    Deallocate (dum4)
    Deallocate (integral_vector)
  End Function

  Subroutine calc_arrest_temperature_sh_sh_SW(vw_option,spinodal_temperature,SW_param,Arrest_temp)
    Use Variables
    Use Structure
    Implicit None
    Real * 8, Intent(in) :: spinodal_temperature,SW_param
    Real * 8 :: tempi,tempf,temp_try,gam_value,error
    Real * 8, Parameter :: precision=1.d0-7
    Real * 8 :: dtemp
    Real * 8, Intent(out) :: Arrest_temp
    Logical, Intent (in) :: vw_option
    dtemp=5.d-1
    tempf=spinodal_temperature+dtemp
    Call calc_ck_sharma_sharma_HS_SW(vw_option,tempf,SW_param)
    gam_value=Gamma_mono()
    If (gam_value>gamma_max)Then
      Print*, "Error, dynamical arrest not avaiable for this volume fraction"
      Stop
    Else
      Do While(gam_value<gamma_max)
        dtemp=dtemp*2.d0
        tempi=tempf
        tempf=tempf+dtemp
        Call calc_ck_sharma_sharma_HS_SW(vw_option,tempf,SW_param)
        gam_value=Gamma_mono()
      End Do
      error=1.d0
      Do While (error>precision)
        dtemp=abs(tempf-tempi)/2.d0
        temp_try=tempi+dtemp
        Call calc_ck_sharma_sharma_HS_SW(vw_option,temp_try,SW_param)
        gam_value=Gamma_mono()
        If(gam_value<gamma_max)Then
          tempi=temp_try
        Else
          tempf=temp_try
        End If
        error=abs(tempf-tempi)/tempf
      End Do
    End If
    Arrest_temp=tempi
  End Subroutine

  Subroutine Gamma_mixture(gamma_values)
    Use Matrix_operations
    Use Variables
    Use Integration_methods
    Implicit None
    Integer :: i1,i2,i3
    Real * 8 :: gamma_test,dimen_dummy,gamma
    Real * 8, Dimension(:,:,:), Allocatable :: dum1_m_k, dum2_m_k, dum3_m_k, dum4_m_k,Integrand_m_k
    Real * 8, Dimension(:,:), Allocatable :: m_dum1,m_dum2,m_dum3,m_dum4,I_m,m_rho_dum
    Real * 8, Dimension(:,:), Intent(out) :: gamma_values
    Real * 8, Dimension(:,:),Allocatable :: gamma_values_test
    Real * 8, Dimension (:), Allocatable :: integral_vector
    Real * 8, Dimension(:), Allocatable :: error
    Real * 8, Parameter :: tolerance=1.d-7
    Logical :: gam_condition,dum_condition
    Allocate(integral_vector(kpoints))
    Allocate(dum1_m_k(kpoints,species,species))
    Allocate(dum2_m_k(kpoints,species,species))
    Allocate(dum3_m_k(kpoints,species,species))
    Allocate(dum4_m_k(kpoints,species,species))
    Allocate(Integrand_m_k(kpoints,species,species))
    Allocate(m_dum1(species,species))
    Allocate(m_dum2(species,species))
    Allocate(m_dum3(species,species))
    Allocate(m_dum4(species,species))
    Allocate(I_m(species,species))
    Allocate(m_rho_dum(species,species))
    Allocate(gamma_values_test(species,species))
    Allocate(error(species))
    If (Dimen==3)Then
      dimen_dummy=6.d0*pi*pi !/* 3D */
    Else If (Dimen==2)Then
      dimen_dummy=4.d0*pi !/* 2D */
    End If
    Call identity_matrix(species,I_m)
    !Initial values of gamma and m_rho_dum calculation
    Do i1=1,species
      Do i2=1,species
        If (i1==i2)Then
          gamma_values(i1,i2)=1.d-5
          If (Dimen==3) Then
            m_rho_dum(i1,i2)=6.d0*eta(i1)/(pi*(sigma(i1)**3))
          Else If (Dimen==2) Then
            m_rho_dum(i1,i2)=4.d0*eta(i1)/(pi*(sigma(i1)**2))
          End if
          m_rho_dum(i1,i2)=sqrt(m_rho_dum(i1,i2))
        Else
          gamma_values(i1,i2)=0.d0
          m_rho_dum(i1,i2)=0.d0
        End If
        Print*, i1,i2,gamma_values(i1,i2)
      End Do
    End Do
    !!
    Do i1=1,kpoints
      Call matrix_conversion_3_to_2(i1,species,lambdak,m_dum1)
      Call inverse_matrix(species,m_dum1,m_dum2)
      Call matrix_constant_multiplication(species,m_dum2,k(i1)**2,m_dum1)
      Call matrix_conversion_2_to_3(i1,species,m_dum1,dum1_m_k)!k^2*lambda^{-1}(k)
      Call matrix_conversion_3_to_2(i1,species,Skinv,m_dum2)
      Call matrix_multiplication(species,m_dum1,m_dum2,m_dum3)
      !Call matrix_constant_multiplication(species,m_dum3,k(i1)**2,m_dum1)
      Call matrix_conversion_2_to_3(i1,species,m_dum3,dum2_m_k) !k^2*lambda^{-1}(k)*S^{-1}(k)
      Call matrix_conversion_3_to_2(i1,species,ck,m_dum1)
      Call matrix_multiplication(species,m_dum1,m_rho_dum,m_dum2)
      Call matrix_conversion_2_to_3(i1,species,m_dum2,dum3_m_k) !c(k)*rho^{1/2}
      Call matrix_conversion_3_to_2(i1,species,sk,m_dum2)
      m_dum3=m_dum2-I_m
      Call inverse_matrix(species,m_rho_dum,m_dum2)
      Call matrix_multiplication(species,m_dum3,m_dum2,m_dum1)
      Call matrix_conversion_2_to_3(i1,species,m_dum1,dum4_m_k) !rho^{1/2}*h(k)
    End Do
    gam_condition=.False. !Inicializador de condicion de gamma
    Do While (gam_condition.Eqv..False.)
      Do i1=1,kpoints
        Call matrix_conversion_3_to_2(i1,species,dum1_m_k,m_dum1)
        Call matrix_multiplication(species,gamma_values,m_dum1,m_dum2)
        m_dum1=I_m+m_dum2
        Call inverse_matrix(species,m_dum1,m_dum2)!First TERM
        Call matrix_conversion_3_to_2(i1,species,dum3_m_k,m_dum1)
        Call matrix_conversion_3_to_2(i1,species,dum2_m_k,m_dum3)
        Call matrix_multiplication(species,gamma_values,m_dum3,m_dum4)
        m_dum3=I_m+m_dum4
        Call inverse_matrix(species,m_dum3,m_dum4)
        Call matrix_multiplication(species,m_dum1,m_dum4,m_dum3)
        Call matrix_conversion_3_to_2(i1,species,dum4_m_k,m_dum1)
        Call matrix_multiplication(species,m_dum3,m_dum1,m_dum4)!Second term
        m_dum3=m_dum2*m_dum4
        If(Dimen==3)Then
          Call matrix_constant_multiplication(species,m_dum3,k(i1)**4,m_dum1)!3D
        Else If (Dimen==2)Then
          Call matrix_constant_multiplication(species,m_dum3,k(i1)**3,m_dum1)!2D
        End If
        Call matrix_conversion_2_to_3(i1,species,m_dum1,Integrand_m_k)
      End Do
      gam_condition=.True.
      Do i1=1,species
        Do i2=1,species
          gamma_values_test(i1,i2)=0.d0
          If(i1==i2)Then
            Call matrix_conversion_3_to_1(i1,i2,kpoints,Integrand_m_k,integral_vector)
            gamma_values_test(i1,i2)=dimen_dummy/rectangle_fixed_differential(integral_vector,dk)
            error(i1)=abs(gamma_values_test(i1,i2)-gamma_values(i1,i2))/gamma_values(i1,i2)
            gamma_values(i1,i2)=gamma_values_test(i1,i2)
            Print*, i1,i2,gamma_values(i1,i2),error(i1),"Holis"
            If (gam_condition.Eqv..True.)Then
              If (error(i1)>tolerance)Then
                gam_condition=.False.
                If (gamma_values(i1,i2)>gamma_max)Then
                  gam_condition=.True.
                End If
              End If
            End If
          End If
        End Do
      End Do
    End Do
    Do i1=1,species
      Do i2=1,species
        If (i1==i2)Then
          Print*, i1,i2, "Gamma value=", gamma_values(i1,i2), "error=",error(i1)
        End If
      End Do
    End Do
    DEAllocate(integral_vector)
    DEAllocate(dum1_m_k)
    DEAllocate(dum2_m_k)
    DEAllocate(dum3_m_k)
    DEAllocate(dum4_m_k)
    DEAllocate(Integrand_m_k)
    DEAllocate(m_dum1)
    DEAllocate(m_dum2)
    DEAllocate(m_dum3)
    DEAllocate(m_dum4)
    DEAllocate(I_m)
    DEAllocate(m_rho_dum)
  End Subroutine

  Subroutine calc_arrest_temperature_sh_sh_SW_mono_system_out_spinodal()
    Use Variables
    Implicit None
    Real * 8 :: Tempi,Tempf,dTemp
    Real * 8, Parameter :: Tolerance=1.0d-7

  End Subroutine

End Module
