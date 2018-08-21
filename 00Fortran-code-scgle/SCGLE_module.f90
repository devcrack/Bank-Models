Module SCGLE
  Use Parameter_Variables
  Use System_Variables
  Use Dimensional_Variables
  Use Static_Variables
  Use SCGLE_Variables
  Use Structure_module
  Use Matrix_operations
  Use Integration_methods
  Use Dynamics_Variables
  Use Writting
  Implicit None
  Real * 8, Dimension(:), Allocatable :: del_z_dimen_dummy

  Contains
    Real * 8 Function Gamma_mono()
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
      If (SDimen==3)Then
        dimen_dummy=36.0*pi*eta(1) !/* 3D */
      Else If (SDimen==2)Then
        dimen_dummy=16.0*eta(1) !/* 2D */
      End If

      Do i1=1,kpoints
        dum1(i1)=(k(i1)*(sk(i1,1,1)-1.d0)*lambdak(i1,1,1))**2
        If (SDimen==3)Then
          dum1(i1)=k(i1)*k(i1)*dum1(i1) !/* 3D */
        Else If (SDimen==2)Then
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
        gamma= dimen_dummy/(gamma*distu*distu)  ! /* Note, IT IS A HARD SPHERE TERM be care of eta */
        error=abs(gamma_test-gamma)/gamma
        !print*, gamma,gamma_test,error
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

    Subroutine Gamma_mixture(gamma_values)
      !NEEDED VARIABLES: Sk, Ski, ck, hk, lambdak, k, m_rho
      Implicit None
      Integer :: i1,i2,i3
      Real * 8 :: dimen_dummy
      Real * 8, Dimension(:,:,:), Allocatable :: dum1_m_k, dum2_m_k, dum3_m_k, dum4_m_k,Integrand_m_k
      Real * 8, Dimension(:,:), Allocatable :: m_dum1,m_dum2,m_dum3,m_dum4,I_m
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
      Allocate(gamma_values_test(species,species))
      Allocate(error(species))
      If (SDimen==3)Then
        dimen_dummy=(4.d0*pi/(3.d0*((2.d0*pi)**3.d0))) !/* 3D */
      Else If (SDimen==2)Then
        dimen_dummy=(2.d0*pi/(2.d0*((2.d0*pi)**2.d0))) !/* 2D */
      End If
      dimen_dummy=1.d0/dimen_dummy
      Call identity_matrix(species,I_m)
      !Initial values of gamma calculation
      Do i1=1,Species
        Do i2=1,Species
          If (i1==i2)Then
            gamma_values(i1,i2)=1.d-5
          Else
            gamma_values(i1,i2)=0.d0
          End If
          Print*, i1,i2,gamma_values(i1,i2)
        End Do
      End Do
      !!
      !$OMP DO Private(m_dum1,m_dum2,m_dum3,i1)
      Do i1=1,kpoints
        Call matrix_conversion_3_to_2(i1,species,lambdak,m_dum1)
        Call inverse_matrix(species,m_dum1,m_dum2)
        Call matrix_constant_multiplication(species,m_dum2,k(i1)**2,m_dum1)
        Call matrix_conversion_2_to_3(i1,species,m_dum1,dum1_m_k)!k^2*lambda^{-1}(k)
        Call matrix_conversion_3_to_2(i1,species,Ski,m_dum2)
        Call matrix_multiplication(species,m_dum1,m_dum2,m_dum3)
        !Call matrix_constant_multiplication(species,m_dum3,k(i1)**2,m_dum1)
        Call matrix_conversion_2_to_3(i1,species,m_dum3,dum2_m_k) !k^2*lambda^{-1}(k)*S^{-1}(k)
        Call matrix_conversion_3_to_2(i1,species,ck,m_dum1)
        !m_dum1=I_m-m_dum2
        !Call inverse_matrix(species,m_rho,m_dum2)
        !Call matrix_multiplication(species,m_dum1,m_dum2,m_dum3)
        Call matrix_multiplication(Species,m_dum1,m_rho,m_dum3)
        Call matrix_conversion_2_to_3(i1,species,m_dum3,dum3_m_k) !c(k)*rho^{1/2}
        !Call matrix_conversion_3_to_2(i1,species,sk,m_dum2)
        !m_dum3=m_dum2-I_m
        !Call inverse_matrix(species,m_rho_dum,m_dum2)
        Call matrix_conversion_3_to_2(i1,Species,hk,m_dum2)
        Call matrix_multiplication(species,m_rho,m_dum2,m_dum1)
        Call matrix_conversion_2_to_3(i1,species,m_dum1,dum4_m_k) !rho^{1/2}*h(k)
      End Do
      !$OMP END DO
      gam_condition=.False. !Inicializador de condicion de gamma
      Do While (gam_condition.Eqv..False.)
        !$OMP DO Private(m_dum1,m_dum2,m_dum3,m_dum4,i1)
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
          If(SDimen==3)Then
            Call matrix_constant_multiplication(species,m_dum3,k(i1)**4,m_dum1)!3D
          Else If (SDimen==2)Then
            Call matrix_constant_multiplication(species,m_dum3,k(i1)**3,m_dum1)!2D
          End If
          Call matrix_conversion_2_to_3(i1,species,m_dum1,Integrand_m_k)
        End Do
        !$OMP END DO
        gam_condition=.True.
        Do i1=1,species
          Do i2=1,species
            gamma_values_test(i1,i2)=0.d0
            If(i1==i2)Then
              Call matrix_conversion_3_to_1(i1,i2,kpoints,Integrand_m_k,integral_vector)
              gamma_values_test(i1,i2)=dimen_dummy/rectangle_fixed_differential(integral_vector,dk)
              error(i1)=abs(gamma_values_test(i1,i2)-gamma_values(i1,i2))/gamma_values_test(i1,i2)
              !Print*, i1,i2,gamma_values(i1,i2),error(i1),gamma_values_test(i1,i2),"Holis"
              gamma_values(i1,i2)=gamma_values_test(i1,i2)
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
      DEAllocate(gamma_values_test)
      DEAllocate(error)
    End Subroutine

    Subroutine Calc_m_times(delta_tau,i_times)
      Implicit None
      Real * 8, Intent(in) :: delta_tau
      Integer, Intent(in) :: i_times
      Integer :: i1
      Do i1=1, i_times
        m_tau(i1)=delta_tau*i1
      End Do
    End Subroutine

    Subroutine calc_Del_z_dimen_dummy()
      Implicit None
      Integer :: i1
      Allocate(del_z_dimen_dummy(kpoints))
      If(SDimen==3) Then
        Do i1=1,kpoints
          Del_z_dimen_dummy(i1)=(4.d0*pi*(k(i1)**4.d0))/(3.d0*((2.d0*pi)**3.d0))
        End Do
      Else If(SDimen==2) Then
        Do i1=1,kpoints
          Del_z_dimen_dummy(i1)=(2.d0*pi*(k(i1)**3.d0))/(2.d0*((2.d0*pi)**2.d0))
        End Do
      End If
    End Subroutine

    Subroutine calc_Del_z_i(i_time,m_dum_Del_z)
      Implicit None
      Integer, Intent(in) :: i_time
      Real * 8, Dimension (:,:), Intent(out) :: m_dum_Del_z
      Real * 8, Dimension (:,:), Allocatable :: m_dum1,m_dum2,m_dum3
      Real * 8, Dimension (:,:), Allocatable :: m_dumS,m_dumSi,m_dumc,m_dumh,m_dumFc,m_dumFs
      Real * 8, Dimension (:,:,:), Allocatable :: m_Integrand
      Real * 8, Dimension (:), Allocatable :: v_Integrand
      Integer :: i1,i2
      Allocate(m_dum1(Species,Species))
      Allocate(m_dum2(Species,Species))
      Allocate(m_dum3(Species,Species))
      Allocate(m_dumS(Species,Species))
      Allocate(m_dumSi(Species,Species))
      Allocate(m_dumh(Species,Species))
      Allocate(m_dumc(Species,Species))
      Allocate(m_dumFs(Species,Species))
      Allocate(m_dumFc(Species,Species))
      Allocate(m_Integrand(kpoints,Species,Species))
      Allocate(v_Integrand(kpoints))
      Do i1=1,kpoints
        Call matrix_conversion_4_to_2(i1,i_time,Species,Fc,m_dumFc)
        Call matrix_conversion_4_to_2(i1,i_time,Species,Fs,m_dumFs)
        Call matrix_conversion_3_to_2(i1,Species,Sk,m_dumS)
        Call matrix_conversion_3_to_2(i1,Species,ck,m_dumc)
        Call matrix_conversion_3_to_2(i1,Species,hk,m_dumh)
        Call inverse_matrix(Species,m_dumS,m_dumSi)
        Call matrix_multiplication(Species,m_dumc,m_rho,m_dum1)
        Call matrix_multiplication(Species,m_dum1,m_dumFc,m_dum2)
        Call matrix_multiplication(Species,m_dum2,m_dumSi,m_dum1)
        Call matrix_multiplication(Species,m_dum1,m_rho,m_dum2)
        Call matrix_multiplication(Species,m_dum2,m_dumh,m_dum1)
        m_dum2=D0M*m_dumFs*m_dum1
        Call matrix_constant_multiplication(Species,m_dum2,del_z_dimen_dummy(i1),m_dum1)
        Call matrix_conversion_2_to_3(i1,Species,m_dum1,m_Integrand)
      End Do
      Do i1=1,Species
        Do i2=1,Species
          If (i1==i2) Then
            Call matrix_conversion_3_to_1(i1,i2,kpoints,m_Integrand,v_Integrand)
            m_dum_Del_z(i1,i2)=rectangle_fixed_differential(v_Integrand,dk)
          Else
            m_dum_Del_z(i1,i2)=0.d0
          End If
        End Do
      End Do
      DEAllocate(m_Integrand)
      DEAllocate(v_Integrand)
      DEAllocate(m_dumFs)
      DEAllocate(m_dumFc)
      DEAllocate(m_dumh)
      DEAllocate(m_dumc)
      DEAllocate(m_dumS)
      DEAllocate(m_dumSi)
      DEAllocate(m_dum1)
      DEAllocate(m_dum2)
      DEAllocate(m_dum3)
    End Subroutine

    Subroutine Short_times_dynamics()
      Implicit None
      Real * 8, Dimension(:,:), Allocatable :: m_dum1,m_dum2,m_dum3,m_dum4,m_dumSi,m_dumFc,m_dumFs,m_dumS
      Real * 8, Dimension(:,:,:), Allocatable :: m_Integrand
      Real * 8, Dimension(:), Allocatable :: v_Integrand
      Real * 8, Dimension(:), Allocatable  :: v_dum1
      Integer :: i1,i2,i3,i4
      Real * 8, Dimension(:,:), Allocatable :: I_m
      Allocate(m_dum1(Species,Species))
      Allocate(m_dum2(Species,Species))
      Allocate(m_dum3(Species,Species))
      Allocate(m_dum4(Species,Species))
      Allocate(m_dumSi(Species,Species))
      Allocate(I_m(Species,Species))
      Allocate(m_Integrand(kpoints,Species,Species))
      Allocate(v_Integrand(kpoints))
      Allocate(m_dumFc(Species,Species))
      Allocate(m_dumFs(Species,Species))
      Allocate(m_dumS(Species,Species))
      Allocate(v_dum1(Species))
      Call identity_matrix(Species,I_m)
      v_dum1=0.d0
      Do i1=1, short_times
        Do i2=1,kpoints
          If (Species>=1)Then
            Call matrix_conversion_3_to_2(i2,Species,Sk,m_dumS)
            Call inverse_matrix(Species,m_dumS,m_dumSi)
            Call matrix_constant_multiplication(Species,D0M,-k(i2)*k(i2)*m_tau(i1),m_dum4)
            m_dumFc=m_dumS+m_dum4
            Call matrix_conversion_2_to_4(i2,i1,Species,m_dumFc,Fc)!Fc
            m_dumFs=I_m+m_dum4
            Call matrix_conversion_2_to_4(i2,i1,Species,m_dumFs,Fs)!Fs
          !Else If(Species==1)Then
          !  Fc(i2,i1,1,1)=Sk(i2,1,1)*exp(-k(i2)*k(i2)*D0M(1,1)*m_tau(i1)/Sk(i2,1,1))
          !  Fs(i2,i1,1,1)=exp(-k(i2)*k(i2)*D0M(1,1)*m_tau(i1)/Sk(i2,1,1))
          End If
        End Do
        Call calc_Del_z_i(i1,m_dum1)
        Call matrix_conversion_2_to_3(i1,Species,m_dum1,Del_z)
        Do i2=1,Species
          v_dum1(i2)=m_dum1(i2,i2)+v_dum1(i2)
          DiffCoeff(i1,i2)=D0M(i2,i2)-(dtau*v_dum1(i2))
          If (i1==1) Then
            MSD(i1,i2)=dtau*DiffCoeff(i1,i2)
            Dl(i2)=dtau*m_dum1(i2,i2)
          Else
            MSD(i1,i2)=MSD(i1-1,i2)+(dtau*DiffCoeff(i1,i2))
            Dl(i2)=Dl(i2)+(dtau*m_dum1(i2,i2))
          End If
        End Do
      End Do
      DEAllocate(v_dum1)
      DEAllocate(m_dum1)
      DEAllocate(m_dum2)
      DEAllocate(m_dum3)
      DEAllocate(m_dum4)
      DEAllocate(m_dumSi)
      DEAllocate(I_m)
      DEAllocate(m_Integrand)
      DEAllocate(v_Integrand)
      DEAllocate(m_dumFc)
      DEAllocate(m_dumFs)
      DEAllocate(m_dumS)
    End Subroutine

    Subroutine calc_medium_times_diffusion_coefficient(i_time)
      Implicit None
      Integer, Intent(in) :: i_time
      Integer :: i1,i2,n2
      Real * 8 :: dum1,dum2,dum3
      n2=i_time/2
      Do i1=1,Species
        dum1=0.d0
        dum2=0.d0
        Do i2=1,n2
          dum1=dum1+(DiffCoeff(i2,i1)*(Del_z(i_time+1-i2,i1,i1)-Del_z(i_time-i2,i1,i1)))
          If (i2>1) Then
            dum2=dum2+((DiffCoeff(i_time+1-i2,i1)-DiffCoeff(i_time-i2,i1))*Del_z(i2,i1,i1))
          End If
        End Do
        dum3=dum1+dum2+(DiffCoeff(n2,i1)*Del_z(n2,i1,i1))-(DiffCoeff(i_time-1,i1)*((1.d0/dtau)+Del_z(1,i1,i1)))
        dum3=-dum3/((1.d0/dtau)+Del_z(1,i1,i1))
        DiffCoeff(i_time,i1)=dum3
      End Do
    End Subroutine

    Subroutine Medium_times_dynamics(Decimation_number,i_time_begin)
      Implicit None
      Real * 8, Dimension(:,:), Allocatable :: m_dum1,m_dum2,m_dum3,m_dum4,m_dum_lambdak
      Real * 8, Dimension(:,:), Allocatable :: m_dum_Fc1,m_dum_Fs1,m_dum_Del_z1
      Real * 8, Dimension(:,:), Allocatable :: m_dum_Fcn1,m_dum_Fsn1,m_dum_Del_zn1
      Real * 8, Dimension(:,:), Allocatable :: m_dum_Fcn2,m_dum_Fsn2,m_dum_Del_zn2
      Real * 8, Dimension(:,:), Allocatable :: m_dum_Fcni,m_dum_Fsni,m_dum_Del_zni
      Real * 8, Dimension(:,:), Allocatable :: m_dum_Fcni1,m_dum_Fsni1,m_dum_Del_zni1
      Real * 8, Dimension(:,:), Allocatable :: m_dum_Fci,m_dum_Fsi,m_dum_Del_zi
      Real * 8, Dimension(:,:), Allocatable :: m_dumSi,m_dumFc,m_dumFs,m_dumS,m_dum_alph_c,m_dum_alph_s
      Real * 8, Dimension(:,:), Allocatable :: m_dumFc_test,m_dumFs_test,m_dum_ck,m_dum_hk
      Real * 8, Dimension(:,:,:), Allocatable :: m_Integrand,alpha_c,alpha_s
      Real * 8, Dimension(:,:,:), Allocatable :: Fcn_dum,Fsn_dum
      Real * 8, Dimension(:), Allocatable :: v_Integrand
      Real * 8, Dimension(:), Allocatable :: error
      Real * 8, Dimension(:,:), Allocatable :: Del_z_test
      Real * 8, Parameter :: tolerance=1.0d-6
      Logical :: Error_condition
      Integer :: i1,i2,i3,i4,n2
      Integer, Intent(in) :: i_time_begin
      Integer, Intent(in) :: Decimation_number
      Real * 8, Dimension(:,:), Allocatable :: I_m
      Allocate(m_dumFc_test(Species,Species))
      Allocate(m_dumFs_test(Species,Species))
      Allocate(error(Species))
      Allocate(Del_z_test(Species,Species))
      Allocate(Fcn_dum(kpoints,Species,Species))
      Allocate(Fsn_dum(kpoints,Species,Species))
      Allocate(m_dum1(Species,Species))
      Allocate(m_dum2(Species,Species))
      Allocate(m_dum3(Species,Species))
      Allocate(m_dum4(Species,Species))
      Allocate(m_dumSi(Species,Species))
      Allocate(I_m(Species,Species))
      Allocate(m_Integrand(kpoints,Species,Species))
      Allocate(v_Integrand(kpoints))
      Allocate(m_dumFc(Species,Species))
      Allocate(m_dumFs(Species,Species))
      Allocate(m_dumS(Species,Species))
      Allocate(m_dum_alph_c(Species,Species))
      Allocate(m_dum_alph_s(Species,Species))
      Allocate(m_dum_Fc1(Species,Species))
      Allocate(m_dum_Fs1(Species,Species))
      Allocate(m_dum_Del_z1(Species,Species))
      Allocate(m_dum_Fcn1(Species,Species))
      Allocate(m_dum_Fsn1(Species,Species))
      Allocate(m_dum_Del_zn1(Species,Species))
      Allocate(m_dum_Fcn2(Species,Species))
      Allocate(m_dum_Fsn2(Species,Species))
      Allocate(m_dum_Del_zn2(Species,Species))
      Allocate(m_dum_Fcni1(Species,Species))
      Allocate(m_dum_Fsni1(Species,Species))
      Allocate(m_dum_Del_zni1(Species,Species))
      Allocate(m_dum_Fci(Species,Species))
      Allocate(m_dum_Fsi(Species,Species))
      Allocate(m_dum_Del_zi(Species,Species))
      Allocate(m_dum_Fcni(Species,Species))
      Allocate(m_dum_Fsni(Species,Species))
      Allocate(m_dum_Del_zni(Species,Species))
      Allocate(m_dum_lambdak(Species,Species))
      Call identity_matrix(Species,I_m) !Identity Matrix
      !Alpha c and Alpha s, dummy matrixes calculation
      Allocate(alpha_c(kpoints,Species,Species))
      Allocate(alpha_s(kpoints,Species,Species))
      Call matrix_constant_multiplication(Species,I_m,1.d0/dtau,m_dum1)
      !$OMP DO Private(m_dum1,m_dum2,m_dum3,m_dum4,m_dum_alph_c,m_dum_alph_s,m_dumS,m_dumSi)
      Do i1=1,kpoints
        Call matrix_conversion_3_to_2(i1,Species,lambdak,m_dum2)
        Call matrix_conversion_3_to_2(1,Species,Del_z,m_dum3)
        Call matrix_conversion_3_to_2(i1,Species,Sk,m_dumS)
        Call inverse_matrix(Species,m_dumS,m_dumSi)
        Call matrix_multiplication(Species,m_dum2,m_dum3,m_dum4)
        m_dum_alph_c=m_dum4+m_dum1
        m_dum_alph_s=m_dum4+m_dum1
        Call matrix_constant_multiplication(Species,D0M,k(i1)*k(i1),m_dum2)
        Call matrix_multiplication(Species,m_dum2,m_dumSi,m_dum3)
        m_dum_alph_c=m_dum_alph_c+m_dum3
        m_dum_alph_s=m_dum_alph_s+m_dum2
        Call inverse_matrix(Species,m_dum_alph_c,m_dum2)
        Call inverse_matrix(Species,m_dum_alph_s,m_dum3)
        Call matrix_conversion_2_to_3(i1,Species,m_dum2,alpha_c)
        Call matrix_conversion_2_to_3(i1,Species,m_dum3,alpha_s)
      End Do
      !$OMP END DO
      Do i1=i_time_begin, medium_times
        n2=i1/2
        !$OMP DO Private(m_dum_lambdak,m_dumS,m_dum_alph_c,m_dum_alph_s,m_dumSi,&
        !$OMP& m_dumS,m_dum_Fc1,m_dum_Fs1,m_dum_Del_z1,m_dum_Fcn1,m_dum_Fsn1,m_dum_Del_zn1,&
        !$OMP& m_dum_Fcn2,m_dum_Fsn2,m_dum_Del_zn2,m_dum1,m_dum2,m_dum3,m_dum4,&
        !$OMP& m_dum_Fci,m_dum_Fsi,m_dum_Del_zi,m_dum_Fcni,m_dum_Fsni,m_dum_Del_zni,&
        !$OMP& m_dum_Fcni1,m_dum_Fsni1,m_dum_Del_zni1,m_dumFc,m_dumFs,i2,i3)
        Do i2=1,kpoints
          Call matrix_conversion_3_to_2(i2,Species,lambdak,m_dum_lambdak)
          Call matrix_conversion_3_to_2(i2,Species,Sk,m_dumS)
          Call matrix_conversion_3_to_2(i2,Species,alpha_c,m_dum_alph_c)
          Call matrix_conversion_3_to_2(i2,Species,alpha_s,m_dum_alph_s)
          Call inverse_matrix(Species,m_dumS,m_dumSi)
          !tau(1) variables
          Call matrix_conversion_4_to_2(i2,1,Species,Fc,m_dum_Fc1)
          Call matrix_conversion_4_to_2(i2,1,Species,Fs,m_dum_Fs1)
          Call matrix_conversion_3_to_2(1,Species,Del_z,m_dum_Del_z1)
          !tau(n-1) variables
          Call matrix_conversion_4_to_2(i2,i1-1,Species,Fc,m_dum_Fcn1)
          Call matrix_conversion_4_to_2(i2,i1-1,Species,Fs,m_dum_Fsn1)
          Call matrix_conversion_3_to_2(i1-1,Species,Del_z,m_dum_Del_zn1)
          !tau(n2) variables
          Call matrix_conversion_4_to_2(i2,n2,Species,Fc,m_dum_Fcn2)
          Call matrix_conversion_4_to_2(i2,n2,Species,Fs,m_dum_Fsn2)
          Call matrix_conversion_3_to_2(n2,Species,Del_z,m_dum_Del_zn2)
          !Operations
          Call matrix_multiplication(Species,m_dum_Del_zn1,m_dum_Fc1,m_dum1)
          Call matrix_multiplication(Species,m_dum_Del_z1,m_dum_Fcn1,m_dum2)
          Call matrix_multiplication(Species,m_dum_Del_zn2,m_dum_Fcn2,m_dum3)
          m_dumFc=m_dum1+m_dum2-m_dum3!inner brackets terms of Fc without the sum term
          Call matrix_multiplication(Species,m_dum_Del_zn1,m_dum_Fs1,m_dum1)
          Call matrix_multiplication(Species,m_dum_Del_z1,m_dum_Fsn1,m_dum2)
          Call matrix_multiplication(Species,m_dum_Del_zn2,m_dum_Fsn2,m_dum3)
          m_dumFs=m_dum1+m_dum2-m_dum3!inner brackets terms of Fs without the sum term
          Do i3=2,n2
            !tau(i3) variables
            Call matrix_conversion_4_to_2(i2,i3,Species,Fc,m_dum_Fci)
            Call matrix_conversion_4_to_2(i2,i3,Species,Fs,m_dum_Fsi)
            Call matrix_conversion_3_to_2(i3,Species,Del_z,m_dum_Del_zi)
            !tau(i1-i3) variables
            Call matrix_conversion_4_to_2(i2,i1-i3,Species,Fc,m_dum_Fcni)
            Call matrix_conversion_4_to_2(i2,i1-i3,Species,Fs,m_dum_Fsni)
            Call matrix_conversion_3_to_2(i1-i3,Species,Del_z,m_dum_Del_zni)
            !tau(i1-i3+1)
            Call matrix_conversion_4_to_2(i2,i1-i3+1,Species,Fc,m_dum_Fcni1)
            Call matrix_conversion_4_to_2(i2,i1-i3+1,Species,Fs,m_dum_Fsni1)
            Call matrix_conversion_3_to_2(i1-i3+1,Species,Del_z,m_dum_Del_zni1)
            !Operations
            m_dum1=m_dum_Del_zni1-m_dum_Del_zni
            m_dum2=m_dum_Fcni1-m_dum_Fcni
            Call matrix_multiplication(Species,m_dum1,m_dum_Fci,m_dum3)
            Call matrix_multiplication(Species,m_dum_Del_zi,m_dum2,m_dum4)
            m_dumFc=m_dumFc-m_dum3-m_dum4!All inner brackets terms of Fc
            m_dum2=m_dum_Fsni1-m_dum_Fsni
            Call matrix_multiplication(Species,m_dum1,m_dum_Fsi,m_dum3)
            Call matrix_multiplication(Species,m_dum_Del_zi,m_dum2,m_dum4)
            m_dumFs=m_dumFs-m_dum3-m_dum4!All inner brackets terms of Fs
          End Do
          Call matrix_multiplication(Species,m_dum_lambdak,m_dumFc,m_dum1)
          Call matrix_multiplication(Species,m_dum_lambdak,m_dumFs,m_dum2)
          Call matrix_constant_multiplication(Species,m_dum_Fcn1,1.d0/dtau,m_dum3)
          Call matrix_constant_multiplication(Species,m_dum_Fsn1,1.d0/dtau,m_dum4)
          m_dumFc=m_dum1+m_dum3
          m_dumFs=m_dum2+m_dum4
          Call matrix_multiplication(Species,m_dum_alph_c,m_dumFc,m_dum1)!Fc(k=i2,tau=i1) all non iterative terms
          Call matrix_multiplication(Species,m_dum_alph_s,m_dumFs,m_dum2)!Fs(k=i2,tau=i1) all non iterative terms
          Call matrix_conversion_2_to_3(i2,Species,m_dum1,Fcn_dum)
          Call matrix_conversion_2_to_3(i2,Species,m_dum2,Fsn_dum)
        End Do
        !$OMP END DO
        !Print *, "holis 2"
        !Begin of iteration
        !$OMP DO Private(m_dum_lambdak,m_dumS,m_dum_alph_c,m_dum_alph_s,m_dumSi,&
        !$OMP& m_dum_Fc1,m_dum_Fs1,m_dumFc,m_dumFs,m_dum1,m_dum2,m_dum3,i2)
        Error_condition=.False.!Initialization of condition
        Call matrix_conversion_3_to_2(i1-1,Species,Del_z,Del_z_test)
        Do While (Error_condition.Eqv..False.)
          Do i2=1, kpoints
            Call matrix_conversion_3_to_2(i2,Species,lambdak,m_dum_lambdak)
            Call matrix_conversion_3_to_2(i2,Species,Sk,m_dumS)
            Call matrix_conversion_3_to_2(i2,Species,alpha_c,m_dum_alph_c)
            Call matrix_conversion_3_to_2(i2,Species,alpha_s,m_dum_alph_s)
            Call inverse_matrix(Species,m_dumS,m_dumSi)
            Call matrix_conversion_4_to_2(i2,1,Species,Fc,m_dum_Fc1)
            Call matrix_conversion_4_to_2(i2,1,Species,Fs,m_dum_Fs1)
            Call matrix_conversion_3_to_2(i2,Species,Fcn_dum,m_dumFc)
            Call matrix_conversion_3_to_2(i2,Species,Fsn_dum,m_dumFs)
            Call matrix_multiplication(Species,m_dum_alph_c,m_dum_lambdak,m_dum1)
            Call matrix_multiplication(Species,m_dum1,Del_z_test,m_dum2)
            m_dum1=m_dumS-m_dum_Fc1
            Call matrix_multiplication(Species,m_dum2,m_dum1,m_dum3)
            m_dumFc_test=m_dum3+m_dumFc!Fc(k(i2),t(i1)) test
            Call matrix_multiplication(Species,m_dum_alph_s,m_dum_lambdak,m_dum1)
            Call matrix_multiplication(Species,m_dum1,Del_z_test,m_dum2)
            m_dum1=I_m-m_dum_Fs1
            Call matrix_multiplication(Species,m_dum2,m_dum1,m_dum3)
            m_dumFs_test=m_dum3+m_dumFs!Fs(k(i2),t(i1)) test
            Call matrix_conversion_2_to_4(i2,i1,Species,m_dumFc_test,Fc)
            Call matrix_conversion_2_to_4(i2,i1,Species,m_dumFs_test,Fs)
          End Do
          !$OMP END DO
          Error_condition=.True.
          Call calc_Del_z_i(i1,m_dum1)
          Call matrix_conversion_2_to_3(i1,Species,m_dum1,Del_z)
          Do i2=1,Species
            error(i2)=Abs((Del_z_test(i2,i2)-Del_z(i1,i2,i2))/Del_z(i1,i2,i2))
            !Print*, Del_z(i1,i2,i2),i1, error(i2),tolerance
            If (error(i2)>tolerance) Then
              Error_condition=.False.
            End If
            Del_z_test(i2,i2)=Del_z(i1,i2,i2)
          End Do
        End Do
        Call calc_medium_times_diffusion_coefficient(i1)
        Do i2=1,Species
          MSD(i1,i2)=MSD(i1-1,i2)+(dtau*DiffCoeff(i1,i2))
          Dl(i2)=Dl(i2)+(dtau*m_dum1(i2,i2))
          !Print*, i1,MSD(i1,i2)
        End Do
      End Do
      DEAllocate(m_dumFc_test)
      DEAllocate(m_dumFs_test)
      DEAllocate(error)
      DEAllocate(Del_z_test)
      DEAllocate(Fcn_dum)
      DEAllocate(Fsn_dum)
      DEAllocate(m_dum1)
      DEAllocate(m_dum2)
      DEAllocate(m_dum3)
      DEAllocate(m_dum4)
      DEAllocate(m_dumSi)
      DEAllocate(I_m)
      DEAllocate(m_Integrand)
      DEAllocate(v_Integrand)
      DEAllocate(m_dumFc)
      DEAllocate(m_dumFs)
      DEAllocate(m_dumS)
      DEAllocate(m_dum_alph_c)
      DEAllocate(m_dum_alph_s)
      DEAllocate(m_dum_Fc1)
      DEAllocate(m_dum_Fs1)
      DEAllocate(m_dum_Del_z1)
      DEAllocate(m_dum_Fcn1)
      DEAllocate(m_dum_Fsn1)
      DEAllocate(m_dum_Del_zn1)
      DEAllocate(m_dum_Fcn2)
      DEAllocate(m_dum_Fsn2)
      DEAllocate(m_dum_Del_zn2)
      DEAllocate(m_dum_Fcni1)
      DEAllocate(m_dum_Fsni1)
      DEAllocate(m_dum_Del_zni1)
      DEAllocate(m_dum_Fci)
      DEAllocate(m_dum_Fsi)
      DEAllocate(m_dum_Del_zi)
      DEAllocate(m_dum_Fcni)
      DEAllocate(m_dum_Fsni)
      DEAllocate(m_dum_Del_zni)
      DEAllocate(m_dum_lambdak)
    End Subroutine

    Subroutine half_times_save()
      Implicit None
      Integer :: i1,i2,i3,i4
      dtau=2.d0*dtau
      Call Calc_m_times(dtau,medium_times)
      Do i1=1, Species
        Do i2=1, Species
          Do i3=1, medium_times/2
            Do i4=1, kpoints
              Fc(i4,i3,i1,i2)=Fc(i4,i3*2,i1,i2)
              Fs(i4,i3,i1,i2)=Fs(i4,i3*2,i1,i2)
            End Do
            Del_z(i3,i1,i2)=Del_z(i3*2,i1,i2)
          End Do
        End Do
      End Do
      Do i1=1, medium_times/2
        Do i2=1, Species
          DiffCoeff(i1,i2)=DiffCoeff(i1*2,i2)
          MSD(i1,i2)=MSD(i1*2,i2)
        End Do
      End Do
    End Subroutine

    Subroutine  Long_time_dynamics(i_k_test,number_i_k_test,i_run)
      Implicit None
      Integer, Intent(in) :: number_i_k_test,i_run
      Integer, Dimension(:), Intent(in) :: i_k_test
      Integer :: i1,i2
      Character (LEN=2) :: char_dum,char_dum2
      Call Dynamics_Variables_mem_allocation(medium_times,kpoints,Species) !Allocation of Dynamic Variables
      Call calc_Del_z_dimen_dummy()!Calculation of Del_z_dimen_dummy variable, a dummy variable of integration, dependant of wave vector k and dimensionality
      Call Calc_m_times(dtau,medium_times)
      Call Short_times_dynamics()
      Call Medium_times_dynamics(0,short_times+1)
      write(char_dum2,"(I2.1)") i_run
      !Call Dynamic_writting_tau_dependant ("Dynamics.dat",i_k_test,1,.True.,.False.)
      Call Dynamic_writting_Diff_msd ("Diff_MSD"//char_dum2 //".dat",1,.True.,.False.)
      !Call Dynamic_writting_Fc ("Fc.dat",i_k_test,1,.True.,.False.)
      !Call Dynamic_writting_Fs ("Fs.dat",i_k_test,1,.True.,.False.)
      Do i2=1, number_i_k_test
        Write(char_dum,"(I2.1)") i2
        Call Dynamic_writting_tau_dependant_multi_k_tests("Fc_"//char_dum2//"_"//&
        &char_dum // ".dat","Fs_"//char_dum2//"_"// char_dum // ".dat",i_k_test(i2),1,.True.,.False.,300+i2,400+i2)
        print*, i_k_test(i2)
      End Do
      Call Dynamic_writting_Delz ("Delz_"//char_dum2 //".dat",1,.True.,.False.)
      Do i1=1, decimations
        Call half_times_save()
        Print*, "decimation:",i1
        Call Medium_times_dynamics(1,(medium_times/2)+1)
        If (i1<decimations) Then
          !Call Dynamic_writting_tau_dependant ("Dynamics.dat",i_k_test,medium_times/2+1,.False.,.False.)
          Call Dynamic_writting_Diff_msd ("Diff_MSD"//char_dum2 //".dat",medium_times/2+1,.False.,.False.)
          !Call Dynamic_writting_Fc ("Fc.dat",i_k_test,medium_times/2+1,.False.,.False.)
          !Call Dynamic_writting_Fs ("Fs.dat",i_k_test,medium_times/2+1,.False.,.False.)
          Do i2=1, number_i_k_test
            Write(char_dum,"(I2.1)") i2
            Call Dynamic_writting_tau_dependant_multi_k_tests("Fc_"//char_dum2//"_"//&
            &char_dum // ".dat","Fs_"//char_dum2//"_"// char_dum // ".dat",i_k_test(i2),&
            &medium_times/2+1,.False.,.False.,300+i2,400+i2)
          End Do
          Call Dynamic_writting_Delz ("Delz_"//char_dum2 //".dat",medium_times/2+1,.False.,.False.)
        Else if (i1==decimations) Then
          !Call Dynamic_writting_tau_dependant ("Dynamics.dat",i_k_test,medium_times/2+1,.False.,.True.)
          Call Dynamic_writting_Diff_msd ("Diff_MSD"//char_dum2 //".dat",medium_times/2+1,.False.,.True.)
          !Call Dynamic_writting_Fc ("Fc.dat",i_k_test,medium_times/2+1,.False.,.True.)
          !Call Dynamic_writting_Fs ("Fs.dat",i_k_test,medium_times/2+1,.False.,.True.)
          Do i2=1, number_i_k_test
            Write(char_dum,"(I2.1)") i2
            Call Dynamic_writting_tau_dependant_multi_k_tests("Fc_"//char_dum2//"_"//&
            &char_dum // ".dat","Fs_"//char_dum2//"_"// char_dum // ".dat",&
            &i_k_test(i2),medium_times/2+1,.False.,.True.,300+i2,400+i2)
          End Do
          Call Dynamic_writting_Delz ("Delz_"//char_dum2 //".dat",medium_times/2+1,.False.,.True.)
        End If
      End Do
      Do i1=1,Species
        Dl(i1)=1.d0/(1.d0+Dl(i1))
        Dl(i1)=Dl(i1)*D0M(i1,i1)
        Print*, "Dl",i1,"=",Dl(i1)
      End Do
      DEAllocate(Del_z_dimen_dummy)
    End Subroutine

End Module
