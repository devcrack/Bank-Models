Module NESCGLE
  Use Parameter_Variables
  Use System_Variables
  Use Dimensional_Variables
  Use Static_Variables
  Use SCGLE_Variables
  Use Structure_module
  Use NESCGLE_variables
  Use Matrix_operations
  Use Integration_methods
  Use Dynamics_Variables
  Use Writting
  Use SCGLE

  Implicit None

Contains

  Subroutine Calc_sku_quench_mono(First_time,Last_time,T0,Tf,u_value)
    Implicit None
    Logical, Intent (in) :: First_time,Last_time
    Real * 8, Intent(in):: T0,Tf,u_value
    Integer :: i1
    !Calculation of the initial and final structure factor
    If (First_time.Eqv..True.)Then
      Call NESCGLE_Static_Variables_mem_alloc(kpoints,Species)
      !Selection of the system
      Call sharma_sharma_attractive_yukawa_2D(T0)
      Call Sk_writting("sk_0.dat")
      Sk0=Sk
      Call sharma_sharma_attractive_yukawa_2D(Tf)
      Call Sk_writting("sk_f.dat")
      Skf=Sk
    End If
    !NE-SCGLE Interpolating function
    Do i1=1,kpoints
      Sk(i1,1,1)=Skf(i1,1,1)+((Sk0(i1,1,1)-Skf(i1,1,1))*exp(-2.d0*(k(i1)**2.d0)*u_value/Skf(i1,1,1)))
    End Do
    If (Last_time.Eqv..True.)Then
      Call NESCGLE_Static_Variables_mem_dealloc()
    End If
  End Subroutine


  !Subroutine that finds out if there is such a thing as a finite value of u(t) in which the system becoms arrested
  Subroutine Calc_ua_mono(ua,T0,Tf,gamma_ua)
    Implicit None
    Real * 8, Intent(out) :: ua
    Real * 8, Intent(in) :: T0, Tf
    Real * 8 :: error, dua, ua_max,ua_min
    Real * 8, Parameter :: Tolerance=1.d-6
    Real * 8, Dimension(1,1) :: gamma_value
    Real * 8, Dimension(1,1), Intent(out) :: gamma_ua
    Logical :: arrested_condition
    ua=1.d0
    Call  Calc_sku_quench_mono(.True.,.False.,T0,Tf,ua)
    Print*, "Holis1"
    Call rho_ini() !Calculation of m_rho and m_rhoi matrixes
    Print*, "Holis2"
    Call calc_ck_hk_ski()!Calculation for the direct correlation function and total correlation function
    Print*, "Holis3"
    Call Gamma_mixture(gamma_value)  !Call of gamma subroutine
    arrested_condition=.True.
    print*, "Holis"
    If (gamma_value(1,1)<gamma_max) Then
      arrested_condition=.True.
      ua_max=ua
      Do While (gamma_value(1,1)<gamma_max)
        ua=ua/2.d0
        Call  Calc_sku_quench_mono(.False.,.False.,T0,Tf,ua)
        Call rho_ini() !Calculation of m_rho and m_rhoi matrixes
        Call calc_ck_hk_ski()!Calculation for the direct correlation function and total correlation function
        Call Gamma_mixture(gamma_value)  !Call of gamma subroutine
        If (gamma_value(1,1)<gamma_max) Then
          ua_max=ua
          gamma_ua=gamma_value
        End If
      End Do
      ua_min=ua
    Else
      ua_min=ua
      Do While (gamma_value(1,1)>gamma_max.and.ua<gamma_max)
        ua=ua*2.d0
        Call  Calc_sku_quench_mono(.False.,.False.,T0,Tf,ua)
        Call rho_ini() !Calculation of m_rho and m_rhoi matrixes
        Call calc_ck_hk_ski()!Calculation for the direct correlation function and total correlation function
        Call Gamma_mixture(gamma_value)  !Call of gamma subroutine
        If (gamma_value(1,1)>gamma_max) Then
          ua_min=ua
        End If
        If (ua>gamma_max) Then
          arrested_condition=.False.
        End If
      End Do
      ua_max=ua
      gamma_ua=gamma_value
    End If
    If (arrested_condition.Eqv..True.) Then
      error=Abs(ua_max-ua_min)/ua_max
      Do while (error>Tolerance)
        ua=(ua_min+ua_max)/2.d0
        Call  Calc_sku_quench_mono(.False.,.False.,T0,Tf,ua)
        Call rho_ini() !Calculation of m_rho and m_rhoi matrixes
        Call calc_ck_hk_ski()!Calculation for the direct correlation function and total correlation function
        Call Gamma_mixture(gamma_value)  !Call of gamma subroutine
        If (gamma_value(1,1)>gamma_max) Then
          ua_min=ua
        Else If(gamma_value(1,1)<gamma_max) Then
          ua_max=ua
          gamma_ua=gamma_value
        End If
        error=Abs(ua_max-ua_min)/ua_max
      End Do
      ua=ua_max
      Call NESCGLE_Static_Variables_mem_dealloc()
    End If
  End Subroutine

  Subroutine Calc_isochoric_ua_Tf_diagram_mono()
    Implicit None
    Real * 8:: T0,Tf,dTf,ua
    Real * 8, Dimension(1:1) :: gamma_ua
    Integer :: i1
    Open(unit=99, File="isochoric_arrest_eta_0,40.dat", Status="Replace")
    Call z_yuk_mono(2.d0,1)
    kpoints=2**12 !Wave vector points
    T0=1.d0
    ua=1.d0
    dTf=1.d-3
    Tf=0.001d0
    dk=1.0d-2    !Parameter System variables assignement and allocation
    SDimen=2  !Space dimension
    Species=1 !Number of species
    Call Sys_Variables_Alloc(species) !Memory Allocation for System Variables
    sigma(1)=distu  !System diameter
    Call D0M_ini(species,sigma) !Diffusion coefficient of the system particle
    Call Static_Variables_mem_alloc(kpoints,species)
    Call Calc_static_k(dk,kpoints)
    Call kc_ini(species,1.305d0)
    Call Calc_lambdak(Species,kpoints)  !Calculation of lambda(k), SCGLE parameter function
    eta(1)=0.40d0
    Call sharma_sharma_attractive_yukawa_2D(T0)
    !Call Sk_writting("sk_test.dat")
    Do while (ua<gamma_max)
      Call Calc_ua_mono(ua,T0,Tf,gamma_ua)
      print*, "holis"
      Write(99,*) Tf,gamma_ua,ua
      Call Sk_writting("sk.dat")
      flush(99)
      Tf=Tf+dTf
    End Do
    Close(99)

  End Subroutine

End Module
