Module NESCGLE
  Implicit None

  Contains
!!!!
  Subroutine calc_sku_instant_isochoric_quench_mono(ski,skf,u)
    Use Variables
    Implicit None
    Integer :: i1
    Real * 8, Dimension(:,:,:) ::ski,skf
    Real * 8 :: u
    Do i1=1,kpoints
      sk(i1,1,1)=skf(i1,1,1)+((ski(i1,1,1)-skf(i1,1,1))*exp(-2.d0*k(i1)*k(i1)*D0(1)*u/skf(i1,1,1)))
    End Do
  End Subroutine
!!!!
  Subroutine calc_ua_instant_isochoric_quench_mono(ski,skf,ua,gamma_value)
    Use NESCGLE_variables
    Use Gamma
    Implicit None
    Real * 8 :: u_change,initial_u,try_u,final_u
    Real * 8 :: error
    Real * 8, Parameter :: Tolerance=1.d-4
    Real * 8, Dimension(:,:,:), Intent(in) :: ski,skf
    Real * 8, Intent (out) :: ua,gamma_value
    u_change=du
    initial_u=0
    final_u=du
    Call calc_sku_instant_isochoric_quench_mono(ski,skf,final_u)
    gamma_value=Gamma_mono()
    Do While (gamma_value>gamma_max.And.final_u<gamma_max)
      u_change=2.d0*u_change
      initial_u=final_u
      final_u=final_u+u_change
      Call calc_sku_instant_isochoric_quench_mono(ski,skf,final_u)
      gamma_value=Gamma_mono()
    End Do
    error=1.d0
    Do While (error>Tolerance)
      u_change=abs(final_u-initial_u)/2.d0
      try_u=initial_u+u_change
      Call calc_sku_instant_isochoric_quench_mono(ski,skf,try_u)
      gamma_value=Gamma_mono()
      If(gamma_value<gamma_max)Then
        final_u=try_u
      Else
        initial_u=try_u
      End If
      error=abs(final_u-initial_u)/final_u
      Print*, "u_0=",initial_u, "u_f=",final_u
    End Do
    ua=final_u
    Call calc_sku_instant_isochoric_quench_mono(ski,skf,ua)
    gamma_value=Gamma_mono()
    Print*, "ua=",ua,"gamma=",gamma_value
  End Subroutine
!!!!
  Subroutine calc_quench_ski_skf_parameters_sw_mono (writting_option,vw_option,initial_Temp,final_Temp,ski,skf,SW_param)
    Use Variables
    Use Structure
    Use Writting
    Implicit None
    Real * 8, Intent(in) :: initial_Temp,final_Temp,SW_param
    Real * 8, Dimension(:,:,:),Intent(out) :: ski,skf
    Logical, Intent(in) :: vw_option,writting_option
    Integer :: i1,i2,i3
    Call calc_ck_sharma_sharma_HS_SW(vw_option,initial_Temp,SW_param)
    Do i1=1, kpoints
      ski(i1,1,1)=sk(i1,1,1)
    End Do
    If (writting_option.Eqv..True.)Then
      Call  Sk_writting ("ski.dat")
    End If
    Call calc_ck_sharma_sharma_HS_SW(vw_option,final_Temp,SW_param)
    Do i1=1, kpoints
      skf(i1,1,1)=sk(i1,1,1)
    End Do
    If (writting_option.Eqv..True.)Then
      Call  Sk_writting ("skf.dat")
    End If
  End Subroutine
!!!!
  Subroutine calc_gamma_vs_Tf_at_ua_mono_sh_sh_SW(vw_option,eta_value,Tempi,minor_Tempf,major_Tempf,Temp_steps,SW_param,File_name)
    Use Variables
    Implicit None
    Integer :: i1,i2
    Integer, Intent(in) :: Temp_steps
    Real * 8 :: dTemp,Temp
    Real * 8 :: gam_value,ua
    Real * 8, Intent (in) :: eta_value,Tempi,minor_Tempf,major_Tempf,SW_param
    Real * 8, Dimension(:,:,:), Allocatable :: ski,skf
    Logical, Intent(in) :: vw_option
    Character(len=*), Intent (in) :: File_name
    Allocate(ski(kpoints,1,1))
    Allocate(skf(kpoints,1,1))
    dTemp=abs(major_Tempf-minor_Tempf)/(1.d0*Temp_steps)
    Open (unit=12, File=File_name,Status="Replace")
    Do i1=0,Temp_steps
      Temp=minor_Tempf+dTemp*i1
      Call calc_quench_ski_skf_parameters_sw_mono (.FALSE.,vw_option,Tempi,Temp,ski,skf,SW_param)
      Call calc_ua_instant_isochoric_quench_mono(ski,skf,ua,gam_value)
      Write (12,*) Temp,gam_value,ua
      Flush(12)
    End Do
    Close(12)
  End Subroutine

End Module
