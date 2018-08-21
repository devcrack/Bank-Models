Module Structure_module
  Use Parameter_Variables
  Use System_Variables
  Use Static_Variables
  Use Matrix_operations
  Use Potential_Variables
  Implicit None

Contains
  !Subroutine for the calculation of a static number of equally spaced wave vectors
  Subroutine  Calc_static_k(difk,knumber)
    Implicit None
    Real * 8, Intent (in) :: difk
    Integer, Intent (in) :: knumber
    Integer :: i1
    Do i1=1, knumber
       k(i1)=i1*difk            !Thi vector is defined in Module Static Variables
       print *, 'k(', i1,') = ',k(i1)
                                ! in the file variables_modules.f90
    End Do
  End Subroutine Calc_static_k
  !Subroutine for the calculation of the Structure factor of a Hard Sphere System
  !Using the Percus-Yevick Closure Relationship
  !vw_option Logical value, True for Verlet-Weiss Correction, False Otherwise
  Subroutine Calc_Sk_hs_py_mono(vw_option)
    Implicit None
    Integer :: i1
    Logical :: vw_option
    Real * 8, Dimension(:), Allocatable :: k_vw
    Real * 8 :: eta_vw,dum1,dum2,dum3,dum4,dumsin,dumcos
    !Start of parameter checks
    If (SDimen/= 3 .or. Species /=1 ) Then
       print*, ('Error, cannot compute the system direct correlation function')
       Stop
    End If
    !End of parameter checks
    Allocate(k_vw(kpoints))
    !Start of VW Correction
    If (vw_option .EQV. .TRUE.) Then
       eta_vw=eta(1)*(1.d0-(eta(1)/16.d0))
       ! print *, 'What the fuck contains k(:)', k(:)*((eta_vw/eta(1))**(1.d0/3.d0)) 
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
       print *, 'SK(', i1, ') = ',Sk(i1,1,1)
       Ski(i1,1,1)=1.d0-Ck(i1,1,1)
    End Do
  End Subroutine Calc_Sk_hs_py_mono

  Subroutine calc_ck_hk_ski()
    Integer :: i1,i2
    Real * 8, Dimension(:,:), Allocatable ::  m_ck_dum,m_hk_dum
    Real * 8, Dimension(:,:), Allocatable :: m_dumS,m_dum1,m_dum2,m_dum3,I_m,m_dumSi
    Allocate(m_ck_dum(Species,Species))
    Allocate(m_hk_dum(Species,Species))
    Allocate(m_dumS(Species,Species))
    Allocate(m_dum1(Species,Species))
    Allocate(m_dum2(Species,Species))
    Allocate(m_dum3(Species,Species))
    Allocate(I_m(Species,Species))
    Allocate(m_dumSi(Species,Species))
    Call identity_matrix(Species,I_m)
    Do i1=1,kpoints
       Call matrix_conversion_3_to_2(i1,Species,Sk,m_dumS)
       Call inverse_matrix(Species,m_dumS,m_dumSi)
       Call matrix_conversion_2_to_3(i1,Species,m_dumSi,Ski)
       m_dum1=m_dumS-I_m
       Call matrix_multiplication(Species,m_rhoi,m_dum1,m_dum2)
       Call matrix_multiplication(Species,m_dum2,m_rhoi,m_hk_dum)
       Call matrix_conversion_2_to_3(i1,Species,m_hk_dum,hk)
       !Call matrix_conversion_2_to_3(i1,Species,m_dum1,hk)
       m_dum1=I_m-m_dumSi
       Call matrix_multiplication(Species,m_rhoi,m_dum1,m_dum2)
       Call matrix_multiplication(Species,m_dum2,m_rhoi,m_ck_dum)
       Call matrix_conversion_2_to_3(i1,Species,m_ck_dum,ck)
       !Call matrix_conversion_2_to_3(i1,Species,m_dum1,ck)
       !Print*, ck(i1,1,1), hk(i1,1,1), k(i1)
    End Do
    DEAllocate(m_ck_dum)
    DEAllocate(m_hk_dum)
    DEAllocate(m_dumS)
    DEAllocate(m_dum1)
    DEAllocate(m_dum2)
    DEAllocate(m_dum3)
    DEAllocate(I_m)
    DEAllocate(m_dumSi)
  End Subroutine calc_ck_hk_ski

  !Santos SW S(k)
  Subroutine calc_sk_hs_sw_mono_santos (Temp)
    Implicit None
    Integer :: i1
    Real *8, Intent(in) :: Temp
    Real *8 :: A0,A01,L1,L2,S1,S2,S3,sw_tau !F(t) Parameters
    Real *8 :: alpha1,alpha2,alpha3,beta1,beta2,beta3,beta4,gamma0,gamma1,gamma2 !L2 Parameters
    Real *8 :: Lam1,Lam2,Lam3,Lam4,Lam5,Lam6,tau1,tau2 !convenient dummies variables
    Complex *16, Dimension (:), Allocatable :: Ftp,Ftm,Gtp,Gtm !F(t),F(-t),G(t),G(-t)
    Real *8 :: Ftdum1,Ftdum2,Ftdum3
    real *8 dum1,dum2,dum3
    Complex*16 :: im
    Real *8 dumeta
    If (SDimen/= 3 .or. Species /=1 ) Then
       print*, ('Error, cannot compute the system direct correlation function')
       Stop
    End If
    dumeta=eta(1)
    im=(0.d0,1.d0)
    Allocate (sw_lambda(1))
    Allocate (Ftp(kpoints))
    Allocate (Ftm(kpoints))
    Allocate (Gtp(kpoints))
    Allocate (Gtm(kpoints))
    sw_lambda(1)=1.1
    Lam1=sw_lambda(1)
    Lam2=Lam1**2
    Lam3=Lam1**3
    Lam4=Lam1**4
    Lam5=Lam1**5
    Lam6=Lam1**6
    !A0
    A0=exp(1.d0/Temp)-1.d0
    !A0 prima
    A01=A0*((Lam1-1.d0)**2)
    print*, "A01=",A01
    !square well tau
    sw_tau=1.d0/(12.d0*A0*Lam1*(Lam1-1.d0))
    tau1=1.d0/sw_tau
    print*, "tau1=",tau1
    tau2=1.d0/(sw_tau**2)
    print*, "tau2=",tau2
    !alpha 1
    alpha1=2.d0*A01*(2.d0+Lam1)+((1.d0+(4.d0*Lam1)+Lam2-3.d0*A01)*tau1/6.d0)
    print*, "alpha1", alpha1
    !alpha2
    alpha2=3.d0+A01*(7.d0+2.d0*Lam1-3.d0*Lam2)-((7.d0+Lam1+16.d0*Lam2+A01*&
         &(23.d0+15.d0*Lam1+15.d0*Lam2+7.d0*Lam3))*tau1/12.d0)
    print*, "alpha2", alpha2
    !alpha3
    alpha3=-2.d0-2.d0*A01*(1.d0+2.d0*Lam1+3.d0*Lam2)+((7.d0+Lam1-2.d0*Lam2+A01*&
         &(7.d0+15.d0*Lam1+21.d0*Lam2+11.d0*Lam3+6.d0*Lam4))*tau1/6.d0)
    print*, "alpha3", alpha3
    !beta1
    beta1=-4.d0-4.d0*A01*(2.d0+Lam1)-(tau2/3.d0)+((5.d0+2.d0*Lam1-Lam2+3.d0*A01)*tau1/3.d0)
    print*, "beta1", beta1
    !beta2
    beta2=6.d0+6.d0*A01*(3.d0+2.d0*Lam1+Lam2)+4.d0*(A01**2)*((2.d0+Lam1)**2)&
         &-((9.d0*(3.d0+Lam1)+A01*(59.d0+51.d0*Lam1+3.d0*Lam2-5.d0*Lam3)&
         &+12.d0*(A01**2)*(2.d0+Lam1))*tau1/6.d0)+((31.d0+8.d0*Lam1+18.d0*Lam2&
         &-4.d0*Lam3+Lam4+12.d0*A01*(5.d0+Lam1)+9.d0*(A01**2))*tau2/36.d0)
    print*, "beta2", beta2
    !beta3
    beta3=-4.d0-12.d0*A01*(1.d0+Lam1+Lam2)-4.d0*(A01**2)*(2.d0+Lam1)*(1.d0&
         &+2.d0*Lam1+3.d0*Lam2)+((3.d0*(4.d0+Lam1+Lam2)+A01*(29.d0&
         &+36.d0*Lam1+30.d0*Lam2+7.d0*Lam3-3.d0*Lam4)+(A01**2)&
         &*(17.d0+43.d0*Lam1+30.d0*Lam2+Lam3-Lam4))*tau1/3.d0)-((29.d0+13.d0*Lam1&
         &+27.d0*Lam2+7.d0*Lam3-4.d0*Lam4+A01*(68.d0+80.d0*Lam1+80.d0*Lam2+16.d0*Lam3&
         &+7.d0*Lam4+Lam5)+3.d0*(A01**2)*(13.d0+13.d0*Lam1+Lam2-3.d0*Lam3))*tau2/36.d0)
    print*, "beta3", beta3
    !beta4 BE WARE, MODIFIED FROM THE ARTICLE
    dum1=(1.d0+A01*(1.d0+2.d0*Lam1+3.d0*Lam2))**2
    dum2=2.d0*A01*(7.d0+15.d0*Lam1+15.d0*Lam2+5.d0*Lam3+6.d0*Lam4)
    dum3= A01*A01*(1.d0+2.d0*Lam1+3.d0*Lam2)*(7.d0+15.d0*Lam1+3.d0*Lam2-Lam3)
    dum2=7.d0+Lam1+4.d0*Lam2+dum2+dum3
    beta4=dum1-(dum2*tau1/6.d0)
    dum1=49.d0+38.d0*Lam1+9.d0*Lam2+8.d0*Lam3+16.d0*Lam4
    dum2=2.d0*A01*(49.d0+88.d0*Lam1+112.d0*Lam2+80.d0*Lam3+35.d0*Lam4-4.d0*Lam5)
    dum3=A01*A01*(49.d0+138.d0*Lam1+219.d0*Lam2+124.d0*Lam3+27.d0*Lam4+18.d0*Lam5+Lam6)
    dum1=dum1+dum2+dum3
    dum1=dum1*tau2/144.d0
    beta4=beta4+dum1
    print*, "beta4",beta4
    !gamma0
    gamma0=1.d0+Lam1-(tau1/6.d0)
    print*, "gamma0",gamma0
    !gamma1
    gamma1=2.d0*(1.d0+Lam1-Lam2)-((3.d0+Lam3)*tau1/6.d0)
    print*, "gamma1",gamma1
    !gamma2
    gamma2=Lam1*(-4.d0*Lam1+((2.d0+3.d0*Lam1+Lam2+Lam3)*tau1/3.d0))
    print*, "gamma2",gamma2
    !L2
    L2=-1.d0+(alpha1*dumeta)+(alpha2*dumeta**2)+(alpha3*dumeta**3)+((1.d0+2.d0*dumeta)*&
         &dsqrt(1.d0+(beta1*dumeta)+(beta2*dumeta**2)+(beta3*dumeta**3)+(beta4*dumeta**4)))
    L2=L2/(12.d0*dumeta*(gamma0+gamma1*dumeta+gamma2*dumeta**2))
    print*, "L2",L2
    !L1
    L1=1.d0+0.5d0*dumeta+2.d0*dumeta*(1.d0+Lam1+Lam2)*L2-0.5d0*dumeta*A01*(3.d0+2.d0*Lam1+Lam2)
    L1=L1/(1.d0+2.d0*dumeta)
    print*, "L1",L1
    !S1
    S1=(-1.5d0+2.d0*(1.d0+Lam1+Lam2)*L2-0.5d0*(3.d0+2.d0*Lam1+Lam2)*A01)!Care inconsistency
    S1=S1*(dumeta/(1.d0+2.d0*dumeta))
    print*, "S1",S1
    !S2
    S2=(-1.d0+dumeta+2.d0*(1.d0-2.d0*dumeta*Lam1*(1.d0+Lam1))*L2-(1.d0-dumeta*((1.d0+Lam1)**2))*A01)
    S2=S2/(2.d0*(1.d0+2.d0*dumeta))
    print*, "S2",S2
    !S3
    S3=(-(((1.d0-dumeta)**2)/(12.d0*dumeta))-(0.5d0*(Lam1+1.d0)-dumeta*Lam2)*L2)
    S3=S3+(4.d0+2.d0*Lam1-dumeta*(3.d0*Lam2+2.d0*Lam1+1.d0))*(A01/12.d0)
    S3=S3/(1.d0+2.d0*dumeta)
    print*, "S3",S3
    !F(t)->F(ik)
    Ftdum1=L1+(L2/(Lam1-1.d0))-A0*(Lam1-1.d0)
    Ftdum2=L2/(Lam1-1.d0)
    Ftdum3=Lam1-1.d0
    Do i1=1, kpoints
       !Ftp(i1)=CMPLX(1.d0+A0,0)+CMPLX(0.D0,Ftdum1*k(i1))-(CMPLX(A0,Ftdum2*k(i1))*exp(CMPLX(0.d0,-Ftdum3*k(i1))))
       !Ftp(i1)=-Ftp(i1)/CMPLX((1.d0-S2*(k(i1)**2)),(S1*k(i1)-S3*(k(i1)**3)))/(12.d0*eta)
       !Ftm(i1)=CMPLX(1.d0+A0,0)+CMPLX(0.D0,-Ftdum1*k(i1))-(CMPLX(A0,-Ftdum2*k(i1))*exp(CMPLX(0.d0,Ftdum3*k(i1))))
       !Ftm(i1)=-Ftp(i1)/CMPLX(12.d0*eta*(1.d0-S2*(k(i1)**2)),-12.d0*eta*(S1*k(i1)-S3*(k(i1)**3)))
       Ftp(i1)=(1.d0+A0+Ftdum1*(im*k(i1)))-(A0+Ftdum2*(im*k(i1)))*exp(-Ftdum3*(im*k(i1)))
       Ftp(i1)=Ftp(i1)/(1.d0+(S1*(im*k(i1)))+(S2*(im*k(i1))**2)+(S3*(im*k(i1))**3))
       Ftp(i1)=-Ftp(i1)/(12.d0*dumeta)
       Ftm(i1)=(1.d0+A0+Ftdum1*(-im*k(i1)))-(A0+Ftdum2*(-im*k(i1)))*exp(-Ftdum3*(-im*k(i1)))
       Ftm(i1)=Ftm(i1)/(1.d0+(S1*(-im*k(i1)))+(S2*(-im*k(i1))**2)+(S3*(-im*k(i1))**3))
       Ftm(i1)=-Ftm(i1)/(12.d0*dumeta)
       !G(t)->G(ik)
       !Gtp(i1)=CMPLX(0.d0,k(i1))*Ftp(i1)*exp(CMPLX(0.d0,-k(i1)))/(CMPLX(1.d0,0.d0)+12.d0*eta*Ftp(i1)*exp(CMPLX(0.d0,-k(i1))))
       Gtp(i1)=(im*k(i1))*Ftp(i1)*exp(-im*k(i1))/(1.d0+12.d0*dumeta*Ftp(i1)*exp(-im*k(i1)))
       !G(-t)->G(-ik)
       !Gtm(i1)=CMPLX(0.d0,-k(i1))*Ftm(i1)*exp(CMPLX(0.d0,k(i1)))/(CMPLX(1.d0,0.d0)+12.d0*eta*Ftm(i1)*exp(CMPLX(0.d0,k(i1))))
       Gtm(i1)=(-im*k(i1))*Ftm(i1)*exp(im*k(i1))/(1.d0+12.d0*dumeta*Ftm(i1)*exp(im*k(i1)))

       sk(i1,1,1)=Realpart((Gtp(i1)-Gtm(i1))/(im*k(i1)))
       sk(i1,1,1)=1.d0-(12.d0*dumeta*sk(i1,1,1))
       !sk(i1,1)=Realpart(Gtm(i1))
    End Do
  End Subroutine calc_sk_hs_sw_mono_santos

  Subroutine Calc_S0_HS_SW_santos(Temp,compressibility)
    Implicit None
    Integer :: i1
    Real *8, Intent(in) :: Temp
    Real *8 :: A0,A01,L1,L2,S1,S2,S3,sw_tau !F(t) Parameters
    Real *8 :: alpha1,alpha2,alpha3,beta1,beta2,beta3,beta4,gamma0,gamma1,gamma2 !L2 Parameters
    Real *8 :: Lam1,Lam2,Lam3,Lam4,Lam5,Lam6,tau1,tau2 !convenient dummies variables
    Real *8 dum1,dum2,dum3
    Real *8, Intent(out) :: compressibility
    Real *8 dumeta,eta1,eta2,eta3,eta4
    If (SDimen/= 3 .or. Species /=1 ) Then
       print*, ('Error, cannot compute the system direct correlation function')
       Stop
    End If
    dumeta=eta(1)
    sw_lambda(1)=1.1
    Lam1=sw_lambda(1)
    Lam2=Lam1**2
    Lam3=Lam1**3
    Lam4=Lam1**4
    Lam5=Lam1**5
    Lam6=Lam1**6
    eta1=eta(1)
    eta2=eta1*eta1
    eta3=eta2*eta1
    eta4=eta3*eta1
    !A0
    A0=exp(1.d0/Temp)-1.d0
    !A0 prima
    A01=A0*((Lam1-1.d0)**2)
    print*, "A01=",A01
    !square well tau
    sw_tau=1.d0/(12.d0*A0*Lam1*(Lam1-1.d0))
    tau1=1.d0/sw_tau
    print*, "tau1=",tau1
    tau2=1.d0/(sw_tau**2)
    print*, "tau2=",tau2
    !alpha 1
    alpha1=2.d0*A01*(2.d0+Lam1)+((1.d0+(4.d0*Lam1)+Lam2-3.d0*A01)*tau1/6.d0)
    print*, "alpha1", alpha1
    !alpha2
    alpha2=3.d0+A01*(7.d0+2.d0*Lam1-3.d0*Lam2)-((7.d0+Lam1+16.d0*Lam2+A01*&
         &(23.d0+15.d0*Lam1+15.d0*Lam2+7.d0*Lam3))*tau1/12.d0)
    print*, "alpha2", alpha2
    !alpha3
    alpha3=-2.d0-2.d0*A01*(1.d0+2.d0*Lam1+3.d0*Lam2)+((7.d0+Lam1-2.d0*Lam2+A01*&
         &(7.d0+15.d0*Lam1+21.d0*Lam2+11.d0*Lam3+6.d0*Lam4))*tau1/6.d0)
    print*, "alpha3", alpha3
    !beta1
    beta1=-4.d0-4.d0*A01*(2.d0+Lam1)-(tau2/3.d0)+((5.d0+2.d0*Lam1-Lam2+3.d0*A01)*tau1/3.d0)
    print*, "beta1", beta1
    !beta2
    beta2=6.d0+6.d0*A01*(3.d0+2.d0*Lam1+Lam2)+4.d0*(A01**2)*((2.d0+Lam1)**2)&
         &-((9.d0*(3.d0+Lam1)+A01*(59.d0+51.d0*Lam1+3.d0*Lam2-5.d0*Lam3)&
         &+12.d0*(A01**2)*(2.d0+Lam1))*tau1/6.d0)+((31.d0+8.d0*Lam1+18.d0*Lam2&
         &-4.d0*Lam3+Lam4+12.d0*A01*(5.d0+Lam1)+9.d0*(A01**2))*tau2/36.d0)
    print*, "beta2", beta2
    !beta3
    beta3=-4.d0-12.d0*A01*(1.d0+Lam1+Lam2)-4.d0*(A01**2)*(2.d0+Lam1)*(1.d0&
         &+2.d0*Lam1+3.d0*Lam2)+((3.d0*(4.d0+Lam1+Lam2)+A01*(29.d0&
         &+36.d0*Lam1+30.d0*Lam2+7.d0*Lam3-3.d0*Lam4)+(A01**2)&
         &*(17.d0+43.d0*Lam1+30.d0*Lam2+Lam3-Lam4))*tau1/3.d0)-((29.d0+13.d0*Lam1&
         &+27.d0*Lam2+7.d0*Lam3-4.d0*Lam4+A01*(68.d0+80.d0*Lam1+80.d0*Lam2+16.d0*Lam3&
         &+7.d0*Lam4+Lam5)+3.d0*(A01**2)*(13.d0+13.d0*Lam1+Lam2-3.d0*Lam3))*tau2/36.d0)
    print*, "beta3", beta3
    !beta4 BE WARE, MODIFIED FROM THE ARTICLE
    dum1=(1.d0+A01*(1.d0+2.d0*Lam1+3.d0*Lam2))**2
    dum2=2.d0*A01*(7.d0+15.d0*Lam1+15.d0*Lam2+5.d0*Lam3+6.d0*Lam4)
    dum3= A01*A01*(1.d0+2.d0*Lam1+3.d0*Lam2)*(7.d0+15.d0*Lam1+3.d0*Lam2-Lam3)
    dum2=7.d0+Lam1+4.d0*Lam2+dum2+dum3
    beta4=dum1-(dum2*tau1/6.d0)
    dum1=49.d0+38.d0*Lam1+9.d0*Lam2+8.d0*Lam3+16.d0*Lam4
    dum2=2.d0*A01*(49.d0+88.d0*Lam1+112.d0*Lam2+80.d0*Lam3+35.d0*Lam4-4.d0*Lam5)
    dum3=A01*A01*(49.d0+138.d0*Lam1+219.d0*Lam2+124.d0*Lam3+27.d0*Lam4+18.d0*Lam5+Lam6)
    dum1=dum1+dum2+dum3
    dum1=dum1*tau2/144.d0
    beta4=beta4+dum1
    print*, "beta4",beta4
    !gamma0
    gamma0=1.d0+Lam1-(tau1/6.d0)
    print*, "gamma0",gamma0
    !gamma1
    gamma1=2.d0*(1.d0+Lam1-Lam2)-((3.d0+Lam3)*tau1/6.d0)
    print*, "gamma1",gamma1
    !gamma2
    gamma2=Lam1*(-4.d0*Lam1+((2.d0+3.d0*Lam1+Lam2+Lam3)*tau1/3.d0))
    print*, "gamma2",gamma2
    !L2
    L2=-1.d0+(alpha1*dumeta)+(alpha2*dumeta**2)+(alpha3*dumeta**3)+((1.d0+2.d0*dumeta)*&
         &dsqrt(1.d0+(beta1*dumeta)+(beta2*dumeta**2)+(beta3*dumeta**3)+(beta4*dumeta**4)))
    L2=L2/(12.d0*dumeta*(gamma0+gamma1*dumeta+gamma2*dumeta**2))
    print*, "L2",L2
    compressibility=5.d0-20.d0*eta1*(1.d0+(2.d0+Lam1)*A01-3.d0*(1.d0+Lam1)*L2)+&
                                !eta**2 terms
         (2.d0*(eta2)*(15.d0-6.d0*(14.d0-Lam1+19.d0*Lam2-Lam3-Lam4)*L2+120.d0*&
         (1.d0+Lam1+Lam2)*(L2**2.d0)+(50.d0+16.d0*Lam1+27.d0*Lam2-2.d0*Lam3-Lam4)*A01-&
         30.d0*(5.d0+4.d0*Lam1+3.d0*Lam2)*(A01*L2)+15.d0*(3.d0+2.d0*Lam1+Lam2)*(A01**2.d0)))+&
                                !eta**3 terms
         (-2.d0*(eta3)*(10.d0+3.d0*(7.d0-53.d0*Lam1-13.d0*Lam2+7.d0*Lam3-8.d0*Lam4)*&
         L2-60.d0*(1.d0+Lam1)*(1.d0-4.d0*Lam1+Lam2)*(1.d0+Lam1+Lam2)*(L2**2.d0)+(19.d0+59.d0*Lam1&
         +9.d0*Lam2-Lam3+4.d0*Lam4)*A01+3.d0*(11.d0-67.d0*Lam1-114.d0*Lam2-66.d0*Lam3-13.d0*&
         Lam4+9.d0*Lam5)*(A01*L2)+3.d0*(3.d0+2.d0*Lam1+Lam2)*(1.d0+7.d0*Lam1+3.d0*Lam2-Lam3)*(A01**2.d0)))+&
                                !eta**4 terms
         ((eta4)*(5.d0-12.d0*(1.d0+Lam1+Lam2+11.d0*Lam3-4.d0*Lam4)*L2+240.d0*Lam3*(1.d0+Lam1+Lam2)&
         *(L2**2.d0)+2.d0*(7.d0+8.d0*Lam1+9.d0*Lam2+10.d0*Lam3-4.d0*Lam4)*A01-12.d0*(1.d0+3.d0*Lam1+6.d0*&
         Lam2+24.d0*Lam3+17.d0*Lam4+9.d0*Lam5)*(A01*L2)+3.d0*(3.d0+2.d0*Lam1+Lam2)*(1.d0+2.d0*Lam1+&
         3.d0*Lam2+4.d0*Lam3)*(A01**2.d0)))
    compressibility=compressibility/(5.d0*((1.d0+2.d0*eta(1))**2.d0))
  End Subroutine Calc_S0_HS_SW_santos

  !Hard Sphere Mixture, Percus-Yevick approximation, Baxter paper
  !Complement Function
  Complex * 16 Function Q_Integral_function(kdum,xdum,adum,bdum,cdum,rhodum)
    Real * 8, Intent(in) :: kdum,xdum,adum,bdum,cdum,rhodum
    Complex * 16 :: I0,I1,I2,Im,dummy1,dummy2
    Im=(0.d0,1.d0)
    dummy1=exp(Im*kdum*xdum)
    I0=-Im/kdum
    I1=(-Im*xdum/kdum)+(1.d0/(kdum**2.d0))
    I2=(-Im*xdum*xdum/kdum)+(2.d0*xdum/(kdum**2.d0))+(2.d0*Im/(kdum**3.d0))
    dummy2=((cdum*I0)+(bdum*I1)+((adum/2.d0)*I2))*dummy1
    dummy2=(2.d0*pi*rhodum*dummy2)
    Q_Integral_function=dummy2
  End Function Q_Integral_function

  Subroutine Calc_sk_HS_py_mixture()
    Implicit None
    Integer :: i1,i2,i3,i4
    Real * 8 :: xi1,xi2,xi3,dumk1,dumk2,dumk3,dumk1n,dumk2n,dumk3n,rhodum
    Real * 8, Dimension(:,:), Allocatable :: Rab, Sab, I_m,cij,m_dumSi,m_dumS,m_dum1
    Real * 8, Dimension(:), Allocatable :: k_vw,eta_vw,ai,bi
    Complex * 16, Dimension(:,:), Allocatable :: cm_dum1,cm_dum2,cm_dum3
    Complex *16 :: cdum1,cdum2,cdum3
    Complex * 16, Dimension(:,:,:), Allocatable :: Qk,QkT
    Complex * 16, Parameter :: Im=(0.d0,1.d0)
    Allocate(m_dum1(Species,Species))
    Allocate (Rab(Species,Species))
    Allocate (Sab(Species,Species))
    Allocate (I_m(Species,Species))
    Allocate(k_vw(kpoints))
    Allocate(eta_vw(Species))
    Allocate(Qk(kpoints,Species,Species))
    Allocate(QkT(kpoints,Species,Species))
    Allocate(cm_dum1(Species,Species))
    Allocate(cm_dum2(Species,Species))
    Allocate(cm_dum3(Species,Species))
    Allocate(ai(Species))
    Allocate(bi(Species))
    Allocate(cij(Species,Species))
    Allocate(m_dumSi(Species,Species))
    Allocate(m_dumS(Species,Species))
    Call identity_matrix(Species,I_m)
    !Calculation of Rab and Sab matrixes, and xi dummies
    xi1=0.d0
    xi2=0.d0
    xi3=0.d0
    Do i1=1,Species
       xi1=(pi/6.d0)*(m_rho(i1,i1)**2.d0)*sigma(i1)+xi1
       xi2=(pi/6.d0)*(m_rho(i1,i1)**2.d0)*(sigma(i1)**2.d0)+xi2
       xi3=(pi/6.d0)*(m_rho(i1,i1)**2.d0)*(sigma(i1)**3.d0)+xi3
       Do i2=1,Species
          Rab(i1,i2)=(sigma(i1)+sigma(i2))/2.d0
          Sab(i1,i2)=(sigma(i1)-sigma(i2))/2.d0
       End Do
    End Do
    !Calculation of ai, bi and cij
    Do i1=1,Species
       ai(i1)=(1.d0-xi3+(3.d0*sigma(i1)*xi2))/((1.d0-xi3)**2.d0)
       bi(i1)=-(3.d0*(sigma(i1)**2.d0)*xi2)/(2.d0*((1.d0-xi3)**2.d0))
       Do i2=1,Species
          cij(i1,i2)=-((ai(i1)/2.d0)*(Rab(i1,i2)**2.d0))-(bi(i1)*Rab(i1,i2))
       End Do
    End Do
    k_vw=k
    eta_vw=eta
    Do i3=1,kpoints
       dumk1=k_vw(i3)
       dumk2=k_vw(i3)**2.d0
       dumk3=k_vw(i3)**3.d0
       dumk1n=-dumk1
       dumk2n=dumk2
       dumk3n=-dumk3
       Do i1=1,Species
          Do i2=1,Species
             rhodum=m_rho(i1,i1)*m_rho(i2,i2)
             cdum1=Q_Integral_function(k_vw(i3),Rab(i1,i2),ai(i1),bi(i1),cij(i1,i2),rhodum)
             cdum2=Q_Integral_function(k_vw(i3),Sab(i1,i2),ai(i1),bi(i1),cij(i1,i2),rhodum)
             Qk(i3,i1,i2)=I_m(i1,i2)-(cdum1-cdum2)
             cdum1=Q_Integral_function(-k_vw(i3),Rab(i1,i2),ai(i1),bi(i1),cij(i1,i2),rhodum)
             cdum2=Q_Integral_function(-k_vw(i3),Sab(i1,i2),ai(i1),bi(i1),cij(i1,i2),rhodum)
             QkT(i3,i2,i1)=I_m(i1,i2)-(cdum1-cdum2)
          End Do
       End Do
       Do i1=1,Species
          Do i2=1,Species
             m_dumSi(i1,i2)=0.d0
             Do i4=1,Species
                m_dumSi(i1,i2)=m_dumSi(i1,i2)+(QkT(i3,i1,i4)*Qk(i3,i4,i2))
             End Do
          End Do
       End Do
       Call inverse_matrix(Species,m_dumSi,m_dumS)
       Call matrix_multiplication(Species,m_dumSi,m_dumS,m_dum1)
       Call matrix_conversion_2_to_3(i3,Species,m_dumS,Sk)
    End Do
    Call identity_matrix(Species,I_m)
    DEAllocate(m_dum1)
    DEAllocate(m_dumS)
    DEAllocate(m_dumSi)
    DEAllocate(cij)
    DEAllocate(ai)
    DEAllocate(bi)
    DEAllocate(cm_dum1)
    DEAllocate(cm_dum2)
    DEAllocate(cm_dum3)
    DEAllocate(Qk)
    DEAllocate(QkT)
    DEAllocate(k_vw)
    DEAllocate(eta_vw)
    DEAllocate(Rab)
    DEAllocate(Sab)
    DEAllocate(I_m)
  End Subroutine Calc_sk_HS_py_mixture

  !Roth Monocomponent HD S(k)
  Subroutine calc_ck_hd_mono_roth()
    Implicit None
    Integer :: i1
    Real * 8 :: dummy1,dummy2,dummy3,ck_dum,J0,J1
    Do i1=1, kpoints
       J0=BESSEL_JN(0, k(i1)*sigma(1)/2.d0)
       J1=BESSEL_JN(1, k(i1)*sigma(1)/2.d0)
       ck_dum=pi/(6.d0*((1.d0-eta(1))**3.d0)*(k(i1)**2.d0))
       dummy1=(-5.d0*((1.d0-eta(1))*(k(i1)*sigma(1)/2.d0)*J0)**2.d0)
       dummy2=(4.d0*((eta(1)-20.d0)*eta(1)+7.d0))+5.d0*(((1.d0-eta(1))*(k(i1)*sigma(1)/2.d0))**2.d0)
       dummy2=dummy2*(J1**2.d0)
       dummy3=4.d0*(eta(1)-13.d0)*(1.d0-eta(1))*(k(i1)*sigma(1)/2.d0)*J1*J0
       ck_dum=ck_dum*(dummy1+dummy2+dummy3);
       Sk(i1,1,1)=1.d0/(1.d0-(4.d0*eta(1)*ck_dum/(pi*sigma(1)*sigma(1))))
    End Do
  End Subroutine calc_ck_hd_mono_roth

  !Roth Binary HD S(k)
  Subroutine calc_ck_hd_binary_mixture_roth()
    Implicit None
    Real * 8 n1,n2,n0,ki
    Real * 8, Dimension(2) :: x_frac,rho
    Real * 8, Dimension(:,:,:), Allocatable :: Nij
    Real * 8, Dimension(:), Allocatable :: dummy_Dk
    Real * 8, Dimension(:,:), Allocatable :: dum_ck,Sk_inv_dum,Ident,Sk_dum
    Integer :: i1,i2,i3
    Real * 8 :: dum1,dum2,dum3,ri,rj
    If (SDimen/= 2 .or. Species /=2 ) Then
       print*, ('Error, cannot compute the system direct correlation function')
       Stop
    End If
    Allocate(Nij(kpoints,species,species))
    Allocate(dummy_Dk(kpoints))
    Allocate (Sk_inv_dum(species,species))
    Allocate (Sk_dum(species,species))
    Allocate (dum_ck(species,species))
    Allocate (Ident(species,species))
    Do i1=1, species
       rho(i1)=4.d0*eta(i1)/(pi*sigma(i1)**2)
       x_frac(i1)=eta(i1)/((sigma(i1)**2)*sum(eta/(sigma**2)))
    End Do
    n2=sum(eta)
    n1=pi*sum(sigma*rho) !4.d0*sum(eta/sigma)
    n0=sum(rho)
    !print*, n2,n1,n0
    !ck
    !n1=n1*dimen_dist
    !n0=n0*dimen_dist*dimen_dist
    Do i1=1,species
       ri=sigma(i1)/(2.d0)
       Do i2=1,species
          rj=sigma(i2)/(2.d0)
          Do i3=1,kpoints
             ki=k(i3)
             !CHECAR EXPRESION
             dum1=ri*BESSEL_JN(1, ki*ri)*(rj*BESSEL_JN(1, ki*rj)*(5.d0*(ki**2)*((1.d0-n2)**2)&
                  &-12.d0*(n1**2)-24.d0*pi*n0*(1.d0-n2))-(12.d0*ki*(1.d0-n2)&
                  &*(n1*rj-n2+1.d0)*BESSEL_JN(0, ki*rj)))
             !
             dum2=-ki*(1.d0-n2)*rj*(12.d0*BESSEL_JN(0, ki*ri)*&
                  &(ki*(1.d0-n2)*ri*BESSEL_JN(0, ki*rj)&
                  &+(n1*ri+1.d0-n2)*BESSEL_JN(1, ki*rj))-7.d0*ki*(1.d0-n2)*ri&
                  &*BESSEL_JN(2, ki*ri)*BESSEL_JN(2,ki*rj))
             !
             dum3=pi/(6.d0*((1.d0-n2)**3)*(ki**2))

             ck(i3,i1,i2)=dum3*(dum1+dum2)
             !ck(i3,i1,i2)=sqrt(x_frac(i1)*x_frac(i2))*ck(i3,i1,i2)
          End Do
       End Do
    End Do
    Call identity_matrix(species,Ident)
    Do i3=1,kpoints
       Do i1=1,species
          Do i2=1,species
             Sk_inv_dum(i1,i2)=Ident(i1,i2)-sum(rho)*sqrt(x_frac(i1)*x_frac(i2))*ck(i3,i1,i2)
             Ski(i3,i1,i2)=Sk_inv_dum(i1,i2)
          End Do
       End Do
       Call inverse_matrix(species,Sk_inv_dum,Sk_dum)
       Do i1=1,species
          Do i2=1,species
             sk(i3,i1,i2)=Sk_dum(i1,i2)
          End Do
       End Do
    End Do
  End Subroutine calc_ck_hd_binary_mixture_roth

  !Rosenfeld Structure Factor
  Subroutine calc_sk_hd_monocomponent_Rosenfeld()
    Implicit None
    Integer :: i1
    Real * 8:: A,B,G,J0,J1,J12,Z1,Z2,X
    If (SDimen/= 2 .or. Species /=1 ) Then
       print*, ('Error, cannot compute the system Structure Factor')
       Stop
    End If
    Z1=(1.d0+0.128d0*(eta(1)**2d0)+0.027d0*(eta(1)**3.d0)+0.06d0*(eta(1)**4.d0))/(1.d0-eta(1))**2.d0
    Z2=((2.d0*0.128d0*eta(1))+(3.d0*0.027d0*(eta(1)**2d0))+(4.d0*0.06d0*(eta(1)**3.d0)))*((1.d0-eta(1))**2.d0)
    Z2=(Z2/(1.d0-eta(1))**4.d0)+(2.d0*(Z1/(1.d0-eta(1))))
    !Z1=1.0d0/(1.d0-eta(1))**2.d0
    !Z2=2.0d0/(1.d0-eta(1))**3.d0
    X=Z1+eta(1)*Z2
    G=Sqrt(Z2/2.d0)
    A=(1.d0+((2.d0*eta(1)-1.d0)*X)+(2.d0*eta(1)*G))/eta(1)
    B=(((1.d0-eta(1))*X)-1.d0-(3.d0*eta(1)*G))/eta(1)
!!!!!!!!!!!!!!!!!!!!CALCULO DE SK
    Do i1=1, kpoints
       J0=BESSEL_JN(0, k(i1)/2.d0)
       J1=BESSEL_JN(1, k(i1)/2.d0)
       J12=BESSEL_JN(1, k(i1))
       Sk(i1,1,1)=A*((2.d0*J1/k(i1))**2.d0)+(B*2.d0*J0*J1/k(i1))+G*2.d0*J12/k(i1)
       Sk(i1,1,1)=-Sk(i1,1,1)*4.d0*eta(1)
       Sk(i1,1,1)=1.d0/(1.d0-Sk(i1,1,1))
    End Do
  End Subroutine calc_sk_hd_monocomponent_Rosenfeld

  !Sharma Sharma approximation for Attractive Potential System
  Subroutine sharma_sharma_attractive_yukawa_2D(Temp)
    Implicit None
    Real * 8, Intent(in) :: Temp
    Real * 8:: J0,Z1,Z2,X,y,uyuk,ck_dum
    Real * 8, Parameter :: dy=1.d-4
    Integer :: i1
    If (SDimen==2) Then
       Call calc_sk_hd_monocomponent_Rosenfeld()
       Do i1=1, kpoints
          J0=BESSEL_JN(0, k(i1)/2.d0)
          ck_dum=1.d0-(1.d0/sk(i1,1,1))	!Ck HD!!!!
          !FOURIER TRANSFORM OF beta u(r)
          y=sigma(1)
          uyuk=0.d0
          Do While (exp(-z_yuk(1)*(y)/sigma(1))>dy)
             J0=BESSEl_JN(0,k(i1)*(y))
             uyuk=uyuk+J0*exp(-(y/sigma(1))*z_yuk(1))
             y=y+dy
          End Do
          uyuk=-8.d0*eta(1)*uyuk*(dy/sigma(1))*exp(z_yuk(1))/Temp!!!!!!!!!!!!!!!!!!!!!!!!!!!rho beta u
          ck_dum=ck_dum-uyuk
          sk(i1,1,1)=1.d0/(1.d0-ck_dum)
       End Do
    Else If (SDimen==3) Then       
       Call Calc_Sk_hs_py_mono(.True.)
       ! print *, 'Reprint SK'
       ! print *, '-------------------------'
       Do i1=1, kpoints          
          ck_dum=(1.d0/sk(i1,1,1))	!Ck HS!!!!
          ! print *, 'ck = ', ck_dum          
          uyuk=24.d0*eta(1)*((k(i1)*cos(k(i1)))+(z_yuk(1)*sin(k(i1)))) / (Temp*k(i1)&
               &*((k(i1)**2.d0)+(z_yuk(1)**2.d0)))
          ! print *, 'uyuk = ', uyuk          
          ! uyuk=24.d0*eta(1)*((k(i1)*cos(k(i1)))+(z_yuk(1)*sin(k(i1))))/(Temp*k(i1)&
          !      &*((k(i1)**2.d0)+(z_yuk(1)**2.d0)))
          ck_dum=ck_dum-uyuk
          Sk(i1,1,1)=1.d0/(ck_dum)
          ! print *, 'SK_Yukis(', i1, ') = ',Sk(i1,1,1)
       End Do
    End If

  End Subroutine sharma_sharma_attractive_yukawa_2D

End Module Structure_module
