MODULE Structure_module
  Use Parameter_Variables
  Use System_Variables
  Use Static_Variables
  Use Potential_Variables
  implicit None
CONTAINS
  !Subroutine for the calculation of a static number of equally spaced wave vectors
  SUBROUTINE  Calc_static_k(difk,knumber)
    implicit none
    real * 8, intent (in) :: difk
    integer,  intent (in) :: knumber
    integer :: i1
    do i1=1, knumber
       k(i1)=i1*difk            !Thi vector is defined in Module Static Variables       
       print *,'K(',i1,') = ',k(i1)
       ! in the file variables_modules.f90
    end do
  END SUBROUTINE Calc_static_k

  SUBROUTINE Calc_Sk_hs_py_mono(vw_option)
    implicit None
    integer :: i1
    logical :: vw_option
    real * 8, Dimension(:), Allocatable :: k_vw
    real * 8 :: eta_vw,dum1,dum2,dum3,dum4,dumsin,dumcos
    !Start of parameter checks
    if (SDimen/= 3 .or. Species /=1 ) Then
       print*, ('Error, cannot compute the system direct correlation function')
       Stop
    End If

    
    !End of parameter checks
    Allocate(k_vw(kpoints))
    ! print *, k(:)               
    !Start of VW Correction
    If (vw_option .EQV. .TRUE.) Then
       eta_vw=eta(1)*(1.d0-(eta(1)/16.d0)) ![x]
       print *, 'eta(1) = ', eta(1)
       ! print *, 'eta_vw = ', eta_vw
       ! print *, 'eta = ', eta
       ! print *, 'What the fuck contains k(:)', k(:)*((eta_vw/eta(1))**(1.d0/3.d0)) 
       k_vw(:)=k(:)*((eta_vw/eta(1))**(1.d0/3.d0)) ![X]
       ! do i1 = 1, kpoints
       !    print *, 'Value of k_vw[',i1,'] = ', k_vw(i1)
       ! end do
    Else
       eta_vw=eta(1)            ![X]
       k_vw(:)=k(:)             ![X]
    End If
    !End of VW correction
    ! print *,'calc=', -(1.0+2.0*eta_vw)**2
    dum1=-((1.0+2.0*eta_vw)**2)
    ! print *, 'Value of Dum1 =', dum1
    dum2=6.0*eta_vw*((1.0+(eta_vw/2.0))**2) ![X]
    ! print *, 'Value of Dum2 =', dum2
    dum3=-eta_vw*((1.0+(2.0*eta_vw))**2)/2.0
    ! print *, 'Value of Dum3 =', dum3
    dum4=((1.0-eta_vw)**4)
    print *, 'Value of Dum4 =', dum4
    dum1=dum1/dum4
    print *, 'Value of Dum1 =', dum1
    dum2=dum2/dum4
    print *, 'Value of Dum2 =', dum2
    dum3=dum3/dum4
    print *, 'Value of Dum3 =', dum3
    do i1=1, kpoints
       dumsin=sin(k_vw(i1))     ![X]
       ! print *, 'Calc Sin[',i1,']',dumsin
       dumcos=cos(k_vw(i1))     ![X]
       ! print *, 'Calc Cos[',i1,']',dumcos

       ! Super Calculate
       ! ####################################################################################
       ! print *, 'Iter #', i1,' Part_1 Super Calculate = ', dum1*(dumsin-k_vw(i1)*dumcos)/(k_vw(i1)**2.d0)
       ck(i1,1,1)= (dum1*(dumsin-k_vw(i1)*dumcos)/(k_vw(i1)**2.d0))+&
            &(dum2*(((2.d0*k_vw(i1))*dumsin)+((-(k_vw(i1)**2.d0)+2.d0)*dumcos)-2.d0)/&
            &(k_vw(i1)**3.d0))+(dum3*(((4.d0*(k_vw(i1)**3.d0)-24.0*k_vw(i1))*dumsin)&
            &+((-(k_vw(i1)**4.d0)+12.d0*(k_vw(i1)**2.d0)-24.d0)*dumcos)+24.d0)/(k_vw(i1)**5))
       Ck(i1,1,1)=  24.d0*eta_vw*ck(i1,1,1)/k_vw(i1)
       ! print *, 'Iter #', i1,'Ck = ',Ck(i1,1,1)
       Sk(i1,1,1)=1.d0/(1.d0-ck(i1,1,1))
       ! print *, 'Iter #', i1,'Sk = ',Sk(i1,1,1)
       ! print *, Sk(i1,1,1)      
       Ski(i1,1,1)=1.d0-Ck(i1,1,1)
       ! print *, 'Iter #', i1,'Ski = ',Ski(i1,1,1)
    end do
  END SUBROUTINE Calc_Sk_hs_py_mono


END MODULE Structure_module
