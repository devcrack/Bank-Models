Module Writting
  Use Dynamics_Variables
  Use System_Variables
  Use Static_Variables
  Implicit None
Contains
  Subroutine Sk_writting (File_name)
    Use Static_Variables
    Use System_Variables
    Implicit None
    Integer :: i1,i2,i3
    Character(len=*), Intent (in) :: File_name
    Character(len=20) :: FMT1, FMT2, FMT3, char_dum
    Character(len=2) :: char_dum_i,char_dum_j
    Open (unit=12, File=File_name,Status="Replace")
    FMT1="(E20.13, 2X)"         ! This format is used for pritn Real numbers
    FMT2="(A20, 2X)"            ! This format is used for print Characters 
    Write (12,FMT2,advance="no") "#k         "
    ! This couple of circles execute one time for  Hard_sphere_dynamics
    Do i1=1,species
       write(char_dum_i,"(I2.1)") i1
       print *, 'fkin species =', species
       print *, 'i1 = ', i1
       Do i2=1,i1
          print *, 'i2 = ', i2
          print *, "What the Fuck have char_dum_j", char_dum_j
          ! char_dum_j: no tiene nada  ._.
          Write(char_dum_j,"(I2.1)") i2
          char_dum="0S(k)_"//char_dum_i//"_"//char_dum_j
          Write (12,FMT2,advance="no") char_dum
       End Do
    End Do
    Write(12,*) ''
    ! Here we write inside the file
    Do i1=1,kpoints
       Write (12,FMT1,advance='no') k(i1)
       print*,'K(',i1,') = ',k(i1) 
       Do i2=1,species
          Do i3=1, i2
             Write (12,FMT1,advance='no')sk(i1,i2,i3)
             print*,'SK(',i1,') = ',sk(i1,i2,i3) 
          End Do
       End Do
       Write(12,*) ''
    End Do
    Close (12)
  End Subroutine Sk_writting

  Subroutine Dynamic_writting_tau_dependant (File_name,ik,it,first_time,last_time)
    Use Dynamics_Variables
    Use System_Variables
    Implicit None
    Integer :: i1,i2,i3
    Integer, Intent(in) :: ik,it
    Logical, Intent(in) :: first_time, last_time
    Character(len=*), Intent (in) :: File_name
    Character(len=20) :: FMT1, FMT2, FMT3, char_dum
    Character(len=2) :: char_dum_i,char_dum_j
    If(first_time.Eqv..True.) Then
       Open (unit=12, File=File_name,Status="Replace")
    End If
    FMT1="(E20.13, 2X)"
    FMT2="(A20, 2X)"
    If(first_time.Eqv..True.) Then
       Write (12,FMT2,advance="no") "#tau       "
       Do i1=1,species
          write(char_dum_i,"(I2.1)") i1
          Do i2=1,i1
             Write(char_dum_j,"(I2.1)") i2
             char_dum="Del(tau)_"//char_dum_i//"_"//char_dum_j
             Write (12,FMT2,advance="no") char_dum
          End Do
       End Do
       Do i1=1,species
          write(char_dum_i,"(I2.1)") i1
          Do i2=1,i1
             Write(char_dum_j,"(I2.1)") i2
             char_dum="Fc(ik,tau)_"//char_dum_i//"_"//char_dum_j
             Write (12,FMT2,advance="no") char_dum
          End Do
       End Do
       Do i1=1,species
          write(char_dum_i,"(I2.1)") i1
          Do i2=1,i1
             Write(char_dum_j,"(I2.1)") i2
             char_dum="Fs(ik,tau)_"//char_dum_i//"_"//char_dum_j
             Write (12,FMT2,advance="no") char_dum
          End Do
       End Do
       Write(12,*) ''
    End If
    Do i1=it,medium_times
       Write (12,FMT1,advance='no') m_tau(i1)
       Do i2=1,species
          Do i3=1, i2
             Write (12,FMT1,advance='no')Del_z(i1,i2,i3)
          End Do
       End Do
       Do i2=1,species
          Do i3=1, i2
             Write (12,FMT1,advance='no')Fc(ik,i1,i2,i3)
          End Do
       End Do
       Do i2=1,species
          Do i3=1, i2
             Write (12,FMT1,advance='no')Fs(ik,i1,i2,i3)
          End Do
       End Do
       Write(12,*) ''
    End Do
    If(last_time.Eqv..True.) Then
       Close (12)
    End If
  End Subroutine Dynamic_writting_tau_dependant

  Subroutine Dynamic_writting_tau_dependant_multi_k_tests (File_name_Fc,File_name_Fs,ik,it,first_time,last_time,Fc_unit,Fs_unit)
    Implicit None
    Integer :: i1,i2,i3
    Integer, Intent(in) :: it,ik,Fc_unit,Fs_unit
    Logical, Intent(in) :: first_time, last_time
    Character(len=*), Intent (in) :: File_name_Fc,File_name_Fs
    Character(len=20) :: FMT1, FMT2, FMT3, char_dum
    Character(len=2) :: char_dum_i,char_dum_j
    print*, ik, k(ik)
    If(first_time.Eqv..True.) Then
       Open (unit=Fc_unit, File=File_name_Fc,Status="Replace")
    End If
    FMT1="(E20.13, 2X)"
    FMT2="(A20, 2X)"
    If(first_time.Eqv..True.) Then
       Write (Fc_unit,FMT2,advance="no") "#tau       "
       Do i1=1,species
          write(char_dum_i,"(I2.1)") i1
          Do i2=1,i1
             Write(char_dum_j,"(I2.1)") i2
             char_dum="Fc(ik,tau)_"//char_dum_i//"_"//char_dum_j
             Write (Fc_unit,FMT2,advance="no") char_dum
          End Do
       End Do
       Write(Fc_unit,FMT2,advance="no") "k_test="
       Write (Fc_unit,"(E20.13)",advance="no") k(ik)
       Write(Fc_unit,*) ''
    End If
    Do i1=it,medium_times
       Write (Fc_unit,FMT1,advance='no') m_tau(i1)
       Do i2=1,species
          Do i3=1, i2
             Write (Fc_unit,FMT1,advance='no')Fc(ik,i1,i2,i3)
          End Do
       End Do
       Write(Fc_unit,*) ''
    End Do
    If(last_time.Eqv..True.) Then
       Close (Fc_unit)
    End If
    !
    If(first_time.Eqv..True.) Then
       Open (unit=Fs_unit, File=File_name_Fs,Status="Replace")
    End If
    If(first_time.Eqv..True.) Then
       Write (Fs_unit,FMT2,advance="no") "#tau       "
       Do i1=1,species
          write(char_dum_i,"(I2.1)") i1
          Do i2=1,i1
             Write(char_dum_j,"(I2.1)") i2
             char_dum="Fs(ik,tau)_"//char_dum_i//"_"//char_dum_j
             Write (Fs_unit,FMT2,advance="no") char_dum
          End Do
       End Do
       Write(Fs_unit,FMT2,advance="no") "k_test="
       Write (Fs_unit,"(E20.13)",advance="no") k(ik)
       Write(Fs_unit,*) ''
    End If
    Do i1=it,medium_times
       Write (Fs_unit,FMT1,advance='no') m_tau(i1)
       Do i2=1,species
          Do i3=1, i2
             Write (Fs_unit,FMT1,advance='no')Fs(ik,i1,i2,i3)
          End Do
       End Do
       Write(Fs_unit,*) ''
    End Do
    If(last_time.Eqv..True.) Then
       Close (Fs_unit)
    End If
  End Subroutine Dynamic_writting_tau_dependant_multi_k_tests

  Subroutine Dynamic_writting_Delz (File_name,it,first_time,last_time)
    Use Dynamics_Variables
    Use System_Variables
    Implicit None
    Integer :: i1,i2,i3
    Integer, Intent(in) :: it
    Logical, Intent(in) :: first_time, last_time
    Character(len=*), Intent (in) :: File_name
    Character(len=20) :: FMT1, FMT2, FMT3, char_dum
    Character(len=2) :: char_dum_i,char_dum_j
    If(first_time.Eqv..True.) Then
       Open (unit=20, File=File_name,Status="Replace")
    End If
    FMT1="(E20.13, 2X)"
    FMT2="(A20, 2X)"
    If(first_time.Eqv..True.) Then
       Write (20,FMT2,advance="no") "#tau       "
       Do i1=1,species
          write(char_dum_i,"(I2.1)") i1
          Do i2=1,i1
             Write(char_dum_j,"(I2.1)") i2
             char_dum="Del(tau)_"//char_dum_i//"_"//char_dum_j
             Write (20,FMT2,advance="no") char_dum
          End Do
       End Do
       Write(20,*) ''
    End If
    Do i1=it,medium_times
       Write (20,FMT1,advance='no') m_tau(i1)
       Do i2=1,species
          Do i3=1, i2
             Write (20,FMT1,advance='no')Del_z(i1,i2,i3)
          End Do
       End Do
       Write(20,*) ''
    End Do
    If(last_time.Eqv..True.) Then
       Close (20)
    End If
  End Subroutine Dynamic_writting_Delz

  Subroutine Dynamic_writting_Fc (File_name,ik,it,first_time,last_time)
    Use Dynamics_Variables
    Use System_Variables
    Implicit None
    Integer :: i1,i2,i3
    Integer, Intent(in) :: ik,it
    Logical, Intent(in) :: first_time, last_time
    Character(len=*), Intent (in) :: File_name
    Character(len=20) :: FMT1, FMT2, FMT3, char_dum
    Character(len=2) :: char_dum_i,char_dum_j
    If(first_time.Eqv..True.) Then
       Open (unit=21, File=File_name,Status="Replace")
    End If
    FMT1="(E20.13, 2X)"
    FMT2="(A20, 2X)"
    If(first_time.Eqv..True.) Then
       Write (21,FMT2,advance="no") "#tau       "
       Do i1=1,species
          write(char_dum_i,"(I2.1)") i1
          Do i2=1,i1
             Write(char_dum_j,"(I2.1)") i2
             char_dum="Fc(ik,tau)_"//char_dum_i//"_"//char_dum_j
             Write (21,FMT2,advance="no") char_dum
          End Do
       End Do
       Write(21,*) ''
    End If
    Do i1=it,medium_times
       Write (21,FMT1,advance='no') m_tau(i1)
       Do i2=1,species
          Do i3=1, i2
             Write (21,FMT1,advance='no')Fc(ik,i1,i2,i3)
          End Do
       End Do
       Write(21,*) ''
    End Do
    If(last_time.Eqv..True.) Then
       Close (21)
    End If
  End Subroutine Dynamic_writting_Fc

  Subroutine Dynamic_writting_Fs (File_name,ik,it,first_time,last_time)
    Use Dynamics_Variables
    Use System_Variables
    Implicit None
    Integer :: i1,i2,i3
    Integer, Intent(in) :: ik,it
    Logical, Intent(in) :: first_time, last_time
    Character(len=*), Intent (in) :: File_name
    Character(len=20) :: FMT1, FMT2, FMT3, char_dum
    Character(len=2) :: char_dum_i,char_dum_j
    If(first_time.Eqv..True.) Then
       Open (unit=22, File=File_name,Status="Replace")
    End If
    FMT1="(E20.13, 2X)"
    FMT2="(A20, 2X)"
    If(first_time.Eqv..True.) Then
       Write (22,FMT2,advance="no") "#tau       "
       Do i1=1,species
          write(char_dum_i,"(I2.1)") i1
          Do i2=1,i1
             Write(char_dum_j,"(I2.1)") i2
             char_dum="Fs(ik,tau)_"//char_dum_i//"_"//char_dum_j
             Write (22,FMT2,advance="no") char_dum
          End Do
       End Do
       Write(22,*) ''
    End If
    Do i1=it,medium_times
       Write (22,FMT1,advance='no') m_tau(i1)
       Do i2=1,species
          Do i3=1, i2
             Write (22,FMT1,advance='no')Fs(ik,i1,i2,i3)
          End Do
       End Do
       Write(22,*) ''
    End Do
    If(last_time.Eqv..True.) Then
       Close (22)
    End If
  End Subroutine Dynamic_writting_Fs

  Subroutine Dynamic_writting_Diff_msd (File_name,it,first_time,last_time)
    Use Dynamics_Variables
    Use System_Variables
    Implicit None
    Integer :: i1,i2,i3
    Integer, Intent(in) :: it
    Logical, Intent(in) :: first_time, last_time
    Character(len=*), Intent (in) :: File_name
    Character(len=20) :: FMT1, FMT2, FMT3, char_dum
    Character(len=2) :: char_dum_i,char_dum_j
    If(first_time.Eqv..True.) Then
       Open (unit=23, File=File_name,Status="Replace")
    End If
    FMT1="(E20.13, 2X)"
    FMT2="(A20, 2X)"
    If(first_time.Eqv..True.) Then
       Write (23,FMT2,advance="no") "#tau       "
       Do i1=1,species
          write(char_dum_i,"(I2.1)") i1
          i2=i1
          Write(char_dum_j,"(I2.1)") i2
          char_dum="DiffCoeff(tau)_"//char_dum_i//"_"//char_dum_j
          Write (23,FMT2,advance="no") char_dum
       End Do
       Do i1=1,species
          write(char_dum_i,"(I2.1)") i1
          i2=i1
          Write(char_dum_j,"(I2.1)") i2
          char_dum="MSD(tau)_"//char_dum_i//"_"//char_dum_j
          Write (23,FMT2,advance="no") char_dum
       End Do
       Write(23,*) ''
    End If
    Do i1=it,medium_times
       Write (23,FMT1,advance='no') m_tau(i1)
       Do i2=1,species
          Write (23,FMT1,advance='no')DiffCoeff(i1,i2)
       End Do
       Do i2=1,species
          Write (23,FMT1,advance='no')MSD(i1,i2)
       End Do
       Write(23,*) ''
    End Do
    If(last_time.Eqv..True.) Then
       Close (23)
    End If
  End Subroutine Dynamic_writting_Diff_msd

  Subroutine SCGLE_run_info_file(File_name,gamma_values,System_kind,i_k_test)
    Use Dynamics_Variables
    Use Static_Variables
    Use SCGLE_Variables
    Use System_Variables
    Implicit None
    Integer, Intent(in) :: i_k_test
    Real * 8, Dimension(:,:), Intent(in) :: gamma_values
    Character(len=*), Intent (in) :: File_name,System_kind
    Integer :: i1
    Open (unit=30, File=File_name,Status="Replace")
    Write(30,*) System_kind, " System"
    print*, System_kind, " System"
    Write(30,*) "Species #=", Species
    print*,"Species #=", Species
    Write(30,*) "System Dimensionality=", SDimen
    print*, "System Dimensionality=", SDimen
    Write(30,*) "k points #=", kpoints
    print*, "k points #=", kpoints
    Write(30,*) "k differential=", dk
    print*, "k differential=", dk
    Write(30,*) "integer for k test=", i_k_test
    print*, "integer for k test=", i_k_test
    Write(30,*) "k test=",k(i_k_test)
    print*,
    Write(30,*) "tau points #=", medium_times
    print*, "k test=",k(i_k_test)
    Write(30,*) "small time points #=", short_times
    print*, "small time points #=", short_times
    Write(30,*) "decimation #=", decimations
    print*, "decimation #=", decimations
    Do i1=1,Species
       Write(30,*) "sigma",i1,"=",sigma(i1)
       print*, "sigma",i1,"=",sigma(i1)
    End Do
    Do i1=1,Species
       Write(30,*) "eta",i1,"=",eta(i1)
       print*, "eta",i1,"=",eta(i1)
    End Do
    Do i1=1,Species
       Write(30,*) "D0",i1,"=",D0M(i1,i1)
       print*, "D0",i1,"=",D0M(i1,i1)
    End Do
    Do i1=1,Species
       Write(30,*) "(limit MSD) Gamma",i1,"=",gamma_values(i1,i1)
       print*,  "(limit MSD) Gamma",i1,"=",gamma_values(i1,i1)
    End Do
    Do i1=1,Species
       Write(30,*) "Dl",i1,"=",Dl(i1)
       print*, "Dl",i1,"=",Dl(i1)
    End Do
    Do i1=1,Species
       Write(30,*) "kc",i1,"=",kc(i1)
       print*, "kc",i1,"=",kc(i1)
    End Do
    Close (30)
  End Subroutine SCGLE_run_info_file

End Module Writting
