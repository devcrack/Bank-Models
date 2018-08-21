Module Data_ask
  Implicit None

  Contains
    !User ask parameters!
    Subroutine user_kpoint_ask
      Use Variables
      Implicit None
      print*, 'Please give me the number of wave vectors that you wish to use'
      read(*, "(I10)") kpoints
    End Subroutine

    Subroutine user_species_ask
      Use Variables
      Implicit None
      print*, 'Please give me the number of species that you wish to use'
      read(*, "(I10)") species
    End Subroutine

    Subroutine user_dimension_ask
      Use Variables
      Implicit None
      print*, 'Please give me the dimension that you wish to use'
      read(*, "(I10)") Dimen
    End Subroutine

    Subroutine user_eta_ask
      Use Variables
      Implicit None
      print*, 'Please give me the value of the space fraction that you wish to use'
      read(*, "(F10.9)") eta
    End Subroutine
    !!
End Module
