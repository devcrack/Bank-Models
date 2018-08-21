Module Integration_methods
  Implicit None
  Contains

  Real * 8 Function rectangle_fixed_differential(integrate_vector,differential)
    Real * 8, Intent(in) :: differential
    Real * 8, Dimension (:), Intent (in) :: integrate_vector
    rectangle_fixed_differential=sum(integrate_vector)*differential
  End Function
End Module
