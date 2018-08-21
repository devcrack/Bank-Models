Module Matrix_operations
Implicit None
Contains
Subroutine identity_matrix(nxn,matrix_out)
  Implicit None
  Integer, Intent (in) :: nxn
  Integer :: i1,i2
  Real * 8, Dimension (:,:), Intent (out) :: matrix_out
  Do i1=1, nxn
    Do i2=1, nxn
      matrix_out(i1,i2)=0.d0
      If(i1==i2)Then
        matrix_out(i1,i2)=1.d0
      End If
    End do
  End do
End Subroutine

Subroutine null_matrix(nxn,matrix_out)
  Implicit None
  Integer, Intent (in) :: nxn
  Integer :: i1,i2
  Real * 8, Dimension (:,:), Intent (out) :: matrix_out
  Do i1=1, nxn
    Do i2=1, nxn
      matrix_out(i1,i2)=0.d0
    End do
  End do
End Subroutine
!
Subroutine minor_matrix(nxn,matrix_in,matrix_out,i,j)
  Implicit None
  Integer, Intent (in) :: nxn,i,j
  Integer :: dumi,dumj
  Real * 8, Dimension (:,:), Intent (in) :: matrix_in
  Real * 8, Dimension (:,:), Intent (out) :: matrix_out
  Integer :: i1,i2
  Do i1=1,nxn
    Do i2=1, nxn
      If (i1/=i.And.i2/=j) Then
        dumi=i1
        dumj=i2
        If (i1>i) Then
          dumi=dumi-1
        End If
        If (i2>j) Then
          dumj=dumj-1
        End If
        matrix_out(dumi,dumj)=matrix_in(i1,i2)
      End If
    End Do
  End Do
End Subroutine
!
Recursive Subroutine Determinant(nxn,matrix_in,Det)
  Implicit None
  Integer, Intent (in) :: nxn
  Real * 8, Dimension (:,:), Intent (in) :: matrix_in
  Real * 8, Intent (out) :: Det
  Real * 8, Dimension (:,:), Allocatable :: m_matrix
  Real * 8 :: dummy_Det
  Integer :: i1,i2
  If (nxn==1) Then
    Det=matrix_in(1,1)
  Else If(nxn==2) Then
    Det=matrix_in(1,1)*matrix_in(2,2)-matrix_in(1,2)*matrix_in(2,1)
  Else
    Det=0
    Allocate(m_matrix(nxn-1,nxn-1))
    Do i1=1, nxn
      Call minor_matrix(nxn,matrix_in,m_matrix,1,i1)
      Call Determinant(nxn-1,m_matrix,dummy_Det)
      Det=Det+((-1.d0)**(i1+1))*dummy_Det
    End Do
  End If
End Subroutine
!
Subroutine cofactor_matrix(nxn,matrix_in,matrix_out)
    Implicit None
    Integer, Intent (in) :: nxn
    Integer i1,i2
    Real * 8, Dimension (:,:), Intent (in) :: matrix_in
    Real * 8, Dimension (:,:), Intent (out) :: matrix_out
    Real * 8, Dimension (:,:), Allocatable :: m_matrix
    Real * 8 :: Determ
    Allocate(m_matrix(nxn-1,nxn-1))
    Do i1=1, nxn
      Do i2=1, nxn
        Call minor_matrix(nxn,matrix_in,m_matrix,i1,i2)
        Call Determinant(nxn-1,m_matrix,Determ)
        matrix_out(i1,i2)=((-1.d0)**(i1+i2))*Determ
      End Do
    End Do
End Subroutine
!
Subroutine transpose_matrix(nxn,matrix_in,matrix_out)
  Implicit None
  Integer, Intent (in) :: nxn
  Integer i1,i2
  Real * 8, Dimension (:,:), Intent (in) :: matrix_in
  Real * 8, Dimension (:,:), Intent (out) :: matrix_out
  Do i1=1, nxn
    Do i2=1, nxn
      matrix_out(i1,i2)=matrix_in(i2,i1)
    End Do
  End Do
End Subroutine
!
Subroutine inverse_matrix(nxn,matrix_in,matrix_out)
  Implicit None
  Integer, Intent (in) :: nxn
  Integer i1,i2
  Real * 8, Dimension (:,:), Intent (in) :: matrix_in
  Real * 8, Dimension (:,:), Intent (out) :: matrix_out
  Real * 8, Dimension (:,:), Allocatable :: cof_matrix,trans_cof_matrix
  Real * 8 :: Det
  If(nxn>1) Then
    Allocate (cof_matrix(nxn,nxn))
    Allocate (trans_cof_matrix(nxn,nxn))
    Call cofactor_matrix(nxn,matrix_in,cof_matrix)
    Call transpose_matrix(nxn,cof_matrix,trans_cof_matrix)
    Call Determinant(nxn,matrix_in,Det)
    Do i1=1, nxn
      Do i2=1, nxn
        matrix_out(i1,i2)=trans_cof_matrix(i1,i2)/Det
      End Do
    End Do
  Else If(nxn==1) Then
    matrix_out(1,1)=1.d0/matrix_in(1,1)
  End If
End Subroutine
!
Subroutine matrix_multiplication(nxn,matrix_in_1,matrix_in_2,matrix_out)
  Implicit None
  Integer, Intent (in) :: nxn
  Integer i1,i2,i3
  Real * 8, Dimension (:,:), Intent (in) :: matrix_in_1,matrix_in_2
  Real * 8, Dimension (:,:), Intent (out) :: matrix_out
  Real * 8 ::dummy1
  Do i1=1, nxn
    Do i2=1, nxn
      dummy1=0.d0
      Do i3=1,nxn
        dummy1=dummy1+matrix_in_1(i1,i3)*matrix_in_2(i3,i2)
      End Do
      matrix_out(i1,i2)=dummy1
    End Do
  End Do
End Subroutine
!
Subroutine matrix_constant_multiplication(nxn,matrix_in,constant_in,matrix_out)
  Implicit None
  Integer, Intent (in) :: nxn
  Integer i1,i2,i3
  Real * 8, Dimension (:,:), Intent (in) :: matrix_in
  Real * 8, Dimension (:,:), Intent (out) :: matrix_out
  Real * 8 :: constant_in
  Real * 8 ::dummy1
  Do i1=1, nxn
    Do i2=1, nxn
      matrix_out(i1,i2)=constant_in*matrix_in(i1,i2)
    End Do
  End Do
End Subroutine
!
Subroutine matrix_conversion_3_to_2(idum,nxn,matrix_in,matrix_out)
  Implicit None
  Integer, Intent (in) :: nxn,idum
  Integer i1,i2,i3
  Real * 8, Dimension (:,:,:), Intent (in) :: matrix_in
  Real * 8, Dimension (:,:), Intent (out) :: matrix_out
  Real * 8 ::dummy1
  Do i1=1, nxn
    Do i2=1, nxn
      matrix_out(i1,i2)=matrix_in(idum,i1,i2)
    End Do
  End Do
End Subroutine


Subroutine matrix_conversion_2_to_3(idum,nxn,matrix_in,matrix_out)
  Implicit None
  Integer, Intent (in) :: nxn,idum
  Integer i1,i2,i3
  Real * 8, Dimension (:,:), Intent (in) :: matrix_in
  Real * 8, Dimension (:,:,:), Intent (out) :: matrix_out
  Real * 8 ::dummy1
  Do i1=1, nxn
    Do i2=1, nxn
      matrix_out(idum,i1,i2)=matrix_in(i1,i2)
    End Do
  End Do
End Subroutine

Subroutine matrix_conversion_3_to_1(idum1,idum2,points,matrix_in,vector_out)
  Implicit None
  Integer, Intent (in) :: points,idum1,idum2
  Integer i1,i2,i3
  Real * 8, Dimension (:,:,:), Intent (in) :: matrix_in
  Real * 8, Dimension (:), Intent (out) :: vector_out
  Real * 8 ::dummy1
  Do i1=1,points
      vector_out(i1)=matrix_in(i1,idum1,idum2)
  End Do
End Subroutine

Subroutine matrix_conversion_2_to_4(idum1,idum2,nxn,matrix_in,matrix_out)
  Implicit None
  Integer, Intent (in) :: idum1,idum2,nxn
  Integer i1,i2
  Real * 8, Dimension (:,:), Intent (in) :: matrix_in
  Real * 8, Dimension (:,:,:,:), Intent (out) :: matrix_out
  Real * 8 ::dummy1
  Do i1=1,nxn
    Do i2=1,nxn
      matrix_out(idum1,idum2,i1,i2)=matrix_in(i1,i2)
    End Do
  End Do
End Subroutine

Subroutine matrix_conversion_4_to_2(idum1,idum2,nxn,matrix_in,matrix_out)
  Implicit None
  Integer, Intent (in) :: idum1,idum2,nxn
  Integer i1,i2
  Real * 8, Dimension (:,:,:,:), Intent (in) :: matrix_in
  Real * 8, Dimension (:,:), Intent (out) :: matrix_out
  Real * 8 ::dummy1
  Do i1=1,nxn
    Do i2=1,nxn
      matrix_out(i1,i2)=matrix_in(idum1,idum2,i1,i2)
    End Do
  End Do
End Subroutine

End Module
