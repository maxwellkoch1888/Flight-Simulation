module math_m 
    implicit none 

    contains 
  !=========================
  ! LU Decomposition from MachLine GitHub
    subroutine lu_solve(N, A, b, x)
    ! Solves a general [A]x=b on an nxn matrix
    ! This replaces A (in place) with its LU decomposition (permuted row-wise)

      implicit none

      integer,intent(in) :: N
      real,dimension(N,N),intent(inout) :: A
      real,dimension(N),intent(in) :: b
      real,dimension(:),allocatable,intent(out) :: x

      integer,allocatable,dimension(:) :: indx
      integer :: D, info

      allocate(indx(N))

      ! Compute decomposition
      call lu_decomp(A, N, indx, D, info)

      ! if the matrix is nonsingular, then backsolve to find X
      if (info == 1) then
          write(*,*) 'Subroutine lu_decomp() failed. The given matrix is singular (i.e. no unique solution). Quitting...'
          stop
      else
          call lu_back_sub(A, N, indx, b, x)
      end if

      ! Cleanup
      deallocate(indx)

  end subroutine lu_solve

  !=========================
  ! Matrix Inverse
    function matrix_inv(A) result(A_inv)
      implicit none
      real, dimension(3,3), intent(in) :: A
      real, dimension(3,3) :: A_inv, A_adj
      real :: det

      ! CALCULATE THE DETERMINANT
      det = A(1,1) * A(2,2) * A(3,3) + A(1,2) * A(2,3) * A(3,1) + A(1,3) * A(2,1) * A(3,2) &
          - A(1,3) * A(2,2) * A(3,1) - A(1,2) * A(2,1) * A(3,3) - A(1,1) * A(2,3) * A(3,2)

      ! CALCULATE THE ADJUGATE MATRIX
      A_adj(1,1) = A(2,2) * A(3,3) - A(2,3) * A(3,2)
      A_adj(1,2) = A(1,3) * A(3,2) - A(1,2) * A(3,3)
      A_adj(1,3) = A(1,2) * A(2,3) - A(1,3) * A(2,2)

      A_adj(2,1) = A(2,3) * A(3,1) - A(2,1) * A(3,3)
      A_adj(2,2) = A(1,1) * A(3,3) - A(1,3) * A(3,1)
      A_adj(2,3) = A(1,3) * A(2,1) - A(1,1) * A(2,3)

      A_adj(3,1) = A(2,1) * A(3,2) - A(2,2) * A(3,1)
      A_adj(3,2) = A(1,2) * A(3,1) - A(1,1) * A(3,2)
      A_adj(3,3) = A(1,1) * A(2,2) - A(1,2) * A(2,1)

      ! CALCULATE THE INVERSE MATRIX
      A_inv = 1/det * A_adj
    end function matrix_inv

  !=========================
  ! CROSS PRODUCT
    function cross_product(vector_a, vector_b) result(vector_c)
      real, intent(in) :: vector_a(3), vector_b(3)
      real :: vector_c(3)
      real :: a1, a2, a3, b1, b2, b3

      ! BREAK DOWN VECTOR 1 AND 2
      a1 = vector_a(1)
      a2 = vector_a(2)
      a3 = vector_a(3)

      b1 = vector_b(1)
      b2 = vector_b(2)
      b3 = vector_b(3)

      ! CALCULATE ORTHOGONAL VECTOR
      vector_c(1) = a2*b3 - a3*b2
      vector_c(2) = a3*b1 - a1*b3
      vector_c(3) = a1*b2 - a2*b1
    end function cross_product
end module math_m