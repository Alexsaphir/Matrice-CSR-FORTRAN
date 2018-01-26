! MODULE: LU
!> @author
!> Mechineau Alexandre
! DESCRIPTION:
!> Ce module permet d'obtenir la décomposition LU d'une matrice.

MODULE LU
USE CSR
IMPLICIT NONE
CONTAINS

        !> Retourne la décomposition LU d'une matrice selon l'algorithme Doolittle sans permutation.
        !! \param M Matrice CSR
        !! \return Tableau de 2 matrice CSR (L, U)
        FUNCTION LU_CSR_MATRIX_SIMPLE(M)
                ! Use the Doolittle algorithm
                
                TYPE(CSR_MATRIX), intent(in) :: M
                
                TYPE(CSR_MATRIX) :: LU_CSR_MATRIX_SIMPLE(2)
                !Assert Matrix of size n*n
                
                TYPE(CSR_MATRIX)        :: U
                TYPE(CSR_MATRIX)        :: Ln
                TYPE(CSR_MATRIX)        :: L
                INTEGER                         :: i
                INTEGER                         :: j
                
                L = CREATE_CSR_MATRIX_SQUARE(M%m_m)
                
                U = M ! A0
                DO i = 1, M%m_m-1, 1
                        Ln = LN_CSR_MATRIX(U, i)
                        ! Compute U
                        U = CSR_MATRIX_PROD(Ln, U)
                        !Fill L with the i column
                        DO j = i, M%m_m, 1
                                IF (j == i) THEN
                                        CALL SET(L, j, i, 1.D0)
                                ELSE
                                        CALL SET(L, j, i, -GET(Ln, j, i))
                                ENDIF
                        ENDDO
                ENDDO
                CALL SET(L, M%m_m, M%m_m, 1.D0)
                LU_CSR_MATRIX_SIMPLE = [ L, U ]
        END FUNCTION LU_CSR_MATRIX_SIMPLE
        
        !> Fonction intérmédiare dans le calcul de LU
        FUNCTION LN_CSR_MATRIX(A, n)
                TYPE(CSR_MATRIX), intent(in)    :: A
                INTEGER, intent(in)             :: n
                
                TYPE(CSR_MATRIX)        :: LN_CSR_MATRIX
                !ASSERT GET(A, n, n) !=0
                
                INTEGER :: i
                
                !Init LN
                LN_CSR_MATRIX = CREATE_CSR_MATRIX(A%m_m, A%m_m)
                !Compute LN
                !Set 1 on the diag
                DO i = 1, A%m_m, 1
                        CALL SET(LN_CSR_MATRIX, i, i, 1.D0)
                ENDDO
                DO i = n + 1, A%m_m, 1
                        CALL SET( LN_CSR_MATRIX, i, n, -GET(A, i, n) / GET(A, n, n) )
                ENDDO
        END FUNCTION LN_CSR_MATRIX

END MODULE LU
