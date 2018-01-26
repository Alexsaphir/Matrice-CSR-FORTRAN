! MODULE: SOLVER
!> @author
!> Mechineau Alexandre
! DESCRIPTION:
!> Ce module permet de résoudre des systèmes linéaires en utilisant soit la décomposition LU, soit la méthode du gradient conjugué.

MODULE SOLVER
USE CSR
USE LU
USE CONJUGATE_GRADIENT
IMPLICIT NONE
CONTAINS
        
        !> Résouds le système \f$AX=B\f$ en utilisant la décomposition LU de A
        !! \param A Matrice CSR
        !! \param B Matrice CSR représentant un vecteur colonne
        !! \return Approximation, dans le pire des cas de la solution du système
        FUNCTION SOLVER_LU(A , B) RESULT(X)
                TYPE(CSR_MATRIX), intent(in) :: A
                TYPE(CSR_MATRIX), intent(in) :: B
                
                TYPE(CSR_MATRIX) :: X
                
                TYPE(CSR_MATRIX) :: LU(2)
                
                LU = LU_CSR_MATRIX_SIMPLE(A)
                X = SOLVER_FORWARD_SUBSTITUTION(LU(1), B)
                X = SOLVER_BACKWARD_SUBSTITUTION(LU(2), X)
        END FUNCTION SOLVER_LU
        
        !> Résouds le système \f$AX=B\f$ en utilisant la méthode du gradiant conjugué
        !! \param A Matrice CSR
        !! \param B Matrice CSR représentant un vecteur colonne
        !! \param eps Précison requise (Optionel)
        !! \return Approximation, dans le pire des cas, de la solution du système
        FUNCTION SOLVER_GC(A, B, eps) RESULT(X)
                TYPE(CSR_MATRIX), INTENT(in)    :: A
                TYPE(CSR_MATRIX), INTENT(in)    :: B
                REAl, OPTIONAL, INTENT(in)      :: eps
                
                TYPE(CSR_MATRIX) :: X
                
                IF (PRESENT(eps)) THEN
                        X = CONJUGATE_GRADIENT_DIRECT(A, B, eps)
                ELSE
                        X = CONJUGATE_GRADIENT_DIRECT(A, B)
                ENDIF
        END FUNCTION SOLVER_GC
        
        !> Résouds le système \f$LX=B\f$ un utilisant la méthode substitution
        !! \param L Matrice CSR triangulaire inférieur
        !! \param B Matrice CSR représentant un vecteur colonne
        !! \returns Solution du système
        FUNCTION SOLVER_FORWARD_SUBSTITUTION(L, B) RESULT(X)
                TYPE(CSR_MATRIX), intent(in) :: L
                TYPE(CSR_MATRIX), intent(in) :: B
                
                TYPE(CSR_MATRIX) :: X
                
                INTEGER         :: i
                INTEGER         :: m
                REAL(KIND=selected_real_kind(15, 307))  :: tmp
                
                X = CREATE_CSR_MATRIX(L%m_m, 1)
                
                DO m = 1, B%m_m, 1
                        tmp = 0.D0
                        DO i = 1, m - 1, 1
                                tmp = tmp + GET(L, m, i)*GET(X, i, 1)
                        ENDDO
                        tmp = (GET(B, m, 1) - tmp) / GET(L, m, m)
                        CALL SET(X, m, 1, tmp)
                ENDDO
        END FUNCTION SOLVER_FORWARD_SUBSTITUTION
        
        !> Résouds le système \f$LX=B\f$ un utilisant la méthode substitution inversée
        !! \param L Matrice CSR triangulaire supérieur
        !! \param B Matrice CSR représentant un vecteur colonne
        !! \returns Solution du système
        FUNCTION SOLVER_BACKWARD_SUBSTITUTION(L, B) RESULT(X)
                TYPE(CSR_MATRIX), intent(in) :: L
                TYPE(CSR_MATRIX), intent(in) :: B
                
                TYPE(CSR_MATRIX) :: X
                
                INTEGER         :: i
                INTEGER         :: m
                REAL(KIND=selected_real_kind(15, 307))  :: tmp
                
                X = CREATE_CSR_MATRIX(L%m_m, 1)
                
                DO m = B%m_m, 1, -1
                        tmp = 0.D0
                        DO i = B%m_m, m + 1, -1
                                tmp = tmp + GET(L, m, i)*GET(X, i, 1)
                        ENDDO
                        tmp = (GET(B, m, 1) - tmp) / GET(L, m, m)
                        CALL SET(X, m, 1, tmp)
                ENDDO
        END FUNCTION SOLVER_BACKWARD_SUBSTITUTION

END MODULE SOLVER
