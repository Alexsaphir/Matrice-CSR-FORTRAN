! MODULE: CONJUGATE_GRADIENT
!> @author
!> Mechineau Alexandre
! DESCRIPTION:
!> Ce module permet d'utiliser l'algorithme du gradient conjugué. Deux implémentations sont fournies, la première effectue le nombre maximal d'étape, la seconde s'arrête à une précision donnée. 

MODULE CONJUGATE_GRADIENT
USE CSR
IMPLICIT NONE
CONTAINS
        
        !> Résouds le système \f$AX=B\f$ en utilisant la méthode du gradiant conjugué
        !! Source https://www.ljll.math.upmc.fr/hecht/ftp/InfoSci/doc-pdf/Master_book_GC.pdf
        !! \param A Matrice CSR
        !! \param B Matrice CSR représentant un vecteur colonne
        !! \param eps Précison requise (Optionel)
        !! \return Approximation, dans le pire des cas, de la solution du système
        FUNCTION CONJUGATE_GRADIENT_DIRECT(A, B, eps) RESULT(X)
                TYPE(CSR_MATRIX), INTENT(in)    :: A
                TYPE(CSR_MATRIX), INTENT(in)    :: B
                REAl, OPTIONAL, INTENT(in)      :: eps
                
                TYPE(CSR_MATRIX) :: X
                
                REAL(KIND=selected_real_kind(15, 307))                  :: alpha
                REAL(KIND=selected_real_kind(15, 307))                  :: beta
                REAL(KIND=selected_real_kind(15, 307))                  :: ps_R_R
                REAL(KIND=selected_real_kind(15, 307))                  :: ps_Rn_Rn
                TYPE(CSR_MATRIX)        :: R
                TYPE(CSR_MATRIX)        :: P
                INTEGER                         :: i
                
                X = CREATE_CSR_MATRIX(A%m_m, 1)! Initial guess set to 0
                R = CREATE_CSR_MATRIX(A%m_m, 1)
                P = CREATE_CSR_MATRIX(A%m_m, 1)
                
                ! Initialize R_0
                R = CSR_MATRIX_VECTOR_SUB(B, CSR_MATRIX_PROD(A, X))
                !Initialize P_0
                P = CREATE_CSR_MATRIX_RECOPY(R)
                DO i = 0, B%m_m, 1
                        alpha = CSR_MATRIX_INNER_PRODUCT(R, R) / CSR_MATRIX_INNER_PRODUCT(P, CSR_MATRIX_PROD(A, P))
                        X = CSR_MATRIX_VECTOR_ADD(X, CSR_MATRIX_SCALAR_PROD(P, alpha))
                        ps_R_R = CSR_MATRIX_INNER_PRODUCT(R, R)
                        R = CSR_MATRIX_VECTOR_SUB(R, CSR_MATRIX_SCALAR_PROD(CSR_MATRIX_PROD(A, P), alpha))
                        ps_Rn_Rn = CSR_MATRIX_INNER_PRODUCT(R, R)
                        beta = ps_Rn_Rn / ps_R_R
                        P = CSR_MATRIX_VECTOR_ADD(R, CSR_MATRIX_SCALAR_PROD(P, beta) )
                        IF (PRESENT(eps)) THEN
                                IF (ps_Rn_Rn < eps) THEN
                                        RETURN
                                ENDIF
                        ENDIF
                ENDDO
        END FUNCTION CONJUGATE_GRADIENT_DIRECT
        

END MODULE CONJUGATE_GRADIENT
