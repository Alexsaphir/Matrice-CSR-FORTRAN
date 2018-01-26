! MODULE: UNIT_TEST
!> @author
!> Mechineau Alexandre
! DESCRIPTION:
!> Ce module contient plusieurs test permettant d'assurer le bon fonctionnement des modules CSR, LU, CONJUGATE_GRADIENT et SOLVER.

MODULE UNIT_TEST
USE CSR
USE LU
USE CONJUGATE_GRADIENT
USE SOLVER

IMPLICIT NONE
CONTAINS
        !> Cette subroutine appelle l'ensemble des tests implémentés et affiche chaqu'un des résultats.
        SUBROUTINE TEST()
                WRITE(*,*) "SET test 1", TEST_SET_1()
                WRITE(*,*) "LU test 1", TEST_LU_1()
                WRITE(*,*) "LU test 2", TEST_LU_2()
                WRITE(*,*) "SOLVER Forward test 1", TEST_SOLVER_FORWARD_1()
                WRITE(*,*) "SOLVER Backward test 1", TEST_SOLVER_BACKWARD_1()
                WRITE(*,*) "SOLVER LU test 1", TEST_SOLVER_LU_1()
                WRITE(*,*) "SOLVER LU test 2", TEST_SOLVER_LU_2()
                WRITE(*,*) "Transpose test 1", TEST_TRANSPOSE_1()
                WRITE(*,*) "Transpose test 2", TEST_TRANSPOSE_2()
                WRITE(*,*) "Inner Product test 1", TEST_INNER_PRODUCT_1()
                WRITE(*,*) "Recopy test 1", TEST_RECOPY()
        END SUBROUTINE

        FUNCTION TEST_SET_1()
                LOGICAL :: TEST_SET_1
                
                TYPE(CSR_MATRIX) :: A
                TYPE(CSR_MATRIX) :: B
                
                A = CREATE_CSR_MATRIX(3,3)
                B = CREATE_CSR_MATRIX(3,3)

                CALL SET(A, 1, 1, 2.D0)
                
                CALL SET(A, 1, 3, -1.D0)
                CALL SET(A, 1, 2, -1.D0)
                CALL SET(A, 2, 1, 2.D0)
                CALL SET(A, 2, 2, -1.D0)
                CALL SET(A, 2, 3, -1.D0)
                CALL SET(A, 3, 1, -7.D0)
                CALL SET(A, 3, 2, 5.D0)
                CALL SET(A, 3, 3, 1.D0)
                
                CALL SET(A, 1, 1, 1.D0)
                CALL SET(A, 1, 2, 2.D0)
                CALL SET(A, 1, 3, 3.D0)
                CALL SET(A, 2, 1, 4.D0)
                CALL SET(A, 2, 2, 5.D0)
                CALL SET(A, 2, 3, 6.D0)
                CALL SET(A, 3, 1, 7.D0)
                CALL SET(A, 3, 2, 8.D0)
                CALL SET(A, 3, 3, 9.D0)
        
                CALL SET(B, 1, 1, 1.D0)
                CALL SET(B, 1, 2, 2.D0)
                CALL SET(B, 1, 3, 3.D0)
                CALL SET(B, 2, 1, 4.D0)
                CALL SET(B, 2, 2, 5.D0)
                CALL SET(B, 2, 3, 6.D0)
                CALL SET(B, 3, 1, 7.D0)
                CALL SET(B, 3, 2, 8.D0)
                CALL SET(B, 3, 3, 9.D0)
                TEST_SET_1 = CSR_MATRIX_EQUAL(A,B)
        
                CALL FREE_CSR_MATRIX(A)
                CALL FREE_CSR_MATRIX(B)
        END FUNCTION

        FUNCTION TEST_LU_1()
                LOGICAL :: TEST_LU_1
        
                TYPE(CSR_MATRIX) :: M
                TYPE(CSR_MATRIX) :: LU(2)
                
                M = CREATE_CSR_MATRIX(3,3)
                CALL SET(M, 1, 1, 2.D0)
                CALL SET(M, 1, 2, -1.D0)
                CALL SET(M, 2, 1, -1.D0)
                CALL SET(M, 2, 2, 2.D0)
                CALL SET(M, 2, 3, -1.D0)
                CALL SET(M, 3, 2, -1.D0)
                CALL SET(M, 3, 3, 2.D0)
        
                LU = LU_CSR_MATRIX_SIMPLE(M)
                
                TEST_LU_1 = CSR_MATRIX_EQUAL(M, CSR_MATRIX_PROD(LU(1), LU(2)))
        
                CALL FREE_CSR_MATRIX(M)
                CALL FREE_CSR_MATRIX(LU(1))
                CALL FREE_CSR_MATRIX(LU(2))
        END FUNCTION
        
        FUNCTION TEST_LU_2()
                LOGICAL :: TEST_LU_2
                TYPE(CSR_MATRIX) :: M
                TYPE(CSR_MATRIX) :: LU(2)
                
                M = CREATE_CSR_MATRIX(2,2)
                CALL SET(M, 1, 1, 4.D0)
                CALL SET(M, 1, 2, 3.D0)
                CALL SET(M, 2, 1, 6.D0)
                CALL SET(M, 2, 2, 3.D0)
        
                LU = LU_CSR_MATRIX_SIMPLE(M)
                
                TEST_LU_2 = CSR_MATRIX_EQUAL(M, CSR_MATRIX_PROD(LU(1), LU(2)))
        
                CALL FREE_CSR_MATRIX(M)
                CALL FREE_CSR_MATRIX(LU(1))
                CALL FREE_CSR_MATRIX(LU(2))
        END FUNCTION
        
        FUNCTION TEST_SOLVER_FORWARD_1()
                LOGICAL :: TEST_SOLVER_FORWARD_1
                
                TYPE(CSR_MATRIX) A
                TYPE(CSR_MATRIX) B
                TYPE(CSR_MATRIX) X
                
                TEST_SOLVER_FORWARD_1 = .FALSE.
                
                A = CREATE_CSR_MATRIX_SQUARE(2)
                B = CREATE_CSR_MATRIX(2, 1)
                
                CALL SET(A, 1, 1 , 1.D0)
                CALL SET(A, 2, 2, 1.D0)
                CALL SET(A, 2, 1 , 1.D0)
                
                CALL SET(B, 1, 1, 1.D0)
                CALL SET(B,2,1,2.D0)
                
                X = SOLVER_FORWARD_SUBSTITUTION(A,B)
                
                IF (CSR_MATRIX_EQUAL(CSR_MATRIX_PROD(A, X), B)) THEN
                        TEST_SOLVER_FORWARD_1 = .TRUE.
                ENDIF
                
                CALL FREE_CSR_MATRIX(A)
                CALL FREE_CSR_MATRIX(B)
                CALL FREE_CSR_MATRIX(X)
        END FUNCTION
        
        FUNCTION TEST_SOLVER_BACKWARD_1()
                LOGICAL :: TEST_SOLVER_BACKWARD_1
                
                TYPE(CSR_MATRIX) A
                TYPE(CSR_MATRIX) B
                TYPE(CSR_MATRIX) X
                
                TEST_SOLVER_BACKWARD_1 = .FALSE.
                
                A = CREATE_CSR_MATRIX_SQUARE(2)
                B = CREATE_CSR_MATRIX(2, 1)
                
                CALL SET(A, 1, 1 , 1.D0)
                CALL SET(A, 2, 2, 1.D0)
                CALL SET(A, 1, 2 , 1.D0)

                CALL SET(B, 1, 1, 1.D0)
                CALL SET(B,2,1,2.D0)
        
                X = SOLVER_BACKWARD_SUBSTITUTION(A,B)
        
                IF (CSR_MATRIX_EQUAL(CSR_MATRIX_PROD(A, X), B)) THEN
                        TEST_SOLVER_BACKWARD_1 = .TRUE.
                ENDIF
                
                CALL FREE_CSR_MATRIX(A)
                CALL FREE_CSR_MATRIX(B)
                CALL FREE_CSR_MATRIX(X)
        END FUNCTION
        

        FUNCTION TEST_SOLVER_LU_1()
                LOGICAL :: TEST_SOLVER_LU_1
                
                TYPE(CSR_MATRIX) A
                TYPE(CSR_MATRIX) B
                TYPE(CSR_MATRIX) X
                
                TEST_SOLVER_LU_1 = .FALSE.
                
                A = CREATE_CSR_MATRIX_SQUARE(2)
                B = CREATE_CSR_MATRIX(2, 1)
                
                CALL SET(A, 1, 1 , 1.D0)
                CALL SET(A, 2, 2, 1.D0)
                CALL SET(A, 2, 1 , 1.D0)
                
                CALL SET(B, 1, 1, 1.D0)
                CALL SET(B, 2, 1, 2.D0)
                
                X = SOLVER_LU(A, B)
                
                IF (CSR_MATRIX_EQUAL(CSR_MATRIX_PROD(A, X), B)) THEN
                        TEST_SOLVER_LU_1 = .TRUE.
                ENDIF
                
                CALL FREE_CSR_MATRIX(A)
                CALL FREE_CSR_MATRIX(B)
                CALL FREE_CSR_MATRIX(X)
        END FUNCTION
                
        FUNCTION TEST_SOLVER_LU_2()
                LOGICAL :: TEST_SOLVER_LU_2
                
                TYPE(CSR_MATRIX) A
                TYPE(CSR_MATRIX) B
                TYPE(CSR_MATRIX) X
                
                TEST_SOLVER_LU_2 = .FALSE.
                
                A = CREATE_CSR_MATRIX_SQUARE(2)
                B = CREATE_CSR_MATRIX(2, 1)
                
                CALL SET(A, 1, 1 , 1.D0)
                CALL SET(A, 2, 2, 1.D0)
                CALL SET(A, 1, 2 , 1.D0)
        
                CALL SET(B, 1, 1, 1.D0)
                CALL SET(B,2,1,2.D0)
        
                X = SOLVER_LU(A,B)
        
                IF (CSR_MATRIX_EQUAL(CSR_MATRIX_PROD(A, X), B)) THEN
                        TEST_SOLVER_LU_2 = .TRUE.
                ENDIF
                
                CALL FREE_CSR_MATRIX(A)
                CALL FREE_CSR_MATRIX(B)
                CALL FREE_CSR_MATRIX(X)
        END FUNCTION
        
        FUNCTION TEST_TRANSPOSE_1() 
                LOGICAL :: TEST_TRANSPOSE_1
        
                TYPE(CSR_MATRIX) B
                
                B = CREATE_CSR_MATRIX(2, 1)

                CALL SET(B, 1, 1, 1.D0)
                CALL SET(B,2,1,2.D0)
                
                TEST_TRANSPOSE_1 = CSR_MATRIX_EQUAL(B, TR(TR(B)))
                
                CALL FREE_CSR_MATRIX(B)
        END FUNCTION
        
        FUNCTION TEST_TRANSPOSE_2()
                LOGICAL :: TEST_TRANSPOSE_2
                TYPE(CSR_MATRIX) A
                
                A = CREATE_CSR_MATRIX(3,3)

                CALL SET(A, 1, 1, 2.D0)
                CALL SET(A, 1, 2, -1.D0)
                CALL SET(A, 1, 3, -1.D0)
                CALL SET(A, 2, 1, 2.D0)
                CALL SET(A, 2, 2, -1.D0)
                CALL SET(A, 2, 3, -1.D0)
                CALL SET(A, 3, 1, -7.D0)
                CALL SET(A, 3, 2, 5.D0)
                CALL SET(A, 3, 3, 1.D0)
                
                TEST_TRANSPOSE_2 = CSR_MATRIX_EQUAL(A, TR(TR(A)))
                
                CALL FREE_CSR_MATRIX(A)
        END FUNCTION
        
        FUNCTION TEST_INNER_PRODUCT_1() RESULT(B)
                LOGICAL :: B
                
                TYPE(CSR_MATRIX) :: A
                TYPE(CSR_MATRIX) :: C
                
                A = CREATE_CSR_MATRIX(5, 1)
                C = CREATE_CSR_MATRIX(5, 1)
                
                CALL A%SET(1, 1, 1.D0)
                CALL A%SET(1, 2, 2.D0)
                CALL A%SET(1, 3, -5.D0)
                CALL A%SET(1, 5, 5.D0)
                
                CALL C%SET(1, 1, 2.D0)
                CALL C%SET(1, 4, 2.D0)
                CALL C%SET(1, 5, 1.D0)
                
                B = CSR_MATRIX_INNER_PRODUCT(A, C) == GET(CSR_MATRIX_PROD(TR(A), C), 1, 1)
        END FUNCTION
        
        FUNCTION TEST_RECOPY() RESULT(B)
                LOGICAL :: B
                
                TYPE(CSR_MATRIX) :: A
                TYPE(CSR_MATRIX) :: C
                
                A = CREATE_CSR_MATRIX(3,3)

                CALL SET(A, 1, 1, 2.D0)
                CALL SET(A, 1, 3, -1.D0)
                CALL SET(A, 1, 2, -1.D0)
                CALL SET(A, 2, 1, 2.D0)
                CALL SET(A, 2, 2, -1.D0)
                CALL SET(A, 2, 3, -1.D0)
                CALL SET(A, 3, 1, -7.D0)
                CALL SET(A, 3, 2, 5.D0)
                CALL SET(A, 3, 3, 1.D0)
                
                C = CREATE_CSR_MATRIX_RECOPY(A)
                
                B = CSR_MATRIX_EQUAL(A, C)
        END FUNCTION
        
        END MODULE UNIT_TEST
