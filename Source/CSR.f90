! MODULE: CSR_MATRIX
!> @author
!> Mechineau Alexandre
! DESCRIPTION:
!> Ce module permet de créer,stocker et manipuler une matrice sous la forme CSR. Plusieurs opérations courantes sont implémentés, avec certains rafinement pour le cas de vecteur.

MODULE CSR
IMPLICIT NONE

PRIVATE :: INCREASE_ROWIDX, DECREASE_ROWIDX, FIND_BEGIN_ROW, FIND_SIZE_ROW
PRIVATE :: DELETE_IN_CSR, REPLACE_IN_CSR, INSERT_IN_CSR

!> Type décrivant une matrice CSR
!! Attention m_col stocke l'index des colonnes en commencant par 1.
!! [CSR]
TYPE CSR_MATRIX
        !> Tableau stockant l'index des lignes
        INTEGER , DIMENSION(:), ALLOCATABLE     :: m_rowIdx
        !> Tableau stockant les colonnes non nulles
        INTEGER , DIMENSION(:), ALLOCATABLE     :: m_col
        !> Tableau stockant les valeurs non nulles
        REAL(KIND=selected_real_kind(15, 307)) , DIMENSION(:), ALLOCATABLE      :: m_val
        
        !> Nombre de colonnes de la matrice
        INTEGER                                 :: m_m ! Number of row
        !> Nombre de lignes de la matrice
        INTEGER                                 :: m_n ! Number of column
        !! [CSR]
        
        
        
        CONTAINS
        PROCEDURE       ::  IS_COLUMN_VECTOR, TR, CSR_MATRIX_SCALAR_PROD, &
                                   IS_IN_MATRIX, GET, GET_COLUMN_INDEX, &
                                   GET_COLUMN_VALUE, SET, &
                                   DISPLAY, DISPLAY_DEBUG, &
                                   CSR_MATRIX_ADDITIVE_INVERSE,&
                                    CSR_MATRIX_SAVE_TO_FILE
!! [CSR1]
END TYPE CSR_MATRIX
!! [CSR1]

CONTAINS
! Method to create/destroy a CSR_MATRIX 
        !> Cette fonction permet de créer une matrice depuis un fichier.
        !! La matrice dans le fichier doit utilisé le style C-CSR modifié pour que les colonnes commence par être indexé par 1.
        !! **On admet que le fichier passé est valide**
        !! \param f Fichier stockant une matrice CSR
        !! \returns Matrice CSR créé depuis le fichier
        
        !! Struct of CSR file
        !! m_m
        !! m_n
        !! m_rowIdx 
        !! SIZE(m_col)
        !! m_col
        !! m_val
        
        FUNCTION CREATE_CSR_MATRIX_FROM_FILE(f) RESULT(M)
                character(len=*), intent(in)    :: f
                
                TYPE (CSR_MATRIX) M
                
                INTEGER, PARAMETER      :: out_unit = 20
                INTEGER                 :: tmp
                
                
                open (unit = out_unit, file = f)
                
                READ(out_unit,*) M%m_m
                READ(out_unit,*) M%m_n
                
                ALLOCATE( M%m_rowIdx(M%m_m + 1) )
                READ(out_unit,*) M%m_rowIdx
                
                READ(out_unit,*) tmp
                ALLOCATE( M%m_col(tmp) )
                ALLOCATE( M%m_val(tmp) )
                
                READ(out_unit,*) M%m_col
                READ(out_unit,*) M%m_val
                CLOSE(out_unit)
        END FUNCTION CREATE_CSR_MATRIX_FROM_FILE
        
        !> Cette subroutine permet de sauvegarder une matrice CSR dans un fichier.
        !! \param M Matrice CSR à sauver dans un fichier
        !! \param f Fichier de sortier
        SUBROUTINE CSR_MATRIX_SAVE_TO_FILE(M, f)
                CLASS (CSR_MATRIX) M
                character(len=*), intent(in)    :: f
                
                INTEGER, PARAMETER      :: out_unit = 20
                INTEGER                 :: i
                
                
                open (unit=out_unit,file=f,action="write",status="replace")
                
                WRITE(out_unit,*) M%m_m
                WRITE(out_unit,*) M%m_n
                
                WRITE(out_unit,*) (M%m_rowIdx(i), i=1,SIZE(M%m_rowIdx))
                
                WRITE(out_unit,*) SIZE(M%m_col)
                
                WRITE(out_unit,*) (M%m_col(i), i=1,SIZE(M%m_col))
                WRITE(out_unit,*) (M%m_val(i), i=1,SIZE(M%m_col))
                
                CLOSE (out_unit)
        END SUBROUTINE CSR_MATRIX_SAVE_TO_FILE
        
        !> Cette fonction permet de créer une matrice carrée vide
        !! \param n Taille d'un coté d'une matrice carrée
        !! \returns Matrice carrée de taille \f$n*n\f$
        FUNCTION CREATE_CSR_MATRIX_SQUARE(n)
                INTEGER                         :: n
                TYPE (CSR_MATRIX)       :: CREATE_CSR_MATRIX_SQUARE
                
                CREATE_CSR_MATRIX_SQUARE = CREATE_CSR_MATRIX(n, n)
        END FUNCTION CREATE_CSR_MATRIX_SQUARE
        
        !> Cette fonction permet de dupliquer une matrice
        !! \param M Une matrice CSR que l'on souhaite recopier
        !! \returns Une nouvelle matrice CSR avec les mêmes valeurs que la matrice passée en argument
        FUNCTION CREATE_CSR_MATRIX_RECOPY(M) RESULT(N)
                TYPE (CSR_MATRIX), intent(in)   :: M
                TYPE (CSR_MATRIX)               :: N
                
                N%m_m = M%m_m
                N%m_n = M%m_n
                
                ! Allocate
                ALLOCATE( N%m_rowIdx(M%m_m + 1) )
                ALLOCATE( N%m_col(SIZE(M%m_col)) )
                ALLOCATE( N%m_val(SIZE(M%m_col)) )
                ! Recopy
                N%m_rowIdx(1:SIZE(N%m_rowIdx)) = M%m_rowIdx(1:SIZE(N%m_rowIdx))
                N%m_col( 1:SIZE(N%m_col) ) = M%m_col( 1:SIZE(N%m_col) )
                N%m_val( 1:SIZE(N%m_val) ) = M%m_val( 1:SIZE(N%m_val) )
        END FUNCTION CREATE_CSR_MATRIX_RECOPY
        
        !> Cette fonction crée une matrice de taille quelconque
        !! \param m Nombre de ligne de la matrice
        !! \param n Nombre de colonne de la matrice
        !! \returns Matrice CSR de taille \f$ m*n\f$
        FUNCTION CREATE_CSR_MATRIX(m, n)
                INTEGER, intent(in)                     :: m
                INTEGER, intent(in)                     :: n
                TYPE (CSR_MATRIX)                       :: CREATE_CSR_MATRIX
                
                CREATE_CSR_MATRIX%m_m = m
                CREATE_CSR_MATRIX%m_n = n
                
                ALLOCATE( CREATE_CSR_MATRIX%m_rowIdx(m + 1) )
                !Initialize to zero
                CREATE_CSR_MATRIX%m_rowIdx = 0
                ALLOCATE( CREATE_CSR_MATRIX%m_col(0) ) 
                ALLOCATE( CREATE_CSR_MATRIX%m_val(0) )
        END FUNCTION
        
        !> Cette subroutine permet de détruire une matrice CSR
        !! \param M Matrice à désallouer
        SUBROUTINE FREE_CSR_MATRIX(M)
                TYPE(CSR_MATRIX) ::  M
                
                M%m_n   = 0
                M%m_m   = 0
                
                DEALLOCATE( M%m_rowIdx )
                DEALLOCATE( M%m_col )
                DEALLOCATE( M%m_val )
        END SUBROUTINE FREE_CSR_MATRIX
        
        
! Function/Subroutine for the CSR_MATRIX

        !> Cette fonction permet de vérifier que des indices sont bien valide pour une matrice donnée
        !! \param M Matrice CSR
        !! \param row Ligne 
        !! \param col Colonne
        !! \returns Retourne VRAI si les indices \f$ (row,col)\f$ sont valides pour la matrice donnée
        FUNCTION IS_IN_MATRIX(M,row, col)
                CLASS(CSR_MATRIX), intent(in)           :: M
                INTEGER, intent(in)                     :: row
                INTEGER, intent(in)                     :: col
                
                LOGICAL :: IS_IN_MATRIX
                
                IS_IN_MATRIX = row>=1 .AND. row<=M%m_m .AND. col>=1 &
                                                        .AND. col<=M%m_n
        END FUNCTION IS_IN_MATRIX
        
        !> Cette fonction permet de  vérifier si une matrice est un vecteur colonne
        !! \param V Une matrice CSR
        !! \returns Retourne VRAI si la matrice est un vecteur colonne
        FUNCTION IS_COLUMN_VECTOR(V) RESULT(B)
                CLASS(CSR_MATRIX), intent(in)   :: V
                
                LOGICAL :: B
                
                IF(V%m_m>0 .AND. V%m_n==1) THEN
                        B = .TRUE.
                ELSE
                        B = .FALSE.
                ENDIF
        END FUNCTION IS_COLUMN_VECTOR
        
        !> Cette fonction permet d'obtnir la valeur d'un élément d'une matrice
        !! \param M Une matrice CSR
        !! \param row Ligne de la valeur souhaitée
        !! \param col Colonne de la valeur souhaitée
        !! \return Valeur contenu à la position \f$ (row,col) \f$. Si La position est invalide retourne 0.
        FUNCTION GET(M, row, col)
                CLASS(CSR_MATRIX), intent(in)           ::  M
                INTEGER, intent(in)                     :: row
                INTEGER, intent(in)                     :: col
                
                REAL(KIND=selected_real_kind(15, 307))                  :: GET
                
                INTEGER                         :: val_beg
                INTEGER                 :: val_end
                INTEGER                         :: i

                ! Compute the number of element before the the current row and the next
                val_beg = M%m_rowIdx(row)
                val_end = M%m_rowIdx(row + 1)
                
                ! 0 value in this row 
                IF ( val_beg == val_end ) THEN 
                        
                        GET = 0.
                        RETURN
                ENDIF
                ! There is at least one value
                DO i = val_beg+1, val_end, 1
                        IF ( M%m_col(i) == col ) THEN
                                GET = M%m_val(i)
                                RETURN
                        ENDIF
                ENDDO
                GET = 0.
        END FUNCTION GET

        !> Cette subroutine permet de modifier une valeur dans une matrice CSR
        !! \param M Une matrice CSR
        !! \param row Ligne de la valeur souhaitée
        !! \param col Colonne de la valeur souhaitée
        !! \param val Valeur souhaitée
        SUBROUTINE SET(M, row, col, val)
                CLASS(CSR_MATRIX)       :: M
                INTEGER, intent(in)     :: row
                INTEGER, intent(in)     :: col
                REAL(KIND=selected_real_kind(15, 307)), intent(in)      :: val 
                
                INTEGER :: beg
                INTEGER :: size_row
                INTEGER :: i

                IF( IS_IN_MATRIX(M, row, col) .eqv. .FALSE.) THEN
                        RETURN
                ENDIF
                
                beg = FIND_BEGIN_ROW(M, row)
                size_row = FIND_SIZE_ROW(M, row)
                
                IF(size_row == 0) THEN
                        IF (val == 0.) THEN
                                RETURN
                        ENDIF
                        CALL INSERT_IN_CSR(M, beg, col, val)
                        CALL INCREASE_ROWIDX(M, row)
                        RETURN
                ENDIF
                
                !Search position in the row to do the insertion
                DO i = beg, (beg + size_row - 1), 1
                        IF (M%m_col(i)<col) THEN
                                !Do nothing
                        ELSE
                                IF (M%m_col(i)==col) THEN
                                        !Replace
                                        IF (val == 0.) THEN
                                                CALL DELETE_IN_CSR(M,i)
                                                CALL DECREASE_ROWIDX(M, row)
                                                RETURN
                                        ENDIF
                                        CALL REPLACE_IN_CSR(M,i,val)
                                        RETURN
                                ELSE
                                        IF (val == 0.) THEN
                                                RETURN
                                        ENDIF
                                        CALL INSERT_IN_CSR(M,i,col,val)
                                        CALL INCREASE_ROWIDX(M, row)
                                        RETURN
                                ENDIF
                        ENDIF
                ENDDO
                ! Insertion is at the end.
                i = beg + size_row
                IF (val == 0.) THEN
                        RETURN
                ENDIF
                CALL INSERT_IN_CSR(M,i,col,val)
                CALL INCREASE_ROWIDX(M, row)
        END SUBROUTINE SET
        
        !> Cette subroutine augmente le nombre d'élément compté à partir d'une ligne donnée
        !! \param M Une matrice CSR
        !! \param row Index de la ligne
        SUBROUTINE INCREASE_ROWIDX(M, row)
                TYPE(CSR_MATRIX)        :: M
                INTEGER, intent(in)     :: row
                
                INTEGER :: i
                
                DO i=row+1, SIZE(M%m_rowIdx), 1
                        M%m_rowIdx(i) = M%m_rowIdx(i) + 1
                ENDDO
        END SUBROUTINE INCREASE_ROWIDX
        
        !> Cette subroutine réduit le nombre d'élément compté à partir d'une ligne donnée
        !! \param M Une matrice CSR
        !! \param row Index de la ligne
        SUBROUTINE DECREASE_ROWIDX(M, row)
                TYPE(CSR_MATRIX)        :: M
                INTEGER, intent(in)     :: row
                
                INTEGER :: i
                
                DO i=row+1, SIZE(M%m_rowIdx), 1
                        M%m_rowIdx(i) = M%m_rowIdx(i) - 1
                ENDDO
        END SUBROUTINE DECREASE_ROWIDX
        
        !> Cette fonction permet de trouver l'indice de début d'une ligne donnée dans m_col et m_val
        !! \param M Une matrice CSR
        !! \param row Index de la ligne
        !! \return Position du début de la ligne donnée en paramètre dans M%m_col et M%m_val
        FUNCTION FIND_BEGIN_ROW(M, row)
                TYPE(CSR_MATRIX)        :: M
                INTEGER, intent(in)     :: row
                
                INTEGER :: FIND_BEGIN_ROW
                
                FIND_BEGIN_ROW = M%m_rowIdx(row) + 1
        END FUNCTION FIND_BEGIN_ROW
        
        !> !> Cette fonction permet de trouver le nombre d'élement pour une ligne donnée
        !! \param M Une matrice CSR
        !! \param row Index de la ligne
        !! \return Nombre d'élement présent sur une ligne donnée
        FUNCTION FIND_SIZE_ROW(M, row)
                TYPE(CSR_MATRIX)        :: M
                INTEGER, intent(in)     :: row
                
                INTEGER :: FIND_SIZE_ROW
                
                FIND_SIZE_ROW = M%m_rowIdx(row + 1) - M%m_rowIdx(row)
        END FUNCTION FIND_SIZE_ROW
        
        !> Cette subroutine permet de supprimer un élement de m_col et m_val
        !! \param M Une matrice CSR
        !! \param i Position à supprimer
        SUBROUTINE DELETE_IN_CSR(M, i)
                TYPE(CSR_MATRIX)        :: M
                INTEGER, intent(in)     :: i
                
                INTEGER , DIMENSION(:), ALLOCATABLE     :: tmp
                REAL(KIND=selected_real_kind(15, 307)) , DIMENSION(:), ALLOCATABLE      :: tmp1
                
                ALLOCATE(tmp( SIZE(M%m_col) - 1 ))
                ALLOCATE(tmp1( SIZE(M%m_val) - 1 ))
                
                tmp  = 0
                tmp1 = 0.
                
                tmp( 1:(i-1) )  = M%m_col( 1:(i-1) )
                tmp1( 1:(i-1) ) = M%m_val( 1:(i-1) )
                
                
                tmp( (i):SIZE(tmp) )  = M%m_col( (i+1):SIZE(M%m_col) )
                tmp1( (i):SIZE(tmp) ) = M%m_col( (i+1):SIZE(M%m_col) )
                
                CALL move_alloc(tmp, M%m_col)
                CALL move_alloc(tmp1, M%m_val)
        END SUBROUTINE DELETE_IN_CSR
        
        !> Cette subroutine permet de remplacer la valeur d'un élément non nul déja présent dans une matrice CSR
        !! \param M Une matrice CSR
        !! \param i Position à remplacer
        !! \param val Valeur à remplacer
        SUBROUTINE REPLACE_IN_CSR(M, i, val)
                TYPE(CSR_MATRIX)        :: M
                INTEGER, intent(in)     :: i
                REAL(KIND=selected_real_kind(15, 307)), intent(in)      :: val
                
                M%m_val(i) = val
        END SUBROUTINE REPLACE_IN_CSR
        
        !> Cette subroutine insert un élémet dans une matrice CSR
        !! \param M Une matrice CSR
        !! \param i Position à insérer
        !! \param col Valeur à inserer dans m_col
        !! \param val Valeur à inserer dans m_val
        SUBROUTINE INSERT_IN_CSR(M, i, col, val)
                TYPE(CSR_MATRIX)        :: M
                INTEGER, intent(in)     :: i
                INTEGER, intent(in)     :: col
                REAL(KIND=selected_real_kind(15, 307)), intent(in)      :: val
                
                INTEGER , DIMENSION(:), ALLOCATABLE     :: tmp
                REAL(KIND=selected_real_kind(15, 307)) , DIMENSION(:), ALLOCATABLE      :: tmp1
                
                ALLOCATE(tmp( SIZE(M%m_col) + 1 ))
                ALLOCATE(tmp1( SIZE(M%m_val) + 1 ))
                
                tmp  = 0
                tmp1 = 0.
                
                tmp( 1:(i-1) )  = M%m_col( 1:(i-1) )
                tmp1( 1:(i-1) ) = M%m_val( 1:(i-1) )
                
                
                IF (i<SIZE(tmp)) THEN
                        tmp( (i+1):SIZE(tmp) )  = M%m_col( i:SIZE(M%m_col) )
                        tmp1( (i+1):SIZE(tmp) ) = M%m_val( i:SIZE(M%m_col) )
                ENDIF
                
                tmp(i)  = col
                tmp1(i) = val
                
                CALL MOVE_ALLOC(tmp, M%m_col)
                CALL MOVE_ALLOC(tmp1, M%m_val)
        END SUBROUTINE INSERT_IN_CSR
        
        !> Cette subroutine affiche une matrice CSR avec des informations avancées
        !! \param M Matrice à afficher
        SUBROUTINE DISPLAY_DEBUG(M)
                CLASS(CSR_MATRIX), intent(in) :: M
                
                ! Raw access to CSR_MATRIX
                WRITE(*, *) "DEBUG INFO"
                WRITE(*, *) "m_m  ", M%m_m
                WRITE(*, *) "m_n  ", M%m_n
                WRITE(*, *) "rowIdx", M%m_rowIdx
                WRITE(*, *) "Col    ", M%m_col
                WRITE(*, *) "val     ", M%m_val
                CALL M%DISPLAY()
        END SUBROUTINE DISPLAY_DEBUG
        
        !> Cette subroutine affiche une matrice CSR
        !! \param M Matrice à afficher
        SUBROUTINE DISPLAY(M)
                CLASS(CSR_MATRIX), intent(in) :: M
                
                REAL(KIND=selected_real_kind(15, 307)) , DIMENSION(:), ALLOCATABLE      :: tmp_row 
                INTEGER                                 :: i
                INTEGER                                 :: j

                ALLOCATE(tmp_row(M%m_n))
                
                DO i = 1, M%m_m, 1
                        DO j = 1, M%m_n, 1
                                ! Use a buffer to display a line because WRITE return to line after each call 
                                ! Can be fixed using some weird parameters for WRITE
                                tmp_row(j) = GET(M, i, j)
                        ENDDO
                        WRITE(*,*) tmp_row
                ENDDO
                WRITE(*,*) ""
        END SUBROUTINE DISPLAY
        
        !Return data of the column of a row
        !> Cette fonction permet d'obtenir les colonnes non-nulles d'une ligne donnée
        !! \param M Matrice CSR
        !! \param row index d'une ligne de la matrice
        !! \returns Ensemble des colonnes à valeur non nulle
        FUNCTION GET_COLUMN_INDEX(M, row)
                CLASS(CSR_MATRIX), intent(in)   :: M
                INTEGER, intent(in)             :: row
                
                INTEGER , DIMENSION(:), ALLOCATABLE     :: GET_COLUMN_INDEX
                
                ! Compute the number of non-null element on the given row
                ! Then, allocate an array 
                ALLOCATE( GET_COLUMN_INDEX( M%m_rowIdx(row + 1) - M%m_rowIdx(row) ) )
                
                ! copy the value into the previous array
                GET_COLUMN_INDEX(1 : M%m_rowIdx(row + 1) - M%m_rowIdx(row)) &
                                                = M%m_col( M%m_rowIdx(row) + 1: M%m_rowIdx(row + 1) )
        END FUNCTION GET_COLUMN_INDEX
        
        !> Cette fonction permet d'obtenir les valeurs non-nulles d'une ligne donnée
        !! \param M Matrice CSR
        !! \param row index d'une ligne de la matrice
        !! \returns Ensemble des valeurs non nulle 
        FUNCTION GET_COLUMN_VALUE(M, row)
                CLASS(CSR_MATRIX), intent(in)   :: M
                INTEGER, intent(in)             :: row
                
                REAL(KIND=selected_real_kind(15, 307)) , DIMENSION(:), ALLOCATABLE      :: GET_COLUMN_VALUE
                
                ! Compute the number of non-null element on the given row
                ! Then, allocate an array 
                ALLOCATE( GET_COLUMN_VALUE( M%m_rowIdx(row + 1) - M%m_rowIdx(row) ) )
                
                ! copy the value into the previous array
                GET_COLUMN_VALUE(1 : M%m_rowIdx(row + 1) - M%m_rowIdx(row)) &
                                                = M%m_val( M%m_rowIdx(row) + 1: M%m_rowIdx(row + 1) )
        END FUNCTION GET_COLUMN_VALUE
        
        !> Cette fonction calcule le produit de deux matrices
        !! \param A Matrice CSR
        !! \param B Matrice CSR
        !! \returns \f$ A*B \f$
        FUNCTION CSR_MATRIX_PROD(A, B)
                TYPE(CSR_MATRIX), intent(in) :: A
                TYPE(CSR_MATRIX), intent(in) :: B
                
                !ASSERT compatible size
                TYPE(CSR_MATRIX) :: CSR_MATRIX_PROD
                
                INTEGER , DIMENSION(:), ALLOCATABLE     :: rowTmpIdx
                REAL(KIND=selected_real_kind(15, 307)) , DIMENSION(:), ALLOCATABLE      :: rowTmpVal
                REAL(KIND=selected_real_kind(15, 307))                                  :: cumSum
                INTEGER                                 :: i
                INTEGER                                 :: j
                INTEGER                                 :: k
                
                CSR_MATRIX_PROD = CREATE_CSR_MATRIX(A%m_m, B%m_n)
                
                DO i = 1, A%m_m, 1
                        !get a row: easy
                        rowTmpIdx = GET_COLUMN_INDEX(A, i)
                        rowTmpVal = GET_COLUMN_VALUE(A, i)
                        
                        DO j = 1, B%m_n, 1
                                cumSum = 0.
                                DO k = 1, SIZE(rowTmpIdx), 1
                                        cumSum = cumSum + rowTmpVal(k) * GET(B, rowTmpIdx(k), j)
                                ENDDO
                                IF (cumSum /= 0. ) THEN
                                        CALL SET(CSR_MATRIX_PROD, i, j, cumSum)
                                ENDIF
                        ENDDO
                ENDDO
        END FUNCTION CSR_MATRIX_PROD
        
        !> Cette fonction permet de vérifier si deux matrices sont égales
        !! \param A Matrice CSR
        !! \param B Matrice CSR
        !! \returns \f$ A=B \f$
        FUNCTION CSR_MATRIX_EQUAL(A, B)
                TYPE(CSR_MATRIX), intent(in)    :: A
                TYPE(CSR_MATRIX), intent(in)    :: B
                
                LOGICAL :: CSR_MATRIX_EQUAL
                
                INTEGER :: i
                INTEGER :: j
                
                CSR_MATRIX_EQUAL = .TRUE.
                
                IF(A%m_m/=B%m_m) THEN
                        CSR_MATRIX_EQUAL = .FALSE.
                        RETURN
                ELSE IF(A%m_n/=B%m_n) THEN
                        CSR_MATRIX_EQUAL = .FALSE.
                        RETURN
                ELSE
                        DO i = 1, A%m_m, 1
                                DO j = 1,A%m_n, 1
                                        IF(GET(A, i, j) /= GET(B, i, j)) THEN
                                                CSR_MATRIX_EQUAL = .FALSE.
                                                RETURN
                                        ENDIF
                                ENDDO
                        ENDDO
                ENDIF
        END FUNCTION CSR_MATRIX_EQUAL
        
        !> Cette fonction calcule la transposée d'une matrice
        !! \param M Matrice CSR
        !! \returns \f$ M^T \f$
        FUNCTION TR(M)
                ! transpose M
                CLASS(CSR_MATRIX), intent(in) :: M
                
                TYPE(CSR_MATRIX) :: TR
                
                INTEGER :: i
                INTEGER :: j
                
                TR = CREATE_CSR_MATRIX(M%m_n, M%m_m)
                
                DO i = 1, M%m_m, 1
                        DO j = 1, M%m_n, 1
                                CALL TR%SET(j, i, M%GET(i, j))
                        ENDDO
                ENDDO
        END FUNCTION TR
        
        !> Cette fonction calcule le produit scalaire de deux vecteurs
        !! \param A Matrice CSR définissant un vecteur colonne
        !! \param B Matrice CSR définissant un vecteur colonne
        !! \returns \f$ \langle A | B \rangle \f$
        FUNCTION CSR_MATRIX_INNER_PRODUCT(A, B) RESULT(X)
                TYPE(CSR_MATRIX), intent(in)    :: A
                TYPE(CSR_MATRIX), intent(in)    :: B
                
                REAL(KIND=selected_real_kind(15, 307)) :: X
                
                INTEGER :: i
                
                X = 0.
                IF(A%m_m == B%m_m .AND. A%m_n == B%m_n .AND. A%m_n == 1) THEN
                ELSE
                        RETURN
                ENDIF
                
                DO i = 2, SIZE(A%m_rowIdx), 1
                        IF (A%m_rowIdx(i) - A%m_rowIdx(i - 1) == 0) THEN
                                ! It's a zero :)
                        ELSE
                                IF (B%m_rowIdx(i) - B%m_rowIdx(i - 1) == 0) THEN
                                        ! It's a zero  again :)
                                ELSE
                                        X = X + A%m_val(A%m_rowIdx(i)) * B%m_val(B%m_rowIdx(i))
                                ENDIF
                        ENDIF
                ENDDO
        END FUNCTION CSR_MATRIX_INNER_PRODUCT
        
        !> Cette fonction calcule la somme de deux vecteurs colonnes
        !! \param A Matrice CSR définissant un vecteur colonne
        !! \param B Matrice CSR définissant un vecteur colonne
        !! \returns \f$ A + B \f$
        FUNCTION CSR_MATRIX_VECTOR_ADD(A , B) RESULT(C)
                TYPE(CSR_MATRIX), intent(in)    :: A
                TYPE(CSR_MATRIX), intent(in)    :: B
                
                TYPE(CSR_MATRIX) :: C
                
                INTEGER :: i
                
                C = CREATE_CSR_MATRIX(A%m_m, 1)
                
                IF(A%m_m == B%m_m .AND. A%m_n == B%m_n .AND. A%m_n == 1) THEN
                        ! GOOD SIZE
                ELSE
                        RETURN
                ENDIF
                
                DO i = 2, SIZE(A%m_rowIdx), 1
                        IF (A%m_rowIdx(i) - A%m_rowIdx(i - 1) == 0) THEN
                                ! It's a zero :)
                                ! But maybe it's not a zero in B
                                IF (B%m_rowIdx(i) - B%m_rowIdx(i - 1) == 0) THEN
                                        ! It's a zero  again :)
                                ELSE
                                        CALL C%SET(i-1, 1, B%m_val(B%m_rowIdx(i)))
                                ENDIF
                        ELSE
                                IF (B%m_rowIdx(i) - B%m_rowIdx(i - 1) == 0) THEN
                                        ! It's a zero Here:)
                                        CALL C%SET(i-1, 1, A%m_val(A%m_rowIdx(i)))
                                ELSE
                                        CALL C%SET(i-1, 1, A%m_val(A%m_rowIdx(i)) + B%m_val(B%m_rowIdx(i)))
                                ENDIF
                        ENDIF
                ENDDO
        END FUNCTION CSR_MATRIX_VECTOR_ADD
        
        !> Cette fonction calcule la différence de deux vecteurs colonnes
        !! \param A Matrice CSR définissant un vecteur colonne
        !! \param B Matrice CSR définissant un vecteur colonne
        !! \returns \f$ A-B \f$
        FUNCTION CSR_MATRIX_VECTOR_SUB(A , B) RESULT(C)
                TYPE(CSR_MATRIX), intent(in)    :: A
                TYPE(CSR_MATRIX), intent(in)    :: B
                
                TYPE(CSR_MATRIX) :: C
                
                INTEGER :: i
                
                C = CREATE_CSR_MATRIX(A%m_m, 1)
                
                IF(A%m_m == B%m_m .AND. A%m_n == B%m_n .AND. A%m_n == 1) THEN
                        ! GOOD SIZE
                ELSE
                        RETURN
                ENDIF
                
                DO i = 2, SIZE(A%m_rowIdx), 1
                        IF (A%m_rowIdx(i) - A%m_rowIdx(i - 1) == 0) THEN
                                ! It's a zero :)
                                ! But maybe it's not a zero in B
                                IF (B%m_rowIdx(i) - B%m_rowIdx(i - 1) == 0) THEN
                                        ! It's a zero  again :)
                                ELSE
                                        CALL C%SET(i - 1, 1, -B%m_val(B%m_rowIdx(i)))
                                ENDIF
                        ELSE
                                IF (B%m_rowIdx(i) - B%m_rowIdx(i - 1) == 0) THEN
                                        ! It's a zero Here:)
                                        CALL C%SET(i - 1, 1, A%m_val(A%m_rowIdx(i)))
                                ELSE
                                        CALL C%SET(i - 1, 1, A%m_val(A%m_rowIdx(i)) - B%m_val(B%m_rowIdx(i)))
                                ENDIF
                        ENDIF
                ENDDO
        END FUNCTION CSR_MATRIX_VECTOR_SUB
        
        !> Cette fonction calcule l'inverse associé à la loi additive d'une matrice
        !! \param A Matrice CSR
        !! \returns \f$ -A \f$
        FUNCTION CSR_MATRIX_ADDITIVE_INVERSE(A) RESULT(X)
                CLASS(CSR_MATRIX), intent(in)   :: A
                
                TYPE(CSR_MATRIX)        :: X
                
                INTEGER         :: i
                
                ! Direct copy of the matrix to a new one
                X = CREATE_CSR_MATRIX_RECOPY(A)
                ! Inverse all value of the matrix
                DO i = 1, SIZE(X%m_val), 1
                        X%m_val(i) = -X%m_val(i)
                ENDDO
        END FUNCTION CSR_MATRIX_ADDITIVE_INVERSE
        
        !> Cette fonction calcule le produit entre un réel et une matrice
        !! \param M Matrice CSR
        !! \param K un réel
        !! \returns \f$ k*M \f$
        FUNCTION CSR_MATRIX_SCALAR_PROD(M, K) RESULT(A)
                CLASS(CSR_MATRIX), intent(in)   :: M
                REAL(KIND=selected_real_kind(15, 307)), intent(in)              :: K
                
                TYPE(CSR_MATRIX)        :: A
                
                INTEGER                 :: i
                
                A%m_m = M%m_m
                A%m_n = M%m_n
                
                ALLOCATE(A%m_rowIdx(SIZE(M%m_rowIdx)))
                ALLOCATE(A%m_col(SIZE(M%m_col)))
                ALLOCATE(A%m_val(SIZE(M%m_col)))
                
                A%m_rowIdx(1:SIZE(M%m_rowIdx)) = M%m_rowIdx(1:SIZE(M%m_rowIdx))
                A%m_col(1:SIZE(M%m_col)) = M%m_col(1:SIZE(M%m_col))
                A%m_val(1:SIZE(M%m_col)) = K * M%m_val(1:SIZE(M%m_col))
        END FUNCTION CSR_MATRIX_SCALAR_PROD
        
END MODULE CSR
