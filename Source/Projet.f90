!> 
!! @mainpage  Utilisation de matrice CSR pour la résolution de AX = B
!! \section Présentation
!! Le module CSR implémente un stockage efficient, en terme de mémoire, pour des matrices creuses.
!! Dans un premier temps, on utilise la décomposition LU pour résoudre le problème.
!! Dans un second temps, on utilise la méthode du gradient conjugué pour résoudre le problème.
!! 
!! \section CSR
!! \snippet CSR.f90 CSR
!! \snippet CSR.f90 CSR1
!! L'implémentation du mode de stockage CSR est fait de tel sorte que rowIdx commence par un 0 et que les colonnes sont indéxés à partir de 1.
!! Pour utiliser des matrices avec des matrices ayant des colonnes indéxés à partir de 0, il suffit juste de rajouter 1 à m_col. Ceci peut être mis en place très facilement grâce à une subfunction.
!! 
!! \section LU
!! On utilise la décomposition LU de la matrice A pour obtenir la solution. Cela revient, alors, à résoudre deux système triangulaire.
!! 
!! \section sec Gradient Conjugué
!! La méthode du gradient conjugué peut être utilisé de deux façon. Soit de manière complète, dans ce cas, l'intégralité des calculs sont fait pour obtenir la solution, soit de manière à rechercher une solution avec une précision donnée, dans ce cas, quand la précision de la solution est suffisante la fonction la retourne directement.
!! 
!! \section Solver
!! Ce module regroupe les appelles pour résoudre un système \f$AX=B\f$. Les deux méthodes proposés sont LU et gradient conjugué.
!! 
!! \section TEST
!! Ce module implémente quelque test pour vérifier le fonctionnement basique de l'ensemble des modules. Ceci permet d'évaluer rapidement si l'ajout ou modification de fonction/subroutine aux différents modules modifie le comportement attendu.
!! 
!! \section Compilation
!! Pour compiler le code executer le script *compile.sh*.
!! 
!! \section Exemple
!! On importe, dans un premier temps les modules nécéssaires.
!! \snippet Projet.f90 Module
!! Dans cet exemple, on va créer une matrice **A** et un vecteur-colonne **B** charger depuis un fichier. On crée aussi un vecteur-colonne pour chaque solution X et Y.
!! \snippet Projet.f90 Create
!! On peut alors charger les matrices depuis les fichiers.
!! \snippet Projet.f90 CreateF
!! On va ensuite calculer la solution du sytème \f$ AX = B \f$ en utilisant les deux méthodes implémentés. X est la solution utilisant la décomposition LU tandis que Y est la solution utilisant la méthode du gradient conjugué.
!! \snippet Projet.f90 Solve
!! On peut alors, afficher le résultat des deux méthodes.
!! \snippet Projet.f90 PrintR
!! On peut alors, afficher l'erreur des deux méthodes.
!! \snippet Projet.f90 PrintE
!! On sauvegarde alors la matrice A dans le fichier **svg1**.
!! \snippet Projet.f90 Svg
!! Il est possible d'initialiser à la main les matrices **A** et **B**. On doit alors créer une matrice vide de dimension voulue puis initialiser les valeurs non nulles
!! \snippet Projet.f90 CreateAB
!! Il est aussi possible de tester le fonctionnement basique des modules.
!! \snippet Projet.f90 TEST
!! A la fin, on doit désallouer les matrices utilisés.
!! \snippet Projet.f90 Free

PROGRAM MAIN
        !! [Module]
        USE CSR
        USE LU
        USE CONJUGATE_GRADIENT
        USE SOLVER
        USE UNIT_TEST
        !! [Module]
        
        !! [Create]
        TYPE(CSR_MATRIX) A
        TYPE(CSR_MATRIX) B
        TYPE(CSR_MATRIX) X
        TYPE(CSR_MATRIX) Y
        TYPE(CSR_MATRIX) R
        !! [Create]
        
        INTEGER :: i,j
        
        !! [CreateAB]
        A = CREATE_CSR_MATRIX_SQUARE(3)
        B = CREATE_CSR_MATRIX(3, 1)
                
        CALL SET(A, 1, 1, 5.D0)
        CALL SET(A, 2, 1, 1.D0)
        CALL SET(A, 2, 2, 1.D0)
        CALL SET(A, 2, 3, 8.D0)
        CALL SET(A, 3, 1, 7.D0)
        CALL SET(A, 3, 2, 1.D0)
        CALL SET(A, 3, 3, 2.D0)
        CALL SET(A, 3, 2, 10.D0)
        
        CALL SET(B, 1, 1, 1.D0)
        CALL SET(B, 2, 1, 1.D0)
        CALL SET(B, 3, 1, 1.D0)
        !! [CreateAB]
        
        CALL FREE_CSR_MATRIX(A)
        CALL FREE_CSR_MATRIX(B)
        
        !! [CreateF]
        A = CREATE_CSR_MATRIX_FROM_FILE("A.csr")
        B = CREATE_CSR_MATRIX_FROM_FILE("B.csr")
        !! [CreateF]
        
        !! [Solve]
        X = SOLVER_LU(A, B)
        Y = SOLVER_GC(A, B)
        !! [Solve]
        
        WRITE(*,*) "Matrice Initiale"
        CALL A%DISPLAY
        CALL B%DISPLAY
        
        !! [PrintR]
        WRITE (*,*) "Solution"
        WRITE(*,*) "LU :"
        CALL X%DISPLAY
        WRITE(*,*) "GC :"
        CALL Y%DISPLAY
        !! [PrintR]
        
        !! [PrintE]
        WRITE (*,*) "Erreur"
        R = CSR_MATRIX_VECTOR_SUB(CSR_MATRIX_PROD(A,X), B)
        WRITE(*,*) "LU :", SQRT(CSR_MATRIX_INNER_PRODUCT(R,R))
        R = CSR_MATRIX_VECTOR_SUB(CSR_MATRIX_PROD(A,Y), B)
        WRITE(*,*) "GC :", SQRT(CSR_MATRIX_INNER_PRODUCT(R,R))
        !! [PrintE]
        
        !! [TEST]
        CALL TEST
        !! [TEST]
        
        !! [Svg]
        CALL A%CSR_MATRIX_SAVE_TO_FILE("A.csr")
        CALL B%CSR_MATRIX_SAVE_TO_FILE("B.csr")
        !! [Svg]
        
        !! [Free]
        CALL FREE_CSR_MATRIX(A)
        CALL FREE_CSR_MATRIX(B)
        CALL FREE_CSR_MATRIX(X)
        CALL FREE_CSR_MATRIX(Y)
        CALL FREE_CSR_MATRIX(R)
        !! [Free]
END PROGRAM MAIN





