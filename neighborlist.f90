MODULE nlist
    
    IMPLICIT NONE
    PUBLIC
    INTEGER :: tnAtoms

    REAL, DIMENSION(:,:), ALLOCATABLE :: pool_pos_car
    INTEGER, DIMENSION(:), ALLOCATABLE :: pool_ids
    INTEGER, DIMENSION(:), ALLOCATABLE :: supersymbols

    INTEGER, DIMENSION(:, :), ALLOCATABLE :: tneighs
    INTEGER, DIMENSION(:, :), ALLOCATABLE :: tneighs_incell
    INTEGER, DIMENSION(:), ALLOCATABLE :: tnum_neigh

    CONTAINS

!-----------------------------------------------------------------------------------!
! calc:  calculate the neighborlist based upon positions, pos, and cell, cell
!-----------------------------------------------------------------------------------!

    SUBROUTINE update_nlist(nAtoms, nelement, atomicNrs, uniqueNrs, pos_car, cell, max_rcut)
!f2py   INTENT(IN) :: nelement
!f2py   INTENT(IN) :: pos_car
!f2py   INTENT(IN) :: cell
!f2py   INTENT(IN) :: max_rcut
       ! input values
        INTEGER :: nAtoms, nelement
        REAL, DIMENSION(nAtoms,3) :: pos_car
        REAL, DIMENSION(3,3) :: cell
        REAL :: max_rcut
        INTEGER, DIMENSION(nelement) :: uniqueNrs
        INTEGER, DIMENSION(nAtoms) :: atomicNrs

        !Variables
        INTEGER :: i, j
        INTEGER, PARAMETER :: MAX_NEIGHS = 100
        INTEGER :: tncells
        INTEGER*4, DIMENSION(nAtoms) :: symbols
        INTEGER, DIMENSION(3) :: ncells
        REAL, DIMENSION(3, 3) :: supercell
        REAL, DIMENSION(nelement, nelement) :: rmins

        DO i = 1, nAtoms
            DO j = 1, nelement
                IF (atomicNrs(i) == uniqueNrs(j)) THEN 
                    symbols(i) = j
                END IF
            END DO
        END DO
        CALL calcCellNum(nAtoms, pos_car, cell, max_rcut, ncells)
        print *, 'calcCellNum done'
        
        tncells = ncells(1) * ncells(2) * ncells(3)
        tnAtoms = nAtoms * tncells
        ALLOCATE(pool_pos_car(tnAtoms,3))
        ALLOCATE(pool_ids(tnAtoms))
        ALLOCATE(supersymbols(tnAtoms))
        ALLOCATE(tneighs(tnAtoms, MAX_NEIGHS))
        ALLOCATE(tneighs_incell(tnAtoms, MAX_NEIGHS))
        ALLOCATE(tnum_neigh(tnAtoms))
    
        CALL genSupercell(nAtoms, pos_car, symbols, cell, ncells, tncells, &
        max_rcut, supercell)
        print *, 'genSupercell done'

        CALL calcNlist(tnAtoms, pool_pos_car, supercell, symbols,&
        max_rcut, nelement, rmins, tncells,& 
        tnum_neigh, tneighs, &
        tneighs_incell)
        print *, 'calcNlist done'

    END SUBROUTINE

    SUBROUTINE calcCellNum(nAtoms, pos_car, cell, rcut, ncells)
        INTEGER :: nAtoms
        REAL :: rcut
        REAL,DIMENSION(nAtoms,3) :: pos_car
        REAL,DIMENSION(3,3) :: cell
        
!! f2py   INTENT(IN) :: nAtoms
!! f2py   INTENT(IN) :: pos_car
!! f2py   INTENT(IN) :: cell
!! f2py   INTENT(IN) :: rcut
        !Variables
        INTEGER :: i, j, k
        INTEGER, DIMENSION(3) :: sort_ids
        REAL, DIMENSION(3) :: latt_const, orth_latt, cross_ab,cross_bc
        REAL :: tmp, cross_ab_norm, cross_bc_norm, sin_3

        !Output
        INTEGER, DIMENSION(3) :: ncells
!f2py   INTENT(OUT) :: ncells

        latt_const(1) = norm(cell(1, :))
        latt_const(2) = norm(cell(2, :))
        latt_const(3) = norm(cell(3, :))

        !Rank lattice constant
        sort_ids = (/1, 2, 3/)
        DO i =1, 3
          DO j = 1, 3-i
             IF (latt_const(j) > latt_const(j+1)) THEN
               k = sort_ids(j)
               sort_ids(j) = sort_ids(j+1)
               sort_ids(j+1) = k
               tmp = latt_const(j) 
               latt_const(j) = latt_const(j+1) 
               latt_const(j+1) = tmp
             END IF
          END DO
        END DO
        orth_latt(sort_ids(1)) = latt_const(1)
 
        cross_ab = cross_product(cell(sort_ids(1), :), cell(sort_ids(2),:))
        cross_ab_norm = norm(cross_ab)
        IF (cross_ab_norm ==0) THEN
           print*, 'Cell is not 3-D'
           CALL EXIT(1)
        END IF
        orth_latt(sort_ids(2)) = cross_ab_norm/latt_const(1) 
        cross_bc = cross_product(cell(sort_ids(3), :),cross_ab)
        cross_bc_norm = norm(cross_bc)
  
        sin_3 = cross_bc_norm/(cross_ab_norm*latt_const(3))
   
        IF (sin_3 <= 0.05) THEN   !sin(4_degree)=0.07
          orth_latt(sort_ids(3)) = latt_const(3)
        ELSE
          orth_latt(sort_ids(3)) = cross_bc_norm/cross_ab_norm
        END IF
        ncells = INT(2*rcut /orth_latt+0.5) + 1
    
    END SUBROUTINE

    SUBROUTINE genSupercell(nAtoms, pos_car, symbols, cell, ncells, tncells, rcut, supercell)

        INTEGER :: nAtoms, tncells
        REAL :: rcut
        REAL,DIMENSION(nAtoms,3) :: pos_car
        INTEGER, DIMENSION(nAtoms) :: symbols
        REAL,DIMENSION(3,3) :: cell
        
        !Variables
        INTEGER :: i, j, k, i_offset, id
        INTEGER, DIMENSION(3) :: ncells
        INTEGER,DIMENSION(tncells,3) :: offset_list  !
        REAL,DIMENSION(tncells,3) :: offset_vector
        
        !Output
        REAL,DIMENSION(3,3) :: supercell

        i_offset = 1
        DO i = 1, ncells(1)
        DO j = 1, ncells(2)
            DO k = 1, ncells(3)
                offset_list(i_offset, :) = (/i-1, j-1, k-1/)
                i_offset = i_offset + 1
            END DO
        END DO
        END DO
        offset_vector = MATMUL(offset_list, cell)
        supercell(1,:) = ncells(1) * cell(1,:)
        supercell(2,:) = ncells(2) * cell(2,:)
        supercell(3,:) = ncells(3) * cell(3,:)
        
        supersymbols = 0
        DO i = 1, tncells
           DO j = 1, nAtoms
                id = (i-1)*nAtoms+j
                pool_pos_car(id,:) = offset_vector(i,:) + pos_car(j,:)
                pool_ids(id) = j
                supersymbols(id) = symbols(j)
           END DO
        END DO

    END SUBROUTINE

    SUBROUTINE calcNlist(nAtoms, pos_car, cell, symbols, rcut, nelement, &
                         rmins, tncells, num_neigh, neighs, neighs_incell)
        IMPLICIT NONE
        REAL,PARAMETER :: PI = 4.*ATAN(1.0_8)
        REAL,PARAMETER :: RAD2DEG = 180.0_8/PI

        ! input values
        INTEGER :: nAtoms, nelement, tncells
        REAL :: rcut
        REAL,DIMENSION(nAtoms,3) :: pos_car
        REAL,DIMENSION(3,3) :: cell
        INTEGER, DIMENSION(nAtoms) :: symbols

        ! variables
        INTEGER, PARAMETER :: MAX_NEIGHS = 100
        INTEGER, DIMENSION(nAtoms) :: ghost_indices
        REAL,DIMENSION(nAtoms,3) :: pos_dir
        REAL,DIMENSION(3,3) :: car2dir, dir2car
        REAL,DIMENSION(3) :: v_dir, v_car, v_unit
        REAL,DIMENSION(3) :: v_car_o, v_dir_o
        REAL :: dist, dist2, rcut2, dist2_o
        INTEGER :: i, j, check, natoms_percell, n_ghost

        !output valeus
        INTEGER, DIMENSION(nAtoms) :: num_neigh
        INTEGER,DIMENSION(nAtoms,MAX_NEIGHS) :: neighs
        INTEGER,DIMENSION(nAtoms,MAX_NEIGHS) :: neighs_incell
        REAL, DIMENSION(nelement, nelement) :: rmins

        rcut2 = rcut**2  ! should check the rcut is defined
        natoms_percell = nAtoms/tncells
        dir2car = TRANSPOSE(cell)
        car2dir = inverse(dir2car)
        ! to direct coordinates
        DO i = 1, nAtoms
            pos_dir(i,:) = MATMUL(car2dir, pos_car(i,:))
        END DO

        ! pair terms
        num_neigh = 0
        neighs = 0

        DO i = 1, nAtoms        
            DO j = i+1, nAtoms
                v_dir = pos_dir(j,:) - pos_dir(i,:)
                v_dir_o = v_dir
                v_dir = MOD(v_dir + 100.5, 1.0) - 0.5
                v_car = MATMUL(dir2car, v_dir)
                dist2 = SUM(v_car*v_car)
                IF (dist2<rcut2) THEN
                    num_neigh(i) = num_neigh(i) + 1
                    num_neigh(j) = num_neigh(j) + 1
                    dist = SQRT(dist2)
                    v_unit = v_car/dist
                    !Record minimum distance of two atoms
                    IF (dist<rmins(symbols(pool_ids(i)), symbols(pool_ids(j)))) THEN
                        rmins(symbols(pool_ids(i)), symbols(pool_ids(j))) = dist
                        rmins(symbols(pool_ids(j)), symbols(pool_ids(i))) = dist
                    END IF

                    ! Lei added
                    ! make neighs to record the index in pairs not origin structure
                    neighs(i, num_neigh(i)) = j
                    neighs(j, num_neigh(j)) = i
                    neighs_incell(i, num_neigh(i)) = pool_ids(j)
                    neighs_incell(j, num_neigh(j)) = pool_ids(i)

                END IF
            END DO
        END DO
    END SUBROUTINE

!-----------------------------------------------------------------------------------!
! norm:  return the norm of A(3) vector
!-----------------------------------------------------------------------------------!

    FUNCTION norm(A)

        REAL,DIMENSION(3) :: A
        REAL :: norm
        norm = SQRT(SUM(A * A))
        RETURN

    END FUNCTION

!-----------------------------------------------------------------------------------!
! cross_product: return the cross product of a(3) and b(3)
!-----------------------------------------------------------------------------------!

    FUNCTION cross_product(a,b)

        REAL,DIMENSION(3) :: a,b
        REAL,DIMENSION(3) :: cross_product

        cross_product(1) = a(2)*b(3) - a(3)*b(2)
        cross_product(2) = a(3)*b(1) - a(1)*b(3)
        cross_product(3) = a(1)*b(2) - a(2)*b(1)

        RETURN
    END FUNCTION cross_product

!-----------------------------------------------------------------------------------!
! inverse:  return the inverse of A(3,3)
!-----------------------------------------------------------------------------------!

    FUNCTION inverse(A)

        REAL,DIMENSION(3,3) :: A
        REAL,DIMENSION(3,3) :: inverse
        REAL :: det

        det = determinant(A)
        IF (det == 0) STOP 'Divide by zero in matrix inverse'
        inverse = adjoint(A)/det

        RETURN
    END FUNCTION inverse

!-----------------------------------------------------------------------------------!
! adjoint: return the adjoint of A(3,3)
!-----------------------------------------------------------------------------------!

    FUNCTION adjoint(A)

        REAL,DIMENSION(3,3) :: A
        REAL,DIMENSION(3,3) :: adjoint

        adjoint(1,1) = A(2,2)*A(3,3) - A(3,2)*A(2,3)
        adjoint(1,2) = A(1,3)*A(3,2) - A(1,2)*A(3,3)
        adjoint(1,3) = A(1,2)*A(2,3) - A(2,2)*A(1,3)

        adjoint(2,1) = A(2,3)*A(3,1) - A(2,1)*A(3,3)
        adjoint(2,2) = A(1,1)*A(3,3) - A(1,3)*A(3,1)
        adjoint(2,3) = A(1,3)*A(2,1) - A(1,1)*A(2,3)

        adjoint(3,1) = A(2,1)*A(3,2) - A(2,2)*A(3,1)
        adjoint(3,2) = A(1,2)*A(3,1) - A(1,1)*A(3,2)
        adjoint(3,3) = A(1,1)*A(2,2) - A(1,2)*A(2,1)

        RETURN
    END FUNCTION adjoint

!-----------------------------------------------------------------------------------!
! determinant: of a 3x3 matrix 
!-----------------------------------------------------------------------------------!

    FUNCTION determinant(A)

        REAL,DIMENSION(3,3) :: A
        REAL :: determinant

        determinant = A(1,1)*A(2,2)*A(3,3) &
                    - A(1,1)*A(2,3)*A(3,2) &
                    - A(1,2)*A(2,1)*A(3,3) &
                    + A(1,2)*A(2,3)*A(3,1) &
                    + A(1,3)*A(2,1)*A(3,2) &
                    - A(1,3)*A(2,2)*A(3,1)

        RETURN
    END FUNCTION

END MODULE
