! This project is licensed under the terms of the MIT license,
! read the LICENSE.txt file
! (c) 2018 Christian Henriksen


! subroutine NNCOPot computes an approximation of the interaction energy
! between NN' and CO .
! ARGUMENTS:
!   Let CoM_N be the center of mass of NN' and CoM_C the center of mass of CO.
! Then the arguments are --
! -- R (Bohr) the distance between CoM_N and CoM_C 
! -- Theta_A (degrees) the angle between the vector going from CoM_N to CoM_C
!       and the vector going from N to N' 
! -- Theta_B (degrees) the angle between the vector going from CoM_N to CoM_C
!       and the vector going from C to O 
! -- Phi (degrees) the dihedral angle N' CoM_N CoM_C O
! -- U (micro Hartrees) contains the interaction energy at the geometry given
!       by R, Theta_A, Theta_B and Phi

! Written in Fortran 95 .

subroutine NNCOPot(R, ThetaA, ThetaB, Phi, U)
    implicit none
    integer, parameter :: dp = 8 

    ! Declaration of subroutine parameters
    real(dp), intent(in) :: R, ThetaA, ThetaB, Phi
    real(dp), intent(out) :: U

    ! Declaration of constants that are global to the subroutine
    real(dp), parameter :: pi = 3.1415926535898d0
    integer, parameter :: num_fixed_atoms = 2
    integer, parameter :: num_moving_atoms = 2
    real(dp), parameter :: fixed_mol(num_fixed_atoms, 3) = &
        transpose(reshape ( &
             (/ -1.03716, 0.0, 0.0, 1.03716, 0.0, 0.0 /), &
             (/ 3, num_fixed_atoms /) ))
    real(dp), parameter :: moving_mol(num_moving_atoms, 3) = &
        transpose(reshape ( &
     (/ -1.218246, 0.0, 0.0, 0.913975, 0.0, 0.0 /), &
     (/ 3, num_fixed_atoms /) ))
    real(dp), parameter :: dil = 4.379907928794446 
    real(dp), parameter :: ref_dist =  8.0d0
    integer, parameter :: num_terms = 79
    integer, parameter :: sym_order = 2
    integer, parameter :: exponents(num_terms, sym_order, num_fixed_atoms, num_moving_atoms) = &
    reshape ( & 
      (/ & 
        0, 0, 0, 2, 0, 0, 1, 0, 1, 0, 3, 0, 1, 0, 0, 1, 1, 0, 0, 2, 0, 4, &
        0, 0, 2, 3, 0, 2, 1, 3, 1, 2, 2, 0, 0, 0, 5, 0, 1, 0, 3, 0, 1, 2, &
        0, 2, 1, 1, 3, 1, 2, 1, 2, 2, 3, 0, 0, 2, 0, 6, 0, 2, 0, 5, 2, 0, &
        4, 1, 0, 3, 0, 1, 1, 4, 1, 1, 1, 2, 3, 0, 0, 1, 0, 0, 1, 1, 0, 0, &
        2, 0, 0, 0, 2, 1, 1, 2, 0, 1, 0, 3, 0, 0, 1, 0, 0, 1, 0, 2, 1, 1, &
        1, 2, 0, 3, 4, 0, 0, 0, 2, 0, 0, 0, 0, 3, 0, 1, 2, 1, 1, 1, 3, 1, &
        2, 2, 0, 1, 0, 5, 0, 0, 0, 3, 0, 0, 3, 0, 0, 2, 0, 1, 1, 2, 1, 5, &
        1, 2, 4, 3, 0, 0, 1, 0, 0, 1, 1, 0, 0, 2, 0, 0, 0, 2, 1, 1, 2, 0, &
        1, 0, 3, 0, 0, 1, 0, 0, 1, 0, 2, 1, 1, 1, 2, 0, 3, 4, 0, 0, 0, 2, &
        0, 0, 0, 0, 3, 0, 1, 2, 1, 1, 1, 3, 1, 2, 2, 0, 1, 0, 5, 0, 0, 0, &
        3, 0, 0, 3, 0, 0, 2, 0, 1, 1, 2, 1, 5, 1, 2, 4, 3, 0, 0, 0, 2, 0, &
        0, 1, 0, 1, 0, 3, 0, 1, 0, 0, 1, 1, 0, 0, 2, 0, 4, 0, 0, 2, 3, 0, &
        2, 1, 3, 1, 2, 2, 0, 0, 0, 5, 0, 1, 0, 3, 0, 1, 2, 0, 2, 1, 1, 3, &
        1, 2, 1, 2, 2, 3, 0, 0, 2, 0, 6, 0, 2, 0, 5, 2, 0, 4, 1, 0, 3, 0, &
        1, 1, 4, 1, 1, 1, 2, 3, 1, 2, 0, 0, 1, 1, 0, 0, 2, 0, 0, 1, 1, 1, &
        2, 1, 0, 4, 0, 2, 0, 0, 3, 1, 1, 0, 2, 0, 0, 0, 1, 0, 0, 5, 0, 0, &
        0, 1, 3, 1, 1, 2, 2, 1, 2, 0, 3, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 4, &
        0, 0, 1, 3, 1, 0, 2, 2, 0, 2, 3, 0, 4, 4, 0, 1, 0, 2, 3, 0, 0, 0, &
        0, 1, 0, 1, 0, 0, 3, 0, 1, 0, 2, 1, 0, 0, 0, 0, 0, 3, 0, 1, 0, 1, &
        2, 1, 1, 1, 2, 1, 0, 1, 1, 0, 0, 2, 1, 0, 4, 1, 2, 1, 3, 2, 2, 0, &
        3, 0, 2, 0, 2, 1, 0, 2, 1, 0, 6, 5, 0, 1, 0, 5, 1, 2, 1, 2, 1, 2, &
        3, 1, 3, 1, 0, 3, 0, 0, 2, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 3, 0, 1, &
        0, 2, 1, 0, 0, 0, 0, 0, 3, 0, 1, 0, 1, 2, 1, 1, 1, 2, 1, 0, 1, 1, &
        0, 0, 2, 1, 0, 4, 1, 2, 1, 3, 2, 2, 0, 3, 0, 2, 0, 2, 1, 0, 2, 1, &
        0, 6, 5, 0, 1, 0, 5, 1, 2, 1, 2, 1, 2, 3, 1, 3, 1, 0, 3, 0, 0, 2, &
        0, 0, 0, 1, 2, 0, 0, 1, 1, 0, 0, 2, 0, 0, 1, 1, 1, 2, 1, 0, 4, 0, &
        2, 0, 0, 3, 1, 1, 0, 2, 0, 0, 0, 1, 0, 0, 5, 0, 0, 0, 1, 3, 1, 1, &
        2, 2, 1, 2, 0, 3, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 4, 0, 0, 1, 3, 1, &
        0, 2, 2, 0, 2, 3, 0, 4, 4, 0, 1, 0, 2, 3, 0, 0 &
      /), &
      (/  79, 2, 2, 2 /) &
    )
    real(dp), parameter :: coeffs(num_terms) = &
    (/ & 
        -32.119688495524052d0, &
        42713.232812234251d0, &
        -90922.334703143715d0, &
        45622.212649960755d0, &
        -43921.70928461336d0, &
        94044.828399882448d0, &
        -47238.943176444518d0, &
        -94238.96876316992d0, &
        222526.29231139383d0, &
        -151150.54159964123d0, &
        23820.699284501588d0, &
        99515.714637618454d0, &
        -151182.26225973005d0, &
        49958.187490857126d0, &
        -87836.876571440938d0, &
        113802.41696291082d0, &
        -27029.722171112306d0, &
        134244.42966897646d0, &
        -261135.26655306265d0, &
        270883.96863537468d0, &
        -309531.19025995804d0, &
        161306.1940339125d0, &
        -143021.61083325851d0, &
        -57524.723068940082d0, &
        27878.686752324153d0, &
        179528.31595040386d0, &
        344483.45010893181d0, &
        -226068.07973029485d0, &
        553579.66313241527d0, &
        -534808.10150174715d0, &
        -86760.136692276414d0, &
        -427291.73099021491d0, &
        377662.85735381918d0, &
        -38254.788107340122d0, &
        128635.64581783509d0, &
        -68503.925029206861d0, &
        -26919.286113998351d0, &
        -28092.37771039395d0, &
        332771.35646888334d0, &
        -243575.9915137917d0, &
        -56721.865322333128d0, &
        69469.778215854807d0, &
        -339359.59562560031d0, &
        346112.63491049834d0, &
        -93145.795373528032d0, &
        -33678.037329836938d0, &
        114684.9292733636d0, &
        -583925.1373044221d0, &
        469398.51013242028d0, &
        -189641.52695470807d0, &
        124011.67780773289d0, &
        86557.422330149304d0, &
        486416.22347594227d0, &
        -480243.41008052812d0, &
        20583.827190640601d0, &
        12990.256705130523d0, &
        -12025.493001092429d0, &
        6503.8933326982078d0, &
        -27982.019937435631d0, &
        27639.182998547651d0, &
        -12989.923393271689d0, &
        -62247.856388897999d0, &
        56617.245072513528d0, &
        6922.3801089256303d0, &
        92503.82184496909d0, &
        -76139.614727798093d0, &
        -4820.5730715874815d0, &
        -8329.7100928079162d0, &
        -75468.23164996953d0, &
        30006.874727631137d0, &
        20848.494313051575d0, &
        -7453.0700775478181d0, &
        34335.277550389961d0, &
        18736.812259763072d0, &
        -59868.759779880173d0, &
        44846.433409049117d0, &
        -38456.681133122474d0, &
        57833.681508431517d0, &
        -22747.360629756819d0  &
    /)

    ! declaration of local variables
    real(dp) :: dist_matrix(num_fixed_atoms, num_moving_atoms)
    real(dp) :: post_matrix(num_fixed_atoms, num_moving_atoms)
    integer :: sym_term_count
    integer :: term_count

    call lin_lin_to_distance_matrix(R, ThetaA, ThetaB, Phi, dist_matrix)
    post_matrix = exp( -(dist_matrix - ref_dist) / dil )

    ! run through all the symmetrized terms accumulating the result in result in U
    U = 0d0
    do sym_term_count = 1, num_terms
        do term_count = 1, sym_order
            U = U + coeffs(sym_term_count) / sym_order &
                * product(post_matrix ** exponents(sym_term_count, term_count, :, :))
        end do
    end do

    ! Interpolate between constant (short range) and U (medium and long range)
    U = interpolate(dist_matrix, U)

contains

! We need to convert lin_lin_coordinates (R, Theta_fixed, Theta_moving, Phi)
! to a distance matrix
subroutine lin_lin_to_distance_matrix(R, Theta_fixed, Theta_moving, Phi, dist_matrix)
    implicit none
    ! uses the global variables
    !  fixed_mol, num_fixed_atoms
    !  moving_mol, num_moving_atoms
    ! Parameters
    real(dp), intent(in) :: R, Theta_fixed, Theta_moving, Phi
    real(dp), intent(out) :: dist_matrix(num_fixed_atoms, num_moving_atoms)
    ! local variables
    real(dp) :: Theta_fixed_rad, Theta_moving_rad, Phi_rad
    real(dp) :: moved_position(num_moving_atoms, 4) 
    integer :: row, col

    Theta_fixed_rad = deg_to_rad(Theta_fixed)
    Theta_moving_rad = deg_to_rad(Theta_moving)
    Phi_rad = deg_to_rad(Phi)

    call homogenize_matrix(moving_mol, num_moving_atoms, moved_position)
    moved_position = matmul(moved_position, &
        transpose(reshape ( &
          (/ cos(Theta_moving_rad), sin(Theta_moving_rad), 0.0d0, 0d0, &
            -sin(Theta_moving_rad), cos(Theta_moving_rad), 0.0d0, 0d0, &
                             0.0d0,                 0.0d0, 1.0d0, 0d0, &
                             0.0d0,                 0.0d0, 0.0d0, 1d0 /), &
          (/ 4, 4 /) )) )

    moved_position = matmul(moved_position, &
        transpose(reshape ( &
        (/ 1d0, 0d0, 0d0, 0d0, &
           0d0, 1d0, 0d0, 0d0, &
           0d0, 0d0, 1d0, 0d0, &
             R, 0d0, 0d0, 1d0 /), &
        (/ 4, 4 /) )) )

    moved_position = matmul(moved_position, &
        transpose(reshape (                              &
          (/ 1.0d0,        0.0d0,        0.0d0,  0d0,    &
             0.0d0, cos(Phi_rad), -sin(Phi_rad), 0d0,    &
             0.0d0, sin(Phi_rad),  cos(Phi_rad), 0d0,    &
             0.0d0,        0.0d0,         0.0d0, 1d0 /), &
          (/ 4, 4 /) )) )

    moved_position = matmul(moved_position, &
        transpose(reshape ( &
          (/ cos(Theta_fixed_rad), -sin(Theta_fixed_rad), 0.0d0, 0.0d0,    &
             sin(Theta_fixed_rad),  cos(Theta_fixed_rad), 0.0d0, 0.0d0,    &
                            0.0d0,                  0.0d0, 1.0d0, 0.0d0,    &
                            0.0d0,                  0.0d0, 0.0d0, 1.0d0 /), &
          (/ 4, 4 /) )) )
    ! moved_position is now in the correct position
    ! we compute cross distances
    do row = 1, num_fixed_atoms
        do col = 1, num_moving_atoms
            dist_matrix(row, col) = sqrt( &
                sum( (fixed_mol(row, :) - moved_position(col, 1:3)) ** 2 ) )
        end do
    end do
end subroutine lin_lin_to_distance_matrix

! convert to radians
function deg_to_rad(angle_deg) result(angle_rad)
    implicit none
    real(dp) :: angle_rad
    real(dp), intent(in) :: angle_deg
    angle_rad = (pi / 180.0d0) * angle_deg
end function deg_to_rad

! homogenize
subroutine homogenize_matrix(A, num_col, A_hom)
    ! A must be num_col by 3 matrix
    implicit none
    real(dp), intent(in) :: A(num_col, 3)
    integer, intent(in) :: num_col
    real(dp), intent(out) :: A_hom(num_col, 4)
    A_hom(1:num_col, 1:3) = A
    A_hom(1:num_col, 4) = 1d0
end subroutine homogenize_matrix

! interpolate
function interpolate(dist_matrix, U) result(Umodified) 
    ! interpolate between a constant (short distances) and
    ! E (mid-and long range)
    implicit none
    real(dp) :: Umodified
    real(dp), intent(in), dimension(num_fixed_atoms, num_moving_atoms) :: dist_matrix
    real(dp), intent(in) :: U

    real(dp), parameter :: potential_zero_distance = 1d9 ! hartrees
    real(dp), parameter :: dist_short = 2.0 ! bohr

    real(dp) :: t 

    t = exp( &
        -product(dist_matrix) / &
        dist_short ** (num_fixed_atoms * num_moving_atoms) )
    Umodified = t * potential_zero_distance + (1d0 - t) * U
end function interpolate

end subroutine NNCOPot
