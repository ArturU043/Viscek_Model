module var
  implicit none

  integer :: a,b, n, L
  real, dimension (2,n) :: r, v, theta
  real :: m, dt

end module var

! a sert a repérer la a-ième particules parmi les n (b1 = x et b2 = y pour une particule a données)
! L est la taille de la boîte
! r est la matrices des positions, initialement random, voir première boucle do in Viscek
! Au début, toutes les particules ont la même vitesse v dont l'opérateur doit pouvoir contrôler la norme
! Initialement, les particules ont des orientations aléatoires theta mais la vnorme doit être conservée, d'où la variable bla

program Viscek
  implicit none
  use var

  real :: bla

  write (*,*) "Entrez un nombre entier de particules que vous voulez modéliser :"
  read (*,*) n
  write (*,*) "Entrez le nombre entier d'étapes (la durée totale qui sera discrétisée):"
  read (*,*) m
  write (*,*) "Entrez le pas de discrétisation du temps :"
  read (*,*) dt
  write (*,*) "Entrez la norme de la vitesse des particules :"
  read (*,*) vnorme

!Initialisation des coposantes des positions et des vitesses:

il y a un probleme dans le theta:

  do a = 1, n
    r(a,b) = L*rand()
    bla = rand()
    theta =
    v(a,b) = vnorme*bla
    do b = 1, 2
      r(a,b) = L*rand()
      v(a,b) = v*(1-bla)
    end do
  end do

  call subroutine positions

end program Viscek


! Calcul de l'évolution temporelles des positions:
subroutine positions
  implicit none
  use var

  integer :: i

  do i = 1, m
    do a=1, n
      r(a,b) =
      do b = 1, 2
        r(a,b) =
      end do
    end do
  end do
