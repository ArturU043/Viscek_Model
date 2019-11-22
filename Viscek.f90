module var
  implicit none

  integer :: a, b, L=5
  integer , parameter :: N = 300
  real, dimension (2,N) :: r, v
  real :: m, dt, vnorme=0.03 , bruit = 0.1

end module var

! a sert a repérer la a-ième particules parmi les n (b1 = x et b2 = y pour une particule a données)
! L est la taille de la boîte
! r est la matrices des positions, initialement random, voir première boucle do in Viscek
! Au début, toutes les particules ont la même vitesse v dont l'opérateur doit pouvoir contrôler la norme
! Initialement, les particules ont des orientations aléatoires theta mais la vnorme doit être conservée, d'où la variable bla

program Viscek
  implicit none
  use var

  real ::o

  !write (*,*) "Entrez un nombre entier de particules que vous voulez modéliser :"
  !read (*,*) N
  !write (*,*) "Entrez le nombre entier d'étapes (la durée totale qui sera discrétisée):"
  !read (*,*) m
  !write (*,*) "Entrez le pas de discrétisation du temps :"
  !read (*,*) dt
  !write (*,*) "Entrez la norme de la vitesse des particules :"
  !read (*,*) vnorme

  ! Initialisation des coposantes des positions et des vitesses:
  do b = 1, N
    r(a,b) = L*rand()
    o = rand()
   
    do a = 1, 2
       r(a,b) = L*rand()
       v(a,b) = vnorme*(1-o)
    end do
  end do

 ! call subroutine positions

end program Viscek


! Calcul de l'évolution temporelles des positions:
subroutine positions
  implicit none
  use var

  integer :: i

 ! Calcul des vraies composantes de theta avec le paramêtre de controle

 ! Calcul des positions au cours du temps:
  do i = 1, m
    do a=1, n
      r(a,b) =
      do b = 1, 2
        r(a,b) =
      end do
    end do
  end do
