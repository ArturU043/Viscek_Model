module var
  implicit none

  integer :: a, b, L=5
  integer , parameter :: N = 300
  real :: m, dt, vnorme=0.03 , eta = 0.1

end module var

! a sert a repérer la a-ième particules parmi les n (b1 = x et b2 = y pour une particule a données)
! L est la taille de la boîte
! r est la matrices des positions, initialement random, voir première boucle do in Viscek
! Au début, toutes les particules ont la même vitesse v dont l'opérateur doit pouvoir contrôler la norme
! Initialement, les particules ont des orientations aléatoires theta mais la vnorme doit être conservée, d'où la variable bla

program Viscek
  use var
  implicit none
  real ,dimension (N, 2) :: r , v
  real ,dimension (N) :: theta

  call init(r,v)

endprogram Viscek

function distance (x1,y1,x2,y2)
  use var
  implicit none
  real ,intent(in) :: x1,x2,y1,y2
  real :: d
  d=sqrt((x1-x2)**2+(y1-y2)**2)
endfunction distance

subroutine evolution (r,v,theta)
  use var
  implicit none

!_v siginifit voisin ; _n siginifit new list
  real ,dimension (N,2), intent(in) :: r
  real ,dimension (N), intent(in) :: theta
  real ,dimension (N), intent(out) :: theta_n
  integer :: n_v
  real :: thet_av , thet_v , D_thet , d

  do i= 1, N                                            ! on fixe la particule i
    n_vu=0                                              ! on initialise le nb de voisins , 0 avant la boucle j
    thet_v=0                                           ! de même pour la somme des angles voisins
    do j=1 , N
      if (i /= j) then                                 ! on balaie
        call distance(r(i,1),r(i,2),r(j,1),r(j,2))
      endif

      if (d <=1) then                                  ! i et j sont voisins
        n_v=n_v+1                                      ! ajoute 1 au nb de voisins
        thet_v=0+theta(j)
      endif                             ! somme des thetas voisins
    enddo

    thet_av=(thet_v+theta(i))/(n_v+1)                  ! moyenne des thetas en incluant la particule i
    D_thet= (eta/2) * (2*rand()-1)
    theta_n(i)=thet_av+D_thet
  enddo
end subroutine evolution


subroutine init(r,v)
  use var
  implicit none
  real ,dimension (N,2), intent(inout) :: r, v
  real ,dimension (N) :: theta

  open (11,file='/users/etu/3670788/PhysNum/viscek/Viscek_Model/positions/pos_0.txt')
  open (12,file='/users/etu/3670788/PhysNum/viscek/Viscek_Model/vitesses/vit_0.txt')
  open (13,file='/users/etu/3670788/PhysNum/viscek/Viscek_Model/angles/tet_0.txt')

   ! Initialisation des coposantes des positions et des vitesses 1:x 2:y
  do b = 1, N
    r(b,1) = L*rand()
    r(b,2) = L*rand()
    theta(b)= 2*3.151592*rand()
    v(b,1) = vnorme*cos(theta(b))
    v(b,2) = sqrt(vnorme**2-v(b,1)**2)
    write(11,*) r(b,1), r(b,2)
    write(12,*) v(b,1), v(b,2)
    write(13,*) theta(b)
  enddo
end subroutine init
