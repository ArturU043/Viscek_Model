module var
  implicit none

  integer :: i, j, a, b, L=5
  integer , parameter :: N = 300
  real :: m, vnorme=0.03 , eta = 0.1 ,dt=1.0

end module var

! a sert a repérer la a-ième particules parmi les n (b1 = x et b2 = y pour une particule a données)
! L est la taille de la boîte
! r est la matrices des positions, initialement random, voir première boucle do in Viscek
! Au début, toutes les particules ont la même vitesse v dont l'opérateur doit pouvoir contrôler la norme
! Initialement, les particules ont des orientations aléatoires theta mais la vnorme doit être conservée, d'où la variable bla

program Viscek
  use var
  implicit none
  real ,dimension (N, 2) :: r , v , r_n , v_n
  real ,dimension (N) :: theta , theta_n


  call init(r,v,theta)
  do a = 1 , 200                          ! 200 pas de temps dt
    call evolution(r,theta,theta_n)       ! influence des voisins sur l'angle
    call deplacement(theta_n,r,r_n,v_n)
    r=r_n ; v=v_n ; theta=theta_n

    do i=1,N
      write(11,*) r(i,1), r(i,2) ,v(i,1), v(i,2)
    end do
  write(11,*)
  write(11,*)
end do
endprogram Viscek




subroutine evolution (r,theta,theta_n)
  use var
  implicit none

                                                        !_v siginifit voisin ; _n siginifit new list
  real ,dimension (N,2), intent(in) :: r
  real ,dimension (N), intent(in) :: theta
  real ,dimension (N), intent(out) :: theta_n
  integer :: n_v
  real :: thet_av , thet_v , D_thet , d

  do i= 1, N                                           ! on fixe la particule i
    n_v=0                                              ! on initialise le nb de voisins , 0 avant la boucle j
    thet_v=0                                           ! de même pour la somme des angles voisins
    do j=1 , N
      if (i /= j) then                                 ! on balaie
        call distance(r(i,1),r(i,2),r(j,1),r(j,2),d)
      endif

      if (d <=1) then                                  ! i et j sont voisins
        n_v=n_v+1                                      ! ajoute 1 au nb de voisins
        thet_v=thet_v+theta(j)
      endif                                            ! somme des thetas voisins
    enddo

    thet_av=(thet_v+theta(i))/(n_v+1)                  ! moyenne des thetas en incluant la particule i
    D_thet= (eta/2) * (2*rand()-1)
    theta_n(i)=thet_av+D_thet
  enddo
end subroutine evolution



subroutine deplacement (theta_n,r,r_n,v_n)
  use var
  implicit none
  real ,dimension (N,2), intent(in) :: r
  real ,dimension (N), intent(in) :: theta_n
  real ,dimension (N,2), intent(out) :: r_n, v_n       ! Utilisation du nouvel angle pour calculer la nouvelle vitesse et position

  do i = 1,N
    v_n(i,1)=vnorme*cos(theta_n(i))
    v_n(i,2)=vnorme*sin(theta_n(i))
    r_n(i,1)=r(i,1)+v_n(i,1)*dt
    r_n(i,2)=r(i,2)+v_n(i,2)*dt
  end do

endsubroutine deplacement



subroutine distance (x1,y1,x2,y2,d)
  use var
  implicit none
  real ,intent(in) :: x1,x2,y1,y2                      ! Norme de la distance entre 2 particules
  real , intent(out):: d
  d=sqrt((x1-x2)**2+(y1-y2)**2)
endsubroutine distance



subroutine init(r,v,theta)
  use var
  implicit none
  real ,dimension (N,2), intent(inout) :: r, v
  real ,dimension (N), intent (inout):: theta

  open (11,file='/users/etu/3670788/PhysNum/viscek/Viscek_Model/positions/pos_v.txt')
  ! Initialisation des coposantes des positions et des vitesses 1:x 2:y
  do b = 1, N
    r(b,1) = L*rand()
    r(b,2) = L*rand()
    theta(b)= 2*3.151592*rand()
    v(b,1) = vnorme*cos(theta(b))
    v(b,2) = vnorme*sin(theta(b))
    write(11,*) r(b,1), r(b,2) ,v(b,1), v(b,2)
  enddo
  write(11,*)
  write(11,*)
end subroutine init
