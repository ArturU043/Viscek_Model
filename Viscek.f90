module var
  implicit none

  integer :: i, k,j, a, b, L=25
  integer , parameter :: N = 300
  real :: vnorme=0.03 , eta=0.1 ,dt=1.0

end module var

! a sert a repérer la a-ième particules parmi les n (b1 = x et b2 = y pour une particule a données)
! L est la taille de la boîte
! r est la matrices des positions, initialement random, voir première boucle do in Viscek
! Au début, toutes les particules ont la même vitesse v dont l'opérateur doit pouvoir contrôler la norme
! _n = new

program Viscek
  use var
  implicit none

  real ,dimension (N, 2) :: r , v , r_n , v_n
  real ,dimension (N) :: theta , theta_n
  real :: v_ordre_x , v_ordre_y , v_ordre

  open (12,file='data/eta_v_ordre.dat')


  !do k=1, 1000
    v_ordre=0.0
    v_ordre_x=0.0
    v_ordre_y=0.0
    call init(r,v,theta)
    do a=1, 100                          ! nombre de pas de temps dt
      call evolution(r,theta,theta_n)       ! influence des voisins sur l'angle
      call deplacement(theta_n,r,r_n,v_n)
      r=r_n ; v=v_n ; theta=theta_n
      do i=1, N
        write(11,*) r(i,1), r(i,2) ,v(i,1), v(i,2)
        v_ordre_x = v_ordre_x + v(i,1)
        v_ordre_y = v_ordre_y + v(i,2)
      end do
      write(11,*)
      write(11,*)
    end do
    v_ordre_x = abs(v_ordre_x /(N*100*vnorme))
    v_ordre_y = abs(v_ordre_y /(N*100*vnorme))
    v_ordre= sqrt(v_ordre_x**2+v_ordre_y**2)
    write(12,*) eta, v_ordre
    write(*,*) v_ordre, v_ordre_x, v_ordre_y
!    eta = eta + 0.05
!   end do
endprogram Viscek



subroutine init(r,v,theta)
  use var
  implicit none

  real ,dimension (N,2), intent(inout) :: r, v
  real ,dimension (N), intent (inout):: theta

  open (11,file='data/pos_v.dat')
  ! Initialisation des coposantes des positions et des vitesses 1:x 2:y
  do b=1, N
    r(b,1) = L*rand()
    r(b,2) = L*rand()
    theta(b)= 2*3.151592654*rand()
    v(b,1) = vnorme*cos(theta(b))
    v(b,2) = vnorme*sin(theta(b))
    write(11,*) r(b,1), r(b,2) ,v(b,1), v(b,2)
  enddo
  write(11,*)
  write(11,*)
end subroutine init



subroutine evolution (r,theta,theta_n)
  use var
  implicit none
                                                  !_v siginifie voisin ; _n siginifie new list
  real ,dimension (N,2), intent(inout) :: r
  real ,dimension (N), intent(in) :: theta
  real ,dimension (N), intent(out) :: theta_n
  integer :: n_v
  real :: thet_av , thet_v , D_thet , d

  do i=1, N                                           ! on fixe la particule i
    n_v=0                                              ! on initialise le nb de voisins , 0 avant la boucle j
    thet_v=0                                             ! de même pour la somme des angles voisins
    do j=1 , N
      if (i /= j) then                                 ! on balaie
        call distance(r(i,1),r(i,2),r(j,1),r(j,2),d)
      endif
      if (d<=1) then                                  ! i et j sont voisins
        n_v=n_v+1                                      ! ajoute 1 au nb de voisins
        thet_v=thet_v+theta(j)
      endif                                            ! somme des thetas voisins
    enddo
    thet_av=(thet_v+theta(i))/(n_v+1)                  ! moyenne des thetas en incluant la particule i
    D_thet= (eta/2)*(2*rand()-1)
    theta_n(i)=thet_av+D_thet
  enddo
end subroutine evolution


subroutine distance (x1,y1,x2,y2,d)
  use var
  implicit none

  real ,intent(inout) :: x1,x2,y1,y2
  real :: x ,y                           ! Norme des composantes de la distance entre 2 particules
  real , intent(out):: d

  x=abs (x1 - x2)
  y=abs (y1 - y2)

  if (x>0.5*L) then
    x = x - anint(x/L)*L
  endif
  if (y>0.5*L) then
    y= y - anint(y/L)*L
  endif
  d=sqrt(x**2+y**2)
endsubroutine distance





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

    if (r_n(i,1)>L) then                              ! Si on depasse le cõté gauche de la boite (L) on revient à droite
      r_n(i,1) = r_n(i,1) - L
    endif
    if (r_n(i,1)<0) then                              ! Si on depasse le cõté droit de la boite (0) on revient à gauche
      r_n(i,1) = r_n(i,1) + L
    endif

   if (r_n(i,2)> L) then
    r_n(i,2) = r_n(i,2) - L
   endif
   if (r_n(i,2)<0) then                              ! Si on depasse le cõté droit de la boite (0) on revient à gauche
      r_n(i,1) = r_n(i,2) + L
   endif
  end do
endsubroutine deplacement
