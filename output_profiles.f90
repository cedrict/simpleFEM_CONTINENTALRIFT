!==============================================!
!
!==============================================!

subroutine output_profiles(np,nel,x,y,T,mat,density,viscosity,exx,eyy,exy,press,vpmode,icon)

implicit none

integer, intent(in) :: np,nel
real(8), dimension(np), intent(in) :: x,y,T
integer, dimension(nel), intent(in) :: mat
real(8), dimension(nel), intent(in) :: density,viscosity,exx,eyy,exy,press,vpmode 
integer, dimension(4,nel) :: icon

real(8) yc,E2,Tc
integer iel

!==============================================!

open(unit=123,file='OUT/profiles.dat')
write(123,*) '#1 y'
write(123,*) '#2 mat'
write(123,*) '#3 density'
write(123,*) '#4 viscosity'
write(123,*) '#5 pressure'
write(123,*) '#6 dev stress 2nd inv.'
write(123,*) '#7 Temperature'
write(123,*) '#8 vpmode'

do iel=1,nel

   yc=0.25*sum(y(icon(1:4,iel)))
   Tc=0.25*sum(T(icon(1:4,iel)))

   E2= sqrt( 0.5d0*(exx(iel)**2+eyy(iel)**2)+exy(iel)**2 )

   write(123,*) yc,mat(iel),density(iel),viscosity(iel),press(iel),2.d0*E2*viscosity(iel),Tc,vpmode(iel)

end do

close(123)

end subroutine

!==============================================!
