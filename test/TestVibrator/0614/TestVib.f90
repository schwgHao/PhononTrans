program test

use iso_c_binding
use libbindC

implicit none

character slabel*20
real(c_double) :: Engi_c, Engf_c, delta_c !ios_binding_c
integer(c_int) :: Neng_c
logical (c_bool) :: isbulk
    
slabel='Cchain'
isbulk=.true.
Engi_c=0.0
Engf_c=1.0
Neng_c=200
delta_c=0.01

call vibrator_c(slabel, isbulk, Engi_c, Engf_c, Neng_c, delta_c)
!call vibrator_c('Cchain', .true., 0.0, 1.0, 200, 0.000001)

end
