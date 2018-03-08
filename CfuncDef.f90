module libbindC
!C++ functions declaration
    use iso_c_binding
    
    interface
        subroutine writensc(slabel, nsc, isbulk, ia1, ia2) bind(c, name="writensc")
            use iso_c_binding
            implicit none
            character(c_char) :: slabel(*)
            integer(c_int) :: nsc(3)
            integer(c_int),value  :: ia1, ia2
            logical(c_bool),value :: isbulk
        end subroutine
        
        subroutine vibrator_c(slabel, isbulk, vibDecayRate, Enli, Enlf, &
        Nenl, Delta, bathTemp1, bathTemp2, Ntemp) bind(c, name="vibrator_c")
            use iso_c_binding
            implicit none
            character(c_char) :: slabel(*)
            logical(c_bool),value :: isbulk, vibDecayRate
            real(c_double),value :: Enli, Enlf, Delta, bathTemp1, bathTemp2
            integer(c_int),value :: Nenl, Ntemp
        end subroutine
    end interface

end module
