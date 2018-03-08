module libbindC
!C++ functions declaration
    use iso_c_binding
    
    interface
        subroutine writensc(slabel, nsc, isbulk) bind(c, name="writensc")
            use iso_c_binding
            implicit none
            character(c_char) :: slabel(*)
            integer(c_int) :: nsc(3)
            logical(c_bool) :: isbulk
        end subroutine
        
        subroutine vibrator_c(slabel, isbulk, Enli, Enlf, Nenl, Delta) bind(c, name="vibrator_c")
            use iso_c_binding
            implicit none
            character(c_char) :: slabel(*)
            logical(c_bool), value :: isbulk
            real(c_double), value :: Enli, Enlf, Delta
            integer(c_int), value :: Nenl
        end subroutine
    end interface

end module
