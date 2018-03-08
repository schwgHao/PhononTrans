module libbindC
!C++ functions declaration
    use iso_c_binding
    
    interface
        subroutine writensc(slabel, nsc) bind(c, name="writensc")
            use iso_c_binding
            implicit none
            character(c_char) :: slabel(*)
            integer(c_int) :: nsc(3)
        end subroutine
        
        subroutine vibrator_c(slabel, isbulk) bind(c, name="vibrator_c")
            use iso_c_binding
            implicit none
            character(c_char) :: slabel(*)
            logical(c_bool) :: isbulk
        end subroutine
    end interface

end module
