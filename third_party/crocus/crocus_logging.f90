! SPDX-License-Identifier: GPL-3.0-or-later
module chm_crocus_logging
    use iso_c_binding, only: c_char, c_int, c_size_t
    use iso_fortran_env, only: error_unit, output_unit
    implicit none
    private
    public :: cro_log, cro_values, cro_nl, cro_debug, cro_info, cro_warn, cro_error, cro_fatal

    integer(c_int), parameter :: cro_debug=0, cro_info=1, cro_warn=2, cro_error=3, cro_fatal=4
    character, parameter :: cro_nl = new_line('a')
    interface
        function chm_crocus_write_log(level, text, length) bind(C) result(handled)
            import c_int, c_char, c_size_t
            integer(c_int), value :: level
            character(kind=c_char), intent(in) :: text(*)
            integer(c_size_t), value :: length
            integer(c_int) :: handled
        end function
    end interface
contains
    subroutine cro_log(level, text)
        integer(c_int), intent(in) :: level
        character(kind=c_char, len=*), intent(in) :: text
        integer :: unit
        if (chm_crocus_write_log(level, text, len_trim(text, kind=c_size_t)) /= 0) return
        unit = output_unit
        if (level >= cro_warn) unit = error_unit
        write(unit, '(a)') trim(text)
        if (level == cro_fatal) flush(unit)
    end subroutine

    ! A full single-precision profile, including its signs and exponents.
    function cro_values(values) result(text)
        real, intent(in) :: values(:)
        character(len=17*size(values)) :: text
        if (size(values) > 0) write(text, '(*(1x,es16.8e3))') values
    end function
end module
