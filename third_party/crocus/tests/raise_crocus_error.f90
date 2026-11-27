subroutine raise_crocus_error() bind(C)
    use modi_abor1_sfx
    call abor1_sfx('Crocus fatal-exit test')
end subroutine
