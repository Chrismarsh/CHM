subroutine native_warning() bind(C)
    use chm_crocus_logging, only: cro_log, cro_warn
    call cro_log(cro_warn, 'test native warning')
end subroutine
