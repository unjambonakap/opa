#include <opa_asm.h>
    ENTRY(__rt_sigsuspend)


    movl    $__NR_rt_sigsuspend, %eax
    int     $0x80


    ret
END(__rt_sigsuspend)

