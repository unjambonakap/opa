#include <opa_asm.h>
    ENTRY(__rt_sigaction)


    movl    $__NR_rt_sigaction, %eax
    int     $0x80


    ret
END(__rt_sigaction)

