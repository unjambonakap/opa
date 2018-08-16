#include <opa_asm.h>
    ENTRY(__rt_sigprocmask)


    movl    $__NR_rt_sigprocmask, %eax
    int     $0x80


    ret
END(__rt_sigprocmask)

