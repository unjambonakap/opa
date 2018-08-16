#include <opa_asm.h>
    ENTRY(__rt_sigsuspend)


    movl    $__NR_rt_sigsuspend, %eax
    syscall


    ret
END(__rt_sigsuspend)

