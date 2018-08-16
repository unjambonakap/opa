#include <opa_asm.h>
    ENTRY(__rt_sigpending)


    movl    $__NR_rt_sigpending, %eax
    syscall


    ret
END(__rt_sigpending)

