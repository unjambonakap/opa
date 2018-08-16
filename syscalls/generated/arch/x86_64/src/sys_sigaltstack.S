#include <opa_asm.h>
    ENTRY(sigaltstack)


    movl    $__NR_sigaltstack, %eax
    syscall


    ret
END(sigaltstack)

