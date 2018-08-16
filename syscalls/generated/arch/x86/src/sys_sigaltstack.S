#include <opa_asm.h>
    ENTRY(sigaltstack)


    movl    $__NR_sigaltstack, %eax
    int     $0x80


    ret
END(sigaltstack)

