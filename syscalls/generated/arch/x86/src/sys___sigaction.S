#include <opa_asm.h>
    ENTRY(__sigaction)


    movl    $__NR_sigaction, %eax
    int     $0x80


    ret
END(__sigaction)

