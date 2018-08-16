#include <opa_asm.h>
    ENTRY(__rt_sigpending)


    movl    $__NR_rt_sigpending, %eax
    int     $0x80


    ret
END(__rt_sigpending)

