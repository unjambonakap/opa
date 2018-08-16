#include <opa_asm.h>
    ENTRY(__rt_sigtimedwait)


    movl    $__NR_rt_sigtimedwait, %eax
    int     $0x80


    ret
END(__rt_sigtimedwait)

