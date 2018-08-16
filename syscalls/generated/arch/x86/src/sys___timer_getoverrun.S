#include <opa_asm.h>
    ENTRY(__timer_getoverrun)


    movl    $__NR_timer_getoverrun, %eax
    int     $0x80


    ret
END(__timer_getoverrun)

