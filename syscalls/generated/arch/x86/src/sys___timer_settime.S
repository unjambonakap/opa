#include <opa_asm.h>
    ENTRY(__timer_settime)


    movl    $__NR_timer_settime, %eax
    int     $0x80


    ret
END(__timer_settime)

