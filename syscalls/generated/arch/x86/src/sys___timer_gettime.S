#include <opa_asm.h>
    ENTRY(__timer_gettime)


    movl    $__NR_timer_gettime, %eax
    int     $0x80


    ret
END(__timer_gettime)

