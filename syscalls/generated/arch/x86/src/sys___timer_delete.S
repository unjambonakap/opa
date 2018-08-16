#include <opa_asm.h>
    ENTRY(__timer_delete)


    movl    $__NR_timer_delete, %eax
    int     $0x80


    ret
END(__timer_delete)

