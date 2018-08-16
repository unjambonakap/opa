#include <opa_asm.h>
    ENTRY(__timer_delete)


    movl    $__NR_timer_delete, %eax
    syscall


    ret
END(__timer_delete)

