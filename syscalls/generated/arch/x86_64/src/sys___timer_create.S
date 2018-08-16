#include <opa_asm.h>
    ENTRY(__timer_create)


    movl    $__NR_timer_create, %eax
    syscall


    ret
END(__timer_create)

