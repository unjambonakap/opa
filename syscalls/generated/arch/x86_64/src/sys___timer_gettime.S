#include <opa_asm.h>
    ENTRY(__timer_gettime)


    movl    $__NR_timer_gettime, %eax
    syscall


    ret
END(__timer_gettime)

