#include <opa_asm.h>
    ENTRY(eventfd)


    movl    $__NR_eventfd2, %eax
    syscall


    ret
END(eventfd)

