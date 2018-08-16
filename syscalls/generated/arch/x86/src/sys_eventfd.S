#include <opa_asm.h>
    ENTRY(eventfd)


    movl    $__NR_eventfd2, %eax
    int     $0x80


    ret
END(eventfd)

