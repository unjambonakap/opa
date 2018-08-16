#include <opa_asm.h>
    ENTRY(getrusage)


    movl    $__NR_getrusage, %eax
    int     $0x80


    ret
END(getrusage)

