#include <opa_asm.h>
    ENTRY(nanosleep)


    movl    $__NR_nanosleep, %eax
    int     $0x80


    ret
END(nanosleep)

