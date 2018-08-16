#include <opa_asm.h>
    ENTRY(getitimer)


    movl    $__NR_getitimer, %eax
    int     $0x80


    ret
END(getitimer)

