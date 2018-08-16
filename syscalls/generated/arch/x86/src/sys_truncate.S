#include <opa_asm.h>
    ENTRY(truncate)


    movl    $__NR_truncate, %eax
    int     $0x80


    ret
END(truncate)

