#include <opa_asm.h>
    ENTRY(truncate64)


    movl    $__NR_truncate64, %eax
    int     $0x80


    ret
END(truncate64)

