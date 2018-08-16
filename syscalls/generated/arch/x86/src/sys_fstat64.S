#include <opa_asm.h>
    ENTRY(fstat64)


    movl    $__NR_fstat64, %eax
    int     $0x80


    ret
END(fstat64)

