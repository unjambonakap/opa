#include <opa_asm.h>
    ENTRY(fallocate64)


    movl    $__NR_fallocate, %eax
    int     $0x80


    ret
END(fallocate64)

