#include <opa_asm.h>
    ENTRY(pwrite64)


    movl    $__NR_pwrite64, %eax
    int     $0x80


    ret
END(pwrite64)

