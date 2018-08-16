#include <opa_asm.h>
    ENTRY(msync)


    movl    $__NR_msync, %eax
    int     $0x80


    ret
END(msync)

