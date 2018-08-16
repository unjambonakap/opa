#include <opa_asm.h>
    ENTRY(readv)


    movl    $__NR_readv, %eax
    int     $0x80


    ret
END(readv)

