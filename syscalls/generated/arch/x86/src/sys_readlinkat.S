#include <opa_asm.h>
    ENTRY(readlinkat)


    movl    $__NR_readlinkat, %eax
    int     $0x80


    ret
END(readlinkat)

