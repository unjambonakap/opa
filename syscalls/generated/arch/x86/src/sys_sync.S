#include <opa_asm.h>
    ENTRY(sync)


    movl    $__NR_sync, %eax
    int     $0x80


    ret
END(sync)

