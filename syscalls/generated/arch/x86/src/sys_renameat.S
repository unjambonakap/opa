#include <opa_asm.h>
    ENTRY(renameat)


    movl    $__NR_renameat, %eax
    int     $0x80


    ret
END(renameat)

