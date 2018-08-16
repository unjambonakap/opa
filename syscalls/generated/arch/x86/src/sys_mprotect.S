#include <opa_asm.h>
    ENTRY(mprotect)


    movl    $__NR_mprotect, %eax
    int     $0x80


    ret
END(mprotect)

