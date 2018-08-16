#include <opa_asm.h>
    ENTRY(___fchmod)


    movl    $__NR_fchmod, %eax
    int     $0x80


    ret
END(___fchmod)

