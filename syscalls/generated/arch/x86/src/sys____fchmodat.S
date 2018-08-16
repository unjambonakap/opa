#include <opa_asm.h>
    ENTRY(___fchmodat)


    movl    $__NR_fchmodat, %eax
    int     $0x80


    ret
END(___fchmodat)

