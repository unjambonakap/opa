#include <opa_asm.h>
    ENTRY(getsid)


    movl    $__NR_getsid, %eax
    int     $0x80


    ret
END(getsid)

