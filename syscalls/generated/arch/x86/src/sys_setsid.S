#include <opa_asm.h>
    ENTRY(setsid)


    movl    $__NR_setsid, %eax
    int     $0x80


    ret
END(setsid)

