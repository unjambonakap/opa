#include <opa_asm.h>
    ENTRY(vmsplice)


    movl    $__NR_vmsplice, %eax
    int     $0x80


    ret
END(vmsplice)

