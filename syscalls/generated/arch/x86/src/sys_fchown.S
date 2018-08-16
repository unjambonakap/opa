#include <opa_asm.h>
    ENTRY(fchown)


    movl    $__NR_fchown32, %eax
    int     $0x80


    ret
END(fchown)

