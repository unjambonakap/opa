#include <opa_asm.h>
    ENTRY(tee)


    movl    $__NR_tee, %eax
    int     $0x80


    ret
END(tee)

