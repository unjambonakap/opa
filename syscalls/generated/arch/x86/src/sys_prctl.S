#include <opa_asm.h>
    ENTRY(prctl)


    movl    $__NR_prctl, %eax
    int     $0x80


    ret
END(prctl)

