#include <opa_asm.h>
    ENTRY(setgroups)


    movl    $__NR_setgroups32, %eax
    int     $0x80


    ret
END(setgroups)

