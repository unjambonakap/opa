#include <opa_asm.h>
    ENTRY(kill)


    movl    $__NR_kill, %eax
    int     $0x80


    ret
END(kill)

