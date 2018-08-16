#include <opa_asm.h>
    ENTRY(setitimer)


    movl    $__NR_setitimer, %eax
    int     $0x80


    ret
END(setitimer)

