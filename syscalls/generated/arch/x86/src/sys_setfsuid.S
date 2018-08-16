#include <opa_asm.h>
    ENTRY(setfsuid)


    movl    $__NR_setfsuid, %eax
    int     $0x80


    ret
END(setfsuid)

