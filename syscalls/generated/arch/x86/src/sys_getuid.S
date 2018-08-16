#include <opa_asm.h>
    ENTRY(getuid)


    movl    $__NR_getuid32, %eax
    int     $0x80


    ret
END(getuid)

