#include <opa_asm.h>
    ENTRY(geteuid)


    movl    $__NR_geteuid32, %eax
    int     $0x80


    ret
END(geteuid)

