#include <opa_asm.h>
    ENTRY(__openat)


    movl    $__NR_openat, %eax
    int     $0x80


    ret
END(__openat)

