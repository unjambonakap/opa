#include <opa_asm.h>
    ENTRY(utimensat)


    movl    $__NR_utimensat, %eax
    int     $0x80


    ret
END(utimensat)

