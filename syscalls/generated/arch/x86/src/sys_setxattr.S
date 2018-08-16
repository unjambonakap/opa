#include <opa_asm.h>
    ENTRY(setxattr)


    movl    $__NR_setxattr, %eax
    int     $0x80


    ret
END(setxattr)

