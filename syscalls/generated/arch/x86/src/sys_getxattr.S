#include <opa_asm.h>
    ENTRY(getxattr)


    movl    $__NR_getxattr, %eax
    int     $0x80


    ret
END(getxattr)

