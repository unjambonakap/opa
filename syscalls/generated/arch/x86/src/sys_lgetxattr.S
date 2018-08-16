#include <opa_asm.h>
    ENTRY(lgetxattr)


    movl    $__NR_lgetxattr, %eax
    int     $0x80


    ret
END(lgetxattr)

