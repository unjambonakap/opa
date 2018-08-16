#include <opa_asm.h>
    ENTRY(llistxattr)


    movl    $__NR_llistxattr, %eax
    int     $0x80


    ret
END(llistxattr)

