#include <opa_asm.h>
    ENTRY(flistxattr)


    movl    $__NR_flistxattr, %eax
    int     $0x80


    ret
END(flistxattr)

