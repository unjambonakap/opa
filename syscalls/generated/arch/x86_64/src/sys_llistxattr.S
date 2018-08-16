#include <opa_asm.h>
    ENTRY(llistxattr)


    movl    $__NR_llistxattr, %eax
    syscall


    ret
END(llistxattr)

