#include <opa_asm.h>
    ENTRY(flistxattr)


    movl    $__NR_flistxattr, %eax
    syscall


    ret
END(flistxattr)

