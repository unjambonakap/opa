#include <opa_asm.h>
    ENTRY(listxattr)


    movl    $__NR_listxattr, %eax
    syscall


    ret
END(listxattr)

