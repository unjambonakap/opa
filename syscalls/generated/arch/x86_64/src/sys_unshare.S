#include <opa_asm.h>
    ENTRY(unshare)


    movl    $__NR_unshare, %eax
    syscall


    ret
END(unshare)

