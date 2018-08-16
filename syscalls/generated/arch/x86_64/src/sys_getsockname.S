#include <opa_asm.h>
    ENTRY(getsockname)


    movl    $__NR_getsockname, %eax
    syscall


    ret
END(getsockname)

