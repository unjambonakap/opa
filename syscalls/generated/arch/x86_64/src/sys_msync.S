#include <opa_asm.h>
    ENTRY(msync)


    movl    $__NR_msync, %eax
    syscall


    ret
END(msync)

