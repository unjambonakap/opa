#include <opa_asm.h>
    ENTRY(fdatasync)


    movl    $__NR_fdatasync, %eax
    syscall


    ret
END(fdatasync)

