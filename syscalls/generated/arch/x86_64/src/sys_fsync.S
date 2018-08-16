#include <opa_asm.h>
    ENTRY(fsync)


    movl    $__NR_fsync, %eax
    syscall


    ret
END(fsync)

