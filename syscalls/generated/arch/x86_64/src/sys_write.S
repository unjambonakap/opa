#include <opa_asm.h>
    ENTRY(write)


    movl    $__NR_write, %eax
    syscall


    ret
END(write)

