#include <opa_asm.h>
    ENTRY(statfs64)


    movl    $__NR_statfs, %eax
    syscall


    ret
END(statfs64)

