#include <opa_asm.h>
    ENTRY(truncate)


    movl    $__NR_truncate, %eax
    syscall


    ret
END(truncate)

